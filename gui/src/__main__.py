import sys, os, math
from itertools import combinations

from PyQt5 import QtGui
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QPushButton, QLabel, QVBoxLayout, QFileDialog,
    QWidget, QListWidget, QFrame, QGroupBox, QMessageBox, QHBoxLayout, QLineEdit,
    QCheckBox, QSizePolicy, QProgressBar, QComboBox, QAbstractItemView, QGridLayout,
    QTableWidget, QTableWidgetItem
)
from PyQt5.QtCore import Qt, QThread, pyqtSignal, QTimer

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas, NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import networkx as nx


import physquirrel
from physquirrel import Split, SplitQuarnet, FourCycle, CycleQuarnet, QuartetSplitSet, DenseQuarnetSet
from physquirrel.utils import CircularOrdering

def longest_distance_to_root(G, root):
    # Initialize distances with -inf
    distances = {node: 0 for node in G.nodes}
    distances[root] = 0  # Distance to itself is 0

    # Iterate over all nodes to compute the longest path to the root
    for node in G.nodes:
        if node == root:
            continue  # Skip the root node, its distance is already set
        
        # Get all simple paths from the node to the root
        paths = list(nx.all_simple_paths(G, source=root, target=node))
        
        # Find the longest path
        longest_path_length = max((len(path) - 1 for path in paths))#, default=0)
        
        # Update the distance for this node
        distances[node] = longest_path_length

    return distances

def hierarchy_pos(G, root=None, width=5., vert_gap=0.5, xcenter=0.5):
    if root is None:
        root = next(iter(nx.topological_sort(G)))  # Start at a topologically sorted root if none is provided

    # Compute longest distances to the root
    distances = longest_distance_to_root(G, root)

    # Create a dictionary to store the level for each node based on distance
    level_mapping = {node: distances[node] for node in G.nodes}

    # Calculate the maximum distance for inversion
    max_distance = max(level_mapping.values())

    def _hierarchy_pos(G, node, width=2., xcenter=0.5, pos=None, parent=None):
        if pos is None:
            pos = {node: (xcenter, max_distance - level_mapping[node])}  # Invert y-coordinate
        else:
            pos[node] = (xcenter, max_distance - level_mapping[node])  # Invert y-coordinate
            
        neighbors = list(G.successors(node))  # Use successors for directed graph

        # Check if the parent exists in neighbors before attempting to remove
        if parent is not None and parent in neighbors:
            neighbors.remove(parent)

        if len(neighbors) != 0:
            dx = width / len(neighbors)  # Calculate space between children
            nextx = xcenter - width / 2 - dx / 2  # Start placing children

            for neighbor in neighbors:
                nextx += dx  # Move to the right for the next child
                pos = _hierarchy_pos(G, neighbor, width=dx, xcenter=nextx, pos=pos, parent=node)

        return pos

    return _hierarchy_pos(G, root, width=width, xcenter=xcenter)


class GUIMSA(physquirrel.msa.MSA):
    def __init__(self, arg):
        self.quarnets = None
        self.sequences = None
        if type(arg) is str:
            if arg.endswith(('.fa', '.fasta')):
                self.sequences = self._read_fasta(arg)
                self.taxa = set(self.sequences.keys())
            elif arg.endswith(('.nex', '.nexus')):
                self.sequences = self._read_nexus(arg)
                self.taxa = set(self.sequences.keys())
            elif arg.endswith(('.txt')):
                quarnet_set = DenseQuarnetSet(arg)
                self.quarnets = quarnet_set.quarnets
                self.taxa = set(quarnet_set.leaves)
            else:
                raise ValueError("Wrong filetype.")

        self.taxa_order = list(self.taxa)
        
        self.alg_thread = None

    def delta_heuristic(self, lam=0.3, weight=True, alg_thread=None):
        """Create a DenseQuarnetSet from the MSA using a similar method as the trinet
        creation algorithm from TriLoNet. The value lambda in [0,1] needs to be tuned with experiments.
        If weight is True, the quarnets will be weighted according to their eta-value."""
        
        if not self.quarnets is None:
            if weight:
                return GUIDenseQuarnetSet(self.quarnets)
            else:
                for q in self.quarnets:
                    q.set_weight(1.0)
                return GUIDenseQuarnetSet(self.quarnets)

        self.alg_thread = alg_thread
        self.alg_thread.update_status.emit("computing δ-values")

        distance = self._compute_distance_matrix()
        delta_sum = {taxum:0 for taxum in self.taxa}
        
        self.alg_thread.update_status.emit("generating quarnets")

        quarnets = []
        index_combinations = combinations(range(len(self.taxa)), 4)
        
        quarnet_numbers = math.comb(len(self.taxa), 4)
        quarnet_it = 0
        for comb in index_combinations:
            quarnet_it += 1
            prev_perc = 0
            
            i, j, k, l = list(comb)
            a = self.taxa_order[i]; b = self.taxa_order[j]; c = self.taxa_order[k]; d = self.taxa_order[l]

            split_distances = {
                Split({a,b},{c,d}): distance[i][j]+distance[k][l], 
                Split({a,c},{b,d}): distance[i][k]+distance[j][l], 
                Split({a,d},{b,c}): distance[i][l]+distance[j][k]
                }
            # sort the split distances in increasing order
            sorted_splits = sorted(split_distances.items(), key=lambda item: item[1])
            
            # compute the delta value
            if sorted_splits[0][1] == sorted_splits[1][1] and sorted_splits[1][1] == sorted_splits[2][1]:
                delta = 0
            else:
                numerator = sorted_splits[2][1] - sorted_splits[1][1]
                denominator = sorted_splits[2][1] - sorted_splits[0][1]
                delta = numerator / denominator
            
            delta_sum[a] += delta; delta_sum[b] += delta; delta_sum[c] += delta; delta_sum[d] += delta
            
            if delta < lam:
                best_split = sorted_splits[0][0]
                q = SplitQuarnet(best_split)
            else:
                wrong_split = sorted_splits[-1][0]
                if wrong_split == Split({a,b}, {c,d}):
                    circ = CircularOrdering([a,c,b,d])
                elif wrong_split == Split({a,c}, {b,d}):
                    circ = CircularOrdering([a,b,c,d])
                elif wrong_split == Split({a,d}, {b,c}):
                    circ = CircularOrdering([a,b,d,c])
                
                q = CycleQuarnet(circ)
            
            if weight == True:
                if delta < lam:
                    w = abs(lam - delta) / lam
                else:
                    w = abs(lam-delta) / (1 - lam)
                q.set_weight(w)
                
            quarnets.append(q)
            perc = int((quarnet_it/quarnet_numbers)*100)
            if perc > prev_perc + 1:
                #time.sleep(0.02)
                prev_perc = perc
                self.alg_thread.update_progress1.emit(int((quarnet_it/quarnet_numbers)*100))

        
        self.alg_thread.update_status.emit("assigning reticulations")

        # Assign reticulations
        # order taxa according to (mean) delta values
        reticulation_order = sorted(delta_sum, key=lambda x: delta_sum[x], reverse=True)

        ret_quarnets = []
        for q in quarnets:
            if isinstance(q, CycleQuarnet):
                circ = q.circular_order
                ret = next((element for element in reticulation_order if element in q.leaves), None)
                
                if weight == True:
                    w = q.weight
                    q = FourCycle(circ, ret, weight=w)
                else:
                    q = FourCycle(circ, ret)

            ret_quarnets.append(q)          
        
        self.alg_thread.update_progress1.emit(100)

        return GUIDenseQuarnetSet(ret_quarnets)

class GUIDenseQuarnetSet(physquirrel.DenseQuarnetSet):   
    def __init__(self, quarnets):
        super().__init__(quarnets)
        self.alg_thread = None
        
    def squirrel(self, outgroup=None, tsp_threshold=13, alg_thread=None, **kwargs):
        """Returns a semi-directed level-1 triangle-free network built up from the
        dense quarnet set. Repeatedly contracts the least supported edge in the
        quartet-joining tree, and then builds the network. If method="best", the
        network with the highest similaraity score is returned, if method="first",
        the last network for which the score did not drop is returned. 
            include_score: whether to return the similarity measure
            verbose: whether to print intermediate info
            visualize: whether to plot intermediate networks
            all_networks: whether to return all networks instead of the best one (ordered according to score)
            tsp_threshold: threshold up to where to solve tsp exactly
            kwargs are passed to the similarity function
        """
        
        self.alg_thread = alg_thread
        
        self.alg_thread.update_status.emit("generating blobtrees")

        tstar = self.tstar_tree()
        qj_tree = self.quartet_joining(starting_tree=tstar)
        split_supports = dict()
        for split in qj_tree.splits(include_trivial=False):
            split_supports[split] = self.split_support(split)
        sorted_splits = [k for k, v in sorted(split_supports.items(), key=lambda item: item[1])]
        
        self.alg_thread.update_status.emit("generating networks")

        tree0 = qj_tree
        quarnets0, qsplitsystem = tree0.quarnets(return_4splits=True)
        score0 = self.similarity(quarnets0, **kwargs)
        
        net_numbers = len(self.leaves) - 3
        
        trees = [tree0]; networks = [tree0]; scores = [score0]
        
        if outgroup is not None:
            net = networks[-1].to_directed_network(outgroup)
            self.alg_thread.intermediate_networks.emit(net, scores[-1])
        else:
            self.alg_thread.intermediate_networks.emit(networks[-1], scores[-1])

        self.alg_thread.update_progress2.emit(int((1/net_numbers)*100))

        for i, split in enumerate(sorted_splits):
            u, v = trees[-1].cutedge_from_split(split)
            split_induced_4splits = trees[-1].quartetsplits_from_cutedge(u, v)

            tree_new = trees[-1].copy()
            tree_new.contract_split(split)
            net_new = self.reconstruct_network_from_tree(tree_new, outgroup=outgroup, tsp_threshold=tsp_threshold)
            
            qsplitsystem = QuartetSplitSet({s for s in qsplitsystem if s not in split_induced_4splits})
            quarnets_new = net_new.quarnets(induced_4splits=qsplitsystem)

            score_new = self.similarity(quarnets_new, **kwargs)

            trees.append(tree_new); networks.append(net_new); scores.append(score_new)
            #time.sleep(0.75)
            self.alg_thread.update_progress2.emit(int(((i+2)/net_numbers)*100))


            if outgroup is not None:
                net = networks[-1].to_directed_network(outgroup)
                self.alg_thread.intermediate_networks.emit(net, scores[-1])
            else:
                self.alg_thread.intermediate_networks.emit(networks[-1], scores[-1])


        paired_lists = list(zip(scores, networks))
        paired_lists.sort(reverse=True, key=lambda x: x[0])
        sorted_scores, sorted_networks = zip(*paired_lists)
        
        self.alg_thread.update_progress2.emit(100)

        return list(sorted_networks), list(sorted_scores)
            

class AlgorithmWorker(QThread):
    # Define signals to communicate with the main thread
    update_progress1 = pyqtSignal(int)  # Signal for updating a progress bar
    update_progress2 = pyqtSignal(int)  # Signal for updating a progress bar
    algorithm_finished = pyqtSignal()  # Signal when the algorithm is finished
    update_status = pyqtSignal(str)
    intermediate_networks = pyqtSignal(object, object)


    def __init__(self, msa_frame, plot_frame, alg_frame, lambd, weights, selected_leaf, tsp):
        super().__init__()
        self.msa_frame = msa_frame
        self.plot_frame = plot_frame
        self.alg_frame = alg_frame
        self.lambd = lambd
        self.weights = weights
        self.selected_leaf = selected_leaf
        self.tsp = tsp

    def run(self):
        self.update_status.emit("initializing algorithm")

        # Run the algorithm
        DQ = self.msa_frame.msa.delta_heuristic(lam=self.lambd, weight=self.weights, alg_thread=self)

        DQ.squirrel(triangles=False, outgroup=self.selected_leaf, tsp_threshold=self.tsp, alg_thread=self)
        
        # Emit the signal when the algorithm is done
        self.algorithm_finished.emit()



class AlgorithmFrame(QFrame):
    def __init__(self, parent=None, msa_frame=None, plot_frame=None):
        super().__init__(parent)
        self.initUI()
        self.worker = None  # Store the worker thread
        self.msa_frame = msa_frame
        self.plot_frame = plot_frame
        
        self.dot_count = 0


    def initUI(self):
        alg_layout = QVBoxLayout(self)
                
        # Create the "Run Algorithm" button and disable it initially
        self.run_button = QPushButton("Reconstruct networks", self)
        self.run_button.setEnabled(False)
        self.run_button.clicked.connect(self.run_algorithm)
        
        alg_layout.addWidget(self.run_button)
        
        # Create a QGroupBox
        group_box = QGroupBox("Algorithm settings", self)
        group_box.setStyleSheet("QGroupBox { font-weight: bold; }")

        # Create a grid layout for the group box
        grid_layout = QGridLayout(group_box)
        
        # Leaf list box section
        self.leaf_combobox_label = QLabel("outgroup:", self)
        self.leaf_combobox_label.setToolTip("Outgroup to root the network. If None, the network is semi-directed.")
        self.leaf_combobox_label.setStyleSheet("color: grey;")
        self.leaf_combobox = QComboBox(self)
        self.leaf_combobox_label.setAlignment(Qt.AlignRight | Qt.AlignVCenter)  # Align right
        self.leaf_combobox.setStyleSheet("QComboBox {combobox-popup: 0;}")
        self.leaf_combobox.setMaxVisibleItems(10)
        # Add to grid layout
        grid_layout.addWidget(self.leaf_combobox_label, 0, 0)  # Row 0, Column 0
        grid_layout.addWidget(self.leaf_combobox, 0, 1)          # Row 0, Column 1
        
        # Lambda value section
        self.lambda_label = QLabel("lambda:", self)
        self.lambda_label.setToolTip("lambda-value to determine if a quarnet has a split.")
        self.lambda_label.setStyleSheet("color: grey;")
        self.lambda_input = QLineEdit(self)
        self.lambda_input.setPlaceholderText("default=0.3")
        self.lambda_input.setText("0.3")  # Set default value
        self.lambda_input.setEnabled(False)
        self.lambda_input.textChanged.connect(self.check_input)
        self.lambda_label.setAlignment(Qt.AlignRight | Qt.AlignVCenter)  # Align right

        # Add to grid layout
        grid_layout.addWidget(self.lambda_label, 1, 0)           # Row 1, Column 0
        grid_layout.addWidget(self.lambda_input, 1, 1)            # Row 1, Column 1
        
        # TSP threshold section
        self.tsp_label = QLabel("TSP threshold:", self)
        self.tsp_label.setToolTip("Upper bound on the number of leaves for which TSP is solved exactly.")
        self.tsp_label.setStyleSheet("color: grey;")
        self.tsp_input = QLineEdit(self)
        self.tsp_input.setPlaceholderText("default=13")
        self.tsp_input.setText("13")  # Set default value
        self.tsp_input.setEnabled(False)
        self.tsp_input.textChanged.connect(self.check_input)
        self.tsp_label.setAlignment(Qt.AlignRight | Qt.AlignVCenter)  # Align right

        # Add to grid layout
        grid_layout.addWidget(self.tsp_label, 2, 0)              # Row 2, Column 0
        grid_layout.addWidget(self.tsp_input, 2, 1)               # Row 2, Column 1
        
        # Checkbox for weights
        self.weights_checkbox = QCheckBox("Weighted quarnets", self)
        self.weights_checkbox.setToolTip("Whether to use weights for the quarnets.")
        self.weights_checkbox.setChecked(True)  
        self.weights_checkbox.setEnabled(False)
        
        # Add to grid layout
        grid_layout.addWidget(self.weights_checkbox, 3, 0, 1, 2)  # Row 3, Column 0 (spanning 2 columns)
        
        # Set the layout for the group box
        group_box.setLayout(grid_layout)
        
        # Finally, add the group box to your main layout (e.g., alg_layout)
        alg_layout.addWidget(group_box)

        # Create the group box for progress bars and status
        self.progress_groupbox = QGroupBox("Algorithm status", self)
        self.progress_groupbox.setStyleSheet("QGroupBox { font-weight: bold; }")
        self.progress_layout = QGridLayout(self.progress_groupbox)


        self.delta_label = QLabel("δ-heuristic:", self)
        self.progress_bar1 = QProgressBar(self)
        self.progress_bar1.setMinimum(0)
        self.progress_bar1.setMaximum(100)
        self.progress_bar1.setValue(0)  # Start at 0
        self.delta_label.setStyleSheet("color: grey;")


        # Add to grid layout
        self.progress_layout.addWidget(self.delta_label, 0, 0, alignment=Qt.AlignRight)
        self.progress_layout.addWidget(self.progress_bar1, 0, 1)

        # SQuIrReL progress bar
        self.squirrel_label = QLabel("Squirrel:", self)
        self.progress_bar2 = QProgressBar(self)
        self.progress_bar2.setMinimum(0)
        self.progress_bar2.setMaximum(100)
        self.progress_bar2.setValue(0)  # Start at 0
        self.squirrel_label.setStyleSheet("color: grey;")

        # Add to grid layout
        self.progress_layout.addWidget(self.squirrel_label, 1, 0, alignment=Qt.AlignRight)
        self.progress_layout.addWidget(self.progress_bar2, 1, 1)

        # Status label
        self.status_name = QLabel("status:", self)
        self.status_label = QLabel("algorithm not started", self)
        self.status_label.setStyleSheet("color: grey;")

        # Add to grid layout
        self.progress_layout.addWidget(self.status_name, 2, 0, alignment=Qt.AlignRight)
        self.progress_layout.addWidget(self.status_label, 2, 1)

        alg_layout.addWidget(self.progress_groupbox)
        


    def check_input(self, file_input='sequences'):
        """File input is either 'sequences' or 'quarnets'."""
        if file_input == 'sequences':
            self.lambda_input.setEnabled(True)
            self.lambda_label.setStyleSheet("color: black;")
            self.delta_label.setStyleSheet("color: black;")
            self.squirrel_label.setStyleSheet("color: black;")
        elif file_input == 'quarnets':
            self.lambda_input.setEnabled(False)
            self.lambda_label.setStyleSheet("color: grey;")
            self.delta_label.setStyleSheet("color: grey;")
            self.squirrel_label.setStyleSheet("color: black;")

        self.weights_checkbox.setEnabled(True)
        self.tsp_input.setEnabled(True)
        
        self.tsp_label.setStyleSheet("color: black;")
        self.leaf_combobox_label.setStyleSheet("color: black;")


        try:
            lambda_value = float(self.lambda_input.text())
            valid = 0 <= lambda_value <= 1
                
            tsp_value = int(self.tsp_input.text())

            valid = 0 <= lambda_value <= 1
            self.run_button.setEnabled(valid and self.msa_frame.selected_file_path is not None)

        except ValueError:
            self.run_button.setEnabled(False)

    def run_algorithm(self):
        
        self.plot_frame.clear_networks()
        
        self.plot_frame.is_directed = bool(self.leaf_combobox.currentText() != "-- None --")

        self.progress_bar1.setValue(0)
        self.progress_bar2.setValue(0)
        self.update_status_text("starting algorithm")
        
        lambd = float(self.lambda_input.text())
        weights = self.weights_checkbox.isChecked()
        tsp = int(self.tsp_input.text())

        selected_leaf = self.leaf_combobox.currentText()
        if selected_leaf == "-- None --":
            selected_leaf = None
        
        # Create and start the worker thread
        self.worker = AlgorithmWorker(self.msa_frame, self.plot_frame, self, lambd, weights, selected_leaf, tsp)
        self.worker.update_progress1.connect(self.progress_bar1.setValue)  # Connect progress updates to the progress bar
        self.worker.update_progress2.connect(self.progress_bar2.setValue)  # Connect progress updates to the progress bar
        self.worker.update_status.connect(self.update_status_text)
        self.worker.intermediate_networks.connect(self.update_networks)
        
        self.worker.algorithm_finished.connect(self.on_algorithm_finished)
        self.worker.start()

    def update_status_text(self, text):
        self.start_dot_animation(self.status_label, text)

    def update_networks(self, network, score):
        self.plot_frame.add_network(network, score)
    
    def on_algorithm_finished(self):
        
        # Stop dot animation
        self.stop_dot_animation()
        
        # Update the label text and make it grey
        self.status_label.setText("algorithm finished")

        # Update the plot frame with the new networks and scores
        #self.plot_frame.is_directed = bool(self.msa_frame.leaf_combobox.currentItem().text() != "-- None --")
        #self.plot_frame.scores = scores
        #self.plot_frame.networks = nets

        self.plot_frame.save_enewick_button.setEnabled(True)
        self.plot_frame.save_all_enewicks_button.setEnabled(True)

        # Update the plot frame with the new networks and scores
        best_score_index = self.plot_frame.scores.index(max(self.plot_frame.scores))
        self.plot_frame.network_table_widget.setCurrentCell(best_score_index, 0)
        
    def set_msa_frame(self, msa_frame):
        self.msa_frame = msa_frame
        
    def start_dot_animation(self, label, base_text):
        # Stop any previous animation if active
        self.stop_dot_animation()

        # Create a timer to animate the dots with any base text
        self.dot_count = 0  # Reset the dot count for new text
        self.dot_timer = QTimer(self)
        self.dot_timer.timeout.connect(lambda: self.update_dots(label, base_text))
        self.dot_timer.start(500)  # 500ms interval

    def stop_dot_animation(self):
        try:
            self.dot_timer.stop()
            self.dot_timer.deleteLater()  # Safely delete the timer
            self.dot_timer = None  # Reset the reference
        except:
            NameError

    def update_dots(self, label, base_text):
        # Generalized function to update the dots for any label and base text
        dots = '.' * ((self.dot_count % 3) + 1)  # Cycles through ".", "..", "..."
        label.setText(f"{base_text}{dots}")
        self.dot_count += 1

    def populate_leaves_combobox(self):
        # Clear and populate the leaves combobox
        self.leaf_combobox.clear()
        self.leaf_combobox.addItem("-- None --")  # Add "None" option as the first item
        leaves = self.msa_frame.msa.taxa
        self.leaf_combobox.addItems(sorted(list(leaves)))

        # Set "None" as the selected item by default
        self.leaf_combobox.setCurrentIndex(0)  # Select the "None" option




class MSAFrame(QFrame):
    def __init__(self, parent=None, alg_frame=None, main_window=None):
        super().__init__(parent)
        self.selected_file_path = None
        self.file_name = None
        self.msa = None
        self.alg_frame = alg_frame
        self.main_window = main_window
        self.initUI()

    def initUI(self):
        main_layout = QVBoxLayout(self)
        
        button_layout = QHBoxLayout(self)
        self.select_file_button = QPushButton("Select File", self)
        self.select_file_button.clicked.connect(self.open_file)
        button_layout.addWidget(self.select_file_button)
        
        self.clear_button = QPushButton("Clear all", self)
        self.clear_button.clicked.connect(self.clear_all)
        self.clear_button.setToolTip("Clear the whole application.")

        button_layout.addWidget(self.clear_button)
        
        main_layout.addLayout(button_layout)

        # MSA information section
        self.msa_info_frame = QGroupBox("File information", self)
        self.msa_info_frame.setStyleSheet("QGroupBox { font-weight: bold; }")
        self.msa_info_layout = QGridLayout(self.msa_info_frame)
        self.msa_info_layout.setColumnStretch(1, 1)  # Stretch the second column

        main_layout.addWidget(self.msa_info_frame)
        self.msa_info_layout.addWidget(QLabel("No file selected", self), 0, 0)


    def open_file(self, clear=True):
        if clear: self.clear_all()
        self.alg_frame.progress_bar1.setValue(0)  # Start at 0
        self.alg_frame.progress_bar2.setValue(0)  # Start at 0

        
        file_path, _ = QFileDialog.getOpenFileName(
            self.parent(), "Select a File", "", "All Files (*);;Fasta Files (*.fasta);;Nexus Files (*.nexus);;Text Files (*.txt)"
        )
        if file_path:
            self.selected_file_path = file_path
            self.file_name = os.path.basename(file_path)
            if self.file_name.endswith(('.fasta', '.fa', '.nexus', '.nex')):
                self.msa = GUIMSA(file_path)
                self.display_msa_info(file_input = 'sequences')
                self.alg_frame.populate_leaves_combobox()
                self.alg_frame.check_input(file_input = 'sequences')  # Check if the run button can be enabled
            elif self.file_name.endswith(('.txt')):
                self.msa = GUIMSA(file_path)
                self.display_msa_info(file_input = 'quarnets')
                self.alg_frame.populate_leaves_combobox()
                self.alg_frame.check_input(file_input = 'quarnets')  # Check if the run button can be enabled

            else:
                QMessageBox.warning(self.parent(), "Invalid filetype", "Please select a valid file.")
        else:
            QMessageBox.warning(self.parent(), "No File", "Please select a valid file.")
    
    def clear_all(self):
        self.main_window.initUI(fileopen=True)
            
    def display_msa_info(self, file_input='sequences'):
        """Displays the info box. file_input is either 'sequences' or 'quarnets'."""
        
        if file_input == 'sequences':
            sequence_length = self.msa.seq_length()
            num_taxa = len(self.msa)
            
        elif file_input == 'quarnets':
            quarnet_nr = len(self.msa.quarnets)
            num_taxa = len(self.msa.taxa)
            
        # Clear previous MSA info
        for i in reversed(range(self.msa_info_layout.count())): 
            self.msa_info_layout.itemAt(i).widget().deleteLater()
        
        # Add MSA info to the grid layout
        type_label = QLabel("input type:", self.msa_info_frame)
        type_label.setAlignment(Qt.AlignRight | Qt.AlignVCenter)  # Right align
        self.msa_info_layout.addWidget(type_label, 0, 0)
        if file_input == 'sequences':
            type_value = QLabel("<b>sequence alignment</b>", self.msa_info_frame)
        elif file_input == 'quarnets':
            type_value = QLabel("<b>dense quarnetset</b>", self.msa_info_frame)
        self.msa_info_layout.addWidget(type_value, 0, 1)
        
        
        filename_label = QLabel("filename:", self.msa_info_frame)
        filename_label.setAlignment(Qt.AlignRight | Qt.AlignVCenter)  # Right align
        self.msa_info_layout.addWidget(filename_label, 1, 0)
        
        filename_value = QLabel(f"<b>{self.file_name}</b>", self.msa_info_frame)
        self.msa_info_layout.addWidget(filename_value, 1, 1)
        
        if file_input == 'sequences':
            sequence_length_label = QLabel("seq. length:", self.msa_info_frame)
            sequence_length_label.setAlignment(Qt.AlignRight | Qt.AlignVCenter)  # Right align
            self.msa_info_layout.addWidget(sequence_length_label, 2, 0)
        
            sequence_length_value = QLabel(f"<b>{sequence_length:,} bp</b>", self.msa_info_frame)
            self.msa_info_layout.addWidget(sequence_length_value, 2, 1)
        
        elif file_input == 'quarnets':
            sequence_length_label = QLabel("nr. of quarnets:", self.msa_info_frame)
            sequence_length_label.setAlignment(Qt.AlignRight | Qt.AlignVCenter)  # Right align
            self.msa_info_layout.addWidget(sequence_length_label, 2, 0)
        
            sequence_length_value = QLabel(f"<b>{quarnet_nr:,}</b>", self.msa_info_frame)
            self.msa_info_layout.addWidget(sequence_length_value, 2, 1)
        
        
        num_taxa_label = QLabel("nr. of taxa:", self.msa_info_frame)
        num_taxa_label.setAlignment(Qt.AlignRight | Qt.AlignVCenter)  # Right align
        self.msa_info_layout.addWidget(num_taxa_label, 3, 0)
        
        num_taxa_value = QLabel(f"<b>{num_taxa}</b>", self.msa_info_frame)
        self.msa_info_layout.addWidget(num_taxa_value, 3, 1)
        
        taxa_list_label = QLabel("taxa list:", self.msa_info_frame)
        taxa_list_label.setAlignment(Qt.AlignRight | Qt.AlignHCenter)  # Right align
        self.msa_info_layout.addWidget(taxa_list_label, 4, 0)
        
        # If you want to add the list box for taxa
        taxa_list = sorted(list(self.msa.taxa))
        
        self.taxa_listbox = QListWidget(self.msa_info_frame)
        self.taxa_listbox.setSelectionMode(QAbstractItemView.NoSelection)  # No selection allowed
        self.taxa_listbox.addItems(taxa_list)
        self.msa_info_layout.addWidget(self.taxa_listbox, 4, 1)  # Add the list box to the layout


    

class PlotFrame(QFrame):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.networks = []  # Store networks after running the algorithm
        self.scores = []    # Store scores corresponding to networks
        self.enewicks = []
        self.is_directed = None
        self.initUI()

    def initUI(self):
        # Main horizontal layout
        self.main_layout = QHBoxLayout(self)  # Main layout that splits left and right sections
        
        # *** LEFT SECTION (Network Table and Save Buttons) ***
        
        # Left layout (vertical) to contain the network list and save buttons
        left_layout = QVBoxLayout()
        
        # Title for the network list section
        network_list_title = QLabel("<b>Reconstructed networks<b>", self)
        left_layout.addWidget(network_list_title)
        
        # Create a QTableWidget to show networks and scores in a grid
        self.network_table_widget = QTableWidget(self)
        self.network_table_widget.setColumnCount(2)
        self.network_table_widget.setHorizontalHeaderLabels(['Network', 'Score'])
        self.network_table_widget.currentCellChanged.connect(self.on_network_selected)
        self.network_table_widget.setSortingEnabled(True)  # Enable sorting
        self.network_table_widget.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.network_table_widget.setFixedWidth(220)  # You can adjust this value to your needs

        # Add the network table to the left layout
        left_layout.addWidget(self.network_table_widget)
        
        # Add save buttons below the table
        self.save_enewick_button = QPushButton("Save selected network", self)
        self.save_enewick_button.setToolTip("Save the selected network as a eNewick string.")

        self.save_all_enewicks_button = QPushButton("Save all networks", self)
        self.save_all_enewicks_button.setToolTip("Save all networks (in the current order) as a eNewick string.")

        self.save_enewick_button.clicked.connect(self.save_enewick_to_file)
        self.save_all_enewicks_button.clicked.connect(self.save_all_enewicks_to_file)
        
        # Initially disable the buttons
        self.save_enewick_button.setEnabled(False)
        self.save_all_enewicks_button.setEnabled(False)
        
        # Create a QGroupBox for the save buttons
        save_group_box = QGroupBox("eNewick save options:", self)
        save_group_box.setStyleSheet("QGroupBox { font-weight: bold; }")

        # Create a layout for the group box
        save_layout = QVBoxLayout()
        
        # Add the save buttons to the group box layout
        save_layout.addWidget(self.save_enewick_button)
        save_layout.addWidget(self.save_all_enewicks_button)
        
        # Set the layout to the group box
        save_group_box.setLayout(save_layout)
        
        # Add the group box to the layout below the table
        left_layout.addWidget(save_group_box)

        
        # Add left layout to the main layout
        self.main_layout.addLayout(left_layout)
        
        # *** RIGHT SECTION (Plot and ENewick String) ***
        
        # Add a vertical separator (QFrame) between left and right sections
        self.main_layout.addSpacing(10)  # Add 10 pixels of space after the separator
        separator = QFrame()
        separator.setFrameShape(QFrame.VLine)  # Create a vertical line
        separator.setFrameShadow(QFrame.Sunken)  # Optional: make it appear sunken
        self.main_layout.addWidget(separator)
        self.main_layout.addSpacing(10)  # Add 10 pixels of space after the separator

        # Right layout (vertical) for plot and ENewick display
        right_layout = QVBoxLayout()
        
        # Title for the plot section
        plot_title = QLabel("<b>Network visualizer<b>", self)
        right_layout.addWidget(plot_title)
        
        # Horizontal layout for the ENewick label and line edit (displayed below the plot)
        enewick_layout = QHBoxLayout()
        
        # ENewick label
        self.enewick_label = QLabel("eNewick:", self)
        enewick_layout.addWidget(self.enewick_label)
        
        # ENewick string display as a selectable QLineEdit
        self.enewick_line_edit = QLineEdit(self)
        self.enewick_line_edit.setReadOnly(True)  # Set to read-only to allow text selection
        enewick_layout.addWidget(self.enewick_line_edit)
        
        # Add the ENewick layout below the plot
        right_layout.addLayout(enewick_layout)
        
        # Add warning for semi-directed network (displayed below the ENewick string)
        self.enewick_warning = QLabel("", self)
        self.enewick_warning.setStyleSheet("color: grey;")
        self.enewick_warning.setWordWrap(True)
        right_layout.addWidget(self.enewick_warning)
        
        # Matplotlib canvas for plotting
        self.figure = Figure()
        self.canvas = FigureCanvas(self.figure)
        
        # Add the canvas to the right layout
        right_layout.addWidget(self.canvas)
        
        # Navigation toolbar for interactivity
        self.toolbar = NavigationToolbar(self.canvas, self)
        right_layout.addWidget(self.toolbar)
        
        
        # Set size policy to make the canvas fill the available space
        self.canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.canvas.updateGeometry()  # Ensure it updates when resized
        
        # Add the right layout to the main layout
        self.main_layout.addLayout(right_layout)
        
        self.plot_network()

    def add_network(self, network, score):
        self.networks.append(network)
        self.scores.append(score)
        i = len(self.networks)

        # Add a new row to the table
        row_position = self.network_table_widget.rowCount()
        self.network_table_widget.insertRow(row_position)
        self.network_table_widget.setItem(row_position, 0, QTableWidgetItem(f"N{i}"))
        self.network_table_widget.setItem(row_position, 1, QTableWidgetItem(f"{score:.5f}"))

        if isinstance(network, physquirrel.SemiDirectedNetwork):
            # Display "NO enwick" in grey if the network is SemiDirected
            for root in network.leaves:
                if network.is_rootable_at(root):
                    rooted_network = network.to_directed_network(root)
                    enewick_string = rooted_network.create_enewick()  
                    break
        else:
            enewick_string = network.create_enewick()

        self.enewicks.append(enewick_string)

        if self.networks:
            self.network_table_widget.setCurrentCell(row_position, 0)  # Select the latest added network
            self.on_network_selected(row_position, 0)  # Call to display the network and its ENewick string

    def clear_networks(self):
        self.network_table_widget.clearContents()
        self.network_table_widget.setRowCount(0)
        self.networks = []
        self.scores = []
        self.plot_network()
        self.enewick_line_edit.setText("")

    def on_network_selected(self, current_row, current_column):
        """Handle network selection from the table and update the plot and ENewick string."""
        # Get the selected network name from column 1
        selected_network_name = self.network_table_widget.item(current_row, 0).text() # Assuming networks is a list of lists or tuples
        # Extract the number from the network name in the format N{xxxx}
        number_str = selected_network_name[1]  # Extract 'xxxx'
        network_index = int(number_str)-1  # Convert to an integer index
        
        # Check if the index is within the valid range
        selected_network = self.networks[network_index]  # Get network using the index
        self.plot_network(selected_network)

        self.enewick_line_edit.setText(self.enewicks[network_index])
        self.enewick_line_edit.setStyleSheet("color: black;")  # Reset text color to black
        
        # Update the ENewick string display
        if isinstance(selected_network, physquirrel.SemiDirectedNetwork):
            self.enewick_warning.setText("*warning*: an arbitrary rooting is used for the eNewick string of a semi-directed network.")
        else:
            self.enewick_warning.setText("")


    def save_enewick_to_file(self):
        """Save the currently selected ENewick string to a text file."""
        current_row = self.network_table_widget.currentRow()
        if current_row == -1:
            QMessageBox.warning(self, "Warning", "No network selected!")
            return
        
        # Get the ENewick string of the selected network
        enewick_string = self.enewicks[current_row]  # Adjust to your method of getting the ENewick string
    
        # Open a file dialog to select the save location
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        file_name, _ = QFileDialog.getSaveFileName(self, "Save ENewick", "", "Text Files (*.txt);;All Files (*)", options=options)
    
        if file_name:
            # Check if the file_name has an extension
            if '.' not in file_name.split('/')[-1]:
                file_name += '.txt'
            with open(file_name, 'w') as file:
                file.write(enewick_string)
            QMessageBox.information(self, "Success", "ENewick string saved successfully!")
    
    def save_all_enewicks_to_file(self):
        """Save all ENewick strings of the networks to a text file."""
        if not self.networks:
            QMessageBox.warning(self, "Warning", "No networks available to save!")
            return

        current_order = []
        for row in range(0, len(self.networks)):
            selected_network_name = self.network_table_widget.item(row, 0).text() # Assuming networks is a list of lists or tuples
            number_str = selected_network_name[1]  # Extract 'xxxx'
            network_index = int(number_str)-1 
            current_order.append(network_index)# Convert to an integer index

        # Prepare to save all ENewick strings
        all_enewicks = [self.enewicks[i] for i in current_order]  # Adjust as necessary
    
        # Open a file dialog to select the save location
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        file_name, _ = QFileDialog.getSaveFileName(self, "Save All ENewicks", "", "Text Files (*.txt);;All Files (*)", options=options)
    
        if file_name:
            # Check if the file_name has an extension
            if '.' not in file_name.split('/')[-1]:
                file_name += '.txt'
            with open(file_name, 'w') as file:
                for enewick in all_enewicks:
                    file.write(enewick + "\n")
            QMessageBox.information(self, "Success", "All ENewick strings saved successfully!")
    

    def plot_network(self, network=None):
        """Plot the selected network in the matplotlib canvas."""
        self.figure.clear()
        #self.figure.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05)
        ax = self.figure.add_subplot(111)
        
        ax.spines['bottom'].set_color('gray')
        ax.spines['top'].set_color('gray') 
        ax.spines['right'].set_color('gray')
        ax.spines['left'].set_color('gray')
        ax.get_xaxis().set_ticks([])
        ax.get_yaxis().set_ticks([])
        ax.format_coord = lambda x, y: ""
        self.figure.tight_layout()
        
        if network is None:
           # Show placeholder when no network is selected
           ax.set_facecolor('whitesmoke')
           ax.text(0.5, 0.5, 'No network selected', fontsize=14, color='dimgrey',
                   verticalalignment='center', horizontalalignment='center', transform=ax.transAxes)

        else:
            font_size = 12
    
            if isinstance(network, physquirrel.DirectedNetwork):
                #pos = nx.drawing.nx_agraph.graphviz_layout(network, prog='dot')
                pos = hierarchy_pos(network)#, vert_gap=0.5, width=5)
                pos = {node: (-y, x) for node, (x, y) in pos.items()}

                nx.draw_networkx_nodes(network, pos, ax=ax, nodelist=network.internal_nodes(), node_size=10, node_color='white',
                                   edgecolors='black', alpha=1)
                nx.draw_networkx_nodes(network, pos, ax=ax, nodelist=network.leaves, node_size=300, node_color='white',
                                   edgecolors='white', alpha=1)
                nx.draw_networkx_edges(network, pos, ax=ax, edgelist=network.non_reticulation_edges(), edge_color='black', width=1,
                                       node_size=200, alpha=0.85, arrows=True, arrowstyle='->', arrowsize=11)    
                nx.draw_networkx_edges(network, pos, ax=ax, edgelist=network.reticulation_edges(), edge_color='red', width=1,
                                       node_size=200, alpha=0.85, arrows=True, arrowstyle='->', arrowsize=11)    
                nx.draw_networkx_labels(network, pos, ax=ax, labels={leaf:leaf for leaf in network.leaves}, font_size=font_size)   

            elif isinstance(network, physquirrel.SemiDirectedNetwork):
                
                pos = nx.kamada_kawai_layout(network)

                nx.draw_networkx_nodes(network, pos, ax=ax, nodelist=network.internal_nodes(), node_size=10, node_color='white',
                                   edgecolors='black', alpha=1)
                nx.draw_networkx_nodes(network, pos, ax=ax, nodelist=network.leaves, node_size=300, node_color='white',
                                   edgecolors='white', alpha=1)
                nx.draw_networkx_edges(network, pos, ax=ax, edgelist=network.undirected_edges(), edge_color='black', width=1,
                                       node_size=200, alpha=0.85)
                nx.draw_networkx_edges(network, pos, ax=ax, edgelist=network.directed_edges, edge_color='red', width=1,
                                       node_size=200, alpha=0.85, arrows=True, arrowstyle='->', arrowsize=11)    
                nx.draw_networkx_labels(network, pos, ax=ax, labels={leaf:leaf for leaf in network.leaves}, font_size=font_size)   

        # Update the plot
        self.canvas.draw()
        self.canvas.updateGeometry()  # Ensure it updates when resized


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.initUI()
    
    def initUI(self, fileopen=False):
        self.setWindowIcon(QtGui.QIcon('squirrel_logo.ico'))
        self.setWindowTitle("Squirrel")
        self.setGeometry(100, 100, 1200, 600)

        # Main container widget
        self.central_widget = QWidget(self)
        self.setCentralWidget(self.central_widget)
        
        # Create a horizontal layout for the main layout
        main_layout = QHBoxLayout(self.central_widget)

        # Create a vertical layout for MSA and Algorithm frames
        left_layout = QVBoxLayout()
        
        # Create frame instances
        self.plot_frame = PlotFrame(self)
        self.alg_frame = AlgorithmFrame(self, plot_frame=self.plot_frame)
        self.msa_frame = MSAFrame(self, self.alg_frame, self)
        self.alg_frame.set_msa_frame(self.msa_frame)

        left_layout.addWidget(self.msa_frame)
        
        # Add a vertical separator (QFrame) between left and right sections
        left_layout.addSpacing(7)  # Add 10 pixels of space after the separator
        separator = QFrame()
        separator.setFrameShape(QFrame.HLine)  # Create a vertical line
        separator.setFrameShadow(QFrame.Sunken)  # Optional: make it appear sunken
        left_layout.addWidget(separator)
        left_layout.addSpacing(7)  # Add 10 pixels of space after the separator


        left_layout.addWidget(self.alg_frame)
        
        main_layout.addLayout(left_layout, 35)
        
        # Add a vertical separator (QFrame) between left and right sections
        main_layout.addSpacing(10)  # Add 10 pixels of space after the separator
        separator = QFrame()
        separator.setFrameShape(QFrame.VLine)  # Create a vertical line
        separator.setFrameShadow(QFrame.Sunken)  # Optional: make it appear sunken
        main_layout.addWidget(separator)
        main_layout.addSpacing(10)  # Add 10 pixels of space after the separator

        main_layout.addWidget(self.plot_frame, 100)
        
        self.central_widget.setLayout(main_layout)

        if fileopen:
            self.msa_frame.open_file(clear=False)

if __name__ == "__main__":
    app = QApplication(sys.argv)
    #app.setStyle("Fusion")  # You can use "Windows", "Macintosh", "Fusion", etc.
    
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())
    
