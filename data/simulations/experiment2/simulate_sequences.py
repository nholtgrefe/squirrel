import physquirrel as psq
import networkx as nx
import random, itertools, os, subprocess, re

# This file contains all the funtions to root semi-directed networks, assign branch lengths, and simulate alignments.
# The functions are listed first, and the script to simulate the alignments from a folder with semi-directed networks is
# at the end of the file.

#%% New class for Directed Networks with branch lengths

class WDirectedNetwork(psq.DirectedNetwork):
    """Subclass of psq.DirectedNetwork that allows for branch lengths."""
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.branch_lengths = {edge: 1.0 for edge in self.edges}
    
    def remove_edge(self, u, v):
        """Removes edge u,v from the network."""
        super().remove_edge(u, v)
        del self.branch_lengths[(u,v)]
        
    def remove_edges_from(self, edge_list):
        """Removes edges from the network."""
        for edge in edge_list:
            self.remove_edge(edge[0], edge[1])
        
    def add_edge(self, u, v, branch_length=1.0):
        """Adds edge u,v to the network."""
        super().add_edge(u, v)
        self.branch_lengths[(u,v)] = branch_length
        
    def add_edges_from(self, edge_list):
        """Add edges to the network."""
        for edge in edge_list:
            self.add_edge(edge[0], edge[1])
    
    def clear(self):
        """Clear the whole network."""
        super().clear()
        self.branch_lengths = dict()
        
    def non_reticulation_edges(self):
        """Returns a list of all non-reticulation edges of the network."""
        ret = []
        for v in self.reticulation_nodes():
            parents = self.predecessors(v)
            for p in parents:
                ret.append((p, v))
        
        return [edge for edge in self.edges if edge not in ret]
    
    def assign_new_branch_lengths(self, new_path_length=0.3):
        """
        Assigns every edge a new branch length as follows. Given an edge (u,v), let
        s be the average number of edges of all root-leaf paths containing (u,v).
        Then, the edge (u,v) is assigned a length of new_path_length * 1/s."""
        
        root = self.root_node()
        all_root_leaf_paths = []
        for leaf in self.leaves:
            all_paths = list(nx.all_simple_paths(self, source=root, target=leaf))
            for path in all_paths:
                all_root_leaf_paths.append(path)
            
        for edge in self.edges:
            # All paths containg the edge
            edge_paths = [path for path in all_root_leaf_paths if edge in zip(path, path[1:])]
            
            # Average length of these paths (in terms of edges)
            average_length = sum(len(path) for path in edge_paths) / len(edge_paths) - 1
            
            # Assign new length
            self.branch_lengths[edge] = new_path_length / average_length
            
    def displayed_trees(self):
        """Returns a list of all trees (as WDirectedNetworks) displayed by the network,
        accounting for mulitplicity. Also readjusts branch lengths when degree-2 nodes are contracted.
        Supresses roots with degree-1."""
        
        trees = []
        
        # Group directed edges by their reticulation nodes
        reticulation_edges = {}
        for edge in self.reticulation_edges():
            source, target = edge
            if target not in reticulation_edges:
                reticulation_edges[target] = []
            reticulation_edges[target].append(edge)
        
        # Get all groups of edges
        edge_groups = list(reticulation_edges.values())
        
        # Generate all combinations with one edge per reticulation node
        combinations = [list(comb) for comb in itertools.product(*edge_groups)]
        
        # Generate displayed trees
        for comb in combinations:
            # Create tree with degree-2 nodes
            tree = WDirectedNetwork(branch_lengths=self.branch_lengths)
            tree.add_edges_from(self.non_reticulation_edges())
            tree.add_edges_from(comb)
            tree.add_leaves_from(self.leaves)
            
            # Assigne branch lengths
            for edge in tree.edges:
                tree.branch_lengths[edge] = self.branch_lengths[edge]
        
            # Contract degree-2 nodes
            tree_degree2_nodes = [v for v in tree.nodes if tree.in_degree(v) == 1 and tree.out_degree(v) == 1]
            for v in tree_degree2_nodes:
                tree.supress_degree2_node(v)
                
            tree_root = tree.root_node()
            if tree.out_degree(tree_root) == 1:
                out_node = list(tree.successors(tree_root))[0]
                tree.remove_edge(tree_root, out_node)
                tree.remove_node(tree_root)
                
            trees.append(tree)    

        return trees
    
    def supress_degree2_node(self, v):
        """Supresses a degree-2 node and sums the two corresponding branch lengths."""
        if not self.in_degree(v) == 1 or not self.out_degree(v) == 1:
            raise ValueError
        u = list(self.predecessors(v))[0]
        w = list(self.successors(v))[0]
        new_branch_length = self.branch_lengths[(u,v)] + self.branch_lengths[(v,w)]
        self.remove_edges_from([(u,v), (v,w)])
        self.add_edge(u, w, new_branch_length)
        self.remove_node(v)
        
    def create_enewick(self, internal_nodes=True):
        """Returns the eNewick string of the network with branch lengths."""

        edge_list = list(self.edges)
        digraph = {}
        self._visited = {}
        inDegree = {}
    
        for arc in edge_list:
            if arc[0] not in digraph:
                digraph[arc[0]] = []
            self._visited[arc[0]] = 0
            self._visited[arc[1]] = 0
    
            # Update the indegree of the head of the arc
            if arc[1] not in inDegree:
                inDegree[arc[1]] = 1
            else:
                inDegree[arc[1]] += 1
            # Add indegree 0 to the tail of the arc if it has not been visited yet
            if arc[0] not in inDegree:
                inDegree[arc[0]] = 0
            # Add the arc to the digraph array
            digraph[arc[0]].append(arc[1])
         
        # Find the root
        root = "r"
        for node, degree in inDegree.items():
            if degree == 0:
                root = node
         
        self._reticulationVerticesFound = 0
        self._reticulationVertices = {}
        
        # Generate eNewick string starting at root
        converted = self._eNewick(root, digraph, inDegree, internal_nodes)
        
        return converted + ";"
            
    def _eNewick(self, source, digraph, inDegree, internalNodesDisplay):
        """Helper function to generate eNewick string with branch lengths."""

        eNewickString = ""
        # If source is a reticulation vertex, compute its number
        if inDegree[source] > 1:
            if source not in self._reticulationVertices:
                self._reticulationVerticesFound += 1
                reticulationNumber = self._reticulationVerticesFound
                self._reticulationVertices[source] = self._reticulationVerticesFound
            else:
                reticulationNumber = self._reticulationVertices[source]
        
        if self._visited[source] == 0:
            # If source was not visited yet, recursively visit its children
            self._visited[source] = 1
            if source in digraph:
                eNewickString = "("
                i = 0
                for child in digraph[source]:
                    if i > 0:
                        eNewickString += ","
                    # Recursively build the string for child nodes
                    eNewickString += self._eNewick(child, digraph, inDegree, internalNodesDisplay)
                    # Append branch length
                    length = self.branch_lengths.get((source, child), 0.0)
                    eNewickString += f":{length}"
                    i += 1
                eNewickString += ")"
        if internalNodesDisplay == 1 or not(source in digraph):
            eNewickString += str(source)
        # If source is a reticulation vertex, label it with its number
        if inDegree[source] > 1:
            eNewickString += f"#H{reticulationNumber}"
        return eNewickString  
    
    def save_to_file(self, filename, overwrite=False):
        """Saves the network to the given file as an enewick string."""
        if overwrite == False and os.path.exists(filename):
            raise ValueError("File already exists.")

        with open(filename, 'w') as file:
            enewick = self.create_enewick()
            file.write(enewick)
            
    def save_displayed_trees_to_file(self, filename, overwrite=False, sequence_length=None):
        """Saves the displayed trees of the network as enewick strings to a single file.
        Takes as optional argument an integer 'sequence length', which is then divided equally over the 
        displayed trees for input to seq-gen. In particular every line gets [len] in front."""
        if overwrite == False and os.path.exists(filename):
            raise ValueError("File already exists.")

        with open(filename, 'w') as file:
            trees = self.displayed_trees()
            nr_trees = len(trees)                
            
            if sequence_length is not None:
                base_sub_length = sequence_length // nr_trees
                remainder = sequence_length % nr_trees
                        
            for i, tree in enumerate(trees):
                if sequence_length is not None:
                    length = base_sub_length
                    if i+1 <= remainder:
                        length += 1
                    if length == 0:
                        continue
                    file.write(f"[{length}]")

                enewick = tree.create_enewick(internal_nodes=False)
                file.write(enewick + "\n")
                
    def load_from_enewick(self, enewick_string):
        """Clears the network and loads the network from the given enewick_string (including branch lengths).
        Requires 'phylox'. Download with 'pip install phylox' first."""
        
        import phylox
        from phylox.newick_parser import dinetwork_to_extended_newick, extended_newick_to_dinetwork

        self.clear()
        net = extended_newick_to_dinetwork(enewick_string, internal_labels=True)
        node_mapping = {u: net.nodes[u].get("label") for u in net.nodes}
        leaf_list = [node_mapping[leaf] for leaf in net.leaves]
        
        for(net_u, net_v) in net.edges:
            u = node_mapping[net_u]
            v = node_mapping[net_v]
            length = net[net_u][net_v]['length']
            self.add_edge(u, v, branch_length=length)
        
        self.add_leaves_from(leaf_list)
        
    def load_from_file(self, filename):
        """Loads networks from a file with enewick string."""
        
        if not os.path.dirname(filename):
            filename = os.path.join(os.getcwd(), filename)
    
        with open(filename, 'r') as file:            
            lines = file.readlines()[0]
            self.load_from_enewick(lines)
            
            
#%% Functions to turn given SemiDirectedNetworks into DirectedNetworks with branch lengths

def is_valid_root_location(sdnetwork, edge):
    """Takes as input an edge of a SemiDirectedNetwork and returns whether the network
    can be rooted on that edge."""
    
    if edge not in sdnetwork.edges:
        raise ValueError
        
    u, v = edge
    if v in sdnetwork.leaves:
        start = u
    elif (u,v) in sdnetwork.directed_edges:
        start = u
    else:
        start = v
    
    # Check if valid root location
    for leaf in sdnetwork.leaves:
        reachable = False
        for path in nx.all_simple_paths(sdnetwork, source=start, target=leaf):
            good_path = True
            for i in range(len(path) - 1):
                # Check if there is a directed edge from sequence[i] to sequence[i+1]
                if (path[i+1], path[i]) in sdnetwork.directed_edges:
                    good_path = False
                    break
            if good_path == True:
                reachable = True
                break
        if reachable == False:
            return False
    
    return True
    

def sdnetwork_to_wdnetwork(sdnetwork):
    """Chooses an arbitrary root location for a SemiDirectedNetwork and turns it into
    a WDirectedNetwork with unit branch lengths."""
    
    root_locations = random.sample(list(sdnetwork.edges), len(sdnetwork.edges))
    
    for root_location in root_locations:
        valid = is_valid_root_location(sdnetwork, root_location)
        if valid == False:
            continue
        else:
            break
        
    u, v = root_location
    if v in sdnetwork.leaves:
        start = u
    elif (u,v) in sdnetwork.directed_edges:
        start = u
    else:
        start = v
    
    good_paths = []
    for target in sdnetwork.leaves:
        for path in nx.all_simple_paths(sdnetwork, source=start, target=target):
            good_path = True
            for i in range(len(path) - 1):
                # Check if there is a directed edge from sequence[i] to sequence[i+1]
                if (path[i+1], path[i]) in sdnetwork.directed_edges:
                    good_path = False
                    break
            if good_path == True:
                good_paths.append(path)
    
    directed_edges = set()
    
    for path in good_paths:
        for i in range(len(path) - 1):
            directed_edges = directed_edges | {(path[i], path[i+1])}
            
    directed_edges = directed_edges - {(u,v), (v,u)}
    
    new_root = psq.utils.id_generator()
    directed_edges = directed_edges | {(new_root, u),(new_root, v)}
    
    D = WDirectedNetwork()
    D.add_edges_from(directed_edges)
    D.add_leaves_from(sdnetwork.leaves)
    
    return D

def create_weighted_networks(input_folder, output_folder):
    """Chooses an arbitrary root for each SemiDirectedNetwork in the input folder. Then, turns this into
    a WDirectedNetwork with branch lengths and saves it."""
    
    files = os.listdir(input_folder)
    files.sort()
    
    for filename in files:
        SDN = psq.SemiDirectedNetwork()
        SDN.load_from_file(input_folder+filename)
        
        WDN = sdnetwork_to_wdnetwork(SDN)
        WDN.assign_new_branch_lengths()  
                
        new_filename = "W" + filename[1:]
        WDN.save_to_file(output_folder + new_filename)



#%% Function to save the displayed trees of the created WDirectedNetworks in seq-gen format

def save_displayed_trees(input_folder, output_folder, network_alignment_lengths):
    """Takes as input a folder with WDNetworks. For every network in this folder and every
    sequence length in network_alignment_lengths, this function saves a file with all displayed
    trees of the network preceded by a partial sequence length for that tree (used as input for seq-gen)."""
    
    files = os.listdir(input_folder)
    files.sort()
    
    for filename in files:
        WDN = WDirectedNetwork()
        WDN.load_from_file(input_folder+filename)
        for length in network_alignment_lengths:
            len_string = str(length).zfill(12)
            new_file = "DT" + filename[3:-4] + f"_LEN{len_string}.txt"
            WDN.save_displayed_trees_to_file(output_folder + new_file, sequence_length=length)



#%% Functions to simulate an alignment for each displayed tree file

def run_seq_gen(input_file, output_file, length, nr_trees, seqgen_path, model="HKY", tstv=4.0):
    """
    Runs Seq-Gen using a file of Newick strings as input and saves the output in output_file.
    Args:
        model (str): Substitution model (default is "HKY").
        length (int): Sequence length
        seqgen_path (str): path to seqgen binary
        nr_trees (int): number of trees in the input_file to simulate along
        tstv (float): TS/TV ratio (default = 4.0)
    """

    # Prepare the Seq-Gen command
    command = [
        seqgen_path,
        "-m" + model,    # Substitution model
        "-l" + str(length),  # Sequence length
        "-t" + str(tstv),    # TS/TV parameter
        "-p" + str(nr_trees),    # Number of displayed trees
        "-on" # Nexus output format
    ]

    # Open input and output files
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        # Run Seq-Gen with input file as stdin and capture output
        result = subprocess.run(command, stdin=infile, stdout=outfile, text=True, check=True)
    
    return result.stdout

def simulate_alignments(input_folder, output_folder, seqgen_path):
    """Simulates a sequence alignment (in .fasta format) along all the displayed tree files in the 'input_folder'."""
    
    files = os.listdir(input_folder)
    files.sort()
    
    for filename in files:
        match = re.match(r"DT(\d+)_L(\d+)_R(\d+)_LEN(\d+)", filename)
        id_string, leaf_string, ret_string, length_string = match.groups()
        
        new_filename = 'SEQ' + filename[2:-4] + '.nexus'
        length = int(length_string)
        nr_trees = min(2**int(ret_string), length)
        
        run_seq_gen(input_folder+filename, output_folder+new_filename, length=length, nr_trees=nr_trees, seqgen_path=seqgen_path, model="HKY", tstv=4.0)

#%% Code to create simulated sequence alignments

# Folder containing all randomly generated SemiDirectedNetworks (from Experiment 1)
input_folder_sdnetworks = "/path/to/folder/with/generated_sdnetworks/"

# Folder to save all newly created DirectedNetworks with branch lengths
output_folder_wdnetworks = "/path/to/folder/with/generated_dnetworks/"

# Turn all SD-networks into WD-networks
create_weighted_networks(input_folder_sdnetworks, output_folder_wdnetworks)

# Save the displayed trees of each WD-network in seq-gen format (one file per sequence length)
output_folder_displayed_trees = "/path/to/temporary/folder/with/seqgen/input/displayed_trees/"
network_alignment_lengths = [1000, 10000, 100000, 1000000]
save_displayed_trees(output_folder_wdnetworks, output_folder_displayed_trees, network_alignment_lengths)

# Simulate an alignment for every displayed trees file, assumes the SeqGen program is installed.
output_folder_alignments = "/path/to/folder/with/simulated_alignments/"
seq_gen_binary = "/path/to/seqgen/binary/Seq-Gen-master/source/seq-gen"
simulate_alignments(output_folder_displayed_trees, output_folder_alignments, seq_gen_binary)
