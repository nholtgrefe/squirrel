import networkx as nx
import matplotlib.pyplot as plt

from .utils import validate

############################################################

class DirectedNetwork(nx.DiGraph):
    """
    Class for directed networks, subclass of nx.DiGraph. The argument 
    'leaves' is optional input for nodes that need to be assigned as leaves.
        self.leaves: set of leaves of the network.
    """
    def __init__(self, leaves=None, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.leaves = set(leaves) if leaves else set()
            
        if validate():
            if not self.leaves.issubset(self.nodes):
                raise ValueError("Leaves not in node set.")
        
        # Used to generate eNewick strings
        self._reticulationVerticesFound = None
        self._visited = None
        self._reticulationVertices = None
        
    def __repr__(self):
        if len(self.leaves) < 15:
            return f"Directed network on {self.leaves}"
        else:
            return "Directed network on >15 leaves"
    
    def remove_node(self, v):
        """Removes node v from the network."""
        super().remove_node(v)
        if v in self.leaves:
            self.leaves.remove(v)
            
    def remove_nodes_from(self, nodes):
        """Removes all nodes in 'nodes' from the network."""
        for v in nodes:
            self.remove_node(v)
        
    def remove_leaf(self, leaf):
        """Removes the leaf from the network."""
        if leaf not in self.leaves:
            raise ValueError("Vertex is not a leaf in the network.")
        self.remove_node(leaf)
    
    def remove_leaves_from(self, leaves):
        """Removes all leaves in 'leaves' from the network."""
        for leaf in leaves:
            self.remove_leaf(leaf)

    def add_leaf(self, leaf):
        """Adds a leaf to the network."""
        self.add_node(leaf)
        self.leaves.add(leaf)

    def add_leaves_from(self, leaves):
        """Adds all leaves in 'leaves' to the network."""
        for leaf in leaves:
            self.add_leaf(leaf)
    
    def copy(self):
        """Returns a copy of the network."""
        N = DirectedNetwork()
        N.add_edges_from(self.edges)
        N.add_nodes_from(self.nodes)
        N.add_leaves_from(self.leaves)
        return N
    
    def clear(self):
        """Clear the whole network."""
        super().clear()
        self.leaves = set()
    
    def reticulation_nodes(self):
        """Returns a list of all reticulation nodes of the network."""
        return [v for v in self.nodes() if self.in_degree(v) == 2 and self.degree(v) == 3]
    
    def reticulation_edges(self):
        """Returns a list of all reticulation edges of the network."""
        res = []
        for v in self.reticulation_nodes():
            parents = self.predecessors(v)
            for p in parents:
                res.append((p, v))
        return res
    
    def non_reticulation_edges(self):
        """Returns a list of all non_reticulation edges of the network."""
        ret_edges = self.reticulation_edges()
        return [(u,v) for (u,v) in self.edges if (u,v) not in ret_edges]
    
    def internal_nodes(self):
        """Returns a list of all non-leaf/internal nodes of the network."""
        return [v for v in self.nodes() if v not in self.leaves]
    
    def root_node(self):
        """Returns the root node of the network."""
        return [v for v in self.nodes if self.in_degree(v)==0][0]
    
    def is_tree(self):
        """Checks if the network is a tree."""
        return len(self.reticulation_nodes()) == 0

    def visualize(self, layout='custom', title=None, leaflabels=True, internal_labels=False, font_size=12):
        """Visualization function with several layout-options: ['custom', 'dot', 'kamada', 'neato', 'twopi', 'circo'].
        If graphviz is not installed, use 'kamada' or 'custom'. Optional title can be given to the plot.
        Label printing can be turned off."""
        if layout == 'kamada':
            pos = nx.kamada_kawai_layout(self)
        elif layout == 'custom':
            pos = _hierarchy_pos(self)#, vert_gap=0.5, width=5)
            pos = {node: (-y, x) for node, (x, y) in pos.items()}
        else:
            pos = nx.drawing.nx_agraph.graphviz_layout(self, prog=layout)
        
        fig, ax = plt.subplots(figsize=(9, 7), dpi=200)
        
        if internal_labels == True:
            int_node_size = 300
        else:
            int_node_size = 50
            
        nx.draw_networkx_nodes(self, pos, nodelist=self.internal_nodes(), node_size=int_node_size, node_color='white',
                           edgecolors='black', alpha=1) #, label="Internal")
        nx.draw_networkx_nodes(self, pos, nodelist=self.leaves, node_size=300, node_color='white',
                           edgecolors='white', alpha=1) #, label="Leaf")
        nx.draw_networkx_edges(self, pos, edgelist=self.edges, edge_color='black', width=1,
                               node_size=200, alpha=0.85, arrows=True, arrowstyle='->', arrowsize=11)    
        if leaflabels == True:
            nx.draw_networkx_labels(self, pos, labels={leaf:leaf for leaf in self.leaves}, font_size=font_size)   
        if internal_labels == True:
            nx.draw_networkx_labels(self, pos, labels={node:node for node in self.internal_nodes()}, font_size=font_size)   


        if title == None: title = "Directed Network"
        plt.title(title)
        plt.show()
        return ax

    def load_from_enewick(self, enewick_string):
        """Clears the network and loads the network from the given enewick_string.
        Requires 'phylox'. Run 'pip install phylox' first."""
        
        import phylox
        from phylox.newick_parser import dinetwork_to_extended_newick, extended_newick_to_dinetwork

        self.clear()
        net = extended_newick_to_dinetwork(enewick_string)
        mapping = {leaf: net.nodes[leaf].get("label") for leaf in net.leaves}

        self.add_edges_from(net.edges)
        self = nx.relabel_nodes(self, mapping, copy=False)
        self.add_leaves_from(list(mapping.values()))

    def create_enewick(self):
        """Returns the enewick string of the network.
        Code adapted from Phillipe Gambette:
            https://github.com/PhilippeGambette/networkDraw/tree/main
        """       

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
        for node,degree in inDegree.items():
           if degree == 0:
              root = node
         
        self._reticulationVerticesFound = 0
        self._reticulationVertices = {}
        
        converted = self._eNewick(root, digraph, inDegree, 0)
        
        return converted + ";"
            
    def _eNewick(self, source, digraph, inDegree, internalNodesDisplay):
        """Helper function to generate eNewick string.
        Code adapted from Phillipe Gambette:
            https://github.com/PhilippeGambette/networkDraw/tree/main
        """

        eNewickString = ""
        # if source is a reticulation vertex, compute its number
        if inDegree[source] > 1:
            if source not in self._reticulationVertices:
                self._reticulationVerticesFound += 1
                reticulationNumber = self._reticulationVerticesFound
                self._reticulationVertices[source] = self._reticulationVerticesFound
            else:
                reticulationNumber = self._reticulationVertices[source]
        
        if self._visited[source] == 0:
            # if source was not visited yet, recursively visit its children
            self._visited[source] = 1
            if source in digraph:
                eNewickString = "("
                i = 0
                for child in digraph[source]:
                    if i > 0:
                        eNewickString += ","
                    eNewickString += self._eNewick(child, digraph, inDegree, internalNodesDisplay)
                    i += 1
                eNewickString += ")"
        if internalNodesDisplay == 1 or not(source in digraph):
            eNewickString += source
        # if source is a reticulation vertex, label it with its number
        if inDegree[source] > 1:
            eNewickString += "#H" + str(reticulationNumber)
        return eNewickString
    
    def load_from_edgelist(self, el_string):
        """Clears the network and loads the network from the given edge_list string."""
        pass
    
    def create_edgelist(self):
        """Returns the edgelist string of the network."""
        
    def save_to_file(self, filename, filetype='enewick', overwrite=False):
        """Saves the directed network in 'enewick' or 'edgelist' type."""
        pass
    
    def load_from_file(self, filename):
        """Loads the network from a file in enewick or edgelist format."""
        pass

############################################################


def _longest_distance_to_root(G, root):
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

def _hierarchy_pos(G, root=None, width=5., vert_gap=0.5, xcenter=0.5):
    if root is None:
        root = next(iter(nx.topological_sort(G)))  # Start at a topologically sorted root if none is provided

    # Compute longest distances to the root
    distances = _longest_distance_to_root(G, root)

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
