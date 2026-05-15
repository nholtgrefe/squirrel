import networkx as nx
import os

from ._config import validate

############################################################

class MixedGraph(nx.Graph):
    """
    Class for partially directed graphs (i.e. mixed graphs). Does not allow for 
    multiple edges between two nodes, even if they have a different direction.
    'directed_edges' is optional input for the directed edges of the graph.
        self.directed_edges: set of tuples of all directed_edges.
    """
    def __init__(self, directed_edges=None, *args, **kwargs):
        super().__init__(*args, **kwargs)        
        self.directed_edges = set(directed_edges) if directed_edges else set()

        if validate():
            if not self.directed_edges.issubset(self.edges):
                raise ValueError("Directed edges not in the edge set.")
    
    def _clear_cache(self):
       pass
       
    @classmethod
    def from_graph(cls, graph, directed_edges=None):
        """Turns an nx.Graph object into a MixedGraph object."""
        return cls(incoming_graph_data=graph.edges, directed_edges=directed_edges)

    def remove_node(self, v):
        """Removes node v from the graph."""
        edges_to_remove = []
        for (s,t) in self.directed_edges:
            if v == s or v == t:
                edges_to_remove.append((s,t))
        self.remove_directed_edges_from(edges_to_remove)
        super().remove_node(v)
        self._clear_cache()
    
    def remove_nodes_from(self, nodes):
        """Removes all nodes in 'nodes' from the graph."""
        for v in nodes:
            self.remove_node(v)
    
    def remove_edge(self, u,v):
        """Removes edge (u,v) from the graph."""
        super().remove_edge(u,v)
        if (u,v) in self.directed_edges:
            self.directed_edges.remove((u,v))
        elif (v,u) in self.directed_edges:
            self.directed_edges.remove((v,u))
        self._clear_cache()
    
    def remove_edges_from(self, edges):
        """Removes all edges in 'edges' from the graph."""
        for (u,v) in edges:
            self.remove_edge(u,v)
        
    def remove_directed_edge(self, u,v):
        """Removes directed edge (u,v) from the graph."""
        if validate():
            if (u,v) not in self.directed_edges:
                raise ValueError("Edge does not exist or has different direction.")
        self.remove_edge(u,v)
    
    def remove_directed_edges_from(self, edges):
        """Removes all directed edges in 'edges' from the graph."""
        for (u,v) in edges:
            self.remove_directed_edge(u,v)
            
    def add_directed_edge(self, u, v):
        """Adds edge directed edge (u,v) to the graph."""
        self.add_edge(u, v)
        if (v,u) in self.directed_edges:
            self.directed_edges.remove((v,u))
        self.directed_edges.add((u,v))
        
    def add_directed_edges_from(self, edges):
        """Adds all directed edges in 'edges' to the graph."""
        for (u,v) in edges:
            self.add_directed_edge(u,v)
     
    def undirected_edges(self):
        """Returns all undirected edges of the graph."""
        return [(u,v) for (u,v) in self.edges() if (u,v) not in self.directed_edges and (v,u) not in self.directed_edges]
    
    def indegree(self, v):
        """Returns the indegree of the vertex v."""
        return len([u for u in self.neighbors(v) if (u,v) in self.directed_edges])

    def outdegree(self, v):
        """Returns the outdegree of the vertex v."""
        return len([u for u in self.neighbors(v) if (v,u) in self.directed_edges])
    
    def is_cutedge(self, u, v):
        """Checks if (u,v) is a cut-edge. Raises an error if it is not an edge in the graph."""
        if validate():
            if (u,v) not in self.edges:
                raise ValueError("Not an edge in the graph.")
        modified_graph = self.copy()
        modified_graph.remove_edge(u,v)
        return nx.number_connected_components(modified_graph) != nx.number_connected_components(self)
            
    def is_cutvertex(self, v):
        """Checks if v is a cut-vertex. Raises an error if it is not a vertex in the graph."""
        if validate():
            if v not in self.nodes:
                raise ValueError("Not a vertex in the graph.")
        modified_graph = self.copy()
        modified_graph.remove_node(v)
        return nx.number_connected_components(modified_graph) != nx.number_connected_components(self)
            
    def identify_two_nodes(self, u, v):
        """Identifies two nodes u and v by the node u."""
        nx.identified_nodes(self, u, v, self_loops=False, copy=False)
        if (u,v) in self.directed_edges:
            self.directed_edges.remove((u,v))
        elif (v,u) in self.directed_edges:
            self.directed_edges.remove((v,u))
        
    def identify_node_set(self, nodes):
        """Identifies all nodes by the first node in the list/set."""
        nodes = list(nodes)
        for i in range(1, len(nodes)):
            self.identify_two_nodes(nodes[0], nodes[i])
    
    def clear(self):
        """Clear the whole mixed graph."""
        super().clear()
        self.directed_edges = set()
        self._clear_cache()

    def save_to_file(self, filename, overwrite=False):
        """Saves the mixed graph to the given file as a list of undirected
        edges and a list of directed edges. If overwrite = False, an error
        is raised when the file already exists."""
        if overwrite == False and os.path.exists(filename):
            raise ValueError("File already exists.")

        with open(filename, 'w') as file:
            file.write("#Undirected\n")
            for u, v in self.undirected_edges():
                file.write(f"{u} {v}\n")
            
            file.write("#Directed\n")
            for u, v in self.directed_edges:
                file.write(f"{u} {v}\n")
    
    def load_from_file(self, filename):
        """Loads the mixed graph from the given file which should contain a list 
        of undirected edges and a list of directed edges. If the filename does not
        specify a folder, the working directory is assumed. Raises an error
        if the file is in an incorrect format."""
        
        if not os.path.dirname(filename):
            filename = os.path.join(os.getcwd(), filename)
    
        with open(filename, 'r') as file:
            self.clear()
            
            lines = file.readlines()
            section = None
                        
            for line in lines:
                line = line.strip()
                if line.startswith("#Undirected"):
                    section = "undirected"
                    continue
                elif line.startswith("#Directed"):
                    section = "directed"
                    continue
                
                if section == "undirected":
                    line_lst = line.split()
                    if len(line_lst) != 2:
                        raise ValueError("Incorrect format.")
                    u, v = line_lst
                    self.add_edge(u, v)
                elif section == "directed":
                    line_lst = line.split("")
                    if len(line_lst) != 2:
                        raise ValueError("Incorrect format.")
                    u, v = line_lst
                    self.add_directed_edge(u, v)

