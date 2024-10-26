import networkx as nx
import matplotlib.pyplot as plt
import os, itertools, random
from functools import lru_cache

from .utils import MixedGraph, Partition, CircularOrdering, CircularSetOrdering
from .utils import id_generator, validate

from .dnetwork import DirectedNetwork
from .splits import Split, SplitSystem, QuartetSplitSet
from .quarnet import QuartetTree, FourCycle, SingleTriangle, DoubleTriangle, CycleQuarnet, SplitQuarnet
from .quarnetset import DenseQuarnetSet
from .trinet import Triangle, ThreeStar
from .trinetset import DenseTrinetSet

############################################################

class SemiDirectedNetwork(MixedGraph):
    """
    Class for semi-directed networks, subclass of MixedGraph. The argument 
    'leaves' is optional input for nodes that need to be assigned as leaves.
        self.leaves: set of leaves of the network.
    """
    def __init__(self, leaves=None, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.leaves = set(leaves) if leaves else set()
            
        if validate():
            if not self.leaves.issubset(self.nodes):
                raise ValueError("Leaves not in node set.")
    
    def _clear_cache(self):
        """Clear the cache for certain methods."""
        super()._clear_cache()
        self.level.cache_clear()
        self.blobs.cache_clear()
        self.reticulation_nodes.cache_clear()
        self.internal_nodes.cache_clear()
        self.splits.cache_clear()
        self.trinets.cache_clear()
        self.quarnets.cache_clear()
       
    @classmethod
    def from_graph(cls, graph, leaves=None, directed_edges=None):
        """Turns an nx.Graph object into a SemiDirectedNetwork object."""
        return cls(incoming_graph_data=graph.edges, directed_edges=directed_edges, leaves=leaves)

    @classmethod
    def from_mixedgraph(cls, mixedgraph, leaves=None):
        """Turns a MixedGraph object into a SemiDirectedNetwork object."""
        return cls(incoming_graph_data=mixedgraph.edges, directed_edges=mixedgraph.directed_edges, leaves=leaves)

    def __repr__(self):
        if len(self.leaves) < 15:
            return f"Semi-directed network on {self.leaves}"
        else:
            return "Semi-directed network on >15 leaves"
    
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
        N = SemiDirectedNetwork()
        N.add_edges_from(self.edges)
        N.add_nodes_from(self.nodes)
        N.add_leaves_from(self.leaves)
        N.add_directed_edges_from(self.directed_edges)
        return N
    
    def clear(self):
        """Clear the whole network."""
        super().clear()
        self.leaves = set()
    
    @lru_cache(maxsize=1)
    def reticulation_nodes(self):
        """Returns a list of all reticulation nodes of the network."""
        return [v for v in self.nodes() if self.indegree(v) == 2 and self.degree(v) == 3]
    
    @lru_cache(maxsize=1)
    def internal_nodes(self):
        """Returns a list of all non-leaf/internal nodes of the network."""
        return [v for v in self.nodes() if v not in self.leaves]
    
    def is_tree(self):
        """Checks if the network is a tree."""
        return len(self.reticulation_nodes()) == 0
    
    def is_rootable(self):
        """Checks if the network has a possible root-location. Only implemented for level-1 networks."""
        # Check for valid root location, by checking whether a path between two
        # reticulations contains at least one reticulatione edge of one of the reticulations
        for r1, r2 in itertools.combinations(self.reticulation_nodes(), 2):
            ret_parents = [u for u in self.nodes if (u,r1) in self.directed_edges or (u,r2) in self.directed_edges]

            path = nx.shortest_path(self, source=r1, target=r2,method='dijkstra')
            if len(set(path) & set(ret_parents)) == 0:
                return False

        return True
    
    def is_rootable_at(self, leaf):
        """Checks if the network can be rooted at the edge connecting 'leaf' to the rest of the network."""
        if validate():
            if not leaf in self.leaves:
                raise ValueError("Leaf not a leaf of the network")
        
        for target in self.leaves:
            if target == leaf: continue
            reachable = False
            for path in nx.all_simple_paths(self, source=leaf, target=target):
                good_path = True
                for i in range(len(path) - 1):
                    # Check if there is a directed edge from sequence[i] to sequence[i+1]
                    if (path[i+1], path[i]) in self.directed_edges:
                        good_path = False
                        break
                if good_path == True:
                    reachable = True
                    break
            if reachable == False:
                return False
        return True
                
    def is_semidirected(self):
        """Checks if the network is indeed a true semi-directed network (with a
        possible root-location). Only implemented for level-1 networks."""
        if validate():
            if not self.level() == 1:
                raise ValueError("Only implemented for level-1 networks.")
                
        # Check connectivity
        if not nx.is_connected(self):
            return False
        
        # Check degrees
        for v in self.nodes:
            if v in self.leaves:
                if self.degree(v) != 1:
                    return False
            else:
                if self.degree(v) != 3:
                    return False
                if self.indegree(v) not in [0,2]:
                    return False
                
        # Check if reticulation edges are indeed not cut-edges
        for (u,v) in self.directed_edges:
            if self.is_cutedge(u, v):
                return False

        return self.is_rootable()
    
    def to_directed_network(self, leaf):
        """Generates a DirectedNetwork from the SemiDirectedNetwork, by choosing
        some valid leaf as the outgroup. Raises an error if the network is not semidirected."""
        if validate():
            if not self.is_semidirected():
                raise ValueError("Not semi directed.")
            if not self.is_rootable_at(leaf):
                raise ValueError("Not rootable at given leaf")
        
        good_paths = []
        for target in self.leaves:
            if target == leaf: continue
            for path in nx.all_simple_paths(self, source=leaf, target=target):
                good_path = True
                for i in range(len(path) - 1):
                    # Check if there is a directed edge from sequence[i] to sequence[i+1]
                    if (path[i+1], path[i]) in self.directed_edges:
                        good_path = False
                        break
                if good_path == True:
                    good_paths.append(path)
        
        directed_edges = set()
        
        for path in good_paths:
            for i in range(len(path) - 1):
                directed_edges = directed_edges | {(path[i], path[i+1])}
                
        x = list(self.neighbors(leaf))[0]
        directed_edges = directed_edges - {(leaf, x)}
        
        new_root = id_generator()
        directed_edges = directed_edges | {(new_root, x),(new_root, leaf)}
        
        D = DirectedNetwork()
        D.add_edges_from(directed_edges)
        D.add_leaves_from(self.leaves)

        return D
        
    
    @lru_cache(maxsize=1)
    def blobs(self, include_trivial=True, include_leaves=True):
        """Returns a list of node-sets, each of which make up a blob of the network.
        include_trivial/include_leaves indicate whether the trivial/leaf blobs should
        be included."""
        if include_trivial == False: include_leaves = False
        blobs = []
        visited = set()
        
        for component in list(nx.biconnected_components(self)):
            if len(component) > 2:
                blobs.append(component)
                visited = visited | component
                
        for v in self.nodes:
            if v not in visited:
                if v not in self.leaves and include_trivial == True:
                    blobs.append({v})
                elif v in self.leaves and include_leaves == True:
                    blobs.append({v})

        return blobs
    
    @lru_cache(maxsize=1)
    def level(self):
        """Returns the level of the network."""
        max_lev = 0
        reticulations = self.reticulation_nodes()
        for blob in self.blobs(include_trivial=False):
            nr_blob_rets = len([v for v in blob if v in reticulations])
            max_lev = max(max_lev, nr_blob_rets)
        return max_lev

    def cutedges(self, include_trivial=True):
        """Returns a list containing all cutedges of the network. The trivial
        cutedges can optionally be excluded."""
        cutedges = []
        for component in nx.biconnected_components(self):
            if len(component) != 2: continue
            u, v = tuple(component)
            if u not in self.leaves and v not in self.leaves:
                cutedges.append((u,v))
            else:
                if include_trivial == True: cutedges.append((u,v))
        return cutedges
    
    def similarity(self, network, measure="QC", reference=None, **kwargs):
        """Returns the similarity (in [0,1]) of the network with another network.
        Different possible measures:
            QC: quarnet consistency
            QCprime: quarnet consistency w.r.t. to reference DenseQuarnetSet
            QS: 1 - normalized quarnet symmetric difference
            TC: trinet consistency (C-measure)
            TS: 1 - normalized trinet symmetric difference (S-measure)
        Raises an error if reference is incompatible with given measure
        kwargs are passed to the corresponding function.
        See: 'A Practical Algorithm for Reconstructing Level-1 Phylogenetic Networks' by Huber et al."""
        
        if validate():
            if not measure in ['QC', 'QCprime', 'QS', 'TC', 'TS']:
                raise ValueError
            if measure == "QCprime" and not isinstance(reference, DenseQuarnetSet):
                raise ValueError

        if measure == "QC":
            return self.quarnets().consistency(network.quarnets(), **kwargs)
        elif measure == "QCprime":
            return self.quarnets().consistency_prime(network.quarnets(), reference, **kwargs)
        elif measure == "QS":
            return 1 - self.quarnets().distance(network.quarnets(), **kwargs)
        elif measure == "TC":
            return self.trinets().consistency(network.trinets())
        elif measure == "TS":
            return 1 - self.trinets().distance(network.trinets())

    def split_from_cutedge(self, u, v):
        """Returns the split induced by the cut-edge (u,v), with set1 on the 
        side of u and set2 on the side of v. Raises an error if (u,v) is not a
        cut-edge."""
        if validate():
            if not self.is_cutedge(u,v):
                raise ValueError("Not a cut-edge")

        modified_graph = self.copy()
        modified_graph.remove_edge(u, v)
        components = list(nx.connected_components(modified_graph))
        assert len(components) == 2
        
        if u in components[0]:
            set1 = set(components[0]) & self.leaves
            set2 = set(components[1]) & self.leaves
        elif u in components[1]:
            set1 = set(components[1]) & self.leaves
            set2 = set(components[0]) & self.leaves
            
        return Split(set1, set2)
    
    def cutedge_from_split(self, split):
        """Returns the cut-edge (u,v) that induces the given split, with u on the side
        of set1 and v on the side of set2. Raises an error if the split is not
        a split in the network."""
        if validate():
            if not split in self.splits():
                raise ValueError("Not a split")
        
        for (u,v) in self.cutedges():
            tested_split = self.split_from_cutedge(u, v)
            if tested_split == split:
                if tested_split.set1 == split.set1:
                    return (u, v)
                else:
                    return (v, u)

    def partition_from_cutvertex(self, v, return_cutedges=False):
        """Returns the partition induced by a cut-vertex v. If return_cutedges=True,
        the method also returns a dictionary where each key is part of the partition,
        and the value is the corresponding cut-edge. Raises an error 
        if v is not a cut-vertex or if it is a leaf."""
        if validate():
            if v in self.leaves:
                raise ValueError("Not an internal cut-vertex")
            if not self.is_cutvertex(v):
                raise ValueError("Not a cut-vertex")
        
        d = dict()
        for u in self.nodes:
            if (u,v) in self.edges:
                s = self.split_from_cutedge(u,v)
                d[frozenset(s.set1)] = (u,v)
        
        if return_cutedges == False:
            return Partition(list(d.keys()))
        else:
            return Partition(list(d.keys())), d
    
    def quartetsplits_from_cutedge(self, u, v):
        """Needs (u,v) to be a non-trivial cutedge and u,v to be cut-vertices. If A1|...|Ak and B1|...|Bl are
        the partitions corresponding to u and v (with split A|B), the method
        returns a QuartetSplitSet with all split a_i a_j | b_l b_k from different sets.
        Raises an error if (u,v) is not a cutedge or u and v are not cutvertices."""
        if validate():
            if not self.is_cutedge(u, v):
                raise ValueError
            if not self.is_cutvertex(u) or not self.is_cutvertex(v):
                raise ValueError
                
        s = self.split_from_cutedge(u, v)
        if validate():
            if s.is_trivial():
                raise ValueError
                
        p1 = self.partition_from_cutvertex(u)
        p2 = self.partition_from_cutvertex(v)
        
        p1_new = Partition([p for p in p1 if p != s.set2])
        p2_new = Partition([p for p in p2 if p != s.set1])
        
        splits = []
        for pp1 in p1_new.subpartitions(2):
            for pp2 in p2_new.subpartitions(2):
                for rp1 in pp1.representative_partitions():
                    for rp2 in pp2.representative_partitions():
                        splits.append(Split(rp1.elements, rp2.elements))

        return QuartetSplitSet(set(splits))
        
    def order_from_cycle(self, cycle):
        """Takes as input a circular ordering of vertices that form a cycle.
        Returns the circular setordering induced by the cycle.
        Also returns the part of the partition that is below the reticulation of the cycle.        
        Raises an error if cycle is not a cycle or the network is not level-1."""
        if validate():
            if self.level() < 1:
                raise ValueError("Not level-1.")
            if cycle not in [CircularOrdering(c) for c in nx.simple_cycles(self)]:
                raise ValueError("Cycle is not a cycle.")
        
        reticulations = self.reticulation_nodes()
        cutedges = self.cutedges(include_trivial=True)
        
        # Find the reticulation in the cycle and move it in the first place of the ordering
        ret = list(set(reticulations) & set(cycle.order))[0]
        cycle = cycle.carrousel(ret)
        # Find all cut-edges incident to the cycle (in the right order)
        cycle_cutedges = [(u,v) for u in cycle.order for v in self.nodes if ((u,v) in cutedges or (v,u) in cutedges)]
        # Create the corresponding set-ordering of the leaves
        setordering = CircularSetOrdering([self.split_from_cutedge(u,v).set2 for (u,v) in cycle_cutedges])
        # Set of leaves below the reticulation
        ret_set = setordering.setorder[0]
        
        return setordering, ret_set
        
    @lru_cache(maxsize=1)
    def splits(self, include_trivial=True):
        """Returns a SplitSystem induced by the cutedges of the network. If
        include_trivial is False, the trivial splits are not included."""
        cutedges = self.cutedges(include_trivial=include_trivial)
        return SplitSystem([self.split_from_cutedge(u,v) for (u,v) in cutedges])
           
    def cutvertex_to_cycle(self, v, circular_setorder, reticulation_set):
        """Replaces a given cut-vertex by a cycle. The structure of the cycle
        is based on circular_setorder (which determines the cyclic ordering of the leafsets
        around the cycle) and reticulation_set (which determines which subset of leaves is 
        below the reticulation). Raises an error if (i) v is not a cut-vertex;
        (ii) circular_setorder does not match the partition induced by v; 
        (iii) reticulation_set is not one of the parts of the circular_setorder."""
        
        induced_partition, cutedge_dict = self.partition_from_cutvertex(v, return_cutedges=True)
        if validate():
            if induced_partition != Partition(circular_setorder.parts):
                raise ValueError("Circular setordering does not match partition induced by cut-vertex.")
            if not reticulation_set in circular_setorder:
                raise ValueError('Reticulation set is not part of the circular setordering.')
        circular_setorder = circular_setorder.carrousel(reticulation_set)
        self.remove_node(v)
        
        cycle_nodes = id_generator(circular_setorder.size())
        cycle_edges = [(cycle_nodes[i], cycle_nodes[i+1]) for i in range(1, len(cycle_nodes)-1)]
        ret_edges = [(cycle_nodes[1], cycle_nodes[0]), (cycle_nodes[-1],cycle_nodes[0])]
        self.add_nodes_from(cycle_nodes)
        self.add_edges_from(cycle_edges)
        self.add_directed_edges_from(ret_edges)

        for i, leaf_subset in enumerate(circular_setorder):
            (u,v) = cutedge_dict[leaf_subset]
            self.add_edge(u, cycle_nodes[i])
    
    def create_split(self, split):
        """Creates the given split in the network by splitting a cut-vertex into a 
        cutedge. Raises an error if the split does not exactly cover the leafset, 
        if the split is not compatible with the current network, or if the split 
        already exists."""
        if validate():
            if not split.elements == self.leaves:
                raise ValueError("Not a valid split for the network.")
            if split in self.splits():
                raise ValueError("Split already exists.")
        is_tree = self.is_tree()

        for v in self.internal_nodes():
            if not is_tree:
                if not self.is_cutvertex(v):
                    continue
            partition, cutedges = self.partition_from_cutvertex(v, return_cutedges=True)
            if partition.is_refinement(split):
                self.remove_node(v)
                u, w = id_generator(2)
                self.add_nodes_from([u,w])
                self.add_edge(u,w)
                for part in partition.parts:
                    cutedge = cutedges[part]
                    leaf_side_node = cutedge[0] if cutedge[1] == v else cutedge[0]
                    if part.issubset(split.set1):
                        self.add_edge(leaf_side_node, u)
                    else:
                        self.add_edge(leaf_side_node, w)
                return
        raise ValueError("Split is not compatible with the tree.")

    def contract_split(self, split):
        """Contracts the cut-edge corresponding to the given non-trivial split. Raises an error if the
        given split is not a split of the network, or if it is a trivial split."""
        
        if validate():
            if not isinstance(split, Split) or not split.elements == self.leaves:
                raise ValueError("Not a split on the leaves of the network.")
            if not split in self.splits() or split.is_trivial():
                raise ValueError("Not a non-trivial split of the network.")
            
        u, v = self.cutedge_from_split(split)
        self.identify_two_nodes(u, v)

    def split_support(self, split, quarnets):
        """Takes as input a split A|B of the network and a DenseQuarnetSet on the leaves
        of the network. Returns the ratio of quarnets on a1,a2,b1,b2 in the set that 'support/agree'
        with the given split. Raises an error if split is not a split of the network,
        or quarnets is not a correct dense quarnet set."""
        
        from .quarnetset import DenseQuarnetSet
        
        if validate():
            if not isinstance(split, Split) or split not in self.splits():
                raise ValueError("Incorrect split")
            if not isinstance(quarnets, DenseQuarnetSet) or quarnets.leaves != self.leaves:
                raise ValueError("Wrong DenseQuarnetSet")
            
        induced_qsplits = split.induced_quartetsplits(include_trivial=False)
        given_qsplits = quarnets.quartetsplits()
        
        nr_agreeing_quartet_splits = 0
        for s in induced_qsplits:
            if s in given_qsplits:
                #leaves = s.elements()
                #w = quarnets.quarnet(leaves).weight
                nr_agreeing_quartet_splits += 1
                
        return nr_agreeing_quartet_splits / len(induced_qsplits)
        
    def blobtree(self):
        """Returns the blobtree of the network."""
        blobtree = self.copy()
        for blob in blobtree.blobs():
            blobtree.identify_node_set(blob)
        return blobtree
    
    @lru_cache(maxsize=1)
    def trinets(self):
        """Returns the trinets of the network as a DenseTrinetSet.
        Raises an error if the network is not level-1 since the method is not 
        implemented for higher levels."""
        if validate():
            if self.level() > 1:
                raise ValueError("Level too high: quarnet generation only implemented for level-0 and level-1 networks.")
        
        trinets = []
        visited_leaf_sets = []
        
        for cycle in nx.simple_cycles(self):
            setordering, ret_set = self.order_from_cycle(CircularOrdering(cycle))
            for sub_order in setordering.suborderings(3):
                for repr_order in sub_order.representative_orderings():
                    reticulation_leaf_lst = list(set(repr_order.elements) & set(ret_set))                    
                    if len(reticulation_leaf_lst) == 0: # 3-star
                        trinets.append(ThreeStar(repr_order))
                    else: # 3-cycle
                        trinets.append(Triangle(repr_order, reticulation_leaf_lst[0]))
                    visited_leaf_sets.append(frozenset(repr_order))
                    
        for three_leaves in itertools.combinations(self.leaves, 3):
            if not frozenset(three_leaves) in visited_leaf_sets:
                trinets.append(ThreeStar(three_leaves))
            
        return DenseTrinetSet(trinets)
    
    @lru_cache(maxsize=1)
    def quarnets(self, triangles=True, reticulations=True, induced_4splits=None, return_4splits=False):
        """Returns the quarnets of the network, if triangles=False, the triangles are
        collapsed. If reticulations=False, the 4-cycles dont contain reticulations.
        The quarnets are returned as a DenseQuarnetSet.
        Optionally takes a QuartetSplitSystem 'induced_4splits' that contains all 2|2-splits induced
        by a split of the network. This speeds up the computation. Also allows for returning the 2|2-splits again.
        Raises an error if the network is not level-1 since the method is not 
        implemented for higher levels. Or if the QuartetSplitSystem has a different leafset."""
        if validate():
            if self.level() > 1:
                raise ValueError("Level too high: quarnet generation only implemented for level-0 and level-1 networks.")
        
        quarnets = []
        
        for cycle in nx.simple_cycles(self):
            setordering, ret_set = self.order_from_cycle(CircularOrdering(cycle))
            for sub_order in setordering.suborderings(4):
                for repr_order in sub_order.representative_orderings():
                    reticulation_leaf_lst = list(set(repr_order.elements) & set(ret_set))
                    if len(reticulation_leaf_lst) == 0: # Quartet tree
                        a, b = repr_order.order[:2]; c, d = repr_order.order[2:]
                        split = Split({a,b},{c,d})
                        if triangles:
                            quarnets.append(QuartetTree(split))
                        else:
                            quarnets.append(SplitQuarnet(split))
                    else: # 4-cycle
                        if reticulations:
                            quarnets.append(FourCycle(repr_order, reticulation_leaf_lst[0]))
                        else:
                            quarnets.append(CycleQuarnet(repr_order))

        if induced_4splits:
            if validate():
                if induced_4splits.elements != self.leaves:
                    raise ValueError
            all_4splits = induced_4splits
        else:
        
            all_4splits = self.splits(include_trivial=False).induced_quartetsplits()

        if triangles:
            trinetset = self.trinets()
        
        for subsplit in all_4splits:
            if triangles:
                a, b = subsplit.set1; c, d = subsplit.set2
                trinet1 = trinetset.trinet({a,b,c})
                trinet2 = trinetset.trinet({a,c,d})
                
                reticulations = []
                if isinstance(trinet1, Triangle):
                    ret1 = trinet1.reticulation_leaf
                    if ret1 == c:
                        ret1 = {c, d}
                    else:
                        ret1 = {ret1}
                    reticulations.append(ret1)
                if isinstance(trinet2, Triangle):
                    ret2 = trinet2.reticulation_leaf
                    if ret2 == a:
                        ret2 = {a, b}
                    else:
                        ret2 = {ret2}
                    reticulations.append(ret2)
                                
                if len(reticulations) == 0:
                    quarnets.append(QuartetTree(subsplit))
                elif len(reticulations) == 1:
                    quarnets.append(SingleTriangle(subsplit, reticulations[0]))
                elif len(reticulations) == 2:
                    quarnets.append(DoubleTriangle(subsplit, reticulations[0], reticulations[1]))

            else:
                quarnets.append(SplitQuarnet(subsplit))
        
        if return_4splits:
            return DenseQuarnetSet(quarnets), all_4splits
        else:
            return DenseQuarnetSet(quarnets)
    
    def visualize(self, layout='kamada', title=None, leaflabels=True, internal_labels=False, font_size=12):
        """Visualization function with several layout-options: ['kamada', 'neato', 'twopi', 'circo'].
        If pygraphviz is not installed, use 'kamada'. Optional title can be given to the plot.
        Label printing can be turned off."""
        if layout == 'kamada':
            pos = nx.kamada_kawai_layout(self)
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
        nx.draw_networkx_edges(self, pos, edgelist=self.undirected_edges(), edge_color='black', width=1,
                               node_size=200, alpha=0.85)
        nx.draw_networkx_edges(self, pos, edgelist=self.directed_edges, edge_color='black', width=1,
                               node_size=200, alpha=0.85, style='dashed', arrows=True, arrowstyle='->', arrowsize=11)    
        if leaflabels == True:
            nx.draw_networkx_labels(self, pos, labels={leaf:leaf for leaf in self.leaves}, font_size=font_size)   
        if internal_labels == True:
            nx.draw_networkx_labels(self, pos, labels={node:node for node in self.internal_nodes()}, font_size=font_size)   


        if title == None: title = "Semi-Directed Network"
        plt.title(title)
        plt.show()
        return ax

    def save_to_file(self, filename, overwrite=False, write_info=True):
        """Saves the network to the given file as a list of undirected
        edges, a list of directed edges and a list of leaves. If overwrite = False,
        an error is raised when the file already exists. If write_info is True,
        the first line contains some summary information as a comment."""
        if overwrite == False and os.path.exists(filename):
            raise ValueError("File already exists.")

        with open(filename, 'w') as file:
            if write_info == True:
                file.write(f"###SemiDirectedNetwork_on_{len(self.leaves)}_leaves_with_{len(self.reticulation_nodes())}_reticulations\n")

            file.write("#Undirected\n")
            for u, v in self.undirected_edges():
                file.write(f"{u} {v}\n")
            
            file.write("#Directed\n")
            for u, v in self.directed_edges:
                file.write(f"{u} {v}\n")
            
            file.write("#Leaves\n")
            for v in self.leaves:
                file.write(f"{v}\n")
    
    def create_enewick(self):
        """Creates an enewick string of the network with an arbitrary rooting."""
        for root in self.leaves:
            if self.is_rootable_at(root):
                rooted_network = self.to_directed_network(root)
                enewick_string = rooted_network.create_enewick()  
                return enewick_string
                    
    def load_from_file(self, filename):
        """Loads the network from the given file which should contain a list 
        of undirected edges, a list of directed edges, and a list of leaves. If 
        the filename does not specify a folder, the working directory is assumed. 
        Raises an error if the file is in an incorrect format."""
        
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
                elif line.startswith("#Leaves"):
                    section = "leaves"
                    continue
                
                if section == "undirected":
                    line_lst = line.split()
                    if len(line_lst) != 2:
                        raise ValueError("Incorrect format.")
                    u, v = line_lst
                    self.add_edge(u, v)
                elif section == "directed":
                    line_lst = line.split()
                    if len(line_lst) != 2:
                        raise ValueError("Incorrect format.")
                    u, v = line_lst
                    self.add_directed_edge(u, v)
                elif section == "leaves":
                    line_lst = line.split()
                    if len(line_lst) != 1:
                        raise ValueError("Incorrect format.")
                    v = line_lst[0]
                    self.add_leaf(v)

############################################################






############################################################

def random_semi_directed_network(n, k):
    """Generates a random semi-directed triangle-free level-1 network with
    n leaves and k reticulations."""
    
    if k == 0:
        # Delete a reticulation edge from a random level-1 network
        N = random_semi_directed_network(n, 1)
        ret_edge = random.choice(list(N.directed_edges))
        r1 = ret_edge[0]; r2 = ret_edge[1]
        N.remove_edge(r1, r2)
        
        n1 = list(set(N.neighbors(r1)) - N.leaves)[0]
        n2 = list(set(N.neighbors(r2)) - N.leaves)[0]

        N.identify_two_nodes(n1, r1)
        N.identify_two_nodes(n2, r2)
        
        return N       
    
    # Start with blobtree with degree-2 nodes
    graph = nx.complete_graph(n)
    spanning_tree = nx.random_spanning_tree(graph)
    spanning_tree = nx.relabel_nodes(spanning_tree, {node:id_generator() for node in spanning_tree.nodes})
    blobtree = SemiDirectedNetwork.from_graph(spanning_tree)
    for i, node in enumerate(spanning_tree.nodes):
        blobtree.add_leaf(f"{i}")
        blobtree.add_edge(f"{i}", node)
        
    # Contract degree-2 nodes
    degree2_nodes = [v for v in blobtree.nodes if blobtree.degree(v) == 2]
    for v in degree2_nodes:
        u = [u for u in blobtree.neighbors(v) if u not in blobtree.leaves][0]
        blobtree.identify_two_nodes(u, v)

    # Enforce reticulation number
    while True:
        reticulations = [v for v in blobtree.nodes if blobtree.degree(v) > 3]
        if len(reticulations) < k:
            candidates = [(u,v) for (u,v) in blobtree.edges if blobtree.degree(u) == 3 and blobtree.degree(v) == 3]
            u, v = random.choice(candidates)
            blobtree.identify_two_nodes(u, v)
        elif len(reticulations) > k:
            r1, r2 = random.sample(reticulations, 2)
            path = nx.shortest_path(blobtree, source=r1, target=r2)
            if set(path) & set(reticulations) != {r1, r2}:
                continue
            else:
                blobtree.identify_node_set(path)
        else:
            break

    # Randomly fill in cycles and reticulations
    reticulations = [v for v in blobtree.nodes if blobtree.degree(v) > 3]
    for v in reticulations:
        leaf_partition = blobtree.partition_from_cutvertex(v)
        random_order = leaf_partition.parts
        random.shuffle(random_order)
        random_circular_order = CircularSetOrdering(random_order)
        while True:
            random_reticulation_set = random.choice(random_circular_order.parts)
            new_blobtree = blobtree.copy()
            new_blobtree.cutvertex_to_cycle(v, random_circular_order, random_reticulation_set)
            if new_blobtree.is_rootable():
                break
        blobtree = new_blobtree
    
    return blobtree

