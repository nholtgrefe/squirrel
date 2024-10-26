import networkx as nx
import itertools, random, math, time, os
import numpy as np
from functools import lru_cache
#from python_tsp.exact import solve_tsp_dynamic_programming
from typing import Dict, List, Optional, Tuple


from .utils import CircularOrdering, id_generator, validate

from .quarnet import SingleTriangle, DoubleTriangle, QuartetTree, SplitQuarnet, FourCycle, CycleQuarnet, Quarnet
from .splits import Split, QuartetSplitSet

############################################################

class QuarnetSet:
    """
    Class for a set of quarnets. Takes as optional input a list of Quarnet instances.
    Raises an error if there are multiple quarnets on the same set of leaves or if
    the argument is not a set of quarnets.
        self.leaves: union of all leaves of the quarnets
        self.quarnets: list of the quarnets
    """
    
    def __init__(self, quarnets=None):
        self.quarnets = set(quarnets) if quarnets else set()
        if validate():
            if not self._is_valid():
                raise ValueError("Argument needs to be a list of quarnets with different leafsets.")          
        self.leaves = self._extract_leaf_set()
        self._quarnet_dict = {frozenset(q.leaves):q for q in self.quarnets}
    
    def __repr__(self):
        string = ""
        for q in self.quarnets:
            string = string + str(q) + ",\n"
        return "QuarnetSet[\n" + string + "]"
    
    def __str__(self):
        string = ""
        for q in self.quarnets:
            string = string + str(q) + ", "
        return "QuarnetSet[" + string + "]"
    
    def __contains__(self, quarnet):        
        if not frozenset(quarnet.leaves) in self._quarnet_dict.keys():
            return False
        return quarnet == self._quarnet_dict[frozenset(quarnet.leaves)]
    
    def __iter__(self):
        self._index = 0
        return self
    
    def __next__(self):
        if self._index < len(self.quarnets):
            quarnet_list = list(self.quarnets)
            result = quarnet_list[self._index]
            self._index += 1
            return result
        else:
            raise StopIteration
            
    def __len__(self):
        """Returns the total number quarnets in the set."""
        return len(self.quarnets)
    
    def weight(self):
        """Returns the total weight of the quarnets in the set."""
        return sum(q.weight for q in self.quarnets)
    
    def is_weighted(self):
        """Returns whether the set contains quarnets with a weight unequal to 1."""
        return self.weight() != len(self)
    
    def add_quarnet(self, quarnet):
        """Add a quarnet to the QuarnetSet. Raises an error if there is already a
        quarnet on the same four leaves."""
        if validate():
            if not isinstance(quarnet, Quarnet) or frozenset(quarnet.leaves) in self._quarnet_dict.keys():
                raise ValueError("Argument needs to be a quarnet on new set of four leaves.")
        self.quarnets = self.quarnets | {quarnet}
        self.leaves = self.leaves | quarnet.leaves
        self._quarnet_dict[frozenset(quarnet.leaves)] = quarnet
    
    def add_quarnets_from(self, quarnets):
        """Adds all quarnets in the list to the QuarnetSet. Raises an error if 
        this results in multiple quarnets on the same leafsets."""
        for quarnet in quarnets:
            self.add_quarnet(quarnet)
    
    def quarnet(self, four_leaves):
        """Returns the quarnet in the quarnet set on the given set of four leaves.
        Raises an error if there is no such quarnet in the set."""
        if validate():
            if frozenset(four_leaves) not in self._quarnet_dict.keys():
                raise ValueError("No quarnet on given leaves")
        return self._quarnet_dict[frozenset(four_leaves)]
    
    def to_densequarnetset(self):
        """Returns the QuarnetSet as a DenseQuarnetSet. Raises an error if it is
        not dense."""
        return DenseQuarnetSet(self)
    
    def nr_cyclequarnets(self):
        """Returns the number of CycleQuarnet quarnets in the set."""
        return sum(1 for q in self.quarnets if isinstance(q, CycleQuarnet))
    
    def nr_fourcycles(self):
        """Returns the number of FourCycle quarnets in the set."""
        return sum(1 for q in self.quarnets if isinstance(q, FourCycle))

    def nr_splitquarnets(self):
        """Returns the number of SplitQuarnets in the set."""
        return sum(1 for q in self.quarnets if isinstance(q, SplitQuarnet))
    
    def nr_quartettrees(self):
        """Returns the number of QuartetTrees in the set."""
        return sum(1 for q in self.quarnets if isinstance(q, QuartetTree))
    
    def nr_singletriangles(self):
        """Returns the number of SingleTriangles in the set."""
        return sum(1 for q in self.quarnets if isinstance(q, SingleTriangle))
    
    def nr_doubletriangles(self):
        """Returns the number of DoubleTriangles in the set."""
        return sum(1 for q in self.quarnets if isinstance(q, DoubleTriangle))
        
    def contains_triangles(self):
        """Returns whether the set contains SplitQuarnets (then
        it returns False) or also QT, ST or DT quarnets (then returns True)."""
        return self.nr_quartettrees() + self.nr_singletriangles() + self.nr_doubletriangles() == self.nr_splitquarnets()
    
    def contains_reticulations(self):
        """Returns whether the set contains FourCycle quarnets (then
        it returns True) or only CycleQuarnets (then returns False)."""
        return self.nr_fourcycles() == self.nr_cyclequarnets()
    
    def quartetsplits(self, threshold = 1.0):
        """Returns a QuartetSplitSet object (on the complete leafset) containing 
        all non-trivial quartet-splits induced by the quarnets in the QuarnetSet that have
        a weight of atleast 'threshold'."""
        splitset = []
        for quarnet in self.quarnets:
            if quarnet.weight >= threshold:
                if isinstance(quarnet, SplitQuarnet):
                    splitset.append(quarnet.split)
        return QuartetSplitSet(splitset, elements=self.leaves)
    
    def collapse_triangles(self):
        """Returns a QuarnetSet containing SplitQuarnets. The
        QuartetTrees, SingleTriangles, DoubleTriangles are contracted to SplitQuarnets."""
        new_quarnets = []
        for q in self.quarnets:
            if isinstance(q, FourCycle):
                new_quarnets.append(q)
            else:
                new_quarnets.append(SplitQuarnet(q.split, weight=q.weight))
        return QuarnetSet(new_quarnets)
        
    def delete_reticulations(self):
        """Returns a QuarnetSet containing CycleQuarnets. The FourCycles are
        reduced."""
        new_quarnets = []
        for q in self.quarnets:
            if isinstance(q, FourCycle):
                new_quarnets.append(CycleQuarnet(q.circular_order, weight=q.weight))
            else:
                new_quarnets.append(q)
        return QuarnetSet(new_quarnets)
    
    def tstar_tree(self):
        """Returns the T* tree of the quarnet set. (Note that if the QuarnetSet is
        not dense, the other quarnets are all assumed to have no non-trivial split."""

        quartetsplit_set = self.quartetsplits()
        bstar_splitsystem = quartetsplit_set.bstar()
        tree = bstar_splitsystem.displayed_tree()
        return tree
        
    def _is_valid(self):
        if not all(isinstance(q, Quarnet) for q in self.quarnets):
            return False
        return len(self.quarnets) == len({frozenset(q.leaves) for q in self.quarnets})
    
    def _extract_leaf_set(self):
        """Returns the union of all leaves of the quarnet set."""
        leaf_set = set()
        for quarnet in self.quarnets:
            leaf_set.update(quarnet.leaves)
        return leaf_set

    def _is_dense(self):
        "Check if the quarnet set is dense with respect to the union of leaves."
        # From the _is_valid method we already know that there are no
        # quarnets with the same leafset. So we just need to check the number of quarnets.
        if self._is_valid() == False:
            return False
        return len(self.quarnets) == math.comb(len(self.leaves), 4)
    
############################################################






############################################################

class DenseQuarnetSet(QuarnetSet):
    """
    Class for a dense set of quarnets; immutable child class of QuarnetSet. 
    Takes as input a list of Quarnet instances or a QuarnetSet. Raises an error 
    if the set does not form a dense set for the union of all leaves.
        self.leaves: union of all leaves of the quarnets
        self.quarnets: list of the quarnets
    """
    
    def __init__(self, quarnets):
        if isinstance(quarnets, QuarnetSet):
            quarnets = quarnets.quarnets
        if isinstance(quarnets, str):
            quarnets = self._load_from_file(quarnets)
        super().__init__(quarnets)

        if validate():
            if self._is_dense() == False:
                raise ValueError("Not a dense quarnet set.")
    
    def __repr__(self):
        return "Dense" + super().__repr__()
    
    def __str__(self):
        return "Dense" + super().__str__()
    
    @lru_cache(maxsize=1)
    def contains_triangles(self):
        """Returns whether the set contains SplitQuarnets (then
        it returns False) or also QT, ST or DT quarnets (then returns True)."""
        return super().contains_triangles()
    
    @lru_cache(maxsize=1)
    def contains_reticulations(self):
        """Returns whether the set contains FourCycles quarnets (then
        it returns True) or only CycleQuarnets (then returns False)."""
        return super().contains_reticulations()
    
    @lru_cache(maxsize=1)
    def is_weighted(self):
        """Returns whether the set contains quarnets with a weight unequal to 1."""
        return super().is_weighted()
    
    def _load_from_file(self, filename):
        """Loads a (triangle-free) DenseQuarnetSet from a file. One line per quarnet.
        'SQ: a b c d 0.5' gives a splitquarnet ab|cd with weight 0.5.
        '4C: a b c d 0.5' gives a four-cycle a-b-c-d with a below the reticulation and weight 0.5.
        if no weight is specified the quarnet gets weight 1.0."""
        
        if not os.path.dirname(filename):
            filename = os.path.join(os.getcwd(), filename)
    
        quarnets = []
        with open(filename, 'r') as file:
            for line in file:
                parts = line.split()
                if len(parts) < 5 or len(parts) > 6:
                    raise ValueError('Invalid format.')
                quarnet_type = parts[0]
                a, b, c, d = parts[1], parts[2], parts[3], parts[4]
                w = float(parts[5]) if len(parts) == 6 else 1.0
                if quarnet_type == 'SQ':
                    q = SplitQuarnet(Split({a,b},{c,d}), weight=w)
                elif quarnet_type == '4C':
                    q = FourCycle(CircularOrdering([a,b,c,d]), a, weight=w)
                quarnets.append(q)
                
        return quarnets
    
    def split_support(self, split):
        """Given a split A|B, returns the (weighted) ratio of quarnets on {a1 a2 b1 b2} that
        agree with the split. Raises an error if the split is trivial or on the wrong leafset."""
        
        if validate():
            if not isinstance(split, Split):
                raise ValueError
            if split.is_trivial():
                raise ValueError
            if split.elements != self.leaves:
                raise ValueError
        
        support = 0
        total_weight = 0
        for a1, a2 in itertools.combinations(split.set1, 2):
            for b1, b2 in itertools.combinations(split.set2, 2):
                q = self.quarnet({a1, a2, b1, b2})
                total_weight += q.weight
                if isinstance(q, SplitQuarnet):
                    if q.split == Split({a1, a2}, {b1, b2}):
                        support += q.weight
        if total_weight == 0:
            return 0
        return support / total_weight
        
    def collapse_triangles(self):
        """Returns a DenseQuarnetSet containing only FourCycles and SplitQuarnets. The
        QuartetTrees, SingleTriangles, DoubleTriangles are contracted to SplitQuarnets."""
        new_quarnets = []
        for q in self.quarnets:
            if isinstance(q, CycleQuarnet):
                new_quarnets.append(q)
            else:
                new_quarnets.append(SplitQuarnet(q.split, weight=q.weight))
        return DenseQuarnetSet(new_quarnets)
    
    def delete_reticulations(self):
        """Returns a DenseQuarnetSet containing CycleQuarnets. The FourCycles are
        reduced."""
        new_quarnets = []
        for q in self.quarnets:
            if isinstance(q, FourCycle):
                new_quarnets.append(CycleQuarnet(q.circular_order, weight=q.weight))
            else:
                new_quarnets.append(q)
        return DenseQuarnetSet(new_quarnets)
    
    def similarity(self, quarnetset, measure="QC", triangles=False, reticulations=True):
        """Returns a similarity measure between two DenseQuarnetSets. Either accepts
        QC or QS. Triangles indicates whether triangles should be accounted for, reticulations indicate whether reticulations in four-cycles
        should be accounted for. For QC returns the weighted version."""
        
        if measure not in ["QC", "QS"]:
            raise ValueError
        if measure == "QC":
            return self.consistency(quarnetset, triangles=triangles, reticulations=reticulations)
        elif measure == "QS":
            return 1 - self.distance(quarnetset, triangles=triangles, reticulations=reticulations)        
    
    def consistency(self, quarnetset, triangles=False, reticulations=True, weighted=True):
        """Returns the C-measure of the quarnetset (i.e. ratio of consistent quarnets).  
        Raises an error if the other quarnetsets have a different leafset or are not dense.
        The optional weighted parameter allows for a weighted version of the C-measure. Triangles indicates
        whether triangles should be accounted for, reticulations indicate whether reticulations in four-cycles
        should be accounted for.
        See: 'A Practical Algorithm for Reconstructing Level-1 Phylogenetic Networks' by Huber et al."""
        
        if validate():
            if not isinstance(quarnetset, DenseQuarnetSet):
                raise ValueError
            if self.leaves != quarnetset.leaves:
                raise ValueError("Wrong leafset.")        
                
        if triangles:
            if not self.contains_triangles() or not quarnetset.contains_triangles():
                raise ValueError("Both sets need to contain triangles.")
        if reticulations:
            if not self.contains_reticulations() or not quarnetset.contains_reticulations():
                raise ValueError("Both sets need to contain reticulations in 4-cycles.")
        
        if not self.is_weighted(): weighted = False

        own_quarnets = self
        other_quarnets = quarnetset
        if not triangles:
            own_quarnets = own_quarnets.collapse_triangles()
            other_quarnets = other_quarnets.collapse_triangles()
        if not reticulations:
            own_quarnets = own_quarnets.delete_reticulations()
            other_quarnets = other_quarnets.delete_reticulations()
            
        if weighted:
            r1 = sum(q.weight for q in own_quarnets if q in other_quarnets)
            r2 = self.weight()
            if r2 == 0: return 1
        else:
            r1 = len(own_quarnets.quarnets & other_quarnets.quarnets)
            r2 = len(own_quarnets.quarnets)
            
        return r1 / r2
    
    def consistency_prime(self, quarnetset, reference, triangles=False, reticulations=True, weighted=True):
        """Returns the C'-measure of the quarnetset w.r.t. the quarnetset reference.  
        If in either of the quarnetsets triangles are collapsed, they will be in both.        
        Raises an error if the other quarnetsets have a different leafset or are not dense.
        The optional weighted parameter allows for a weighted version of the C'-measure.
        Triangles indicates whether triangles should be accounted for, reticulations indicate whether reticulations in four-cycles
        should be accounted for.
        See: 'A Practical Algorithm for Reconstructing Level-1 Phylogenetic Networks' by Huber et al."""
        
        if validate():
            if not isinstance(quarnetset, DenseQuarnetSet) or not isinstance(reference, DenseQuarnetSet):
                raise ValueError
            if self.leaves != quarnetset.leaves or self.leaves != reference.leaves:
                raise ValueError("Wrong leafset.")        
        
        if triangles:
            if not self.contains_triangles() or not quarnetset.contains_triangles() or not reference.contains_triangles():
                raise ValueError("All sets need to contain triangles.")
        if reticulations:
            if not self.contains_reticulations() or not quarnetset.contains_reticulations() or not reference.contains_reticulations():
                raise ValueError("All sets need to contain reticulations in 4-cycles.")

        if not self.is_weighted(): weighted = False
        
        own_quarnets = self
        other_quarnets = quarnetset
        reference_quarnets = reference
        if not triangles:
            own_quarnets = own_quarnets.collapse_triangles()
            other_quarnets = other_quarnets.collapse_triangles()
            reference_quarnets = reference_quarnets.collapse_triangles()
        if not reticulations:
            own_quarnets = own_quarnets.delete_reticulations()
            other_quarnets = other_quarnets.delete_reticulations()
            reference_quarnets = reference_quarnets.delete_reticulations()
        
        if weighted:
            r1 = sum(q.weight for q in own_quarnets if q in reference_quarnets and q in other_quarnets)
            r2 = sum(q.weight for q in own_quarnets if q in reference_quarnets)
            if r2 == 0: return 1
        else:
            r1 = len(own_quarnets.quarnets & reference_quarnets.quarnets & other_quarnets.quarnets)
            r2 = len(own_quarnets.quarnets & reference_quarnets.quarnets)
        
        return r1 / r2
    
    def distance(self, quarnetset, triangles=False, reticulations=True, normalize=True):
        """Returns the S-distance (symmetric difference) w.r.t. another quarnetset. The distance can optionally
        be normalized. Raises an error if the other quarnetset has a different leafset or is not dense.
        Triangles indicates whether triangles should be accounted for, reticulations indicate whether reticulations in four-cycles
        should be accounted for.
        See: 'A Practical Algorithm for Reconstructing Level-1 Phylogenetic Networks' by Huber et al."""
        
        if validate():
            if not isinstance(quarnetset, DenseQuarnetSet):
                raise ValueError
            if self.leaves != quarnetset.leaves:
                raise ValueError("Wrong leafset.")
                
        if triangles:
            if not self.contains_triangles() or not quarnetset.contains_triangles():
                raise ValueError("Both sets need to contain triangles.")
        if reticulations:
            if not self.contains_reticulations() or not quarnetset.contains_reticulations():
                raise ValueError("Both sets need to contain reticulations in 4-cycles.")
        
        own_quarnets = self
        other_quarnets = quarnetset
        if not triangles:
            own_quarnets = own_quarnets.collapse_triangles()
            other_quarnets = other_quarnets.collapse_triangles()
        if not reticulations:
            own_quarnets = own_quarnets.delete_reticulations()
            other_quarnets = other_quarnets.delete_reticulations()
            
        if normalize:
            return len(own_quarnets.quarnets ^ other_quarnets.quarnets) / len(own_quarnets.quarnets | other_quarnets.quarnets)
        else:
            return len(own_quarnets.quarnets ^ other_quarnets.quarnets)

    def quartet_joining(self, threshold=0.0, starting_tree=None):
        """Returns the tree from the quartet-joining algorithm for the quarnet-splits. 
        The threshold can be used for early stopping. In particular, if the highest
        score in an iteration is less than threshold times the maximum possible 
        score, the algorithm stops. If threshold is 0.0, the algorithm continues 
        until we have a binary tree, thus recreating the exact QuartetJoining algorithm. 
        (Note that if the QuarnetSet is not dense, the other quarnets are all 
         assumed to have no non-trivial split.) Can also be started with a partially
        resolved tree, instead of the star-graph"""
        
        threshold = min(1, max(0, threshold))
    
        from .sdnetwork import SemiDirectedNetwork
                
        if starting_tree is None:
            tree = SemiDirectedNetwork()
            center_node = id_generator()
            tree.add_node(center_node)
            tree.add_leaves_from(self.leaves)
            tree.add_edges_from([(center_node,leaf) for leaf in self.leaves])
        else:
            if validate():
                if not isinstance(starting_tree, SemiDirectedNetwork):
                    raise ValueError
                if not starting_tree.is_tree():
                    raise ValueError
                if not starting_tree.leaves == self.leaves:
                    raise ValueError
            tree = starting_tree
                
        center_nodes = [u for u in tree.internal_nodes() if tree.degree(u) > 3]
        nr_iterations = len(self.leaves) - 3 - len(tree.splits(include_trivial=False))
        for i in range(nr_iterations):
            score = dict()
            for center_node in center_nodes:
                if tree.degree(center_node) == 3:
                    continue
                C = dict()
                for u in tree.neighbors(center_node):
                    C[u] = tree.split_from_cutedge(u, center_node).set1
                            
                neighbour_pairs = itertools.combinations(tree.neighbors(center_node), 2)
                for (u1 ,u2) in neighbour_pairs:
                    score[(u1,u2, center_node)] = 0
                    for (A1, A2) in itertools.combinations(C.values(), 2):
                        if C[u1] not in [A1, A2] and C[u2] not in [A1, A2]:
                            score[(u1,u2, center_node)] += self._omega_bar(frozenset(C[u1]), frozenset(C[u2]), frozenset(A1), frozenset(A2))
        
            u1star, u2star, cstar = max(score, key=score.get)
            
            max_possible = 2 * math.comb(len(C.keys())-2, 2)
    
            # Check threshold
            if score[(u1star, u2star, cstar)] < threshold * max_possible:
                return tree
            
            # Otherwise create split
            tree.remove_edges_from([(cstar, u1star), (cstar, u2star)])
            w = id_generator()
            tree.add_node(w)
            tree.add_edges_from([(w,u1star), (w,u2star), (w,cstar)])
                        
        return tree

    def TSP_ordering(self, method='SA', threshold=13):
        """Returns a CircularOrdering of all leaves with the TSP method. 
            threshold = maximum size for which to solve TSP optimally (if None, always solved optimally)
            method = heuristic to use if size is larger than threshold.
                one of ['SA', 'greedy', 'christofides']"""
                
        if threshold is None:
            threshold = 99999999999

        if len(self.leaves) <= threshold:    
            leaf_order = list(self.leaves)
            distance_matrix = np.zeros((len(self.leaves),len(self.leaves)))
            for quarnet in self.quarnets:
                for (leaf1, leaf2) in itertools.combinations(quarnet.leaves, 2):
                    i = leaf_order.index(leaf1)
                    j = leaf_order.index(leaf2)
                    delta = quarnet.kalmanson_distance(leaf1, leaf2)
                    distance_matrix[i][j] += delta
                    distance_matrix[j][i] += delta
    
            permutation, _ = _solve_tsp_dynamic_programming(distance_matrix)
            tour = [leaf_order[i] for i in permutation]
            
        else:
            complete_graph = nx.Graph()
            complete_graph.add_nodes_from(self.leaves)
            weight_zero_edges = [(u,v,0) for (u,v) in itertools.combinations(self.leaves, 2)]
            complete_graph.add_weighted_edges_from(weight_zero_edges)
            
            for quarnet in self.quarnets:
                for (leaf1, leaf2) in itertools.combinations(quarnet.leaves, 2):
                    delta = quarnet.kalmanson_distance(leaf1, leaf2)
                    complete_graph[leaf1][leaf2]['weight'] += delta
            
            tsp = nx.approximation.traveling_salesman_problem
            if method == "SA":               
                method = nx.approximation.simulated_annealing_tsp
                tour = tsp(complete_graph, cycle=False, method=method, init_cycle='greedy', seed=0)

            elif method == "greedy":
                method = nx.approximation.greedy_tsp
                tour = tsp(complete_graph, cycle=False, method=method)

            else:
                method = nx.approximation.christofides
                tour = tsp(complete_graph, cycle=False, method=method)

            #weight_matrix = nx.floyd_warshall_numpy(complete_graph)

        return CircularOrdering(tour)
    

    def squirrel(self, outgroup=None, method="best", include_score=False, all_networks=False, verbose=False, visualize=False, tsp_threshold=13, **kwargs):
        """Returns a semi-directed level-1 triangle-free network built up from the
        dense quarnet set. Repeatedly contracts the least supported edge in the
        quartet-joining tree, and then builds the network. If method="best", the
        network with the highest similarity score is returned, if method="first",
        the last network for which the score did not drop is returned. If an outgroup is specified
        the method returns DirectedNetworks rooted at the outgroup.
            include_score: whether to return the similarity measure
            verbose: whether to print intermediate info
            visualize: whether to plot intermediate networks
            all_networks: whether to return all networks instead of the best one (ordered according to score)
            tsp_threshold: threshold up to where to solve tsp exactly
            kwargs are passed to the similarity function
        """
        
        if method not in ["best", "first"]:
            raise ValueError("Incorrect method.")
        

        tstar = self.tstar_tree()
        qj_tree = self.quartet_joining(starting_tree=tstar)
        split_supports = dict()
        for split in qj_tree.splits(include_trivial=False):
            split_supports[split] = self.split_support(split)
        sorted_splits = [k for k, v in sorted(split_supports.items(), key=lambda item: item[1])]

        tree0 = qj_tree
        quarnets0, qsplitsystem = tree0.quarnets(return_4splits=True)
        score0 = self.similarity(quarnets0, **kwargs)

        trees = [tree0]; networks = [tree0]; scores = [score0]

        for i, split in enumerate(sorted_splits):
            if visualize: networks[-1].visualize(title=f"Net {i}")
            u, v = trees[-1].cutedge_from_split(split)
            split_induced_4splits = trees[-1].quartetsplits_from_cutedge(u, v)

            tree_new = trees[-1].copy()
            tree_new.contract_split(split)
            net_new = self.reconstruct_network_from_tree(tree_new, outgroup=outgroup, tsp_threshold=tsp_threshold)
            
            qsplitsystem = QuartetSplitSet({s for s in qsplitsystem if s not in split_induced_4splits})
            quarnets_new = net_new.quarnets(induced_4splits=qsplitsystem)

            score_new = self.similarity(quarnets_new, **kwargs)

            trees.append(tree_new); networks.append(net_new); scores.append(score_new)

            if verbose == True:
                print(f"\rIteration {i+1}: score = {score_new}", end='')
                print("")
                time.sleep(.01)
            
            if method == "first" and score_new < scores[-2]:
                break

        if visualize: networks[-1].visualize(title=f"Net {i}")

        best_network_index = scores.index(max(scores))
        if verbose == True:
            print(f"\rBest network from iteration {best_network_index+1}: score = {scores[best_network_index]}", end='')
        
        if outgroup:
            networks = [network.to_directed_network(outgroup) for network in networks]
            
        if all_networks == False:
            if include_score == False:
                return networks[best_network_index]
            else:
                return networks[best_network_index], scores[best_network_index]
        
        elif all_networks == True:
            paired_lists = list(zip(scores, networks))
            paired_lists.sort(reverse=True, key=lambda x: x[0])
            sorted_scores, sorted_networks = zip(*paired_lists)
            if include_score == False:
                return list(sorted_networks)
            else:
                return list(sorted_networks), list(sorted_scores)
        
    def reconstruct_network_from_tree(self, blobtree=None, outgroup=None, tsp_threshold=13):
        """Returns a semi-directed level-1 triangle-free network built up from the
        dense quarnet set. 'blobtree' is an optional argument specifying the tree
        used to fill in the cycles. If blobtree=None, the T*-tree is used. Raises an 
        error if the supplied blobtree is not a tree."""
        if blobtree == None:
            network = self.tstar_tree()
        else:
            network = blobtree.copy()
        if validate():
            if not network.is_tree():
                raise ValueError("Need a tree as input.")

        large_degree_vertices = sorted(
            [v for v in network.internal_nodes() if network.degree(v) > 3],
            key=lambda v: network.degree(v),
            reverse=True
        )
        for v in large_degree_vertices:
            leaf_partition = network.partition_from_cutvertex(v)
            repr_quarnets, set_mapping = self._vote_quarnets(leaf_partition)
            repr_circular_order = repr_quarnets.TSP_ordering(threshold=tsp_threshold)
            inv_set_mapping = {rep:s for s, rep in set_mapping.items()}
            circular_setorder = repr_circular_order.to_circular_setordering(inv_set_mapping)

            repr_reticulation_leaves = repr_quarnets._vote_reticulation()
            
            for i, ret in enumerate(repr_reticulation_leaves):
                reticulation_set = inv_set_mapping[ret]
                new_network = network.copy()
                new_network.cutvertex_to_cycle(v, circular_setorder, reticulation_set)
                if outgroup is None:
                    if new_network.is_rootable():
                        break
                else:
                    if new_network.is_rootable_at(outgroup):
                        break

            network = new_network
            
        return network
 
    def _vote_reticulation(self):
        """Returns a sorted list of the leaves, where the first leaf appears in most four cycles,
        the second leaf is in second place, etc. If the set only contains a single quarnet 4-cycle, its
        reticulation is placed first, and the rest randomly. If the set contains a single splitquarnet:
        completely random."""
        
        votes = {leaf: 0 for leaf in self.leaves}
        
        if self.__len__() == 1:
            q = self.quarnet(self.leaves)
            if isinstance(q, FourCycle):
                votes[q.reticulation_leaf] += q.weight
                        
        else:
            for quarnet in self.quarnets:
                if isinstance(quarnet, CycleQuarnet):                    
                    for leaf in quarnet.leaves:
                        votes[leaf] += quarnet.weight
                        

        return sorted(votes, key=votes.get)[::-1]
    
    def _vote_quarnets(self, leaf_partition):
        """Given some partition of the leaves, we represent each set in the partition
        by some representative leaf. Then, we create a DenseQuarnetSet on those 
        representative leaves, where each quarnet-structure is determined by voting.
        Also returns a dictionairy mapping each leafset of the partition to its
        representative leaf. Raises an error if the input is not a valid leaf-partition."""
        if validate():
            if not leaf_partition.elements == self.leaves:
                raise ValueError("Not a partition of the leafset")
        set_mapping = {leafset: id_generator() for leafset in leaf_partition}
        mapping = {leaf: value for key, value in set_mapping.items() for leaf in key}
        repr_quarnets = []
        for four_sets in leaf_partition.subpartitions(4):
            repr_quarnet = self._vote_quarnet(four_sets, mapping)
            repr_quarnets.append(repr_quarnet)
        return DenseQuarnetSet(repr_quarnets), set_mapping
 
    def _vote_quarnet(self, leaf_subpartition, mapping):
        """Given are a leaf_subpartition (i.e. their union need not be the complete
        leaf set) of size 4, and a mapping that maps each of the sets to some label.
        Returns a quarnet on those labels, the structure is determined by the most
        appearing quarnet for all quarnets induced by the leaf_subpartition."
        Raises an error if the input is an incorrect leaf-subpartition or the mapping
        does not map every leaf in it."""
        
        if validate():
            if leaf_subpartition.size() != 4:
                raise ValueError("Leaf subpartition is not of size 4.")
            if not leaf_subpartition.elements.issubset(set(mapping.keys())):
                raise ValueError("Not all leaves are mapped.")
            if not leaf_subpartition.elements.issubset(self.leaves):
                raise ValueError("Not a subpartition of the leafset")
        
        votes = []
        ret_counts = {}
        
        # *a*bcd and ab*c*d are the same. The ret only matters for the ret in the end.
        for four_leafset in leaf_subpartition.representative_partitions():
            q = self.quarnet(four_leafset.elements)
            if isinstance(q, FourCycle):
                relabelled_ret = mapping[q.reticulation_leaf]
                q = CycleQuarnet(q.circular_order)
                
                relabelled_q = q.relabel(mapping)
                votes.append(relabelled_q)
                
                if relabelled_q in ret_counts.keys():
                    ret_counts[relabelled_q].append(relabelled_ret)
                else:
                    ret_counts[relabelled_q] = [relabelled_ret]
        
            if isinstance(q, QuartetTree) or isinstance(q, SingleTriangle) or isinstance(q, DoubleTriangle):
                q = SplitQuarnet(q.split)
                relabelled_q = q.relabel(mapping)
                votes.append(relabelled_q)
                
            else:
                relabelled_q = q.relabel(mapping)
                votes.append(relabelled_q)
            
        counts = {candidate: 0 for candidate in set(votes)}
        for vote in votes:
            counts[vote] += vote.weight
        
        final_q = max(counts, key=counts.get)
        sum_count_weight = sum(list(counts.values()))
        
        if sum_count_weight == 0:
            final_weight = 0
        else:
            final_weight = counts[final_q] / sum_count_weight
        
        if isinstance(final_q, CycleQuarnet):
            rets = ret_counts[final_q]
            ret = max(set(rets), key=rets.count)
            final_q = FourCycle(final_q.circular_order, ret)
        
        final_q.weight = final_weight

        return final_q
    
        #### Version where every 4-cycle is counted individually, so where reticulation matters.
        # for four_leafset in leaf_subpartition.representative_partitions():
        #     q = self.quarnet(four_leafset.elements)
        #     if isinstance(q, QuartetTree) or isinstance(q, SingleTriangle) or isinstance(q, DoubleTriangle):
        #         q = SplitQuarnet(q.split)
        #     relabelled_q = q.relabel(mapping)
        #     votes.append(relabelled_q)
        # counts = {candidate: 0 for candidate in set(votes)}
        # for vote in votes:
        #     counts[vote] += vote.weight
        # return max(counts, key=counts.get)
        
    def _shake(self, epsilon=0.5, include_score=False):
        """Returns a perturbed DenseQuarnetSet. 'epsilon' is a parameter between 
        0 and 1, which determines what percentage of quarnets are changed.
        This method is for testing purposes only. If include_score = True the consistency is also returned
        , showing the ratio of quarnets in the perturbed set that are also in the original set."""
        if self.contains_triangles():
            raise ValueError
        
        epsilon = max(0, min(1, epsilon))
        nr_perturbed_quarnets = math.floor(epsilon * len(self))
        perturbed_indices = random.sample(range(len(self)), nr_perturbed_quarnets)        
        
        new_quarnets = []
        for i, q in enumerate(self.quarnets):
            if i in perturbed_indices:
                new_q = self._random_quarnet(q.leaves)
            else:
                new_q = q.copy()
            new_quarnets.append(new_q)
        res = DenseQuarnetSet(new_quarnets)
        
        if include_score == False:
            return res
        else:
            return res, self.consistency(res, triangles=False)
        
    def _random_quarnet(self, four_leaves):
        """Returns a random quarnet on the given set of four_leaves according to 
        a uniform distribution of the six types. The quarnet in the set on the four leaves is excluded.
        Raises an error if the argument is not a size-4 set, the leaves are not a subset of 
        self.leaves, or exclude is not a quarnet on the same leaves. Also raises an error
        if the quarnet set contains triangles."""
        if validate():
            if self.contains_triangles():
                raise ValueError
            if len(four_leaves) != 4:
                raise ValueError
            if not four_leaves.issubset(self.leaves):
                raise ValueError
        
        a, b, c, d = list(four_leaves)
        distr = {
             SplitQuarnet(Split({a,b},{c,d})): 1.0,
             SplitQuarnet(Split({a,c},{b,d})): 1.0,
             SplitQuarnet(Split({a,d},{b,c})): 1.0,
             FourCycle(CircularOrdering([a,b,c,d]), a): 0.25,
             FourCycle(CircularOrdering([a,b,c,d]), b): 0.25,
             FourCycle(CircularOrdering([a,b,c,d]), c): 0.25,
             FourCycle(CircularOrdering([a,b,c,d]), d): 0.25,
             FourCycle(CircularOrdering([a,c,b,d]), a): 0.25,
             FourCycle(CircularOrdering([a,c,b,d]), b): 0.25,
             FourCycle(CircularOrdering([a,c,b,d]), c): 0.25,
             FourCycle(CircularOrdering([a,c,b,d]), d): 0.25,
             FourCycle(CircularOrdering([a,b,d,c]), a): 0.25,
             FourCycle(CircularOrdering([a,b,d,c]), b): 0.25,
             FourCycle(CircularOrdering([a,b,d,c]), c): 0.25,
             FourCycle(CircularOrdering([a,b,d,c]), d): 0.25,
        
             }
        
        q = self.quarnet(four_leaves)
        
        distr_new = dict()
        if isinstance(q, SplitQuarnet):
            distr_new = {key: value for key, value in distr.items() if not key == q}
        
        elif isinstance(q, FourCycle):
            distr_new = dict()
            for key, value in distr.items():
                if isinstance(key, SplitQuarnet):
                    distr_new[key] = value
                elif isinstance(key, FourCycle):
                    if key.circular_order != q.circular_order:
                        distr_new[key] = value
        
        
        choice = random.choices(list(distr_new.keys()), weights=list(distr_new.values()), k=1)[0]
        return choice
    
    @lru_cache(maxsize=1)
    def _omega_bar(self, X1, X2, X3, X4):
        """Returns omega_bar as defined in the paper:
        "Consistency of the QNet algorithm for generating planar split networks
        from weighted quartets" by Grunewald et al."""
        res = 0
        total_sum = 0
        for (x1, x2, x3, x4) in itertools.product(X1, X2, X3, X4):
            q = self.quarnet({x1, x2, x3, x4})
            total_sum += q.weight
            if isinstance(q, SplitQuarnet):
                if (q.split.set1 == {x1, x2} and q.split.set2 == {x3, x4}) \
                  or (q.split.set2 == {x1, x2} and q.split.set1 == {x3, x4}):
                    res += q.weight
        if total_sum == 0:
            return 1
        return res / total_sum
    
    def add_quarnet(self, quarnet):
        raise AttributeError("DenseQuarnetSet is an immutable class.")

    def add_quarnets_from(self, quarnets):
        raise AttributeError("DenseQuarnetSet is an immutable class.")
   
############################################################


def _solve_tsp_dynamic_programming(
    distance_matrix: np.ndarray,
    maxsize: Optional[int] = None,
) -> Tuple[List, float]:
    """
    Taken from the 'python_tsp' package source code (under the MIT license).
    
    Solve TSP to optimality with dynamic programming

    Parameters
    ----------
    distance_matrix
        Distance matrix of shape (n x n) with the (i, j) entry indicating the
        distance from node i to j. It does not need to be symmetric

    maxsize
        Parameter passed to ``lru_cache`` decorator. Used to define the maximum
        size for the recursion tree. Defaults to `None`, which essentially
        means "take as much space as needed".

    Returns
    -------
    permutation
        A permutation of nodes from 0 to n that produces the least total
        distance

    distance
        The total distance the optimal permutation produces

    Notes
    -----
    Algorithm: cost of the optimal path
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Consider a TSP instance with 3 nodes: {0, 1, 2}. Let dist(0, {1, 2}) be the
    distance from 0, visiting all nodes in {1, 2} and going back to 0. This can
    be computed recursively as:

        dist(0, {1, 2}) = min(
            c_{0, 1} + dist(1, {2}),
            c_{0, 2} + dist(2, {1}),
        )

    wherein c_{0, 1} is the cost from going from 0 to 1 in the distance matrix.
    The inner dist(1, {2}) is computed as:

        dist(1, {2}) = min(
            c_{1, 2} + dist(2, {}),
        )

    and similarly for dist(2, {1}). The stopping point in the recursion is:

        dist(2, {}) = c_{2, 0}.

    This process can be generalized as:

        dist(ni, N) =   min   ( c_{ni, nj} + dist(nj, N - {nj}) )
                      nj in N

    and

        dist(ni, {}) = c_{ni, 0}

    With starting point as dist(0, {1, 2, ..., tsp_size}). The notation
    N - {nj} is the difference operator, meaning set N without node nj.


    Algorithm: compute the optimal path
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    The previous process returns the distance of the optimal path. To find the
    actual path, we need to store in a memory the following key/values:

        memo[(ni, N)] = nj_min

    with nj_min the node in N that provided the smallest value of dist(ni, N).
    Then, the process goes backwards starting from
    memo[(0, {1, 2, ..., tsp_size})].

    In the previous example, suppose memo[(0, {1, 2})] = 1.
    Then, look for memo[(1, {2})] = 2.
    Then, since the next step would be memo[2, {}], stop there. The optimal
    path would be 0 -> 1 -> 2 -> 0.

    Reference
    ---------
    https://en.wikipedia.org/wiki/Held%E2%80%93Karp_algorithm#cite_note-5
    """
    # Get initial set {1, 2, ..., tsp_size} as a frozenset because @lru_cache
    # requires a hashable type
    N = frozenset(range(1, distance_matrix.shape[0]))
    memo: Dict[Tuple, int] = {}

    # Step 1: get minimum distance
    @lru_cache(maxsize=maxsize)
    def dist(ni: int, N: frozenset) -> float:
        if not N:
            return distance_matrix[ni, 0]

        # Store the costs in the form (nj, dist(nj, N))
        costs = [
            (nj, distance_matrix[ni, nj] + dist(nj, N.difference({nj})))
            for nj in N
        ]
        nmin, min_cost = min(costs, key=lambda x: x[1])
        memo[(ni, N)] = nmin

        return min_cost

    best_distance = dist(0, N)

    # Step 2: get path with the minimum distance
    ni = 0  # start at the origin
    solution = [0]

    while N:
        ni = memo[(ni, N)]
        solution.append(ni)
        N = N.difference({ni})

    return solution, best_distance
