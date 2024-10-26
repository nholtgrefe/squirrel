from .utils import CircularOrdering, validate

from .splits import Split

############################################################
# Can also be implemented as subclass of the SemiDirectedNetwork class
# but that brings redundant overhead.

class Quarnet:
    """Class for quarnets with optional weight in [0,1].""" 
    def __init__(self, weight=1.0):
        self.weight = weight
        if validate():
            if not self._is_valid():
                raise ValueError("Invalid weight.")
    
    def _is_valid(self):
        if self.weight > 1.0 or self.weight < 0.0:
            raise ValueError
            
    def set_weight(self, weight):
        if validate():
            if weight > 1.0 or weight < 0.0:
                raise ValueError
        self.weight = weight

############################################################





############################################################

class CycleQuarnet(Quarnet):
    """Class for quarnets with a cyclic ordering. This is a subclass of the Quarnet class.
    (Used primarily to disregard the reticulation in 4-cycles.) Takes as input 
    an instance of the CircularOrdering class indicating the order of the four leaves.
        self.circular_order: CircularOrdering of the leaves with the reticulation leaf in
                                first place"""
    
    def __init__(self, circular_order, **kwargs):
        self.circular_order = circular_order
        self.leaves = self.circular_order.elements
        super().__init__(**kwargs)
        
        if validate():
            if not self._is_valid():
                raise ValueError("Invalid cycle quarnet.")
    
    def __repr__(self):
        a, b, c, d = self.circular_order.order
        return f"Cycle quarnet with circular order {a}-{b}-{c}-{d}" 

    def __str__(self):
        a, b, c, d = self.circular_order.order
        return f"(CY: {a}-{b}-{c}-{d})" 
    
    def __hash__(self):
        return hash(self.circular_order)
    
    def __eq__(self, other):
        if not isinstance(other, CycleQuarnet):
            return False
        return self.circular_order == other.circular_order
    
    def copy(self):
        """Returns a copy of the quarnet."""
        return CycleQuarnet(self.circular_order, weight=self.weight)
    
    def kalmanson_distance(self, leaf1, leaf2, weighted=True):
        """Returns the distance between two leaves of the cycle, in light of the
        kalmanson property. Optionally alters the distances if the quarnet is weighted.
        Raises an error if the leaves are equal or are not leaves
        of the quartet."""
        if validate():
            if leaf1 not in self.leaves or leaf2 not in self.leaves or leaf1 == leaf2:
                raise ValueError()
        
        if weighted:
            c = (1 - self.weight) / 2
        else:
            c = 0
                
        if self.circular_order.are_neighbours(leaf1, leaf2):
            return 1 + c
        else:
            return 2 - c
        
    def relabel(self, mapping):
        """Returns a copy of the quarnet with a relabelling of the leaves specified
        by 'mapping': a dictionairy that maps each leaf to a new label. The dictionairy
        is allowed to contain more keys than the leaves of the quarnet. Raises 
        an error if the four leaves are not mapped to four distinct new labels."""
        sub_mapping = {leaf: mapping[leaf] for leaf in mapping.keys() if leaf in self.leaves}
        if validate():
            if len(set(sub_mapping.values())) < 4:
                raise ValueError

        new_circular_order = [sub_mapping[leaf] for leaf in self.circular_order]
        return CycleQuarnet(CircularOrdering(new_circular_order), weight=self.weight)
    
    def _is_valid(self):
        if super()._is_valid() == False:
            return False
        if len(self.circular_order) != 4:
            return False
        return True
    
    
############################################################





############################################################

class FourCycle(CycleQuarnet):
    """
    Class for 4-cycle quarnets. This is a subclass of the CycleQuarnet class.
    Takes as input an instance of the CircularOrdering class indicating the 
    order of the four leaves, and a 'reticulation_leaf', indicating the leaf below
    the reticulation. Raises an error if the order is not of length 4 or the reticulation
    leaf is not in the order. After initialization, the order is 
    rotated to have the reticulation leaf in first place.
        self.reticulation_leaf: leaf below the reticulation
        self.circular_order: CircularOrdering of the leaves with the reticulation leaf in
                                first place
    """

    def __init__(self, circular_order, reticulation_leaf, **kwargs):
        self.reticulation_leaf = reticulation_leaf
        super().__init__(circular_order, **kwargs)
        self.circular_order = self.circular_order.carrousel(reticulation_leaf) # move reticulation to place 1

        if validate():
            if not self._is_valid():
                raise ValueError("Invalid 4-cycle.")
    
    def __repr__(self):
        a, b, c, d = self.circular_order.order
        return f"4-cycle with circular order *{a}*-{b}-{c}-{d}" 

    def __str__(self):
        a, b, c, d = self.circular_order.order
        return f"(4C: *{a}*-{b}-{c}-{d})" 
    
    def __hash__(self):
        return hash((self.circular_order, self.reticulation_leaf))
    
    def __eq__(self, other):
        if not isinstance(other, FourCycle):
            return False
        return self.circular_order == other.circular_order \
            and self.reticulation_leaf == other.reticulation_leaf
    
    def copy(self):
        """Returns a copy of the quarnet."""
        return FourCycle(self.circular_order, self.reticulation_leaf, weight=self.weight)
        
    def relabel(self, mapping):
        """Returns a copy of the quarnet with a relabelling of the leaves specified
        by 'mapping': a dictionairy that maps each leaf to a new label. The dictionairy
        is allowed to contain more keys than the leaves of the quarnet. Raises 
        an error if the four leaves are not mapped to four distinct new labels."""
        sub_mapping = {leaf: mapping[leaf] for leaf in mapping.keys() if leaf in self.leaves}
        if validate():
            if len(set(sub_mapping.values())) < 4:
                raise ValueError
        new_reticulation_leaf = sub_mapping[self.reticulation_leaf]
        new_circular_order = [sub_mapping[leaf] for leaf in self.circular_order]
        return FourCycle(CircularOrdering(new_circular_order), new_reticulation_leaf, weight=self.weight)
    
    def _is_valid(self):
        if super()._is_valid() == False:
            return False
        if self.reticulation_leaf not in self.circular_order.elements:
            return False
        return True
    
    # def _construct(self):
    #     a, b, c, d = self.circular_order.order
    #     int1, int2, int3, int4 = id_generator(4)
    #     self.add_leaves_from([a,b,c,d])
    #     self.add_nodes_from([int1, int2, int3, int4])
    #     self.add_edges_from([(a,int1),(b,int2),(c,int3),(d,int4),(int1,int2),(int2,int3),(int3,int4),(int4,int1)])
    #     self.add_directed_edges_from([(int2,int1),(int4,int1)])
    
############################################################





############################################################

class SplitQuarnet(Quarnet):
    """Class for quarnets that contain a non-trivial split. 
    (Used primarily to disregard triangles in the quarnets.)
        self.split: the split of the quarnet."""
    
    def __init__(self, split, **kwargs):
        self.split = split
        self.leaves = self.split.elements
        super().__init__(**kwargs)
        
        if validate():
            if not self._is_valid():
                raise ValueError("Invalid quarnet.")
    
    def __repr__(self):
        a, b = self.split.set1; c, d = self.split.set2
        return f"Split-quarnet with split {a} {b} | {c} {d}"
    
    def __str__(self):
        a, b = self.split.set1; c, d = self.split.set2
        return f"(SQ: {a} {b} | {c} {d})"
    
    def __eq__(self, other):
        if not isinstance(other, SplitQuarnet):
            return False
        return self.split == other.split
    
    def __hash__(self):
        return hash(self.split)

    def copy(self):
        """Returns a copy of the quarnet."""
        return SplitQuarnet(self.split, weight=self.weight)
    
    def kalmanson_distance(self, leaf1, leaf2, weighted=True):
        """Returns the distance between two leaves of the quartet, in light of the
        kalmanson property. Optionally alters the distances if the quarnet is weighted.
        Raises an error if the leaves are equal or are not leaves of the quartet."""
        if validate():
            if leaf1 not in self.leaves or leaf2 not in self.leaves or leaf1 == leaf2:
                raise ValueError()
        
        if weighted:
            c = (1 - self.weight) / 2
        else:
            c = 0
        
        if not Split({leaf1}, {leaf2}).is_subsplit(self.split):
            return 1 + c
        else:
            return 2 - c
        
    def relabel(self, mapping):
        """Returns a copy of the quarnet with a relabelling of the leaves specified
        by 'mapping': a dictionairy that maps each leaf to a new label. The dictionairy
        is allowed to contain more keys than the leaves of the quarnet. Raises 
        an error if the four leaves are not mapped to four distinct new labels."""
        sub_mapping = {leaf: mapping[leaf] for leaf in mapping.keys() if leaf in self.leaves}
        if validate():
            if len(set(sub_mapping.values())) < 4:
                raise ValueError
        newset1 = {sub_mapping[leaf] for leaf in self.split.set1}
        newset2 = {sub_mapping[leaf] for leaf in self.split.set2}
        return SplitQuarnet(Split(newset1, newset2), weight=self.weight)        

    def _is_valid(self):
        if super()._is_valid() == False:
            return False
        if len(self.split.set1) != 2 or len(self.split.set2) != 2:
            return False
        return True
    
############################################################







############################################################

class QuartetTree(SplitQuarnet):
    """
    Class for (resolved) quartet trees; subclass of the SplitQuarnet class.
    Takes as input a Split indicating the split of the quartet.
    Raises an error if the given split is not a 2-2 split
        self.split: the given split
    """

    def __init__(self, split, **kwargs):
        super().__init__(split, **kwargs)
        
        #if validate():
        #    if not self._is_valid():
        #        raise ValueError("Invalid quartet.")
    
    def __repr__(self):
        a, b = self.split.set1; c, d = self.split.set2
        return f"Quartet tree with split {a} {b} | {c} {d}"
    
    def __str__(self):
        a, b = self.split.set1; c, d = self.split.set2
        return f"(QT: {a} {b} | {c} {d})"
    
    def __eq__(self, other):
        if not isinstance(other, QuartetTree):
            return False
        return self.split == other.split
    
    def __hash__(self):
        return hash(self.split)

    def copy(self):
        """Returns a copy of the quarnet."""
        return QuartetTree(self.split, weight=self.weight)
        
    def relabel(self, mapping):
        """Returns a copy of the quarnet with a relabelling of the leaves specified
        by 'mapping': a dictionairy that maps each leaf to a new label. The dictionairy
        is allowed to contain more keys than the leaves of the quarnet. Raises 
        an error if the four leaves are not mapped to four distinct new labels."""
        sub_mapping = {leaf: mapping[leaf] for leaf in mapping.keys() if leaf in self.leaves}
        if validate():
            if len(set(sub_mapping.values())) < 4:
                raise ValueError
        newset1 = {sub_mapping[leaf] for leaf in self.split.set1}
        newset2 = {sub_mapping[leaf] for leaf in self.split.set2}
        return QuartetTree(Split(newset1, newset2), weight=self.weight)        
    
    # def _construct(self):
    #     a, b = self.split.set1
    #     c, d = self.split.set2
    #     int1, int2 = id_generator(2)
    #     self.add_leaves_from([a,b,c,d])
    #     self.add_nodes_from([int1, int2])
    #     self.add_edges_from([(a,int1),(b, int1),(int1,int2),(c,int2),(d,int2)])

############################################################







############################################################

class SingleTriangle(SplitQuarnet):
    """
    Class for single triangle quarnet; subclass of the SplitQuarnet class.
    Takes as input a Split indicating the split of the quarnet, and a set of
    1 or 2 leaves indicating the leaves below the reticulation of the quarnet.
    Raises an error if the given split is not a 2-2 split, or if the 2 leaves below
    the reticulation are not on the same side of the split
        self.split: the given split
        self.reticulation_leaves: set of leaves below the reticulation
        self.type: 1 if the reticulation points towards a leaf, 2 if points towards the split edge
    """

    def __init__(self, split, reticulation_leaves, **kwargs):
        self.reticulation_leaves = reticulation_leaves
        super().__init__(split, **kwargs)
        if validate():
            if not self._is_valid():
                raise ValueError("Invalid quarnet.")
        self.type = int(len(self.reticulation_leaves))
    
    def __repr__(self):
        a, b = self.split.set1; c, d = self.split.set2
        return f"Single Triangle with split {a} {b} | {c} {d}"
    
    def __str__(self):
        a, b = self.split.set1; c, d = self.split.set2
        return f"(ST: {a} {b} | {c} {d})"
    
    def __eq__(self, other):
        if not isinstance(other, SingleTriangle):
            return False
        return self.split == other.split and self.reticulation_leaves == other.reticulation_leaves
    
    def __hash__(self):
        return hash((self.split, frozenset(self.reticulation_leaves)))

    def copy(self):
        """Returns a copy of the quarnet."""
        return SingleTriangle(self.split, self.reticulation_leaves, weight=self.weight)
        
    def relabel(self, mapping):
        """Returns a copy of the quarnet with a relabelling of the leaves specified
        by 'mapping': a dictionairy that maps each leaf to a new label. The dictionairy
        is allowed to contain more keys than the leaves of the quarnet. Raises 
        an error if the four leaves are not mapped to four distinct new labels."""
        sub_mapping = {leaf: mapping[leaf] for leaf in mapping.keys() if leaf in self.leaves}
        if validate():
            if len(set(sub_mapping.values())) < 4:
                raise ValueError
        newset1 = {sub_mapping[leaf] for leaf in self.split.set1}
        newset2 = {sub_mapping[leaf] for leaf in self.split.set2}
        reticulations = {sub_mapping[leaf] for leaf in self.reticulation_leaves}
        return SingleTriangle(Split(newset1, newset2), reticulations, weight=self.weight)        

    def _is_valid(self):
        if super()._is_valid() == False:
            return False
        if not self.reticulation_leaves.issubset(self.split.set1) and not self.reticulation_leaves.issubset(self.split.set2):
            return False
        return True

############################################################







############################################################

class DoubleTriangle(SplitQuarnet):
    """
    Class for double triangle quarnet; subclass of the SplitQuarnet class.
    Takes as input a Split indicating the split of the quarnet, and two sets of leaves
    indicating the leaves below the reticulations of the quarnet.
    Raises an error if the given split is not a 2-2 split, or if the leaves below the reticulation are
    either not on the same side of the split or the reticulations point towards each other.
        self.split: the given split
        self.reticulation_set1: set of leaves below reticulation 1
        self.reticulation_set2: set of leaves below reticulation 2
        self.type: 1 if both reticulation point towards leaves, 2 if one points towards the split edge
    """

    def __init__(self, split, reticulation_set1, reticulation_set2, **kwargs):
        self.reticulation_set1 = reticulation_set1
        self.reticulation_set2 = reticulation_set2
        super().__init__(split, **kwargs)

        if validate():
            if not self._is_valid():
                raise ValueError("Invalid quarnet.")
                
        if len(self.reticulation_set1) == 1 and len(self.reticulation_set2) == 1:
            self.type = 1
        else:
            self.type = 2

    def __repr__(self):
        a, b = self.split.set1; c, d = self.split.set2
        return f"Double Triangle with split {a} {b} | {c} {d}"
    
    def __str__(self):
        a, b = self.split.set1; c, d = self.split.set2
        return f"(DT: {a} {b} | {c} {d})"
    
    def __eq__(self, other):
        if not isinstance(other, DoubleTriangle):
            return False
        if not self.split == other.split:
            return False
        return (self.reticulation_set1 == other.reticulation_set1 and self.reticulation_set2 == other.reticulation_set2) \
            or (self.reticulation_set1 == other.reticulation_set2 and self.reticulation_set2 == other.reticulation_set1)

    def __hash__(self):
        rets = {frozenset(self.reticulation_set1), frozenset(self.reticulation_set2)}
        return hash((self.split, frozenset(rets)))

    def copy(self):
        """Returns a copy of the quarnet."""
        return DoubleTriangle(self.split, self.reticulation_set1, self.reticulation_set2, weight=self.weight)
        
    def relabel(self, mapping):
        """Returns a copy of the quarnet with a relabelling of the leaves specified
        by 'mapping': a dictionairy that maps each leaf to a new label. The dictionairy
        is allowed to contain more keys than the leaves of the quarnet. Raises 
        an error if the four leaves are not mapped to four distinct new labels."""
        sub_mapping = {leaf: mapping[leaf] for leaf in mapping.keys() if leaf in self.leaves}
        if validate():
            if len(set(sub_mapping.values())) < 4:
                raise ValueError
        newset1 = {sub_mapping[leaf] for leaf in self.split.set1}
        newset2 = {sub_mapping[leaf] for leaf in self.split.set2}
        reticulations1 = {sub_mapping[leaf] for leaf in self.reticulation_set1}
        reticulations2 = {sub_mapping[leaf] for leaf in self.reticulation_set2}
        return DoubleTriangle(Split(newset1, newset2), reticulations1, reticulations2, weight=self.weight)        

    def _is_valid(self):
        if super()._is_valid() == False:
            return False
        if self.reticulation_set1 == self.reticulation_set2:
            return False
        if not (self.reticulation_set1.issubset(self.split.set1) or self.reticulation_set1.issubset(self.split.set2)):
            return False
        if not (self.reticulation_set2.issubset(self.split.set1) or self.reticulation_set2.issubset(self.split.set2)):
            return False
        if len(self.reticulation_set1) == 2 and len(self.reticulation_set2) == 2:
            return False
        return True

############################################################


