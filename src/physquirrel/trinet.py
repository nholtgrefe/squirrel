
from .utils import validate

############################################################
# Can also be implemented as subclass of the SemiDirectedNetwork class
# but that brings redundant overhead.

class Trinet:
    """Class for trinets""" 
    def __init__(self):
        pass
    
############################################################





############################################################

class ThreeStar(Trinet):
    """
    Class for 3-star trinets. This is a subclass of the Trinet class.
    Takes as input a set of three leaves. Raises an error otherwise.
    if the order is not of length 4 or the reticulation.
        self.leaves = leafset
    """
    def __init__(self, leaves):
        self.leaves = set(leaves)
        
        if validate():
            if not self._is_valid():
                raise ValueError("Invalid Trinet.")

    def __repr__(self):
        a, b, c = self.leaves
        return f"3-star on {a}-{b}-{c}" 

    def __str__(self):
        a, b, c = self.leaves
        return f"(3S: {a}-{b}-{c})"
    
    def __hash__(self):
        return hash(frozenset(self.leaves))
    
    def __eq__(self, other):
        if not isinstance(other, ThreeStar):
            return False
        return self.leaves == other.leaves
    
    def copy(self):
        """Returns a copy of the trinet."""
        return ThreeStar(self.leaves)
      
    def _is_valid(self):
        return len(self.leaves) == 3
    
############################################################





############################################################

class Triangle(Trinet):
    """
    Class for triangle trinets. This is a subclass of the Trinet class.
    Takes as input a set of three leaves and the leaf below the reticulation.
    Raises an error if that leaf is not in the leafset or there are not 3 leaves.
        self.leaves = leafset
        self.reticulation_leaf = the leaf below the reticulation
    """
    def __init__(self, leaves, reticulation_leaf):
        self.leaves = set(leaves)
        self.reticulation_leaf = reticulation_leaf
        
        if validate():
            if not self._is_valid():
                raise ValueError("Invalid Trinet.")

    def __repr__(self):
        lst = list(self.leaves)
        i = lst.index(self.reticulation_leaf)
        order = lst[i:] + lst[:i]
        a, b, c = order
        return f"Triangle on *{a}*-{b}-{c}" 

    def __str__(self):
        lst = list(self.leaves)
        i = lst.index(self.reticulation_leaf)
        order = lst[i:] + lst[:i]
        a, b, c = order
        return f"(Tr: *{a}*-{b}-{c})"
    
    def __hash__(self):
        return hash((frozenset(self.leaves), self.reticulation_leaf))
    
    def __eq__(self, other):
        if not isinstance(other, Triangle):
            return False
        return self.leaves == other.leaves and self.reticulation_leaf == other.reticulation_leaf
    
    def copy(self):
        """Returns a copy of the trinet."""
        return Triangle(self.leaves, self.reticulation_leaf)
      
    def _is_valid(self):
        return len(self.leaves) == 3 and self.reticulation_leaf in self.leaves
