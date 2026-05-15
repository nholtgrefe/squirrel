import math

from .utils import validate
from .trinet import Trinet, ThreeStar, Triangle

############################################################

class TrinetSet:
    """
    Class for a set of trinets. Takes as optional input a list of Trinet instances.
    Raises an error if there are multiple trinets on the same set of leaves or if
    the argument is not a set of trinets.
        self.leaves: union of all leaves of the trinets
        self.trinets: list of the trinets
    """
    
    def __init__(self, trinets=None):
        self.trinets = set(trinets) if trinets else set()
        if validate():
            if not self._is_valid():
                raise ValueError("Argument needs to be a list of trinets with different leafsets.")          
        self.leaves = self._extract_leaf_set()
        self._trinet_dict = {frozenset(t.leaves):t for t in self.trinets}
    
    def __repr__(self):
        string = ""
        for t in self.trinets:
            string = string + str(t) + ",\n"
        return "TrinetSet[\n" + string + "]"
    
    def __str__(self):
        string = ""
        for t in self.trinets:
            string = string + str(t) + ", "
        return "TrinetSet[" + string + "]"
    
    def __contains__(self, trinet):        
        if not frozenset(trinet.leaves) in self._trinet_dict.keys():
            return False
        return trinet == self._trinet_dict[frozenset(trinet.leaves)]
    
    def __iter__(self):
        self._index = 0
        return self
    
    def __next__(self):
        if self._index < len(self.trinets):
            trinet_list = list(self.trinets)
            result = trinet_list[self._index]
            self._index += 1
            return result
        else:
            raise StopIteration
            
    def __len__(self):
        """Returns the total number trinets in the set."""
        return len(self.trinets)
    
    def add_trinet(self, trinet):
        """Add a trinet to the trinetSet. Raises an error if there is already a
        trinet on the same four leaves."""
        if validate():
            if not isinstance(trinet, trinet) or frozenset(trinet.leaves) in self._trinet_dict.keys():
                raise ValueError("Argument needs to be a trinet on new set of three leaves.")
        self.trinets = self.trinets | {trinet}
        self.leaves = self.leaves | trinet.leaves
        self._trinet_dict[frozenset(trinet.leaves)] = trinet
    
    def add_trinets_from(self, trinets):
        """Adds all trinets in the list to the trinetSet. Raises an error if 
        this results in multiple trinets on the same leafsets."""
        for trinet in trinets:
            self.add_trinet(trinet)
    
    def trinet(self, three_leaves):
        """Returns the trinet in the trinet set on the given set of three leaves.
        Raises an error if there is no such trinet in the set."""
        if validate():
            if frozenset(three_leaves) not in self._trinet_dict.keys():
                raise ValueError("No trinet on given leaves")
        return self._trinet_dict[frozenset(three_leaves)]
    
    def to_densetrinetset(self):
        """Returns the trinetSet as a DenseTrinetSet. Raises an error if it is
        not dense."""
        return DenseTrinetSet(self)
    
    def nr_threestars(self):
        """Returns the number of ThreeStar trinets in the set."""
        return sum(1 for t in self.trinets if isinstance(t, ThreeStar))

    def nr_triangles(self):
        """Returns the number of Triangles in the set."""
        return sum(1 for t in self.trinets if isinstance(t, Triangle))

    def _is_valid(self):
        if not all(isinstance(t, Trinet) for t in self.trinets):
            return False
        return len(self.trinets) == len({frozenset(t.leaves) for t in self.trinets})
    
    def _extract_leaf_set(self):
        """Returns the union of all leaves of the trinet set."""
        leaf_set = set()
        for trinet in self.trinets:
            leaf_set.update(trinet.leaves)
        return leaf_set
    
    def _is_dense(self):
        "Check if the trinet set is dense with respect to the union of leaves."
        # From the _is_valid method we already know that there are no
        # trinets with the same leafset. So we just need to check the number of quarnets.
        if self._is_valid() == False:
            return False
        return len(self.trinets) == math.comb(len(self.leaves), 3)
    
############################################################






############################################################

class DenseTrinetSet(TrinetSet):
    """
    Class for a dense set of trinet; immutable child class of TrinetSet. 
    Takes as input a list of Trinet instances or a TrinetSet. Raises an error 
    if the set does not form a dense set for the union of all leaves.
        self.leaves: union of all leaves of the trinets
        self.trinets: list of the trinets
    """
    
    def __init__(self, trinets):
        if isinstance(trinets, TrinetSet):
            trinets = trinets.trinets
        super().__init__(trinets)
        if validate():
            if self._is_dense() == False:
                raise ValueError("Not a dense trinet set.")
    
    def __repr__(self):
        return "Dense" + super().__repr__()
    
    def __str__(self):
        return "Dense" + super().__str__()
    
    def consistency(self, trinetset):
        """Returns the ratio of trinets consistent with a trinet in another trinetset (C-measure).
        Raises an error if the other trinetsets has a different leafset or is not dense.
        See: 'A Practical Algorithm for Reconstructing Level-1 Phylogenetic Networks' by Huber et al."""
        if validate():
            if not isinstance(trinetset, DenseTrinetSet):
                raise ValueError
            if self.leaves != trinetset.leaves:
                raise ValueError("Wrong leafset.")
        return len(self.trinets & trinetset.trinets) / math.comb(len(self.leaves),3)

    def distance(self, trinetset, normalize=True):
        """Returns the S-distance (symmetric difference) w.r.t. another trinetset. The distance can optionally
        be normalized. Raises an error if the other trinetset has a different leafset or is not dense.
        See: 'A Practical Algorithm for Reconstructing Level-1 Phylogenetic Networks' by Huber et al."""
        if validate():
            if not isinstance(trinetset, DenseTrinetSet):
                raise ValueError
            if self.leaves != trinetset.leaves:
                raise ValueError("Wrong leafset.")
        if normalize:
            return len(self.trinets ^ trinetset.trinets) / len(self.trinets | trinetset.trinets)
        else:
            return len(self.trinets ^ trinetset.trinets)

    def add_trinet(self, trinet):
        raise AttributeError("DenseTrinetSet is an immutable class.")

    def add_trinets_from(self, trinets):
        raise AttributeError("DenseTrinetSet is an immutable class.")

############################################################
