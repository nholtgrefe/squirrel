import itertools
from ._config import validate

############################################################

class Partition:
    """
    General class for partitions of sets. Takes as input a list of sets of 
    elements. Raises an error if the sets are not disjoint.
        self.elements: set containing the union of all elements
        self.parts: list containing the partition sets
    """
    
    def __init__(self, parts: list[set]):
        self.parts = [frozenset(part) for part in parts]
        self.elements = set().union(*self.parts)

        if validate():
            if not self._is_valid():
                raise ValueError("Invalid partition: sets overlap")
    
    def __repr__(self):
        unfrozen_set = [set(part) for part in self.parts]
        return f"Partition({unfrozen_set})"

    def __eq__(self, other):
        if not isinstance(other, Partition):
            return False
        return set(self.parts) == set(other.parts)
    
    def __hash__(self):
        return hash(frozenset(self.parts))
        
    def __contains__(self, subset):
        return any(subset == part for part in self.parts)
    
    def __iter__(self):
        self._index = 0
        return self
    
    def __next__(self):
        if self._index < len(self.parts):
            result = self.parts[self._index]
            self._index += 1
            return result
        else:
            raise StopIteration
        
    def __len__(self):
        """Returns the total number of elements the partition covers."""
        return len(self.elements)
    
    def size(self):
        """Returns how many parts the partition consists of."""
        return len(self.parts)
        
    def _is_valid(self):
        return sum(len(part) for part in self.parts) == len(self.elements)
    
    def subpartitions(self, size=4):
        """Returns a list of all subpartitions (i.e. containing not all sets)
        of a specified size."""
        subparts = self._fixed_size_sub(self.parts, size)
        return [Partition(part) for part in subparts]
    
    def representative_partitions(self):
        """Generates a list of all partitions that contain exactly one 
        element per set of the current partition"""
        combinations = list(itertools.product(*self.parts))
        partitions = [[{elt} for elt in comb] for comb in combinations]
        return [Partition(part) for part in partitions]
    
    def is_refinement(self, other):
        """Checks if the partition is a refinement of another partition (i.e. each part
        is a subset of a part of the other partition). Raises an error if 'other' is
        not a partition or covers different elements."""
        if not isinstance(other, Partition):
            raise ValueError("The argument must be an instance of Partition")
        if not self.elements == other.elements:
            raise ValueError("Other partition covers different elements.")
        for part in self.parts:
            if not any(part.issubset(other_part) for other_part in other.parts):
                return False
        return True
    
    @staticmethod
    def _fixed_size_sub(structure, k):
        """Helper function: returns all sublists/subsets/other types of size k 
        of a given 'structure'."""
        t = type(structure)
        return [t(comb) for comb in itertools.combinations(structure, k)]
