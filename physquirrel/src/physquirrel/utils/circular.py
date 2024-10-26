import itertools

from ._config import validate
from .partition import Partition

############################################################

class CircularSetOrdering(Partition):
    """
    Class for circular orderings of sets, child-class of the general Partition 
    class. Takes as input a list of sets of elements.
        self.setorder: list with the ordering of the sets
    """
    
    def __init__(self, setorder: list[set], **kwargs):
        parts = [frozenset(s) for s in setorder]
        super().__init__(parts, **kwargs)
        self.setorder = parts

    def __repr__(self):
        return f"CircularSetOrdering({[set(a) for a in self.setorder]})"
    
    def __iter__(self):
        self._index = 0
        return self
    
    def __next__(self):
        if self._index < len(self.setorder):
            result = self.setorder[self._index]
            self._index += 1
            return result
        else:
            raise StopIteration
    
    def __eq__(self, other):
        """Checks if two circular setorderings are equal up to cyclic permutations 
        and taking the reverse."""
        if not isinstance(other, CircularSetOrdering):
            return False
        
        if self.size() != other.size() or len(self) != len(other):
            return False
        
        permutations = [self.setorder[i:] + self.setorder[:i] for i in range(self.size())]
        permutations += [list(reversed(self.setorder[i:] + self.setorder[:i])) for i in range(self.size())]
        return any(permutation == other.setorder for permutation in permutations)
    
    def __hash__(self):
        """Returns a hash value based on the cyclic permutations and reversals of the setorder."""
        # hash value calculated based on the minimum lexicographically ordered 
        # tuple among all cyclic permutations and reversals of the setorder attribute
        order = tuple(tuple(sorted(elt)) for elt in self.setorder)
        cyclic_permutations = [tuple(order[i:] + order[:i]) for i in range(self.size())]
        reversed_permutations = [tuple(reversed(p)) for p in cyclic_permutations]
        return hash(min(cyclic_permutations + reversed_permutations))

    def carrousel(self, first):
        """Returns a new CircularSetOrdering where the set 'first' is moved to 
        the first place in the setorder, while keeping the respective cyclic order,
        making the order equivalent."""
        if first not in self.setorder:
            raise ValueError()
        i = self.setorder.index(first)
        new_setorder = self.setorder[i:] + self.setorder[:i]
        return CircularSetOrdering(new_setorder)
    
    def reverse(self):
        """Returns a new CircularSetOrdering where the order is reversed,
        making the order equivalent."""
        new_setorder = self.setorder[::-1]
        return CircularSetOrdering(new_setorder)
        
    def is_singletons(self):
        """Checks if the ordering consists of just singleton sets."""
        return all((len(s) == 1) for s in self.setorder)
    
    def to_circular_ordering(self):
        """Returns the circular ordering corresponding to the circular setordering, 
        if the sets are all singletons. Raises an error otherwise."""
        if not self.is_singletons():
            raise ValueError("Not only singleton sets.")
        return CircularOrdering([list(elt)[0] for elt in self.setorder])
    
    def are_neighbours(self, set1, set2):
        """Returns boolean indicating whether two sets are neighbours in the 
        circular setordering."""
        if validate():
            if set1 == set2 or set1 not in self.setorder or set2 not in self.setorder:
                raise ValueError
        
        i = self.setorder.index(set1)
        j = self.setorder.index(set2)
        if abs(i-j) == 1 or abs(i-j) == self.size()-1:
            return True
        else:
            return False
    
    def suborderings(self, size=4):
        """Returns a list of all suborderings with a specified size."""
        sublists = self._fixed_size_sub(self.setorder, size)
        return [CircularSetOrdering(lst) for lst in sublists]
    
    def representative_orderings(self):
        """Generates a list of all CircularOrderings that contain exactly one 
        element per set of the CircularSetOrdering"""
        combinations = list(itertools.product(*self.setorder))
        return [CircularOrdering(list(combination)) for combination in combinations]
    
############################################################





############################################################
   
class CircularOrdering(CircularSetOrdering):
    """
    Class for circular orderings of single elements, child-class of the 
    CircularSetOrdering class. Takes as input a list of single elements (e.g. strings, int).
        self.order: list with the ordering of the elements
    """
    
    def __init__(self, order: list, **kwargs):
        self.order = order
        parts = [frozenset({s}) for s in order]
        super().__init__(parts, **kwargs)
    
    def __repr__(self):
        return f"CircularOrdering({self.order})"
    
    def __eq__(self, other):
        if not isinstance(other, CircularOrdering):
            return False
        return super().__eq__(other)
    
    def __hash__(self):
        return super().__hash__()
     
    def __contains__(self, element):
        return element in self.order
    
    def __iter__(self):
        self._index = 0
        return self
    
    def __next__(self):
        if self._index < len(self.order):
            result = self.order[self._index]
            self._index += 1
            return result
        else:
            raise StopIteration
    
    def carrousel(self, first):
        """Returns a new CircularOrdering where the set 'first' is moved to 
        the first place in the order, while keeping the respective cyclic order,
        making the order equivalent."""
        i = self.order.index(first)
        new_order = self.order[i:] + self.order[:i]
        return CircularOrdering(new_order)
    
    def reverse(self):
        """Returns a new CircularOrdering where the order is reversed,
        making the order equivalent."""
        new_order = self.order[::-1]
        return CircularOrdering(new_order)
        
    def are_neighbours(self, elt1, elt2):
        """Returns boolean indicating whether two elements are neighbours in the 
        circular ordering."""
        return super().are_neighbours({elt1}, {elt2})
    
    def suborderings(self, size=4):
        """Returns a list of all suborderings with a specified size."""
        sublists = self._fixed_size_sub(self.order, size)
        return [CircularOrdering(lst) for lst in sublists]
    
    def to_circular_setordering(self, mapping=None):
        """Returns a circular setordering obtained from the circular ordering.
        Mapping is an optional dictionairy that maps each element (and possible not-used elements)
        to a set. If mapping = None, the elements are mapped to their singleton sets."""
        if mapping == None:
            mapping = {elt: frozenset({elt}) for elt in self.order}
        new_setorder = [mapping[elt] for elt in self.order]
        return CircularSetOrdering(new_setorder)

############################################################

