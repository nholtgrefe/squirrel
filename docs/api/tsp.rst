TSP Solvers
===========

Traveling Salesman Problem solvers used internally by the cycle resolution
step to find circular orderings of taxon sets around reticulation cycles.

Both functions return a :class:`~phylozoo.core.primitives.circular_ordering.CircularOrdering`
over the labels of the input :class:`~phylozoo.DistanceMatrix`.

.. automodule:: physquirrel.algorithms.tsp
   :members:
   :exclude-members: _held_karp_numba
