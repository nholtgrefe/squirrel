Internal Algorithm Steps
========================

The individual algorithm components that make up the Squirrel pipeline.
These are exposed at the top level of physquirrel and can be used independently
for custom workflows.

T* tree and B* splits
----------------------

.. automodule:: physquirrel.algorithms.tstar_tree
   :members:

Quartet joining
---------------

.. automodule:: physquirrel.algorithms.qjoining
   :members:
   :exclude-members: _omega_bar

Unresolve tree and split support
---------------------------------

.. automodule:: physquirrel.algorithms.unresolve_tree
   :members:

Cycle resolution
-----------------

.. automodule:: physquirrel.algorithms.cycle_resolution
   :members:
   :exclude-members: _qprofiles_to_circular_ordering, _qprofiles_to_hybrid_ranking, _insert_cycle

Profile similarity and scoring
-------------------------------

.. automodule:: physquirrel.algorithms.qsimilarity
   :members:
   :exclude-members: _circular_orders_from_cycles

Quartet distance
----------------

.. automodule:: physquirrel.algorithms.qdistance
   :members:
