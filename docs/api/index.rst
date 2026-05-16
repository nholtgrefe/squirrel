API Reference
=============

All public symbols are importable directly from ``physquirrel``. This page maps
names to their detailed documentation sections.

.. list-table:: Where to find what
   :widths: 40 60
   :header-rows: 1

   * - Symbol(s)
     - Page
   * - :class:`~physquirrel.datatypes.sqprofile.SqQuartetProfile`,
       :class:`~physquirrel.datatypes.sqprofileset.SqQuartetProfileSet`
     - :doc:`datatypes`
   * - :func:`~physquirrel.datatypes.io.to_psq`,
       :func:`~physquirrel.datatypes.io.from_psq`,
       :func:`~physquirrel.datatypes.io.to_profile_list`,
       :func:`~physquirrel.datatypes.io.from_profile_list`
     - :doc:`io`
   * - :func:`~physquirrel.algorithms.squirrel.squirrel`,
       :func:`~physquirrel.algorithms.delta_heuristic.delta_heuristic`
     - :doc:`pipeline`
   * - :func:`~physquirrel.algorithms.tstar_tree.tstar_tree`,
       :func:`~physquirrel.algorithms.tstar_tree.bstar`,
       :func:`~physquirrel.algorithms.qjoining.quartet_joining`,
       :func:`~physquirrel.algorithms.qjoining.adapted_quartet_joining`,
       :func:`~physquirrel.algorithms.unresolve_tree.unresolve_tree`,
       :func:`~physquirrel.algorithms.unresolve_tree.split_support`,
       :func:`~physquirrel.algorithms.cycle_resolution.resolve_cycles`,
       :func:`~physquirrel.algorithms.qsimilarity.sqprofileset_from_network`,
       :func:`~physquirrel.algorithms.qsimilarity.sqprofileset_similarity`
     - :doc:`algorithms`
   * - :func:`~physquirrel.convenience.delta_heuristic_from_msa`,
       :func:`~physquirrel.convenience.squirrel_from_distances`,
       :func:`~physquirrel.convenience.squirrel_from_msa`
     - :doc:`convenience`
   * - :func:`~physquirrel.algorithms.tsp.optimal_tsp_tour`,
       :func:`~physquirrel.algorithms.tsp.approximate_tsp_tour`
     - :doc:`tsp`

.. toctree::
   :maxdepth: 2
   :caption: Modules

   datatypes
   io
   pipeline
   algorithms
   convenience
   tsp
