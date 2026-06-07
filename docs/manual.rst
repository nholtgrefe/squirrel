Manual
======

`physquirrel` reconstructs level-1 phylogenetic networks from quartet data using
the Squirrel algorithm. This guide covers the full workflow: preparing input
data, running the algorithm, using convenience shortcuts, and working with the
output network. Note that physquirrel is built on top of `PhyloZoo`; see the
`PhyloZoo documentation <https://nholtgrefe.github.io/phylozoo/>`_ for details
on working with the core datatypes (:class:`~phylozoo.core.network.sdnetwork.sd_phynetwork.SemiDirectedPhyNetwork`,
:class:`~phylozoo.core.split.base.Split`, :class:`~phylozoo.core.quartet.base.Quartet`, :class:`~phylozoo.core.distance.base.DistanceMatrix`,
and :class:`~phylozoo.core.sequence.base.MSA`). The main datatypes from phylozoo used in physquirrel are
re-exported in the top-level namespace for convenience.

The standard import convention used throughout is:

.. code-block:: python

   import physquirrel as psq

Input data
----------

physquirrel works with three input types, depending on how far upstream you
want to start:

.. list-table::
   :widths: 30 30 40
   :header-rows: 1

   * - Input
     - Type
     - Entry point
   * - Quarnets (4-leaf subnetworks)
     - :class:`~physquirrel.datatypes.sqprofileset.SqQuartetProfileSet`
     - :func:`~physquirrel.algorithms.squirrel.squirrel`
   * - Pairwise distance matrix
     - :class:`~phylozoo.core.distance.base.DistanceMatrix`
     - :func:`~physquirrel.convenience.squirrel_from_distances`
   * - Multiple-sequence alignment
     - :class:`~phylozoo.core.sequence.base.MSA`
     - :func:`~physquirrel.convenience.squirrel_from_msa`

From a set of quarnets
~~~~~~~~~~~~~~~~~~~~~~~

A quarnet is a 4-leaf subnetwork that can be either a split (tree-like) or a
4-cycle (network-like). In physquirrel, quarnets are represented by
:class:`~physquirrel.datatypes.sqprofile.SqQuartetProfile` objects, which are then collected into a
:class:`~physquirrel.datatypes.sqprofileset.SqQuartetProfileSet`. These objects build on top of
phylozoo's :class:`~phylozoo.core.quartet.qprofile.QuartetProfile`: the
standard object to represent quartets.

For full control, construct an :class:`~physquirrel.datatypes.sqprofileset.SqQuartetProfileSet` from
individual :class:`~physquirrel.datatypes.sqprofile.SqQuartetProfile` objects:

.. code-block:: python

   import physquirrel as psq

   # A split quarnet: AB|CD
   q = psq.Quartet(psq.Split({'A', 'B'}, {'C', 'D'}))
   profile = psq.SqQuartetProfile([q])

   # A cycle quarnet: two resolved quartets sharing four taxa
   q1 = psq.Quartet(psq.Split({'A', 'B'}, {'C', 'D'}))
   q2 = psq.Quartet(psq.Split({'A', 'D'}, {'B', 'C'}))
   cycle_profile = psq.SqQuartetProfile([q1, q2], reticulation_leaf='A')

   profileset = psq.SqQuartetProfileSet(profiles=[profile, cycle_profile])

Each :class:`~physquirrel.datatypes.sqprofile.SqQuartetProfile` holds exactly **1 or 2** resolved
quartets on the same four taxa. A 2-quartet profile encodes a 4-cycle; an
optional ``reticulation_leaf`` identifies which taxon is below the hybrid node.
Quarnets can also be weighted by passing a ``profile_weight`` to the constructor; 
this is used in cycle resolution and scoring to downweight low-confidence profiles.

From a distance matrix
~~~~~~~~~~~~~~~~~~~~~~~

The delta heuristic infers a profile set from pairwise distances, represented
by a :class:`~phylozoo.core.distance.base.DistanceMatrix`. For each 4-taxon subset it computes
a score :math:`\delta \in [0,1]`; below the threshold ``lam`` the quartet is
classified as a split, otherwise as a cycle.

.. code-block:: python

   import numpy as np
   import physquirrel as psq

   matrix = np.array([
       [0.0, 0.1, 0.9, 0.9],
       [0.1, 0.0, 0.9, 0.9],
       [0.9, 0.9, 0.0, 0.1],
       [0.9, 0.9, 0.1, 0.0],
   ])
   dm = psq.DistanceMatrix(matrix, labels=['A', 'B', 'C', 'D'])

   profileset = psq.delta_heuristic(dm, lam=0.3, weight=True)

The ``weight`` parameter controls whether profile confidence (distance to the
threshold) is stored as profile-level weights in the returned
:class:`~physquirrel.datatypes.sqprofileset.SqQuartetProfileSet`.

From an MSA
~~~~~~~~~~~

When starting from aligned sequences, these first need to be converted into a distance matrix.
This can be done through phylozoo's :func:`~phylozoo.core.sequence.distances.hamming_distances` function,
which computes pairwise Hamming distances between sequences in the MSA. The resulting distance matrix
can then be passed to the delta heuristic as described above.
Alternatively, physquirrel provides a convenience function that handles this internally:
it computes Hamming distances and applies the delta heuristic in one call.

Use the ``lam`` parameter to control the split/cycle classification threshold,
and the optional ``weight`` parameter to store profile confidence as weights in
the returned :class:`~physquirrel.datatypes.sqprofileset.SqQuartetProfileSet`.

MSAs are represented by :class:`~phylozoo.core.sequence.base.MSA`. Both FASTA and NEXUS formats
are supported. Load from a file with :meth:`~phylozoo.core.sequence.base.MSA.load`, or construct
directly from a dictionary of sequences:

.. code-block:: python

   import physquirrel as psq

   # From a FASTA file (.fasta, .fa, .fas)
   msa = psq.MSA.load("alignment.fasta")

   # From a NEXUS file (.nexus, .nex, .nxs) — specify format explicitly
   msa = psq.MSA.load("alignment.nex", format="nexus")

   # From a dictionary of sequences
   msa = psq.MSA({"A": "ACGT", "B": "ACGT", "C": "TGCA", "D": "TGCA"})

   profileset = psq.delta_heuristic_from_msa(msa, lam=0.3)

Running Squirrel
----------------

The main entry point is the function :func:`~physquirrel.algorithms.squirrel.squirrel`. It accepts
an :class:`~physquirrel.datatypes.sqprofileset.SqQuartetProfileSet` (see previous section)and returns a
:class:`~phylozoo.core.network.sdnetwork.sd_phynetwork.SemiDirectedPhyNetwork`:

.. code-block:: python

   import physquirrel as psq

   network = psq.squirrel(profileset)

   print(f"eNewick string:  {network.to_string()}")

The pipeline runs internally:

1. **T\* tree** (:func:`~physquirrel.algorithms.tstar_tree.tstar_tree`) — B\*-compatible split system → initial tree.
2. **Adapted quartet joining** (:func:`~physquirrel.algorithms.qjoining.adapted_quartet_joining`) — refines the tree using quartet support.
3. **Unresolve tree** (:func:`~physquirrel.algorithms.unresolve_tree.unresolve_tree`) — iteratively contracts low-support splits, yielding a
   sequence of contracted trees from fully resolved to star.
4. **Cycle resolution** (:func:`~physquirrel.algorithms.cycle_resolution.resolve_cycles`) — inserts reticulation cycles at internal vertices of
   each contracted tree using a TSP-based circular ordering.
5. **Scoring** (:func:`~physquirrel.algorithms.qsimilarity.sqprofileset_similarity`) — each candidate network is scored by profile similarity
   (consistency with the input profile set); the best is returned.

Rooting at an outgroup
~~~~~~~~~~~~~~~~~~~~~~

Pass ``outgroup`` to return a :class:`~phylozoo.core.network.dnetwork.base.DirectedPhyNetwork` rooted
at the edge incident to the specified leaf:

.. code-block:: python

   import physquirrel as psq

   rooted = psq.squirrel(profileset, outgroup='A')
   assert isinstance(rooted, psq.DirectedPhyNetwork)
   assert rooted.root_node is not None

Convenience shortcuts
~~~~~~~~~~~~~~~~~~~~~

For common workflows, physquirrel provides one-call convenience functions that
compose :func:`~physquirrel.algorithms.delta_heuristic.delta_heuristic` and :func:`~physquirrel.algorithms.squirrel.squirrel`
internally.

**From a distance matrix**:

.. code-block:: python

   import physquirrel as psq

   network = psq.squirrel_from_distances(dm, lam=0.3, outgroup='A')

**From an MSA**:

.. code-block:: python

   import physquirrel as psq

   network = psq.squirrel_from_msa(msa, lam=0.3)

All keyword arguments accepted by :func:`~physquirrel.algorithms.squirrel.squirrel` (e.g.
``rho``, ``tsp_threshold``) can be forwarded through both convenience functions.


Parallelization
~~~~~~~~~~~~~~~~

The :func:`~physquirrel.algorithms.squirrel.squirrel` algorithm can be parallelized by 
handling the candidate trees in parallel during cycle resolution and scoring. This is controlled by the optional ``parallel`` argument, which accepts a `ParallelConfig` object from phylozoo. By default, the algorithm runs sequentially.
   
**Threading** (shared memory, lower overhead, good for 2–4 cores):

.. code-block:: python

   from phylozoo.utils.parallel import ParallelConfig, ParallelBackend
   import physquirrel as psq

   # Use threading with 4 worker threads
   parallel_config = ParallelConfig(backend=ParallelBackend.THREADING, n_jobs=4)
   network = psq.squirrel(profileset, parallel=parallel_config)

**Multiprocessing** (true parallelism, higher overhead, good for 4+ cores):

.. code-block:: python

   from phylozoo.utils.parallel import ParallelConfig, ParallelBackend
   import physquirrel as psq

   # Use multiprocessing with 8 worker processes
   parallel_config = ParallelConfig(backend=ParallelBackend.MULTIPROCESSING, n_jobs=8)
   network = psq.squirrel(profileset, parallel=parallel_config)

**Sequential** (default, no parallelization):

.. code-block:: python

   import physquirrel as psq

   # Omit the parallel argument or pass None
   network = psq.squirrel(profileset)
   # or explicitly:
   network = psq.squirrel(profileset, parallel=None)

Notes:

* Threading is recommended for most cases due to lower overhead; it works well for 2–4 cores.
* Multiprocessing enables true parallelism but incurs higher overhead (serialization, process spawning);
  use it for large datasets or 4+ cores.
* The number of jobs (``n_jobs``) should typically match your CPU core count. Use ``n_jobs=-1`` to
  automatically detect the number of available cores.

I/O and plotting
-----------------

Loading and saving quartet profile sets
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:class:`~physquirrel.datatypes.sqprofileset.SqQuartetProfileSet` supports two file formats.
The class uses the :class:`~phylozoo.utils.io.IOMixin` interface from phylozoo, which infers the format from the file extension by default.
Alternatively, the format can be specified explicitly when loading or saving. The four main PhyloZoo functions 
for I/O are: `save`, `to_string`, `load`, and `from_string`.

**psq format** (``.psq``, default)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``psq`` format is a compact JSON file that stores all information
losslessly: the full taxon set, per-profile quartet weights, reticulation leaf
assignments, and profile-level confidence weights. It is the recommended format
for archiving and exchanging profile sets between physquirrel runs.

.. code-block:: python

   import physquirrel as psq

   # Infer a profile set from a distance matrix
   profileset = psq.delta_heuristic(dm, lam=0.3)

   # Save to file — format inferred from .psq extension
   profileset.save("profiles.psq")

   # Load back
   loaded = psq.SqQuartetProfileSet.load("profiles.psq")


A ``profiles.psq`` file looks like:

.. code-block:: json

   {
     "taxa": ["A", "B", "C", "D"],
     "profiles": [
       {
         "taxa": ["A", "B", "C", "D"],
         "quartets": [
           {"quartet": {"type": "resolved",
                        "split": {"set1": ["A", "B"], "set2": ["C", "D"]}},
            "weight": 1.0}
         ],
         "reticulation_leaf": null,
         "profile_weight": 0.85
       }
     ]
   }

**profile_list format** (``.txt``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``profile_list`` format is a plain-text file with one profile per line.
Each line begins with a type tag followed by four taxon names and an optional
weight:

.. code-block:: text

   # physquirrel profile list
   # Lines starting with '#' are comments and are ignored.
   SQ A B C D 0.85    # split quarnet: split AB|CD, confidence weight 0.85
   4C A B C E 0.60    # 4-cycle: circular order A,B,C,E; A is the reticulation leaf

- ``SQ a b c d [weight]`` — a split quarnet with split :math:`ab|cd`.
- ``4C a b c d [weight]`` — a 4-cycle with circular order :math:`a,b,c,d`
  where :math:`a` is the reticulation leaf (below the hybrid node).

The weight field is optional and defaults to 1.0. This format is compatible
with the quarnet file format used in earlier versions of physquirrel.

.. code-block:: python

   import physquirrel as psq

   # Save as profile_list
   profileset.save("profiles.txt", format="profile_list")

   # Or get the string directly
   txt = profileset.to_string(format="profile_list")
   print(txt)

   # Load back — format must be specified explicitly for .txt
   loaded = psq.SqQuartetProfileSet.load("profiles.txt", format="profile_list")
   loaded2 = psq.from_profile_list(txt)

Saving networks
~~~~~~~~~~~~~~~

:class:`~phylozoo.core.network.sdnetwork.sd_phynetwork.SemiDirectedPhyNetwork` and
:class:`~phylozoo.core.network.dnetwork.base.DirectedPhyNetwork` are provided by phylozoo and support
the **eNewick** format (default) via the standard IOMixin interface:

.. code-block:: python

   import physquirrel as psq

   # Serialize to an eNewick string
   enewick_str = network.to_string()              # eNewick by default

   # Write directly to file (.nwk, .newick, .enewick, .enw)
   network.save("network.nwk")

See the `PhyloZoo documentation <https://nholtgrefe.github.io/phylozoo/>`_ for
the full list of supported formats and options.

Visualization
~~~~~~~~~~~~~~

physquirrel networks can be plotted directly using phylozoo's unified
:func:`~phylozoo.viz.plot` function, which dispatches by network type and
uses :mod:`matplotlib`:

.. code-block:: python

   from phylozoo.viz import plot

   plot(network)                          # semi-directed network, twopi layout
   plot(network, layout="pz-radial")      # radial layout (trees only)
   plot(rooted_network)                   # :class:`~phylozoo.core.network.dnetwork.base.DirectedPhyNetwork`, DAG layout

The ``viz`` extra dependency must be installed (``pip install physquirrel[viz]``
or ``pip install physquirrel[graphviz]``); see the :doc:`Installation Guide <installation>`
for details.

Quick Links
-----------

* :doc:`Quickstart <quickstart>` — a quick walkthrough of the main workflow
* :doc:`Installation Guide <installation>` — requirements and install options
* :doc:`API Reference <api/index>` — complete index of all functions and classes
* `PhyloZoo Documentation <https://nholtgrefe.github.io/phylozoo/>`_ — core datatypes and utilities for working with networks, splits, quartets, distance matrices, and MSAs.
