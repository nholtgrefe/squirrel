Changelog
=========

Version 2.0.0
-------------

physquirrel 2.0 is a major release with a new architecture, new data types, and
several algorithmic improvements.

**Breaking changes**

- All core datatypes (phylogenetic networks, splits, quartets, distance matrices,
  MSAs) are now provided by `phylozoo <https://github.com/nholtgrefe/phylozoo>`_.
  The previous self-contained implementations have been removed.
- ``QuarnetSet`` class removed; replaced by :class:`~physquirrel.SqQuartetProfileSet`.
- Runtime dependencies changed: ``matplotlib`` and ``scipy`` removed;
  ``phylozoo`` (≥ 0.1.2) and ``numba`` (≥ 0.56) added.
- The ``[viz]`` optional dependency now relies on PhyloZoo's plotting utilities
  (``physquirrel[viz]``).

**New data types**

- :class:`~physquirrel.SqQuartetProfile` — extends PhyloZoo's
  :class:`~phylozoo.core.quartet.qprofile.QuartetProfile`; holds 1–2 resolved
  quartets and an optional ``reticulation_leaf`` for 4-cycle profiles.
- :class:`~physquirrel.SqQuartetProfileSet` — collection of
  :class:`~physquirrel.SqQuartetProfile` objects, one per 4-leaf subset, with
  optional per-profile weights.

**New I/O formats**

- ``.psq`` (JSON) — lossless, compact serialization via
  :func:`~physquirrel.to_psq` / :func:`~physquirrel.from_psq`.
- Profile list (plain text) — human-readable ``SQ`` / ``4C`` format via
  :func:`~physquirrel.to_profile_list` / :func:`~physquirrel.from_profile_list`.

**New algorithms**

- :func:`~physquirrel.quartet_distance_with_partition` — weighted quartet
  distance between a profile set and a partition, used internally in cycle
  resolution.
- :func:`~physquirrel.split_support` — refactored to correctly handle cycle
  profiles when scoring candidate splits.

**Algorithmic improvements**

- New ``representative_mode`` parameter in :func:`~physquirrel.squirrel`
  (``'best'`` or ``'average'``):

  - ``'best'`` (default) — elects a plurality quarnet per 4-tuple of partition
    sets, matching v1 behaviour.
  - ``'average'`` — averages rho-distances across all representative
    leaf-partitions for a fuller use of the empirical signal.

- Weighted profile support: profiles can carry per-profile confidence weights
  that propagate into cycle resolution and scoring.
- Parallel execution via :class:`~phylozoo.utils.parallel.ParallelConfig`;
  supports both threading and multiprocessing backends.
- JIT-accelerated Held-Karp TSP solver via numba (replaces scipy).

**Other**

- Full Sphinx documentation with API reference, manual, quickstart, and
  installation guide.
- Comprehensive test suite: 242 tests including integration tests on HIV-1,
  alignments.

----

Version 1.0.7
-------------

Last release of the v1.x series.
