"""
Quartet distance module.

Provides functions for computing distance matrices from quartet profiles
using quartet distance metrics.
"""

from __future__ import annotations

import itertools
from typing import TYPE_CHECKING, Literal

import numpy as np

from phylozoo.core.distance import DistanceMatrix
from phylozoo.utils.exceptions import PhyloZooValueError

if TYPE_CHECKING:
    from phylozoo.core.quartet.qprofileset import QuartetProfileSet
    from phylozoo.core.primitives.partition import Partition


def _topo_key(
    profile,
    leaf_to_set: 'dict[str, frozenset]',
    X_indices: 'dict[frozenset, int]',
    four_taxa_set: 'frozenset[str]',
) -> tuple:
    """Canonical topology key for a profile within a 4-subpartition.

    Returns ('s', frozenset) for split topologies and ('c', tuple) for cycle
    topologies. Used by 'best' mode to group profiles by plurality topology.

    For a split profile with split set1|set2 over partition sets at global
    indices {i,j,k,l}: key = ('s', frozenset({frozenset({i,j}), frozenset({k,l})}))
    where {i,j} are the indices of the two sets on the same side.

    For a cycle profile with circular ordering inducing adjacency at global
    indices: key = ('c', canonicalized 4-tuple of indices), with canonicalization
    starting at the minimum index and choosing the lexicographically smaller
    of the forward and backward traversal.
    """
    quartets = profile.quartets
    if len(quartets) == 1:
        set1 = next(iter(quartets)).split.set1
        side_a = frozenset(X_indices[leaf_to_set[leaf]] for leaf in four_taxa_set if leaf in set1)
        side_b = frozenset(X_indices[leaf_to_set[leaf]] for leaf in four_taxa_set if leaf not in set1)
        return ('s', frozenset({side_a, side_b}))
    else:
        ordering = next(iter(profile.circular_orderings))
        idx_tup = tuple(X_indices[leaf_to_set[leaf]] for leaf in ordering.order)
        start = idx_tup.index(min(idx_tup))
        fwd = tuple(idx_tup[(start + k) % 4] for k in range(4))
        bwd = tuple(idx_tup[(start - k) % 4] for k in range(4))
        return ('c', min(fwd, bwd))


def quartet_distance_with_partition(
    profileset: 'QuartetProfileSet',
    partition: 'Partition',
    rho: tuple[float, float, float, float] = (0.5, 1.0, 0.5, 1.0),
    weighted_distance: bool = True,
    representative_mode: Literal['average', 'best'] = 'best',
) -> DistanceMatrix:
    """
    Compute a distance matrix between partition sets based on quartet profiles.

    For each 4-subpartition of the partition, considers all quartets with one
    leaf from each set and aggregates their rho-distance contributions. The
    result is a k×k distance matrix over the k partition sets (TSP input for
    cycle resolution).

    Parameters
    ----------
    profileset : QuartetProfileSet
        Must be dense (a profile for every 4-taxon combination).
    partition : Partition
        A partition of the taxa. Distance is computed between the sets.
    rho : tuple[float, float, float, float]
        (rho_c, rho_s, rho_a, rho_o) — distance contributions for different
        quartet topologies. Must satisfy rho_c ≤ rho_s and rho_a ≤ rho_o.
        Default: (0.5, 1.0, 0.5, 1.0).

        For a **split profile** (1 quartet), the distance between two leaves is:

        - ``rho_c`` — same side of the split
        - ``rho_s`` — opposite sides of the split

        For a **cycle profile** (2 quartets), the distance between two leaves is:

        - ``rho_a`` — adjacent in the circular ordering
        - ``rho_o`` — opposite in the circular ordering
    weighted_distance : bool
        If True, each profile's rho-distance contribution is scaled by its
        profile weight before averaging. High-confidence profiles then pull
        the distance more than low-confidence ones, mirroring the Kalmanson
        weighting used in physquirrel v1. Default: False.
    representative_mode : {'average', 'best'}
        How profiles within each 4-subpartition are aggregated:

        - ``'average'`` (default): average rho-distance contributions over
          **all** representative leaf-partitions (one leaf from each set),
          then shrink the per-pair mean toward 1.5 by ``1 - confidence``
          where ``confidence`` is the weight fraction on the most-supported
          delta value for that pair. Reflects the full empirical distribution
          of quartet signals while still down-weighting ambiguous pairs
          (Kalmanson-style shrink, but over all contributing quarnets rather
          than a single elected representative).

        - ``'best'``: first vote (weighted by profile weight) across all
          representative leaf-partitions to determine the **plurality topology**
          (most common split direction or cycle orientation). Then average
          contributions from only the profiles that match the winning topology.
          This mirrors the v1 behaviour of building a single representative
          quarnet per 4-tuple of partition sets via a plurality vote.

    Returns
    -------
    DistanceMatrix
        k×k symmetric matrix with zero diagonal. Labels are the partition sets
        (frozensets).

    Raises
    ------
    PhyloZooValueError
        If rho constraints are violated, partition elements don't match
        profileset taxa, or profileset is not dense.
    """
    if len(rho) != 4:
        raise PhyloZooValueError(f"Rho vector must have 4 elements, got {len(rho)}")
    rho_c, rho_s, rho_a, rho_o = rho
    if rho_a > rho_o or rho_c > rho_s:
        raise PhyloZooValueError(
            "Rho vector must satisfy: rho_a <= rho_o and rho_c <= rho_s"
        )

    if partition.elements != profileset.taxa:
        raise PhyloZooValueError(
            "Partition elements must match profile set taxa. "
            f"Partition has {partition.elements}, profile set has {profileset.taxa}"
        )

    if not profileset.is_dense:
        raise PhyloZooValueError(
            "Profile set must be dense (have a profile for every 4-taxon combination)"
        )

    if representative_mode not in ('average', 'best'):
        raise PhyloZooValueError(
            f"representative_mode must be 'average' or 'best', got {representative_mode!r}"
        )

    n = len(partition)
    set_order = sorted(partition.parts, key=lambda p: sorted(p))
    set_to_index: dict[frozenset, int] = {part: idx for idx, part in enumerate(set_order)}
    D = np.zeros((n, n), dtype=np.float64)

    # Pre-extract rho scalars so the inner loop doesn't unpack the tuple each iteration.
    delta_same  = 2.0 * rho_c   # split profile, same side
    delta_diff  = 2.0 * rho_s   # split profile, different sides
    delta_adj   = 2.0 * rho_a   # cycle profile, adjacent
    delta_opp   = 2.0 * rho_o   # cycle profile, opposite

    for four_sub_partition in partition.subpartitions(4):
        X_indices: dict[frozenset, int] = {
            X: set_to_index[X] for X in four_sub_partition.parts
        }
        leaf_to_set: dict[str, frozenset] = {
            leaf: part
            for part in four_sub_partition.parts
            for leaf in part
        }
        pair_keys = list(itertools.combinations(sorted(X_indices.values()), 2))

        if representative_mode == 'best':
            _accumulate_best(
                four_sub_partition, profileset, leaf_to_set, X_indices,
                pair_keys, D, weighted_distance,
                delta_same, delta_diff, delta_adj, delta_opp,
            )
        else:
            _accumulate_average(
                four_sub_partition, profileset, leaf_to_set, X_indices,
                pair_keys, D, weighted_distance,
                delta_same, delta_diff, delta_adj, delta_opp,
            )

    constant = 2 * n - 4
    D += constant
    np.fill_diagonal(D, 0)

    return DistanceMatrix(D, labels=set_order)


def _accumulate_average(
    four_sub_partition,
    profileset,
    leaf_to_set: dict,
    X_indices: dict,
    pair_keys: list,
    D: np.ndarray,
    weighted_distance: bool,
    delta_same: float,
    delta_diff: float,
    delta_adj: float,
    delta_opp: float,
) -> None:
    """Average mode: accumulate rho-distance contributions from all representative leaf-partitions.

    Unlike the unshrunk average, every contributing quarnet adds its (weighted)
    delta to a per-pair mass distribution. After averaging, the pair's mean delta
    is shrunk toward 1.5 by ``1 - confidence``, where ``confidence`` is the
    fraction of contributing weight sitting on the single most-supported delta
    value for that pair. This mirrors the Kalmanson confidence shrinkage used by
    'best' mode (``scaled = avg - (avg - 1.5) * (1 - winner_frac)``), but the
    confidence is computed over *all* contributing quarnets — each weighted
    individually — instead of a single elected representative quarnet.
    """
    d_wsum: dict[tuple[int, int], float] = {k: 0.0 for k in pair_keys}
    d_wtot: dict[tuple[int, int], float] = {k: 0.0 for k in pair_keys}
    # Per-pair weight mass per distinct delta value, for the confidence shrink.
    d_mass: dict[tuple[int, int], dict[float, float]] = {k: {} for k in pair_keys}

    for four_leaf_partition in four_sub_partition.representative_partitions():
        four_taxa_set = four_leaf_partition.elements
        profile = profileset.get_profile(four_taxa_set)
        if profile is None:
            raise PhyloZooValueError(
                f"No profile found for 4-taxon set {four_taxa_set}. "
                "Profile set must be dense."
            )

        quartets = profile.quartets
        num_quartets = len(quartets)
        if num_quartets == 0 or num_quartets > 2:
            raise PhyloZooValueError(
                f"Profile for {four_taxa_set} must have 1 or 2 quartets, "
                f"got {num_quartets}"
            )
        for q in quartets:
            if not q.is_resolved():
                raise PhyloZooValueError(
                    f"Profile for {four_taxa_set} contains an unresolved quartet."
                )

        w = profileset._profiles[four_taxa_set][1] if weighted_distance else 1.0

        if num_quartets == 1:
            set1 = next(iter(quartets)).split.set1
            for leaf1, leaf2 in itertools.combinations(four_taxa_set, 2):
                same = (leaf1 in set1) == (leaf2 in set1)
                delta = delta_same if same else delta_diff
                i = X_indices[leaf_to_set[leaf1]]
                j = X_indices[leaf_to_set[leaf2]]
                if i > j:
                    i, j = j, i
                d_wsum[(i, j)] += delta * w
                d_wtot[(i, j)] += w
                m = d_mass[(i, j)]
                m[delta] = m.get(delta, 0.0) + w
        else:
            ordering = next(iter(profile.circular_orderings))
            order = list(ordering.order)
            adj: dict[str, frozenset] = {
                order[k]: frozenset({order[k - 1], order[(k + 1) % 4]})
                for k in range(4)
            }
            for leaf1, leaf2 in itertools.combinations(four_taxa_set, 2):
                delta = delta_adj if leaf2 in adj[leaf1] else delta_opp
                i = X_indices[leaf_to_set[leaf1]]
                j = X_indices[leaf_to_set[leaf2]]
                if i > j:
                    i, j = j, i
                d_wsum[(i, j)] += delta * w
                d_wtot[(i, j)] += w
                m = d_mass[(i, j)]
                m[delta] = m.get(delta, 0.0) + w

    for key in pair_keys:
        wtotal = d_wtot[key]
        if wtotal > 0:
            avg = d_wsum[key] / wtotal
            # Confidence = weight fraction on the single most-supported delta.
            confidence = max(d_mass[key].values()) / wtotal
            # Shrink toward 1.5 when confidence is low (mirrors 'best' mode).
            scaled = avg - (avg - 1.5) * (1.0 - confidence)
            i, j = key
            D[i, j] += scaled
            D[j, i] += scaled


def _tkey_canonical(k: tuple) -> tuple:
    """Deterministic sort key for a topology key, used to break ties in max()."""
    if k[0] == 's':
        return (k[0], tuple(sorted(tuple(sorted(s)) for s in k[1])))
    return k  # ('c', (i,j,k,l)) — already canonical


def _accumulate_best(
    four_sub_partition,
    profileset,
    leaf_to_set: dict,
    X_indices: dict,
    pair_keys: list,
    D: np.ndarray,
    weighted_distance: bool,
    delta_same: float,
    delta_diff: float,
    delta_adj: float,
    delta_opp: float,
) -> None:
    """Best mode: vote for the plurality topology per 4-subpartition, then use only its contributions."""
    # topo_votes: total weight voting for each topology key
    topo_votes: dict[tuple, float] = {}
    # Per-topology weighted sums and totals, keyed by (topology_key, pair_key)
    topo_wsum: dict[tuple, dict[tuple[int, int], float]] = {}
    topo_wtot: dict[tuple, dict[tuple[int, int], float]] = {}

    for four_leaf_partition in four_sub_partition.representative_partitions():
        four_taxa_set = four_leaf_partition.elements
        profile = profileset.get_profile(four_taxa_set)
        if profile is None:
            raise PhyloZooValueError(
                f"No profile found for 4-taxon set {four_taxa_set}. "
                "Profile set must be dense."
            )

        quartets = profile.quartets
        num_quartets = len(quartets)
        if num_quartets == 0 or num_quartets > 2:
            raise PhyloZooValueError(
                f"Profile for {four_taxa_set} must have 1 or 2 quartets, "
                f"got {num_quartets}"
            )
        for q in quartets:
            if not q.is_resolved():
                raise PhyloZooValueError(
                    f"Profile for {four_taxa_set} contains an unresolved quartet."
                )

        w = profileset._profiles[four_taxa_set][1] if weighted_distance else 1.0
        tkey = _topo_key(profile, leaf_to_set, X_indices, four_taxa_set)
        topo_votes[tkey] = topo_votes.get(tkey, 0.0) + w

        if tkey not in topo_wsum:
            topo_wsum[tkey] = {k: 0.0 for k in pair_keys}
            topo_wtot[tkey] = {k: 0.0 for k in pair_keys}

        if num_quartets == 1:
            set1 = next(iter(quartets)).split.set1
            for leaf1, leaf2 in itertools.combinations(four_taxa_set, 2):
                same = (leaf1 in set1) == (leaf2 in set1)
                delta = delta_same if same else delta_diff
                i = X_indices[leaf_to_set[leaf1]]
                j = X_indices[leaf_to_set[leaf2]]
                if i > j:
                    i, j = j, i
                topo_wsum[tkey][(i, j)] += delta * w
                topo_wtot[tkey][(i, j)] += w
        else:
            ordering = next(iter(profile.circular_orderings))
            order = list(ordering.order)
            adj: dict[str, frozenset] = {
                order[k]: frozenset({order[k - 1], order[(k + 1) % 4]})
                for k in range(4)
            }
            for leaf1, leaf2 in itertools.combinations(four_taxa_set, 2):
                delta = delta_adj if leaf2 in adj[leaf1] else delta_opp
                i = X_indices[leaf_to_set[leaf1]]
                j = X_indices[leaf_to_set[leaf2]]
                if i > j:
                    i, j = j, i
                topo_wsum[tkey][(i, j)] += delta * w
                topo_wtot[tkey][(i, j)] += w

    if not topo_votes:
        return

    total_votes = sum(topo_votes.values())
    winner = max(topo_votes, key=lambda k: (topo_votes[k], _tkey_canonical(k)))
    winner_frac = topo_votes[winner] / total_votes if total_votes > 0 else 1.0

    for key in pair_keys:
        wtotal = topo_wtot[winner][key]
        if wtotal > 0:
            avg = topo_wsum[winner][key] / wtotal
            # Scale toward 1.5 when confidence is low — mirrors Kalmanson: delta = 1+(1-wf)/2 or 2-(1-wf)/2
            scaled = avg - (avg - 1.5) * (1.0 - winner_frac)
            i, j = key
            D[i, j] += scaled
            D[j, i] += scaled
