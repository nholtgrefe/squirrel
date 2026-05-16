"""
Quartet distance module.

Provides functions for computing distance matrices from quartet profiles
using quartet distance metrics.
"""

from __future__ import annotations

import itertools
from typing import TYPE_CHECKING

import numpy as np

from phylozoo.core.distance import DistanceMatrix
from phylozoo.utils.exceptions import PhyloZooValueError

if TYPE_CHECKING:
    from phylozoo.core.quartet.qprofileset import QuartetProfileSet
    from phylozoo.core.primitives.partition import Partition


def quartet_distance_with_partition(
    profileset: 'QuartetProfileSet',
    partition: 'Partition',
    rho: tuple[float, float, float, float] = (0.5, 1.0, 0.5, 1.0),
    weighted_distance: bool = False,
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

    n = len(partition)
    set_order = list(partition.parts)
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
        # Two flat dicts instead of one dict-of-lists: one dict lookup per update.
        pair_keys = list(itertools.combinations(X_indices.values(), 2))
        d_wsum: dict[tuple[int, int], float] = {k: 0.0 for k in pair_keys}
        d_wtot: dict[tuple[int, int], float] = {k: 0.0 for k in pair_keys}

        for four_leaf_partition in four_sub_partition.representative_partitions():
            four_taxa_set = four_leaf_partition.elements
            profile = profileset.get_profile(four_taxa_set)
            if profile is None:
                raise PhyloZooValueError(
                    f"No profile found for 4-taxon set {four_taxa_set}. "
                    "Profile set must be dense."
                )

            # Cache the quartets dict once — avoids re-calling the property per pair.
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
                # Split profile: precompute set1 once; membership test is O(1).
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
            else:
                # Cycle profile: precompute adjacency as a dict of neighbour sets.
                # One circular ordering suffices (all are equivalent rotations/reflections).
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

        for key in pair_keys:
            wtotal = d_wtot[key]
            if wtotal > 0:
                avg = d_wsum[key] / wtotal
                i, j = key
                D[i, j] += avg
                D[j, i] += avg

    constant = 2 * n - 4
    D += constant
    np.fill_diagonal(D, 0)

    return DistanceMatrix(D, labels=set_order)
