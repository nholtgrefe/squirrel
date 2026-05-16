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
    from phylozoo.core.quartet.qprofile import QuartetProfile
    from phylozoo.core.quartet.qprofileset import QuartetProfileSet
    from phylozoo.core.primitives.partition import Partition


def _rho_distance(
    profile: 'QuartetProfile',
    leaf1: str,
    leaf2: str,
    rho: tuple[float, float, float, float],
) -> float:
    """
    Compute rho-distance between two leaves in a quartet profile.

    If the profile has 1 quartet (split quarnet), uses split-based logic:
    - Same side of split: rho_c
    - Different sides: rho_s

    If the profile has 2 quartets (four-cycle), uses circular ordering logic:
    - Adjacent in ordering: rho_a
    - Opposite in ordering: rho_o

    Parameters
    ----------
    profile : QuartetProfile
    leaf1 : str
    leaf2 : str
    rho : tuple[float, float, float, float]
        (rho_c, rho_s, rho_a, rho_o)

    Returns
    -------
    float
    """
    rho_c, rho_s, rho_a, rho_o = rho
    num_quartets = len(profile.quartets)

    if num_quartets == 1:
        quartet = next(iter(profile.quartets))
        split = quartet.split
        if split is None:
            raise PhyloZooValueError("Quartet must be resolved (have a split)")
        same_side = (leaf1 in split.set1 and leaf2 in split.set1) or (
            leaf1 in split.set2 and leaf2 in split.set2
        )
        return float(rho_c) if same_side else float(rho_s)

    elif num_quartets == 2:
        circular_orderings = profile.circular_orderings
        if circular_orderings is None or len(circular_orderings) == 0:
            raise PhyloZooValueError(
                "Profile with 2 quartets must have a circular ordering"
            )
        ordering = next(iter(circular_orderings))
        return float(rho_a) if ordering.are_neighbors(leaf1, leaf2) else float(rho_o)

    else:
        raise PhyloZooValueError(
            f"Profile must have 1 or 2 quartets, got {num_quartets}"
        )


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

    for four_sub_partition in partition.subpartitions(4):
        X_indices: dict[frozenset, int] = {
            X: set_to_index[X] for X in four_sub_partition.parts
        }
        leaf_to_set: dict[str, frozenset] = {
            leaf: part
            for part in four_sub_partition.parts
            for leaf in part
        }
        # (i, j) -> [weighted_distance_distance_sum, total_weight]
        d_accum: dict[tuple[int, int], list[float]] = {
            (i, j): [0.0, 0.0]
            for i, j in itertools.combinations(X_indices.values(), 2)
        }

        for four_leaf_partition in four_sub_partition.representative_partitions():
            four_taxa_set = four_leaf_partition.elements
            profile = profileset.get_profile(four_taxa_set)
            if profile is None:
                raise PhyloZooValueError(
                    f"No profile found for 4-taxon set {four_taxa_set}. "
                    "Profile set must be dense."
                )

            num_quartets = len(profile.quartets)
            if num_quartets == 0 or num_quartets > 2:
                raise PhyloZooValueError(
                    f"Profile for {four_taxa_set} must have 1 or 2 quartets, "
                    f"got {num_quartets}"
                )
            for quartet in profile.quartets:
                if not quartet.is_resolved():
                    raise PhyloZooValueError(
                        f"Profile for {four_taxa_set} contains an unresolved quartet."
                    )

            # When weighted_distance, scale by profile confidence; otherwise treat all equally.
            w = profileset._profiles[four_taxa_set][1] if weighted_distance else 1.0

            for leaf1, leaf2 in itertools.combinations(four_taxa_set, 2):
                rho_dist = _rho_distance(profile, leaf1, leaf2, rho)
                Xi = leaf_to_set[leaf1]
                Xj = leaf_to_set[leaf2]
                i = X_indices[Xi]
                j = X_indices[Xj]
                if i > j:
                    i, j = j, i
                d_accum[(i, j)][0] += 2 * rho_dist * w
                d_accum[(i, j)][1] += w

        for (i, j), (wsum, wtotal) in d_accum.items():
            if wtotal > 0:
                avg = wsum / wtotal
                D[i, j] += avg
                D[j, i] += avg

    constant = 2 * n - 4
    D += constant
    np.fill_diagonal(D, 0)

    return DistanceMatrix(D, labels=set_order)
