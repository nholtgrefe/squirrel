from __future__ import annotations

import itertools

from phylozoo import DistanceMatrix
from phylozoo.core.quartet.base import Quartet
from phylozoo.core.split.base import Split
from phylozoo.core.primitives.circular_ordering import CircularOrdering
from phylozoo.utils.exceptions import PhyloZooValueError

from ..datatypes.sqprofile import SqQuartetProfile
from ..datatypes.sqprofileset import SqQuartetProfileSet


def delta_heuristic(
    distance_matrix: DistanceMatrix,
    lam: float = 0.3,
    weight: bool = True,
) -> SqQuartetProfileSet:
    """
    Infer quartet profiles from a distance matrix using the delta heuristic.

    For each set of 4 taxa, computes delta = (d2 - d1) / (d2 - d0) where d0 <= d1 <= d2
    are the three possible split distances. If delta < lam the quartet is a split;
    otherwise a 4-cycle. Reticulation leaves are assigned globally by highest delta sum.

    Parameters
    ----------
    distance_matrix : DistanceMatrix
        Pairwise distances between taxa.
    lam : float
        Threshold in [0, 1]. Below → split, above → cycle. Default 0.3.
    weight : bool
        Whether to weight profiles by delta confidence. Default True.

    Returns
    -------
    SqQuartetProfileSet
    """
    if lam < 0 or lam > 1:
        raise PhyloZooValueError(f"Lambda must be in [0, 1], got {lam}")

    taxa_list = list(distance_matrix.labels)
    n = len(taxa_list)

    if n < 4:
        return SqQuartetProfileSet()

    taxon_to_idx = {taxon: idx for idx, taxon in enumerate(taxa_list)}
    dist_array = distance_matrix.np_array
    delta_sum: dict[str, float] = {taxon: 0.0 for taxon in taxa_list}
    profiles: list[tuple[SqQuartetProfile, float]] = []

    for i, j, k, l in itertools.combinations(range(n), 4):
        a, b, c, d = taxa_list[i], taxa_list[j], taxa_list[k], taxa_list[l]

        split1_dist = dist_array[i, j] + dist_array[k, l]
        split2_dist = dist_array[i, k] + dist_array[j, l]
        split3_dist = dist_array[i, l] + dist_array[j, k]

        dists = [split1_dist, split2_dist, split3_dist]
        split_configs = [({a, b}, {c, d}), ({a, c}, {b, d}), ({a, d}, {b, c})]
        sorted_indices = sorted(range(3), key=lambda idx: dists[idx])

        d0, d1, d2 = dists[sorted_indices[0]], dists[sorted_indices[1]], dists[sorted_indices[2]]

        if d0 == d1 == d2:
            delta = 0.0
        else:
            denominator = d2 - d0
            delta = 0.0 if denominator == 0 else (d2 - d1) / denominator

        delta_sum[a] += delta
        delta_sum[b] += delta
        delta_sum[c] += delta
        delta_sum[d] += delta

        if delta < lam:
            best_idx = sorted_indices[0]
            best_set1, best_set2 = split_configs[best_idx]
            quartet = Quartet(Split(best_set1, best_set2))
            profile = SqQuartetProfile([quartet])
            w = (max(abs(lam - delta) / lam, 1e-9) if lam > 0 else 1.0) if weight else 1.0
        else:
            wrong_idx = sorted_indices[2]
            wrong_set1, wrong_set2 = split_configs[wrong_idx]

            if wrong_set1 == {a, b} and wrong_set2 == {c, d}:
                circ = CircularOrdering([a, c, b, d])
            elif wrong_set1 == {a, c} and wrong_set2 == {b, d}:
                circ = CircularOrdering([a, b, c, d])
            else:
                circ = CircularOrdering([a, b, d, c])

            order_list = list(circ.order)
            a_ord, b_ord, c_ord, d_ord = order_list
            q1 = Quartet(Split({a_ord, b_ord}, {c_ord, d_ord}))
            q2 = Quartet(Split({a_ord, d_ord}, {b_ord, c_ord}))
            # Equal within-profile weights (sum to 1.0); confidence goes to profile weight
            profile = SqQuartetProfile([q1, q2])
            w = (max(abs(lam - delta) / (1.0 - lam), 1e-9) if lam < 1.0 else 1.0) if weight else 1.0

        profiles.append((profile, w))

    reticulation_order = [t for t, _ in sorted(delta_sum.items(), key=lambda x: x[1], reverse=True)]

    final_profiles: list[tuple[SqQuartetProfile, float]] = []
    for profile, w in profiles:
        if len(profile.quartets) == 2:
            ret_leaf = next((t for t in reticulation_order if t in profile.taxa), None)
            if ret_leaf is not None:
                final_profiles.append(
                    (SqQuartetProfile(dict(profile.quartets), reticulation_leaf=ret_leaf), w)
                )
            else:
                final_profiles.append((profile, w))
        else:
            final_profiles.append((profile, w))

    return SqQuartetProfileSet(profiles=final_profiles)
