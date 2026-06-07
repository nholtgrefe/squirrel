from __future__ import annotations

from typing import TYPE_CHECKING

from phylozoo import DistanceMatrix, MSA
from phylozoo.core.sequence.distances import hamming_distances

from .algorithms.delta_heuristic import delta_heuristic
from .algorithms.squirrel import squirrel
from .datatypes.sqprofileset import SqQuartetProfileSet

if TYPE_CHECKING:
    from phylozoo import SemiDirectedPhyNetwork, DirectedPhyNetwork


def delta_heuristic_from_msa(
    msa: MSA,
    lam: float = 0.3,
    weight: bool = True,
) -> SqQuartetProfileSet:
    """
    Compute a SqQuartetProfileSet from an MSA using Hamming distances and the delta heuristic.

    Parameters
    ----------
    msa : MSA
    lam : float
        Delta threshold separating split quartets from cycle quartets. Default 0.3.
    weight : bool
        Weight profiles by delta confidence. Default True.

    Returns
    -------
    SqQuartetProfileSet
    """
    dm: DistanceMatrix = hamming_distances(msa)
    return delta_heuristic(dm, lam=lam, weight=weight)


def squirrel_from_distances(
    distance_matrix: DistanceMatrix,
    lam: float = 0.3,
    weight: bool = True,
    outgroup: str | None = None,
    **kwargs,
) -> 'SemiDirectedPhyNetwork | DirectedPhyNetwork':
    """
    Reconstruct a phylogenetic network directly from a distance matrix.

    Runs delta_heuristic then squirrel in one call.

    Parameters
    ----------
    distance_matrix : DistanceMatrix
    lam : float
        Delta threshold. Default 0.3.
    weight : bool
        Weight profiles by delta confidence. Default True.
    outgroup : str | None
        If provided, returns a rooted DirectedPhyNetwork. Default None.
    **kwargs
        Passed to squirrel (e.g. rho, tsp_threshold, parallel).

    Returns
    -------
    SemiDirectedPhyNetwork | DirectedPhyNetwork
    """
    profileset = delta_heuristic(distance_matrix, lam=lam, weight=weight)
    return squirrel(profileset, outgroup=outgroup, **kwargs)


def squirrel_from_msa(
    msa: MSA,
    lam: float = 0.3,
    weight: bool = True,
    outgroup: str | None = None,
    **kwargs,
) -> 'SemiDirectedPhyNetwork | DirectedPhyNetwork':
    """
    Reconstruct a phylogenetic network directly from an MSA.

    Computes Hamming distances, runs delta_heuristic, then squirrel.

    Parameters
    ----------
    msa : MSA
    lam : float
        Delta threshold. Default 0.3.
    weight : bool
        Weight profiles by delta confidence. Default True.
    outgroup : str | None
        If provided, returns a rooted DirectedPhyNetwork. Default None.
    **kwargs
        Passed to squirrel (e.g. rho, tsp_threshold, parallel).

    Returns
    -------
    SemiDirectedPhyNetwork | DirectedPhyNetwork
    """
    dm: DistanceMatrix = hamming_distances(msa)
    return squirrel_from_distances(dm, lam=lam, weight=weight, outgroup=outgroup, **kwargs)
