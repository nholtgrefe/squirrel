import itertools
from collections.abc import Iterator
from typing import Any

from phylozoo import SemiDirectedPhyNetwork, Split
from phylozoo.core.network.sdnetwork.classifications import is_tree
from phylozoo.core.network.sdnetwork.conversions import sdnetwork_from_graph
from phylozoo.core.network.sdnetwork.derivations import split_from_cutedge
from phylozoo.core.network.sdnetwork.features import cut_edges
from phylozoo.core.primitives.m_multigraph.transformations import identify_vertices
from phylozoo.core.quartet.qprofileset import QuartetProfileSet
from phylozoo.utils.exceptions import PhyloZooValueError


def split_support(profileset: QuartetProfileSet, split: Split) -> float:
    """
    Compute quartet support for a split A|B.

    Returns the weighted ratio of trivial, resolved profiles that agree with the
    split over all 4-taxon sets {a1, a2, b1, b2} with a1, a2 ∈ A, b1, b2 ∈ B.

    Parameters
    ----------
    profileset : QuartetProfileSet
    split : Split

    Returns
    -------
    float in [0, 1], or 0.0 if no matching quartets found.
    """
    if split.is_trivial:
        raise PhyloZooValueError("Split must be non-trivial")
    if split.elements != profileset.taxa:
        raise PhyloZooValueError("Split taxa must match profile set taxa")

    support = 0.0
    total_weight = 0.0

    for a1, a2 in itertools.combinations(split.set1, 2):
        for b1, b2 in itertools.combinations(split.set2, 2):
            four_taxa = frozenset({a1, a2, b1, b2})
            profile = profileset.get_profile(four_taxa)
            if profile is None or not profile.is_trivial():
                continue

            quartet = next(iter(profile.quartets))
            quartet_split = quartet.split
            if quartet_split is None:
                continue

            pw = profileset.get_profile_weight(four_taxa) or 1.0
            qw = profile.get_weight(quartet) or 1.0
            total_weight += pw * qw

            if quartet_split == Split({a1, a2}, {b1, b2}):
                support += pw * qw

    return support / total_weight if total_weight > 0 else 0.0


def unresolve_tree(
    tree: SemiDirectedPhyNetwork,
    profileset: QuartetProfileSet,
) -> Iterator[SemiDirectedPhyNetwork]:
    """
    Yield a sequence of trees by iteratively contracting least-supported splits.

    Computes split support for all non-trivial splits once, sorts them by support
    (lowest first), and contracts them one at a time. Yields after each contraction,
    starting with the original tree.

    Parameters
    ----------
    tree : SemiDirectedPhyNetwork
    profileset : QuartetProfileSet

    Yields
    ------
    SemiDirectedPhyNetwork
    """
    if not is_tree(tree):
        raise PhyloZooValueError("Tree must be a valid tree")
    if tree.taxa != profileset.taxa:
        raise PhyloZooValueError("Tree taxa must match profile set taxa")

    current_tree = tree
    cut_edges_set = cut_edges(current_tree)
    split_to_edge: dict[Split, tuple[Any, Any, int]] = {}

    for u, v, key in cut_edges_set:
        edge_split = split_from_cutedge(current_tree, u, v, key=key)
        split_to_edge[edge_split] = (u, v, key)

    non_trivial_splits = [s for s in split_to_edge if not s.is_trivial]

    if not non_trivial_splits:
        yield current_tree
        return

    split_supports = {s: split_support(profileset, s) for s in non_trivial_splits}
    sorted_splits = sorted(non_trivial_splits, key=lambda x: split_supports[x])

    yield current_tree

    for split in sorted_splits:
        u, v, key = split_to_edge[split]

        graph_copy = current_tree._graph.copy()
        identify_vertices(graph_copy, [u, v])
        current_tree = sdnetwork_from_graph(graph_copy, network_type='semi-directed')

        del split_to_edge[split]

        updated: dict[Split, tuple[Any, Any, int]] = {}
        for s, (eu, ev, ek) in split_to_edge.items():
            new_u = u if eu == v else eu
            new_v = u if ev == v else ev
            if new_u == new_v:
                continue
            updated[s] = (new_u, new_v, ek)
        split_to_edge = updated

        yield current_tree
