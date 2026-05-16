import itertools
from typing import Any

from phylozoo import SemiDirectedPhyNetwork, Split
from phylozoo.core.network.sdnetwork.classifications import is_tree
from phylozoo.core.network.sdnetwork.conversions import sdnetwork_from_graph
from phylozoo.core.network.sdnetwork.derivations import split_from_cutedge
from phylozoo.core.quartet.base import Quartet
from phylozoo.core.quartet.qprofileset import QuartetProfileSet
from phylozoo.utils.exceptions import PhyloZooValueError


def _omega_bar(
    profileset: QuartetProfileSet,
    X1: frozenset[str],
    X2: frozenset[str],
    X3: frozenset[str],
    X4: frozenset[str],
) -> float:
    """
    Quartet support score for the split X1 ∪ X2 | X3 ∪ X4.

    Only considers trivial, resolved profiles.
    Returns 1.0 if no quartets are found.
    """
    res = 0.0
    total_sum = 0.0

    for (x1, x2, x3, x4) in itertools.product(X1, X2, X3, X4):
        four_taxa = frozenset({x1, x2, x3, x4})
        profile = profileset.get_profile(four_taxa)
        if profile is None or not profile.is_trivial():
            continue

        quartet = next(iter(profile.quartets))
        quartet_split = quartet.split
        if quartet_split is None:
            continue

        pw = profileset.get_profile_weight(four_taxa) or 1.0
        qw = profile.get_weight(quartet) or 1.0
        total_weight = pw * qw
        total_sum += total_weight

        if (quartet_split.set1 == {x1, x2} and quartet_split.set2 == {x3, x4}) or \
           (quartet_split.set2 == {x1, x2} and quartet_split.set1 == {x3, x4}):
            res += total_weight

    return res / total_sum if total_sum > 0 else 1.0


def adapted_quartet_joining(
    profileset: QuartetProfileSet,
    starting_tree: SemiDirectedPhyNetwork,
) -> SemiDirectedPhyNetwork:
    """
    Adapted quartet joining algorithm starting from an arbitrary tree.

    Iteratively refines the tree by joining the best-scoring neighbor pair at each
    high-degree internal node until the tree is binary.

    Parameters
    ----------
    profileset : QuartetProfileSet
    starting_tree : SemiDirectedPhyNetwork

    Returns
    -------
    SemiDirectedPhyNetwork
    """
    if not is_tree(starting_tree):
        raise PhyloZooValueError("Starting tree must be a valid tree")
    if starting_tree.taxa != profileset.taxa:
        raise PhyloZooValueError("Starting tree taxa must match profile set taxa")

    current_tree = starting_tree
    n = len(profileset.taxa)
    if n < 4:
        return current_tree

    total_edges = current_tree.number_of_edges()
    nr_iterations = n - 3 - (total_edges - n)

    for _ in range(nr_iterations):
        center_nodes = [
            node for node in current_tree.internal_nodes
            if current_tree.degree(node) > 3
        ]
        if not center_nodes:
            break

        score: dict[tuple[Any, Any, Any], float] = {}

        for center_node in center_nodes:
            if current_tree.degree(center_node) == 3:
                continue

            neighbors_list = list(current_tree.neighbors(center_node))
            C: dict[Any, frozenset[str]] = {}

            for u in neighbors_list:
                _, (_, taxa_u), _ = split_from_cutedge(
                    current_tree, u, center_node, return_node_taxa=True
                )
                C[u] = taxa_u

            for (u1, u2) in itertools.combinations(neighbors_list, 2):
                score[(u1, u2, center_node)] = 0.0
                for (A1, A2) in itertools.combinations(C.values(), 2):
                    if C[u1] not in [A1, A2] and C[u2] not in [A1, A2]:
                        score[(u1, u2, center_node)] += _omega_bar(
                            profileset, C[u1], C[u2], A1, A2,
                        )

        if not score:
            break

        u1star, u2star, cstar = max(score, key=score.get)
        tree_graph = current_tree._graph.copy()

        if tree_graph.has_edge(cstar, u1star):
            tree_graph.remove_edge(cstar, u1star)
        if tree_graph.has_edge(cstar, u2star):
            tree_graph.remove_edge(cstar, u2star)

        w = next(tree_graph.generate_node_ids(1))
        tree_graph.add_node(w)
        tree_graph.add_undirected_edge(w, u1star)
        tree_graph.add_undirected_edge(w, u2star)
        tree_graph.add_undirected_edge(w, cstar)

        current_tree = sdnetwork_from_graph(tree_graph, network_type='semi-directed')

    return current_tree


def quartet_joining(profileset: QuartetProfileSet) -> SemiDirectedPhyNetwork:
    """
    Quartet joining algorithm starting from a star tree.

    Constructs an initial star tree (all taxa connected to a single center node)
    and refines it using adapted_quartet_joining.

    Parameters
    ----------
    profileset : QuartetProfileSet

    Returns
    -------
    SemiDirectedPhyNetwork
    """
    if len(profileset.taxa) < 4:
        raise PhyloZooValueError("Quartet joining requires at least 4 taxa")

    center = "_center"
    taxa_list = sorted(profileset.taxa)
    undirected_edges = [(center, taxon) for taxon in taxa_list]
    nodes = [(taxon, {'label': taxon}) for taxon in taxa_list]

    starting_tree = SemiDirectedPhyNetwork(undirected_edges=undirected_edges, nodes=nodes)
    return adapted_quartet_joining(profileset, starting_tree)
