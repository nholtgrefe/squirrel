from typing import Any, Literal, TYPE_CHECKING

from phylozoo import SemiDirectedPhyNetwork
from .tsp import approximate_tsp_tour, optimal_tsp_tour
from phylozoo.core.primitives.circular_ordering import CircularSetOrdering
from phylozoo.core.primitives.partition import Partition
from phylozoo.core.primitives.m_multigraph.features import source_components
from .qdistance import quartet_distance_with_partition
from phylozoo.core.network.sdnetwork.conversions import sdnetwork_from_graph
from phylozoo.core.network.sdnetwork.derivations import partition_from_blob
from phylozoo.core.network.sdnetwork.features import cut_vertices
from phylozoo.utils.exceptions import PhyloZooValueError, PhyloZooAlgorithmError

if TYPE_CHECKING:
    from phylozoo.core.network.sdnetwork import MixedPhyNetwork
    from phylozoo.core.primitives.m_multigraph import MixedMultiGraph
    from ..datatypes.sqprofileset import SqQuartetProfileSet
    from phylozoo.core.quartet.qprofileset import QuartetProfileSet


def _qprofiles_to_circular_ordering(
    profileset: 'QuartetProfileSet',
    partition: Partition,
    rho: tuple[float, float, float, float] = (0.5, 1.0, 0.5, 1.0),
    tsp_method: Literal['optimal', 'simulated_annealing', 'greedy', 'christofides'] = 'optimal',
    weighted_distance: bool = False,
) -> CircularSetOrdering:
    """
    Compute the optimal circular set ordering from quartet profiles via TSP.

    Parameters
    ----------
    profileset : QuartetProfileSet
    partition : Partition
    rho : tuple of 4 floats
        Distance contributions (rho_c, rho_s, rho_a, rho_o). Default: (0.5, 1.0, 0.5, 1.0).
    tsp_method : str
        One of 'optimal', 'simulated_annealing', 'greedy', 'christofides'.
    weighted_distance : bool
        If True, profile weights are used to scale distance contributions.

    Returns
    -------
    CircularSetOrdering
    """
    if len(partition) < 3:
        raise PhyloZooValueError(
            f"Partition must have at least 3 sets for TSP, got {len(partition)}"
        )

    dist_matrix = quartet_distance_with_partition(
        profileset=profileset, partition=partition, rho=rho, weighted_distance=weighted_distance
    )

    if tsp_method == 'optimal':
        tour = optimal_tsp_tour(dist_matrix)
    elif tsp_method in ('simulated_annealing', 'greedy', 'christofides'):
        tour = approximate_tsp_tour(dist_matrix, method=tsp_method)
    else:
        raise PhyloZooValueError(
            f"Invalid tsp_method: {tsp_method}. "
            "Must be one of ['optimal', 'simulated_annealing', 'greedy', 'christofides']"
        )

    set_order = [set(fs) for fs in tour.order]
    return CircularSetOrdering(set_order)


def _qprofiles_to_hybrid_ranking(
    profileset: 'SqQuartetProfileSet',
    partition: Partition,
    weights: bool = True,
) -> list[frozenset[str]]:
    """
    Rank partition sets by likelihood of being the hybrid (reticulation) set.

    For 4-set partitions, uses reticulation_leaf from 4-cycle profiles.
    For larger partitions, aggregates cycle percentages from 4-subpartitions.

    Parameters
    ----------
    profileset : SqQuartetProfileSet
    partition : Partition
    weights : bool
        Whether to use profile weights in voting. Default True.

    Returns
    -------
    list[frozenset[str]]
        Partition sets ordered from most to least likely to be hybrid.
    """
    set_voting: dict[frozenset[str], float] = {frozenset(part): 0.0 for part in partition.parts}

    if len(partition) == 4:
        for four_leaf_partition in partition.representative_partitions():
            four_taxa_set = four_leaf_partition.elements
            profile = profileset.get_profile(four_taxa_set)
            if profile is None:
                continue
            if len(profile.quartets) == 2 and profile.reticulation_leaf is not None:
                ret_leaf = profile.reticulation_leaf
                for part in partition.parts:
                    if ret_leaf in part:
                        X = frozenset(part)
                        if weights:
                            profile_weight = profileset._profiles[four_taxa_set][1]
                            set_voting[X] += profile_weight
                        else:
                            set_voting[X] += 1.0
                        break
    else:
        for four_sub_partition in partition.subpartitions(4):
            splits_count = 0.0
            cycles_count = 0.0
            for four_leaf_partition in four_sub_partition.representative_partitions():
                four_taxa_set = four_leaf_partition.elements
                profile = profileset.get_profile(four_taxa_set)
                if profile is None:
                    continue
                if weights:
                    pw = profileset._profiles[four_taxa_set][1]
                    if len(profile.quartets) == 1:
                        splits_count += pw
                    elif len(profile.quartets) == 2:
                        cycles_count += pw
                else:
                    if len(profile.quartets) == 1:
                        splits_count += 1.0
                    elif len(profile.quartets) == 2:
                        cycles_count += 1.0

            total = splits_count + cycles_count
            if total > 0:
                cycle_perc = cycles_count / total
                for part in four_sub_partition.parts:
                    set_voting[frozenset(part)] += cycle_perc

    return sorted(set_voting, key=set_voting.get, reverse=True)


def _insert_cycle(
    network: SemiDirectedPhyNetwork,
    vertex: Any,
    circular_setorder: CircularSetOrdering,
    reticulation_ranking: list[frozenset[str]] | None = None,
    outgroup: str | None = None,
) -> 'SemiDirectedPhyNetwork | MixedPhyNetwork':
    """
    Replace a cut-vertex with an undirected cycle, then orient one hybrid node.

    Parameters
    ----------
    network : SemiDirectedPhyNetwork
    vertex : Any
        Must be a cut-vertex.
    circular_setorder : CircularSetOrdering
        Must match the partition induced by vertex.
    reticulation_ranking : list[frozenset[str]] | None
        Ordered list of sets to try as hybrid. If None, returns MixedPhyNetwork.
    outgroup : str | None

    Returns
    -------
    SemiDirectedPhyNetwork or MixedPhyNetwork
    """
    cut_verts = cut_vertices(network)
    if vertex not in cut_verts:
        raise PhyloZooValueError(f"Vertex {vertex} is not a cut-vertex in the network")

    induced_partition, edge_taxa_list = partition_from_blob(
        network, {vertex}, return_edge_taxa=True
    )

    if len(circular_setorder.parts) != len(induced_partition.parts):
        raise PhyloZooValueError(
            f"Circular set ordering ({len(circular_setorder.parts)} sets) does not match "
            f"partition induced by vertex ({len(induced_partition.parts)} sets)."
        )
    for s in circular_setorder.parts:
        if s not in induced_partition.parts:
            raise PhyloZooValueError(f"Set {s} in ordering is not in partition.")
    for s in induced_partition.parts:
        if s not in circular_setorder.parts:
            raise PhyloZooValueError(f"Set {s} in partition is not in ordering.")

    if reticulation_ranking is not None:
        for ret_set in reticulation_ranking:
            ret_frozen = frozenset(ret_set) if not isinstance(ret_set, frozenset) else ret_set
            if ret_frozen not in circular_setorder:
                raise PhyloZooValueError(
                    f"Reticulation set {ret_set} is not part of the circular set ordering"
                )

    graph_copy = network._graph.copy()

    def replace_vertex_with_cycle(
        graph: 'MixedMultiGraph',
        vertex: Any,
        circular_setorder: CircularSetOrdering,
        edge_taxa_list: list[tuple[Any, Any, frozenset[str]]],
        original_graph: 'MixedMultiGraph',
    ) -> tuple['MixedMultiGraph', list[Any], dict[frozenset[str], Any]]:
        graph.remove_node(vertex)
        n_cycle = len(circular_setorder)
        cycle_nodes = list(graph.generate_node_ids(n_cycle))

        for i in range(n_cycle):
            u = cycle_nodes[i]
            v = cycle_nodes[(i + 1) % n_cycle]
            graph.add_undirected_edge(u, v)

        set_to_node: dict[frozenset[str], Any] = {}
        for i, part_set in enumerate(circular_setorder.setorder):
            set_to_node[frozenset(part_set)] = cycle_nodes[i]

        for u, v, taxa_set in edge_taxa_list:
            cycle_node = set_to_node.get(taxa_set)
            if cycle_node is None:
                raise PhyloZooAlgorithmError(
                    f"Could not find cycle node for taxa set {taxa_set}"
                )
            attrs: dict[str, Any] = {}
            if original_graph.has_edge(u, vertex):
                for key in original_graph[u][vertex]:
                    attrs = original_graph[u][vertex][key].copy()
                    break
                attrs.pop('gamma', None)
            elif original_graph.has_edge(vertex, u):
                for key in original_graph[vertex][u]:
                    attrs = original_graph[vertex][u][key].copy()
                    break
                attrs.pop('gamma', None)
            graph.add_undirected_edge(u, cycle_node, **attrs)

        return graph, cycle_nodes, set_to_node

    def try_hybrid(
        graph: 'MixedMultiGraph',
        ret_set: frozenset[str],
        circular_setorder: CircularSetOrdering,
        cycle_nodes: list[Any],
        outgroup: str | None = None,
    ) -> 'MixedMultiGraph | None':
        setorder_list = list(circular_setorder.setorder)
        ret_frozen = frozenset(ret_set) if not isinstance(ret_set, frozenset) else ret_set

        ret_idx = None
        for i, s in enumerate(setorder_list):
            if frozenset(s) == ret_frozen:
                ret_idx = i
                break
        if ret_idx is None:
            return None

        test_graph = graph.copy()
        n_cycle = len(cycle_nodes)
        prev_idx = (ret_idx - 1) % n_cycle
        next_idx = (ret_idx + 1) % n_cycle

        hybrid_node = cycle_nodes[ret_idx]
        parent1 = cycle_nodes[prev_idx]
        parent2 = cycle_nodes[next_idx]

        test_graph.remove_edge(parent1, hybrid_node)
        test_graph.remove_edge(parent2, hybrid_node)
        test_graph.add_directed_edge(parent1, hybrid_node, gamma=0.5)
        test_graph.add_directed_edge(parent2, hybrid_node, gamma=0.5)

        components = source_components(test_graph)
        if len(components) != 1:
            return None

        if outgroup is not None:
            nodes_in_component, _, _ = components[0]
            outgroup_node = network._label_to_node.get(outgroup)
            if outgroup_node is None or outgroup_node not in nodes_in_component:
                return None

        return test_graph

    graph_with_cycle, cycle_nodes, set_to_node = replace_vertex_with_cycle(
        graph_copy, vertex, circular_setorder, edge_taxa_list, network._graph
    )

    if reticulation_ranking is None:
        return sdnetwork_from_graph(graph_with_cycle, network_type='mixed')

    for ret_set in reticulation_ranking:
        ret_frozen = frozenset(ret_set) if not isinstance(ret_set, frozenset) else ret_set
        hybrid_graph = try_hybrid(
            graph_with_cycle, ret_frozen, circular_setorder, cycle_nodes, outgroup
        )
        if hybrid_graph is not None:
            return sdnetwork_from_graph(hybrid_graph, network_type='semi-directed')

    raise PhyloZooAlgorithmError(
        "No valid hybrid configuration found from reticulation_ranking that "
        "maintains exactly one source component"
    )


def resolve_cycles(
    profileset: 'SqQuartetProfileSet',
    tree: SemiDirectedPhyNetwork,
    outgroup: str | None = None,
    rho: tuple[float, float, float, float] = (0.5, 1.0, 0.5, 1.0),
    tsp_threshold: int | None = 13,
    weighted_distance: bool = False,
) -> SemiDirectedPhyNetwork:
    """
    Convert a tree to a level-1 network by replacing high-degree vertices with cycles.

    For each internal vertex with degree > 3 (processed highest-degree first):
    1. Compute the partition induced by that vertex.
    2. Solve TSP to get the circular set ordering.
    3. Get hybrid ranking from quartet profiles.
    4. Insert a directed cycle with the best hybrid configuration.

    Parameters
    ----------
    profileset : SqQuartetProfileSet
    tree : SemiDirectedPhyNetwork
    outgroup : str | None
    rho : tuple of 4 floats
    tsp_threshold : int | None
        Max partition size for optimal TSP; larger uses simulated annealing.
    weighted_distance : bool
        If True, profile weights are used to scale TSP distance contributions
        so that high-confidence profiles pull the circular ordering more than
        low-confidence ones. Default: False.

    Returns
    -------
    SemiDirectedPhyNetwork
    """
    network = tree
    cut_verts = cut_vertices(network)
    if not cut_verts:
        return network

    nodes_by_degree = sorted(cut_verts, key=lambda v: network.degree(v), reverse=True)

    for vertex in nodes_by_degree:
        if network.degree(vertex) <= 3:
            continue

        induced_partition = partition_from_blob(network, {vertex}, return_edge_taxa=False)

        tsp_method = (
            'optimal'
            if (tsp_threshold is None or len(induced_partition) <= tsp_threshold)
            else 'simulated_annealing'
        )
        circular_setorder = _qprofiles_to_circular_ordering(
            profileset, induced_partition, rho=rho, tsp_method=tsp_method, weighted_distance=weighted_distance
        )
        hybrid_ranking = _qprofiles_to_hybrid_ranking(profileset, induced_partition, weights=True)

        network = _insert_cycle(
            network, vertex, circular_setorder,
            reticulation_ranking=hybrid_ranking,
            outgroup=outgroup,
        )

        if not isinstance(network, SemiDirectedPhyNetwork):
            raise PhyloZooAlgorithmError(
                f"Network became MixedPhyNetwork after inserting cycle at vertex {vertex}"
            )

    return network
