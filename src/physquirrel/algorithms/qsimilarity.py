from typing import Any, Iterator, TYPE_CHECKING

from phylozoo.core.quartet.base import Quartet
from phylozoo.core.network.sdnetwork.derivations import induced_splits, partition_from_blob
from phylozoo.core.network.sdnetwork.classifications import level, is_binary, has_parallel_edges
from phylozoo.core.network.sdnetwork.features import blobs
from phylozoo.core.split.base import Split
from phylozoo.core.split.algorithms import induced_quartetsplits
from phylozoo.core.primitives.circular_ordering import CircularSetOrdering
from phylozoo.utils.exceptions import PhyloZooValueError, PhyloZooAlgorithmError

from ..datatypes.sqprofile import SqQuartetProfile
from ..datatypes.sqprofileset import SqQuartetProfileSet

if TYPE_CHECKING:
    from phylozoo.core.network.sdnetwork import MixedPhyNetwork, SemiDirectedPhyNetwork


def _circular_orders_from_cycles(
    network: 'MixedPhyNetwork',
) -> Iterator[tuple[CircularSetOrdering, frozenset[str] | None]]:
    """
    Yield (circular_set_ordering, reticulation_set) for each cycle in the network.

    Requires a binary, level-1, parallel-edge-free network.
    """
    if not is_binary(network):
        raise PhyloZooValueError("Network must be binary for _circular_orders_from_cycles")
    if has_parallel_edges(network):
        raise PhyloZooValueError("Network must not have parallel edges")
    net_level = level(network)
    if net_level > 1:
        raise PhyloZooValueError(
            f"Network must be level-1 or lower, but got level {net_level}"
        )

    def _find_cycle_in_blob(graph: Any, blob: set[Any]) -> list[Any]:
        start = next(iter(blob))
        cycle_nodes: list[Any] = [start]
        visited: set[Any] = {start}
        current = start
        previous = None

        while len(cycle_nodes) < len(blob):
            next_node = None
            for neighbor in graph.neighbors(current):
                if neighbor in blob and neighbor != previous:
                    next_node = neighbor
                    break
            if next_node is None:
                raise PhyloZooAlgorithmError(
                    f"Could not find next node in blob {blob} from {current}."
                )
            cycle_nodes.append(next_node)
            visited.add(next_node)
            previous = current
            current = next_node

        if start not in graph.neighbors(current):
            raise PhyloZooAlgorithmError(
                f"Cycle in blob {blob} is not closed."
            )
        return cycle_nodes

    network_blobs = blobs(network, trivial=False, leaves=False)
    combined_graph = network._graph._combined

    for blob in network_blobs:
        if len(blob) < 3:
            raise PhyloZooAlgorithmError(
                f"Blob {blob} has fewer than 3 nodes."
            )

        cycle_nodes = _find_cycle_in_blob(combined_graph, blob)

        network_hybrid_nodes = network.hybrid_nodes
        cycle_hybrids = [node for node in cycle_nodes if node in network_hybrid_nodes]
        hybrid: Any | None = cycle_hybrids[0] if cycle_hybrids else None

        _, edge_taxa_list = partition_from_blob(network, blob, return_edge_taxa=True)
        cycle_order_dict = {node: i for i, node in enumerate(cycle_nodes)}
        sorted_edge_taxa = sorted(edge_taxa_list, key=lambda x: cycle_order_dict.get(x[1], -1))
        taxa_sets = [taxa_set for _, _, taxa_set in sorted_edge_taxa]

        ret_set: frozenset[str] | None = None
        if hybrid is not None:
            for (_, v, taxa_set) in sorted_edge_taxa:
                if v == hybrid:
                    ret_set = taxa_set
                    break
            if ret_set is None:
                raise PhyloZooAlgorithmError(
                    f"Could not find reticulation set for hybrid {hybrid} in cycle {blob}"
                )

        set_ordering = CircularSetOrdering(taxa_sets)
        yield set_ordering, ret_set


def sqprofileset_similarity(
    profileset1: SqQuartetProfileSet,
    profileset2: SqQuartetProfileSet,
    weighted: bool = True,
) -> float:
    """
    Compute the C-measure (consistency) between two SqQuartetProfileSets.

    Two profiles are consistent if they have the same quartets and reticulation_leaf.

    Parameters
    ----------
    profileset1 : SqQuartetProfileSet
    profileset2 : SqQuartetProfileSet
    weighted : bool
        Use profile weights. Default True.

    Returns
    -------
    float in [0, 1]. Returns 1.0 if profileset1 is empty.
    """
    if profileset1.taxa != profileset2.taxa:
        raise PhyloZooValueError(
            f"Profile sets must have the same taxa. "
            f"Set 1: {profileset1.taxa}, Set 2: {profileset2.taxa}"
        )

    if len(profileset1) == 0:
        return 1.0

    if weighted:
        consistent_weight = 0.0
        total_weight = 0.0
        for taxa_set, (profile1, weight1) in profileset1.profiles.items():
            total_weight += weight1
            profile2 = profileset2.get_profile(taxa_set)
            if profile2 is not None:
                profile2_full = profileset2.profiles[taxa_set][0]
                if (
                    profile1._quartets == profile2_full._quartets
                    and profile1.reticulation_leaf == profile2_full.reticulation_leaf
                ):
                    consistent_weight += weight1
        return consistent_weight / total_weight if total_weight > 0 else 1.0
    else:
        consistent_count = 0
        total_count = len(profileset1)
        for taxa_set, (profile1, _) in profileset1.profiles.items():
            profile2 = profileset2.get_profile(taxa_set)
            if profile2 is not None and (
                profile1._quartets == profile2._quartets
                and profile1.reticulation_leaf == profile2.reticulation_leaf
            ):
                consistent_count += 1
        return consistent_count / total_count if total_count > 0 else 1.0


def sqprofileset_from_network(
    network: 'SemiDirectedPhyNetwork',
) -> SqQuartetProfileSet:
    """
    Compute the SqQuartetProfileSet displayed by a semi-directed level-1 network.

    Requires the network to be binary, level-1, and parallel-edge-free.

    Parameters
    ----------
    network : SemiDirectedPhyNetwork

    Returns
    -------
    SqQuartetProfileSet
    """
    taxa_list = sorted(network.taxa)
    if len(taxa_list) < 4:
        return SqQuartetProfileSet()

    if not is_binary(network):
        raise PhyloZooValueError("Network must be binary for sqprofileset_from_network")
    if has_parallel_edges(network):
        raise PhyloZooValueError("Network must not have parallel edges")
    net_level = level(network)
    if net_level > 1:
        raise PhyloZooValueError(
            f"Network must be level-1 or lower, but got level {net_level}"
        )

    sq_quartet_profiles: list[SqQuartetProfile] = []

    # Step 1: cycle-based quartets
    for set_ordering, ret_set in _circular_orders_from_cycles(network):
        if len(set_ordering) < 4:
            continue
        for sub_order in set_ordering.suborderings(4):
            if ret_set is None:
                for repr_order in sub_order.representative_orderings():
                    a, b, c, d = list(repr_order.order)
                    q1 = Quartet(Split({a, b}, {c, d}))
                    q2 = Quartet(Split({a, d}, {b, c}))
                    sq_quartet_profiles.append(SqQuartetProfile([q1, q2]))
            elif any(ret_set == part for part in sub_order.parts):
                for repr_order in sub_order.representative_orderings():
                    a, b, c, d = list(repr_order.order)
                    ret_leaf = next(iter(ret_set & repr_order.elements))
                    q1 = Quartet(Split({a, b}, {c, d}))
                    q2 = Quartet(Split({a, d}, {b, c}))
                    sq_quartet_profiles.append(
                        SqQuartetProfile([q1, q2], reticulation_leaf=ret_leaf)
                    )
            else:
                for repr_order in sub_order.representative_orderings():
                    a, b, c, d = list(repr_order.order)
                    quartet = Quartet(Split({a, b}, {c, d}))
                    sq_quartet_profiles.append(SqQuartetProfile([quartet]))

    # Step 2: split-based quartets
    split_system = induced_splits(network)
    seen_quartet_splits: set[Split] = set()
    for split in split_system.splits:
        if split.is_trivial:
            continue
        for quartet_split in induced_quartetsplits(split):
            if quartet_split in seen_quartet_splits:
                continue
            seen_quartet_splits.add(quartet_split)
            sq_quartet_profiles.append(SqQuartetProfile([Quartet(quartet_split)]))

    return SqQuartetProfileSet(profiles=sq_quartet_profiles)
