from typing import Any

from phylozoo import SemiDirectedPhyNetwork, DirectedPhyNetwork
from phylozoo.core.network.sdnetwork.derivations import root_at_outgroup
from phylozoo.utils.parallel import ParallelConfig, ParallelBackend

from .cycle_resolution import resolve_cycles
from .qjoining import adapted_quartet_joining
from .qsimilarity import sqprofileset_from_network, sqprofileset_similarity
from ..datatypes.sqprofileset import SqQuartetProfileSet
from .tstar_tree import tstar_tree
from .unresolve_tree import unresolve_tree


def _process_contracted_tree(
    args: tuple[SemiDirectedPhyNetwork | str, SqQuartetProfileSet | str, str | None, dict[str, Any]],
) -> tuple[SemiDirectedPhyNetwork | str, float]:
    """
    Resolve cycles on a contracted tree and compute similarity.

    Accepts either direct objects (threading/sequential) or serialized strings
    (multiprocessing).
    """
    tree_input, profileset_input, outgroup, kwargs = args
    is_serialized = isinstance(tree_input, str)

    if is_serialized:
        from phylozoo.core.network.sdnetwork import io as _sdnetwork_io  # noqa: F401
        from phylozoo import SemiDirectedPhyNetwork as _SdNet
        from ..datatypes import io  # noqa: F401 — registers pz format in worker
        from ..datatypes.sqprofileset import SqQuartetProfileSet as _SqPS

        contracted_tree = _SdNet.from_string(tree_input, format='phylozoo-dot')
        profileset = _SqPS.from_string(profileset_input, format='pz')  # type: ignore
    else:
        contracted_tree = tree_input  # type: ignore
        profileset = profileset_input  # type: ignore

    net_new = resolve_cycles(profileset, contracted_tree, outgroup=outgroup, **kwargs)
    quarnets_new = sqprofileset_from_network(net_new)
    score_new = sqprofileset_similarity(profileset, quarnets_new)

    if is_serialized:
        return (net_new.to_string(format='phylozoo-dot'), score_new)
    return (net_new, score_new)


def squirrel(
    profileset: SqQuartetProfileSet,
    outgroup: str | None = None,
    parallel: ParallelConfig | None = None,
    **kwargs: Any,
) -> SemiDirectedPhyNetwork | DirectedPhyNetwork:
    """
    Reconstruct a phylogenetic network from a squirrel quartet profile set.

    Pipeline:
    1. Compute T* tree.
    2. Apply adapted quartet joining.
    3. Iteratively contract least-supported splits (unresolve_tree).
    4. For each contracted tree: resolve cycles, score with C-measure.
    5. Return best network.

    Parameters
    ----------
    profileset : SqQuartetProfileSet
    outgroup : str | None
        If provided, returns a DirectedPhyNetwork rooted at the outgroup edge.
    parallel : ParallelConfig | None
        Parallelization config. None → sequential.
    **kwargs
        Passed to resolve_cycles (rho, tsp_threshold).

    Returns
    -------
    SemiDirectedPhyNetwork | DirectedPhyNetwork
    """
    tstar = tstar_tree(profileset)
    qj_tree = adapted_quartet_joining(profileset, starting_tree=tstar)

    if parallel is None or parallel.backend == ParallelBackend.SEQUENTIAL:
        networks: list[SemiDirectedPhyNetwork] = []
        scores: list[float] = []
        for contracted_tree in unresolve_tree(qj_tree, profileset):
            net_new, score_new = _process_contracted_tree(
                (contracted_tree, profileset, outgroup, kwargs)
            )
            networks.append(net_new)
            scores.append(score_new)

    elif parallel.backend == ParallelBackend.MULTIPROCESSING:
        contracted_trees = list(unresolve_tree(qj_tree, profileset))
        profileset_pz = profileset.to_string(format='pz')
        process_args = [
            (tree.to_string(format='phylozoo-dot'), profileset_pz, outgroup, kwargs)
            for tree in contracted_trees
        ]
        executor = parallel.get_executor()
        try:
            results = list(executor.map(_process_contracted_tree, process_args))
        finally:
            if hasattr(executor, '_pool') and executor._pool is not None:
                executor._pool.close()
                executor._pool.join()
        networks = [
            SemiDirectedPhyNetwork.from_string(net_pzdot, format='phylozoo-dot')
            for net_pzdot, _ in results
        ]
        scores = [score for _, score in results]

    else:  # THREADING
        contracted_trees = list(unresolve_tree(qj_tree, profileset))
        process_args = [
            (tree, profileset, outgroup, kwargs) for tree in contracted_trees
        ]
        executor = parallel.get_executor()
        try:
            results = list(executor.map(_process_contracted_tree, process_args))
        finally:
            if hasattr(executor, '_executor'):
                executor._executor.shutdown(wait=True)
        networks = [net for net, _ in results]
        scores = [score for _, score in results]

    best_idx = scores.index(max(scores))
    best_network = networks[best_idx]

    if outgroup is not None:
        best_network = root_at_outgroup(best_network, outgroup)

    return best_network
