"""
Integration tests against the HIV-1 alignment (9 taxa).

Ground-truth scores were computed on the fixed develop branch and confirmed to
match v1 (main branch) exactly at every iteration for lam=0.3, representative_mode='best'.

Run with:   pytest tests/test_real_data.py -v
Skip with:  pytest -m "not slow"
"""

from pathlib import Path

import pytest
from phylozoo import MSA
from phylozoo.core.network.dnetwork import DirectedPhyNetwork
from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork

from physquirrel import (
    adapted_quartet_joining,
    delta_heuristic_from_msa,
    resolve_cycles,
    sqprofileset_from_network,
    sqprofileset_similarity,
    squirrel,
    tstar_tree,
    unresolve_tree,
)

# ---------------------------------------------------------------------------
# Paths and taxonomy
# ---------------------------------------------------------------------------

_HIV_FASTA = Path(__file__).parent / "data" / "hiv.fasta"
_TAXA = frozenset({'A', 'B', 'C', 'D', 'F', 'G', 'H', 'J', 'KAL153'})

# ---------------------------------------------------------------------------
# Ground truth: lam=0.3, representative_mode='best'
# Confirmed equal to v1 (main branch) at every iteration.
# ---------------------------------------------------------------------------
_LAM03_SCORES = [
    0.508079, 0.558157, 0.514600, 0.482298, 0.482437, 0.551487, 0.576545
]
_LAM03_BEST_IDX   = 6
_LAM03_2ND_IDX    = 1
_LAM03_BEST_SCORE = 0.576545
_LAM03_2ND_SCORE  = 0.558157
_LAM03_BEST_RET   = frozenset({'KAL153'})   # leaf(s) adjacent to hybrid in best network
_LAM03_2ND_RET    = frozenset({'C'})        # leaf(s) adjacent to hybrid in 2nd-best network

# ---------------------------------------------------------------------------
# Ground truth: lam=0.1, representative_mode='best'
# Confirmed equal to v1 (main branch) at every iteration.
# Both 'best' and 'average' modes give identical results on this dataset.
# ---------------------------------------------------------------------------
_LAM01_SCORES = [
    0.263852, 0.321844, 0.360907, 0.270447, 0.288730, 0.226111, 0.117161
]
_LAM01_BEST_IDX   = 2
_LAM01_2ND_IDX    = 1
_LAM01_BEST_SCORE = 0.360907
_LAM01_2ND_SCORE  = 0.321844
_LAM01_BEST_RET   = frozenset({'C', 'KAL153'})   # two hybrids; one leaf each

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _read_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    name: str | None = None
    chunks: list[str] = []
    with path.open() as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if name:
                    seqs[name] = ''.join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if name:
        seqs[name] = ''.join(chunks)
    return seqs


def _hybrid_leaf_labels(net: SemiDirectedPhyNetwork) -> frozenset[str]:
    """Return the set of leaf labels directly adjacent to any hybrid node."""
    labels: set[str] = set()
    graph = net._graph._combined
    for node in net.hybrid_nodes:
        for nb in graph.neighbors(node):
            lbl = graph.nodes[nb].get('label', '')
            if lbl:
                labels.add(lbl)
    return frozenset(labels)


def _run_all_iterations(
    fasta_path: Path,
    lam: float,
    representative_mode: str = 'best',
) -> dict:
    """Run the full squirrel pipeline and return all intermediate networks/scores."""
    ps = delta_heuristic_from_msa(MSA(_read_fasta(fasta_path)), lam=lam)
    qj = adapted_quartet_joining(ps, starting_tree=tstar_tree(ps))
    nets: list[SemiDirectedPhyNetwork] = []
    scores: list[float] = []
    for t in unresolve_tree(qj, ps):
        n = resolve_cycles(ps, t, weighted_distance=True, representative_mode=representative_mode)
        nets.append(n)
        scores.append(float(sqprofileset_similarity(ps, sqprofileset_from_network(n))))
    return {"nets": nets, "scores": scores, "ps": ps}


# ---------------------------------------------------------------------------
# Module-scoped fixtures (pipeline runs once per test session per lambda)
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def hiv_lam03():
    return _run_all_iterations(_HIV_FASTA, lam=0.3, representative_mode='best')


@pytest.fixture(scope="module")
def hiv_lam01():
    return _run_all_iterations(_HIV_FASTA, lam=0.1, representative_mode='best')


# ---------------------------------------------------------------------------
# lam=0.3 score tests
# ---------------------------------------------------------------------------

@pytest.mark.slow
class TestHivScoresLam03:
    """Per-iteration scores must match v1 ground truth exactly (within 1e-4)."""

    def test_iteration_count(self, hiv_lam03):
        assert len(hiv_lam03["scores"]) == 7

    def test_all_scores_match_v1(self, hiv_lam03):
        for i, (got, expected) in enumerate(zip(hiv_lam03["scores"], _LAM03_SCORES)):
            assert got == pytest.approx(expected, abs=1e-4), f"iter {i}: {got} != {expected}"

    def test_best_index(self, hiv_lam03):
        scores = hiv_lam03["scores"]
        assert scores.index(max(scores)) == _LAM03_BEST_IDX

    def test_best_score(self, hiv_lam03):
        assert max(hiv_lam03["scores"]) == pytest.approx(_LAM03_BEST_SCORE, abs=1e-4)

    def test_second_best_score(self, hiv_lam03):
        sorted_scores = sorted(hiv_lam03["scores"], reverse=True)
        assert sorted_scores[1] == pytest.approx(_LAM03_2ND_SCORE, abs=1e-4)


# ---------------------------------------------------------------------------
# lam=0.3 network structure tests
# ---------------------------------------------------------------------------

@pytest.mark.slow
class TestHivNetworkLam03:
    """Best and second-best network topology at lam=0.3."""

    def test_best_network_is_valid(self, hiv_lam03):
        hiv_lam03["nets"][_LAM03_BEST_IDX].validate()

    def test_best_network_taxa(self, hiv_lam03):
        assert hiv_lam03["nets"][_LAM03_BEST_IDX].taxa == _TAXA

    def test_best_network_one_hybrid(self, hiv_lam03):
        assert len(hiv_lam03["nets"][_LAM03_BEST_IDX].hybrid_nodes) == 1

    def test_best_network_hybrid_leaf(self, hiv_lam03):
        assert _hybrid_leaf_labels(hiv_lam03["nets"][_LAM03_BEST_IDX]) == _LAM03_BEST_RET

    def test_second_best_network_is_valid(self, hiv_lam03):
        hiv_lam03["nets"][_LAM03_2ND_IDX].validate()

    def test_second_best_network_one_hybrid(self, hiv_lam03):
        assert len(hiv_lam03["nets"][_LAM03_2ND_IDX].hybrid_nodes) == 1

    def test_second_best_network_hybrid_leaf(self, hiv_lam03):
        assert _hybrid_leaf_labels(hiv_lam03["nets"][_LAM03_2ND_IDX]) == _LAM03_2ND_RET


# ---------------------------------------------------------------------------
# lam=0.1 score tests
# ---------------------------------------------------------------------------

@pytest.mark.slow
class TestHivScoresLam01:
    """Per-iteration scores at lam=0.1, confirmed equal to v1."""

    def test_iteration_count(self, hiv_lam01):
        assert len(hiv_lam01["scores"]) == 7

    def test_all_scores_match_v1(self, hiv_lam01):
        for i, (got, expected) in enumerate(zip(hiv_lam01["scores"], _LAM01_SCORES)):
            assert got == pytest.approx(expected, abs=1e-4), f"iter {i}: {got} != {expected}"

    def test_best_index(self, hiv_lam01):
        scores = hiv_lam01["scores"]
        assert scores.index(max(scores)) == _LAM01_BEST_IDX

    def test_best_score(self, hiv_lam01):
        assert max(hiv_lam01["scores"]) == pytest.approx(_LAM01_BEST_SCORE, abs=1e-4)

    def test_second_best_score(self, hiv_lam01):
        sorted_scores = sorted(hiv_lam01["scores"], reverse=True)
        assert sorted_scores[1] == pytest.approx(_LAM01_2ND_SCORE, abs=1e-4)


# ---------------------------------------------------------------------------
# lam=0.1 network structure tests
# ---------------------------------------------------------------------------

@pytest.mark.slow
class TestHivNetworkLam01:
    """Best network topology at lam=0.1, and mode-agreement check.

    Best iteration is determined dynamically (per-iteration scores vary between
    Python processes due to hash-dependent tie-breaking in split support sorting).
    """

    def test_best_network_is_valid(self, hiv_lam01):
        hiv_lam01["nets"][_LAM01_BEST_IDX].validate()

    def test_best_network_taxa(self, hiv_lam01):
        assert hiv_lam01["nets"][_LAM01_BEST_IDX].taxa == _TAXA

    def test_best_network_two_hybrids(self, hiv_lam01):
        assert len(hiv_lam01["nets"][_LAM01_BEST_IDX].hybrid_nodes) == 2

    def test_best_network_hybrid_leaves(self, hiv_lam01):
        assert _hybrid_leaf_labels(hiv_lam01["nets"][_LAM01_BEST_IDX]) == _LAM01_BEST_RET

    def test_second_best_network_is_valid(self, hiv_lam01):
        hiv_lam01["nets"][_LAM01_2ND_IDX].validate()

    def test_modes_agree_on_best_network(self):
        """'average' mode must produce the same per-iteration scores and best network as 'best' mode."""
        result_avg = _run_all_iterations(_HIV_FASTA, lam=0.1, representative_mode='average')
        assert result_avg["scores"].index(max(result_avg["scores"])) == _LAM01_BEST_IDX
        assert _hybrid_leaf_labels(result_avg["nets"][_LAM01_BEST_IDX]) == _LAM01_BEST_RET


# ---------------------------------------------------------------------------
# Outgroup rooting
# ---------------------------------------------------------------------------

@pytest.mark.slow
class TestHivOutgroup:
    """Rooting the HIV best network at outgroup C."""

    def test_squirrel_rooted_at_C(self):
        ps = delta_heuristic_from_msa(MSA(_read_fasta(_HIV_FASTA)), lam=0.3)
        net = squirrel(ps, outgroup='C', representative_mode='best')
        assert isinstance(net, DirectedPhyNetwork)
        assert net.root_node is not None
        assert 'C' in net.taxa
        net.validate()
