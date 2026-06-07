"""
Tests for the main squirrel algorithm module.
"""

import pytest
import numpy as np

from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
from phylozoo.core.network.dnetwork import DirectedPhyNetwork
from phylozoo.core.quartet.base import Quartet
from phylozoo.core.split.base import Split
from phylozoo.core.network.sdnetwork.derivations import root_at_outgroup
from phylozoo import DistanceMatrix
from phylozoo.utils.parallel import ParallelConfig, ParallelBackend
from physquirrel import (
    squirrel,
    SqQuartetProfile,
    SqQuartetProfileSet,
    sqprofileset_from_network,
    sqprofileset_similarity,
    delta_heuristic,
)


class TestSquirrelMain:
    """Test cases for the main squirrel() function."""

    def test_squirrel_basic(self) -> None:
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        profile1 = SqQuartetProfile([q1])
        profileset = SqQuartetProfileSet(profiles=[profile1])

        network = squirrel(profileset)

        assert isinstance(network, SemiDirectedPhyNetwork)
        assert network.number_of_nodes() > 0
        assert len(network.taxa) == 4
        network.validate()

    def test_squirrel_with_outgroup(self) -> None:
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        profile1 = SqQuartetProfile([q1])
        profileset = SqQuartetProfileSet(profiles=[profile1])

        network = squirrel(profileset, outgroup='A')

        assert isinstance(network, DirectedPhyNetwork)
        assert network.root_node is not None
        network.validate()

    def test_squirrel_with_cycle_profiles(self) -> None:
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        q2 = Quartet(Split({'A', 'C'}, {'B', 'D'}))
        profile1 = SqQuartetProfile([q1, q2], reticulation_leaf='A')
        profileset = SqQuartetProfileSet(profiles=[profile1])

        network = squirrel(profileset)

        assert isinstance(network, SemiDirectedPhyNetwork)
        assert len(network.taxa) == 4
        network.validate()

    def test_squirrel_with_kwargs(self) -> None:
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        profile1 = SqQuartetProfile([q1])
        profileset = SqQuartetProfileSet(profiles=[profile1])

        network = squirrel(
            profileset,
            rho=(0.4, 0.9, 0.6, 1.1),
            tsp_threshold=10
        )

        assert isinstance(network, SemiDirectedPhyNetwork)
        network.validate()

    def test_squirrel_larger_profileset(self) -> None:
        quartets = [
            Quartet(Split({'A', 'B'}, {'C', 'D'})),
            Quartet(Split({'A', 'B'}, {'C', 'E'})),
            Quartet(Split({'A', 'B'}, {'D', 'E'})),
            Quartet(Split({'A', 'C'}, {'D', 'E'})),
            Quartet(Split({'B', 'C'}, {'D', 'E'}))
        ]
        profiles = [SqQuartetProfile([q]) for q in quartets]
        profileset = SqQuartetProfileSet(profiles=profiles)

        network = squirrel(profileset)

        assert isinstance(network, SemiDirectedPhyNetwork)
        assert len(network.taxa) == 5
        network.validate()

    def test_squirrel_returns_best_network(self) -> None:
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        profile1 = SqQuartetProfile([q1])
        profileset = SqQuartetProfileSet(profiles=[profile1])

        network = squirrel(profileset)

        reconstructed_profileset = sqprofileset_from_network(network)
        similarity = sqprofileset_similarity(profileset, reconstructed_profileset)

        assert 0.0 <= similarity <= 1.0


class TestRootAtOutgroup:
    """Test cases for root_at_outgroup function."""

    def test_root_at_outgroup_basic(self) -> None:
        network = SemiDirectedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[
                (1, {'label': 'A'}),
                (2, {'label': 'B'}),
                (4, {'label': 'C'})
            ]
        )

        rooted = root_at_outgroup(network, 'A')

        assert isinstance(rooted, DirectedPhyNetwork)
        assert rooted.root_node is not None
        assert 'A' in rooted.taxa
        rooted.validate()

    def test_root_at_outgroup_different_leaf(self) -> None:
        network = SemiDirectedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[
                (1, {'label': 'A'}),
                (2, {'label': 'B'}),
                (4, {'label': 'C'})
            ]
        )

        rooted = root_at_outgroup(network, 'B')

        assert isinstance(rooted, DirectedPhyNetwork)
        assert rooted.root_node is not None
        assert 'B' in rooted.taxa
        rooted.validate()

    def test_root_at_outgroup_with_hybrid(self) -> None:
        network = SemiDirectedPhyNetwork(
            directed_edges=[(5, 4), (6, 4)],
            undirected_edges=[
                (5, 1), (5, 6),
                (6, 2),
                (4, 3)
            ],
            nodes=[
                (1, {'label': 'A'}),
                (2, {'label': 'B'}),
                (3, {'label': 'C'})
            ]
        )

        rooted = root_at_outgroup(network, 'A')

        assert isinstance(rooted, DirectedPhyNetwork)
        assert rooted.root_node is not None
        assert 'A' in rooted.taxa
        rooted.validate()

    def test_root_at_outgroup_invalid_taxon(self) -> None:
        network = SemiDirectedPhyNetwork(
            undirected_edges=[(5, 1), (5, 2), (5, 3)],
            nodes=[
                (1, {'label': 'A'}),
                (2, {'label': 'B'}),
                (3, {'label': 'C'})
            ]
        )

        with pytest.raises(ValueError, match="not found"):
            root_at_outgroup(network, 'X')


class TestSqProfilesetFromNetwork:
    """Test cases for sqprofileset_from_network edge cases."""

    def test_sqprofileset_from_network_fewer_than_four_taxa(self) -> None:
        network = SemiDirectedPhyNetwork(
            undirected_edges=[(5, 1), (5, 2), (5, 3)],
            nodes=[
                (1, {'label': 'A'}),
                (2, {'label': 'B'}),
                (3, {'label': 'C'})
            ]
        )

        profileset = sqprofileset_from_network(network)
        assert len(profileset) == 0

    def test_sqprofileset_from_network_non_binary_error(self) -> None:
        network = SemiDirectedPhyNetwork(
            undirected_edges=[
                (5, 1), (5, 2), (5, 3), (5, 4), (5, 6)
            ],
            nodes=[
                (1, {'label': 'A'}),
                (2, {'label': 'B'}),
                (3, {'label': 'C'}),
                (4, {'label': 'D'}),
                (6, {'label': 'E'})
            ]
        )

        with pytest.raises(ValueError, match="binary"):
            sqprofileset_from_network(network)

    def test_sqprofileset_from_network_parallel_edges_error(self) -> None:
        from tests.fixtures.sd_networks import LEVEL_1_SDNETWORK_PARALLEL_EDGES

        with pytest.raises(ValueError, match="parallel edges"):
            sqprofileset_from_network(LEVEL_1_SDNETWORK_PARALLEL_EDGES)

    def test_sqprofileset_from_network_level_greater_than_one_error(self) -> None:
        from tests.fixtures.sd_networks import LEVEL_2_SDNETWORK_DIAMOND_HYBRID

        with pytest.raises(ValueError, match="(level-1|binary)"):
            sqprofileset_from_network(LEVEL_2_SDNETWORK_DIAMOND_HYBRID)


class TestSqProfilesetSimilarity:
    """Test cases for sqprofileset_similarity function."""

    def test_sqprofileset_similarity_identical(self) -> None:
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        profile1 = SqQuartetProfile([q1])
        profileset = SqQuartetProfileSet(profiles=[profile1])

        similarity = sqprofileset_similarity(profileset, profileset)
        assert similarity == 1.0

    def test_sqprofileset_similarity_different(self) -> None:
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        profile1 = SqQuartetProfile([q1])
        profileset1 = SqQuartetProfileSet(profiles=[profile1])

        q2 = Quartet(Split({'A', 'C'}, {'B', 'D'}))
        profile2 = SqQuartetProfile([q2])
        profileset2 = SqQuartetProfileSet(profiles=[profile2])

        similarity = sqprofileset_similarity(profileset1, profileset2)
        assert 0.0 <= similarity <= 1.0
        assert similarity < 1.0

    def test_sqprofileset_similarity_empty(self) -> None:
        profileset1 = SqQuartetProfileSet()
        profileset2 = SqQuartetProfileSet()

        similarity = sqprofileset_similarity(profileset1, profileset2)
        assert similarity == 1.0

    def test_sqprofileset_similarity_weighted(self) -> None:
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        profile1 = SqQuartetProfile([q1])
        profileset1 = SqQuartetProfileSet(profiles=[profile1])
        profileset2 = SqQuartetProfileSet(profiles=[profile1])

        similarity_weighted = sqprofileset_similarity(profileset1, profileset2, weighted=True)
        similarity_unweighted = sqprofileset_similarity(profileset1, profileset2, weighted=False)

        assert 0.0 <= similarity_weighted <= 1.0
        assert 0.0 <= similarity_unweighted <= 1.0

    def test_sqprofileset_similarity_partial_overlap(self) -> None:
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        q2 = Quartet(Split({'A', 'C'}, {'B', 'D'}))

        profile1 = SqQuartetProfile([q1])
        profile2 = SqQuartetProfile([q2])

        profileset1 = SqQuartetProfileSet(profiles=[profile1])
        profileset2 = SqQuartetProfileSet(profiles=[profile2])

        similarity = sqprofileset_similarity(profileset1, profileset2)
        assert 0.0 <= similarity <= 1.0


class TestSquirrelIntegration:
    """Integration tests for the full squirrel pipeline."""

    def test_squirrel_pipeline_from_distance_matrix(self) -> None:
        matrix = np.array([
            [0.0, 0.1, 0.9, 0.9],
            [0.1, 0.0, 0.9, 0.9],
            [0.9, 0.9, 0.0, 0.1],
            [0.9, 0.9, 0.1, 0.0]
        ])
        dm = DistanceMatrix(matrix, labels=['A', 'B', 'C', 'D'])

        profileset = delta_heuristic(dm, lam=0.3)
        assert len(profileset) > 0

        network = squirrel(profileset)

        assert isinstance(network, SemiDirectedPhyNetwork)
        assert len(network.taxa) == 4
        network.validate()

    def test_squirrel_pipeline_with_outgroup(self) -> None:
        matrix = np.array([
            [0.0, 0.1, 0.9, 0.9],
            [0.1, 0.0, 0.9, 0.9],
            [0.9, 0.9, 0.0, 0.1],
            [0.9, 0.9, 0.1, 0.0]
        ])
        dm = DistanceMatrix(matrix, labels=['A', 'B', 'C', 'D'])

        profileset = delta_heuristic(dm, lam=0.3)
        network = squirrel(profileset, outgroup='A')

        assert isinstance(network, DirectedPhyNetwork)
        assert network.root_node is not None
        network.validate()

    def test_squirrel_handles_empty_profileset(self) -> None:
        empty_profileset = SqQuartetProfileSet()

        with pytest.raises((ValueError, RuntimeError)):
            squirrel(empty_profileset)

    def test_squirrel_handles_single_quartet(self) -> None:
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        profile1 = SqQuartetProfile([q1])
        profileset = SqQuartetProfileSet(profiles=[profile1])

        network = squirrel(profileset)

        assert isinstance(network, SemiDirectedPhyNetwork)
        assert len(network.taxa) == 4
        network.validate()

    def test_squirrel_consistency(self) -> None:
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        profile1 = SqQuartetProfile([q1])
        profileset = SqQuartetProfileSet(profiles=[profile1])

        network1 = squirrel(profileset)
        network2 = squirrel(profileset)

        assert isinstance(network1, SemiDirectedPhyNetwork)
        assert isinstance(network2, SemiDirectedPhyNetwork)
        network1.validate()
        network2.validate()

        assert network1.taxa == network2.taxa


class TestRepresentativeMode:
    """Tests for the representative_mode parameter of squirrel()."""

    def test_default_mode_is_best(self, five_taxon_profileset) -> None:
        net = squirrel(five_taxon_profileset)
        assert isinstance(net, SemiDirectedPhyNetwork)
        net.validate()

    def test_best_mode_works(self, five_taxon_profileset) -> None:
        net = squirrel(five_taxon_profileset, representative_mode='best')
        assert isinstance(net, SemiDirectedPhyNetwork)
        net.validate()

    def test_cycle_profileset_both_modes_return_valid_network(
        self, cycle_four_taxon_profileset
    ) -> None:
        for mode in ('average', 'best'):
            net = squirrel(cycle_four_taxon_profileset, representative_mode=mode)
            assert isinstance(net, SemiDirectedPhyNetwork)
            net.validate()

    def test_tree_profileset_no_hybrid_either_mode(
        self, simple_four_taxon_profileset
    ) -> None:
        for mode in ('average', 'best'):
            net = squirrel(simple_four_taxon_profileset, representative_mode=mode)
            assert len(net.hybrid_nodes) == 0

    def test_tree_profileset_perfect_score_either_mode(
        self, simple_four_taxon_profileset
    ) -> None:
        ps = simple_four_taxon_profileset
        for mode in ('average', 'best'):
            net = squirrel(ps, representative_mode=mode)
            score = sqprofileset_similarity(ps, sqprofileset_from_network(net))
            assert score == pytest.approx(1.0, abs=1e-6)


class TestParallelExecution:
    """Tests for the parallel parameter of squirrel()."""

    def test_parallel_none_works(self, five_taxon_profileset) -> None:
        net = squirrel(five_taxon_profileset, parallel=None)
        assert isinstance(net, SemiDirectedPhyNetwork)
        net.validate()

    def test_threading_returns_valid_network(self, five_taxon_profileset) -> None:
        cfg = ParallelConfig(backend=ParallelBackend.THREADING, n_jobs=2)
        net = squirrel(five_taxon_profileset, parallel=cfg)
        assert isinstance(net, SemiDirectedPhyNetwork)
        net.validate()

    def test_threading_matches_sequential_taxa(self, five_taxon_profileset) -> None:
        net_seq = squirrel(five_taxon_profileset, parallel=None, representative_mode='best')
        cfg = ParallelConfig(backend=ParallelBackend.THREADING, n_jobs=2)
        net_par = squirrel(five_taxon_profileset, parallel=cfg, representative_mode='best')
        assert net_seq.taxa == net_par.taxa

    def test_multiprocessing_returns_valid_network(self, five_taxon_profileset) -> None:
        cfg = ParallelConfig(backend=ParallelBackend.MULTIPROCESSING, n_jobs=2)
        net = squirrel(five_taxon_profileset, parallel=cfg)
        assert isinstance(net, SemiDirectedPhyNetwork)
        net.validate()

    def test_multiprocessing_matches_sequential_taxa(self, five_taxon_profileset) -> None:
        net_seq = squirrel(five_taxon_profileset, parallel=None, representative_mode='best')
        cfg = ParallelConfig(backend=ParallelBackend.MULTIPROCESSING, n_jobs=2)
        net_par = squirrel(five_taxon_profileset, parallel=cfg, representative_mode='best')
        assert net_seq.taxa == net_par.taxa
