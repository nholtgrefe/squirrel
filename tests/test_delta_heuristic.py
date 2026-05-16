"""
Tests for the delta heuristic module.
"""

import pytest
import numpy as np

from phylozoo import DistanceMatrix
from phylozoo.core.split.base import Split
from phylozoo.core.quartet.base import Quartet
from physquirrel import delta_heuristic, SqQuartetProfileSet, squirrel
from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork


class TestDeltaHeuristic:
    """Test cases for delta heuristic algorithm."""

    def test_fewer_than_four_taxa_returns_empty(self) -> None:
        matrix = np.array([[0.0, 0.1, 0.2], [0.1, 0.0, 0.3], [0.2, 0.3, 0.0]])
        dm = DistanceMatrix(matrix, labels=['A', 'B', 'C'])
        profileset = delta_heuristic(dm)
        assert len(profileset) == 0

    def test_four_taxa_creates_one_profile(self) -> None:
        matrix = np.array([
            [0.0, 0.1, 0.9, 0.9],
            [0.1, 0.0, 0.9, 0.9],
            [0.9, 0.9, 0.0, 0.1],
            [0.9, 0.9, 0.1, 0.0]
        ])
        dm = DistanceMatrix(matrix, labels=['A', 'B', 'C', 'D'])
        profileset = delta_heuristic(dm, lam=0.3)

        assert len(profileset) == 1
        profile = list(profileset.profiles.values())[0][0]
        assert len(profile.quartets) == 1

        quartet = next(iter(profile.quartets))
        assert quartet.split == Split({'A', 'B'}, {'C', 'D'})

    def test_delta_threshold_split_vs_cycle(self) -> None:
        matrix = np.array([
            [0.0, 0.5, 0.5, 0.5],
            [0.5, 0.0, 0.5, 0.5],
            [0.5, 0.5, 0.0, 0.5],
            [0.5, 0.5, 0.5, 0.0]
        ])
        dm = DistanceMatrix(matrix, labels=['A', 'B', 'C', 'D'])

        profileset = delta_heuristic(dm, lam=0.3)
        assert len(profileset) == 1
        profile = list(profileset.profiles.values())[0][0]
        assert len(profile.quartets) == 1

    def test_weighted_profiles(self) -> None:
        matrix = np.array([
            [0.0, 0.1, 0.9, 0.9],
            [0.1, 0.0, 0.9, 0.9],
            [0.9, 0.9, 0.0, 0.1],
            [0.9, 0.9, 0.1, 0.0]
        ])
        dm = DistanceMatrix(matrix, labels=['A', 'B', 'C', 'D'])

        profileset_weighted = delta_heuristic(dm, lam=0.3, weight=True)
        profile_weighted = list(profileset_weighted.profiles.values())[0][0]
        quartet = next(iter(profile_weighted.quartets))
        weight = profile_weighted.get_weight(quartet)
        assert weight is not None
        assert weight > 0

        profileset_unweighted = delta_heuristic(dm, lam=0.3, weight=False)
        profile_unweighted = list(profileset_unweighted.profiles.values())[0][0]
        quartet_unweighted = next(iter(profile_unweighted.quartets))
        weight_unweighted = profile_unweighted.get_weight(quartet_unweighted)
        assert weight_unweighted == 1.0

    def test_reticulation_leaf_assignment(self) -> None:
        matrix = np.array([
            [0.0, 0.2, 0.3, 0.4],
            [0.2, 0.0, 0.4, 0.3],
            [0.3, 0.4, 0.0, 0.2],
            [0.4, 0.3, 0.2, 0.0]
        ])
        dm = DistanceMatrix(matrix, labels=['A', 'B', 'C', 'D'])

        profileset = delta_heuristic(dm, lam=0.01, weight=False)

        has_cycle = any(len(profile.quartets) == 2 for profile, _ in profileset.profiles.values())

        if has_cycle:
            cycle_profile = next(
                (profile for profile, _ in profileset.profiles.values() if len(profile.quartets) == 2),
                None
            )
            assert cycle_profile is not None
            assert cycle_profile.reticulation_leaf is not None
            assert cycle_profile.reticulation_leaf in cycle_profile.taxa

    def test_lambda_validation(self) -> None:
        matrix = np.array([
            [0.0, 0.1, 0.9, 0.9],
            [0.1, 0.0, 0.9, 0.9],
            [0.9, 0.9, 0.0, 0.1],
            [0.9, 0.9, 0.1, 0.0]
        ])
        dm = DistanceMatrix(matrix, labels=['A', 'B', 'C', 'D'])

        with pytest.raises(ValueError, match="Lambda must be in"):
            delta_heuristic(dm, lam=-0.1)

        with pytest.raises(ValueError, match="Lambda must be in"):
            delta_heuristic(dm, lam=1.5)

    def test_five_taxa_creates_five_profiles(self) -> None:
        matrix = np.array([
            [0.0, 0.1, 0.2, 0.3, 0.4],
            [0.1, 0.0, 0.2, 0.3, 0.4],
            [0.2, 0.2, 0.0, 0.3, 0.4],
            [0.3, 0.3, 0.3, 0.0, 0.4],
            [0.4, 0.4, 0.4, 0.4, 0.0]
        ])
        dm = DistanceMatrix(matrix, labels=['A', 'B', 'C', 'D', 'E'])
        profileset = delta_heuristic(dm, lam=0.3)

        assert len(profileset) == 5

    def test_all_profiles_have_four_taxa(self) -> None:
        matrix = np.array([
            [0.0, 0.1, 0.2, 0.3, 0.4],
            [0.1, 0.0, 0.2, 0.3, 0.4],
            [0.2, 0.2, 0.0, 0.3, 0.4],
            [0.3, 0.3, 0.3, 0.0, 0.4],
            [0.4, 0.4, 0.4, 0.4, 0.0]
        ])
        dm = DistanceMatrix(matrix, labels=['A', 'B', 'C', 'D', 'E'])
        profileset = delta_heuristic(dm, lam=0.3)

        for profile, _ in profileset.profiles.values():
            assert len(profile.taxa) == 4

    def test_profiles_are_valid_sqprofiles(self) -> None:
        matrix = np.array([
            [0.0, 0.1, 0.9, 0.9],
            [0.1, 0.0, 0.9, 0.9],
            [0.9, 0.9, 0.0, 0.1],
            [0.9, 0.9, 0.1, 0.0]
        ])
        dm = DistanceMatrix(matrix, labels=['A', 'B', 'C', 'D'])
        profileset = delta_heuristic(dm, lam=0.3)

        for profile, _ in profileset.profiles.values():
            assert len(profile.quartets) in (1, 2)
            for quartet in profile.quartets:
                assert quartet.is_resolved()


class TestDeltaHeuristicToNetwork:
    """Test cases for delta heuristic to network pipeline."""

    def test_delta_heuristic_to_network_no_errors(self) -> None:
        matrix = np.array([
            [0.0, 0.1, 0.9, 0.9, 0.9],
            [0.1, 0.0, 0.9, 0.9, 0.9],
            [0.9, 0.9, 0.0, 0.1, 0.9],
            [0.9, 0.9, 0.1, 0.0, 0.9],
            [0.9, 0.9, 0.9, 0.9, 0.0]
        ])
        dm = DistanceMatrix(matrix, labels=['A', 'B', 'C', 'D', 'E'])

        profileset = delta_heuristic(dm, lam=0.3)
        assert len(profileset) > 0

        network = squirrel(profileset)

        assert isinstance(network, SemiDirectedPhyNetwork)
        assert network.number_of_nodes() > 0
        assert len(network.taxa) == 5
        network.validate()

    def test_delta_heuristic_to_network_multiple_taxa(self) -> None:
        matrix = np.array([
            [0.0, 0.1, 0.9, 0.9, 0.9, 0.9],
            [0.1, 0.0, 0.9, 0.9, 0.9, 0.9],
            [0.9, 0.9, 0.0, 0.1, 0.9, 0.9],
            [0.9, 0.9, 0.1, 0.0, 0.9, 0.9],
            [0.9, 0.9, 0.9, 0.9, 0.0, 0.1],
            [0.9, 0.9, 0.9, 0.9, 0.1, 0.0]
        ])
        dm = DistanceMatrix(matrix, labels=['A', 'B', 'C', 'D', 'E', 'F'])

        profileset = delta_heuristic(dm, lam=0.3)
        assert len(profileset) > 0

        network = squirrel(profileset)
        assert isinstance(network, SemiDirectedPhyNetwork)
        assert len(network.taxa) == 6
        network.validate()


class TestNetworkRoundTrip:
    """Test cases for network -> profileset -> network round-trip."""

    def test_network_round_trip_simple_tree(self) -> None:
        from physquirrel import sqprofileset_from_network
        from tests.fixtures.sd_networks import SDTREE_MEDIUM_BINARY
        from phylozoo.core.network.sdnetwork.features import blobs, k_blobs
        from phylozoo.core.network.sdnetwork.classifications import has_parallel_edges
        from phylozoo.core.network.sdnetwork.isomorphism import is_isomorphic

        network = SDTREE_MEDIUM_BINARY

        assert not has_parallel_edges(network)

        profileset = sqprofileset_from_network(network)
        assert len(profileset) > 0

        reconstructed = squirrel(profileset)

        assert isinstance(reconstructed, SemiDirectedPhyNetwork)
        reconstructed.validate()

        assert is_isomorphic(network, reconstructed)
