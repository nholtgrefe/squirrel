"""
Tests for qsimilarity module.
"""

import pytest

from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
from physquirrel.algorithms.qsimilarity import _circular_orders_from_cycles
from tests.fixtures.sd_networks import LEVEL_1_SDNETWORK_SINGLE_HYBRID


class TestCircularOrdersFromCycles:
    """Tests for _circular_orders_from_cycles function."""

    def test_non_binary_network_raises_error(self) -> None:
        network = SemiDirectedPhyNetwork(
            undirected_edges=[
                (3, 1), (3, 2), (3, 4), (3, 5)
            ],
            nodes=[
                (1, {'label': 'A'}),
                (2, {'label': 'B'}),
                (4, {'label': 'C'}),
                (5, {'label': 'D'}),
            ]
        )

        with pytest.raises(ValueError, match="Network must be binary"):
            list(_circular_orders_from_cycles(network))

    def test_single_hybrid_network_with_cycle(self) -> None:
        network = LEVEL_1_SDNETWORK_SINGLE_HYBRID

        results = list(_circular_orders_from_cycles(network))

        assert len(results) == 1

        set_ordering, ret_set = results[0]

        assert set_ordering is not None
        assert len(set_ordering) == 3

        assert ret_set is not None

        assert ret_set == {'A', 'B'} or ret_set == frozenset({'A', 'B'})

        all_taxa_in_ordering = set()
        for taxa_set in set_ordering:
            all_taxa_in_ordering.update(taxa_set)

        assert all_taxa_in_ordering == {'A', 'B', 'C', 'D'}

        assert ret_set in set_ordering

    def test_tree_network_no_cycles(self) -> None:
        network = SemiDirectedPhyNetwork(
            undirected_edges=[
                (3, 1),
                (3, 2),
                (3, 4),
            ],
            nodes=[
                (1, {'label': 'A'}),
                (2, {'label': 'B'}),
                (4, {'label': 'C'}),
            ]
        )

        results = list(_circular_orders_from_cycles(network))

        assert len(results) == 0

    def test_network_with_cycle_yields_tuples(self) -> None:
        network = LEVEL_1_SDNETWORK_SINGLE_HYBRID

        results = list(_circular_orders_from_cycles(network))

        assert len(results) >= 1

        for set_ordering, ret_set in results:
            assert set_ordering is not None
            assert isinstance(ret_set, (frozenset, type(None)))


class TestSqProfilesetFromNetwork:
    """Tests for sqprofileset_from_network function."""

    def test_two_hybrids_network(self) -> None:
        from phylozoo.core.quartet.base import Quartet
        from phylozoo.core.split.base import Split
        from physquirrel import SqQuartetProfile, SqQuartetProfileSet, sqprofileset_from_network
        from tests.fixtures.sd_networks import LEVEL_1_SDNETWORK_TWO_HYBRIDS_SEPARATE

        network = LEVEL_1_SDNETWORK_TWO_HYBRIDS_SEPARATE

        result_profileset = sqprofileset_from_network(network)

        profile_1 = result_profileset.get_profile(frozenset({'A', 'B', 'C', 'D'}))
        profile_2 = result_profileset.get_profile(frozenset({'A', 'B', 'C', 'E'}))
        profile_3 = result_profileset.get_profile(frozenset({'A', 'B', 'D', 'E'}))
        profile_4 = result_profileset.get_profile(frozenset({'A', 'C', 'D', 'E'}))
        profile_5 = result_profileset.get_profile(frozenset({'B', 'C', 'D', 'E'}))

        assert profile_1 is not None
        assert profile_2 is not None
        assert profile_3 is not None
        assert profile_4 is not None
        assert profile_5 is not None

        assert profile_4.taxa == frozenset({'A', 'C', 'D', 'E'})
        assert len(profile_4) >= 1
        assert len(profile_4) <= 2
        if len(profile_4) == 2:
            assert profile_4.reticulation_leaf == 'A'

        assert profile_5.taxa == frozenset({'B', 'C', 'D', 'E'})
        assert len(profile_5) >= 1
        assert len(profile_5) <= 2
        if len(profile_5) == 2:
            assert profile_5.reticulation_leaf == 'B'

        expected_profileset = SqQuartetProfileSet(profiles=[
            profile_1, profile_2, profile_3, profile_4, profile_5,
        ])

        assert len(result_profileset) == len(expected_profileset)
        assert len(result_profileset) == 5

        assert result_profileset.taxa == expected_profileset.taxa

        for taxa_set in [frozenset({'A', 'B', 'C', 'D'}),
                         frozenset({'A', 'B', 'C', 'E'}),
                         frozenset({'A', 'B', 'D', 'E'}),
                         frozenset({'A', 'C', 'D', 'E'}),
                         frozenset({'B', 'C', 'D', 'E'})]:
            result_profile = result_profileset.get_profile(taxa_set)
            expected_profile = expected_profileset.get_profile(taxa_set)
            assert result_profile is not None
            assert expected_profile is not None
            assert result_profile._quartets == expected_profile._quartets
            assert result_profile.reticulation_leaf == expected_profile.reticulation_leaf
            assert result_profileset.get_profile_weight(taxa_set) == expected_profileset.get_profile_weight(taxa_set)
