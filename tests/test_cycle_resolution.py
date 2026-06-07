"""
Tests for cycle resolution module.
"""

import pytest

from phylozoo.core.quartet.base import Quartet
from phylozoo.core.split.base import Split
from phylozoo.core.primitives.partition import Partition
from physquirrel import SqQuartetProfile, SqQuartetProfileSet
from physquirrel.algorithms.cycle_resolution import _qprofiles_to_hybrid_ranking


class TestQProfilesToHybridRanking:
    """Tests for the _qprofiles_to_hybrid_ranking function."""

    def test_four_set_partition_with_cycle(self) -> None:
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        q2 = Quartet(Split({'A', 'C'}, {'B', 'D'}))
        profileset = SqQuartetProfileSet(
            profiles=[(q1, 1.0), (q2, 1.0)],
            taxa={'A', 'B', 'C', 'D'}
        )
        profile = profileset.get_profile(frozenset({'A', 'B', 'C', 'D'}))
        assert profile is not None
        assert len(profile.quartets) == 2

        sq_profile = SqQuartetProfile(
            {q1: 0.5, q2: 0.5},
            reticulation_leaf='A'
        )
        profileset_with_ret = SqQuartetProfileSet(
            profiles=[sq_profile],
            taxa={'A', 'B', 'C', 'D'}
        )

        partition = Partition([{'A'}, {'B'}, {'C'}, {'D'}])
        order = _qprofiles_to_hybrid_ranking(profileset_with_ret, partition, weights=True)

        assert isinstance(order, list)
        assert len(order) == 4
        assert all(isinstance(s, frozenset) for s in order)
        assert order[0] == frozenset({'A'})

    def test_four_set_partition_without_cycle(self) -> None:
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        q2 = Quartet(Split({'A', 'C'}, {'B', 'D'}))
        profileset = SqQuartetProfileSet(
            profiles=[q1, q2],
            taxa={'A', 'B', 'C', 'D'}
        )

        partition = Partition([{'A'}, {'B'}, {'C'}, {'D'}])
        order = _qprofiles_to_hybrid_ranking(profileset, partition, weights=True)

        assert isinstance(order, list)
        assert len(order) == 4
        assert all(isinstance(s, frozenset) for s in order)

    def test_larger_partition_with_cycles(self) -> None:
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        q2 = Quartet(Split({'A', 'C'}, {'B', 'D'}))
        q3 = Quartet(Split({'A', 'B'}, {'E', 'F'}))
        q4 = Quartet(Split({'A', 'E'}, {'B', 'F'}))

        q5 = Quartet(Split({'C', 'D'}, {'E', 'F'}))
        q6 = Quartet(Split({'C', 'E'}, {'D', 'F'}))
        cycle_profile = SqQuartetProfile(
            {q5: 0.5, q6: 0.5},
            reticulation_leaf='C'
        )

        profileset = SqQuartetProfileSet(
            profiles=[q1, q2, q3, q4, cycle_profile],
            taxa={'A', 'B', 'C', 'D', 'E', 'F'}
        )

        partition = Partition([{'A'}, {'B'}, {'C'}, {'D'}, {'E'}, {'F'}])
        order = _qprofiles_to_hybrid_ranking(profileset, partition, weights=True)

        assert isinstance(order, list)
        assert len(order) == 6
        assert all(isinstance(s, frozenset) for s in order)

    def test_weights_parameter(self) -> None:
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        q2 = Quartet(Split({'A', 'C'}, {'B', 'D'}))
        sq_profile = SqQuartetProfile(
            {q1: 0.5, q2: 0.5},
            reticulation_leaf='A'
        )
        profileset = SqQuartetProfileSet(
            profiles=[sq_profile],
            taxa={'A', 'B', 'C', 'D'}
        )

        partition = Partition([{'A'}, {'B'}, {'C'}, {'D'}])

        order_with_weights = _qprofiles_to_hybrid_ranking(profileset, partition, weights=True)
        order_without_weights = _qprofiles_to_hybrid_ranking(profileset, partition, weights=False)

        assert order_with_weights[0] == frozenset({'A'})
        assert order_without_weights[0] == frozenset({'A'})

    def test_empty_profileset(self) -> None:
        profileset = SqQuartetProfileSet(profiles=[], taxa={'A', 'B', 'C', 'D'})
        partition = Partition([{'A'}, {'B'}, {'C'}, {'D'}])
        order = _qprofiles_to_hybrid_ranking(profileset, partition, weights=True)

        assert isinstance(order, list)
        assert len(order) == 4

    def test_partition_with_missing_profiles(self) -> None:
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        profileset = SqQuartetProfileSet(
            profiles=[q1],
            taxa={'A', 'B', 'C', 'D', 'E', 'F'}
        )

        partition = Partition([{'A'}, {'B'}, {'C'}, {'D'}, {'E'}, {'F'}])
        order = _qprofiles_to_hybrid_ranking(profileset, partition, weights=True)

        assert isinstance(order, list)
        assert len(order) == 6


class TestInsertCycle:
    """Tests for the _insert_cycle function."""

    def test_basic_insert_cycle(self) -> None:
        from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
        from phylozoo.core.primitives.circular_ordering import CircularSetOrdering
        from physquirrel.algorithms.cycle_resolution import _insert_cycle

        net = SemiDirectedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
        )

        ordering = CircularSetOrdering([{'A'}, {'B'}, {'C'}])
        ranking = [frozenset({'A'}), frozenset({'B'}), frozenset({'C'})]
        new_net = _insert_cycle(net, 3, ordering, reticulation_ranking=ranking)

        assert isinstance(new_net, SemiDirectedPhyNetwork)
        assert 3 not in new_net._graph.nodes
        assert len(new_net._graph.nodes) > len(net._graph.nodes)
        assert new_net.number_of_edges() > net.number_of_edges()

    def test_insert_cycle_with_four_leaves(self) -> None:
        from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
        from phylozoo.core.primitives.circular_ordering import CircularSetOrdering
        from physquirrel.algorithms.cycle_resolution import _insert_cycle

        net = SemiDirectedPhyNetwork(
            undirected_edges=[(5, 1), (5, 2), (5, 3), (5, 4)],
            nodes=[
                (1, {'label': 'A'}), (2, {'label': 'B'}),
                (3, {'label': 'C'}), (4, {'label': 'D'})
            ]
        )

        ordering = CircularSetOrdering([{'A'}, {'B'}, {'C'}, {'D'}])
        ranking = [frozenset({'A'}), frozenset({'B'}), frozenset({'C'}), frozenset({'D'})]
        new_net = _insert_cycle(net, 5, ordering, reticulation_ranking=ranking)

        assert isinstance(new_net, SemiDirectedPhyNetwork)
        original_leaf_nodes = {1, 2, 3, 4}
        new_nodes = set(new_net._graph.nodes)
        cycle_nodes = new_nodes - original_leaf_nodes
        assert len(cycle_nodes) == 4

    def test_insert_cycle_validation(self) -> None:
        from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
        from phylozoo.core.primitives.circular_ordering import CircularSetOrdering
        from physquirrel.algorithms.cycle_resolution import _insert_cycle

        net = SemiDirectedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
        )
        ordering = CircularSetOrdering([{'A'}, {'B'}, {'C'}])
        ranking = [frozenset({'A'})]

        with pytest.raises(ValueError, match="is not a cut-vertex"):
            _insert_cycle(net, 1, ordering, reticulation_ranking=ranking)

        wrong_ordering = CircularSetOrdering([{'A'}, {'B'}])
        with pytest.raises(ValueError, match="does not match partition"):
            _insert_cycle(net, 3, wrong_ordering, reticulation_ranking=ranking)

        bad_ranking = [frozenset({'X'})]
        with pytest.raises(ValueError, match="not part of the circular set ordering"):
            _insert_cycle(net, 3, ordering, reticulation_ranking=bad_ranking)

    def test_insert_cycle_preserves_taxa(self) -> None:
        from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
        from phylozoo.core.primitives.circular_ordering import CircularSetOrdering
        from physquirrel.algorithms.cycle_resolution import _insert_cycle

        net = SemiDirectedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
        )

        ordering = CircularSetOrdering([{'A'}, {'B'}, {'C'}])
        ranking = [frozenset({'A'}), frozenset({'B'}), frozenset({'C'})]
        new_net = _insert_cycle(net, 3, ordering, reticulation_ranking=ranking)

        assert new_net.taxa == net.taxa
        assert 'A' in new_net.taxa
        assert 'B' in new_net.taxa
        assert 'C' in new_net.taxa

    def test_insert_cycle_returns_semidirected_with_ranking(self) -> None:
        from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
        from phylozoo.core.primitives.circular_ordering import CircularSetOrdering
        from physquirrel.algorithms.cycle_resolution import _insert_cycle

        net = SemiDirectedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
        )

        ordering = CircularSetOrdering([{'A'}, {'B'}, {'C'}])
        ranking = [frozenset({'A'}), frozenset({'B'}), frozenset({'C'})]
        new_net = _insert_cycle(net, 3, ordering, reticulation_ranking=ranking)

        assert isinstance(new_net, SemiDirectedPhyNetwork)
        assert 'A' in new_net.taxa

    def test_insert_cycle_without_ranking_returns_mixed(self) -> None:
        from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork, MixedPhyNetwork
        from phylozoo.core.primitives.circular_ordering import CircularSetOrdering
        from physquirrel.algorithms.cycle_resolution import _insert_cycle

        net = SemiDirectedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
        )

        ordering = CircularSetOrdering([{'A'}, {'B'}, {'C'}])
        new_net = _insert_cycle(net, 3, ordering, reticulation_ranking=None)
        assert isinstance(new_net, MixedPhyNetwork)
        assert not isinstance(new_net, SemiDirectedPhyNetwork)
