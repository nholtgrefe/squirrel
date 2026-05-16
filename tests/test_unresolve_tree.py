"""
Tests for the unresolve_tree module.
"""

import pytest

from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
from phylozoo.core.network.sdnetwork.classifications import is_tree
from phylozoo.core.quartet.base import Quartet
from phylozoo.core.quartet.qprofileset import QuartetProfileSet
from phylozoo.core.split.base import Split
from physquirrel import unresolve_tree, split_support


class TestSplitSupport:
    """Tests for the split_support function."""

    def test_split_support_single_matching_quartet(self) -> None:
        split = Split({'A', 'B'}, {'C', 'D'})
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        profileset = QuartetProfileSet(profiles=[q1])
        support = split_support(profileset, split)
        assert support == 1.0

    def test_split_support_no_matching_quartets(self) -> None:
        split = Split({'A', 'B'}, {'C', 'D'})
        q1 = Quartet(Split({'A', 'C'}, {'B', 'D'}))
        profileset = QuartetProfileSet(profiles=[q1])
        support = split_support(profileset, split)
        assert support == 0.0

    def test_split_support_partial_match(self) -> None:
        split = Split({'A', 'B'}, {'C', 'D', 'E'})
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        q2 = Quartet(Split({'A', 'B'}, {'C', 'E'}))
        q3 = Quartet(Split({'A', 'C'}, {'B', 'D'}))
        profileset = QuartetProfileSet(profiles=[q1, q2, q3], taxa={'A', 'B', 'C', 'D', 'E'})
        support = split_support(profileset, split)
        assert support == 1.0

    def test_split_support_trivial_split_raises_error(self) -> None:
        split = Split({'A'}, {'B', 'C', 'D'})
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        profileset = QuartetProfileSet(profiles=[q1])
        with pytest.raises(ValueError, match="Split must be non-trivial"):
            split_support(profileset, split)

    def test_split_support_taxa_mismatch_raises_error(self) -> None:
        split = Split({'A', 'B'}, {'C', 'D', 'E'})
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        profileset = QuartetProfileSet(profiles=[q1])
        with pytest.raises(ValueError, match="Split taxa must match profile set taxa"):
            split_support(profileset, split)

    def test_split_support_non_trivial_profile_ignored(self) -> None:
        split = Split({'A', 'B'}, {'C', 'D'})
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        q2 = Quartet(Split({'A', 'C'}, {'B', 'D'}))
        from phylozoo.core.quartet.qprofile import QuartetProfile
        profile = QuartetProfile([q1, q2])
        profileset = QuartetProfileSet(profiles=[profile])
        support = split_support(profileset, split)
        assert support == 0.0

    def test_split_support_unresolved_quartet_ignored(self) -> None:
        split = Split({'A', 'B'}, {'C', 'D', 'E'})
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        q_star = Quartet({'A', 'B', 'C', 'E'})
        profileset = QuartetProfileSet(profiles=[q1, q_star], taxa={'A', 'B', 'C', 'D', 'E'})
        support = split_support(profileset, split)
        assert support == 1.0


class TestUnresolveTree:
    """Tests for the unresolve_tree generator function."""

    def test_unresolve_tree_single_split(self) -> None:
        tree = SemiDirectedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (4, 3), (4, 5), (4, 6)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (5, {'label': 'C'}), (6, {'label': 'D'})]
        )
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        profileset = QuartetProfileSet(profiles=[q1])
        trees = list(unresolve_tree(tree, profileset))
        assert len(trees) >= 2
        assert all(is_tree(t) for t in trees)
        assert trees[0].number_of_nodes() == tree.number_of_nodes()
        assert trees[-1].number_of_nodes() < trees[0].number_of_nodes()

    def test_unresolve_tree_multiple_splits(self) -> None:
        tree = SemiDirectedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (4, 3), (4, 5), (4, 6)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (5, {'label': 'C'}), (6, {'label': 'D'})]
        )
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        q2 = Quartet(Split({'A', 'C'}, {'B', 'D'}))
        profileset = QuartetProfileSet(profiles=[q1, q2])
        trees = list(unresolve_tree(tree, profileset))
        assert len(trees) >= 2
        assert all(is_tree(t) for t in trees)

    def test_unresolve_tree_no_non_trivial_splits(self) -> None:
        tree = SemiDirectedPhyNetwork(
            undirected_edges=[(0, 1), (0, 2), (0, 3)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'})]
        )
        profileset = QuartetProfileSet(profiles=[], taxa={'A', 'B', 'C'})
        trees = list(unresolve_tree(tree, profileset))
        assert len(trees) == 1
        assert trees[0].number_of_nodes() == tree.number_of_nodes()

    def test_unresolve_tree_invalid_tree_raises_error(self) -> None:
        network = SemiDirectedPhyNetwork(
            directed_edges=[(5, 4), (6, 4)],
            undirected_edges=[(4, 1), (5, 6), (5, 8), (6, 9)],
            nodes=[(1, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'})]
        )
        profileset = QuartetProfileSet(profiles=[], taxa={'A', 'B', 'C'})
        with pytest.raises(ValueError, match="Tree must be a valid tree"):
            list(unresolve_tree(network, profileset))

    def test_unresolve_tree_taxa_mismatch_raises_error(self) -> None:
        tree = SemiDirectedPhyNetwork(
            undirected_edges=[(1, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        profileset = QuartetProfileSet(profiles=[q1])
        with pytest.raises(ValueError, match="Tree taxa must match profile set taxa"):
            list(unresolve_tree(tree, profileset))

    def test_unresolve_tree_all_trees_valid(self) -> None:
        tree = SemiDirectedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (4, 3), (4, 5), (4, 6)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (5, {'label': 'C'}), (6, {'label': 'D'})]
        )
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        q2 = Quartet(Split({'A', 'C'}, {'B', 'D'}))
        profileset = QuartetProfileSet(profiles=[q1, q2])
        trees = list(unresolve_tree(tree, profileset))
        for t in trees:
            assert is_tree(t)
            assert t.taxa == tree.taxa

    def test_unresolve_tree_decreasing_nodes(self) -> None:
        tree = SemiDirectedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (4, 3), (4, 5), (4, 6)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (5, {'label': 'C'}), (6, {'label': 'D'})]
        )
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        profileset = QuartetProfileSet(profiles=[q1])
        trees = list(unresolve_tree(tree, profileset))
        if len(trees) > 1:
            for i in range(1, len(trees)):
                assert trees[i].number_of_nodes() <= trees[i - 1].number_of_nodes()
