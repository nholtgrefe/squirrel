"""
Tests for the qjoining module.
"""

import pytest

from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
from phylozoo.core.network.sdnetwork.classifications import is_tree
from phylozoo.core.quartet.base import Quartet
from phylozoo.core.quartet.qprofileset import QuartetProfileSet
from phylozoo.core.split.base import Split
from physquirrel import quartet_joining, adapted_quartet_joining


class TestQuartetJoining:
    """Tests for the quartet_joining function."""

    def test_quartet_joining_simple_four_taxa(self) -> None:
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        profileset = QuartetProfileSet(profiles=[q1])

        tree = quartet_joining(profileset)

        assert isinstance(tree, SemiDirectedPhyNetwork)
        assert is_tree(tree)
        assert tree.taxa == {'A', 'B', 'C', 'D'}

    def test_quartet_joining_five_taxa(self) -> None:
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        q2 = Quartet(Split({'A', 'B'}, {'C', 'E'}))
        q3 = Quartet(Split({'A', 'B'}, {'D', 'E'}))
        profileset = QuartetProfileSet(profiles=[q1, q2, q3])

        tree = quartet_joining(profileset)

        assert isinstance(tree, SemiDirectedPhyNetwork)
        assert is_tree(tree)
        assert tree.taxa == {'A', 'B', 'C', 'D', 'E'}

    def test_quartet_joining_fewer_than_four_taxa_raises_error(self) -> None:
        profileset_empty = QuartetProfileSet()

        with pytest.raises(ValueError, match="Quartet joining requires at least 4 taxa"):
            quartet_joining(profileset_empty)

        profileset_no_profiles = QuartetProfileSet(profiles=[])

        with pytest.raises(ValueError, match="Quartet joining requires at least 4 taxa"):
            quartet_joining(profileset_no_profiles)

    def test_quartet_joining_produces_binary_tree(self) -> None:
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        q2 = Quartet(Split({'A', 'B'}, {'C', 'E'}))
        q3 = Quartet(Split({'A', 'B'}, {'D', 'E'}))
        profileset = QuartetProfileSet(profiles=[q1, q2, q3])

        tree = quartet_joining(profileset)

        for node in tree.internal_nodes:
            assert tree.degree(node) == 3

    def test_quartet_joining_larger_set(self) -> None:
        quartets = [
            Quartet(Split({'A', 'B'}, {'C', 'D'})),
            Quartet(Split({'A', 'B'}, {'C', 'E'})),
            Quartet(Split({'A', 'B'}, {'C', 'F'})),
            Quartet(Split({'A', 'B'}, {'D', 'E'})),
            Quartet(Split({'A', 'B'}, {'D', 'F'})),
            Quartet(Split({'A', 'B'}, {'E', 'F'})),
        ]
        profileset = QuartetProfileSet(profiles=quartets)

        tree = quartet_joining(profileset)

        assert isinstance(tree, SemiDirectedPhyNetwork)
        assert is_tree(tree)
        assert tree.taxa == {'A', 'B', 'C', 'D', 'E', 'F'}
        assert tree.number_of_edges() == tree.number_of_nodes() - 1


class TestAdaptedQuartetJoining:
    """Tests for the adapted_quartet_joining function."""

    def test_adapted_quartet_joining_with_star_tree(self) -> None:
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        profileset = QuartetProfileSet(profiles=[q1])

        center = '_center'
        starting = SemiDirectedPhyNetwork(
            undirected_edges=[(center, 'A'), (center, 'B'), (center, 'C'), (center, 'D')],
            nodes=[('A', {'label': 'A'}), ('B', {'label': 'B'}), ('C', {'label': 'C'}), ('D', {'label': 'D'})]
        )

        tree = adapted_quartet_joining(profileset, starting)

        assert isinstance(tree, SemiDirectedPhyNetwork)
        assert is_tree(tree)
        assert tree.taxa == {'A', 'B', 'C', 'D'}

    def test_adapted_quartet_joining_with_non_star_tree(self) -> None:
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        q2 = Quartet(Split({'A', 'B'}, {'C', 'E'}))
        profileset = QuartetProfileSet(profiles=[q1, q2])

        starting = SemiDirectedPhyNetwork(
            undirected_edges=[
                (3, 'A'), (3, 'B'), (3, 4),
                (4, 5), (4, 'E'),
                (5, 'C'), (5, 'D')
            ],
            nodes=[('A', {'label': 'A'}), ('B', {'label': 'B'}), ('C', {'label': 'C'}),
                   ('D', {'label': 'D'}), ('E', {'label': 'E'})]
        )

        tree = adapted_quartet_joining(profileset, starting)

        assert isinstance(tree, SemiDirectedPhyNetwork)
        assert is_tree(tree)
        assert tree.taxa == {'A', 'B', 'C', 'D', 'E'}

    def test_adapted_quartet_joining_taxa_mismatch_raises_error(self) -> None:
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        profileset = QuartetProfileSet(profiles=[q1])

        starting = SemiDirectedPhyNetwork(
            undirected_edges=[('_center', 'A'), ('_center', 'B'), ('_center', 'C'), ('_center', 'X')],
            nodes=[('A', {'label': 'A'}), ('B', {'label': 'B'}), ('C', {'label': 'C'}), ('X', {'label': 'X'})]
        )

        with pytest.raises(ValueError, match="Starting tree taxa must match profile set taxa"):
            adapted_quartet_joining(profileset, starting)

    def test_adapted_quartet_joining_produces_binary_tree(self) -> None:
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        q2 = Quartet(Split({'A', 'B'}, {'C', 'E'}))
        profileset = QuartetProfileSet(profiles=[q1, q2])

        center = '_center'
        starting = SemiDirectedPhyNetwork(
            undirected_edges=[(center, 'A'), (center, 'B'), (center, 'C'), (center, 'D'), (center, 'E')],
            nodes=[('A', {'label': 'A'}), ('B', {'label': 'B'}), ('C', {'label': 'C'}),
                   ('D', {'label': 'D'}), ('E', {'label': 'E'})]
        )

        tree = adapted_quartet_joining(profileset, starting)

        for node in tree.internal_nodes:
            assert tree.degree(node) == 3


class TestQuartetJoiningIntegration:
    """Integration tests comparing quartet_joining and adapted_quartet_joining."""

    def test_quartet_joining_equivalent_to_adapted_with_star_tree(self) -> None:
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        q2 = Quartet(Split({'A', 'B'}, {'C', 'E'}))
        profileset = QuartetProfileSet(profiles=[q1, q2])

        tree1 = quartet_joining(profileset)

        center = '_center'
        starting = SemiDirectedPhyNetwork(
            undirected_edges=[(center, 'A'), (center, 'B'), (center, 'C'), (center, 'D'), (center, 'E')],
            nodes=[('A', {'label': 'A'}), ('B', {'label': 'B'}), ('C', {'label': 'C'}),
                   ('D', {'label': 'D'}), ('E', {'label': 'E'})]
        )
        tree2 = adapted_quartet_joining(profileset, starting)

        assert is_tree(tree1)
        assert is_tree(tree2)
        assert tree1.taxa == tree2.taxa
        assert tree1.taxa == {'A', 'B', 'C', 'D', 'E'}


class TestQuartetJoiningEdgeCases:
    """Edge case tests for quartet joining functions."""

    def test_quartet_joining_with_single_quartet(self) -> None:
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        profileset = QuartetProfileSet(profiles=[q1])

        tree = quartet_joining(profileset)

        assert is_tree(tree)
        assert tree.taxa == {'A', 'B', 'C', 'D'}
        assert len(tree.internal_nodes) == 2

    def test_quartet_joining_with_multiple_quartets_same_taxa(self) -> None:
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        q2 = Quartet(Split({'A', 'C'}, {'B', 'D'}))
        q3 = Quartet(Split({'A', 'D'}, {'B', 'C'}))
        profileset = QuartetProfileSet(profiles=[q1, q2, q3])

        tree = quartet_joining(profileset)

        assert is_tree(tree)
        assert tree.taxa == {'A', 'B', 'C', 'D'}

    def test_quartet_joining_already_binary_tree(self) -> None:
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        profileset = QuartetProfileSet(profiles=[q1])

        tree = quartet_joining(profileset)

        for node in tree.internal_nodes:
            assert tree.degree(node) == 3
