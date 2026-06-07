"""
Tests for the tstar_tree module.
"""

import pytest

from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
from phylozoo.core.network.sdnetwork.classifications import is_tree
from phylozoo.core.quartet.base import Quartet
from phylozoo.core.quartet.qprofile import QuartetProfile
from phylozoo.core.quartet.qprofileset import QuartetProfileSet
from phylozoo.core.split.base import Split
from phylozoo.core.split.classifications import is_tree_compatible
from phylozoo.core.split.splitsystem import SplitSystem
from physquirrel import bstar, tstar_tree


class TestBstar:
    """Tests for the bstar function."""

    def test_bstar_simple_four_taxa(self) -> None:
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        q2 = Quartet(Split({'A', 'C'}, {'B', 'D'}))
        profileset = QuartetProfileSet(profiles=[q1, q2])

        result = bstar(profileset)

        assert isinstance(result, SplitSystem)
        assert len(result) >= 4
        assert result.elements == {'A', 'B', 'C', 'D'}

    def test_bstar_fewer_than_four_taxa_raises_error(self) -> None:
        profileset_empty = QuartetProfileSet()

        with pytest.raises(ValueError, match="B\\* algorithm requires at least 4 taxa"):
            bstar(profileset_empty)

    def test_bstar_only_trivial_profiles(self) -> None:
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        profile1 = QuartetProfile({q1: 1.0})

        q2 = Quartet(Split({'A', 'C'}, {'B', 'E'}))
        q3 = Quartet(Split({'A', 'E'}, {'B', 'C'}))
        profile2 = QuartetProfile({q2: 0.6, q3: 0.4})

        profileset = QuartetProfileSet(profiles=[profile1, profile2])

        result = bstar(profileset)

        assert isinstance(result, SplitSystem)
        assert result.elements == {'A', 'B', 'C', 'D', 'E'}

    def test_bstar_five_taxa(self) -> None:
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        q2 = Quartet(Split({'A', 'B'}, {'C', 'E'}))
        q3 = Quartet(Split({'A', 'C'}, {'D', 'E'}))

        profileset = QuartetProfileSet(profiles=[q1, q2, q3])

        result = bstar(profileset)

        assert isinstance(result, SplitSystem)
        assert result.elements == {'A', 'B', 'C', 'D', 'E'}
        assert len(result) >= 5

    def test_bstar_larger_set(self) -> None:
        quartets = [
            Quartet(Split({'A', 'B'}, {'C', 'D'})),
            Quartet(Split({'A', 'B'}, {'E', 'F'})),
            Quartet(Split({'C', 'D'}, {'E', 'F'})),
            Quartet(Split({'A', 'C'}, {'B', 'D'})),
        ]

        profileset = QuartetProfileSet(profiles=quartets)

        result = bstar(profileset)

        assert isinstance(result, SplitSystem)
        assert result.elements == {'A', 'B', 'C', 'D', 'E', 'F'}
        assert len(result) >= 6

    def test_bstar_returns_compatible_splits(self) -> None:
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        q2 = Quartet(Split({'A', 'C'}, {'B', 'D'}))
        profileset = QuartetProfileSet(profiles=[q1, q2])

        result = bstar(profileset)

        from phylozoo.core.split.algorithms import is_tree_compatible

        assert is_tree_compatible(result), "B* splits are not pairwise compatible"


class TestTstarTree:
    """Tests for the tstar_tree function."""

    def test_tstar_tree_simple_four_taxa(self) -> None:
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        q2 = Quartet(Split({'A', 'C'}, {'B', 'D'}))
        profileset = QuartetProfileSet(profiles=[q1, q2])

        tree = tstar_tree(profileset)

        assert isinstance(tree, SemiDirectedPhyNetwork)
        assert is_tree(tree)
        assert tree.taxa == {'A', 'B', 'C', 'D'}

    def test_tstar_tree_fewer_than_four_taxa_raises_error(self) -> None:
        profileset = QuartetProfileSet()

        with pytest.raises(ValueError, match="B\\* algorithm requires at least 4 taxa"):
            tstar_tree(profileset)

    def test_tstar_tree_five_taxa(self) -> None:
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        q2 = Quartet(Split({'A', 'B'}, {'C', 'E'}))
        q3 = Quartet(Split({'A', 'C'}, {'D', 'E'}))

        profileset = QuartetProfileSet(profiles=[q1, q2, q3])

        tree = tstar_tree(profileset)

        assert isinstance(tree, SemiDirectedPhyNetwork)
        assert is_tree(tree)
        assert tree.taxa == {'A', 'B', 'C', 'D', 'E'}

    def test_tstar_tree_uses_bstar_splits(self) -> None:
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        q2 = Quartet(Split({'A', 'C'}, {'B', 'D'}))
        profileset = QuartetProfileSet(profiles=[q1, q2])

        bstar_splits = bstar(profileset)
        tree = tstar_tree(profileset)

        from phylozoo.core.network.sdnetwork.derivations import induced_splits

        tree_splits = induced_splits(tree)

        for split in bstar_splits.splits:
            assert split in tree_splits.splits, \
                f"Split {split} from bstar is not in the tree splits"

    def test_tstar_tree_consistent_with_bstar(self) -> None:
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        q2 = Quartet(Split({'A', 'B'}, {'C', 'E'}))
        q3 = Quartet(Split({'A', 'B'}, {'D', 'E'}))
        q4 = Quartet(Split({'A', 'C'}, {'D', 'E'}))
        q5 = Quartet(Split({'B', 'C'}, {'D', 'E'}))

        profileset = QuartetProfileSet(profiles=[q1, q2, q3, q4, q5])

        bstar_splits = bstar(profileset)
        tree = tstar_tree(profileset)

        assert is_tree(tree)
        assert tree.taxa == profileset.taxa

        from phylozoo.core.network.sdnetwork.derivations import induced_splits

        tree_splits = induced_splits(tree)

        for split in bstar_splits.splits:
            assert split in tree_splits.splits, \
                f"Split {split} from bstar is not in the tree splits"

    def test_tstar_tree_empty_profileset(self) -> None:
        profileset = QuartetProfileSet()

        with pytest.raises(ValueError, match="B\\* algorithm requires at least 4 taxa"):
            tstar_tree(profileset)

    def test_tstar_tree_single_quartet(self) -> None:
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        profileset = QuartetProfileSet(profiles=[q1])

        tree = tstar_tree(profileset)

        assert isinstance(tree, SemiDirectedPhyNetwork)
        assert is_tree(tree)
        assert tree.taxa == {'A', 'B', 'C', 'D'}
        assert tree.number_of_nodes() > 0
        assert tree.number_of_edges() > 0
