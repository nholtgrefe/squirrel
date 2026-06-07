"""
Tests for SqQuartetProfile class.
"""

import pytest

from phylozoo.core.quartet.base import Quartet
from phylozoo.core.split.base import Split
from physquirrel import SqQuartetProfile


class TestSqQuartetProfileInit:
    """Tests for SqQuartetProfile initialization."""

    def test_init_single_quartet(self) -> None:
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        profile = SqQuartetProfile([q1])

        assert profile.taxa == frozenset({'A', 'B', 'C', 'D'})
        assert len(profile) == 1
        assert profile.reticulation_leaf is None
        assert profile.circular_ordering is None
        assert profile.get_weight(q1) == 1.0

    def test_init_single_quartet_with_weight(self) -> None:
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        profile = SqQuartetProfile({q1: 1.0})

        assert len(profile) == 1
        assert profile.get_weight(q1) == 1.0
        assert profile.reticulation_leaf is None

    def test_init_two_quartets_without_reticulation_leaf(self) -> None:
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        q2 = Quartet(Split({'A', 'C'}, {'B', 'D'}))
        profile = SqQuartetProfile([q1, q2])

        assert len(profile) == 2
        assert profile.reticulation_leaf is None
        assert profile.circular_ordering is not None
        assert profile.get_weight(q1) == pytest.approx(0.5)
        assert profile.get_weight(q2) == pytest.approx(0.5)

    def test_init_two_quartets_with_reticulation_leaf(self) -> None:
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        q2 = Quartet(Split({'A', 'C'}, {'B', 'D'}))
        profile = SqQuartetProfile([q1, q2], reticulation_leaf='A')

        assert len(profile) == 2
        assert profile.reticulation_leaf == 'A'
        assert profile.circular_ordering is not None
        assert 'A' in profile.circular_ordering.elements

    def test_init_two_quartets_with_weights(self) -> None:
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        q2 = Quartet(Split({'A', 'C'}, {'B', 'D'}))
        profile = SqQuartetProfile({q1: 0.7, q2: 0.3}, reticulation_leaf='B')

        assert len(profile) == 2
        assert profile.reticulation_leaf == 'B'
        assert profile.get_weight(q1) == 0.7
        assert profile.get_weight(q2) == 0.3

    def test_init_unresolved_quartet_error(self) -> None:
        star_q = Quartet({'A', 'B', 'C', 'D'})

        with pytest.raises(ValueError, match="SqQuartetProfile cannot contain unresolved quartets"):
            SqQuartetProfile([star_q])

        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        with pytest.raises(ValueError, match="SqQuartetProfile cannot contain unresolved quartets"):
            SqQuartetProfile([q1, star_q])

    def test_init_zero_quartets_error(self) -> None:
        with pytest.raises(ValueError, match="QuartetProfile must have at least one quartet"):
            SqQuartetProfile([])

        with pytest.raises(ValueError, match="QuartetProfile must have at least one quartet"):
            SqQuartetProfile({})

    def test_init_three_quartets_error(self) -> None:
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        q2 = Quartet(Split({'A', 'C'}, {'B', 'D'}))
        q3 = Quartet(Split({'A', 'D'}, {'B', 'C'}))

        with pytest.raises(ValueError, match="SqQuartetProfile must contain exactly 1 or 2 resolved quartets"):
            SqQuartetProfile([q1, q2, q3])

    def test_init_reticulation_leaf_with_one_quartet_error(self) -> None:
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))

        with pytest.raises(ValueError, match="reticulation_leaf can only be provided when profile has 2 quartets"):
            SqQuartetProfile([q1], reticulation_leaf='A')

    def test_init_reticulation_leaf_not_in_taxa_error(self) -> None:
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        q2 = Quartet(Split({'A', 'C'}, {'B', 'D'}))

        with pytest.raises(ValueError, match="Reticulation leaf 'X' must be one of the taxa"):
            SqQuartetProfile([q1, q2], reticulation_leaf='X')

    def test_init_different_taxa_error(self) -> None:
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        q2 = Quartet(Split({'E', 'F'}, {'G', 'H'}))

        with pytest.raises(ValueError, match="All quartets must have the same taxa"):
            SqQuartetProfile([q1, q2])

    def test_init_non_positive_weight_error(self) -> None:
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))

        with pytest.raises(ValueError, match="Weight must be positive"):
            SqQuartetProfile({q1: 0.0})

        with pytest.raises(ValueError, match="Weight must be positive"):
            SqQuartetProfile({q1: -0.5})


class TestSqQuartetProfileProperties:
    """Tests for SqQuartetProfile properties."""

    def test_reticulation_leaf_property(self) -> None:
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        q2 = Quartet(Split({'A', 'C'}, {'B', 'D'}))

        profile1 = SqQuartetProfile([q1])
        assert profile1.reticulation_leaf is None

        profile2 = SqQuartetProfile([q1, q2], reticulation_leaf='A')
        assert profile2.reticulation_leaf == 'A'

        profile3 = SqQuartetProfile([q1, q2])
        assert profile3.reticulation_leaf is None

    def test_circular_ordering_property(self) -> None:
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        q2 = Quartet(Split({'A', 'C'}, {'B', 'D'}))

        profile1 = SqQuartetProfile([q1])
        assert profile1.circular_ordering is None

        profile2 = SqQuartetProfile([q1, q2])
        assert profile2.circular_ordering is not None
        assert profile2.circular_ordering.elements == frozenset({'A', 'B', 'C', 'D'})

    def test_inherited_properties(self) -> None:
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        q2 = Quartet(Split({'A', 'C'}, {'B', 'D'}))
        profile = SqQuartetProfile({q1: 0.6, q2: 0.4})

        assert profile.taxa == frozenset({'A', 'B', 'C', 'D'})
        assert profile.is_resolved() is True
        assert profile.is_trivial() is False
        assert q1 in profile
        assert q2 in profile


class TestSqQuartetProfileEdgeCases:
    """Tests for edge cases and special behaviors."""

    def test_all_three_possible_quartets_error(self) -> None:
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        q2 = Quartet(Split({'A', 'C'}, {'B', 'D'}))
        q3 = Quartet(Split({'A', 'D'}, {'B', 'C'}))

        with pytest.raises(ValueError, match="SqQuartetProfile must contain exactly 1 or 2 resolved quartets"):
            SqQuartetProfile([q1, q2, q3])

    def test_reticulation_leaf_all_four_taxa(self) -> None:
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        q2 = Quartet(Split({'A', 'C'}, {'B', 'D'}))

        for leaf in ['A', 'B', 'C', 'D']:
            profile = SqQuartetProfile([q1, q2], reticulation_leaf=leaf)
            assert profile.reticulation_leaf == leaf

    def test_circular_ordering_consistency(self) -> None:
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        q2 = Quartet(Split({'A', 'C'}, {'B', 'D'}))
        profile = SqQuartetProfile([q1, q2])

        ordering = profile.circular_ordering
        assert ordering is not None
        assert len(ordering.elements) == 4
        assert ordering.elements == frozenset({'A', 'B', 'C', 'D'})
        assert len(ordering.order) == 4

    def test_repr(self) -> None:
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        profile1 = SqQuartetProfile([q1])
        repr1 = repr(profile1)
        assert 'SqQuartetProfile' in repr1
        assert 'reticulation_leaf=None' in repr1

        q2 = Quartet(Split({'A', 'C'}, {'B', 'D'}))
        profile2 = SqQuartetProfile([q1, q2], reticulation_leaf='A')
        repr2 = repr(profile2)
        assert 'SqQuartetProfile' in repr2
        assert "reticulation_leaf='A'" in repr2
