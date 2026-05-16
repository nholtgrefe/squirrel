"""
Tests for SqQuartetProfileSet class.
"""

import pytest

from phylozoo.core.quartet.base import Quartet
from phylozoo.core.quartet.qprofile import QuartetProfile
from phylozoo.core.split.base import Split
from physquirrel import SqQuartetProfile, SqQuartetProfileSet
from phylozoo.utils.exceptions import PhyloZooParseError, PhyloZooFormatError


class TestSqQuartetProfileSetInit:
    """Tests for SqQuartetProfileSet initialization."""

    def test_init_empty(self) -> None:
        profileset = SqQuartetProfileSet()
        assert len(profileset) == 0
        assert len(profileset.taxa) == 0

    def test_init_from_sq_profiles(self) -> None:
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({5, 6}, {7, 8}))
        sq_profile1 = SqQuartetProfile([q1])
        sq_profile2 = SqQuartetProfile([q2])

        profileset = SqQuartetProfileSet(profiles=[sq_profile1, sq_profile2])

        assert len(profileset) == 2
        assert profileset.get_profile_weight(frozenset({1, 2, 3, 4})) == 1.0
        assert profileset.get_profile_weight(frozenset({5, 6, 7, 8})) == 1.0

    def test_init_from_sq_profiles_with_weights(self) -> None:
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({5, 6}, {7, 8}))
        sq_profile1 = SqQuartetProfile([q1])
        sq_profile2 = SqQuartetProfile([q2])

        profileset = SqQuartetProfileSet(profiles=[(sq_profile1, 2.0), (sq_profile2, 1.5)])

        assert len(profileset) == 2
        assert profileset.get_profile_weight(frozenset({1, 2, 3, 4})) == 2.0
        assert profileset.get_profile_weight(frozenset({5, 6, 7, 8})) == 1.5

    def test_init_from_quartets_auto_conversion(self) -> None:
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({5, 6}, {7, 8}))

        profileset = SqQuartetProfileSet(profiles=[q1, q2])

        assert len(profileset) == 2
        assert profileset.get_profile_weight(frozenset({1, 2, 3, 4})) == 1.0
        assert profileset.get_profile_weight(frozenset({5, 6, 7, 8})) == 1.0

        profile1 = profileset.get_profile(frozenset({1, 2, 3, 4}))
        assert isinstance(profile1, SqQuartetProfile)
        assert len(profile1) == 1

    def test_init_from_quartets_with_weights(self) -> None:
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({5, 6}, {7, 8}))

        profileset = SqQuartetProfileSet(profiles=[(q1, 0.8), (q2, 1.2)])

        assert len(profileset) == 2
        assert profileset.get_profile_weight(frozenset({1, 2, 3, 4})) == 0.8
        assert profileset.get_profile_weight(frozenset({5, 6, 7, 8})) == 1.2

    def test_init_from_quartets_grouped_by_taxa(self) -> None:
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({1, 3}, {2, 4}))

        profileset = SqQuartetProfileSet(profiles=[q1, q2])

        assert len(profileset) == 1
        profile = profileset.get_profile(frozenset({1, 2, 3, 4}))
        assert profile is not None
        assert len(profile) == 2
        assert q1 in profile
        assert q2 in profile
        assert profileset.get_profile_weight(frozenset({1, 2, 3, 4})) == 2.0

    def test_init_from_two_quartets_same_taxa_no_reticulation_leaf(self) -> None:
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({1, 3}, {2, 4}))

        profileset = SqQuartetProfileSet(profiles=[q1, q2])

        assert len(profileset) == 1
        profile = profileset.get_profile(frozenset({1, 2, 3, 4}))
        assert profile is not None
        assert len(profile) == 2
        assert q1 in profile
        assert q2 in profile
        assert profile.reticulation_leaf is None

    def test_init_mixed_sq_profiles_and_quartets(self) -> None:
        q1 = Quartet(Split({1, 2}, {3, 4}))
        sq_profile1 = SqQuartetProfile([q1])
        q2 = Quartet(Split({5, 6}, {7, 8}))

        profileset = SqQuartetProfileSet(profiles=[sq_profile1, q2])
        assert len(profileset) == 2

    def test_init_non_sq_profile_error(self) -> None:
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({1, 3}, {2, 4}))
        regular_profile = QuartetProfile([q1, q2])

        with pytest.raises(ValueError, match="Expected SqQuartetProfile or Quartet"):
            SqQuartetProfileSet(profiles=[regular_profile])

    def test_init_non_positive_weight_error(self) -> None:
        q1 = Quartet(Split({1, 2}, {3, 4}))
        sq_profile1 = SqQuartetProfile([q1])

        with pytest.raises(ValueError, match="Profile weight must be positive"):
            SqQuartetProfileSet(profiles=[(sq_profile1, 0.0)])

        with pytest.raises(ValueError, match="Profile weight must be positive"):
            SqQuartetProfileSet(profiles=[(sq_profile1, -0.5)])

    def test_init_duplicate_quartet_raises(self) -> None:
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({1, 3}, {2, 4}))

        with pytest.raises(ValueError, match="appears multiple times in the input"):
            SqQuartetProfileSet(profiles=[q1, q2, q1])

    def test_init_three_quartets_same_taxa_raises(self) -> None:
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({1, 3}, {2, 4}))
        q3 = Quartet(Split({1, 4}, {2, 3}))

        with pytest.raises(ValueError, match="SqQuartetProfile must contain exactly 1 or 2 resolved quartets"):
            SqQuartetProfileSet(profiles=[q1, q2, q3])


class TestSqQuartetProfileSetProperties:
    """Tests for SqQuartetProfileSet properties."""

    def test_profiles_property(self) -> None:
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({5, 6}, {7, 8}))
        sq_profile1 = SqQuartetProfile([q1])
        sq_profile2 = SqQuartetProfile([q2])

        profileset = SqQuartetProfileSet(profiles=[sq_profile1, sq_profile2])

        for taxa, (profile, weight) in profileset.profiles.items():
            assert isinstance(profile, SqQuartetProfile)
            assert weight > 0

    def test_get_profile_returns_sq_profile(self) -> None:
        q1 = Quartet(Split({1, 2}, {3, 4}))
        sq_profile1 = SqQuartetProfile([q1])

        profileset = SqQuartetProfileSet(profiles=[sq_profile1])

        profile = profileset.get_profile(frozenset({1, 2, 3, 4}))
        assert isinstance(profile, SqQuartetProfile)
        assert len(profile) == 1

        assert profileset.get_profile(frozenset({5, 6, 7, 8})) is None

    def test_inherited_properties(self) -> None:
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({5, 6}, {7, 8}))
        sq_profile1 = SqQuartetProfile([q1])
        sq_profile2 = SqQuartetProfile([q2])

        profileset = SqQuartetProfileSet(profiles=[sq_profile1, sq_profile2])

        assert profileset.taxa == frozenset({1, 2, 3, 4, 5, 6, 7, 8})
        assert len(profileset) == 2
        assert profileset.has_profile(frozenset({1, 2, 3, 4})) is True
        assert profileset.has_profile(frozenset({9, 10, 11, 12})) is False


class TestSqQuartetProfileSetEdgeCases:
    """Tests for edge cases and special behaviors."""

    def test_profiles_with_reticulation_leaf(self) -> None:
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({1, 3}, {2, 4}))
        sq_profile = SqQuartetProfile([q1, q2], reticulation_leaf=1)

        profileset = SqQuartetProfileSet(profiles=[sq_profile])

        profile = profileset.get_profile(frozenset({1, 2, 3, 4}))
        assert profile is not None
        assert profile.reticulation_leaf == 1

    def test_auto_conversion_preserves_quartet_weights(self) -> None:
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({1, 3}, {2, 4}))

        profileset = SqQuartetProfileSet(profiles=[(q1, 0.6), (q2, 0.4)])

        profile = profileset.get_profile(frozenset({1, 2, 3, 4}))
        assert profile is not None
        assert profileset.get_profile_weight(frozenset({1, 2, 3, 4})) == 1.0
        assert profile.get_weight(q1) == 0.6
        assert profile.get_weight(q2) == 0.4

    def test_taxa_parameter(self) -> None:
        q1 = Quartet(Split({1, 2}, {3, 4}))
        sq_profile1 = SqQuartetProfile([q1])

        profileset = SqQuartetProfileSet(
            profiles=[sq_profile1],
            taxa=frozenset({1, 2, 3, 4, 5, 6, 7, 8})
        )

        assert profileset.taxa == frozenset({1, 2, 3, 4, 5, 6, 7, 8})
        assert len(profileset) == 1

    def test_taxa_parameter_not_superset_error(self) -> None:
        q1 = Quartet(Split({1, 2}, {3, 4}))
        sq_profile1 = SqQuartetProfile([q1])

        with pytest.raises(ValueError, match="Provided taxa must be a superset"):
            SqQuartetProfileSet(
                profiles=[sq_profile1],
                taxa=frozenset({1, 2, 3})
            )


class TestSqQuartetProfileSetIO:
    """Tests for SqQuartetProfileSet I/O operations."""

    def test_to_string_psq_format(self) -> None:
        q1 = Quartet(Split({1, 2}, {3, 4}))
        sq_profile1 = SqQuartetProfile([q1])
        profileset = SqQuartetProfileSet(profiles=[sq_profile1])

        psq_string = profileset.to_string(format='psq')

        assert isinstance(psq_string, str)
        assert 'profiles' in psq_string
        assert 'taxa' in psq_string

    def test_from_string_psq_format(self) -> None:
        q1 = Quartet(Split({1, 2}, {3, 4}))
        sq_profile1 = SqQuartetProfile([q1])
        original_profileset = SqQuartetProfileSet(profiles=[sq_profile1])

        psq_string = original_profileset.to_string(format='psq')
        reconstructed_profileset = SqQuartetProfileSet.from_string(psq_string, format='psq')

        assert len(reconstructed_profileset) == len(original_profileset)
        assert reconstructed_profileset.taxa == original_profileset.taxa
        assert reconstructed_profileset.has_profile(frozenset({1, 2, 3, 4}))

    def test_to_string_from_string_roundtrip(self) -> None:
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({1, 3}, {2, 4}))
        sq_profile = SqQuartetProfile([q1, q2], reticulation_leaf=1)
        original_profileset = SqQuartetProfileSet(profiles=[sq_profile])

        psq_string = original_profileset.to_string(format='psq')
        reconstructed_profileset = SqQuartetProfileSet.from_string(psq_string, format='psq')

        assert len(reconstructed_profileset) == len(original_profileset)
        original_profile = original_profileset.get_profile(frozenset({1, 2, 3, 4}))
        reconstructed_profile = reconstructed_profileset.get_profile(frozenset({1, 2, 3, 4}))

        assert original_profile is not None
        assert reconstructed_profile is not None
        assert len(original_profile) == len(reconstructed_profile)
        assert original_profile.reticulation_leaf == reconstructed_profile.reticulation_leaf

        for q in original_profile:
            assert q in reconstructed_profile
            assert original_profile.get_weight(q) == reconstructed_profile.get_weight(q)

    def test_to_string_from_string_with_weights(self) -> None:
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({5, 6}, {7, 8}))
        sq_profile1 = SqQuartetProfile([q1])
        sq_profile2 = SqQuartetProfile([q2])
        original_profileset = SqQuartetProfileSet(
            profiles=[(sq_profile1, 2.0), (sq_profile2, 1.5)]
        )

        psq_string = original_profileset.to_string(format='psq')
        reconstructed_profileset = SqQuartetProfileSet.from_string(psq_string, format='psq')

        assert len(reconstructed_profileset) == 2
        assert reconstructed_profileset.get_profile_weight(frozenset({1, 2, 3, 4})) == 2.0
        assert reconstructed_profileset.get_profile_weight(frozenset({5, 6, 7, 8})) == 1.5

    def test_from_string_invalid_json(self) -> None:
        with pytest.raises(PhyloZooParseError):
            SqQuartetProfileSet.from_string("not valid json {", format='psq')

    def test_from_string_missing_keys(self) -> None:
        with pytest.raises(PhyloZooFormatError):
            SqQuartetProfileSet.from_string('{"profiles": []}', format='psq')

    def test_from_string_invalid_quartet_type(self) -> None:
        invalid_data = '''{
            "taxa": [1, 2, 3, 4],
            "profiles": [{
                "taxa": [1, 2, 3, 4],
                "quartets": [{
                    "quartet": {"type": "invalid"},
                    "weight": 1.0
                }],
                "reticulation_leaf": null,
                "profile_weight": 1.0
            }]
        }'''
        with pytest.raises(PhyloZooFormatError):
            SqQuartetProfileSet.from_string(invalid_data, format='psq')
