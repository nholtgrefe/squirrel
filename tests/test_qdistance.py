"""
Tests for the qdistance module (quartet_distance_with_partition).
"""

import pytest

from phylozoo import Quartet, Split
from phylozoo.core.distance import DistanceMatrix
from phylozoo.core.primitives.partition import Partition
from phylozoo.core.quartet.qprofile import QuartetProfile
from phylozoo.core.quartet.qprofileset import QuartetProfileSet
from physquirrel import SqQuartetProfile, SqQuartetProfileSet
from physquirrel.algorithms.qdistance import quartet_distance_with_partition


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _split_profileset(*taxa_splits, extra_profiles=None):
    """Build a SqQuartetProfileSet from (taxa_frozenset, Split, weight) triples.

    Pass extra_profiles as a list of (SqQuartetProfile, weight) tuples for
    profiles that need non-default construction (e.g. cycle profiles).
    """
    profiles = []
    for taxa_set, split, weight in taxa_splits:
        q = Quartet(split)
        profiles.append((SqQuartetProfile([q]), weight))
    if extra_profiles:
        profiles.extend(extra_profiles)
    all_taxa = frozenset().union(*(t for t, _, __ in taxa_splits))
    if extra_profiles:
        for p, _ in extra_profiles:
            all_taxa = all_taxa | frozenset().union(*p.quartets)
    return SqQuartetProfileSet(profiles=profiles, taxa=all_taxa)


def _dense_five_taxon_profileset(weights=None):
    """Return a dense SqQuartetProfileSet on {A,B,C,D,E} with AB|CD-style splits."""
    four_sets = [
        frozenset({'A', 'B', 'C', 'D'}),
        frozenset({'A', 'B', 'C', 'E'}),
        frozenset({'A', 'B', 'D', 'E'}),
        frozenset({'A', 'C', 'D', 'E'}),
        frozenset({'B', 'C', 'D', 'E'}),
    ]
    splits = [
        Split({'A', 'B'}, {'C', 'D'}),
        Split({'A', 'B'}, {'C', 'E'}),
        Split({'A', 'B'}, {'D', 'E'}),
        Split({'A', 'C'}, {'D', 'E'}),
        Split({'B', 'C'}, {'D', 'E'}),
    ]
    if weights is None:
        weights = [1.0] * 5
    profiles = [
        (SqQuartetProfile([Quartet(s)]), w)
        for s, w in zip(splits, weights)
    ]
    return SqQuartetProfileSet(
        profiles=profiles,
        taxa=frozenset({'A', 'B', 'C', 'D', 'E'}),
    )


# ---------------------------------------------------------------------------
# TestQuartetDistanceWithPartition
# ---------------------------------------------------------------------------

class TestQuartetDistanceWithPartition:
    """Tests for the unweighted behaviour of quartet_distance_with_partition."""

    def test_returns_distance_matrix(self) -> None:
        q = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        ps = SqQuartetProfileSet(profiles=[SqQuartetProfile([q])], taxa=frozenset({'A', 'B', 'C', 'D'}))
        dm = quartet_distance_with_partition(ps, Partition([{'A'}, {'B'}, {'C'}, {'D'}]))
        assert isinstance(dm, DistanceMatrix)

    def test_size_equals_number_of_partition_sets(self) -> None:
        q = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        ps = SqQuartetProfileSet(profiles=[SqQuartetProfile([q])], taxa=frozenset({'A', 'B', 'C', 'D'}))
        dm = quartet_distance_with_partition(ps, Partition([{'A'}, {'B'}, {'C'}, {'D'}]))
        assert len(dm) == 4

    def test_diagonal_is_zero(self) -> None:
        q = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        ps = SqQuartetProfileSet(profiles=[SqQuartetProfile([q])], taxa=frozenset({'A', 'B', 'C', 'D'}))
        dm = quartet_distance_with_partition(ps, Partition([{'A'}, {'B'}, {'C'}, {'D'}]))
        for label in dm.labels:
            assert dm.get_distance(label, label) == pytest.approx(0.0)

    def test_symmetric(self) -> None:
        q = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        ps = SqQuartetProfileSet(profiles=[SqQuartetProfile([q])], taxa=frozenset({'A', 'B', 'C', 'D'}))
        dm = quartet_distance_with_partition(ps, Partition([{'A'}, {'B'}, {'C'}, {'D'}]))
        for s1 in dm.labels:
            for s2 in dm.labels:
                assert dm.get_distance(s1, s2) == pytest.approx(dm.get_distance(s2, s1))

    def test_constant_added_correctly(self) -> None:
        # With rho=(0,0,0,0) only the constant 2n-4 remains in off-diagonal entries.
        q = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        ps = SqQuartetProfileSet(profiles=[SqQuartetProfile([q])], taxa=frozenset({'A', 'B', 'C', 'D'}))
        dm = quartet_distance_with_partition(ps, Partition([{'A'}, {'B'}, {'C'}, {'D'}]),
                                             rho=(0.0, 0.0, 0.0, 0.0))
        # n=4 sets → constant = 2*4-4 = 4
        for s1 in dm.labels:
            for s2 in dm.labels:
                expected = 0.0 if s1 == s2 else 4.0
                assert dm.get_distance(s1, s2) == pytest.approx(expected)

    def test_singleton_partition_same_side(self) -> None:
        # A,B on same side of split → delta = 2*rho_c = 1.0; total = 1.0 + 4 = 5.0
        q = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        ps = SqQuartetProfileSet(profiles=[SqQuartetProfile([q])], taxa=frozenset({'A', 'B', 'C', 'D'}))
        dm = quartet_distance_with_partition(ps, Partition([{'A'}, {'B'}, {'C'}, {'D'}]))
        assert dm.get_distance(frozenset({'A'}), frozenset({'B'})) == pytest.approx(5.0)

    def test_singleton_partition_opposite_sides(self) -> None:
        # A,C on opposite sides → delta = 2*rho_s = 2.0; total = 2.0 + 4 = 6.0
        q = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        ps = SqQuartetProfileSet(profiles=[SqQuartetProfile([q])], taxa=frozenset({'A', 'B', 'C', 'D'}))
        dm = quartet_distance_with_partition(ps, Partition([{'A'}, {'B'}, {'C'}, {'D'}]))
        assert dm.get_distance(frozenset({'A'}), frozenset({'C'})) == pytest.approx(6.0)

    def test_three_set_partition_five_taxa(self) -> None:
        ps = _dense_five_taxon_profileset()
        dm = quartet_distance_with_partition(ps, Partition([{'A', 'B'}, {'C'}, {'D', 'E'}]))
        assert len(dm) == 3
        # constant for 3 sets: 2*3-4 = 2; all off-diagonal ≥ 2
        for s1 in dm.labels:
            for s2 in dm.labels:
                if s1 != s2:
                    assert dm.get_distance(s1, s2) >= 2.0

    def test_five_singleton_partition(self) -> None:
        ps = _dense_five_taxon_profileset()
        dm = quartet_distance_with_partition(
            ps, Partition([{'A'}, {'B'}, {'C'}, {'D'}, {'E'}])
        )
        assert len(dm) == 5
        # constant for 5 sets: 2*5-4 = 6; all off-diagonal ≥ 6
        for s1 in dm.labels:
            for s2 in dm.labels:
                if s1 != s2:
                    assert dm.get_distance(s1, s2) >= 6.0

    def test_cycle_profile_distances(self) -> None:
        # Circular order A,B,C,D → adjacent (A,B) closer than opposite (A,C)
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        q2 = Quartet(Split({'A', 'D'}, {'B', 'C'}))
        profile = SqQuartetProfile([q1, q2])
        ps = SqQuartetProfileSet(profiles=[profile], taxa=frozenset({'A', 'B', 'C', 'D'}))
        dm = quartet_distance_with_partition(ps, Partition([{'A'}, {'B'}, {'C'}, {'D'}]))
        d_adjacent = dm.get_distance(frozenset({'A'}), frozenset({'B'}))
        d_opposite = dm.get_distance(frozenset({'A'}), frozenset({'C'}))
        assert d_adjacent < d_opposite

    def test_nanuq_rho(self) -> None:
        # rho=(0,1,0.5,1) → same-side delta=0, constant=4 → d(A,B)=4.0
        q = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        ps = SqQuartetProfileSet(profiles=[SqQuartetProfile([q])], taxa=frozenset({'A', 'B', 'C', 'D'}))
        dm = quartet_distance_with_partition(ps, Partition([{'A'}, {'B'}, {'C'}, {'D'}]),
                                             rho=(0.0, 1.0, 0.5, 1.0))
        assert dm.get_distance(frozenset({'A'}), frozenset({'B'})) == pytest.approx(4.0)
        assert dm.get_distance(frozenset({'A'}), frozenset({'C'})) == pytest.approx(6.0)

    def test_averaging_across_multiple_representative_quartets(self) -> None:
        # Partition {A,B}|{C}|{D}|{E}: 2 representative quartets per set pair.
        # Both {A,C,D,E} and {B,C,D,E} have AB-like structure → consistent averages.
        ps = _dense_five_taxon_profileset()
        dm = quartet_distance_with_partition(
            ps, Partition([{'A', 'B'}, {'C'}, {'D'}, {'E'}])
        )
        assert isinstance(dm, DistanceMatrix)
        assert len(dm) == 4

    # --- Validation errors ---

    def test_invalid_rho_length_raises(self) -> None:
        q = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        ps = SqQuartetProfileSet(profiles=[SqQuartetProfile([q])], taxa=frozenset({'A', 'B', 'C', 'D'}))
        with pytest.raises(Exception, match="4 elements"):
            quartet_distance_with_partition(ps, Partition([{'A'}, {'B'}, {'C'}, {'D'}]),
                                            rho=(0.5, 1.0, 0.5))

    @pytest.mark.parametrize('rho', [
        (0.5, 1.0, 1.5, 1.0),  # rho_a > rho_o
        (1.5, 1.0, 0.5, 1.0),  # rho_c > rho_s
    ])
    def test_invalid_rho_constraints_raises(self, rho) -> None:
        q = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        ps = SqQuartetProfileSet(profiles=[SqQuartetProfile([q])], taxa=frozenset({'A', 'B', 'C', 'D'}))
        with pytest.raises(Exception, match="satisfy"):
            quartet_distance_with_partition(ps, Partition([{'A'}, {'B'}, {'C'}, {'D'}]), rho=rho)

    def test_partition_taxon_mismatch_raises(self) -> None:
        q = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        ps = SqQuartetProfileSet(profiles=[SqQuartetProfile([q])], taxa=frozenset({'A', 'B', 'C', 'D'}))
        with pytest.raises(Exception, match="match"):
            quartet_distance_with_partition(ps, Partition([{'A'}, {'B'}, {'C'}, {'X'}]))

    def test_non_dense_profileset_raises(self) -> None:
        q = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        ps = SqQuartetProfileSet(profiles=[SqQuartetProfile([q])],
                                 taxa=frozenset({'A', 'B', 'C', 'D', 'E'}))
        with pytest.raises(Exception, match="dense"):
            quartet_distance_with_partition(ps, Partition([{'A'}, {'B'}, {'C'}, {'D'}, {'E'}]))

    def test_unresolved_quartet_in_profile_raises(self) -> None:
        # Use base QuartetProfileSet — SqQuartetProfileSet rejects unresolved quartets itself.
        star = Quartet(frozenset({'A', 'B', 'C', 'D'}))
        ps = QuartetProfileSet(profiles=[QuartetProfile([star])],
                               taxa=frozenset({'A', 'B', 'C', 'D'}))
        with pytest.raises(Exception):
            quartet_distance_with_partition(ps, Partition([{'A'}, {'B'}, {'C'}, {'D'}]))


# ---------------------------------------------------------------------------
# TestWeightedDistance
# ---------------------------------------------------------------------------

class TestWeightedDistance:
    """Tests for the weighted_distance parameter."""

    def test_uniform_weights_matches_unweighted(self) -> None:
        """When all profile weights are 1.0, weighted_distance gives the same result."""
        ps = _dense_five_taxon_profileset(weights=[1.0, 1.0, 1.0, 1.0, 1.0])
        partition = Partition([{'A'}, {'B'}, {'C'}, {'D'}, {'E'}])
        dm_uw = quartet_distance_with_partition(ps, partition, weighted_distance=False)
        dm_w = quartet_distance_with_partition(ps, partition, weighted_distance=True)
        for s1 in dm_uw.labels:
            for s2 in dm_uw.labels:
                assert dm_w.get_distance(s1, s2) == pytest.approx(dm_uw.get_distance(s1, s2))

    def test_low_weight_profile_contributes_less_than_high_weight(self) -> None:
        """In a two-representative scenario, the lower-weight profile is down-weighted.

        Partition {A,B}|{C}|{D}|{E}: representatives are {A,C,D,E} (weight=1.0, same-side
        A,C → delta=1.0) and {B,C,D,E} (weight=w_b, opposite-side B,C → delta=2.0).

        Weighted average = (1.0 + 2.0×w_b) / (1.0 + w_b).
        As w_b → 0, average → 1.0 (all same-side signal).
        As w_b → 1, average → 1.5 (equal mix, same as unweighted).
        """
        taxa = frozenset({'A', 'B', 'C', 'D', 'E'})
        partition = Partition([{'A', 'B'}, {'C'}, {'D'}, {'E'}])
        set_ab, set_c = frozenset({'A', 'B'}), frozenset({'C'})

        def make_ps(w_b: float) -> SqQuartetProfileSet:
            return SqQuartetProfileSet(
                profiles=[
                    (SqQuartetProfile([Quartet(Split({'A', 'B'}, {'C', 'D'}))]), 1.0),
                    (SqQuartetProfile([Quartet(Split({'A', 'B'}, {'C', 'E'}))]), 1.0),
                    (SqQuartetProfile([Quartet(Split({'A', 'B'}, {'D', 'E'}))]), 1.0),
                    (SqQuartetProfile([Quartet(Split({'A', 'C'}, {'D', 'E'}))]), 1.0),
                    (SqQuartetProfile([Quartet(Split({'B', 'D'}, {'C', 'E'}))]), w_b),
                ],
                taxa=taxa,
            )

        # With w_b=0.01 (near zero), weighted avg ≈ 4 + 1.0 = 5.0 (same-side dominates)
        dm_low = quartet_distance_with_partition(make_ps(0.01), partition, weighted_distance=True)
        # With w_b=1.0 (equal), weighted avg = unweighted avg = 4 + 1.5 = 5.5
        dm_equal = quartet_distance_with_partition(make_ps(1.0), partition, weighted_distance=True)

        assert dm_low.get_distance(set_ab, set_c) < dm_equal.get_distance(set_ab, set_c)
        assert dm_low.get_distance(set_ab, set_c) == pytest.approx(5.0, abs=0.05)

    def test_weighted_average_between_component_deltas(self) -> None:
        """Weighted distance is always between the min and max individual rho distances."""
        taxa = frozenset({'A', 'B', 'C', 'D', 'E'})
        profiles = [
            (SqQuartetProfile([Quartet(Split({'A', 'B'}, {'C', 'D'}))]), 1.0),
            (SqQuartetProfile([Quartet(Split({'A', 'B'}, {'C', 'E'}))]), 1.0),
            (SqQuartetProfile([Quartet(Split({'A', 'B'}, {'D', 'E'}))]), 1.0),
            # same-side (delta=1.0), weight=0.7
            (SqQuartetProfile([Quartet(Split({'A', 'C'}, {'D', 'E'}))]), 0.7),
            # opposite-side (delta=2.0), weight=0.3
            (SqQuartetProfile([Quartet(Split({'B', 'D'}, {'C', 'E'}))]), 0.3),
        ]
        ps = SqQuartetProfileSet(profiles=profiles, taxa=taxa)
        partition = Partition([{'A', 'B'}, {'C'}, {'D'}, {'E'}])
        set_ab, set_c = frozenset({'A', 'B'}), frozenset({'C'})

        dm = quartet_distance_with_partition(ps, partition, weighted_distance=True)
        d = dm.get_distance(set_ab, set_c)

        # Distance must be between 4+1.0=5.0 (all same-side) and 4+2.0=6.0 (all opposite)
        assert 5.0 <= d <= 6.0

    def test_non_uniform_weights_change_distance(self) -> None:
        """With non-uniform weights and different rho distances per representative
        quartet, weighted_distance=True and False produce distinct results.

        Setup: partition {A,B}|{C}|{D}|{E}, 5-taxon dense profileset.
          - Profile {A,C,D,E}: split AC|DE  (A and C same side → rho_c), weight=1.0
          - Profile {B,C,D,E}: split BD|CE  (B and C diff sides → rho_s), weight=0.5
          - Three filler profiles (not used by this partition), weight=1.0

        For pair ({A,B}, {C}) the two representative quartets contribute:
          rho_dist(A,C) = 0.5 from {A,C,D,E} and rho_dist(B,C) = 1.0 from {B,C,D,E}.
          delta_1 = 1.0 (w=1.0), delta_2 = 2.0 (w=0.5).

          Unweighted: avg=1.5, confidence=0.5, shrink=0 → scaled=1.5 → d=5.5
          Weighted:   avg=2/1.5=4/3, confidence=1.0/1.5=2/3,
                      scaled = 4/3 - (4/3-3/2)*(1-2/3) = 4/3+1/18 = 25/18 → d=4+25/18≈5.389
        """
        taxa = frozenset({'A', 'B', 'C', 'D', 'E'})
        profiles = [
            # Filler: not reached by the partition's representative partitions
            (SqQuartetProfile([Quartet(Split({'A', 'B'}, {'C', 'D'}))]), 1.0),
            (SqQuartetProfile([Quartet(Split({'A', 'B'}, {'C', 'E'}))]), 1.0),
            (SqQuartetProfile([Quartet(Split({'A', 'B'}, {'D', 'E'}))]), 1.0),
            # Used: A from {A,B} → profile {A,C,D,E}
            (SqQuartetProfile([Quartet(Split({'A', 'C'}, {'D', 'E'}))]), 1.0),
            # Used: B from {A,B} → profile {B,C,D,E}
            (SqQuartetProfile([Quartet(Split({'B', 'D'}, {'C', 'E'}))]), 0.5),
        ]
        ps = SqQuartetProfileSet(profiles=profiles, taxa=taxa)
        partition = Partition([{'A', 'B'}, {'C'}, {'D'}, {'E'}])

        set_ab = frozenset({'A', 'B'})
        set_c = frozenset({'C'})

        dm_uw = quartet_distance_with_partition(ps, partition, weighted_distance=False, representative_mode='average')
        dm_w = quartet_distance_with_partition(ps, partition, weighted_distance=True, representative_mode='average')

        assert dm_uw.get_distance(set_ab, set_c) == pytest.approx(5.5)
        assert dm_w.get_distance(set_ab, set_c) == pytest.approx(4 + 2.0 / 1.5)
        assert dm_uw.get_distance(set_ab, set_c) != pytest.approx(dm_w.get_distance(set_ab, set_c))

    def test_higher_weight_pulls_distance_toward_its_value(self) -> None:
        """The higher-weight profile dominates the weighted average.

        Two representatives for pair ({A,B}, {C}):
          rho_dist(A,C)=0.5, w_high;  rho_dist(B,C)=1.0, w_low.
        As w_high increases relative to w_low, weighted distance → 4 + 1.0 (same-side limit).
        As w_low increases, weighted distance → 4 + 2.0 (opposite-side limit).
        """
        taxa = frozenset({'A', 'B', 'C', 'D', 'E'})

        def make_ps(w_acde: float, w_bcde: float) -> SqQuartetProfileSet:
            return SqQuartetProfileSet(
                profiles=[
                    (SqQuartetProfile([Quartet(Split({'A', 'B'}, {'C', 'D'}))]), 1.0),
                    (SqQuartetProfile([Quartet(Split({'A', 'B'}, {'C', 'E'}))]), 1.0),
                    (SqQuartetProfile([Quartet(Split({'A', 'B'}, {'D', 'E'}))]), 1.0),
                    (SqQuartetProfile([Quartet(Split({'A', 'C'}, {'D', 'E'}))]), w_acde),
                    (SqQuartetProfile([Quartet(Split({'B', 'D'}, {'C', 'E'}))]), w_bcde),
                ],
                taxa=taxa,
            )

        partition = Partition([{'A', 'B'}, {'C'}, {'D'}, {'E'}])
        set_ab, set_c = frozenset({'A', 'B'}), frozenset({'C'})

        # High weight on same-side profile → distance closer to 4 + 1.0 = 5.0
        dm_near = quartet_distance_with_partition(
            make_ps(1.0, 0.01), partition, weighted_distance=True
        )
        # High weight on opposite-side profile → distance closer to 4 + 2.0 = 6.0
        dm_far = quartet_distance_with_partition(
            make_ps(0.01, 1.0), partition, weighted_distance=True
        )

        assert dm_near.get_distance(set_ab, set_c) < dm_far.get_distance(set_ab, set_c)
        assert dm_near.get_distance(set_ab, set_c) == pytest.approx(5.0, abs=0.1)
        assert dm_far.get_distance(set_ab, set_c) == pytest.approx(6.0, abs=0.1)

    def test_weighted_distance_is_symmetric(self) -> None:
        """weighted_distance=True still produces a symmetric matrix."""
        ps = _dense_five_taxon_profileset(weights=[1.0, 0.8, 0.6, 0.4, 0.2])
        partition = Partition([{'A'}, {'B'}, {'C'}, {'D'}, {'E'}])
        dm = quartet_distance_with_partition(ps, partition, weighted_distance=True)
        for s1 in dm.labels:
            for s2 in dm.labels:
                assert dm.get_distance(s1, s2) == pytest.approx(dm.get_distance(s2, s1))
