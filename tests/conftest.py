"""
Shared pytest fixtures for the physquirrel test suite.
"""

import pytest
import numpy as np

from phylozoo import DistanceMatrix, Quartet, Split
from physquirrel import SqQuartetProfile, SqQuartetProfileSet, delta_heuristic


@pytest.fixture
def simple_four_taxon_profileset():
    """A minimal SqQuartetProfileSet with a single split quartet on {A, B, C, D}."""
    q = Quartet(Split({'A', 'B'}, {'C', 'D'}))
    return SqQuartetProfileSet(profiles=[SqQuartetProfile([q])])


@pytest.fixture
def cycle_four_taxon_profileset():
    """A SqQuartetProfileSet with a cycle profile on {A, B, C, D}."""
    q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
    q2 = Quartet(Split({'A', 'C'}, {'B', 'D'}))
    profile = SqQuartetProfile([q1, q2], reticulation_leaf='A')
    return SqQuartetProfileSet(profiles=[profile])


@pytest.fixture
def five_taxon_profileset():
    """A SqQuartetProfileSet with 5 taxa (A–E)."""
    quartets = [
        Quartet(Split({'A', 'B'}, {'C', 'D'})),
        Quartet(Split({'A', 'B'}, {'C', 'E'})),
        Quartet(Split({'A', 'B'}, {'D', 'E'})),
        Quartet(Split({'A', 'C'}, {'D', 'E'})),
        Quartet(Split({'B', 'C'}, {'D', 'E'})),
    ]
    return SqQuartetProfileSet(profiles=[SqQuartetProfile([q]) for q in quartets])


@pytest.fixture
def tree_distance_matrix():
    """Distance matrix producing a clean tree signal (AB close, CD close)."""
    matrix = np.array([
        [0.0, 0.1, 0.9, 0.9],
        [0.1, 0.0, 0.9, 0.9],
        [0.9, 0.9, 0.0, 0.1],
        [0.9, 0.9, 0.1, 0.0],
    ])
    return DistanceMatrix(matrix, labels=['A', 'B', 'C', 'D'])
