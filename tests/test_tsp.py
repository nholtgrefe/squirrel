"""
Tests for TSP solvers (physquirrel.algorithms.tsp).

Ported from phylozoo/tests/core/distance/test_operations.py.
"""

import pytest
import numpy as np

from phylozoo import DistanceMatrix
from phylozoo.core.primitives.circular_ordering import CircularOrdering
from physquirrel.algorithms.tsp import approximate_tsp_tour, optimal_tsp_tour


class TestOptimalTspTour:
    """Test optimal_tsp_tour function."""

    def test_small_tsp_tour(self) -> None:
        matrix = np.array([
            [0, 1, 2, 3],
            [1, 0, 1, 2],
            [2, 1, 0, 1],
            [3, 2, 1, 0]
        ])
        dm = DistanceMatrix(matrix, labels=['A', 'B', 'C', 'D'])
        tour = optimal_tsp_tour(dm)

        assert isinstance(tour, CircularOrdering)
        assert len(tour) == len(dm.labels)
        assert set(tour.order) == set(dm.labels)

        expected_tour = CircularOrdering(['A', 'B', 'C', 'D'])
        assert tour == expected_tour

        tour_distance = sum(
            dm.get_distance(tour.order[i], tour.order[(i + 1) % len(tour)])
            for i in range(len(tour))
        )
        assert abs(tour_distance - 6.0) < 1e-10

    def test_single_label_tour(self) -> None:
        matrix = np.array([[0]])
        dm = DistanceMatrix(matrix, labels=['A'])
        tour = optimal_tsp_tour(dm)
        assert isinstance(tour, CircularOrdering)
        assert tour.order == ('A',)

    def test_two_label_tour(self) -> None:
        matrix = np.array([[0, 1], [1, 0]])
        dm = DistanceMatrix(matrix, labels=['A', 'B'])
        tour = optimal_tsp_tour(dm)
        assert isinstance(tour, CircularOrdering)
        assert len(tour) == 2
        assert set(tour.order) == {'A', 'B'}

        tour_distance = sum(
            dm.get_distance(tour.order[i], tour.order[(i + 1) % len(tour)])
            for i in range(len(tour))
        )
        assert abs(tour_distance - 2.0) < 1e-10

    def test_optimal_tour_property(self) -> None:
        matrix = np.array([
            [0, 2, 3, 4],
            [2, 0, 1, 2],
            [3, 1, 0, 1],
            [4, 2, 1, 0]
        ])
        dm = DistanceMatrix(matrix, labels=['A', 'B', 'C', 'D'])
        tour = optimal_tsp_tour(dm)

        assert isinstance(tour, CircularOrdering)
        assert len(tour) == len(set(tour.order))
        assert set(tour.order) == set(dm.labels)

        expected_tour = CircularOrdering(['A', 'B', 'C', 'D'])
        assert tour == expected_tour

        tour_distance = sum(
            dm.get_distance(tour.order[i], tour.order[(i + 1) % len(tour)])
            for i in range(len(tour))
        )
        assert abs(tour_distance - 8.0) < 1e-10


class TestApproximateTspTour:
    """Test approximate_tsp_tour function."""

    def test_simulated_annealing(self) -> None:
        matrix = np.array([
            [0, 1, 2, 3],
            [1, 0, 1, 2],
            [2, 1, 0, 1],
            [3, 2, 1, 0]
        ])
        dm = DistanceMatrix(matrix, labels=['A', 'B', 'C', 'D'])
        tour = approximate_tsp_tour(dm, method='simulated_annealing')

        assert isinstance(tour, CircularOrdering)
        assert len(tour) == len(dm.labels)
        assert set(tour.order) == set(dm.labels)

    def test_greedy_method(self) -> None:
        matrix = np.array([
            [0, 1, 2, 3],
            [1, 0, 1, 2],
            [2, 1, 0, 1],
            [3, 2, 1, 0]
        ])
        dm = DistanceMatrix(matrix, labels=['A', 'B', 'C', 'D'])
        tour = approximate_tsp_tour(dm, method='greedy')

        assert isinstance(tour, CircularOrdering)
        assert len(tour) == len(dm.labels)
        assert set(tour.order) == set(dm.labels)

    def test_christofides_method(self) -> None:
        matrix = np.array([
            [0, 1, 2, 3],
            [1, 0, 1, 2],
            [2, 1, 0, 1],
            [3, 2, 1, 0]
        ])
        dm = DistanceMatrix(matrix, labels=['A', 'B', 'C', 'D'])
        tour = approximate_tsp_tour(dm, method='christofides')

        assert isinstance(tour, CircularOrdering)
        assert len(tour) == len(dm.labels)
        assert set(tour.order) == set(dm.labels)

    def test_invalid_method_raises_error(self) -> None:
        matrix = np.array([[0, 1], [1, 0]])
        dm = DistanceMatrix(matrix, labels=['A', 'B'])
        with pytest.raises(ValueError, match="Method must be one of"):
            approximate_tsp_tour(dm, method='invalid_method')

    def test_single_label_approximate(self) -> None:
        matrix = np.array([[0]])
        dm = DistanceMatrix(matrix, labels=['A'])
        tour = approximate_tsp_tour(dm, method='greedy')
        assert isinstance(tour, CircularOrdering)
        assert tour.order == ('A',)

    def test_empty_matrix(self) -> None:
        matrix = np.zeros((0, 0))
        dm = DistanceMatrix(matrix, labels=[])
        tour1 = approximate_tsp_tour(dm, method='greedy')
        tour2 = optimal_tsp_tour(dm)
        assert isinstance(tour1, CircularOrdering)
        assert isinstance(tour2, CircularOrdering)
        assert len(tour1) == 0
        assert len(tour2) == 0


class TestTspTourProperties:
    """Test properties that hold for all TSP solvers."""

    def test_tour_contains_all_labels(self) -> None:
        matrix = np.array([
            [0, 1, 2, 3, 4],
            [1, 0, 1, 2, 3],
            [2, 1, 0, 1, 2],
            [3, 2, 1, 0, 1],
            [4, 3, 2, 1, 0]
        ])
        dm = DistanceMatrix(matrix, labels=['A', 'B', 'C', 'D', 'E'])

        for method in ('simulated_annealing', 'greedy', 'christofides'):
            tour = approximate_tsp_tour(dm, method=method)
            assert isinstance(tour, CircularOrdering)
            assert set(tour.order) == set(dm.labels)
            assert len(tour) == len(dm.labels)

    def test_optimal_vs_approximate_same_labels(self) -> None:
        matrix = np.array([
            [0, 1, 2, 3],
            [1, 0, 1, 2],
            [2, 1, 0, 1],
            [3, 2, 1, 0]
        ])
        dm = DistanceMatrix(matrix, labels=['A', 'B', 'C', 'D'])

        optimal = optimal_tsp_tour(dm)
        approximate = approximate_tsp_tour(dm, method='greedy')

        assert set(optimal.order) == set(approximate.order) == set(dm.labels)

    def test_deterministic_optimal(self) -> None:
        matrix = np.array([
            [0, 1, 2],
            [1, 0, 1],
            [2, 1, 0]
        ])
        dm = DistanceMatrix(matrix, labels=['A', 'B', 'C'])

        tour1 = optimal_tsp_tour(dm)
        tour2 = optimal_tsp_tour(dm)

        assert tour1 == tour2
        assert tour1 == CircularOrdering(['A', 'B', 'C'])

        tour_distance = sum(
            dm.get_distance(tour1.order[i], tour1.order[(i + 1) % len(tour1)])
            for i in range(len(tour1))
        )
        assert abs(tour_distance - 4.0) < 1e-10

    def test_larger_matrix_approximate(self) -> None:
        n = 8
        matrix = np.array([[abs(i - j) for j in range(n)] for i in range(n)], dtype=float)
        dm = DistanceMatrix(matrix)

        tour = approximate_tsp_tour(dm, method='greedy')
        assert isinstance(tour, CircularOrdering)
        assert len(tour) == n
        assert set(tour.order) == set(range(n))
