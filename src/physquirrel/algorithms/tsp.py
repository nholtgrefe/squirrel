"""
TSP solvers used by the cycle resolution algorithm.

Implements:
  - optimal_tsp_tour   — Held-Karp DP, O(n^2 * 2^n), Numba-accelerated
  - approximate_tsp_tour — networkx heuristics (simulated annealing / greedy / Christofides)
"""

from __future__ import annotations

import itertools

import networkx as nx
import numpy as np
from numba import njit

from phylozoo import DistanceMatrix
from phylozoo.core.primitives.circular_ordering import CircularOrdering
from phylozoo.utils.exceptions import PhyloZooValueError


@njit(cache=True)
def _held_karp_numba(
    matrix: np.ndarray, n: int, memo: np.ndarray, next_node: np.ndarray
) -> float:
    """Numba-accelerated Held-Karp algorithm using bitmasks."""
    for i in range(n):
        memo[i, 0] = matrix[i, 0]

    max_mask = 1 << (n - 1)

    for mask in range(1, max_mask):
        for ni in range(n):
            if ni == 0:
                continue
            if mask & (1 << (ni - 1)):
                continue
            if mask == 0:
                min_cost = matrix[ni, 0]
                best_nj = 0
            else:
                min_cost = np.inf
                best_nj = -1
                for nj in range(1, n):
                    if nj == ni:
                        continue
                    if not (mask & (1 << (nj - 1))):
                        continue
                    mask_without_nj = mask & ~(1 << (nj - 1))
                    cost = matrix[ni, nj] + memo[nj, mask_without_nj]
                    if cost < min_cost:
                        min_cost = cost
                        best_nj = nj
            memo[ni, mask] = min_cost
            next_node[ni, mask] = best_nj

    full_mask = max_mask - 1
    min_cost = np.inf
    best_first = -1
    for nj in range(1, n):
        cost = matrix[0, nj] + memo[nj, full_mask & ~(1 << (nj - 1))]
        if cost < min_cost:
            min_cost = cost
            best_first = nj

    return min_cost


def optimal_tsp_tour(distance_matrix: DistanceMatrix) -> CircularOrdering:
    """
    Solve TSP to optimality using the Held-Karp dynamic programming algorithm.

    Time complexity O(n^2 * 2^n) — practical for n ≤ ~20.

    Parameters
    ----------
    distance_matrix : DistanceMatrix

    Returns
    -------
    CircularOrdering
        Canonical circular ordering of all labels.
    """
    n = len(distance_matrix)
    if n == 0:
        return CircularOrdering([])
    if n == 1:
        return CircularOrdering([distance_matrix.labels[0]])

    matrix = np.ascontiguousarray(distance_matrix._matrix, dtype=np.float64)

    max_mask = 1 << (n - 1)
    memo = np.full((n, max_mask), np.inf, dtype=np.float64)
    next_node = np.full((n, max_mask), -1, dtype=np.int32)

    _held_karp_numba(matrix, n, memo, next_node)

    solution = [0]
    full_mask = max_mask - 1

    min_cost = np.inf
    best_first = -1
    for nj in range(1, n):
        mask_without_nj = full_mask & ~(1 << (nj - 1))
        cost = matrix[0, nj] + memo[nj, mask_without_nj]
        if cost < min_cost:
            min_cost = cost
            best_first = nj

    ni = best_first
    current_mask = full_mask & ~(1 << (ni - 1))
    solution.append(ni)

    while current_mask != 0:
        next_ni = int(next_node[ni, current_mask])
        if next_ni == 0:
            break
        solution.append(next_ni)
        ni = next_ni
        current_mask = current_mask & ~(1 << (ni - 1))

    label_tour = [distance_matrix.labels[i] for i in solution]
    return CircularOrdering(label_tour)


def approximate_tsp_tour(
    distance_matrix: DistanceMatrix,
    method: str = 'simulated_annealing',
) -> CircularOrdering:
    """
    Find an approximate TSP tour using a networkx heuristic.

    Parameters
    ----------
    distance_matrix : DistanceMatrix
    method : str
        One of ``'simulated_annealing'`` (default), ``'greedy'``,
        or ``'christofides'``.

    Returns
    -------
    CircularOrdering
        Canonical circular ordering of all labels.

    Raises
    ------
    PhyloZooValueError
        If *method* is not one of the supported options.
    """
    if method not in ('simulated_annealing', 'greedy', 'christofides'):
        raise PhyloZooValueError(
            f"Method must be one of ['simulated_annealing', 'greedy', 'christofides'], "
            f"got '{method}'"
        )

    n = len(distance_matrix)
    if n == 0:
        return CircularOrdering([])
    if n == 1:
        return CircularOrdering([distance_matrix.labels[0]])

    complete_graph = nx.Graph()
    complete_graph.add_nodes_from(distance_matrix.indices)
    complete_graph.add_weighted_edges_from(
        (i, j, distance_matrix._matrix[i, j])
        for i, j in itertools.combinations(distance_matrix.indices, 2)
    )

    tsp_func = nx.approximation.traveling_salesman_problem

    if method == 'simulated_annealing':
        tour = tsp_func(
            complete_graph,
            cycle=False,
            method=nx.approximation.simulated_annealing_tsp,
            init_cycle='greedy',
            seed=0,
        )
    elif method == 'greedy':
        tour = tsp_func(
            complete_graph, cycle=False, method=nx.approximation.greedy_tsp
        )
    else:  # christofides
        tour = tsp_func(
            complete_graph, cycle=False, method=nx.approximation.christofides
        )

    label_tour = [distance_matrix.labels[i] for i in tour]
    return CircularOrdering(label_tour)
