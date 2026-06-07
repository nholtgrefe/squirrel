"""
Fixture SemiDirectedPhyNetwork instances for use in physquirrel tests.

Adapted from sandbox/phylozoo_source/tests/fixtures/sd_networks.py.
Only the fixtures actually used in the test suite are included.
"""

from phylozoo import SemiDirectedPhyNetwork

# ============================================================================
# TREES
# ============================================================================

SDTREE_SMALL_BINARY = SemiDirectedPhyNetwork(
    undirected_edges=[
        (3, 1),
        (3, 2),
        (3, 4)
    ],
    nodes=[
        (1, {'label': 'A'}),
        (2, {'label': 'B'}),
        (4, {'label': 'C'})
    ]
)
"""Small binary tree with 3 taxa (A, B, C)."""

SDTREE_MEDIUM_BINARY = SemiDirectedPhyNetwork(
    undirected_edges=[
        (11, 1),
        (11, 2),
        (11, 16),
        (16, 12),
        (16, 19),
        (12, 3),
        (12, 4),
        (13, 5),
        (13, 6),
        (13, 17),
        (17, 14),
        (17, 19),
        (14, 7),
        (14, 8),
        (15, 9),
        (15, 10),
        (15, 18),
        (18, 1018),
        (18, 20),
        (19, 20),
        (20, 1020)
    ],
    nodes=[
        (1, {'label': 'L1'}),
        (2, {'label': 'L2'}),
        (3, {'label': 'L3'}),
        (4, {'label': 'L4'}),
        (5, {'label': 'L5'}),
        (6, {'label': 'L6'}),
        (7, {'label': 'L7'}),
        (8, {'label': 'L8'}),
        (9, {'label': 'L9'}),
        (10, {'label': 'L10'}),
        (1018, {'label': 'Dummy1018'}),
        (1020, {'label': 'Dummy1020'})
    ]
)
"""Binary tree with 12 leaves (L1–L10, Dummy1018, Dummy1020)."""

# ============================================================================
# SIMPLE HYBRIDS
# ============================================================================

LEVEL_1_SDNETWORK_SINGLE_HYBRID = SemiDirectedPhyNetwork(
    directed_edges=[
        (5, 4),
        (6, 4)
    ],
    undirected_edges=[
        (5, 3),
        (5, 6),
        (6, 7),
        (4, 8),
        (8, 1),
        (8, 2)
    ],
    nodes=[
        (3, {'label': 'C'}),
        (7, {'label': 'D'}),
        (1, {'label': 'A'}),
        (2, {'label': 'B'})
    ]
)
"""Level-1 network with single hybrid node. Taxa: A, B, C, D. Hybrid node 4 has parents 5 and 6."""

LEVEL_1_SDNETWORK_TWO_HYBRIDS_SEPARATE = SemiDirectedPhyNetwork(
    directed_edges=[
        (5, 4),
        (6, 4),
        (7, 9),
        (8, 9)
    ],
    undirected_edges=[
        (5, 11),
        (5, 6),
        (11, 7),
        (11, 8),
        (6, 3),
        (7, 12),
        (8, 13),
        (4, 1),
        (9, 2)
    ],
    nodes=[
        (3, {'label': 'C'}),
        (1, {'label': 'A'}),
        (12, {'label': 'D'}),
        (13, {'label': 'E'}),
        (2, {'label': 'B'})
    ]
)
"""Level-1 network with two separate hybrid nodes. Taxa: A, B, C, D, E."""

LEVEL_2_SDNETWORK_DIAMOND_HYBRID = SemiDirectedPhyNetwork(
    directed_edges=[
        (7, 5),
        (7, 6),
        (8, 5),
        (8, 6),
        (11, 6)
    ],
    undirected_edges=[
        (10, 7),
        (10, 8),
        (10, 11),
        (11, 14),
        (11, 15),
        (5, 12),
        (12, 1),
        (12, 2),
        (6, 13),
        (13, 3),
        (13, 4)
    ],
    nodes=[
        (14, {'label': 'E'}),
        (15, {'label': 'F'}),
        (1, {'label': 'A'}),
        (2, {'label': 'B'}),
        (3, {'label': 'C'}),
        (4, {'label': 'D'})
    ]
)
"""Level-2 diamond hybrid network. Taxa: A, B, C, D, E, F."""

# ============================================================================
# PARALLEL EDGES
# ============================================================================

LEVEL_1_SDNETWORK_PARALLEL_EDGES = SemiDirectedPhyNetwork(
    directed_edges=[
        (10, 5),
        (10, 5, 1)
    ],
    undirected_edges=[
        (10, 6),
        (6, 3),
        (6, 4),
        (5, 11),
        (11, 1),
        (11, 2)
    ],
    nodes=[
        (3, {'label': 'C'}),
        (4, {'label': 'D'}),
        (1, {'label': 'A'}),
        (2, {'label': 'B'})
    ]
)
"""Level-1 network with parallel edges. Taxa: A, B, C, D."""
