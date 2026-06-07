"""
physquirrel 2.0 — SQUIRREL inference package for phylogenetic networks.

Core datatypes are provided by phylozoo; this package owns the inference pipeline.
"""

__version__ = "2.0.0"

# Datatypes
from .datatypes import SqQuartetProfile, SqQuartetProfileSet, to_psq, from_psq, to_profile_list, from_profile_list

# Inference pipeline
from .algorithms import (
    squirrel,
    delta_heuristic,
    tstar_tree,
    bstar,
    quartet_joining,
    adapted_quartet_joining,
    unresolve_tree,
    split_support,
    resolve_cycles,
    sqprofileset_from_network,
    sqprofileset_similarity,
    quartet_distance_with_partition,
)

# Convenience helpers
from .convenience import delta_heuristic_from_msa, squirrel_from_distances, squirrel_from_msa

# Re-export commonly needed phylozoo types so users can do `from physquirrel import ...`
from phylozoo import (
    SemiDirectedPhyNetwork,
    DirectedPhyNetwork,
    Split,
    SplitSystem,
    DistanceMatrix,
    MSA,
    Quartet,
    QuartetProfile,
    QuartetProfileSet,
)

__all__ = [
    '__version__',
    # Inference
    'squirrel',
    'delta_heuristic',
    'SqQuartetProfile',
    'SqQuartetProfileSet',
    'tstar_tree',
    'bstar',
    'quartet_joining',
    'adapted_quartet_joining',
    'unresolve_tree',
    'split_support',
    'resolve_cycles',
    'sqprofileset_from_network',
    'sqprofileset_similarity',
    'quartet_distance_with_partition',
    'to_psq',
    'from_psq',
    'to_profile_list',
    'from_profile_list',
    'delta_heuristic_from_msa',
    'squirrel_from_distances',
    'squirrel_from_msa',
    # phylozoo re-exports
    'SemiDirectedPhyNetwork',
    'DirectedPhyNetwork',
    'Split',
    'SplitSystem',
    'DistanceMatrix',
    'MSA',
    'Quartet',
    'QuartetProfile',
    'QuartetProfileSet',
]
