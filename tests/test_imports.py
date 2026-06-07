"""
Smoke tests for the physquirrel public API surface.
Verifies that all re-exports from physquirrel are importable.
"""


class TestTopLevelImports:
    """Test that all public symbols are importable from physquirrel."""

    def test_version(self) -> None:
        from physquirrel import __version__
        assert __version__ == "2.0.0"

    def test_squirrel(self) -> None:
        from physquirrel import squirrel
        assert callable(squirrel)

    def test_delta_heuristic(self) -> None:
        from physquirrel import delta_heuristic
        assert callable(delta_heuristic)

    def test_sq_quartet_profile(self) -> None:
        from physquirrel import SqQuartetProfile
        assert SqQuartetProfile is not None

    def test_sq_quartet_profileset(self) -> None:
        from physquirrel import SqQuartetProfileSet
        assert SqQuartetProfileSet is not None

    def test_tstar_tree(self) -> None:
        from physquirrel import tstar_tree
        assert callable(tstar_tree)

    def test_bstar(self) -> None:
        from physquirrel import bstar
        assert callable(bstar)

    def test_quartet_joining(self) -> None:
        from physquirrel import quartet_joining
        assert callable(quartet_joining)

    def test_adapted_quartet_joining(self) -> None:
        from physquirrel import adapted_quartet_joining
        assert callable(adapted_quartet_joining)

    def test_unresolve_tree(self) -> None:
        from physquirrel import unresolve_tree
        assert callable(unresolve_tree)

    def test_split_support(self) -> None:
        from physquirrel import split_support
        assert callable(split_support)

    def test_resolve_cycles(self) -> None:
        from physquirrel import resolve_cycles
        assert callable(resolve_cycles)

    def test_sqprofileset_from_network(self) -> None:
        from physquirrel import sqprofileset_from_network
        assert callable(sqprofileset_from_network)

    def test_sqprofileset_similarity(self) -> None:
        from physquirrel import sqprofileset_similarity
        assert callable(sqprofileset_similarity)

    def test_to_psq(self) -> None:
        from physquirrel import to_psq
        assert callable(to_psq)

    def test_from_psq(self) -> None:
        from physquirrel import from_psq
        assert callable(from_psq)

    def test_to_profile_list(self) -> None:
        from physquirrel import to_profile_list
        assert callable(to_profile_list)

    def test_from_profile_list(self) -> None:
        from physquirrel import from_profile_list
        assert callable(from_profile_list)

    def test_delta_heuristic_from_msa(self) -> None:
        from physquirrel import delta_heuristic_from_msa
        assert callable(delta_heuristic_from_msa)

    def test_squirrel_from_distances(self) -> None:
        from physquirrel import squirrel_from_distances
        assert callable(squirrel_from_distances)

    def test_squirrel_from_msa(self) -> None:
        from physquirrel import squirrel_from_msa
        assert callable(squirrel_from_msa)


class TestPhylozooReexports:
    """Test that phylozoo types are re-exported from physquirrel."""

    def test_semi_directed_phy_network(self) -> None:
        from physquirrel import SemiDirectedPhyNetwork
        assert SemiDirectedPhyNetwork is not None

    def test_directed_phy_network(self) -> None:
        from physquirrel import DirectedPhyNetwork
        assert DirectedPhyNetwork is not None

    def test_split(self) -> None:
        from physquirrel import Split
        assert Split is not None

    def test_split_system(self) -> None:
        from physquirrel import SplitSystem
        assert SplitSystem is not None

    def test_distance_matrix(self) -> None:
        from physquirrel import DistanceMatrix
        assert DistanceMatrix is not None

    def test_msa(self) -> None:
        from physquirrel import MSA
        assert MSA is not None

    def test_quartet(self) -> None:
        from physquirrel import Quartet
        assert Quartet is not None

    def test_quartet_profile(self) -> None:
        from physquirrel import QuartetProfile
        assert QuartetProfile is not None

    def test_quartet_profileset(self) -> None:
        from physquirrel import QuartetProfileSet
        assert QuartetProfileSet is not None


class TestModuleStructure:
    """Test that the submodule structure is correct."""

    def test_datatypes_module(self) -> None:
        import physquirrel.datatypes as dt
        assert hasattr(dt, 'SqQuartetProfile')
        assert hasattr(dt, 'SqQuartetProfileSet')
        assert hasattr(dt, 'to_psq')
        assert hasattr(dt, 'from_psq')
        assert hasattr(dt, 'to_profile_list')
        assert hasattr(dt, 'from_profile_list')

    def test_algorithms_module(self) -> None:
        import physquirrel.algorithms as alg
        assert hasattr(alg, 'squirrel')
        assert hasattr(alg, 'delta_heuristic')
        assert hasattr(alg, 'tstar_tree')
        assert hasattr(alg, 'bstar')
        assert hasattr(alg, 'quartet_joining')
        assert hasattr(alg, 'adapted_quartet_joining')
        assert hasattr(alg, 'unresolve_tree')
        assert hasattr(alg, 'split_support')
        assert hasattr(alg, 'resolve_cycles')
        assert hasattr(alg, 'sqprofileset_from_network')
        assert hasattr(alg, 'sqprofileset_similarity')
