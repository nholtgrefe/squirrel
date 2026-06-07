from typing import Mapping

from phylozoo.core.quartet.base import Quartet
from phylozoo.core.quartet.qprofileset import QuartetProfileSet
from phylozoo.utils.exceptions import PhyloZooValueError
from phylozoo.utils.io import IOMixin

from .sqprofile import SqQuartetProfile


class SqQuartetProfileSet(QuartetProfileSet, IOMixin):
    """
    Collection of SqQuartetProfile objects with .psq/.txt I/O support.

    Accepts SqQuartetProfile objects, Quartet objects (auto-grouped by taxa),
    or tuples of (profile/quartet, weight).
    """

    def __init__(
        self,
        profiles: list[SqQuartetProfile | Quartet | tuple[SqQuartetProfile, float] | tuple[Quartet, float]] | None = None,
        taxa: frozenset[str] | None = None,
    ) -> None:
        if profiles is None:
            profiles = []

        sq_profiles: list[tuple[SqQuartetProfile, float]] = []
        quartet_data: dict[frozenset[str], dict[Quartet, float]] = {}

        for item in profiles:
            if isinstance(item, tuple):
                obj, weight = item
            else:
                obj = item
                weight = 1.0

            if isinstance(obj, SqQuartetProfile):
                sq_profiles.append((obj, weight))
            elif isinstance(obj, Quartet):
                qtaxa = obj.taxa
                if qtaxa not in quartet_data:
                    quartet_data[qtaxa] = {}
                if obj in quartet_data[qtaxa]:
                    raise PhyloZooValueError(
                        f"Quartet {obj} appears multiple times in the input."
                    )
                quartet_data[qtaxa][obj] = weight
            else:
                raise PhyloZooValueError(f"Expected SqQuartetProfile or Quartet, got {type(obj)}")

        converted: list[tuple[SqQuartetProfile, float]] = list(sq_profiles)
        for taxa_set, quartets_dict in quartet_data.items():
            if not quartets_dict:
                raise PhyloZooValueError(f"Cannot have empty profile for taxa {taxa_set}")
            profile_weight = sum(quartets_dict.values())
            # Normalize within-profile weights to sum to 1.0 (required by QuartetProfile)
            normalized = {q: w / profile_weight for q, w in quartets_dict.items()}
            sq_profile = SqQuartetProfile(normalized)
            converted.append((sq_profile, profile_weight))

        super().__init__(profiles=converted, taxa=taxa)

        for profile, _ in self._profiles.values():
            if not isinstance(profile, SqQuartetProfile):
                raise PhyloZooValueError(
                    f"All profiles must be SqQuartetProfile instances, got {type(profile)}"
                )

    @property
    def profiles(self) -> Mapping[frozenset[str], tuple[SqQuartetProfile, float]]:
        return self._profiles  # type: ignore

    def get_profile(self, taxa: frozenset[str]) -> SqQuartetProfile | None:
        result = self._profiles.get(taxa)
        return result[0] if result else None  # type: ignore

    def __repr__(self) -> str:
        if not self._profiles:
            return "SqQuartetProfileSet(profiles={})"
        items = []
        for taxa, (profile, weight) in self._profiles.items():
            qd = dict(profile.quartets)
            ret = f", reticulation_leaf='{profile.reticulation_leaf}'" if profile.reticulation_leaf else ""
            items.append(f"(SqQuartetProfile({repr(qd)}{ret}), {weight})")
        return f"SqQuartetProfileSet(profiles=[{', '.join(items)}])"

    _default_format: str = 'psq'
    _supported_formats: list[str] = ['psq', 'profile_list']
