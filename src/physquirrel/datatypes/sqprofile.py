from typing import TYPE_CHECKING, Mapping

from phylozoo.core.quartet.qprofile import QuartetProfile
from phylozoo.core.primitives.circular_ordering import CircularOrdering
from phylozoo.utils.exceptions import PhyloZooValueError

if TYPE_CHECKING:
    from phylozoo.core.quartet.base import Quartet


class SqQuartetProfile(QuartetProfile):
    """
    Squirrel quartet profile: exactly 1 or 2 resolved quartets.

    If it contains 2 quartets (forming a 4-cycle), it may have an optional
    reticulation_leaf identifying the leaf below the hybrid node.
    """

    __slots__ = ('_reticulation_leaf',)

    def __init__(
        self,
        quartets: (
            dict['Quartet', float]
            | Mapping['Quartet', float]
            | list['Quartet']
            | list[tuple['Quartet', float]]
        ),
        reticulation_leaf: str | None = None,
    ) -> None:
        super().__init__(quartets)

        resolved = [q for q in self._quartets if q.is_resolved()]
        unresolved = [q for q in self._quartets if not q.is_resolved()]

        if unresolved:
            raise PhyloZooValueError(
                f"SqQuartetProfile cannot contain unresolved quartets, "
                f"got {len(unresolved)} unresolved"
            )

        if len(resolved) not in (1, 2):
            raise PhyloZooValueError(
                f"SqQuartetProfile must contain exactly 1 or 2 resolved quartets, "
                f"got {len(resolved)}"
            )

        if len(resolved) == 2:
            orderings = self.circular_orderings
            if not orderings:
                raise PhyloZooValueError(
                    "SqQuartetProfile with 2 quartets must form a circular ordering."
                )
            if reticulation_leaf is not None and reticulation_leaf not in self._taxa:
                raise PhyloZooValueError(
                    f"Reticulation leaf '{reticulation_leaf}' must be one of the taxa: {set(self._taxa)}"
                )
        else:
            if reticulation_leaf is not None:
                raise PhyloZooValueError(
                    "reticulation_leaf can only be provided when profile has 2 quartets"
                )

        object.__setattr__(self, '_reticulation_leaf', reticulation_leaf)

    @property
    def reticulation_leaf(self) -> str | None:
        return self._reticulation_leaf

    @property
    def circular_ordering(self) -> CircularOrdering | None:
        if len(self._quartets) == 2:
            orderings = self.circular_orderings
            if not orderings:
                return None
            return next(iter(orderings))
        return None

    def __repr__(self) -> str:
        ret = f"reticulation_leaf='{self._reticulation_leaf}'" if self._reticulation_leaf else "reticulation_leaf=None"
        return (
            f"SqQuartetProfile(taxa={set(self._taxa)}, {ret}, quartets={dict(self._quartets)})"
        )
