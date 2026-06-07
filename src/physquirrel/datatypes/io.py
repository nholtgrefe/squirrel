"""
I/O module for SqQuartetProfileSet.

Two formats are supported:

**psq** (default, ``.psq``) — compact JSON, lossless round-trip::

    {
        "taxa": ["A", "B", "C", "D"],
        "profiles": [
            {
                "taxa": ["A", "B", "C", "D"],
                "quartets": [
                    {"quartet": {"type": "resolved",
                                 "split": {"set1": ["A", "B"], "set2": ["C", "D"]}},
                     "weight": 1.0}
                ],
                "reticulation_leaf": null,
                "profile_weight": 0.85
            }
        ]
    }

**profile_list** (``.txt``) — plain-text, one profile per line::

    # Lines starting with '#' are comments.
    SQ a b c d [weight]   # split quarnet: split ab|cd
    4C a b c d [weight]   # 4-cycle: circular order a,b,c,d; a is the reticulation leaf
"""

from __future__ import annotations

import json
from typing import Any, TYPE_CHECKING

from phylozoo.core.quartet.base import Quartet
from phylozoo.core.split.base import Split
from phylozoo.utils.exceptions import PhyloZooFormatError, PhyloZooParseError
from phylozoo.utils.io import FormatRegistry

from .sqprofile import SqQuartetProfile

if TYPE_CHECKING:
    from .sqprofileset import SqQuartetProfileSet


# ---------------------------------------------------------------------------
# psq (JSON) format
# ---------------------------------------------------------------------------

def to_psq(profileset: 'SqQuartetProfileSet', **kwargs: Any) -> str:
    """Serialize a :class:`~physquirrel.SqQuartetProfileSet` to a JSON string.

    The resulting string uses the ``.psq`` format, which stores the full taxon
    set, per-profile quartet weights, reticulation leaf assignments, and
    profile-level confidence weights. Round-tripping via :func:`from_psq` is
    lossless.

    Parameters
    ----------
    profileset : SqQuartetProfileSet
        The profile set to serialize.

    Returns
    -------
    str
        A compact JSON string. The top-level object has two keys:

        - ``"taxa"`` — sorted list of all taxon labels.
        - ``"profiles"`` — list of profile objects, each with:

          - ``"taxa"`` — the four taxon labels for this profile.
          - ``"quartets"`` — list of ``{"quartet": ..., "weight": float}``
            entries. A resolved quartet has type ``"resolved"`` with a
            ``"split"`` sub-object; a star quartet has type ``"star"`` with a
            ``"taxa"`` list.
          - ``"reticulation_leaf"`` — taxon label or ``null``.
          - ``"profile_weight"`` — confidence weight for this profile.
    """
    profiles_data = []
    for four_taxa, (profile, weight) in profileset.profiles.items():
        quartets_data = []
        for q, q_weight in profile.quartets.items():
            if q.is_resolved():
                split = q.split
                if split is None:
                    raise PhyloZooFormatError("Resolved quartet should have a split")
                quartet_data = {
                    'type': 'resolved',
                    'split': {
                        'set1': sorted(split.set1),
                        'set2': sorted(split.set2),
                    },
                }
            else:
                quartet_data = {
                    'type': 'star',
                    'taxa': sorted(q.taxa),
                }
            quartets_data.append({'quartet': quartet_data, 'weight': q_weight})

        profiles_data.append({
            'taxa': sorted(four_taxa),
            'quartets': quartets_data,
            'reticulation_leaf': profile.reticulation_leaf,
            'profile_weight': weight,
        })

    data = {
        'profiles': profiles_data,
        'taxa': sorted(profileset.taxa),
    }
    return json.dumps(data, indent=None)


def from_psq(psq_string: str, **kwargs: Any) -> 'SqQuartetProfileSet':
    """Parse a JSON string in ``.psq`` format into a :class:`~physquirrel.SqQuartetProfileSet`.

    Parameters
    ----------
    psq_string : str
        A JSON string as produced by :func:`to_psq` or by saving a profile
        set with ``profileset.save("file.psq")``.

    Returns
    -------
    SqQuartetProfileSet

    Raises
    ------
    PhyloZooParseError
        If *psq_string* is not valid JSON.
    PhyloZooFormatError
        If the JSON structure does not match the expected ``.psq`` schema.
    """
    try:
        data = json.loads(psq_string)
    except json.JSONDecodeError as e:
        raise PhyloZooParseError(f"Invalid JSON in psq format: {e}") from e

    if not isinstance(data, dict):
        raise PhyloZooFormatError("psq format must be a JSON object")
    if 'profiles' not in data or 'taxa' not in data:
        raise PhyloZooFormatError("psq format must contain 'profiles' and 'taxa' keys")

    reconstructed_profiles = []
    for profile_data in data['profiles']:
        if 'taxa' not in profile_data or 'quartets' not in profile_data:
            raise PhyloZooFormatError("Profile must contain 'taxa' and 'quartets' keys")

        quartets_list = []
        for q_data in profile_data['quartets']:
            if 'quartet' not in q_data or 'weight' not in q_data:
                raise PhyloZooFormatError("Quartet data must contain 'quartet' and 'weight' keys")

            qi = q_data['quartet']
            q_weight = q_data['weight']

            if qi['type'] == 'resolved':
                sd = qi['split']
                quartet = Quartet(Split(set(sd['set1']), set(sd['set2'])))
            elif qi['type'] == 'star':
                quartet = Quartet(frozenset(qi['taxa']))
            else:
                raise PhyloZooFormatError(f"Unknown quartet type: {qi['type']}")

            quartets_list.append((quartet, q_weight))

        sq_profile = SqQuartetProfile(
            quartets=quartets_list,
            reticulation_leaf=profile_data.get('reticulation_leaf'),
        )
        reconstructed_profiles.append((sq_profile, profile_data['profile_weight']))

    from .sqprofileset import SqQuartetProfileSet
    taxa = frozenset(data['taxa'])
    return SqQuartetProfileSet(profiles=reconstructed_profiles, taxa=taxa)


# ---------------------------------------------------------------------------
# profile_list (plain text) format
# ---------------------------------------------------------------------------

def to_profile_list(profileset: 'SqQuartetProfileSet', **kwargs: Any) -> str:
    """Serialize a :class:`~physquirrel.SqQuartetProfileSet` to a plain-text string.

    The ``.txt`` profile list format has one profile per line:

    - ``SQ a b c d weight`` — split quarnet with split :math:`ab|cd`.
    - ``4C a b c d weight`` — 4-cycle with circular order :math:`a,b,c,d`;
      :math:`a` is the reticulation leaf (placed first because it is below the
      hybrid node).

    This format is compatible with the quarnet file format used in earlier
    versions of physquirrel.

    Parameters
    ----------
    profileset : SqQuartetProfileSet
        The profile set to serialize.

    Returns
    -------
    str
        A plain-text string. Example output for a profile set with one split
        quarnet and one 4-cycle::

            SQ A B C D 0.85
            4C A B C E 0.60
    """
    lines = []
    for _four_taxa, (profile, weight) in profileset.profiles.items():
        quartets = list(profile.quartets.keys())
        if len(quartets) == 1:
            split = quartets[0].split
            if split is None:
                raise PhyloZooFormatError("Resolved quartet in profile_list must have a split")
            s1 = sorted(split.set1)
            s2 = sorted(split.set2)
            lines.append(f"SQ {s1[0]} {s1[1]} {s2[0]} {s2[1]} {weight}")
        else:
            ordering = profile.circular_ordering
            if ordering is None:
                raise PhyloZooFormatError(
                    "Cannot serialize 2-quartet profile without a valid circular ordering"
                )
            order = list(ordering.order)
            if profile.reticulation_leaf is not None and profile.reticulation_leaf in order:
                idx = order.index(profile.reticulation_leaf)
                order = order[idx:] + order[:idx]
            lines.append(f"4C {order[0]} {order[1]} {order[2]} {order[3]} {weight}")
    return '\n'.join(lines)


def from_profile_list(profile_list_string: str, **kwargs: Any) -> 'SqQuartetProfileSet':
    """Parse a plain-text profile list string into a :class:`~physquirrel.SqQuartetProfileSet`.

    Each non-blank, non-comment line must follow one of two formats:

    - ``SQ a b c d [weight]`` — a split quarnet with split :math:`ab|cd`.
      The first two taxon names go on one side of the split, the last two on
      the other.
    - ``4C a b c d [weight]`` — a 4-cycle with circular order
      :math:`a,b,c,d`. The first taxon (:math:`a`) is treated as the
      reticulation leaf.

    Lines starting with ``#`` are treated as comments and ignored. The weight
    field is optional and defaults to 1.0.

    Parameters
    ----------
    profile_list_string : str
        A string in the profile list format, as produced by
        :func:`to_profile_list` or written by hand. Example::

            # My quartet profiles
            SQ Human Chimp Gorilla Orang 0.9
            4C Human Chimp Gorilla Macaque 0.6

    Returns
    -------
    SqQuartetProfileSet

    Raises
    ------
    PhyloZooParseError
        If a line has too few tokens.
    PhyloZooFormatError
        If a line begins with an unrecognised type tag (not ``SQ`` or ``4C``).
    """
    profiles = []
    for lineno, raw in enumerate(profile_list_string.splitlines(), start=1):
        line = raw.strip()
        if not line or line.startswith('#'):
            continue
        tokens = line.split()
        kind = tokens[0]
        if kind == 'SQ':
            if len(tokens) < 5:
                raise PhyloZooParseError(
                    f"Line {lineno}: SQ requires 4 taxa, got {len(tokens) - 1}: {line!r}"
                )
            a, b, c, d = tokens[1], tokens[2], tokens[3], tokens[4]
            weight = float(tokens[5]) if len(tokens) > 5 else 1.0
            quartet = Quartet(Split({a, b}, {c, d}))
            profile = SqQuartetProfile([quartet])
            profiles.append((profile, weight))
        elif kind == '4C':
            if len(tokens) < 5:
                raise PhyloZooParseError(
                    f"Line {lineno}: 4C requires 4 taxa, got {len(tokens) - 1}: {line!r}"
                )
            a, b, c, d = tokens[1], tokens[2], tokens[3], tokens[4]
            weight = float(tokens[5]) if len(tokens) > 5 else 1.0
            q1 = Quartet(Split({a, b}, {c, d}))
            q2 = Quartet(Split({a, d}, {b, c}))
            profile = SqQuartetProfile([q1, q2], reticulation_leaf=a)
            profiles.append((profile, weight))
        else:
            raise PhyloZooFormatError(
                f"Line {lineno}: unknown profile type {kind!r}, expected 'SQ' or '4C'"
            )

    from .sqprofileset import SqQuartetProfileSet
    return SqQuartetProfileSet(profiles=profiles)


# ---------------------------------------------------------------------------
# Format registration
# ---------------------------------------------------------------------------

def _register_formats() -> None:
    from .sqprofileset import SqQuartetProfileSet
    FormatRegistry.register(
        SqQuartetProfileSet,
        'psq',
        reader=from_psq,
        writer=to_psq,
        extensions=['.psq'],
        default=True,
    )
    FormatRegistry.register(
        SqQuartetProfileSet,
        'profile_list',
        reader=from_profile_list,
        writer=to_profile_list,
        extensions=['.txt'],
        default=False,
    )


_register_formats()
