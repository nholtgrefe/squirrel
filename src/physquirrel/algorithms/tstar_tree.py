import itertools

from phylozoo import SemiDirectedPhyNetwork, Split, SplitSystem
from phylozoo.core.quartet.qprofileset import QuartetProfileSet
from phylozoo.core.split.algorithms import tree_from_splitsystem
from phylozoo.utils.exceptions import PhyloZooValueError


def bstar(profileset: QuartetProfileSet) -> SplitSystem:
    """
    Compute the B*-set of compatible splits from quartet profiles.

    Uses the incremental O(n^5) algorithm from Berry & Gascuel (2000).
    Only considers trivial (single-quartet), resolved profiles.

    Parameters
    ----------
    profileset : QuartetProfileSet

    Returns
    -------
    SplitSystem
    """
    all_taxa = list(profileset.taxa)
    if len(all_taxa) < 4:
        raise PhyloZooValueError("B* algorithm requires at least 4 taxa")

    quartet_splits: dict[frozenset[str], Split] = {}
    for taxa_set, (profile, _) in profileset.profiles.items():
        if not profile.is_trivial():
            continue
        split = profile.split
        if split is not None:
            quartet_splits[taxa_set] = split

    def obtain_split(four_elements: frozenset[str]) -> Split | None:
        return quartet_splits.get(four_elements)

    order = all_taxa.copy()
    a, b, c, d = order[0:4]

    bstar_splits = [
        Split({a}, {b, c, d}),
        Split({b}, {a, c, d}),
        Split({c}, {a, b, d}),
        Split({d}, {a, b, c}),
    ]
    abcd_split = obtain_split(frozenset({a, b, c, d}))
    if abcd_split is not None:
        bstar_splits.append(abcd_split)

    for i, element in enumerate(order):
        if i < 4:
            continue

        previous_taxa = frozenset(order[0:i])
        new_bstar = [Split({element}, previous_taxa)]

        for split in bstar_splits:
            add1 = True
            for x in split.set1:
                for y, z in itertools.combinations(split.set2, 2):
                    qs = obtain_split(frozenset({element, x, y, z}))
                    if qs is not None and qs != Split({x, element}, {y, z}):
                        add1 = False
                        break
                if not add1:
                    break

            add2 = True
            for x, y in itertools.combinations(split.set1, 2):
                for z in split.set2:
                    qs = obtain_split(frozenset({element, x, y, z}))
                    if qs is not None and qs != Split({x, y}, {z, element}):
                        add2 = False
                        break
                if not add2:
                    break

            if add1:
                new_bstar.append(Split(split.set1 | {element}, split.set2))
            if add2:
                new_bstar.append(Split(split.set1, split.set2 | {element}))

        bstar_splits = new_bstar

    return SplitSystem(bstar_splits)


def tstar_tree(profileset: QuartetProfileSet) -> SemiDirectedPhyNetwork:
    """
    Compute the T* tree from a quartet profile set.

    Builds the B*-set of splits then converts to a SemiDirectedPhyNetwork tree.

    Parameters
    ----------
    profileset : QuartetProfileSet

    Returns
    -------
    SemiDirectedPhyNetwork
    """
    bstar_splits = bstar(profileset)
    return tree_from_splitsystem(bstar_splits, check_compatibility=False)
