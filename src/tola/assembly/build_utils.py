"""
Utility objects used by BuildAssembly
"""

import logging
import textwrap
from abc import ABC, abstractmethod

from tola.assembly.fragment import Fragment
from tola.assembly.overlap_result import OverlapResult
from tola.assembly.scaffold import Scaffold

log = logging.getLogger(__name__)


class FoundFragment:
    """
    Little object to store a Fragment found and the list of Scaffolds it was
    found in.
    """

    __slots__ = "fragment", "scaffolds"

    def __init__(self, fragment: Fragment):
        self.fragment = fragment
        self.scaffolds = []

    @property
    def scaffold_count(self):
        return len(self.scaffolds)

    def add_scaffold(self, scaffold: Scaffold):
        self.scaffolds.append(scaffold)

    def remove_scaffold(self, scaffold: Scaffold):
        self.scaffolds.remove(scaffold)


class OverhangPremise(ABC):
    """
    Stores a "what-if" for removal of a terminal (start or end) Fragment. Used
    to decide which OverlapResult to remove a Fragment from, where the
    Fragment is present in more than one OverlapResult.
    """

    __slots__ = "scaffold", "fragment"

    def __init__(self, scaffold: OverlapResult, fragment: Fragment):
        self.scaffold = scaffold
        self.fragment = fragment

    def __str__(self):
        return (
            f"{self.__class__.__name__}\n"
            f"  bait overlap: {self.bait_overlap:12_d}\n  if applied:\n"
            f"      overhang: {self.overhang_if_applied:12_d}\n"
            f"   error delta: {self.overhang_error_delta_if_applied:12_d}\n\n"
            + textwrap.indent(f"{self.scaffold}\n", "  ")
        )

    @property
    @abstractmethod
    def bait_overlap(self) -> int:
        """
        The overlap between the bait and the `Fragment` on this end of the
        `OverlapResult`.  (`0` if they don't overlap.)
        """

    @property
    @abstractmethod
    def overhang_if_applied(self) -> int:
        """
        Overhang of the `Fragment` at this end of the `OverlapResult` beyond
        the `bait` that would be left if `Fragment` at this end was removed.
        Value is negative if removing the `Fragment` would leave an
        underhang.
        """

    @property
    @abstractmethod
    def overhang_error_delta_if_applied(self) -> int:
        """
        Change (positive or negative) in the absolute size of the overhang
        (or underhang) if this `OverhangPremise` were applied.
        """

    @abstractmethod
    def apply(self):
        """
        Carry out this premise's action by removing the scaffold at this end.
        """

    def improves(self, err_length) -> bool:
        if len(self.scaffold.rows) == 1:
            return False
        return self.overhang_error_delta_if_applied < 0 and (
            # Guard against removing fragments which would produce a large
            # negative overhang - they should be cut instead.
            self.overhang_if_applied > -3 * err_length
        )

    def makes_worse(self, err_length) -> bool:
        return not self.improves(err_length)


class StartOverhangPremise(OverhangPremise):
    @property
    def bait_overlap(self) -> int:
        return self.scaffold.start_row_bait_overlap

    @property
    def overhang_if_applied(self) -> int:
        return self.scaffold.overhang_if_start_removed()

    @property
    def overhang_error_delta_if_applied(self) -> int:
        return abs(self.scaffold.overhang_if_start_removed()) - abs(
            self.scaffold.start_overhang
        )

    def apply(self):
        self.scaffold.discard_start()


class EndOverhangPremise(OverhangPremise):
    @property
    def bait_overlap(self) -> int:
        return self.scaffold.end_row_bait_overlap

    @property
    def overhang_if_applied(self) -> int:
        return self.scaffold.overhang_if_end_removed()

    @property
    def overhang_error_delta_if_applied(self) -> int:
        return abs(self.scaffold.overhang_if_end_removed()) - abs(
            self.scaffold.end_overhang
        )

    def apply(self):
        self.scaffold.discard_end()


class OverhangResolver:
    """
    Takes in a list of "problem" OverlapResults which share a Fragment.
    Performs one round of comparing OverlapResult pairs, choosing which of
    the two to remove the shared, terminal Fragment from. Returns a list of
    the OverlapPremises which were applied.
    """

    def __init__(self, error_length: int):
        self.premises_by_fragment_key: dict[
            tuple[str, int, int], list[OverhangPremise]
        ] = {}
        self.error_length = error_length

    def add_overhang_premise(self, fragment: Fragment, scffld: OverlapResult):
        if scffld.rows[0] is fragment:
            premise = StartOverhangPremise(scffld, fragment)
        elif scffld.rows[-1] is fragment:
            premise = EndOverhangPremise(scffld, fragment)
        else:
            return

        fk = fragment.key_tuple
        self.premises_by_fragment_key.setdefault(fk, []).append(premise)

    def make_fixes(self) -> list[OverhangPremise]:
        fixes_made = []
        err_length = self.error_length

        for prem_list in self.premises_by_fragment_key.values():
            prem_count = len(prem_list)

            log.debug(
                f"\n{prem_count} OverhangPremises for {prem_list[0].fragment}:\n"
                + textwrap.indent("".join(f"\n{prem}" for prem in prem_list), "  ")
            )

            if prem_count == 2:
                # To prevent cuts being made which result in Fragments smaller
                # than a Pretext pixel, remove Fragment from the OverlapResult
                # scaffold with the shortest overlap to the bait Fragment.
                frst, scnd = prem_list
                if frst.bait_overlap < err_length and scnd.bait_overlap < err_length:
                    if frst.bait_overlap < scnd.bait_overlap:
                        frst.apply()
                        fixes_made.append(frst)
                    else:
                        scnd.apply()
                        fixes_made.append(scnd)
                    continue

            if prem_count > 1:
                # Can only discard overhanging fragments present in more than
                # one Scaffold, or we would be removing sequence data from
                # the assembly.
                best_to_worst = sorted(
                    prem_list, key=lambda x: x.overhang_error_delta_if_applied
                )
                bst = best_to_worst[0]
                nxt = best_to_worst[1]
                if bst.improves(err_length) and nxt.makes_worse(err_length):
                    bst.apply()  # Remove the overhanging fragment
                    fixes_made.append(bst)

        return fixes_made
