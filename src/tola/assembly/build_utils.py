"""
Utility objects used by BuildAssembly
"""

import re

from tola.assembly.fragment import Fragment
from tola.assembly.overlap_result import OverlapResult
from tola.assembly.scaffold import Scaffold


class ChrNamer:
    """
    Tracks naming of chromosomes as Pretext assembly is processed
    """

    def __init__(self):
        self.chr_name_n = 0
        self.current_chr_name = None
        self.haplotig_n = 0
        self.current_haplotig = None
        self.haplotig_scaffolds = []
        self.unloc_n = 0
        self.unloc_scaffolds = []

    def make_chr_name(self, scaffold: Scaffold) -> str:
        tag_set = scaffold.fragment_tags()

        chr_name = None
        is_painted = False  # Has HiC contacts
        for tags in tag_set:
            for t in tags:
                if t == "Painted":
                    is_painted = True
                elif m := re.match(r"[A-Z]\d*$", t):
                    cn = m.group(0)
                    if chr_name and cn != chr_name:
                        msg = (
                            f"Found more than one chr_name name: '{chr_name}' and '{cn}'"
                            f" in scaffold:\n\n{scaffold}"
                        )
                        raise ValueError(msg)
                    chr_name = cn

        if not chr_name:
            if is_painted:
                self.chr_name_n += 1
                chr_name = f"RL_{self.chr_name_n}"
            else:
                chr_name = scaffold.rows[0].name

        self.current_chr_name = chr_name
        self.unloc_n = 0
        self.unloc_scaffolds = []

        return chr_name

    def unloc_name(self) -> str:
        self.unloc_n += 1
        return f"{self.current_chr_name}_unloc_{self.unloc_n}"

    def haplotig_name(self) -> str:
        self.haplotig_n += 1
        return f"H_{self.haplotig_n}"

    def name_scaffold(self, scaffold: Scaffold, fragment: Fragment) -> None:
        name = None
        if "Unloc" in fragment.tags:
            name = self.unloc_name()
            self.unloc_scaffolds.append(scaffold)
        elif "Haplotig" in fragment.tags:
            name = self.haplotig_name()
            self.haplotig_scaffolds.append(scaffold)
        else:
            name = self.current_chr_name
        scaffold.name = name

    def rename_unlocs_by_size(self) -> None:
        self.rename_by_size(self.unloc_scaffolds)

    def rename_haplotigs_by_size(self) -> None:
        self.rename_by_size(self.haplotig_scaffolds)

    def rename_by_size(self, scaffolds: list[Scaffold]) -> None:
        if not scaffolds:
            return
        names = [s.name for s in scaffolds]
        by_size = sorted(scaffolds, key=lambda s: s.length, reverse=True)
        for s, n in zip(by_size, names, strict=True):
            s.name = n


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

    def add_scaffold(self, scaffold: Scaffold) -> None:
        self.scaffolds.append(scaffold)

    def remove_scaffold(self, scaffold: Scaffold) -> None:
        self.scaffolds.remove(scaffold)


class OverhangPremise:
    """
    Stores a "what-if" for removal of a terminal (start or end) Fragment. Used
    to decide which OverlapResult to remove a Fragment from, where the
    Fragment is present in more than one OverlapResult.
    """

    __slots__ = "scaffold", "fragment", "position", "error_increase"

    def __init__(self, scaffold: OverlapResult, fragment: Fragment, position: int):
        self.scaffold = scaffold
        self.fragment = fragment
        if position == 1:
            self.error_increase = scaffold.error_increase_if_start_removed()
        elif position == -1:
            self.error_increase = scaffold.error_increase_if_end_removed()
        else:
            msg = f"position must be '1' (start) or '-1' (end) not '{position}'"
            raise ValueError(msg)
        self.position = position

    @property
    def improves(self) -> bool:
        if len(self.scaffold.rows) == 1:
            return False
        return True if self.error_increase < 0 else False

    @property
    def makes_worse(self) -> bool:
        if len(self.scaffold.rows) == 1:
            return True
        return True if self.error_increase > 0 else False

    def apply(self) -> None:
        if self.position == 1:
            self.scaffold.discard_start()
        else:
            self.scaffold.discard_end()


class OverhangResolver:
    """
    Takes in a list of "problem" OverlapResults which share a Fragment.
    Performs one round of comparing OverlapResult pairs, choosing which of
    the two to remove the shared, terminal Fragment from. Returns a list of
    the OverlapPremises which were applied.
    """

    def __init__(self):
        self.premises_by_fragment_key = {}

    def add_overhang_premise(self, fragment: Fragment, scffld: OverlapResult) -> None:
        if scffld.rows[0] is fragment:
            premise = OverhangPremise(scffld, fragment, 1)
        elif scffld.rows[-1] is fragment:
            premise = OverhangPremise(scffld, fragment, -1)
        else:
            return

        fk = fragment.key_tuple
        self.premises_by_fragment_key.setdefault(fk, []).append(premise)

    def make_fixes(self) -> list[OverlapResult]:
        fixes_made = []
        for prem_list in self.premises_by_fragment_key.values():
            # Can only discard overhanging fragments present in more than one
            # Scaffold, or we would be removing sequence from the assembly.
            best_to_worst = sorted(prem_list, key=lambda x: x.error_increase)
            if len(prem_list) > 1:
                bst = best_to_worst[0]
                nxt = best_to_worst[1]
                if bst.improves and nxt.makes_worse:
                    bst.apply()  # Remove the overhanging fragment
                    fixes_made.append(bst)

        return fixes_made
