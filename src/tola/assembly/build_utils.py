"""
Utility objects used by BuildAssembly
"""

import logging
import re
import textwrap

from tola.assembly.fragment import Fragment
from tola.assembly.overlap_result import OverlapResult
from tola.assembly.scaffold import Scaffold


class ScaffoldNamer:
    """
    Tracks naming of chromosomes as Pretext assembly is processed
    """

    def __init__(self, autosome_prefix="RL_"):
        self.autosome_prefix = autosome_prefix
        self.current_scaffold_name = None
        self.current_rank = None
        self.current_haplotype = None
        self.haplotig_n = 0
        self.haplotig_scaffolds = []
        self.target_tags = False
        self.unloc_n = 0
        self.unloc_scaffolds = []

        # Halplotype names stored under their lower case names
        self.haplotype_lc_dict = {}

    OTHER_KNOWN_TAGS = {
        "Contaminant",
        "Cut",
        "Haplotig",
        "Unloc",
    }

    def make_scaffold_name(self, scaffold: Scaffold) -> None:
        """
        Using the tags in the Scaffold from Pretext, work out what the
        chromosome name should be.
        """
        scaffold_name = None
        haplotype = None
        is_painted = False  # Has HiC contacts
        rank = None

        for tag in scaffold.fragment_tags():
            if tag == "Painted":
                is_painted = True
            elif tag == "Target":
                self.target_tags = True
            elif re.fullmatch(r"([A-Z]\d*|[IVX_]+)", tag):
                # This tag looks like a chromosome name
                if scaffold_name and tag != scaffold_name:
                    msg = (
                        f"Found more than one scaffold_name name: '{scaffold_name}'"
                        f" and '{tag}' in scaffold:\n\n{scaffold}"
                    )
                    raise ValueError(msg)
                scaffold_name = tag
                rank = 2
            elif tag not in self.OTHER_KNOWN_TAGS:
                # Any tag that doesn't look like a chromosome name is assumed
                # to be a haplotype, and we only expect to find one within
                # each Pretext Scaffold
                if haplotype:
                    msg = (
                        f"Found both '{haplotype}' and '{tag}', when only one'"
                        f" is expected, in scaffold:\n\n{scaffold}"
                    )
                    raise ValueError(msg)
                else:
                    # The first occurance of a haplotype tag sets its case.
                    # i.e.  "Hap1" will be used if it is seen before "HAP1".
                    haplotype = self.haplotype_lc_dict.setdefault(tag.lower(), tag)

        if not scaffold_name:
            if is_painted:
                scaffold_name = scaffold.name
                if not rank:
                    rank = 1  # Rank for autosomes
            else:
                # Unpainted scaffolds keep the name they have in the input
                # assembly
                scaffold_name = scaffold.rows[0].name
                rank = 3

                # If we don't have a haplotype from a tag, does its name begin
                # with the name of a haplotype?  (This will fail if unplaced
                # contigs from a haplotype appear before the first Scaffold
                # assigned to that haplotype in the Pretext Assembly.)
                if not haplotype and (m := re.match(r"([^_]+)_", scaffold_name)):
                    lc_prefix = m.group(1).lower()
                    haplotype = self.haplotype_lc_dict.get(lc_prefix)

        self.current_scaffold_name = scaffold_name
        self.current_rank = rank
        self.current_haplotype = haplotype
        self.unloc_n = 0
        self.unloc_scaffolds = []

    def label_scaffold(
        self,
        scaffold: Scaffold,
        fragment: Fragment,
        scaffold_tags: set[str],
        original_name: str,
    ) -> None:
        name = self.current_scaffold_name
        if "Contaminant" in fragment.tags:
            scaffold.tag = "Contaminant"
        elif "Haplotig" in fragment.tags:
            name = self.haplotig_name()
            scaffold.tag = "Haplotig"
            self.haplotig_scaffolds.append(scaffold)
        elif "Unloc" in fragment.tags:
            if "Painted" not in scaffold_tags:
                msg = f"Unloc in unpainted scaffold {original_name!r}: {fragment}"
                raise ValueError(msg)
            name = self.unloc_name()
            self.unloc_scaffolds.append(scaffold)
        elif self.target_tags and "Target" not in scaffold_tags:
            scaffold.tag = "Contaminant"

        scaffold.name = name
        scaffold.haplotype = self.current_haplotype
        scaffold.rank = self.current_rank
        scaffold.original_name = original_name

    def haplotig_name(self) -> str:
        self.haplotig_n += 1
        return f"H_{self.haplotig_n}"

    def unloc_name(self) -> str:
        self.unloc_n += 1
        return f"{self.current_scaffold_name}_unloc_{self.unloc_n}"

    def rename_haplotigs_by_size(self) -> None:
        self.rename_by_size(self.haplotig_scaffolds)

    def rename_unlocs_by_size(self) -> None:
        self.rename_by_size(self.unloc_scaffolds)

    def rename_by_size(self, scaffolds: list[Scaffold]) -> None:
        if not scaffolds:
            return
        names = [s.name for s in scaffolds]
        by_size = sorted(scaffolds, key=lambda s: s.length, reverse=True)
        for s, n in zip(by_size, names, strict=True):
            s.name = n


class ChrNamer:
    def __init__(self, chr_prefix="SUPER_"):
        self.chr_prefix = chr_prefix
        self.scaffolds = []
        self.haplotypes_seen = {}
        self.groups = []

    def add_scaffold(self, haplotype, scffld):
        self.haplotypes_seen[haplotype] = True
        self.scaffolds.append((haplotype, scffld))

    def new_group(self):
        grp = {}
        for hap in self.haplotypes_seen:
            grp[hap] = {}
        self.groups.append(grp)
        return grp

    def group_length(self, group, first_haplotype):
        pass

    def name_autosomes(self):
        first_haplotype, *_ = self.haplotypes_seen.keys()
        group = self.new_group()
        last_haplotype = None
        for haplotype, scffld in self.scaffolds:
            orig = scffld.original_name
            if not orig:
                msg = "Missing original_name value in Scaffold:\n{scffld}"
                raise ValueError(msg)
            if haplotype != last_haplotype and group.get(haplotype):
                # New haplotype which is already filled in this group
                group = self.new_group()
            last_haplotype = haplotype
            group[haplotype].setdefault(orig, []).append(scffld)


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

    def apply(self) -> None:
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

    def apply(self) -> None:
        self.scaffold.discard_end()


class OverhangResolver:
    """
    Takes in a list of "problem" OverlapResults which share a Fragment.
    Performs one round of comparing OverlapResult pairs, choosing which of
    the two to remove the shared, terminal Fragment from. Returns a list of
    the OverlapPremises which were applied.
    """

    def __init__(self, error_length=None):
        self.premises_by_fragment_key = {}
        self.error_length = error_length

    def add_overhang_premise(self, fragment: Fragment, scffld: OverlapResult) -> None:
        if scffld.rows[0] is fragment:
            premise = StartOverhangPremise(scffld, fragment)
        elif scffld.rows[-1] is fragment:
            premise = EndOverhangPremise(scffld, fragment)
        else:
            return

        fk = fragment.key_tuple
        self.premises_by_fragment_key.setdefault(fk, []).append(premise)

    def make_fixes(self) -> list[OverlapResult]:
        fixes_made = []
        err_length = self.error_length

        for prem_list in self.premises_by_fragment_key.values():
            prem_count = len(prem_list)

            logging.debug(
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
