import logging
import math
from collections.abc import Iterator

from tola.assembly.assembly import Assembly
from tola.assembly.assembly_stats import AssemblyStats
from tola.assembly.build_utils import (
    ChrNamer,
    FoundFragment,
    OverhangResolver,
    ScaffoldNamer,
)
from tola.assembly.fragment import Fragment
from tola.assembly.gap import Gap
from tola.assembly.indexed_assembly import IndexedAssembly
from tola.assembly.overlap_result import OverlapResult
from tola.assembly.scaffold import Scaffold


class BuildAssembly(Assembly):
    """
    Class for building an Assembly from a Pretext Assembly and the
    IndexedAssembly source. Stores a list of mutable OverlapResults rather
    than Scaffolds, which are fused into Scaffolds by name and returned in
    new Assembly object(s) from the assemblies_with_scaffolds_fused()
    method.
    """

    def __init__(
        self,
        name,
        header=None,
        scaffolds=None,
        default_gap=None,
        bp_per_texel=None,
        autosome_prefix=None,
    ):
        super().__init__(name, header, scaffolds, bp_per_texel)
        self.default_gap = default_gap
        self.found_fragments = {}
        self.fragments_found_more_than_once = {}
        self.scaffold_namer = ScaffoldNamer()
        self.assembly_stats = AssemblyStats()
        if autosome_prefix:
            self.autosome_prefix = autosome_prefix

    @property
    def autosome_prefix(self):
        return self.scaffold_namer.autosome_prefix

    @autosome_prefix.setter
    def autosome_prefix(self, prefix: str):
        self.scaffold_namer.autosome_prefix = prefix
        self.assembly_stats.autosome_prefix = prefix

    @property
    def error_length(self) -> int:
        """
        Expected maximum resolution from bp_per_texel as an integer which is
        guaranteed to be larger than the smallest length from Pretext, even
        if the resolution's floating point value after the decimal point is
        zero.  i.e. 2300.000000 becomes 2301
        """
        return 1 + math.floor(self.bp_per_texel)

    def remap_to_input_assembly(
        self, prtxt_asm: Assembly, input_asm: IndexedAssembly
    ) -> None:
        if not self.bp_per_texel:
            self.bp_per_texel = prtxt_asm.bp_per_texel
        self.assembly_stats.input_assembly = input_asm
        self.find_assembly_overlaps(prtxt_asm, input_asm)
        self.discard_overhanging_fragments()
        self.cut_remaining_overhangs()
        self.scaffold_namer.rename_haplotigs_by_size()
        self.add_missing_scaffolds_from_input(input_asm)

    def find_assembly_overlaps(
        self, prtxt_asm: Assembly, input_asm: IndexedAssembly
    ) -> None:
        logging.info(f"Pretext resolution = {self.bp_per_texel:,.0f} bp per texel\n")
        scaffold_namer = self.scaffold_namer
        err_length = self.error_length
        for prtxt_scffld in prtxt_asm.scaffolds:
            scaffold_namer.make_scaffold_name(prtxt_scffld)
            prtxt_scffld_tags = prtxt_scffld.fragment_tags()
            for prtxt_frag in prtxt_scffld.fragments():
                if found := input_asm.find_overlaps(prtxt_frag):
                    scaffold_namer.label_scaffold(
                        found, prtxt_frag, prtxt_scffld_tags, prtxt_scffld.name
                    )
                    found.trim_large_overhangs(err_length)
                    if found.rows:
                        self.add_scaffold(found)
                        self.store_fragments_found(found)
                else:
                    logging.warning(f"No overlaps found for: {prtxt_frag}")
            scaffold_namer.rename_unlocs_by_size()

    def discard_overhanging_fragments(self) -> None:
        multi = self.fragments_found_more_than_once

        while multi:
            ovr_resolver = OverhangResolver(self.error_length)
            for fnd in multi.values():
                for scffld in fnd.scaffolds:
                    ovr_resolver.add_overhang_premise(fnd.fragment, scffld)
            fixes_made = ovr_resolver.make_fixes()
            if fixes_made:
                for premise in fixes_made:
                    # Remove the Scaffold we fixed
                    fk = premise.fragment.key_tuple
                    if fxd := multi.get(fk):
                        fxd.remove_scaffold(premise.scaffold)
                        if fxd.scaffold_count <= 1:
                            # Fragment is no longer in more than one Scaffold,
                            # so remove it from fragments_found_more_than_once
                            del multi[fk]
            else:
                break

    def cut_remaining_overhangs(self) -> None:
        multi = self.fragments_found_more_than_once

        for fnd in multi.values():
            self.cut_fragments(fnd)

        self.fragments_found_more_than_once = {}

    def cut_fragments(self, fnd: FoundFragment) -> None:
        """
        Make a new Fragment for each region of the fragment found in each
        OverlapResult
        """
        frgmnt = fnd.fragment
        ordered_scaffolds = sorted(
            fnd.scaffolds, key=lambda s: s.fragment_start_if_trimmed(frgmnt)
        )

        sub_fragments = []
        last_i = len(ordered_scaffolds) - 1
        for i, scffld in enumerate(ordered_scaffolds):
            keep_start = i == 0
            keep_end = i == last_i
            sub_fragments.append(scffld.trim_fragment(frgmnt, keep_start, keep_end))
        self.qc_sub_fragments(fnd, sub_fragments)

        self.assembly_stats.cuts += len(sub_fragments) - 1

        logging.info(
            f"Contig:\n  {frgmnt.length:15,d}  {frgmnt}\ncut into:\n"
            + "".join(f"  {sub.length:15,d}  {sub}\n" for sub in sub_fragments)
        )

    def qc_sub_fragments(
        self, fnd: FoundFragment, sub_fragments: list[Fragment]
    ) -> None:
        """
        Check that sub fragments abut each other and do not overlap, and that
        no sequence from the cut fragment has been lost.
        """
        abut_count = 0
        overlap_count = 0
        pairs_with_gaps = []
        srtd_frags = sorted(sub_fragments, key=lambda frag: (frag.start, frag.end))
        for i, frag_a in enumerate(srtd_frags[:-1]):
            frag_b = srtd_frags[i + 1]
            if frag_a.abuts(frag_b):
                abut_count += 1
            if frag_a.overlaps(frag_b):
                overlap_count += 1
            if g := frag_a.gap_between(frag_b):
                pairs_with_gaps.append((frag_a, frag_b, g))

        sub_frags_length = sum(f.length for f in sub_fragments)

        msg = ""
        if fnd.fragment.length != sub_frags_length:
            msg += (
                f"Sum of fragment lengths {sub_frags_length:_d} does not"
                f" match orginal fragment length {fnd.fragment.length:_d}\n"
            )
        if overlap_count != 0:
            msg += (
                f"Expecting 0 but got {overlap_count} overlaps in new sub fragments\n"
            )
        lgth = len(sub_fragments)
        if abut_count != lgth - 1:
            msg += f"Expecting {lgth - 1} abutting sub fragments but got {abut_count}\n"
        for frag_a, frag_b, g in pairs_with_gaps:
            pixels = g / self.bp_per_texel
            msg += (
                f"Gap of length {g} ({pixels:.1f} pixels)"
                f" between:\n  {frag_a}\nand:\n  {frag_b}\n"
            )

        if msg:
            msg += "\n" + "\n\n".join(str(s) for s in fnd.scaffolds)
            raise ValueError(msg)

    def log_multi_scaffolds(self) -> None:
        multi = self.fragments_found_more_than_once

        for fnd in multi.values():
            ff = fnd.fragment
            logging.warning(
                f"\nFragment {ff} ({ff.length}) found in:\n"
                + "\n".join(
                    (
                        f"{scffld.start_overhang:9d} {scffld.end_overhang:9d}"
                        f"  {scffld.bait} ({scffld.bait.length})"
                    )
                    for scffld in fnd.scaffolds
                )
            )

    def store_fragments_found(self, scffld: Scaffold) -> None:
        store = self.found_fragments
        multi = self.fragments_found_more_than_once
        for ff in scffld.fragments():
            ff_tuple = ff.key_tuple
            if fnd := store.get(ff_tuple):
                # Already have it, so record that we've found it more than
                # once
                multi[ff_tuple] = fnd
            else:
                fnd = FoundFragment(ff)
                store[ff_tuple] = fnd
            fnd.add_scaffold(scffld)

    def add_missing_scaffolds_from_input(self, input_asm: Assembly) -> None:
        scaffold_namer = self.scaffold_namer
        found_frags = self.found_fragments
        for scffld in input_asm.scaffolds:
            new_scffld = None
            last_added_i = None
            for i, frag in scffld.idx_fragments():
                if not found_frags.get(frag.key_tuple):
                    if not new_scffld:
                        new_scffld = Scaffold(scffld.name)
                        new_scffld.rank = 3
                    if last_added_i is not None and last_added_i != i - 1:
                        # Last added row was not the previous row in the
                        # scaffold
                        prev_row = scffld.rows[i - 1]
                        if isinstance(prev_row, Gap):
                            new_scffld.add_row(prev_row)
                        else:
                            new_scffld.add_row(self.default_gap)
                    new_scffld.add_row(frag)
                    last_added_i = i

            if new_scffld:
                scaffold_namer.make_scaffold_name(new_scffld)
                if (
                    scaffold_namer.target_tags
                    and "Target" not in scffld.fragment_tags()
                ):
                    new_scffld.tag = "Contaminant"
                new_scffld.haplotype = scaffold_namer.current_haplotype
                self.add_scaffold(new_scffld)

    def assemblies_with_scaffolds_fused(self) -> list[Assembly]:
        chr_namer = ChrNamer(chr_prefix=self.autosome_prefix)
        assemblies = {}
        for scffld in self.scaffolds_fused_by_name():
            if tag := scffld.tag:
                asm_key = tag
                asm_name = f"{self.name}.{tag.lower()}s"
            elif hap := scffld.haplotype:
                asm_key = hap
                asm_name = f"{self.name}.{hap.lower()}"
            else:
                asm_key = None
                asm_name = self.name
            new_asm = assemblies.setdefault(asm_key, Assembly(asm_name))
            new_asm.add_scaffold(scffld)
            if scffld.rank == 1:
                # Add autosome to the ChrNamer
                chr_namer.add_scaffold(asm_key, scffld)
            elif scffld.rank == 2:
                chr_namer.add_chr_prefix(scffld)

        # ChrNamer names autosome chromosomes by size
        chr_namer.name_chromosomes()

        for asm in assemblies.values():
            # Sort scaffolds by name
            asm.smart_sort_scaffolds()

        self.assembly_stats.make_stats(assemblies)

        return assemblies

    def scaffolds_fused_by_name(self) -> Iterator[Scaffold]:
        gap = self.default_gap
        hap_name_scaffold = {}
        for scffld in self.scaffolds:
            if not scffld.rows:
                # discard_overhanging_fragments() may have removed the only
                # row from an OverlapResult
                continue

            build_scffld = hap_name_scaffold.setdefault(
                (scffld.haplotype, scffld.name),
                Scaffold(
                    scffld.name,
                    tag=scffld.tag,
                    haplotype=scffld.haplotype,
                    rank=scffld.rank,
                    original_name=scffld.original_name,
                ),
            )
            if isinstance(scffld, OverlapResult):
                build_scffld.append_scaffold(scffld.to_scaffold(), gap)
            else:
                build_scffld.append_scaffold(scffld)

        for scffld in hap_name_scaffold.values():
            yield scffld
