import logging
import math

from tola.assembly.assembly import Assembly, AssemblyDict
from tola.assembly.assembly_stats import AssemblyStats
from tola.assembly.build_utils import FoundFragment, OverhangResolver
from tola.assembly.fragment import Fragment
from tola.assembly.gap import Gap
from tola.assembly.indexed_assembly import IndexedAssembly
from tola.assembly.naming_utils import ChrNamer, ScaffoldNamer
from tola.assembly.overlap_result import OverlapResult
from tola.assembly.scaffold import Scaffold

log = logging.getLogger(__name__)


class LongScaffoldCuttingError(Exception):
    """
    Error when cutting a `Scaffold` longer than
    `BuildAssembly.max_contig_length`.
    """


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
        max_contig_length=2_000_000_000,
    ):
        super().__init__(name, header, scaffolds, bp_per_texel)
        self.default_gap = Gap(200, "scaffold") if default_gap is None else default_gap
        self.found_fragments = {}
        self.fragments_found_more_than_once = {}
        self.scaffold_namer: ScaffoldNamer = ScaffoldNamer()
        self.assembly_stats: AssemblyStats = AssemblyStats()
        if autosome_prefix:
            self.autosome_prefix = autosome_prefix
        self.max_contig_length = max_contig_length

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
        self.scaffold_namer.rename_unlocs_by_size()
        self.scaffold_namer.rename_haplotigs_by_size()
        self.add_missing_scaffolds_from_input(input_asm)

    def find_assembly_overlaps(
        self, prtxt_asm: Assembly, input_asm: IndexedAssembly
    ) -> None:
        log.info(f"Pretext resolution = {self.bp_per_texel:,.0f} bp per texel\n")
        scaffold_namer = self.scaffold_namer
        err_length = self.error_length
        for prtxt_scffld in prtxt_asm.scaffolds:
            prtxt_scffld_tags = prtxt_scffld.fragment_tags()
            scaffold_namer.make_scaffold_name(prtxt_scffld, prtxt_scffld_tags)
            for prtxt_frag in prtxt_scffld.fragments():
                if found := input_asm.find_overlaps(prtxt_frag):
                    scaffold_namer.label_scaffold(
                        scaffold=found,
                        fragment=prtxt_frag,
                        scaffold_tags=prtxt_scffld_tags,
                        original_name=prtxt_scffld.name,
                    )
                    found.trim_large_overhangs(err_length)
                    if found.rows:
                        self.add_scaffold(found)
                        self.store_fragments_found(found)
                else:
                    log.warning(f"No overlaps found for: {prtxt_frag}")

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

        log.info(
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
            log.warning(
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

    def assembly_with_scaffolds_in_map_order(self) -> AssemblyDict:
        scaffolds, _ = self.__build_name_and_sort_assemblies()
        return {None: Assembly("Pretext", scaffolds=scaffolds)}

    def assemblies_with_scaffolds_fused(self) -> AssemblyDict:
        _, assemblies = self.__build_name_and_sort_assemblies()
        return assemblies

    def __build_name_and_sort_assemblies(
        self,
    ) -> tuple[list[Scaffold], AssemblyDict]:
        chr_namer = ChrNamer(chr_prefix=self.autosome_prefix)

        scaffolds = self.scaffolds_fused_by_name()

        assemblies = {}
        for scffld in scaffolds:
            curated = True
            hap = None
            if tag := scffld.tag:
                # Set assembly name from tag, which will be one of:
                #   Contaminant | FalseDuplicate | Haplotig
                curated = False
                asm_key = tag
            elif hap := scffld.haplotype:
                asm_key = hap
            else:
                asm_key = None

            if not (new_asm := assemblies.get(asm_key)):
                new_asm = Assembly(self.name, curated=curated)
                assemblies[asm_key] = new_asm

            cut_scaffolds: list[Scaffold]
            if self.max_contig_length is not None:
                cut_scaffolds = self.cut_scaffold_if_too_long(scffld)
            else:
                cut_scaffolds = [scffld]

            for cut_scffld in cut_scaffolds:
                new_asm.add_scaffold(cut_scffld)

                if cut_scffld.rank == 1:
                    # Add autosome to the ChrNamer
                    chr_namer.add_scaffold(asm_key, cut_scffld)
                elif cut_scffld.rank == 2:
                    chr_namer.add_chr_prefix(cut_scffld, hap)
                elif hap is not None:
                    chr_namer.add_haplotype_prefix(cut_scffld, hap)

        # ChrNamer names autosome chromosomes by size
        chr_namer.name_chromosomes()

        for asm in assemblies.values():
            # Sort scaffolds by name
            asm.smart_sort_scaffolds()

        self.assembly_stats.make_stats(assemblies)

        return scaffolds, assemblies

    def cut_scaffold_if_too_long(self, scffld: Scaffold) -> list[Scaffold]:
        whole = scffld.length
        pieces = math.ceil(whole / self.max_contig_length)

        if pieces == 1:
            return [scffld]

        # If, for example, we need to cut the scaffold into 3 pieces, this
        # loop will take the first 1/3 off the scaffold, then 1/2 of what's
        # remaining.
        cut_parts = []
        to_cut = scffld
        for div in range(pieces, 1, -1):
            cut_at = whole // div
            gap_i = self.index_of_nearest_gap_to_ideal_cut_site(to_cut, cut_at)
            rows = to_cut.rows

            # First part is everything up to, but not including, the gap
            cut = to_cut.clone_empty()
            cut.rows = rows[:gap_i]
            cut_parts.append(cut)

            # Second part is everything after the gap
            to_cut = to_cut.clone_empty()
            to_cut.rows = rows[gap_i + 1 :]
        cut_parts.append(to_cut)

        # Add suffix "_1", "_2" etc... to cut scaffolds
        for i, part in enumerate(cut_parts):
            part.name = f"{part.name}_{i + 1}"

        # Format report of cuts made
        whole_str = f"{whole:,d}"
        wl = len(whole_str)
        nl = len(scffld.name) + 4
        log.info(
            f"Cut {scffld.name:<{nl}}  {whole:{wl},d} bp (including gaps) into:\n"
            + "".join(
                [f"    {x.name:<{nl}}  {x.length:{wl},d} bp\n" for x in cut_parts]
            )
        )

        for part in cut_parts:
            if part.length > self.max_contig_length:
                msg = (
                    f"Scaffold '{part.name}' length = {part.length:,d} is longer"
                    f" than max contig length {self.max_contig_length:,d}"
                )
                raise LongScaffoldCuttingError(msg)

        return cut_parts

    def index_of_nearest_gap_to_ideal_cut_site(self, to_cut: Scaffold, cut_at: int):
        # Make a temporary `IndexedAssembly` to efficiently search for an row
        # which overlaps the cut coordinate.
        idx_asm = IndexedAssembly(
            f"Temporary Assembly for cutting '{to_cut.name}' at {cut_at:_d}",
            scaffolds=[to_cut],
        )

        # Find the row which overlaps the ideal cut site
        ovr_i_j = idx_asm.overlapping_indices_by_scaffold_start_end(
            to_cut, cut_at, cut_at
        )
        if not ovr_i_j:
            msg = (
                f"Failed to find an element at {cut_at:,d}"
                f" within '{to_cut.name}' of length {to_cut.length:,d}"
            )
            raise LongScaffoldCuttingError(msg)

        ovr_i = ovr_i_j[0]
        rows = to_cut.rows
        ele = rows[ovr_i]
        if not isinstance(ele, Gap):
            # This isn't a gap, so we need to find the nearest
            gap_i_before = None
            for i in range(ovr_i, 0, -1):
                if isinstance(rows[i], Gap):
                    gap_i_before = i
                    break

            gap_i_after = None
            for i in range(ovr_i, len(rows)):
                if isinstance(rows[i], Gap):
                    gap_i_after = i
                    break

            if gap_i_before is None and gap_i_after is None:
                # If this ever happens we could implement cutting here
                msg = (
                    "Failed to find a gap before or"
                    f" after {cut_at:,d} in '{to_cut.name}'\n"
                    "Maybe implement cutting?"
                )
                raise LongScaffoldCuttingError(msg)

            if gap_i_before is None:
                ovr_i = gap_i_after
            elif gap_i_after is None:
                ovr_i = gap_i_before
            else:
                length_before = (
                    idx_asm.start_end_of_row(to_cut.name, gap_i_before)[0] - 1
                )
                length_after = idx_asm.start_end_of_row(to_cut.name, gap_i_after)[0] - 1

                # Choose the gap before or after, whichever is nearest to the
                # ideal cut point.
                ovr_i = (
                    gap_i_before
                    if abs(cut_at - length_before) < abs(length_after - cut_at)
                    else gap_i_after
                )

        return ovr_i

    def scaffolds_fused_by_name(self) -> list[Scaffold]:
        gap = self.default_gap
        hap_name_scaffold: dict[tuple[str | None, str], Scaffold] = {}
        for scffld in self.scaffolds:
            if not scffld.rows:
                # discard_overhanging_fragments() may have removed the only
                # row from an OverlapResult
                continue

            idx = scffld.haplotype, scffld.name
            build_scffld = hap_name_scaffold.get(idx)
            if not build_scffld:
                hap_name_scaffold[idx] = build_scffld = scffld.clone_empty()

            if isinstance(scffld, OverlapResult):
                build_scffld.append_scaffold(scffld.to_scaffold(), gap)
            else:
                build_scffld.append_scaffold(scffld)

        return list(hap_name_scaffold.values())
