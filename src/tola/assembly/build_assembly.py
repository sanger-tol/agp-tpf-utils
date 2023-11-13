import logging

from tola.assembly.assembly import Assembly
from tola.assembly.gap import Gap
from tola.assembly.overlap_result import OverlapResult
from tola.assembly.scaffold import Scaffold


class BuildAssembly(Assembly):
    """
    Class for building an Assembly from a Pretext Assembly and the
    IndexedAssembly source. Stores a list of mutable OverlapResults rather
    than Scaffolds, which are fused into Scaffolds by name and returned in a
    new Assembly object from the assembly_with_scaffolds_fused() method.
    """

    def __init__(
        self, name, header=None, scaffolds=None, default_gap=None, bp_per_texel=None
    ):
        super().__init__(name, header, scaffolds, bp_per_texel)
        self.default_gap = default_gap
        self.problem_scaffolds = []
        self.found_fragments = {}

    def remap_to_input_assembly(self, prtxt_asm, input_asm):
        if not self.bp_per_texel:
            self.bp_per_texel = prtxt_asm.bp_per_texel
        self.find_assembly_overlaps(prtxt_asm, input_asm)
        self.discard_overhanging_fragments(self.bp_per_texel)
        self.add_missing_scaffolds_from_input(input_asm)

    def find_assembly_overlaps(self, prtxt_asm, input_asm):
        scffld_n = 0
        bp_per_texel = self.bp_per_texel
        for prtxt_scffld in prtxt_asm.scaffolds:
            scffld_n += 1
            scffld_name = f"R{scffld_n}"
            for prtxt_frag in prtxt_scffld.fragments():
                if found := input_asm.find_overlaps(prtxt_frag):
                    found.name = scffld_name
                    self.add_scaffold(found)
                    found.trim_large_overhangs(bp_per_texel)
                    self.store_fragments_found(found)

                    if found.has_problem_overhang(bp_per_texel):
                        self.problem_scaffolds.append(found)
                else:
                    logging.warn(f"No overlaps found for: {prtxt_frag}")

    def discard_overhanging_fragments(self, bp_per_texel):
        problems = self.problem_scaffolds
        while prob_count := len(problems):
            ovr_resolver = OverhangResolver(bp_per_texel, problems)
            fix_count, problems = ovr_resolver.make_fixes()
            logging.debug(
                f"Discarded {fix_count} overhanging fragments"
                f" in {prob_count} problem scaffolds"
            )
            if not fix_count:
                break
        self.problem_scaffolds = problems

    def log_problem_scaffolds(self):
        bp_per_texel = self.bp_per_texel
        if probs := self.problem_scaffolds:
            for scffld in probs:
                # Log problem regions
                if logging.root.level < logging.INFO:
                    logging.debug(scffld)
                else:
                    err = scffld.length_error_in_texels(bp_per_texel)
                    logging.info(
                        f"{scffld.start_overhang:9d} {scffld.end_overhang:9d}"
                        + f"  {err:6.2f} pixels  {scffld.bait}",
                    )

    def store_fragments_found(self, scffld):
        store = self.found_fragments
        for ff in scffld.fragments():
            ff_tuple = ff.key_tuple
            store[ff_tuple] = 1 + store.get(ff_tuple, 0)

    def add_missing_scaffolds_from_input(self, input_asm):
        found_frags = self.found_fragments
        for scffld in input_asm.scaffolds:
            new_scffld = None
            last_added_i = None
            for i, frag in scffld.idx_fragments():
                if not found_frags.get(frag.key_tuple):
                    if not new_scffld:
                        new_scffld = Scaffold(scffld.name)
                    if last_added_i is not None and not last_added_i == i - 1:
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
                self.add_scaffold(new_scffld)

    def assembly_with_scaffolds_fused(self):
        new_asm = Assembly(self.name)
        for scffld in self.scaffolds_fused_by_name():
            new_asm.add_scaffold(scffld)
        return new_asm

    def scaffolds_fused_by_name(self):
        gap = self.default_gap
        new_scffld = None
        current_name = ""
        for scffld in self.scaffolds:
            if scffld.name != current_name:
                if new_scffld:
                    yield new_scffld
                current_name = scffld.name
                new_scffld = Scaffold(scffld.name)
            if isinstance(scffld, OverlapResult):
                new_scffld.append_scaffold(scffld.to_scaffold(), gap)
            else:
                new_scffld.append_scaffold(scffld)

        if new_scffld:
            yield new_scffld


class OverhangPremise:
    """
    Stores a "what-if" for removal of a terminal (start or end) Fragment. Used
    to decide which OverlapResult to remove a Fragment from, where the
    Fragment is present in more than one OverlapResult.
    """

    __slots__ = "scaffold", "position", "error_increase"

    def __init__(self, scaffold, position):
        self.scaffold = scaffold
        if position == 1:
            self.error_increase = scaffold.error_increase_if_start_removed()
        elif position == -1:
            self.error_increase = scaffold.error_increase_if_end_removed()
        else:
            msg = f"position must be '1' (start) or '-1' (end) not '{position}'"
            raise ValueError(msg)
        self.position = position

    @property
    def fragment(self):
        return self.scaffold.rows[0 if self.position == 1 else -1]

    @property
    def improves(self):
        return True if self.error_increase < 0 else False

    @property
    def makes_worse(self):
        return True if self.error_increase > 0 else False

    def apply(self):
        if self.position == 1:
            self.scaffold.discard_start()
        else:
            self.scaffold.discard_end()


class OverhangResolver:
    """
    Takes in a list of "problem" OverlapResults, i.e. with start or end
    overhangs longer than the expected bp_per_pixel inaccuracy. Performs one
    round of comparing OverlapResult pairs, choosing which of the two to
    remove the shared, terminal Frament from. Returns a count of the number
    of fixes applied, and a list of the remaining "problem" OverlapResults.
    """

    def __init__(self, bp_per_texel, scaffolds=None):
        self.bp_per_texel = bp_per_texel
        self.premises_by_fragment_key = {}
        if scaffolds:
            self.add_scaffolds(scaffolds)

    def add_scaffolds(self, scaffolds):
        for scffld in scaffolds:
            if scffld.start_overhang > self.bp_per_texel:
                self.add_overhang_premise(scffld, 1)
            if scffld.end_overhang > self.bp_per_texel:
                self.add_overhang_premise(scffld, -1)

    def add_overhang_premise(self, scffld, position):
        premise = OverhangPremise(scffld, position)
        fk = premise.fragment.key_tuple
        self.premises_by_fragment_key.setdefault(fk, []).append(premise)

    def make_fixes(self):
        fixes_made = 0
        still_a_problem = {}
        for prem_list in self.premises_by_fragment_key.values():
            # Can only discard overhanging fragments present in more than one
            # Scaffold, or we would be removing sequence from the assembly.
            if len(prem_list) > 1:
                best_to_worst = sorted(prem_list, key=lambda x: x.error_increase)
                bst = best_to_worst[0]
                nxt = best_to_worst[1]
                if bst.improves and nxt.makes_worse:
                    bst.apply()  # Remove the overhanging fragment
                    fixes_made += 1
            for premise in prem_list:
                scffld = premise.scaffold
                if scffld.has_problem_overhang(self.bp_per_texel):
                    # Store Scaffolds keyed on bait so that Scaffolds with two
                    # bad overhangs are not stored twice.
                    still_a_problem[scffld.bait.key_tuple] = scffld

        return fixes_made, list(still_a_problem.values())
