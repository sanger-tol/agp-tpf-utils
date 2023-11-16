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
        self.found_fragments = {}
        self.fragments_found_more_than_once = {}

    def remap_to_input_assembly(self, prtxt_asm, input_asm):
        if not self.bp_per_texel:
            self.bp_per_texel = prtxt_asm.bp_per_texel
        self.find_assembly_overlaps(prtxt_asm, input_asm)
        self.discard_overhanging_fragments()
        self.cut_remaining_overhangs()
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
                else:
                    logging.warn(f"No overlaps found for: {prtxt_frag}")

    def discard_overhanging_fragments(self):
        multi = self.fragments_found_more_than_once

        while multi:
            ovr_resolver = OverhangResolver()
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

    def cut_remaining_overhangs(self):
        multi = self.fragments_found_more_than_once

        for fk, fnd in multi.items():
            self.cut_fragments(fnd)

        self.fragments_found_more_than_once = {}

    def cut_fragments(self, fnd):
        # Make a new Fragment for each region of the fragment found in each
        # OverlapResult
        frgmnt = fnd.fragment
        sub_fragments = [s.trim_fragment(frgmnt) for s in fnd.scaffolds]

        # Check that sub fragments abut each other and do not overlap
        abut_count = 0
        overlap_count = 0
        lgth = len(sub_fragments)
        for i in range(0, lgth):
            frag_a = sub_fragments[i]
            for j in range(i + 1, lgth):
                frag_b = sub_fragments[j]
                if frag_a.abuts(frag_b):
                    abut_count += 1
                if frag_a.overlaps(frag_b):
                    overlap_count += 1
        msg = ""
        if not overlap_count == 0:
            msg += (
                f"Expecting 0 but got {overlap_count} overlaps in new sub fragments\n"
            )
        if not abut_count == lgth - 1:
            msg += f"Execting {lgth -1} abutting sub fragments but got {abut_count}\n"
        if msg:
            msg += "\n" + "\n\n".join(str(s) for s in fnd.scaffolds)
            raise ValueError(msg)

        sub_fragments.sort(key=lambda f: f.start)
        logging.warn(
            f"Fragment {frgmnt} cut into:\n"
            + "".join(f"  {sub}\n" for sub in sub_fragments)
        )

    def log_multi_scaffolds(self):
        multi = self.fragments_found_more_than_once

        for fnd in multi.values():
            ff = fnd.fragment
            logging.warn(
                f"\nFragment {ff} ({ff.length}) found in:\n"
                + "\n".join(
                    (
                        f"{scffld.start_overhang:9d} {scffld.end_overhang:9d}"
                        + f"  {scffld.bait} ({scffld.bait.length})"
                    )
                    for scffld in fnd.scaffolds
                )
            )

    def store_fragments_found(self, scffld):
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


class FoundFragment:
    """
    Little object to store fragments found and the list of Scaffolds it was
    found in.
    """

    __slots__ = "fragment", "scaffolds"

    def __init__(self, fragment):
        self.fragment = fragment
        self.scaffolds = []

    @property
    def scaffold_count(self):
        return len(self.scaffolds)

    def add_scaffold(self, scaffold):
        self.scaffolds.append(scaffold)

    def remove_scaffold(self, scaffold):
        self.scaffolds.remove(scaffold)


class OverhangPremise:
    """
    Stores a "what-if" for removal of a terminal (start or end) Fragment. Used
    to decide which OverlapResult to remove a Fragment from, where the
    Fragment is present in more than one OverlapResult.
    """

    __slots__ = "scaffold", "fragment", "position", "error_increase"

    def __init__(self, scaffold, fragment, position):
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
    Takes in a list of "problem" OverlapResults which share a Fragment.
    Performs one round of comparing OverlapResult pairs, choosing which of
    the two to remove the shared, terminal Fragment from. Returns a list of
    the OverlapPremises which were applied.
    """

    def __init__(self):
        self.premises_by_fragment_key = {}

    def add_overhang_premise(self, fragment, scffld):
        if scffld.rows[0] is fragment:
            premise = OverhangPremise(scffld, fragment, 1)
        elif scffld.rows[-1] is fragment:
            premise = OverhangPremise(scffld, fragment, -1)
        else:
            return

        fk = fragment.key_tuple
        self.premises_by_fragment_key.setdefault(fk, []).append(premise)

    def make_fixes(self):
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
