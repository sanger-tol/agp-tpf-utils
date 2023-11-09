import logging

from tola.assembly.assembly import Assembly
from tola.assembly.scaffold import Scaffold


class BuildAssembly(Assembly):
    def __init__(self, name, header=None, scaffolds=None, default_gap=None):
        super().__init__(name, header, scaffolds)
        self.default_gap = default_gap
        self.problem_scaffolds = []
        self.found_fragments = {}

    def remap_to_input_assembly(self, prtxt_asm, input_asm):
        self.find_assembly_overlaps(prtxt_asm, input_asm)
        while True:
            if not self.fix_problem_scaffolds(prtxt_asm.bp_per_texel):
                break

    def find_assembly_overlaps(self, prtxt_asm, input_asm):
        scffld_n = 0
        bp_per_texel = prtxt_asm.bp_per_texel
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

    def fix_problem_scaffolds(self, bp_per_texel):
        overhangs = {}
        START = 0
        END = -1

        def add_to_overhang(idx, scffld):
            frag = scffld.rows[idx]
            k = frag.frag_tuple
            i_scaffolds = overhangs.get(k, [])
            discrep = (
                scffld.error_increase_if_start_removed()
                if idx == START
                else scffld.error_increase_if_end_removed()
            )

            ### Turn results into an object to improve readability
            i_scaffolds.append((discrep, idx, scffld))
            overhangs[k] = i_scaffolds

        for scffld in self.problem_scaffolds:
            if scffld.start_overhang > bp_per_texel:
                add_to_overhang(START, scffld)
            if scffld.end_overhang > bp_per_texel:
                add_to_overhang(END, scffld)

        fix_made = False
        still_a_problem = {}
        for i_scffld_list in overhangs.values():
            if len(i_scffld_list) > 1:
                best_to_worst = sorted(i_scffld_list)
                bst = best_to_worst[0]
                nxt = best_to_worst[1]
                ### May get three or more if a large contig is cut into three ###

                # If removing the Fragment would improve "bst" and make "nxt"
                # worse, then remove it from "bst"
                if bst[0] < 0 and nxt[0] > 0:
                    if bst[1] == START:
                        bst[2].discard_start()
                    else:
                        bst[2].discard_end()
                    fix_made = True

            for ele in i_scffld_list:
                scffld = ele[2]
                if scffld.has_problem_overhang(bp_per_texel):
                    still_a_problem[scffld.bait.frag_tuple] = scffld

        self.problem_scaffolds = still_a_problem.values()
        return fix_made

    def log_problem_scaffolds(self):
        if probs := self.problem_scaffolds:
            for scffld in probs:
                # Log problem regions
                if logging.root.level < logging.INFO:
                    logging.debug(scffld)
                else:
                    logging.info(
                        f"overhangs: {scffld.start_overhang:9d} {scffld.end_overhang:9d}"
                        + f"  {scffld.bait}",
                    )

    def store_fragments_found(self, scffld):
        store = self.found_fragments
        for ff in scffld.fragments():
            ff_tuple = ff.frag_tuple
            store[ff_tuple] = 1 + store.get(ff_tuple, 0)

    def assembly_with_scaffolds_fused(self):
        new_asm = Assembly(self.name)
        for scffld in self.scaffolds_fused_by_name():
            new_asm.add_scaffold(scffld)
        return new_asm

    def scaffolds_fused_by_name(self):
        gap = self.default_gap
        new_scffld = None
        current_name = ""
        for ovr_res in self.scaffolds:
            if ovr_res.name != current_name:
                if new_scffld:
                    yield new_scffld
                current_name = ovr_res.name
                new_scffld = Scaffold(ovr_res.name)
            new_scffld.append_scaffold(ovr_res.to_scaffold(), gap)

        if new_scffld:
            yield new_scffld
