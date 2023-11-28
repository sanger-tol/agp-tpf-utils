import logging

from tola.assembly.assembly import Assembly
from tola.assembly.indexed_assembly import IndexedAssembly


class AssemblyStats:
    def __init__(
        self, autosome_prefix: str, assembly: Assembly, input_assembly: IndexedAssembly
    ) -> None:
        self.autosome_prefix = autosome_prefix
        self.assembly_name = assembly.name
        self.breaks = 0
        self.joins = 0
        self._make_stats(assembly, input_assembly)

    def _make_stats(self, asm: Assembly, input_asm: IndexedAssembly) -> None:
        self._count_breaks_and_joins_at_scaffold_ends(asm, input_asm)
        self._count_breaks_and_joins_within_scaffolds(asm)

    def _count_breaks_and_joins_at_scaffold_ends(
        self, asm: Assembly, input_asm: IndexedAssembly
    ) -> None:
        breaks = 0

        for scffld in asm.scaffolds:
            first = scffld.rows[0]
            last = scffld.rows[-1]
            if first.start != 1:
                breaks += 1
            source_scaffold = input_asm.scaffold_by_name(last.name)
            if last.end != source_scaffold.length:
                breaks += 1

        self.breaks += breaks

    def _count_breaks_and_joins_within_scaffolds(self, asm: Assembly) -> None:
        changes = 0

        for scffld in asm.scaffolds:
            idx_frags = tuple(scffld.idx_fragments())

            for idx in range(1, len(idx_frags)):
                i, a_frag = idx_frags[idx - 1]
                j, b_frag = idx_frags[idx]
                if a_frag.name != b_frag.name:
                    # Fragments come from different scaffolds
                    changes += 1
                else:
                    # Calculate gap length from rows between i and j
                    gap_length = sum(g.length for g in scffld.rows[i + 1 : j])

                    # If the strands match, check if the fragments follow each other
                    if a_frag.strand == 1 and b_frag.strand == 1:
                        if a_frag.end + gap_length + 1 != b_frag.start:
                            changes += 1
                    elif a_frag.strand == -1 and b_frag.strand == -1:
                        if b_frag.end + gap_length + 1 != a_frag.start:
                            changes += 1
                    else:
                        # Strands do not match, so one has been flipped
                        changes += 1

        self.breaks += changes
        self.joins += changes

    def log_stats(self):
        logging.info(f"Assembly {self.assembly_name} has {self.breaks} breaks and {self.joins} joins")
