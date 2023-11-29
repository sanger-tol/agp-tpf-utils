import logging

from tola.assembly.assembly import Assembly
from tola.assembly.indexed_assembly import IndexedAssembly


class AssemblyStats:
    def __init__(self, autosome_prefix: str = None) -> None:
        self.autosome_prefix = autosome_prefix
        self.input_assembly = None
        self.cuts = 0
        self.breaks = 0
        self.joins = 0

    def make_stats(self, output_assemblies: list[Assembly]) -> None:
        input_set = self.input_assembly.fragment_junction_set()
        output_set = set()
        for asm in output_assemblies:
            output_set.update(asm.fragment_junction_set())

        # Breaks are junctions in the input that are not in the output
        self.breaks = len(input_set - output_set)

        # Joins are junctions in the output that were not in the input
        self.joins = len(output_set - input_set)

    def log_stats(self):
        cut_plural = "cut in a contig" if self.cuts == 1 else "cuts in contigs"
        break_plural = "one break at a gap" if self.breaks == 1 else "breaks at gaps"
        join_plural = "join" if self.joins == 1 else "joins"
        logging.info(
            f"Curation made {self.cuts} {cut_plural}, {self.breaks}"
            f" {break_plural} and {self.joins} {join_plural}"
        )
