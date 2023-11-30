import logging

from tola.assembly.assembly import Assembly
from tola.assembly.indexed_assembly import IndexedAssembly
from tola.assembly.scaffold import Scaffold


class AssemblyStats:
    def __init__(self, autosome_prefix: str = "RL_") -> None:
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

    def log_curation_stats(self):
        cut_plural = "cut in a contig" if self.cuts == 1 else "cuts in contigs"
        break_plural = "break at a gap" if self.breaks == 1 else "breaks at gaps"
        join_plural = "join" if self.joins == 1 else "joins"
        logging.info(
            f"Curation made {self.cuts} {cut_plural}, {self.breaks}"
            f" {break_plural} and {self.joins} {join_plural}"
        )

    def log_assembly_chromosomes(self, asm: Assembly):
        autosome_prefix = self.autosome_prefix
        ranked_scaffolds = {}
        for scffld in asm.scaffolds:
            rank = scffld.rank(autosome_prefix)
            ranked_scaffolds.setdefault(rank, []).append(scffld)

        logging.info(f"\n{asm.name}")
        logging.info(f"    {asm.fragments_length:15,d}  bp sequence (minus gaps)")
        is_main = False  # Only show rank headings for main assemblies
        for rank, scaffolds in ranked_scaffolds.items():
            rank_i, rank_name = rank
            if rank_i == 0:
                is_main = True
            total = sum(s.fragments_length for s in scaffolds)
            if is_main:
                logging.info(f"  {rank_name}:")

            if rank_i != 2:
                scaffolds = self.merge_unlocs(scaffolds)
            logging.info(f"    n = {len(scaffolds)}")

            if rank_i == 1:
                for scffld in scaffolds:
                    self.log_scaffold_length(scffld)
            else:
                scaffolds = sorted(
                    scaffolds, key=lambda s: s.fragments_length, reverse=True
                )
                self.log_scaffold_length(scaffolds[0])
                if len(scaffolds) > 2:
                    logging.info("                ...  ...")
                if len(scaffolds) > 1:
                    self.log_scaffold_length(scaffolds[-1])

            if len(scaffolds) > 1 and len(ranked_scaffolds) > 1:
                logging.info(f"    {total:15,d}  bp total")

    def log_scaffold_length(self, scffld: Scaffold):
        logging.info(f"    {scffld.fragments_length:15,d}  {scffld.name}")

    def merge_unlocs(self, scaffolds: list[Scaffold]):
        this = Scaffold(scaffolds[0].name)
        merged = [this]
        for scffld in scaffolds:
            if not scffld.name.startswith(this.name):
                this = Scaffold(scffld.name)
                merged.append(this)
            this.append_scaffold(scffld)
        return merged
