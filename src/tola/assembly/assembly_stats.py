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
        self.assembly_scaffold_lengths = {}

    def make_stats(self, output_assemblies: dict[str | None, Assembly]) -> None:
        input_set = self.input_assembly.fragment_junction_set()
        output_set = set()
        for asm in output_assemblies.values():
            output_set |= asm.fragment_junction_set()

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

    def ranked_scaffolds(self, asm: Assembly):
        autosome_prefix = self.autosome_prefix
        ranked_scaffolds = {}
        for scffld in asm.scaffolds:
            rank = scffld.rank(autosome_prefix)
            ranked_scaffolds.setdefault(rank, []).append(scffld)
        return ranked_scaffolds

    def build_assembly_scaffold_lengths(self, asm: Assembly):
        ranked_scaffolds = self.ranked_scaffolds(asm)

        ranked_names_lengths = {}
        for rank, scaffolds in ranked_scaffolds.items():
            rank_i, rank_name = rank
            # If this isn't the Unplaced rank, merge in any Unlocs before
            # storing lengths
            if rank_i != 2:
                scaffolds = self.merge_unlocs(scaffolds)
            name_length = {}
            for scffld in scaffolds:
                name_length[scffld.name] = scffld.fragments_length
            ranked_names_lengths[rank] = name_length

        return ranked_names_lengths

    def get_assembly_scaffold_lengths(self, asm_key: str | None, asm: Assembly):
        return self.assembly_scaffold_lengths.setdefault(
            asm_key, self.build_assembly_scaffold_lengths(asm)
        )

    def chromosome_names(self, asm: Assembly):
        chr_names = []
        for rank, scaffolds in self.ranked_scaffolds(asm).items():
            if rank[0] in (0, 1):
                current_chr = None
                current_chr_list = None
                for scffld in scaffolds:
                    name = scffld.name
                    if current_chr and name.startswith(current_chr):
                        current_chr_list.append(name)
                    else:
                        current_chr = name
                        current_chr_list = [name]
                        chr_names.append(current_chr_list)
        return chr_names if chr_names else None

    def log_assembly_chromosomes(self, asm_key: str | None, asm: Assembly):
        ranked_names_lengths = self.get_assembly_scaffold_lengths(asm_key, asm)

        logging.info(f"\n{asm.name}")
        logging.info(f"    {asm.fragments_length:15,d}  bp sequence (minus gaps)")
        is_main = False
        for rank, name_length in ranked_names_lengths.items():
            # rank_i is used in logic below so that code doesn't have to be
            # changed if we rename a rank
            rank_i, rank_name = rank

            # Does this assembly have autosomes or named chromosomes?
            if rank_i in (0, 1):
                is_main = True

            # Only show rank headings for main assemblies
            if is_main:
                logging.info(f"  {rank_name}:")
            logging.info(f"    n = {len(name_length)}")

            if rank_i == 1:
                # Show all the named scaffolds
                for name, length in name_length.items():
                    self.log_scaffold_length(name, length)
            else:
                # Print a summary of largest ... smallest scaffolds
                scaffolds = sorted(
                    name_length.items(), key=lambda nl: nl[1], reverse=True
                )
                # Show longest scaffold
                self.log_scaffold_length(*scaffolds[0])

                # Show the middle scaffold if there are only three
                if len(scaffolds) == 3:
                    self.log_scaffold_length(*scaffolds[1])
                # Omit the "..." line if there are only one or two scaffolds
                elif len(scaffolds) > 2:
                    logging.info("                ...  ...")

                # Show shortest scaffold
                if len(scaffolds) > 1:
                    self.log_scaffold_length(*scaffolds[-1])

            # Only show total if there's more than one item, or we would show
            # the same number twice.
            if len(name_length) > 1:
                total = sum(name_length.values())
                logging.info(f"    {total:15,d}  bp total")

    def log_scaffold_length(self, name, length):
        logging.info(f"    {length:15,d}  {name}")

    def merge_unlocs(self, scaffolds: list[Scaffold]):
        this = Scaffold(scaffolds[0].name)
        merged = [this]
        for scffld in scaffolds:
            if not scffld.name.startswith(this.name):
                this = Scaffold(scffld.name)
                merged.append(this)
            this.append_scaffold(scffld)
        return merged
