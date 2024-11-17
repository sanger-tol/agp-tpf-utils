import csv
import io
import logging

from tola.assembly.assembly import Assembly
from tola.assembly.scaffold import Scaffold


class AssemblyStats:
    def __init__(self, autosome_prefix: str = "SUPER_") -> None:
        self.autosome_prefix = autosome_prefix
        self.input_assembly = None
        self.cuts = 0
        self.breaks = 0
        self.joins = 0
        self.per_assembly_stats = {}
        self.assembly_scaffold_lengths = {}

    def make_stats(self, output_assemblies: dict[str | None, Assembly]) -> None:
        input_junction_sets = self.input_assembly.fragment_junctions_by_asm_prefix()
        input_set = set()
        for junc_set in input_junction_sets.values():
            input_set |= junc_set

        output_junction_sets = {}
        output_set = set()
        for name, asm in output_assemblies.items():
            junc_set = asm.fragment_junction_set()
            output_set |= junc_set
            output_junction_sets[name] = junc_set

        # Breaks are junctions in the input that are not in the output
        total_breaks = input_set - output_set
        self.breaks = len(total_breaks)

        # Joins are junctions in the output that were not in the input
        total_joins = output_set - input_set
        self.joins = len(total_joins)

        for name, junc_set in output_junction_sets.items():
            junc_key = name.lower() if name else None
            if input_asm_set := input_junction_sets.get(junc_key):
                self.per_assembly_stats[name or "Primary"] = {
                    # Breaks are junctions which were in the input, but are
                    # not in this assembly.  This set is then intersected
                    # with the total set of breaks to avoid counting
                    # junctions in scaffolds which have been moved between
                    # haplotypes.
                    "manual_breaks": len((input_asm_set - junc_set) & total_breaks),
                    # Joins are anything new in this assembly compared to the
                    # input which is also in the total set of joins.
                    "manual_joins": len((junc_set - input_asm_set) & total_joins),
                }

    def log_curation_stats(self):
        cut_plural = "cut in a contig" if self.cuts == 1 else "cuts in contigs"
        break_plural = "break at a gap" if self.breaks == 1 else "breaks at gaps"
        join_plural = "join" if self.joins == 1 else "joins"
        logging.info(
            f"Curation made {self.cuts} {cut_plural}, {self.breaks}"
            f" {break_plural} and {self.joins} {join_plural}"
        )

    def ranked_scaffolds(self, asm: Assembly):
        ranked_scaffolds = {}
        for scffld in asm.scaffolds:
            rank = scffld.rank
            ranked_scaffolds.setdefault(rank, []).append(scffld)
        return ranked_scaffolds

    def build_assembly_scaffold_lengths(self, asm: Assembly):
        ranked_scaffolds = self.ranked_scaffolds(asm)

        ranked_names_lengths = {}
        for rank, scaffolds in ranked_scaffolds.items():
            # If this isn't the Unplaced rank, merge in any Unlocs before
            # storing lengths
            if rank != 3:
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
        current_chr = None
        current_chr_list = None
        for scffld in asm.scaffolds:
            if scffld.rank in (1, 2):
                name = scffld.name
                if current_chr and name.startswith(current_chr):
                    current_chr_list.append(name)
                else:
                    current_chr = name
                    current_chr_list = [name]
                    chr_names.append(current_chr_list)
        return chr_names if chr_names else None

    def chromosomes_report_csv(self, hap_asm: dict[str, Assembly]):
        csv_str = io.StringIO()
        # Would prefer to use quoting=csv.QUOTE_STRINGS but it was introduced
        # only in Python 3.12
        csvr = csv.writer(csv_str, quoting=csv.QUOTE_NONNUMERIC)
        csvr.writerow(
            (
                "assembly",
                "seq_name",
                "chr_name",
                "localised",
                "pretext_scaffold",
                "length",
                "length_minus_gaps",
            )
        )
        head_pos = csv_str.tell()

        prefix = self.autosome_prefix
        current_root = None
        current_chr = None
        for hap, asm in hap_asm.items():
            if not hap:
                hap = "Primary"
            for scffld in asm.scaffolds:
                if scffld.rank in (1, 2):
                    name = scffld.name
                    if current_root and name.startswith(current_root):
                        chr_name = current_chr
                    else:
                        current_root = name
                        chr_name = current_chr = name.replace(prefix, "", 1)
                    csvr.writerow(
                        (
                            hap,
                            name,
                            chr_name,
                            "true" if name == current_root else "false",
                            scffld.original_name,
                            scffld.length,
                            scffld.fragments_length,
                        )
                    )

        return csv_str.getvalue() if csv_str.tell() > head_pos else None

    def log_assembly_chromosomes(self, asm_key: str | None, asm: Assembly):
        ranked_names_lengths = self.get_assembly_scaffold_lengths(asm_key, asm)

        logging.info(f"\n{asm.name}")
        logging.info(f"    {asm.fragments_length:15,d}  bp sequence (minus gaps)")
        is_main = False
        rank_label = {
            1: "Autosomes",
            2: "Named",
            3: "Unplaced",
        }
        for rank, name_length in ranked_names_lengths.items():
            # Does this assembly have autosomes or named chromosomes?
            if rank in (1, 2):
                is_main = True

            # Only show rank headings for main assemblies
            if is_main:
                logging.info(f"  {rank_label[rank]}:")
            logging.info(f"    n = {len(name_length)}")

            if rank == 2:
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
