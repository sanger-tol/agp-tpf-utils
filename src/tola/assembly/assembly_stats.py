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

    def chromosome_name_csv(self, asm: Assembly):
        prefix = self.autosome_prefix
        current_root = None
        current_chr = None

        csv_str = io.StringIO()
        for scffld in asm.scaffolds:
            if scffld.rank in (1, 2):
                name = scffld.name
                if current_root and name.startswith(current_root):
                    chr_name = current_chr
                else:
                    current_root = name
                    chr_name = current_chr = name.replace(prefix, "", 1)
                csv_str.write(
                    ",".join(
                        (
                            name,
                            chr_name,
                            "yes" if name == current_root else "no",
                        )
                    )
                )
                csv_str.write("\n")

        return csv_str.getvalue() if csv_str.tell() else None

    def chromosomes_report_csv(self, hap_asm: dict[str | None, Assembly]):
        csv_str = io.StringIO()
        # Would prefer to use quoting=csv.QUOTE_STRINGS but it was introduced
        # only in Python 3.12
        csvr = csv.writer(csv_str, quoting=csv.QUOTE_NONNUMERIC)
        csvr.writerow(
            (
                "assembly",
                "seq_name",
                "chromosome",
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

    def log_sanity_checks(self, hap_asm: dict[str | None, Assembly]) -> None:
        for check in (
            self.check_consistent_autosome_count,
            self.check_for_large_haplotigs,
        ):
            if msg_list := check(hap_asm):
                for msg in msg_list:
                    logging.warning(msg)

    def check_consistent_autosome_count(
        self, hap_asm: dict[str | None, Assembly]
    ) -> str | None:
        chr_counts = {}
        for hap, asm in hap_asm.items():
            if autosomes := [x for x in asm.scaffolds if x.rank == 1]:
                chr_counts[hap if hap else "Primary"] = len(autosomes)
        if len(chr_counts) > 1:
            distinct_counts = set(chr_counts.values())
            if len(distinct_counts) > 1:
                return [
                    "Mismatch in autosome count between "
                    + (" and ".join([f"{x} = {n}" for x, n in chr_counts.items()]))
                ]
        return None

    def check_for_large_haplotigs(
        self, hap_asm: dict[str | None, Assembly]
    ) -> str | None:
        htigs = hap_asm.get("Haplotig")
        if not htigs:
            return None

        shortest = None
        for hap, asm in hap_asm.items():
            if hap == "Haplotig":
                continue
            ranked_names_lengths = self.get_assembly_scaffold_lengths(hap, asm)
            for rank in (1,2):
                if names_lengths := ranked_names_lengths.get(rank):
                    for frags_len in names_lengths.values():
                        if shortest and frags_len > shortest:
                            continue
                        shortest = frags_len

        if not shortest:
            return None

        msg_list = []
        ht_names_lengths = self.get_assembly_scaffold_lengths("Haplotig", htigs)[3]
        for ht in htigs.scaffolds:
            ht_len = ht_names_lengths[ht.name]
            if ht_len > shortest:
                msg_list.append(
                    f"Haplotig {ht.name} ({ht.original_name}) is {ht_len:,d} bp"
                    f" which is longer than the shortest chromosome ({shortest:,d} bp)"
                )
        return msg_list if msg_list else None

    def merge_unlocs(self, scaffolds: list[Scaffold]):
        this = Scaffold(scaffolds[0].name)
        merged = [this]
        for scffld in scaffolds:
            if not scffld.name.startswith(this.name):
                this = Scaffold(scffld.name)
                merged.append(this)
            this.append_scaffold(scffld)
        return merged
