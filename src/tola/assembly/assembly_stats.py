import csv
import io
import logging
from typing import TypeAlias

from tola.assembly.assembly import Assembly, AssemblyDict

log = logging.getLogger(__name__)


RankedNameLengths: TypeAlias = dict[int, dict[str, int]]


class AssemblyStatsError(Exception):
    """Error from AssemblyStats"""


class AssemblyStats:
    def __init__(self, autosome_prefix: str = "SUPER_"):
        self.autosome_prefix = autosome_prefix
        self.input_assembly: Assembly | None = None
        self.cuts = 0
        self.breaks = None
        self.joins = None
        self.interventions_per_gbp = None
        self.percent_assembly_in_chromosomes = None
        self.per_assembly_stats = {}
        self.assembly_scaffold_lengths = {}

    def make_stats(self, output_assemblies: AssemblyDict):
        if not self.input_assembly:
            msg = "Missing input_assembly attribute"
            raise AssemblyStatsError(msg)
        self.__build_junction_stats(output_assemblies)
        self.__build_length_stats(output_assemblies)

    def __build_junction_stats(self, output_assemblies: AssemblyDict):

        # These stats are going to be wrong if re-curating assemblies with
        # Fragment names beginning with "SUPER_".
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

    def __build_length_stats(self, output_assemblies: AssemblyDict):
        input_asm_length = self.input_assembly.fragments_length

        # Calculate the number of breaks and joins made per Gbp of the input assembly
        self.interventions_per_gbp = round(
            (self.breaks + self.joins) / (input_asm_length / 1e9), 3
        )

        # Calculate the percentage of the assembly that was placed in chromosomes
        tier_1_and_2_length = 0
        for hap, asm in output_assemblies.items():
            ranked_scffld_lengths = self.get_assembly_scaffold_lengths(hap, asm)
            for rank in (1, 2):
                if scffld_lengths := ranked_scffld_lengths.get(rank):
                    tier_1_and_2_length += sum(scffld_lengths.values())
        self.percent_assembly_in_chromosomes = round(
            100 * (tier_1_and_2_length / input_asm_length), 1
        )

    def log_curation_stats(self):
        cut_plural = "cut in a contig" if self.cuts == 1 else "cuts in contigs"
        break_plural = "break at a gap" if self.breaks == 1 else "breaks at gaps"
        join_plural = "join" if self.joins == 1 else "joins"
        log.info(
            f"Curation made {self.cuts} {cut_plural}, {self.breaks}"
            f" {break_plural} and {self.joins} {join_plural}"
        )
        log.info(
            f"Assembly placed in chromosomes = {self.percent_assembly_in_chromosomes}%"
        )
        log.info(f"Interventions per Gbp = {self.interventions_per_gbp}")

    def ranked_scaffolds(self, asm: Assembly):
        ranked_scaffolds = {}
        for scffld in asm.scaffolds:
            rank = scffld.rank
            ranked_scaffolds.setdefault(rank, []).append(scffld)
        return ranked_scaffolds

    def build_assembly_scaffold_lengths(self, asm: Assembly) -> RankedNameLengths:
        ranked_scaffolds = self.ranked_scaffolds(asm)

        ranked_names_lengths = {}
        for rank, scaffolds in ranked_scaffolds.items():
            name_length = {}

            if rank == 3:
                for scffld in scaffolds:
                    name_length[scffld.name] = scffld.fragments_length
            else:
                current_name = None
                last_orig = None
                for scffld in scaffolds:
                    orig = scffld.original_name
                    if orig != last_orig:
                        current_name = scffld.name
                        last_orig = orig
                    name_length[current_name] = (
                        name_length.get(current_name, 0) + scffld.fragments_length
                    )

            ranked_names_lengths[rank] = name_length

        return ranked_names_lengths

    def get_assembly_scaffold_lengths(
        self,
        asm_key: str | None,
        asm: Assembly,
    ) -> RankedNameLengths:
        scaff_lengths = self.assembly_scaffold_lengths.get(asm_key)
        if not scaff_lengths:
            scaff_lengths = self.assembly_scaffold_lengths[asm_key] = (
                self.build_assembly_scaffold_lengths(asm)
            )
        return scaff_lengths

    def chromosome_name_csv(self, haplotype: str, asm: Assembly):
        csv_str = io.StringIO()
        for scffld in asm.scaffolds:
            if scffld.rank in (1, 2):
                csv_str.write(
                    ",".join(
                        (
                            scffld.name,
                            scffld.chr_name,
                            "yes" if scffld.localised else "no",
                        )
                    )
                )
                csv_str.write("\n")

        return csv_str.getvalue() if csv_str.tell() else None

    def chromosomes_report_csv(self, hap_asm: AssemblyDict):
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

        for hap, asm in hap_asm.items():
            for scffld in asm.scaffolds:
                if scffld.rank in (1, 2):
                    csvr.writerow(
                        (
                            scffld.haplotype or hap or "Primary",
                            scffld.name,
                            scffld.chr_name,
                            "true" if scffld.localised else "false",
                            scffld.original_name,
                            scffld.length,
                            scffld.fragments_length,
                        )
                    )

        return csv_str.getvalue() if csv_str.tell() > head_pos else None

    def log_assembly_chromosomes(self, asm_key: str | None, asm: Assembly):
        ranked_names_lengths = self.get_assembly_scaffold_lengths(asm_key, asm)

        log.info(f"\n{asm.name}")
        log.info(f"    {asm.fragments_length:15,d}  bp sequence (minus gaps)")
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
                log.info(f"  {rank_label[rank]}:")
            log.info(f"    n = {len(name_length)}")

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
                    log.info("                ...  ...")

                # Show shortest scaffold
                if len(scaffolds) > 1:
                    self.log_scaffold_length(*scaffolds[-1])

            # Only show total if there's more than one item, or we would show
            # the same number twice.
            if len(name_length) > 1:
                total = sum(name_length.values())
                log.info(f"    {total:15,d}  bp total")

    def log_scaffold_length(self, name, length):
        log.info(f"    {length:15,d}  {name}")

    def log_sanity_checks(self, hap_asm: dict[str | None, Assembly]) -> None:
        for check in (
            self.check_consistent_autosome_count,
            self.check_for_large_haplotigs,
        ):
            if msg_list := check(hap_asm):
                for msg in msg_list:
                    log.warning(msg)

    def check_consistent_autosome_count(
        self, hap_asm: AssemblyDict
    ) -> list[str] | None:
        chr_counts = {}
        for hap, asm in hap_asm.items():
            ranked_names_lengths = self.get_assembly_scaffold_lengths(hap, asm)
            if autosomes := ranked_names_lengths.get(1):
                chr_counts[hap or "Primary"] = len(autosomes)

        if len(chr_counts) > 1:
            distinct_counts = set(chr_counts.values())
            if len(distinct_counts) > 1:
                return [
                    "Mismatch in autosome count between "
                    + (" and ".join([f"{x} = {n}" for x, n in chr_counts.items()]))
                ]
        return None

    def check_for_large_haplotigs(self, hap_asm: AssemblyDict) -> list[str] | None:
        htigs = hap_asm.get("Haplotig")
        if not htigs:
            return None

        shortest = None
        for hap, asm in hap_asm.items():
            if hap == "Haplotig":
                continue
            ranked_names_lengths = self.get_assembly_scaffold_lengths(hap, asm)
            for rank in (1, 2):
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
        return msg_list or None
