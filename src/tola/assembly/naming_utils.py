import logging
import re

from tola.assembly.fragment import Fragment
from tola.assembly.scaffold import Scaffold
from tola.assembly.terminal_table import TerminalTable, bold, bold_red

log = logging.getLogger(__name__)


class TaggingError(Exception):
    """Error in Pretext tags"""


class ScaffoldNamer:
    """
    Labels Scaffolds with named chromosomes (sex chromosomes, B chromosomes),
    ranks (autosomes=1, named=2, unplaced=3), and haplotypes as the Pretext
    assembly is processed. Also saves the original pretext scaffold name
    under the `original_name` attribute.
    """

    def __init__(self, autosome_prefix="SUPER_"):
        self.autosome_prefix = autosome_prefix
        self.current_scaffold_name = None
        self.current_rank = None
        self.current_haplotype = None
        self.haplotig_n = 0
        self.haplotig_scaffolds = []
        self.primary_haplotype = None
        self.target_tags = False
        self.unloc_n = 0
        self.unloc_scaffolds = []

        # Halplotype names stored under their lower case names
        self.haplotype_lc_dict = {}

    OTHER_KNOWN_TAGS = {
        "Contaminant",
        "Cut",
        "FalseDuplicate",
        "Haplotig",
        "Singleton",
        "Unloc",
    }

    def make_scaffold_name(self, scaffold: Scaffold, fragment_tags=None):
        """
        Using the tags from Pretext in the Scaffold, work out what the
        haplotype is, if it has been named, and what its rank is.
        """
        scaffold_name = None
        haplotype = None
        is_painted = False  # Is a curated chromsome
        rank = None
        primary_tag = False

        if not fragment_tags:
            fragment_tags = scaffold.fragment_tags()

        for tag in fragment_tags:
            if tag == "Painted":
                is_painted = True
            elif tag == "Target":
                self.target_tags = True
            elif tag == "Primary":
                primary_tag = True
            elif re.fullmatch(r"([A-Z]\d*[A-Z]?|[IVX_]+|\d+[A-Z]+)", tag):
                # This tag looks like a chromosome name, e.g. "X1", "I_II", "2RL"
                if scaffold_name and tag != scaffold_name:
                    msg = (
                        f"Found more than one scaffold_name name: '{scaffold_name}'"
                        f" and '{tag}' in scaffold:\n\n{scaffold}"
                    )
                    raise TaggingError(msg)
                scaffold_name = tag
                rank = 2
            elif tag not in self.OTHER_KNOWN_TAGS:
                # Any tag that doesn't look like a chromosome name is assumed
                # to be a haplotype, and we only expect to find one within
                # each Pretext Scaffold
                if haplotype:
                    msg = (
                        f"Found both '{haplotype}' and '{tag}', when only one'"
                        f" is expected, in scaffold:\n\n{scaffold}"
                    )
                    raise TaggingError(msg)
                else:
                    haplotype = self.get_set_haplotype(tag)

        # If we don't have a haplotype from a tag, try to get it from the
        # first row of the scaffold
        if not haplotype:
            haplotype = self.haplotype_from_first_row_name(scaffold)

        if primary_tag and not self.primary_haplotype:
            if not haplotype:
                msg = (
                    f"Failed to determine haplotype for Primary"
                    f" from scaffold:\n\n{scaffold}"
                )
                raise TaggingError(msg)
            self.primary_haplotype = self.get_set_haplotype(haplotype)
            log.debug(f"Primary haplotype is '{self.primary_haplotype}'")

        if not scaffold_name:
            if is_painted:
                scaffold_name = scaffold.name
                if not rank:
                    rank = 1  # Rank for autosomes
            else:
                # Unpainted scaffolds keep the name they have in the input
                # assembly
                scaffold_name = scaffold.rows[0].name
                rank = 3

        if prim := self.primary_haplotype:
            self.current_haplotype = "Primary" if haplotype == prim else haplotype
        else:
            self.current_haplotype = haplotype

        self.current_scaffold_name = scaffold_name
        self.current_rank = rank
        self.unloc_n = 0
        self.unloc_scaffolds.append([])

    def label_scaffold(
        self,
        scaffold: Scaffold,
        fragment: Fragment,
        scaffold_tags: set[str],
        original_name: str,
    ):
        name = self.current_scaffold_name
        rank = self.current_rank
        if "Contaminant" in fragment.tags or (
            self.target_tags and "Target" not in scaffold_tags
        ):
            scaffold.tag = "Contaminant"
            rank = 3
        if "FalseDuplicate" in fragment.tags:
            scaffold.tag = "FalseDuplicate"
            rank = 3
        elif "Haplotig" in fragment.tags:
            name = self.haplotig_name()
            scaffold.tag = "Haplotig"
            rank = 3
            self.haplotig_scaffolds.append(scaffold)
        elif "Unloc" in fragment.tags:
            if "Painted" not in scaffold_tags:
                msg = f"Unloc in unpainted scaffold {original_name!r}: {fragment}"
                raise TaggingError(msg)
            scaffold.chr_name = name
            name = self.unloc_name()
            self.unloc_scaffolds[-1].append(scaffold)
        elif "Painted" in scaffold_tags:
            scaffold.localised = True
            scaffold.chr_name = name

        scaffold.name = name
        scaffold.haplotype = self.current_haplotype
        scaffold.rank = rank
        scaffold.original_name = original_name
        scaffold.original_tags = scaffold_tags

    def haplotype_from_first_row_name(self, scaffold):
        if m := re.search(r"^([^_]+)_.+_\d+$", scaffold.rows[0].name):
            return self.get_set_haplotype(m.group(1))
        else:
            return None

    def get_set_haplotype(self, haplotype):
        """
        The first occurance of a haplotype tag sets its case.
        i.e. "Hap1" will be used if it is seen before "HAP1".
        """
        return self.haplotype_lc_dict.setdefault(haplotype.lower(), haplotype)

    def haplotig_name(self) -> str:
        self.haplotig_n += 1
        return f"H_{self.haplotig_n}"

    def unloc_name(self) -> str:
        self.unloc_n += 1
        return f"{self.current_scaffold_name}_unloc_{self.unloc_n}"

    def rename_haplotigs_by_size(self):
        self.rename_by_size(self.haplotig_scaffolds)

    def rename_unlocs_by_size(self):
        for unloc_list in self.unloc_scaffolds:
            self.rename_by_size(unloc_list)

    def rename_by_size(self, scaffolds: list[Scaffold]):
        if not scaffolds:
            return
        names = [s.name for s in scaffolds]
        by_size = sorted(scaffolds, key=lambda s: s.length, reverse=True)
        for s, n in zip(by_size, names, strict=True):
            s.name = n


class ChrGroup:
    """
    Used to group together adjacent chromosomes from each haplotype so that
    they can be kept together when sorting by size in the first haplotype.
    """

    def __init__(self, haplotypes):
        self.data = data = {}
        for hap in haplotypes:
            data[hap] = {}

    def __repr__(self):
        data_summary = {hap: list(self.data[hap].keys()) for hap in self.data}
        return f"{self.__class__.__name__}(\n  data={data_summary})\n"

    def haplotype_dict(self, hap_name):
        return self.data.get(hap_name)

    def add_scaffold_to_haplotype(self, hap_name, scaffold):
        log.debug(f"Adding scaffold to '{hap_name}':\n{scaffold}")
        self.data.get(hap_name).setdefault(scaffold.original_name, []).append(scaffold)

    def original_tags_of_haplotype_scaffold(self, hap_name, scffld_name):
        return self.haplotype_dict(hap_name)[scffld_name][0].original_tags or ()

    def length_of_first_haplotype(self):
        first, *_ = self.data.values()

        length = 0
        for hap_scffld_grp in first.values():
            for scffld in hap_scffld_grp:
                length += scffld.fragments_length
        return length

    @staticmethod
    def multi_chr_list(chr_name, multi_count):
        """
        Adds the suffix "A", "B", "C" etc... to the supplied chromosome name
        for when there are multiple chromosomes in a group.
        """
        if multi_count == 1:
            return [chr_name]
        else:
            chr_list = []
            for ltr in range(ord("A"), ord("A") + multi_count):
                chr_list.append(chr_name + chr(ltr))
            return chr_list

    def max_hap_set_count(self):
        return max(len(hap_set) for hap_set in self.data.values())

    def name_chromosome(self, chr_prefix, chr_n):
        """
        Replace the original Pretext scaffold name with the supplied
        `chr_prefix` and `chr_n` for each haplotype within the group.

        e.g.
              "Scaffold_10"         > "SUPER_9A"
              "Scaffold_10_unloc_1" > "SUPER_9A_unloc_1"
              "Scaffold_10_unloc_2" > "SUPER_9A_unloc_2"
              "Scaffold_11"         > "SUPER_9B"
        """
        for hap_name, hap_set in self.data.items():
            hap_suffix = (
                "" if (hap_name is None or hap_name == "Primary") else f"_{hap_name}"
            )
            chr_names = self.multi_chr_list(str(chr_n), len(hap_set))
            for orig, scffld_list in hap_set.items():
                this_chr = chr_names.pop(0)
                for scffld in scffld_list:
                    scffld.name = scffld.name.replace(
                        orig,
                        chr_prefix + this_chr + hap_suffix,
                    )
                    scffld.chr_name = this_chr


class ChrNamerError(Exception):
    """
    An error in the expected pattern of scaffolds and haplotypes and tags when
    naming the autosomes.
    """


class ChrNamerUsageError(Exception):
    """
    Programming usage error in ChrNamer
    """


class ChrNamer:
    """
    Groups chromosomes across haplotypes and sorts and names them.
    """

    def __init__(self, chr_prefix="SUPER_"):
        self.chr_prefix = chr_prefix
        self.scaffolds = []
        self.haplotypes_seen = {}
        self.groups: list[ChrGroup] | None = None

    def add_scaffold(self, haplotype, scffld):
        # A dict is used to store the haplotypes seen since order is
        # significant and sets do not preserve order.
        # haplotype = str(hap)
        self.haplotypes_seen[haplotype] = True
        self.scaffolds.append((haplotype, scffld))

    def new_group(self):
        grp = ChrGroup(self.haplotypes_seen)
        self.groups.append(grp)
        return grp

    def add_chr_prefix(self, scffld: Scaffold, haplotype: str | None):
        prefix = self.chr_prefix
        if not scffld.name.startswith(prefix):
            scffld.name = prefix + scffld.name
        if haplotype and haplotype != "Primary":
            suffix = "_" + haplotype
            if not scffld.name.lower().endswith(suffix.lower()):
                scffld.name = scffld.name + suffix

    def add_haplotype_prefix(self, scffld: Scaffold, haplotype: str):
        if haplotype == "Primary":
            return
        prefix = haplotype + "_"
        if not scffld.name.lower().startswith(prefix.lower()):
            scffld.name = prefix + scffld.name

    def name_chromosomes(self):
        if not self.haplotypes_seen:
            # No autosomes to name
            return
        self.groups = []
        self.build_groups()
        self.groups.sort(key=lambda x: x.length_of_first_haplotype(), reverse=True)
        chr_prefix = self.chr_prefix
        for i, grp in enumerate(self.groups):
            grp.name_chromosome(chr_prefix, i + 1)

    def check_for_painted_scaffolds_missing_haplotype_tag(self):
        if len(self.haplotypes_seen) > 1 and None in self.haplotypes_seen:
            untagged = "".join(
                [f"  {scffld.name}\n" for hap, scffld in self.scaffolds if hap is None]
            )
            msg = f"Haplotype tag missing from Painted scaffolds:\n{untagged}"
            raise TaggingError(msg)

    def build_groups(self):
        self.check_for_painted_scaffolds_missing_haplotype_tag()

        first_haplotype, *other_haplotypes = self.haplotypes_seen
        group = self.new_group()
        last_haplotype = None
        last_orig = None
        for haplotype, scffld in self.scaffolds:
            orig = scffld.original_name
            if not orig:
                msg = f"Missing original_name value in Scaffold:\n{scffld}"
                raise ChrNamerUsageError(msg)

            # Do we already have a scaffold in this haplotype in the ChrGroup?
            if group.haplotype_dict(haplotype):
                if other_haplotypes:
                    if haplotype != last_haplotype:
                        # New haplotype which already has an entry in this
                        # group, so we must be in a new group.
                        group = self.new_group()
                    elif (
                        orig != last_orig  # i.e. not an Unloc
                        and "Singleton"
                        in group.original_tags_of_haplotype_scaffold(
                            haplotype, last_orig
                        )
                    ):
                        # Previous scaffold is tagged as a Singleton
                        group = self.new_group()
                elif orig != last_orig:
                    # There will be mulitple scaffolds in a row from with the
                    # same original_name when there are Unlocs.
                    # When there's only one haplotype, we make a new ChrGroup
                    # for each original_name, i.e. Pretext scaffold name.
                    group = self.new_group()

            # Append to the list under Haplotype > Pretext Scaffold in the
            # ChrGroup. i.e. Each ChrGroup will be structured like this:
            #
            #   ChrGroup(
            #       data={
            #           "Hap1": {
            #               "Scaffold_9": [
            #                   Scaffold(name="Scaffold_9"),
            #               ]
            #           },
            #           "Hap2": {
            #               "Scaffold_10": [
            #                   Scaffold(name="Scaffold_10"),
            #                   Scaffold(name="Scaffold_10_unloc_1)",
            #               ]
            #           },
            #       }
            #   )

            group.add_scaffold_to_haplotype(haplotype, scffld)
            last_haplotype = haplotype
            last_orig = orig

        table = self.check_groups()
        if table.errors:
            s = "" if len(table.errors) == 1 else "s"
            msg = f"Error{s} naming autosomes:\n"
            raise ChrNamerError(msg, *table.error_render())
        else:
            log.debug("\n" + table.render())

    def check_groups(self):
        tbl = TerminalTable()
        hdr = tbl.new_header()
        for hap in self.haplotypes_seen:
            hdr.new_cell().new_line(str(hap), bold)

        ignore = set(self.haplotypes_seen)
        ignore.add("Cut")
        ignore.add("Painted")
        if self.groups is None:
            msg = "Groups attribute not set"
            raise ChrNamerUsageError(msg)
        for grp in self.groups:
            row_count = grp.max_hap_set_count()
            for row_idx in range(row_count):
                row = tbl.new_row()
                for i, hap in enumerate(self.haplotypes_seen):
                    # Make a new cell, which may be empty
                    cell = row.new_cell()

                    # Are there any scaffolds for this haplotype in this ChrGroup?
                    if scaffolds := grp.data.get(hap):
                        # Is there a scaffold for this row of the ChrGroup?
                        if row_idx < len(scaffolds):
                            # Get the scaffold on this row
                            scffld_name = list(scaffolds)[row_idx]
                            scffld = scaffolds[scffld_name][0]
                            cell.new_line(scffld_name)

                            # # The first haplotype should only have one
                            # # scaffold in the group
                            # if i == 0 and row_idx > 0:
                            #     cell.new_line(f"<Consecutive {hap}>", bold_red)
                            #     tbl.mark_error()

                            # Show the fragments length of this scaffold
                            s_length = sum(
                                x.fragments_length for x in scaffolds[scffld_name]
                            )
                            cell.new_line(f"{s_length:,} bp")

                            # Show any tags on this scaffold
                            for tag in sorted(scffld.original_tags or ()):
                                if tag not in ignore:
                                    cell.new_line(tag)

                    # If this is the first line of the ChrGroup, is the first
                    # haplotype missing?
                    elif row_idx == 0 and i == 0:
                        cell.new_line("<empty>", bold_red)
                        tbl.mark_error()

        return tbl
