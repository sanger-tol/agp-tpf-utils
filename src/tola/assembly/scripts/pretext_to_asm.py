import logging
import os
import re
import sys
from pathlib import Path
from shutil import which

import click

from tola.assembly.assembly import Assembly
from tola.assembly.assembly_set import AssemblySet
from tola.assembly.build_assembly import BuildAssembly
from tola.assembly.file_writers import (
    write_assemblies,
    write_assembly,
    write_assembly_stats,
    write_chr_csv_files,
    write_chr_report_csv,
    write_info_yaml,
    write_sum_chrs,
)
from tola.assembly.gfa_stats import GfaStatsError
from tola.assembly.indexed_assembly import IndexedAssembly
from tola.assembly.naming_utils import ChrNamerError, TaggingError
from tola.assembly.parser import format_from_file_extn, parse_agp, parse_tpf
from tola.assembly.sanger_files import (
    AssemblyYaml,
    AssemblyYamlError,
    AssemblyYamlLocationError,
    find_yaml,
)
from tola.fasta.index import FastaIndex
from tola.fasta.stream import FastaCollection

log = logging.getLogger(__name__)


def bd(txt):
    return click.style(txt, bold=True)


def it(txt):
    return click.style(txt, italic=True)


def ul(txt):
    return click.style(txt, underline=True)


@click.command(
    help=f"""
      Uses fragments in the assembly (AGP) produced by PretextView to find
      matching fragments in the assembly which was fed into Pretext and
      output an assembly made from the input assembly fragments.

      {ul("Named Chromsomes")}

        Upper case letters followed by zero or more digits are assumed to be
      chromosome names. {it("e.g.")} 'X', 'W', 'B1'

      {ul("Known Tags")}

        {bd("Contaminant")} tagged scaffolds are saved in a separate
      'Contaminants' file.

        When there are large numbers of contaminant scaffolds in the assembly,
      {bd("Target")} tags can insted be used to label the non-contaminant
      scaffolds and reduce the amount of labelling necessary in PretextView.
      Any un-tagged scaffolds will then be treated as if they were tagged
      with {it("Contaminant")}. (Any contaminants occurring before the first
      {it("Target")} tag in the PretextView AGP must still be individually
      tagged with{it("Contaminant")}.)

        {bd("FalseDuplicate")} for tagging duplicated regions in multi-haplotype
      Pretext maps which should be removed, not moved to another haplotype.

        {bd("Haplotig")} taggged scaffolds are saved in a separate 'Haplotigs'
      file. The haplotig scaffolds receive names 'H_1' to 'H_{it("n")}',
      sorted and numbered from longest to shortest.

        {bd("Primary")} in a multi-haplotpye Pretext map where only one of the
      haplotypes is being curated, is used to tag the first 'Painted' chromosome
      in the curated haplotype.

        {bd("Singleton")} is used to flag autosomes which were not found in
      any other haplotype.

        {bd("Unloc")} tagged scaffolds receive names '{it("CHR")}_unloc_1'
      to '{it("CHR")}_unloc_{it("n")}', added to the end of their
      chromosome and sorted and numbered from longest to
      shortest.

      {ul("Haplotypes")}

        Any other tags are assumed to be the name of a haplotype, and their
      assemblies are placed in separate files. Unplaced scaffolds for each
      haplotype are identified by their names beginning with the
      haplotype's name followed by an underscore. {it("i.e.")} 'Hap2_' for
      'Hap2'
      """,
)
@click.option(
    "--assembly",
    "-a",
    "assembly_file",
    type=click.Path(
        path_type=Path,
        exists=True,
        readable=True,
        resolve_path=True,
        dir_okay=False,
    ),
    required=True,
    help="""Assembly before curation, usually a FASTA file.
      FASTA files will be indexed, creating a '.fai' and a '.agp' file
      alongside the assembly if they are missing or are older than the
      FASTA.  Supports gzip compresed '.gz' files.""",
)
@click.option(
    "--pretext",
    "-p",
    "pretext_file",
    type=click.Path(
        path_type=Path,
        exists=True,
        readable=True,
        dir_okay=False,
    ),
    required=True,
    help="Assembly file from Pretext, which is usually an AGP.",
)
@click.option(
    "--output",
    "-o",
    "output_file",
    type=click.Path(
        path_type=Path,
        dir_okay=False,
    ),
    help=f"""Output file template, typically: '<ToLID>.<VERSION>.fa'

      {it("e.g.")} --output mVulVul1.2.fa

      for version 2 of the assembly of 'mVulVul1'. If <VERSION> is not
      specified, it defaults to '1'.

      The output file type is determined from its extension. When the outuput
      is FASTA ('.fa'), an AGP format file ('.fa.agp') is also written.  FASTA
      output files can be gzip compressed ('.fa.gz').

      The names of output files created are printed to STDERR.

      If not given, prints to STDOUT in 'STR' format.
      """,
)
@click.option(
    "--autosome-prefix",
    "-c",
    default="SUPER_",
    show_default=True,
    help="Prefix for naming autosomal chromosomes.",
)
@click.option(
    "--clobber/--no-clobber",
    "-f",
    default=True,
    show_default=True,
    help="Overwrite any existing output files.",
)
@click.option(
    "--log-level",
    "-l",
    type=click.Choice(
        ["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        case_sensitive=False,
    ),
    default="INFO",
    show_default=True,
    help="Diagnostic messages to show.",
)
@click.option(
    "--write-log/--no-write-log",
    "-w/-W",
    default=True,
    show_default=True,
    help="Write messages into a '.log' file alongside the output file.",
)
@click.option(
    "--info-yaml",
    "info_yaml_file",
    type=click.Path(
        path_type=Path,
        exists=True,
        readable=True,
        dir_okay=False,
    ),
    help="""

      Location of the YAML information file indicating the location of
      haplotigs and mitochondrion and chloroplast genome assemblies.  These
      are added to the output assemblies.  Files should be in FASTA format,
      optionally gzip compressed. Keys used are:

      \b
        haplotigs  - haplotigs not present in the
                     curated map
        mito       - mitochrondrial genome
        plastid    - chloroplast genome
      \b
    """,
)
@click.option(
    "--auto-find-yaml/--no-auto-find-yaml",
    default=False,
    show_default=True,
    envvar="SANGER_AUTO_FIND_YAML",
    help="""
      Find the YAML draft assembly information file using the directory
      structure expected at the Wellcome Sanger Institute, which is
      any ".yaml" file inside any "assembly/drafts" subdirectory off any
      parent directory of the "--assembly" input file.

      This option will default FASTA output files to be gzip compressed even
      when the '.gz' file extension is not used.
    """,
)
@click.option(
    "--keep-map-order",
    "-k",
    flag_value=True,
    default=False,
    show_default=True,
    help="Output a single assembly in the order of input Pretext map AGP.",
)
@click.option(
    "--default-assembly-name",
    "default_asm_name",
    type=str,
    help="""
        Name for output files which are not split per haplotype. Defaults
        to 'map-order' if the '--keep-map-order' option is given.
        """,
)
@click.option(
    "--max-contig-length",
    "max_contig_length",
    type=int,
    help=f"""
        Maximum length for a single contig.  If a contig exceeds this size, it
        will be broken at gaps into a number of pieces such that each is less
        than this number.  The pieces will have the suffixes '_1', '_2'
        {it("etc...")} added to their names.  [default: 2 Gbp]
        """,
    default=2_000_000_000,
)
@click.option(
    "--no-max-contig-length",
    flag_value=True,
    default=False,
    help="Turns off splitting on --max-contig-length",
)
@click.option(
    "--min-contig-length",
    "min_contig_length",
    type=int,
    help="""
        Minimum length for a single contig.  Contigs shorter that this are discarded.
        """,
    default=1000,
    show_default=True,
)
def cli(
    assembly_file: Path,
    pretext_file: Path,
    output_file: Path | None,
    autosome_prefix: str,
    clobber: bool,
    log_level: str,
    write_log: bool,
    auto_find_yaml: bool,
    info_yaml_file: Path,
    keep_map_order: bool,
    default_asm_name: str,
    max_contig_length: int,
    no_max_contig_length: bool,
    min_contig_length: int,
):
    logfile = setup_logging(log_level, output_file, write_log, clobber)

    # Locate and parse the draft assembly YAML file
    if info_yaml_file:
        info_yaml_file = info_yaml_file.absolute()
    elif auto_find_yaml:
        try:
            info_yaml_file = find_yaml(assembly_file.parent)
        except AssemblyYamlLocationError as ayle:
            for msg in ayle.args:
                log.error(msg)
            sys.exit("Error finding draft assembly YAML file")
    draft_yaml = AssemblyYaml(info_yaml_file) if info_yaml_file else None
    if draft_yaml:
        check_for_executable("gfastats")

    if keep_map_order and not default_asm_name:
        default_asm_name = "map-order"

    asm, fai = parse_assembly_file(assembly_file, "FASTA")
    input_asm = IndexedAssembly.new_from_assembly(asm)
    prtxt_asm, _ = parse_assembly_file(pretext_file, "AGP")

    # Trap "-a" and "-p" arguments being switched
    if not prtxt_asm.bp_per_texel:
        exit(
            f"No bp_per_texel value in the PretextView AGP file '{pretext_file}'\n"
            "(Are the -a/--assembly and -p/--pretext arguments the right way around?)"
        )

    build_asm = BuildAssembly(
        "stdout",
        autosome_prefix=autosome_prefix,
        min_contig_length=min_contig_length,
        max_contig_length=None if no_max_contig_length else max_contig_length,
        assembly_yaml=draft_yaml,
    )
    build_asm.remap_to_input_assembly(prtxt_asm, input_asm)

    try:
        out_assemblies = (
            build_asm.assembly_with_scaffolds_in_map_order()
            if keep_map_order
            else build_asm.assemblies_with_scaffolds_fused()
        )
    except ChrNamerError as cne:
        for msg in cne.args:
            log.info(msg)
        page_messages(cne.args)
        sys.exit("Error naming chromosomes")
    except TaggingError as te:
        for msg in te.args:
            log.warning(msg)
        sys.exit("Error in Pretext tags")

    # Build colletion of FASTA indexes for writing assembly
    fai_coll = FastaCollection(fai) if fai else None
    if fai_coll and draft_yaml:
        draft_yaml.add_indexes_to_collection(fai_coll)

    stats = build_asm.assembly_stats
    if output_file:
        out_fmt, out_dir, out_root, asm_version, suffix, gz_flag = parse_output_file(
            output_file, gz=auto_find_yaml
        )
        out_template = out_dir / f"{out_root}.{asm_version}"
        write_info_yaml(out_template, stats, out_assemblies, clobber)

        # Rename assemblies for output files
        out_assemblies = name_assemblies(
            out_assemblies, out_root, asm_version, default_asm_name
        )

        curated_asm_files = write_assemblies(
            fai_coll, out_fmt, out_dir, suffix, out_assemblies, clobber
        )
        write_chr_csv_files(out_dir, stats, out_assemblies, clobber, gz_flag)
        write_chr_report_csv(out_template, stats, out_assemblies, clobber)
        if draft_yaml and fai_coll:
            try:
                write_assembly_stats(draft_yaml, curated_asm_files, clobber)
                write_sum_chrs(out_dir, stats, out_assemblies, clobber)
            except GfaStatsError as ge:
                for msg in ge.args:
                    log.warning(msg)
                sys.exit("Error running gfastats")
            except AssemblyYamlError as ye:
                for msg in ye.args:
                    log.warning(msg)
                sys.exit(f"Error using draft assembly YAML: '{draft_yaml.file}'")

    else:
        for asm in out_assemblies.values():
            write_assembly(fai_coll, asm, None, None, clobber)

    for asm_key, out_asm in out_assemblies.items():
        stats.log_assembly_chromosomes(asm_key, out_asm)
    log.info("")
    stats.log_curation_stats()
    stats.log_sanity_checks(out_assemblies)
    if logfile:
        click.echo(f"  Log saved: '{logfile}'", err=True)


def setup_logging(log_level, output_file, write_log, clobber):
    conf = {
        "level": getattr(logging, log_level),
        # Leave messages unchanged:
        "format": "%(message)s",
        # Change config if called a second time (e.g. during testing):
        "force": True,
    }
    logfile = None
    if output_file and write_log:
        logfile = output_file.with_suffix(".log")
        conf["filename"] = logfile
        conf["filemode"] = "w" if clobber else "x"

    try:
        logging.basicConfig(**conf)  # ty: ignore[no-matching-overload]
    except FileExistsError:
        click.echo(f"ERROR: log file '{logfile}' already exists", err=True)
        sys.exit(1)

    if logfile:
        # Also print warnings to STDERR if logging to a file
        err_hdlr = logging.StreamHandler(sys.stderr)
        err_hdlr.setLevel(logging.WARNING)
        err_hdlr.setFormatter(logging.Formatter("%(levelname)s: %(message)s"))
        logging.getLogger().addHandler(err_hdlr)

    return logfile


def check_for_executable(cmd: str):
    if not which(cmd):
        log.error(f"No such executable {str!r} in PATH")
        sys.exit(1)


def name_assemblies(
    asm_dict: AssemblySet,
    root: str,
    version: str,
    default_asm_name: str,
) -> AssemblySet:
    """Rename assemblies for their output files"""

    ret_asm = AssemblySet()

    # A combined Pretext map of two or more haplotypes where only one of them
    # has been curated. One of the painted chromosomes in the curated
    # haplotype has been tagged with 'Primary'
    if asm_dict.get("Primary"):
        # <ToLID>.1.primary.curated.fa    <- Sequence from "Hap1" tagged scaffolds
        # <ToLID>.1.primary.chromosome.list.csv
        # <ToLID>.1.haplotigs.curated.fa  <- Sequence from "Hap2" tagged scaffolds
        other_asm = []
        for asm_key, asm in asm_dict.items():
            if asm_key == "Primary":
                asm.name = f"{root}.{version}.primary"
            elif asm.curated:
                other_asm.append(asm)
                continue
            else:
                asm.name = f"{root}.{version}.{asm_key.lower()}s"  # ty: ignore[unresolved-attribute]
            ret_asm[asm_key] = asm
        if other_asm:
            # Join the other haplotypes in the other assemblies for an
            # 'haplotigs' alternate assembly file
            htigs = merge_assemblies(other_asm)
            htigs.curated = True
            new_key = "haplotigs"
            htigs.name = f"{root}.{version}.{new_key}"
            ret_asm[new_key] = htigs

    # A single haplotype Pretext map
    elif asm_dict.get(None):
        # <ToLID>.1.primary.curated.fa
        # <ToLID>.1.primary.chromosome.list.csv
        # <ToLID>.1.haplotigs.curated.fa  <- Sequence from "Haplotig"
        #                                    tagged scaffolds
        other_asm = []
        for asm_key, asm in asm_dict.items():
            if asm_key is None:
                asm.name = (
                    f"{root}.{version}.{default_asm_name}"
                    if default_asm_name
                    else f"{root}.{version}.primary"
                )
                ret_asm[None] = asm
            elif asm_key == "Haplotig":
                new_key = "haplotigs"
                asm.name = f"{root}.{version}.{new_key}"
                asm.curated = True
                ret_asm[new_key] = asm
            else:
                asm.name = f"{root}.{version}.{asm_key.lower()}s"
                ret_asm[asm_key] = asm

    # Two or more haplotypes in a combined map
    else:
        # <ToLID>.hap1.1.primary.curated.fa
        # <ToLID>.hap1.1.primary.chromosome.list.csv
        # <ToLID>.hap2.1.primary.curated.fa
        # <ToLID>.hap2.1.primary.chromosome.list.csv
        for asm_key, asm in asm_dict.items():
            if asm.curated:
                asm.name = f"{root}.{asm_key.lower()}.{version}.primary"  # ty: ignore[unresolved-attribute]
            else:
                asm.name = f"{root}.{version}.{asm_key.lower()}s"  # ty: ignore[unresolved-attribute]
            ret_asm[asm_key] = asm

    return ret_asm


def merge_assemblies(asm_list):
    new = Assembly("merge")
    for asm in asm_list:
        if not new.source_haplotype:
            new.source_haplotype = asm.source_haplotype
        for scffld in asm.scaffolds:
            new.add_scaffold(scffld)
    return new


def parse_output_file(file: Path, gz=False) -> tuple[str, Path, str, str, str, bool]:
    out_fmt = format_from_file_extn(file)
    if out_fmt is None:
        reason = (
            "Missing extension"
            if file.suffix == ""
            else f"Unknown extension '{file.suffix}'"
        )
        click.echo(f"{reason} on file name: '{file}'", err=True)
        sys.exit(1)

    if file.suffix.lower() == ".gz":
        gz = True
        file = file.parent / file.stem

    sfx = f".{out_fmt.lower()}"
    if sfx.startswith(file.suffix.lower()):
        out_root = file.stem
        sfx = file.suffix
    else:
        out_root = file.name

    # Only gzip FASTA, not any other output formats
    if out_fmt == "FASTA" and gz:
        sfx += ".gz"

    # Is there a version suffix?
    if m := re.search(r"\.(\d+)$", out_root):
        version = m.group(1)
        # Clip version suffix off file name root
        out_root = Path(out_root).stem
    else:
        version = "1"

    return out_fmt, file.parent, out_root, version, sfx, gz


def parse_assembly_file(
    path: Path, default_format: str | None = None
) -> tuple[Assembly, FastaIndex | None]:
    fmt = format_from_file_extn(path, default_format)
    if fmt == "AGP":
        return parse_agp(path.open(), path.stem), None
    elif fmt == "TPF":
        return parse_tpf(path.open(), path.stem), None
    elif fmt == "FASTA":
        fai = FastaIndex(path)
        fai.auto_load()
        return fai.assembly, fai
    else:
        msg = f"Unknown assembly file format '{fmt}'"
        raise ValueError(msg)


def page_messages(itr):
    if sys.stdout.isatty():
        os.environ.setdefault(
            "LESS",
            " ".join(
                (
                    "--no-init",
                    "--quit-if-one-screen",
                    "--ignore-case",
                    "--RAW-CONTROL-CHARS",
                )
            ),
        )
        click.echo_via_pager(itr, color=True)
    else:
        for msg in itr:
            click.echo(msg, color=False, err=True)


if __name__ == "__main__":
    cli()
