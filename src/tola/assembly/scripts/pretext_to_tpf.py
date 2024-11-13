import logging
import pathlib
import sys

import click
import yaml

from tola.assembly.build_assembly import BuildAssembly
from tola.assembly.format import format_agp, format_tpf
from tola.assembly.gap import Gap
from tola.assembly.indexed_assembly import IndexedAssembly
from tola.assembly.parser import parse_agp, parse_tpf
from tola.assembly.scripts.asm_format import format_from_file_extn
from tola.fasta.index import FastaIndex
from tola.fasta.stream import FastaStream


def bd(txt):
    return click.style(txt, bold=True)


def it(txt):
    return click.style(txt, italic=True)


def ul(txt):
    return click.style(txt, underline=True)


@click.command(
    help=f"""
      Uses fragments in the assembly (AGP) produced by PretextView to find
      matching fragments in the assembly (TPF) which was fed into Pretext and
      output an assembly made from the input assembly fragments.

      {ul("Named Chromsomes")}

        Upper case letters followed by zero or more digits are assumed to be
      chromosome names. {it("e.g.")} 'X', 'W', 'B1'

      {ul("Known Tags")}

        {bd("Contaminant")} tagged scaffolds are saved in a separate
      'Contaminants' file.

        When there are large numbers of contaminant scaffolds in the
        assembly, {bd("Target")} tags can insted be used to label the
        non-contaminant scaffolds and reduce the amount of labelling
        necessary in PretextView. Any un-tagged scaffolds will then be
        treated as if they were tagged with {it("Contaminant")}.
        (Any contaminants occurring before the first {it("Target")} tag in
        the PretextView AGP must still be individually tagged with
        {it("Contaminant")}.)

        {bd("Haplotig")} taggged scaffolds are saved in a separate 'Haplotigs'
      file. The haplotig scaffolds receive names 'H_1' to 'H_{it("n")}',
      sorted and numbered from longest to shortest.

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
        path_type=pathlib.Path,
        exists=True,
        readable=True,
    ),
    required=True,
    help="Assembly file from before curation, which is usually a TPF.",
)
@click.option(
    "--pretext",
    "-p",
    "pretext_file",
    type=click.Path(
        path_type=pathlib.Path,
        exists=True,
        readable=True,
    ),
    required=True,
    help="Assembly file from Pretext, which is usually an AGP.",
)
@click.option(
    "--output",
    "-o",
    "output_file",
    type=click.Path(
        path_type=pathlib.Path,
        dir_okay=False,
    ),
    help="""Output file, usually a TPF.
      If not given, prints to STDOUT in 'STR' format.""",
)
# @click.option(
#     "--version",
#     "-v",
#     default="1",
#     show_default=True,
#     help="Version of the assembly, used when naming the output files.",
# )
@click.option(
    "--autosome-prefix",
    "-c",
    default="RL_",
    show_default=True,
    help="Prefix for naming autosomal chromosomes.",
)
@click.option(
    "--clobber/--no-clobber",
    "-f",
    default=True,
    show_default=True,
    help="Overwrite an existing output file.",
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
    default=False,
    show_default=True,
    help="Write messages into a '.log' file alongside the output file",
)
def cli(
    assembly_file,
    pretext_file,
    output_file,
    autosome_prefix,
    clobber,
    log_level,
    write_log,
):
    logfile = setup_logging(log_level, output_file, write_log, clobber)

    asm, fai = parse_assembly_file(assembly_file, "TPF")
    input_asm = IndexedAssembly.new_from_assembly(asm)
    prtxt_asm, _ = parse_assembly_file(pretext_file, "AGP")

    # Trap "-a" and "-p" arguments being switched
    if not prtxt_asm.bp_per_texel:
        exit(
            f"No bp_per_texel value in the PretextView AGP file '{pretext_file}'\n"
            "(Are the -a, --assembly and -p, --pretext arguments the right way around?)"
        )

    out_name = output_file.stem if output_file else "stdout"
    build_asm = BuildAssembly(
        out_name,
        default_gap=Gap(200, "scaffold"),
        autosome_prefix=autosome_prefix,
    )
    build_asm.remap_to_input_assembly(prtxt_asm, input_asm)

    out_assemblies = build_asm.assemblies_with_scaffolds_fused()
    for out_asm in out_assemblies.values():
        write_assembly(fai, out_asm, output_file, clobber)
    stats = build_asm.assembly_stats
    if output_file:
        write_chr_csv_files(output_file, stats, out_assemblies, clobber)
        write_info_yaml(output_file, stats, out_assemblies, clobber)
    for asm_key, out_asm in out_assemblies.items():
        stats.log_assembly_chromosomes(asm_key, out_asm)
    logging.info("")
    stats.log_curation_stats()
    if logfile:
        click.echo(f"Log saved to file '{logfile}'", err=True)


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
        logging.basicConfig(**conf)
    except FileExistsError:
        click.echo(f"Error: log file '{logfile}' already exists", err=True)
        sys.exit(1)
    return logfile


"""

<ToLID>.hap1.1.inter.csv                     # Supply better file with all chr

<ToLID>.hap1.1.primary.curated.fa
<ToLID>.hap1.1.all_haplotigs.curated.fa      # Do not make
<ToLID>.hap1.1.primary.chromosome.list.csv

<ToLID>.hap2.1.primary.curated.fa
<ToLID>.hap2.1.all_haplotigs.curated.fa
"""

def write_assembly(fai, out_asm, output_file, clobber):
    if output_file:
        out_fmt = format_from_file_extn(output_file, "TPF")
        mode = "b" if out_fmt == "FASTA" else ""
        if out_asm.name != output_file.stem:
            output_file = output_file.with_stem(out_asm.name)
        out_fh = get_output_filehandle(output_file, clobber, mode)
    else:
        out_fmt = "STR"
        out_fh = sys.stdout

    if out_fmt == "TPF":
        format_tpf(out_asm, out_fh)
    elif out_fmt == "AGP":
        format_agp(out_asm, out_fh)
    elif out_fmt == "FASTA":
        stream = FastaStream(out_fh, fai)
        stream.write_assembly(out_asm)

        # Save a .agp file alongside the .fa / .fasta
        output_agp = output_file.with_suffix(".agp")
        agp_fh = get_output_filehandle(output_agp, clobber)
        format_agp(out_asm, agp_fh)

    elif out_fmt == "STR":
        out_fh.write("\n")
        out_fh.write(str(out_asm))


def write_chr_csv_files(output_file, stats, out_assemblies, clobber):
    for asm_key, asm in out_assemblies.items():
        if chr_names := stats.chromosome_names(asm):
            csv_file = output_file.parent / (
                f"chrs_{asm_key}.csv" if asm_key else "chrs.csv"
            )
            with get_output_filehandle(csv_file, clobber) as csv_fh:
                for cn_list in chr_names:
                    csv_fh.write(",".join(cn_list) + "\n")


def write_info_yaml(output_file, stats, out_assemblies, clobber):
    asm_stats = stats.per_assembly_stats
    info = {"assemblies": asm_stats}
    if len(asm_stats) > 1:
        info["manual_breaks"] = stats.breaks
        info["manual_joins"] = stats.joins

    haplotig_count = 0
    if h_asm := out_assemblies.get("Haplotig"):
        haplotig_count = len(h_asm.scaffolds)
    info["manual_haplotig_removals"] = haplotig_count

    yaml_file = output_file.with_name(output_file.stem + "_info.yaml")
    with get_output_filehandle(yaml_file, clobber) as yaml_fh:
        yaml_fh.write(yaml.safe_dump(info, sort_keys=False))

    op = "Overwrote" if clobber else "Created"
    click.echo(f"{op} file '{yaml_file}'", err=True)


def get_output_filehandle(path, clobber, mode=""):
    op = "Overwrote" if path.exists() else "Created"
    try:
        out_fh = path.open("w" + mode if clobber else "x" + mode)
    except FileExistsError:
        click.echo(f"Error: output file '{path}' already exists", err=True)
        sys.exit(1)
    click.echo(f"{op} file '{path}'", err=True)
    return out_fh


def parse_assembly_file(path, default_format=None):
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


if __name__ == "__main__":
    cli()
