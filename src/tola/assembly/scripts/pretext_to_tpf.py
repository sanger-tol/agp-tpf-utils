import click
import logging
import pathlib
import sys

from tola.assembly.assembly_stats import AssemblyStats
from tola.assembly.build_assembly import BuildAssembly
from tola.assembly.format import format_agp, format_tpf
from tola.assembly.gap import Gap
from tola.assembly.indexed_assembly import IndexedAssembly
from tola.assembly.parser import parse_agp, parse_tpf
from tola.assembly.scripts.asm_format import format_from_file_extn


@click.command(
    help=f"""
      Uses fragments in the assembly (AGP) produced by PretextView to find
      matching fragments in the assembly (TPF) which was fed into Pretext and
      output an assembly made from the input assembly fragments.

      {click.style("Named Chromsomes", underline=True)}

        Upper case letters followed by zero or more digits are assumed to be
      chromosome names. {click.style("e.g.", italic=True)} 'X', 'W', 'B1'

      {click.style("Known Tags", underline=True)}

        {click.style("Contaminant", bold=True)} tagged scaffolds are saved in
      a separate 'Contaminants' file.

        {click.style("Haplotig", bold=True)} taggged scaffolds are saved in a
      separate 'Haplotigs' file. The haplotig scaffolds receive names 'H_1'
      to 'H_{click.style("n", italic=True)}', sorted and numbered from
      longest to shortest.

        {click.style("Unloc", bold=True)} tagged scaffolds receive names
      '{click.style("CHR", italic=True)}_unloc_1' to
      '{click.style("CHR", italic=True)}_unloc_{click.style("n", italic=True)}',
      added to the end of their chromosome and sorted and numbered from
      longest to shortest.

      {click.style("Haplotypes", underline=True)}

        Any other tags are assumed to be the name of a haplotype, and their
      assemblies are placed in separate files. Unplaced scaffolds for each
      haplotype are identified by their names beginning with the
      haplotype's name followed by an underscore.
      {click.style("i.e.", italic=True)} 'HAP2_' for 'HAP2'
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
@click.option(
    "--autosome-prefix",
    "-c",
    default="RL_",
    show_default=True,
    help="Prefix for naming autosomal chromosomes.",
)
@click.option(
    "--clobber/--no-clobber",
    default=False,
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

    input_asm = IndexedAssembly.new_from_assembly(
        parse_assembly_file(assembly_file, "TPF")
    )
    prtxt_asm = parse_assembly_file(pretext_file, "AGP")
    out_name = output_file.stem if output_file else "stdout"
    build_asm = BuildAssembly(
        out_name,
        default_gap=Gap(200, "scaffold"),
        autosome_prefix=autosome_prefix,
    )
    build_asm.remap_to_input_assembly(prtxt_asm, input_asm)

    out_assemblies = build_asm.assemblies_with_scaffolds_fused()
    stats = build_asm.assembly_stats
    for asm_key, out_asm in out_assemblies.items():
        stats.log_assembly_chromosomes(asm_key, out_asm)
    logging.info("")
    stats.log_curation_stats()
    for out_asm in out_assemblies.values():
        write_assembly(out_asm, output_file, clobber)
    if output_file:
        write_chr_csv_files(output_file, stats, out_assemblies, clobber)
    if logfile:
        click.echo(f"Log saved to file '{logfile}'", err=True)


def setup_logging(log_level, output_file, write_log, clobber):
    conf = {
        "level": getattr(logging, log_level),
        "format": "%(message)s",  # Leaves message unchanged
    }
    logfile = None
    if output_file and write_log:
        logfile = output_file.with_suffix(".log")
        conf["filename"] = logfile
        conf["filemode"] = "w" if clobber else "x"
    logging.basicConfig(**conf)
    return logfile


def write_assembly(out_asm, output_file, clobber):
    if output_file:
        out_fmt = format_from_file_extn(output_file, "TPF")
        if out_asm.name != output_file.stem:
            output_file = output_file.with_stem(out_asm.name)
        try:
            out_fh = output_file.open("w" if clobber else "x")
        except FileExistsError:
            click.echo(f"Error: output file '{output_file}' already exists", err=True)
            sys.exit(1)
    else:
        out_fmt = "STR"
        out_fh = sys.stdout

    if out_fmt == "TPF":
        format_tpf(out_asm, out_fh)
    elif out_fmt == "AGP":
        format_agp(out_asm, out_fh)
    elif out_fmt == "STR":
        out_fh.write("\n")
        out_fh.write(str(out_asm))

    if out_fmt != "STR":
        op = "Overwrote" if clobber else "Created"
        click.echo(f"{op} file '{output_file}'", err=True)


def write_chr_csv_files(output_file, stats, out_assemblies, clobber):
    for asm_key, asm in out_assemblies.items():
        if chr_names := stats.chromosome_names(asm_key, asm):
            csv_file = output_file.parent / (
                f"chrs_{asm_key}.csv" if asm_key else "chrs.csv"
            )
            with csv_file.open("w" if clobber else "x") as csv_fh:
                for cn in chr_names:
                    csv_fh.write(cn + "\n")
            op = "Overwrote" if clobber else "Created"
            click.echo(f"{op} file '{csv_file}'")


def report_overlaps(pairs):
    for pr in pairs:
        f1, s1 = pr[0]
        f2, s2 = pr[1]
        logging.warn(f"\nDuplicated component:\n{s1.name} {f1}\n{s2.name} {f2}")


def parse_assembly_file(path, default_format=None):
    fmt = format_from_file_extn(path, default_format)
    path_fh = path.open("r")
    if fmt == "AGP":
        return parse_agp(path_fh, path.stem)
    elif fmt == "TPF":
        return parse_tpf(path_fh, path.stem)
    else:
        msg = f"Unknown assembly file format '{fmt}'"
        raise ValueError(msg)


if __name__ == "__main__":
    cli()
