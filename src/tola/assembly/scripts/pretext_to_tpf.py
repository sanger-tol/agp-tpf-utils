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
    help="""Uses fragments in the assembly (AGP) produced by PretextView to
      find matching fragments in the assembly (TPF) fed into Pretext and
      output an assembly made from the input assembly fragments.""",
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
    help="Assembly file from before curation, which is usually a TPF",
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
    help="Assembly file from Pretext, which is usually an AGP",
)
@click.option(
    "--output",
    "-o",
    "output_file",
    type=click.Path(
        path_type=pathlib.Path,
        dir_okay=False,
    ),
    help="Output file, usually a TPF",
)
@click.option(
    "--clobber/--no-clobber",
    default=False,
    show_default=True,
    help="Overwrite an existing output file",
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
    help="Diagnostic messages to show",
)
def cli(assembly_file, pretext_file, output_file, clobber, log_level):
    logging.basicConfig(
        level=getattr(logging, log_level),
        format="%(message)s",  # Leaves message unchanged
    )

    input_asm = IndexedAssembly.new_from_assembly(
        parse_assembly_file(assembly_file, "TPF")
    )
    prtxt_asm = parse_assembly_file(pretext_file, "AGP")
    out_name = output_file.stem if output_file else "stdout"
    build_asm = BuildAssembly(out_name, default_gap=Gap(200, "scaffold"))
    build_asm.remap_to_input_assembly(prtxt_asm, input_asm)

    out_assemblies = build_asm.assemblies_with_scaffolds_fused()
    for out_asm in out_assemblies:
        write_assembly(out_asm, output_file, clobber)
    logging.info("")
    build_asm.assembly_stats.log_stats()


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
        logging.warn(f"Wrote assembly to: {output_file}")


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
