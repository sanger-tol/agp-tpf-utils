import click
import pathlib
import sys

from tola.assembly.assembly import Assembly
from tola.assembly.gap import Gap
from tola.assembly.indexed_assembly import IndexedAssembly
from tola.assembly.format import format_agp, format_tpf
from tola.assembly.parser import parse_agp, parse_tpf
from tola.assembly.scaffold import Scaffold
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
    help="Input file from assembly, which is usually a TPF",
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
    help="Input file from Pretext, which is usually an AGP",
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
def cli(assembly_file, pretext_file, output_file, clobber):
    input_asm = IndexedAssembly.new_from_assembly(
        parse_assembly_file(assembly_file, "TPF")
    )
    prtxt_asm = parse_assembly_file(pretext_file, "AGP")

    if output_file:
        out_fmt = format_from_file_extn(output_file, "TPF")
        out_name = output_file.stem
        try:
            out_fh = output_file.open("w" if clobber else "x")
        except FileExistsError:
            click.echo(f"Error: output file '{output_file}' already exists", err=True)
            sys.exit(1)
    else:
        out_fmt = "STR"
        out_name = "stdout"
        out_fh = sys.stdout

    default_gap = Gap(200, "scaffold")
    scffld_n = 0
    out_asm = Assembly(out_name)
    for prtxt_scffld in prtxt_asm.scaffolds:
        scffld_n += 1
        out_scffld = Scaffold(f"R{scffld_n}")
        out_asm.add_scaffold(out_scffld)
        for prtxt_frag in prtxt_scffld.fragments():
            ### Deal with negative strand
            if overlaps := input_asm.fetch_overlaps(prtxt_frag):
                if out_scffld.last_row_is_fragment:
                    out_scffld.add_row(default_gap)
                out_scffld.rows.extend(overlaps)
            else:
                click.echo(f"No overaps found for: {prtxt_frag}", err=True)

    if out_fmt == "TPF":
        format_tpf(out_asm, out_fh)
    elif out_fmt == "AGP":
        format_agp(out_asm, out_fh)
    elif out_fmt == "STR":
        out_fh.write(str(out_asm))


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
