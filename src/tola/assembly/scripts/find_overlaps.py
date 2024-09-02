import pathlib
import re
import sys

import click
from tola.assembly.fragment import Fragment
from tola.assembly.indexed_assembly import IndexedAssembly
from tola.assembly.parser import parse_agp, parse_tpf
from tola.assembly.scripts.asm_format import format_from_file_extn


@click.command(
    help="""Find all the overlapping Fragments to the '<name>:<start>-<end>'
      specification(s) provided on the command line in all the AGP or TPF
      formatted assembly files listed"""
)
@click.argument(
    "files_and_specs",
    nargs=-1,
    required=True,
)
def cli(files_and_specs):
    files = []
    specs = []
    for ele in files_and_specs:
        if m := re.search(r"([-\w]+):([\d_]+)-([\d_]+)", ele):
            name, start, end = m.groups()
            frag = Fragment(name, start, end, 1)
            specs.append(frag)
        else:
            path = pathlib.Path(ele)
            if not path.exists():
                error_exit(f"No such file: '{path}'", err=True)
            files.append(path)

    if not specs:
        error_exit(f"No '<name>:<start>-<end>' specifications given", err=True)

    for path in files:
        fmt = format_from_file_extn(path)
        if not fmt:
            error_exit(
                (
                    "Failed to determine file fromat"
                    + f" from file extension of '{path}'"
                ),
                err=True,
            )
        click.echo(f"\nFile: {path}")
        with path.open("r") as asm_fh:
            if fmt == "AGP":
                asm = parse_agp(asm_fh, path.stem)
            elif fmt == "TPF":
                asm = parse_tpf(asm_fh, path.stem)
            for scffld in asm.scaffolds:
                for i, frgmnt in scffld.idx_fragments():
                    for bait in specs:
                        if ovr := bait.overlap_length(frgmnt):
                            report_overlap(scffld, i, frgmnt, bait, ovr)


def report_overlap(scffld, i, frgmnt, bait, ovr):
    if len(scffld.rows) == 1:
        pos = "only row"
    elif i == 0:
        pos = "first row"
    elif scffld.rows[-1] is frgmnt:
        pos = "last row"
    else:
        pos = f"row {1 + i}"

    click.echo(
        f"  {pos} of {scffld.name}: {frgmnt}\n    overlaps {bait} by {ovr:,d} bp"
    )


def error_exit(msg, code=2):
    click.echo(msg, err=True)
    sys.exit(code)


if __name__ == "__main__":
    cli()
