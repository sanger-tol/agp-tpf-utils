import pathlib
import sys
from typing import TextIO

import click

from tola.assembly.format import format_agp, format_tpf
from tola.assembly.parser import format_from_file_extn, parse_agp, parse_tpf
from tola.fasta.index import FastaIndex


@click.command(
    help="""
      Parse and reformat ToL AGP and TPF files. Parses the files provided on
      the comamnd line, or STDIN if none are given.
    """,
)
@click.argument(
    "input_files",
    nargs=-1,
    type=click.Path(
        dir_okay=False,
        exists=True,
        readable=True,
        path_type=pathlib.Path,
    ),
)
@click.option(
    "--input-format",
    "-i",
    type=click.Choice(
        ["AGP", "TPF", "FASTA"],
        case_sensitive=False,
    ),
    help="""
      Format of input. Automatically determined from each input file's
      extension, or defaults to 'AGP'.
    """,
)
@click.option(
    "--output-file",
    "-o",
    type=click.Path(
        path_type=pathlib.Path,
        exists=False,
    ),
    help="""
      Output file. Format is guessed from extension. If no output file is
      given, ouput is printed to STDOUT.
    """,
)
@click.option(
    "--format",
    "-f",
    "output_format",
    type=click.Choice(
        ["AGP", "TPF", "STR", "REPR"],
        case_sensitive=False,
    ),
    help="""
      Format of output. Automatically determined from output file extension,
      or defaults to 'AGP'. 'STR' is a human-readable format, and 'REPR' is
      the parsed assembly object's data structure.
    """,
)
@click.option(
    "--name",
    "-n",
    "assembly_name",
    help="Name of the assembly. Defaults to the file name or 'stdin'.",
)
@click.option(
    "--qc-overlaps/--no-qc-overlaps",
    default=False,
    help="Report to STDERR any fragments within each assembly which overlap",
)
def cli(
    input_files, input_format, output_file, output_format, assembly_name, qc_overlaps
):
    if output_file:
        if not output_format:
            output_format = format_from_file_extn(output_file)
        out_fh = output_file.open("w")
    else:
        out_fh = sys.stdout

    if not output_format:
        output_format = "AGP"

    if input_files:
        for pth in input_files:
            # Select the format of the input file
            in_fmt = input_format or format_from_file_extn(pth, default="AGP")

            # Select the name of the assembly
            asm_name = assembly_name if assembly_name else pth.stem

            try:
                process_file(
                    pth,
                    in_fmt,
                    asm_name,
                    out_fh,
                    output_format,
                    qc_overlaps,
                )
            except Exception as e:
                msg = f"Error processing file '{pth}'"
                raise ValueError(msg) from e
    else:
        # Process STDIN
        in_fmt = input_format if input_format else "AGP"
        if in_fmt == "FASTA":
            sys.exit("Cannot process FASTA on STDIN")
        asm_name = assembly_name if assembly_name else "stdin"
        process_file(
            sys.stdin,
            in_fmt,
            asm_name,
            out_fh,
            output_format,
            qc_overlaps,
        )


def process_file(in_file, in_fmt, asm_name, out_fh, out_fmt, qc_overlaps):
    in_fh: TextIO = in_file if in_file is sys.stdin else in_file.open("r")

    if in_fmt == "AGP":
        asm = parse_agp(in_fh, asm_name)
    elif in_fmt == "TPF":
        asm = parse_tpf(in_fh, asm_name)
    elif in_fmt == "FASTA":
        fai = FastaIndex(in_file)
        fai.auto_load()
        asm = fai.assembly
    else:
        msg = f"Unknown input format: '{in_fmt}'"
        raise ValueError(msg)

    if qc_overlaps and (pairs := asm.find_overlapping_fragments()):
        report_overlaps(asm_name, pairs)

    if out_fmt == "AGP":
        format_agp(asm, out_fh)
    elif out_fmt == "TPF":
        format_tpf(asm, out_fh)
    elif out_fmt == "STR":
        out_fh.write(str(asm))
    elif out_fmt == "REPR":
        out_fh.write(repr(asm))
    else:
        msg = f"Unknown output format: '{out_fmt}'"
        raise ValueError(msg)

    if in_fh is not sys.stdin:
        in_fh.close()


def report_overlaps(asm_name, pairs):
    click.echo(f"\nOverlaps detected in assembly '{asm_name}'", err=True)
    for pr in pairs:
        f1, s1 = pr[0]
        f2, s2 = pr[1]
        click.echo(f"\nOverlap:\n{s1.name} {f1}\n{s2.name} {f2}", err=True)


if __name__ == "__main__":
    cli()
