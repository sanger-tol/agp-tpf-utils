import re
import sys
from io import BytesIO
from pathlib import Path

from tola.assembly.assembly import Assembly
from tola.assembly.fragment import Fragment
from tola.assembly.gap import Gap
from tola.assembly.scaffold import Scaffold
from tola.fasta.info import FastaInfo


def index_fasta_file(
    file: Path, buffer_size: int = 250_000
) -> tuple[dict[str, FastaInfo], Assembly]:
    name = ""
    seq_length = 0
    region_start = None
    region_end = None
    seq_buffer = BytesIO()

    idx_dict = {}
    asm = Assembly(
        file.name,
        header=[f"Built from FASTA file '{file.absolute()}'"],
    )

    def store_info():
        process_seq_buffer()
        if region_end:
            seq_regions.append((region_start, region_end))

        if idx_dict.get(name):
            msg = f"More than one sequence named '{name}' in FASTA file '{file}'"
            raise ValueError(msg)
        idx_dict[name] = FastaInfo(
            seq_length,
            file_offset,
            residues_per_line,
            residues_per_line + line_end_bytes,
        )

        scffld = Scaffold(name)
        prev = (0, 0)
        for region in seq_regions:
            start, end = region
            if start != prev[1]:
                gap_length = start - prev[1]
                scffld.add_row(Gap(gap_length, "scaffold"))
            scffld.add_row(Fragment(name, start + 1, end, 1))
            prev = region
        if rem := seq_length - prev[1]:
            scffld.add_row(Gap(rem, "scaffold"))

        asm.add_scaffold(scffld)

    def process_seq_buffer():
        # Outer scope variables which we "rebind" in this function.
        # See https://peps.python.org/pep-3104/ for explanation.
        nonlocal seq_length, region_start, region_end

        # Take the value from the sequence buffer and empty it
        seq_bytes = seq_buffer.getvalue()
        seq_buffer.seek(0)
        seq_buffer.truncate(0)

        # Treat any non-ACGT character as an "N" (i.e. gap)
        for m in re.finditer(rb"[ACGTacgt]+", seq_bytes):
            start = seq_length + m.start()
            end = seq_length + m.end()
            if start == region_end:
                region_end = end
            else:
                if region_end:
                    seq_regions.append((region_start, region_end))
                region_start = start
                region_end = end

        seq_length += len(seq_bytes)

    # Reading the file in bytes mode is about 10% faster than text mode, which
    # has the overhead of decoding to UTF-8.
    with file.open("rb") as fh:
        for line in fh:
            # ord(">") == 62
            if line[0] == 62:
                # If this isn't the first sequence in the file, store the
                # accumulated data from the previous sequence.
                if name:
                    store_info()

                # Get new name by splitting on whitespace beyond the first
                # character and taking the first element of the array.
                # (This also allows space characters following the ">"
                # character of the header.)
                name = line[1:].split()[0].decode()
                if not name:
                    msg = f"Failed to parse sequence name from line:\n{line}"
                    raise ValueError(msg)

                # Reset variables for new sequence
                seq_length = 0
                residues_per_line = 0
                region_start = 0
                region_end = None
                seq_regions = []

                # The first residue of the sequence will be where the file
                # pointer now is.
                file_offset = fh.tell()

                # We assume each sequence entry will have the same line
                # endings.  Check for Windows "\r\n" line ending where the
                # second to last byte will be ord("\r") == 13
                line_end_bytes = 2 if line[-2] == 13 else 1
            else:
                if not residues_per_line:
                    residues_per_line = len(line) - line_end_bytes

                seq_buffer.write(line[:-line_end_bytes])
                if seq_buffer.tell() > buffer_size:
                    process_seq_buffer()

    # Store info for the last sequence in the file
    if name:
        store_info()

    if idx_dict:
        return idx_dict, asm
    else:
        msg = f"No data in FASTA file '{file.absolute()}'"
        raise ValueError(msg)


if __name__ == "__main__":
    for file in sys.argv[1:]:
        idx_dict, asm = index_fasta_file(Path(file))
        for name, info in idx_dict.items():
            sys.stdout.write(info.fai_row(name))
