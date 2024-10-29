#!/usr/bin/env python3

import io
import re
import sys
from pathlib import Path


class FastaInfo:
    __slots__ = (
        "name",
        "length",
        "file_offset",
        "residues_per_line",
        "max_line_length",
        "seq_regions",
    )

    def __init__(
        self,
        name,
        length,
        file_offset,
        residues_per_line,
        max_line_length,
        seq_regions=None,
    ):
        self.name = name
        self.length = length
        self.file_offset = file_offset
        self.residues_per_line = residues_per_line
        self.max_line_length = max_line_length
        self.seq_regions = seq_regions

    def fai_row(self):
        """Returns a row for a Fasta Index (.fai) file."""
        numbers = "\t".join(
            str(x)
            for x in (
                self.length,
                self.file_offset,
                self.residues_per_line,
                self.max_line_length,
            )
        )
        return f"{self.name}\t{numbers}\n"

    def regions(self):
        s = io.StringIO()
        for start, end in self.seq_regions:
            s.write(f"{end - start + 1:14,d}  {self.name}:{start}-{end}\n")

        return s.getvalue()


def index_fasta_bytes(file: Path, buffer_size: int = 10e6):
    name = None
    seq_length = None
    file_offset = None
    residues_per_line = None
    region_start = None
    region_end = None
    seq_regions = None
    line_end_bytes = None
    seq_buffer = io.BytesIO()

    info = []

    # Opening the file in bytes mode means that Windows ("\r\n") or UNIX
    # ("\n") line endings are preserved.  It is also about 10% faster than
    # decoding to UTF-8.
    with file.open("rb") as fh:
        for line in fh:
            # ord(">") == 62
            if line[0] == 62:
                # If this isn't the first sequence in the file, store the
                # accumulated data from the previous sequence.
                if name:
                    if region_end:
                        seq_regions.append((region_start + 1, region_end))
                    info.append(
                        FastaInfo(
                            name,
                            seq_length,
                            file_offset,
                            residues_per_line,
                            residues_per_line + line_end_bytes,
                            seq_regions,
                        )
                    )

                # Get new name by splitting on whitespace beyond the first
                # character and taking the first element of the array. This
                # also allows space characters following the ">" character of
                # the header.
                name = line[1:].split()[0].decode("utf8")
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
                residues = len(line) - line_end_bytes
                if residues > residues_per_line:
                    residues_per_line = residues

                # Treat any non-ACGT character as an "N" (i.e. gap)
                for m in re.finditer(rb"[ACGTacgt]+", line[:-line_end_bytes]):
                    start = seq_length + m.start()
                    end = seq_length + m.end()
                    if start == region_end:
                        region_end = end
                    else:
                        if region_end:
                            seq_regions.append((region_start + 1, region_end))
                        region_start = start
                        region_end = end

                seq_length += residues
    if name:
        if region_end:
            seq_regions.append((region_start + 1, region_end))
        info.append(
            FastaInfo(
                name,
                seq_length,
                file_offset,
                residues_per_line,
                residues_per_line + line_end_bytes,
                seq_regions,
            )
        )

    return info


if __name__ == "__main__":
    for file in sys.argv[1:]:
        info = index_fasta_bytes(Path(file))
        for fst in info:
            sys.stdout.write("\n")
            sys.stdout.write(fst.fai_row())
            sys.stdout.write(fst.regions())
