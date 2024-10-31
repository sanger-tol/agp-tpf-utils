#!/usr/bin/env python3

import io
import logging
import re
import sys
from functools import cached_property
from pathlib import Path

from tola.assembly.assembly import Assembly
from tola.assembly.format import format_agp
from tola.assembly.fragment import Fragment
from tola.assembly.gap import Gap
from tola.assembly.parser import parse_agp
from tola.assembly.scaffold import Scaffold


class IndexUsageError(Exception):
    """Unexpected usage of FastaIndex"""


IUPAC_COMPLEMENT = bytes.maketrans(
    b"ACGTRYMKSWHBVDNacgtrymkswhbvdn",
    b"TGCAYRKMSWDVBHNtgcayrkmswdvbhn",
)


def reverse_complement(seq: bytes):
    return seq[::-1].translate(IUPAC_COMPLEMENT)


class FastaIndex:
    def __init__(self, fasta_file: Path | str):
        if not isinstance(fasta_file, Path):
            fasta_file = Path(fasta_file)
        if not fasta_file.exists():
            missing = str(fasta_file)
            raise FileNotFoundError(missing)
        self.fasta_file = fasta_file
        self.fai_file = Path(str(fasta_file) + ".fai")
        self.agp_file = Path(str(fasta_file) + ".agp")
        self.index = None
        self.assembly = None

    def auto_load(self):
        if self.check_for_index_files():
            self.load_index()
            self.load_assembly()
        else:
            self.run_indexing()

    def check_for_index_files(self):
        """
        Check that the .agp and fai files exist and are newer than the FASTA
        sequence file.
        """
        fasta_mtime = self.fasta_file.stat().st_mtime
        for idx_file in self.fai_file, self.agp_file:
            if not idx_file.exists():
                return False
            if not idx_file.stat().st_mtime > fasta_mtime:
                logging.warning(
                    f"Index file '{idx_file}' is older than"
                    f" FASTA file '{self.fasta_file}'"
                )
                return False
        return True

    def load_index(self):
        if self.index:
            msg = "Index FAI already loaded"
            raise IndexUsageError(msg)

        idx_dict = {}
        with self.fai_file.open() as idx:
            for line in idx:
                name, length, file_offset, residues_per_line, max_line_length = (
                    line.split()
                )
                idx_dict[name] = FastaInfo(
                    length,
                    file_offset,
                    residues_per_line,
                    max_line_length,
                )
            self.index = idx_dict

    def write_index(self):
        idx_dict = self.index
        if not idx_dict:
            msg = "No index data to write to FAI file"
            raise IndexUsageError(msg)
        if self.fai_file.exists():
            logging.warning(f"Overwriting FAI index file '{self.fai_file}'")
        with self.fai_file.open("w") as idx_fh:
            for name, info in idx_dict.items():
                idx_fh.write(info.fai_row(name))

    def load_assembly(self):
        if self.assembly:
            msg = "Assembly AGP already loaded"
            raise IndexUsageError(msg)
        self.assembly = parse_agp(self.agp_file.open(), self.fasta_file.name)

    def write_assembly(self):
        asm = self.assembly
        if not asm:
            msg = "No assembly data to write to AGP file"
            raise IndexUsageError(msg)
        if self.agp_file.exists():
            logging.warning(f"Overwriting AGP assembly file '{self.agp_file}'")
        with self.agp_file.open("w") as agp_fh:
            format_agp(asm, agp_fh)

    def run_indexing(self):
        idx_dict, assembly = index_fasta_file(self.fasta_file)
        self.index = idx_dict
        self.assembly = assembly
        self.write_index()
        self.write_assembly()

    @cached_property
    def fasta_fileandle(self):
        return self.fasta_file.open("rb")

    def get_sequence(self, frag: Fragment):
        fh = self.fasta_fileandle
        info = self.index.get(frag.name)
        if not info:
            msg = f"No sequence in index named '{frag.name}'"
            raise ValueError(msg)
        rpl = info.residues_per_line
        mll = info.max_line_length
        line_end_bytes = mll - rpl
        seq = io.BytesIO()

        lines_to_seek = mll * ((frag.start - 1) // rpl)
        line_offset = (frag.start - 1) % rpl
        fh.seek(info.file_offset + lines_to_seek + line_offset)
        head = 0
        if line_offset:
            ### Wrong if frag.length < head
            head = rpl - line_offset
            seq.write(fh.read(head))
            fh.seek(line_end_bytes, 1)
        remainder = frag.length - head
        whole_lines = remainder // rpl
        for _ in range(whole_lines):
            seq.write(fh.read(rpl))
            fh.seek(line_end_bytes, 1)
        if tail := remainder % rpl:
            seq.write(fh.read(tail))
        yield seq


class FastaInfo:
    __slots__ = (
        "length",
        "file_offset",
        "residues_per_line",
        "max_line_length",
    )

    def __init__(
        self,
        length,
        file_offset,
        residues_per_line,
        max_line_length,
    ):
        self.length = int(length)
        self.file_offset = int(file_offset)
        self.residues_per_line = int(residues_per_line)
        self.max_line_length = int(max_line_length)

    def __eq__(self, othr):
        for attr in self.__slots__:
            if getattr(self, attr) != getattr(othr, attr):
                return False
        return True

    def __repr__(self):
        return (
            "FastaInfo("
            + (", ".join(f"{attr}={getattr(self, attr)!r}" for attr in self.__slots__))
            + ")"
        )

    def fai_row(self, name):
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
        return f"{name}\t{numbers}\n"


def index_fasta_file(file: Path, buffer_size: int = 250_000):
    name = None
    seq_length = None
    file_offset = None
    residues_per_line = None
    region_start = None
    region_end = None
    seq_regions = None
    line_end_bytes = None
    seq_buffer = io.BytesIO()

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
