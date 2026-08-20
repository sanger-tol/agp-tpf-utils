import gzip
import logging
from functools import cached_property
from io import BytesIO
from pathlib import Path

from tola.assembly.format import format_agp
from tola.assembly.fragment import Fragment
from tola.assembly.gap import Gap
from tola.assembly.parser import parse_agp
from tola.fasta.info import FastaInfo
from tola.fasta.parse import index_fasta_file
from tola.fasta.simple import FastaSeq, revcomp_bytes_io

log = logging.getLogger(__name__)


class IndexUsageError(Exception):
    """Unexpected usage of FastaIndex"""


class FastaIndex:
    def __init__(
        self,
        fasta_file: Path,
        buffer_size: int = 250_000,
        source: str | None = None,
    ):
        if not fasta_file.exists():
            missing = str(fasta_file)
            raise FileNotFoundError(missing)
        self.fasta_file: Path = fasta_file
        self.buffer_size: int = buffer_size
        self.fai_file: Path = fasta_file.with_name(f"{fasta_file.name}.fai")
        self.agp_file: Path = fasta_file.with_name(f"{fasta_file.name}.agp")
        self.source = source

    def auto_load(self):
        if self.check_for_index_files():
            self.load_index()
            self.load_assembly()
        else:
            self.run_indexing()
            self.write_index()
            self.write_assembly()

    def check_for_index_files(self) -> bool:
        """
        Check that the .agp and fai files exist and are newer than the FASTA
        sequence file.
        """
        fasta_mtime = self.fasta_file.stat().st_mtime
        for idx_file in self.fai_file, self.agp_file:
            if not idx_file.exists():
                return False
            if not idx_file.stat().st_mtime > fasta_mtime:
                log.warning(
                    f"Index file '{idx_file}' is older than"
                    f" FASTA file '{self.fasta_file}'"
                )
                return False
        return True

    def load_index(self) -> None:
        if hasattr(self, "index"):
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

    def write_index(self) -> None:
        if not hasattr(self, "index"):
            msg = "No index data to write to FAI file"
            raise IndexUsageError(msg)
        if self.fai_file.exists():
            log.warning(f"Overwriting FAI index file '{self.fai_file}'")
        with self.fai_file.open("w") as idx_fh:
            for name, info in self.index.items():
                idx_fh.write(info.fai_row(name))

    def load_assembly(self) -> None:
        if hasattr(self, "assembly"):
            msg = "Assembly AGP already loaded"
            raise IndexUsageError(msg)
        self.assembly = parse_agp(
            self.agp_file.open(), self.fasta_file.name, source=self.source
        )

    def write_assembly(self) -> None:
        if not hasattr(self, "assembly"):
            msg = "No assembly data to write to AGP file"
            raise IndexUsageError(msg)
        if self.agp_file.exists():
            log.warning(f"Overwriting AGP assembly file '{self.agp_file}'")
        with self.agp_file.open("w") as agp_fh:
            format_agp(self.assembly, agp_fh)

    def run_indexing(self) -> None:
        idx_dict, assembly = index_fasta_file(self.fasta_file, self.buffer_size)
        self.index = idx_dict
        self.assembly = assembly

    @cached_property
    def fasta_fileandle(self):
        ff = self.fasta_file
        return gzip.open(ff, "rb") if ".gz" in ff.suffix.lower() else ff.open("rb")

    def get_info(self, name):
        info = self.index.get(name)
        if not info:
            msg = f"No sequence in index named '{name}'"
            raise ValueError(msg)
        return info

    def get_gap_iter(self, gap: Gap, gap_character=b"N"):
        """
        Returns an iterator of `BytesIO` objects for gap characters for the Gap.
        Keeps memory usage below `buffer_size` for large gaps.
        """
        max_length = self.buffer_size
        length = gap.length
        chunk_count = 1 + (length // max_length)
        for i in range(chunk_count):
            chunk_start = i * max_length
            chunk_end = min(length, chunk_start + max_length)
            yield BytesIO(gap_character * (chunk_end - chunk_start))

    def get_sequence_iter(self, frag: Fragment):
        """
        Returns an iterator of `BytesIO` objects for sequence characters of
        the `Fragment`, keeping memory usage by the sequence data below
        `buffer_size`.
        """
        info = self.get_info(frag.name)

        if frag.strand == -1:
            return self.rev_chunks(info, frag.start, frag.end)
        else:
            return self.fwd_chunks(info, frag.start, frag.end)

    def fwd_chunks(self, info: FastaInfo, start, end):
        max_length = self.buffer_size
        chunk_count = 1 + ((end - start) // max_length)
        for i in range(chunk_count):
            offset = i * max_length
            chunk_start = start + offset
            chunk_end = min(end, chunk_start + max_length - 1)
            yield self.sequence_bytes(info, chunk_start, chunk_end)

    def rev_chunks(self, info: FastaInfo, start, end):
        max_length = self.buffer_size
        chunk_count = (end - start) // max_length

        # Loop backwards from last chunk to the first, yeilding the
        # reverse-complement of each chunk.
        for i in range(chunk_count, -1, -1):
            offset = i * max_length
            chunk_start = start + offset
            chunk_end = min(end, chunk_start + max_length - 1)
            yield revcomp_bytes_io(self.sequence_bytes(info, chunk_start, chunk_end))

    def all_fasta_seq(self):
        for name in self.index:
            yield self.get_fasta_seq(name)

    def get_fasta_seq(self, name) -> FastaSeq:
        info = self.get_info(name)
        seq_bytes = self.sequence_bytes(info, 1, info.length).getvalue()
        return FastaSeq(name, seq_bytes)

    def sequence_bytes(self, info: FastaInfo, start, end) -> BytesIO:
        start -= 1  # Switch to Python coordinates
        rpl = info.residues_per_line
        mll = info.max_line_length
        line_end_bytes = mll - rpl

        frst_line = start // rpl
        last_line = (end - 1) // rpl

        frst_offset = start % rpl
        last_offset = end % rpl

        # Seek to the first residue
        fh = self.fasta_fileandle
        fh.seek(info.file_offset + frst_offset + mll * frst_line)

        seq = BytesIO()
        if frst_line == last_line:
            # Sequence fragment is all on one line of the FASTA file
            seq.write(fh.read(end - start))
            return seq
        else:
            # Read sequence to the end of the first line
            seq.write(fh.read(rpl - frst_offset))
            fh.seek(line_end_bytes, 1)

            # Read all the whole lines
            last_whole_line = last_line if last_offset == 0 else last_line - 1
            for _ in range(last_whole_line - frst_line):
                seq.write(fh.read(rpl))
                fh.seek(line_end_bytes, 1)

            # Read any sequence on the last line
            if last_offset:
                seq.write(fh.read(last_offset))
            return seq
