from io import BufferedIOBase

from tola.assembly.assembly import Assembly
from tola.assembly.gap import Gap
from tola.assembly.scaffold import Scaffold
from tola.fasta.index import FastaIndex


class FastaCollectionError(Exception):
    """Error using `FastaCollection`"""


class FastaCollection:
    """
    A collection of `FastaIndex` objects keyed by `source` to handle
    multiple source FASTA files.
    """

    def __init__(self, fai: FastaIndex | None = None):
        self.__source_idx: dict[str | None, FastaIndex] = {}
        if fai:
            self.add_faidx(fai)

    def add_faidx(self, idx: FastaIndex) -> None:
        """
        Add a new `FastaIndex` under `source`.  Raises a
        `FastaCollectionError` exception if an entry already exists for
        `source`.
        """
        source = idx.source
        if self.__source_idx.get(source):
            msg = f"Already have a FastaIndex for source = '{source}'"
            raise FastaCollectionError(msg)
        self.__source_idx[source] = idx

    def get_faidx(self, source: str | None) -> FastaIndex:
        """
        Return the `FastaIndex` for `source`.
        """
        return self.__source_idx[source]

    def default_faidx(self) -> FastaIndex:
        """
        Returns either the `FastaIndex` stored under the `None` key, or the
        first one stored.
        """
        if not self.__source_idx:
            msg = "Cannot return a default FastaIndex from an empty collection"
            raise FastaCollectionError(msg)
        return self.__source_idx.get(None) or list(self.__source_idx.values())[0]


class FastaStream:
    """
    Streams the sequence data for an `Assembly` to `out` fetching the data for
    sequence `Fragment`s from a `FastaCollection`.
    """

    def __init__(
        self,
        out: BufferedIOBase,
        index_collection: FastaCollection,
        *,
        line_length: int = 60,
        gap_character: bytes = b"N",
    ):
        self.out = out
        self.index_collection = index_collection
        self.line_length = line_length
        self.gap_character = gap_character

    def write_assembly(self, assembly: Assembly):
        for scffld in assembly.scaffolds:
            self.write_scaffold(scffld)

    def write_scaffold(self, scaffold: Scaffold):
        out = self.out
        coll = self.index_collection
        default_fai = coll.default_faidx()
        line_length = self.line_length
        want = line_length

        out.write(f">{scaffold.name}\n".encode())
        for row in scaffold.rows:
            itr = (
                default_fai.get_gap_iter(row, self.gap_character)
                if isinstance(row, Gap)
                else coll.get_faidx(row.source).get_sequence_iter(row)
            )
            for chunk in itr:
                chunk.seek(0)
                while True:
                    if seq := chunk.read(want):
                        out.write(seq)
                        want -= len(seq)
                        if want == 0:
                            out.write(b"\n")
                            want = line_length
                    else:
                        break

        if want != line_length:
            out.write(b"\n")
