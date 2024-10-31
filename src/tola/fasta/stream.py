from io import BufferedIOBase, BytesIO

from tola.assembly.assembly import Assembly
from tola.assembly.gap import Gap
from tola.assembly.scaffold import Scaffold
from tola.fasta.index import FastaIndex


class FastaStream:
    def __init__(
        self,
        out: BufferedIOBase,
        index: FastaIndex,
        line_length=60,
        gap_character=b"N",
    ):
        self.out = out
        self.index = index
        self.line_length = line_length
        self.gap_character = gap_character

    def write_assembly(self, assembly: Assembly):
        for scffld in assembly.scaffolds:
            self.write_scaffold(scffld)

    def write_scaffold(self, scaffold: Scaffold):
        out = self.out
        fai = self.index
        line_length = self.line_length
        want = line_length

        out.write(f">{scaffold.name}\n".encode().lower())
        for row in scaffold.rows:
            itr = (
                self.gap_seq(row)
                if isinstance(row, Gap)
                else fai.get_sequence(row)
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

    def gap_seq(self, gap: Gap):
        yield BytesIO(self.gap_character * gap.length)
