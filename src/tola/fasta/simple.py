import io

IUPAC_COMPLEMENT = bytes.maketrans(
    b"ACGTRYMKSWHBVDNacgtrymkswhbvdn",
    b"TGCAYRKMSWDVBHNtgcayrkmswdvbhn",
)


def reverse_complement(seq: bytes):
    return seq[::-1].translate(IUPAC_COMPLEMENT)


def revcomp_bytes_io(seq: io.BytesIO):
    return io.BytesIO(reverse_complement(seq.getvalue()))


class FastaSeq:
    __slots__ = "name", "description", "sequence"

    def __init__(self, name: str, sequence: bytes, description: str = None):
        self.name = name
        self.sequence = sequence
        self.description = description

    @property
    def length(self):
        return len(self.sequence)

    def __str__(self, line_length=60):
        out = io.StringIO()
        out.write(f">{self.name}")
        if desc := self.description:
            out.write(f" {desc}")
        out.write("\n")
        seq_length = self.length
        seq = self.sequence.decode()
        line_count = 1 + ((seq_length - 1) // line_length)
        for i in range(line_count):
            x = i * line_length
            y = min(seq_length, x + line_length)
            out.write(seq[x:y] + "\n")

        return out.getvalue()

    def fasta_bytes(self, line_length=60):
        out = io.BytesIO()
        out.write(b">" + self.name.encode())
        if desc := self.description:
            out.write(b" " + desc.encode())
        out.write(b"\n")
        seq_length = self.length
        seq = self.sequence
        line_count = 1 + ((seq_length - 1) // line_length)
        for i in range(line_count):
            x = i * line_length
            y = min(seq_length, x + line_length)
            out.write(seq[x:y] + b"\n")

        return out.getvalue()

    def rev_comp(self):
        return FastaSeq(
            self.name,
            reverse_complement(self.sequence),
            self.description,
        )
