class FastaInfo:
    """
    Represents a line on a `.fai` FASTA index file.
    """
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

    def fai_row(self, name) -> str:
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
