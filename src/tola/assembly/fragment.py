class Fragment:
    __slots__ = "_name", "_start", "_end", "_strand", "_tags"

    def __init__(self, name, start, end, strand, tags=()):
        self._name = str(name)
        self._start = int(start)
        self._end = int(end)
        self._strand = int(strand)
        self._tags = tags

        if self.strand not in (0, 1, -1):
            msg = f"strand '{self.strand}' should be one of: 0, 1, -1"
            raise ValueError(msg)

        if self.start > self.end:
            msg = f"start '{self.start}' must be <= end 'self.end'"
            raise ValueError(msg)

    @property
    def name(self):
        return self._name

    @property
    def start(self):
        return self._start

    @property
    def end(self):
        return self._end

    @property
    def strand(self):
        return self._strand

    @property
    def tags(self):
        """tuple of tags data, empty if there is none"""
        return self._tags

    @property
    def length(self):
        return self.end - self.start + 1

    STRAND_STR = ".", "+", "-"

    @property
    def strand_str(self):
        """
        Returns a string depending on the strand:
            0: "."
            1: "+"
           -1: "-"
        """
        return self.STRAND_STR[self.strand]

    def attr_values(self):
        return tuple(self.__getattribute__(x) for x in self.__slots__)

    def __eq__(self, othr):
        if self is othr:
            return True
        else:
            return self.attr_values() == othr.attr_values()

    def __str__(self):
        return f"{self.name}:{self.start}-{self.end}({self.strand_str})" + (
            (" " + " ".join(self.tags)) if self.tags else ""
        )

    def __repr__(self):
        return (
            f"{self.__class__.__name__}(name='{self.name}',"
            f" start={self.start}, end={self.end}, strand={self.strand}"
            + (f", tags={self.tags})" if self.tags else ")")
        )

    def overlaps(self, othr):
        if self.name != othr.name:
            return False
        if self.end >= othr.start and self.start <= othr.end:
            return True
        return False

    def abuts(self, othr):
        if self.name != othr.name:
            return False
        if self.end + 1 == othr.start or othr.end + 1 == self.start:
            return True
        return False

    def reverse(self):
        return self.__class__(self.name, self.start, self.end, -1 * self.strand)

    def rename(self, new_name):
        return self.__class__(new_name, self.start, self.end, self.strand)
