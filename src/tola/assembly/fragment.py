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
            msg = f"start '{self.start}' must be <= end '{self.end}'"
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
        return self._end - self._start + 1

    @property
    def key_tuple(self) -> tuple[str, int, int]:
        return self._name, self._start, self._end

    def junction_tuple(self, othr) -> tuple:
        """
        Encodes the positions of two adjacent Fragments in a Scaffold, with
        reverse strand ends encoded by flipping the order of the name and
        coordinate.
        """
        if self.strand == 1:
            if othr.strand == 1:
                #      fwd >>>              fwd >>>
                return self.name, self.end, othr.name, othr.start
            elif othr.strand == -1:
                #      fwd >>>                          <<< rev
                return self.name, self.end, othr.end, othr.name
        elif self.strand == -1:
            if othr.strand == 1:
                #                    <<< rev  fwd >>>
                return self.start, self.name, othr.name, othr.start
            elif othr.strand == -1:
                # For the rev-rev case, junction should match fwd-fwd
                #      rev >>>              rev >>>
                return othr.name, othr.end, self.name, self.start

        msg = f"strand == 0 not supported:\n  {self}\n  {othr}"
        raise ValueError(msg)

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
        return bool(self.end >= othr.start and self.start <= othr.end)

    def overlap_length(self, othr):
        if self.name != othr.name:
            return None

        ovr_start = max(self.start, othr.start)
        ovr_end = min(self.end, othr.end)
        if ovr_start > ovr_end:
            return None
        else:
            return ovr_end - ovr_start + 1

    def abuts(self, othr):
        if self.name != othr.name:
            return False
        return bool(self.end + 1 == othr.start or othr.end + 1 == self.start)

    def gap_between(self, othr):
        """ Returns `None` if no gap, zero if Fragments abut, and the length
            of the gap otherwise. """
        if self.name != othr.name:
            return None

        gap_start = min(self.end, othr.end)
        gap_end = max(self.start, othr.start)
        if gap_start < gap_end:
            return gap_end - gap_start - 1
        else:
            return None

    def reverse(self):
        return self.__class__(
            self.name,
            self.start,
            self.end,
            -1 * self.strand,
            self.tags,
        )

    def rename(self, new_name):
        return self.__class__(
            new_name,
            self.start,
            self.end,
            self.strand,
            self.tags,
        )
