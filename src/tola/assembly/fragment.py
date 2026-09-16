from typing import TypeAlias

Junction: TypeAlias = tuple[str | int, str | int, str | int, str | int]


class Fragment:
    """
    Represents a contiguous piece of DNA containing no "N" characters.  Often
    called a "contig".  A `Fragment` is an immutable object.  Methods which
    change any of its fields return a new `Fragment` object.
    """
    __slots__ = "__name", "__start", "__end", "__strand", "__tags", "__source"

    def __init__(
        self,
        name: str,
        start: int,
        end: int,
        strand: int,
        tags: tuple[str, ...] = (),
        *,
        source: str | None = None,
    ):
        self.__name = name
        self.__start = start
        self.__end = end
        self.__strand = strand
        self.__tags = tags
        self.__source = source

        if self.__strand not in (0, 1, -1):
            msg = f"strand '{self.__strand}' should be one of: 0, 1, -1"
            raise ValueError(msg)

        if self.__start > self.__end:
            msg = f"start '{self.__start}' must be <= end '{self.__end}'"
            raise ValueError(msg)

    @property
    def name(self):
        return self.__name

    @property
    def start(self):
        return self.__start

    @property
    def end(self):
        return self.__end

    @property
    def strand(self):
        return self.__strand

    @property
    def tags(self):
        """tuple of tags data, empty if there are none"""
        return self.__tags

    @property
    def source(self):
        return self.__source

    @property
    def length(self):
        return self.__end - self.__start + 1

    @property
    def key_tuple(self) -> tuple[str, int, int]:
        return self.__name, self.__start, self.__end

    def junction_tuple(self, othr: "Fragment") -> Junction:
        """
        Encodes the positions of two adjacent Fragments in a Scaffold, with
        reverse strand ends encoded by flipping the order of the name and
        coordinate.
        """
        if self.__strand == 1:
            if othr.__strand == 1:
                #      fwd >>>                  fwd >>>
                return self.__name, self.__end, othr.__name, othr.__start
            elif othr.__strand == -1:
                #      fwd >>>                                  <<< rev
                return self.__name, self.__end, othr.__end, othr.__name
        elif self.__strand == -1:
            if othr.__strand == 1:
                #                        <<< rev  fwd >>>
                return self.__start, self.__name, othr.__name, othr.__start
            elif othr.__strand == -1:
                # For the rev-rev case, junction should match fwd-fwd
                #      rev >>>                  rev >>>
                return othr.__name, othr.__end, self.__name, self.__start

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
        return self.STRAND_STR[self.__strand]

    def attr_values(self):
        return (
            self.__name,
            self.__start,
            self.__end,
            self.__strand,
            self.__tags,
            self.__source,
        )

    def __eq__(self, othr):
        if self is othr:
            return True
        else:
            return self.attr_values() == othr.attr_values()

    def __str__(self):
        return (
            f"{self.__name}:{self.__start}-{self.__end}({self.strand_str})"
            + ((" " + " ".join(self.__tags)) if self.__tags else "")
            + (f" source={self.__source!r}" if self.__source else "")
        )

    def __repr__(self):
        return (
            f"{self.__class__.__name__}(name='{self.__name}',"
            f" start={self.__start}, end={self.__end}, strand={self.__strand}"
            + (f", tags={self.__tags!r})" if self.__tags else ")")
            + (f" source={self.__source!r}" if self.__source else "")
        )

    def overlaps(self, othr: "Fragment"):
        if self.__name != othr.__name:
            return False
        return bool(self.__end >= othr.__start and self.__start <= othr.__end)

    def overlap_length(self, othr: "Fragment"):
        if self.__name != othr.__name:
            return None

        ovr_start = max(self.__start, othr.__start)
        ovr_end = min(self.__end, othr.__end)
        if ovr_start > ovr_end:
            return None
        else:
            return ovr_end - ovr_start + 1

    def abuts(self, othr: "Fragment"):
        if self.__name != othr.__name:
            return False
        return bool(self.__end + 1 == othr.__start or othr.__end + 1 == self.__start)

    def gap_between(self, othr: "Fragment"):
        """Returns `None` if no gap, zero if Fragments abut, and the length
        of the gap otherwise."""
        if self.__name != othr.__name:
            return None

        gap_start = min(self.__end, othr.__end)
        gap_end = max(self.__start, othr.__start)
        if gap_start < gap_end:
            return gap_end - gap_start - 1
        else:
            return None

    def reverse(self):
        return self.__class__(
            self.__name,
            self.__start,
            self.__end,
            -1 * self.__strand,
            self.__tags,
        )

    def rename(self, new_name):
        return self.__class__(
            new_name,
            self.__start,
            self.__end,
            self.__strand,
            self.__tags,
        )
