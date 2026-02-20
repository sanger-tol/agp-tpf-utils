import io
from collections.abc import Generator
from typing import Self

from tola.assembly.fragment import Fragment, Junction
from tola.assembly.gap import Gap


class Scaffold:
    def __init__(
        self,
        name,
        rows: list[Fragment | Gap] | None = None,
        tag=None,
        haplotype=None,
        rank=0,
        original_name=None,
        original_tags: set[str] | None = None,
    ):
        self.name = str(name)
        if rows:
            self.rows = [*rows]
        else:
            self.rows: list[Fragment | Gap] = []
        self.tag: str | None = tag
        self.haplotype: str | None = haplotype
        self.rank: int = rank
        self.original_name: str | None = original_name
        self.original_tags: set[str] | None = original_tags

    def __repr__(self):
        txt = io.StringIO()
        txt.write(f"{self.__class__.__name__}(\n    name='{self.name}',\n")
        if orig := self.original_name:
            txt.write(f"    original_name='{orig}',\n")
        if orig_tags := self.original_tags:
            txt.write(f"    original_tags={sorted(orig_tags)!r},\n")
        if rnk := self.rank:
            txt.write(f"    rank='{rnk}',\n")
        txt.write("    rows=[\n")
        for row in self.rows:
            txt.write(f"        {row!r},\n")
        txt.write("    ],\n)")
        return txt.getvalue()

    def __str__(self):
        txt = io.StringIO()
        txt.write(f"{self.name}")
        if (orig := self.original_name) and orig != self.name:
            txt.write(f" ({orig})")
        if orig_tags := self.original_tags:
            txt.write(f" original_tags={sorted(orig_tags)!r}")
        if rnk := self.rank:
            txt.write(f" rank={rnk}")
        txt.write("\n")
        for row in self.rows:
            txt.write(f"  {row.length:14_d}  {row}\n")
        return txt.getvalue()

    def add_row(self, row: Fragment | Gap):
        self.rows.append(row)

    @property
    def length(self) -> int:
        return sum(r.length for r in self.rows)

    @property
    def fragments_length(self) -> int:
        return sum(f.length for f in self.fragments())

    @property
    def gaps_length(self) -> int:
        return sum(g.length for g in self.gaps())

    @property
    def last_row_is_fragment(self) -> bool:
        if self.rows:
            return isinstance(self.rows[-1], Fragment)
        else:
            return False

    def fragments(self) -> Generator[Fragment]:
        for row in self.rows:
            if isinstance(row, Fragment):
                yield row

    def idx_fragments(self) -> Generator[tuple[int, Fragment]]:
        for i, row in enumerate(self.rows):
            if isinstance(row, Fragment):
                yield i, row

    def gaps(self) -> Generator[Fragment]:
        for row in self.rows:
            if isinstance(row, Gap):
                yield row

    def idx_gaps(self) -> Generator[tuple[int, Gap]]:
        for i, row in enumerate(self.rows):
            if isinstance(row, Gap):
                yield i, row

    def fragment_tags(self) -> set[str]:
        tag_set = set()
        for frag in self.fragments():
            for t in frag.tags:
                tag_set.add(t)
        return tag_set

    def reverse(self) -> Self:
        new = self.__class__(
            self.name,
            tag=self.tag,
            haplotype=self.haplotype,
            rank=self.rank,
            original_name=self.original_name,
            original_tags=self.original_tags,
        )
        new.rows = self.rows[::-1]
        for i, frag in new.idx_fragments():
            new.rows[i] = frag.reverse()
        return new

    def append_scaffold(self, othr, gap=None):
        # Add a gap if it is not the first row
        if gap and self.rows:
            self.add_row(gap)
        self.rows.extend(othr.rows)

    def fragment_junction_set(self) -> set[Junction]:
        junctions = set()
        itr = self.fragments()

        try:
            prev = next(itr)
        except StopIteration:
            return junctions

        for this in itr:
            junctions.add(prev.junction_tuple(this))
            prev = this

        return junctions
