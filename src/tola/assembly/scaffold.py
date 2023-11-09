import io

from tola.assembly.fragment import Fragment
from tola.assembly.gap import Gap


class Scaffold:
    def __init__(self, name, rows=None):
        self.name = str(name)
        if rows:
            self.rows = [*rows]
        else:
            self.rows = []

    def __repr__(self):
        txt = io.StringIO()
        txt.write(
            f"{self.__class__.__name__}(\n"
            + f"    name='{self.name}',\n"
            + "    rows=[\n",
        )
        for row in self.rows:
            txt.write(f"        {row!r},\n")
        txt.write("    ],\n)")
        return txt.getvalue()

    def __str__(self):
        txt = io.StringIO()
        txt.write(f"{self.name}\n")
        p = 0
        for row in self.rows:
            txt.write(f"  {p + 1 :11d} {p + row.length :11d}  {row}\n")
            p += row.length
        return txt.getvalue()

    def add_row(self, row):
        self.rows.append(row)

    @property
    def last_row_is_fragment(self):
        if self.rows:
            return isinstance(self.rows[-1], Fragment)
        else:
            return False

    def fragments(self):
        for row in self.rows:
            if isinstance(row, Fragment):
                yield row

    def idx_fragments(self):
        for i, row in enumerate(self.rows):
            if isinstance(row, Fragment):
                yield i, row

    def gaps(self):
        for row in self.rows:
            if isinstance(row, Gap):
                yield row

    def idx_gaps(self):
        for i, row in enumerate(self.rows):
            if isinstance(row, Gap):
                yield i, row

    def reverse(self):
        new = self.__class__(self.name)
        new.rows = self.rows[::-1]
        for i, frag in new.idx_fragments():
            new.rows[i] = frag.reverse()
        return new

    def append_scaffold(self, othr, gap=None):
        # Add a gap if it is not the first row
        if gap and self.rows:
            self.add_row(gap)
        self.rows.extend(othr.rows)
