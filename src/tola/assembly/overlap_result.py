import io

from tola.assembly.gap import Gap
from tola.assembly.scaffold import Scaffold


class OverlapResult(Scaffold):
    def __init__(self, bait, rows, start, end, name=None):
        if not name:
            name = f"matches to {bait.name} {bait.start} to {bait.end}"
        super().__init__(name, rows)
        self.bait = bait
        self.start = start
        self.end = end

    def __str__(self):
        txt = io.StringIO()
        txt.write(
            f"bait: {self.bait}\nlength: {self.length}\n"
            f"overhang: {self.start_overhang}\n"
        )
        p = self.start - 1
        for row in self.rows:
            txt.write(f"  {p + 1 :11d} {p + row.length :11d}  {row}\n")
            p += row.length
        txt.write(f"overhang: {self.end_overhang}\n")
        return txt.getvalue()

    @property
    def length(self):
        return self.end - self.start + 1

    @property
    def start_overhang(self):
        return self.bait.start - self.start

    @property
    def end_overhang(self):
        return self.end - self.bait.end

    def to_scaffold(self):
        scffld = Scaffold(self.name, self.rows)
        if self.bait.strand == -1:
            return scffld.reverse()
        else:
            return scffld

    def remove_leading_and_trailing_gaps(self):
        while self.rows and isinstance(self.rows[0], Gap):
            discard = self.rows.pop(0)
            self.start += discard.length
        while self.rows and isinstance(self.rows[-1], Gap):
            discard = self.rows.pop(-1)
            self.end -= discard.length

    def trim_large_overhangs(self, err_length):
        if (
            self.start_overhang > err_length
            and self.bait.overlap_length(self.rows[0]) < err_length
        ):
            discard = self.rows.pop(0)
            self.start += discard.length

        if (
            self.end_overhang > err_length
            and self.bait.overlap_length(self.rows[-1]) < err_length
        ):
            discard = self.rows.pop(-1)
            self.end -= discard.length

        self.remove_leading_and_trailing_gaps()
