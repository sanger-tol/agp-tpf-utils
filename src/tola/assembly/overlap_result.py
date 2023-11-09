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
            f"difference: {self.length - self.bait.length}\n"
            f"overhang: {self.start_overhang}\n"
        )
        p = self.start - 1
        for row in self.rows:
            txt.write(f" {p + 1 :11d} {p + row.length :11d} {row.length:9d}  {row}\n")
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

    def has_problem_overhang(self, bp_per_texel):
        if (
            abs(self.start_overhang) > bp_per_texel
            or abs(self.end_overhang) > bp_per_texel
        ):
            return True
        else:
            return False

    def discard_start(self):
        discard = self.rows.pop(0)
        self.start += discard.length
        while self.rows and isinstance(self.rows[0], Gap):
            gap = self.rows.pop(0)
            self.start += gap.length

    def discard_end(self):
        discard = self.rows.pop(-1)
        self.end -= discard.length
        while self.rows and isinstance(self.rows[-1], Gap):
            gap = self.rows.pop(-1)
            self.end -= gap.length

    def error_increase_if_start_removed(self):
        length_if_rem = self.length
        length_if_rem -= self.rows[0].length
        for r in self.rows[1:]:
            if isinstance(r, Gap):
                length_if_rem -= r.length
            else:
                break
        return self._error_delta(length_if_rem)

    def error_increase_if_end_removed(self):
        length_if_rem = self.length
        length_if_rem -= self.rows[-1].length
        for r in self.rows[-2::-1]:  # Step backwards from second to last element
            if isinstance(r, Gap):
                length_if_rem -= r.length
            else:
                break
        return self._error_delta(length_if_rem)

    def _error_delta(self, length_if_rem):
        before = abs(self.length - self.bait.length)
        after = abs(length_if_rem - self.bait.length)
        return after - before  # More negative is better

    def trim_large_overhangs(self, err_length):
        if (
            self.start_overhang > err_length
            and self.bait.overlap_length(self.rows[0]) < err_length
        ):
            self.discard_start()

        if (
            self.end_overhang > err_length
            and self.bait.overlap_length(self.rows[-1]) < err_length
        ):
            self.discard_end()
