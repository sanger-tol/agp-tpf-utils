import io

from tola.assembly.gap import Gap
from tola.assembly.fragment import Fragment
from tola.assembly.scaffold import Scaffold


class OverlapResult(Scaffold):
    def __init__(self, bait, rows, start, end, name=None, tag=None, haplotype=None):
        if not name:
            name = f"matches to {bait.name} {bait.start} to {bait.end}"
        super().__init__(name, rows, tag, haplotype)
        self.bait = bait
        self.start = start
        self.end = end

    def __repr__(self) -> str:
        scffld_repr = super().__repr__()
        rows_start = "    rows=[\n"
        return scffld_repr.replace(
            rows_start,
            (
                f"    bait={self.bait!r},\n"
                + f"    start={self.start},\n"
                + f"    end={self.end},\n"
                + rows_start
            ),
            1,
        )

    def __str__(self) -> str:
        """
        String representation of OverlapResult object useful during
        development.
        """
        txt = io.StringIO()
        txt.write(
            f"{self.name}\n"
            f"  length: {self.length:14_d}\n"
            f"  bait:   {self.bait.length:14_d}  {self.bait}\n"
            f"  diff:   {self.length - self.bait.length:14_d}\n"
            f"  overhang: {self.start_overhang:_d}\n"
        )
        p = self.start - 1
        for row in self.rows:
            txt.write(
                f"    {p + 1 :14_d} {p + row.length :14_d}"
                f" {row.length:11_d}  {row}\n"
            )
            p += row.length
        txt.write(f"  overhang: {self.end_overhang:_d}")
        return txt.getvalue()

    @property
    def length(self) -> int:
        return self.end - self.start + 1

    @property
    def start_overhang(self) -> int:
        return self.bait.start - self.start

    @property
    def end_overhang(self) -> int:
        return self.end - self.bait.end

    @property
    def start_row_bait_overlap(self) -> int:
        """
        Length of the overlap between the bait and the first row.
        """
        start = max(self.bait.start, self.start)
        end = min(
            self.bait.end,
            self.start + self.rows[0].length - 1,  # end of first row
        )

        if end < start:
            return 0
        else:
            return end - start + 1

    @property
    def end_row_bait_overlap(self) -> int:
        """
        Length of the overlap between the bait and the last row.
        """
        start = max(
            self.bait.start,
            self.end - self.rows[-1].length + 1,  # start of last row
        )
        end = min(self.bait.end, self.end)

        if end < start:
            return 0
        else:
            return end - start + 1

    @property
    def length_error(self) -> int:
        return self.length - self.bait.length

    def length_error_in_texels(self, bp_per_texel) -> float:
        return abs(self.length_error) / bp_per_texel

    def fragment_start_if_trimmed(self, frag: Fragment) -> int:
        if frag.strand == 1:
            if self.rows[0] is frag:
                # Forward strand fragment at start of OverlapResult
                return frag.start + self.start_overhang
        else:
            if self.rows[-1] is frag:
                # Reverse strand fragment at end of OverlapResult
                return frag.start + self.end_overhang

        # Start would not be changed
        return frag.start

    def trim_fragment(
        self, trim: Fragment, keep_start=False, keep_end=False
    ) -> Fragment:
        start = trim.start
        end = trim.end

        # Fragment may be at both start and end
        idx = None
        if self.rows[0] is trim:
            # fragment is at the start
            idx = 0
            start_ovr = self.start_overhang
            if start_ovr > 0 and not keep_start:
                if trim.strand == 1:
                    start += start_ovr
                else:
                    end -= start_ovr
                self.start += start_ovr
        if self.rows[-1] is trim:
            # fragment is at the end
            idx = -1
            end_ovr = self.end_overhang
            if end_ovr > 0 and not keep_end:
                if trim.strand == 1:
                    end -= end_ovr
                else:
                    start += end_ovr
                self.end -= end_ovr
        if idx is None:
            msg = f"Fragment {trim} not found in:\n{self}"
            raise ValueError(msg)

        bait_tags = (x for x in self.bait.tags if x != "Painted")
        new = Fragment(trim.name, start, end, trim.strand, ("Cut", *bait_tags))
        self.rows[idx] = new
        return new

    def reverse(self):
        """
        Not implemented - no use case (yet)
        """
        raise NotImplementedError

    def append_scaffold(self):
        """
        Not implemented - no use case (yet)
        """
        raise NotImplementedError

    def to_scaffold(self) -> Scaffold:
        scffld = Scaffold(self.name, self.rows)
        if self.bait.strand == -1:
            return scffld.reverse()
        else:
            return scffld

    def discard_start(self) -> None:
        discard = self.rows.pop(0)
        self.start += discard.length
        while self.rows and isinstance(self.rows[0], Gap):
            gap = self.rows.pop(0)
            self.start += gap.length

    def discard_end(self) -> None:
        discard = self.rows.pop(-1)
        self.end -= discard.length
        while self.rows and isinstance(self.rows[-1], Gap):
            gap = self.rows.pop(-1)
            self.end -= gap.length

    def error_if_start_removed(self) -> int:
        length_if_rem = self.length
        length_if_rem -= self.rows[0].length
        for r in self.rows[1:]:
            if isinstance(r, Gap):
                length_if_rem -= r.length
            else:
                break
        return length_if_rem - self.bait.length

    def error_if_end_removed(self) -> int:
        length_if_rem = self.length
        length_if_rem -= self.rows[-1].length
        for r in self.rows[-2::-1]:  # Step backwards from second to last element
            if isinstance(r, Gap):
                length_if_rem -= r.length
            else:
                break
        return length_if_rem - self.bait.length

    def trim_large_overhangs(self, err_length: int) -> None:
        if len(self.rows) == 1 and self.bait.length > err_length:
            return

        if (
            self.start_overhang > err_length
            and self.start_row_bait_overlap < err_length
        ):
            self.discard_start()

        if self.end_overhang > err_length and self.end_row_bait_overlap < err_length:
            self.discard_end()
