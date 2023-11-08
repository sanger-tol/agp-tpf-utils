from tola.assembly.assembly import Assembly
from tola.assembly.gap import Gap
from tola.assembly.overlap_result import OverlapResult


class IndexedAssembly(Assembly):
    def __init__(self, name, header=None, scaffolds=None):
        self.name = str(name)
        self.header = header if header else []
        self._scaffold_dict = {}
        self._scaffold_index = {}
        if scaffolds:
            for scffld in scaffolds:
                self.add_scaffold(scffld)

    @classmethod
    def new_from_assembly(cls, asm):
        return cls(asm.name, asm.header, asm.scaffolds)

    @property
    def scaffolds(self):
        return self._scaffold_dict.values()

    def add_scaffold(self, scffld):
        if self._scaffold_dict.get(scffld.name):
            msg = f"Already have Scaffold named '{scffld.name}'"
            raise ValueError(msg)

        # Index rows by their end position in the Scaffold
        end = 0
        idx = []
        for row in scffld.rows:
            end += row.length
            idx.append(end)

        # Store the Scaffold and its index
        self._scaffold_dict[scffld.name] = scffld
        self._scaffold_index[scffld.name] = idx

    def scaffold_by_name(self, name):
        if scffld := self._scaffold_dict.get(name):
            return scffld
        else:
            msg = f"No such Scaffold '{name}' in Assembly '{self.name}'"
            raise ValueError(msg)

    def find_overlaps(self, bait):
        """
        Given a Fragment bait, returns an OverlapResult (a subclass of
        Scaffold) with rows from the Scaffold within the IndexedAssembly
        which overlap. Any leading or trailing Gaps in the overlapping rows
        are removed.
        """
        scffld = self.scaffold_by_name(bait.name)
        if not scffld.rows:
            msg = f"Scaffold '{scffld.name}' is empty"
            raise ValueError(msg)

        idx = self._scaffold_index.get(bait.name)
        if not idx:
            msg = f"Scaffold '{scffld.name}' is not indexed."
            raise ValueError(msg)

        bait_start = bait.start
        bait_end = bait.end

        # Binary search: a = left; m = midpoint; z = right
        a = 0
        z = len(idx)
        ovr = None
        while a < z:
            m = a + ((z - a) // 2)
            scffld_start = 1 if m == 0 else 1 + idx[m - 1]
            scffld_end = idx[m]
            if scffld_end < bait_start:
                # Scaffold row at "m" is to the left of the bait
                a = m + 1
            elif scffld_start > bait_end:
                # Scaffold row at "m" is to the right of the bait
                z = m
            else:
                # Must be overlapping
                ovr = m
                break

        if ovr is None:
            return None

        # The span of overlapping entries may extend to the left or right
        # of "ovr"
        i_ovr = j_ovr = ovr
        for i in range(ovr - 1, -1, -1):  # if ovr == 4, range will be [3, 2, 1, 0]
            if idx[i] < bait_start:
                break
            i_ovr = i
        for j in range(ovr + 1, len(idx)):
            s_start = 1 if j == 0 else 1 + idx[j - 1]
            if s_start > bait_end:
                break
            j_ovr = j

        # Walk start and end pointers back to ignore Gaps on the ends
        while isinstance(scffld.rows[i_ovr], Gap):
            i_ovr += 1
        while isinstance(scffld.rows[j_ovr], Gap):
            j_ovr -= 1
        if not i_ovr <= j_ovr:
            return None

        overlaps = scffld.rows[i_ovr : j_ovr + 1]
        overlap_start = 1 if i_ovr == 0 else 1 + idx[i_ovr - 1]
        overlap_end = idx[j_ovr]

        return OverlapResult(
            bait=bait,
            start=overlap_start,
            end=overlap_end,
            rows=overlaps,
        )
