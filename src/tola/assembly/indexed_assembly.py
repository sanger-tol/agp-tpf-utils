from tola.assembly.assembly import Assembly


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

    @staticmethod
    def start_end_at_index(idx, i):
        start = 1 if i == 0 else 1 + idx[i - 1]
        end = idx[i]
        return start, end

    def fetch_overlaps(self, fragment):
        scffld = self.scaffold_by_name(fragment.name)
        if not scffld.rows:
            msg = f"Scaffold '{scffld.name}' is empty"
            raise ValueError(msg)

        idx = self._scaffold_index.get(fragment.name)
        if not idx:
            msg = f"Scaffold '{scffld.name}' is not indexed."
            raise ValueError(msg)

        bait_start = fragment.start
        bait_end = fragment.end

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

        return scffld.rows[i_ovr : j_ovr + 1]
