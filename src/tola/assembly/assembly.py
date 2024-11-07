import io
import re
import textwrap

from tola.assembly.scaffold import Scaffold


class Assembly:
    def __init__(self, name, header=None, scaffolds=None, bp_per_texel=None):
        self.name = str(name)
        self.scaffolds = scaffolds if scaffolds else []
        self.header = header if header else []
        if bp_per_texel:
            self.bp_per_texel = bp_per_texel

    def __repr__(self):
        txt = io.StringIO()
        txt.write(f"{self.__class__.__name__}(\n    name='{self.name}',\n")

        if self.header:
            txt.write("    header=[\n")
            for line in self.header:
                txt.write(f"        '{line}',\n")
            txt.write("    ],\n")

        if self.scaffolds:
            txt.write("    scaffolds=[\n")
            for scffld in self.scaffolds:
                txt.write(textwrap.indent(f"{scffld!r},\n", "        "))
            txt.write("    ],\n)")
        else:
            txt.write(")")

        return txt.getvalue()

    def __str__(self):
        txt = io.StringIO()
        txt.write(f"{self.__class__.__name__}: {self.name}\n")
        for line in self.header:
            txt.write(f"  # {line}\n")
        for scffld in self.scaffolds:
            txt.write("\n")
            txt.write(textwrap.indent(str(scffld), "  "))

        return txt.getvalue()

    def add_header_line(self, txt: str):
        self.header.append(txt)

    def add_scaffold(self, scffld: Scaffold):
        self.scaffolds.append(scffld)

    @property
    def bp_per_texel(self):
        if hasattr(self, "_bp_per_texel"):
            return self._bp_per_texel
        else:
            bpt = None
            for txt in self.header:
                if m := re.match(r"HiC MAP RESOLUTION: ([\d\.]+) bp/texel", txt):
                    bpt = float(m.group(1))
            self._bp_per_texel = bpt
            return bpt

    @bp_per_texel.setter
    def bp_per_texel(self, bp_per_texel: float):
        self._bp_per_texel = bp_per_texel

    @property
    def length(self):
        return sum(s.length for s in self.scaffolds)

    @property
    def fragments_length(self):
        return sum(s.fragments_length for s in self.scaffolds)

    @property
    def gaps_length(self):
        return sum(s.gaps_length for s in self.scaffolds)

    NEMATODE_CHR_INT = {
        "I": 1,
        "II": 2,
        "III": 3,
        "IV": 4,
        # "V" is not included because it is a sex chromosome in the UV system
        #  in plants, and will sort in the correct order within nematode
        #  chromosomes anyway.
    }

    @staticmethod
    def name_natural_key(obj):
        return tuple(
            (Assembly.NEMATODE_CHR_INT.get(x) or int(x)) if i % 2 else x
            for i, x in enumerate(re.split(r"(I+V?|\d+)", obj.name))
        )

    def fragment_junction_set(self):
        junctions = set()
        for scffld in self.scaffolds:
            junctions |= scffld.fragment_junction_set()
        return junctions

    def fragment_junctions_by_asm_prefix(self):
        prefix_junctions = {}
        for scffld in self.scaffolds:
            try:
                first = next(scffld.fragments())
            except StopIteration:
                continue
            m = re.match(r"([A-Za-z]+\d+)_", first.name)
            asm_name = m.group(1).lower() if m else None
            prefix_junctions.setdefault(asm_name, set()).update(
                scffld.fragment_junction_set()
            )
        return prefix_junctions

    def scaffolds_sorted_by_name(self):
        return sorted(self.scaffolds, key=self.name_natural_key)

    def smart_sort_scaffolds(self):
        def smart_sort_key(scaffold: Scaffold):
            return scaffold.rank, self.name_natural_key(scaffold)

        self.scaffolds.sort(key=smart_sort_key)

    def find_overlapping_fragments(self):
        over_pairs = []

        def detect_overlap(v1, v2):
            if v1[0].overlaps(v2[0]):
                over_pairs.append((v1, v2))

        self.all_vs_all_fragments(detect_overlap)
        return over_pairs if over_pairs else None

    def all_vs_all_fragments(self, compare_func):
        frags = []
        for scffld in self.scaffolds:
            frags.extend((x, scffld) for x in scffld.fragments())
        lgth = len(frags)
        for i in range(0, lgth):
            for j in range(i + 1, lgth):
                compare_func(frags[i], frags[j])
