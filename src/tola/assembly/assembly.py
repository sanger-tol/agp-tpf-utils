import io
import re
import textwrap


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

    def add_header_line(self, txt):
        self.header.append(txt)

    def add_scaffold(self, scffld):
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
    def bp_per_texel(self, bp_per_texel):
        self._bp_per_texel = bp_per_texel

    @staticmethod
    def name_natural_key(obj):
        return tuple(
            int(x) if i % 2 else x for i, x in enumerate(re.split(r"(\d+)", obj.name))
        )

    def scaffolds_sorted_by_name(self):
        return sorted(self.scaffolds, key=self.name_natural_key)

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
