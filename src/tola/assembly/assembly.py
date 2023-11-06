import io
import re
import textwrap

from functools import cached_property

class Assembly:
    def __init__(self, name, header=None, scaffolds=None):
        self.name = str(name)
        self.scaffolds = scaffolds if scaffolds else []
        self.header = header if header else []

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

    @cached_property
    def bp_per_texel(self):
        for txt in self.header:
            if m := re.match(r"HiC MAP RESOLUTION: ([\d\.]+) bp/texel", txt):
                return float(m.group(1))
        return None

    @staticmethod
    def name_natural_key(obj):
        return tuple(
            int(x) if i % 2 else x for i, x in enumerate(re.split(r"(\d+)", obj.name))
        )

    def sort_scaffolds_by_name(self):
        self.scaffolds.sort(key=self.name_natural_key)

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
