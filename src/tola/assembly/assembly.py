import io
import re
import textwrap

from tola.assembly.fragment import Fragment
from tola.assembly.gap import Gap
from tola.assembly.scaffold import Scaffold


class Assembly:
    def __init__(self, name, scaffolds=None):
        self.name = str(name)
        if scaffolds:
            self.scaffolds = scaffolds
        else:
            self.scaffolds = []

    def __repr__(self):
        txt = io.StringIO()
        txt.write(
            f"{self.__class__.__name__}(\n"
            + f"    name='{self.name}',\n"
            + f"    scaffolds=[\n"
        )
        for scffld in self.scaffolds:
            txt.write(textwrap.indent(f"{scffld!r},\n", "        "))
        txt.write("    ],\n)")
        return txt.getvalue()

    def __str__(self):
        txt = io.StringIO()
        txt.write(f"{self.__class__.__name__}: {self.name}\n")
        for scffld in self.scaffolds:
            txt.write(textwrap.indent(str(scffld), "  "))
        return txt.getvalue()

    @staticmethod
    def name_natural_key(obj):
        return tuple(
            int(x) if i % 2 else x for i, x in enumerate(re.split(r"(\d+)", obj.name))
        )

    def sort_scaffolds_by_name(self):
        self.scaffolds.sort(key=self.name_natural_key)
