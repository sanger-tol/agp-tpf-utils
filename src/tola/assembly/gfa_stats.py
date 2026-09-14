import re
import subprocess
import sys
from collections.abc import Iterable
from functools import cached_property
from io import StringIO
from pathlib import Path
from typing import Literal

from tola.assembly.file_utils import get_output_filehandle


class GfaStatsError(Exception):
    """Error running `gfastats` process"""


class BeforeAfterStats:
    def __init__(self, before: Path, after: Path):
        self.__before = before
        self.__after = after

        asm_root = re.sub(r"\.fa(\.gz)?$", "", after.name)
        self.__gfa_out = after.parent / f"{asm_root}.gfastats"
        self.__vgp_out = after.parent / f"{asm_root}.gather_vgp_stats"

    def write_stats(self, clobber: bool = False):
        before = GfaStats(self.__before).run()
        after = GfaStats(self.__after).run()

        gfa_io = get_output_filehandle(self.__gfa_out, clobber)
        gfa_io.write(before.short_stats_string())
        gfa_io.write("\n")
        gfa_io.write(after.short_stats_string())
        gfa_io.close()

        vgp_stats = VgpStats()
        vgp_stats.load_before_stats(before)
        vgp_stats.load_after_stats(after)
        vgp_io = get_output_filehandle(self.__vgp_out, clobber)
        vgp_io.write(str(vgp_stats))
        vgp_io.close()


class Stat:
    """
    An object representing a line from `gfastats` output.
    """
    __slots__ = (
        "gfa_name",
        "col_name",
        "category",
        "measure",
        "value",
    )

    def __init__(self, gfa_name: str, raw_value: str):
        """
        Parses a `gfastats` line's name and value into a useful object.
        """
        self.gfa_name = gfa_name
        py_type = int

        first, *rest = gfa_name.lower().split(" ")
        col_name = "_".join(rest)

        match first:
            case "scaffold" | "contig" | "gap":
                self.category = first
                self.measure = rest[0]
                self.col_name = f"{first}_{rest[0]}"
                if self.measure == "aun":
                    py_type = float
            case "#":
                self.category = rest[0].rstrip("s")
                self.measure = "count"
                self.col_name = rest[0]
            case "total":
                self.category = rest[0]
                self.measure = first
                self.col_name = col_name
            case "average":
                py_type = float
                self.category = rest[0]
                self.measure = "mean_length"
                self.col_name = f"{col_name}_mean"
            case "largest":
                self.category = rest[0]
                self.measure = "longest_length"
                self.col_name = f"{col_name}_longest"
            case "smallest":
                self.category = rest[0]
                self.measure = "shortest_length"
                self.col_name = f"{col_name}_shortest"
            case "base":
                self.category = "contig"
                self.measure = "composition"
                self.col_name = "base_counts_a_c_g_t"
                self.value = [int(x) for x in raw_value.split(":")]
                py_type = list
            case "gc":
                self.category = "contig"
                self.measure = "gc_content"
                self.col_name = "gc_content_percent"
                py_type = float
            case _:
                msg = f"Unexpected line: '{gfa_name}\t{raw_value}'"
                raise GfaStatsError(msg)

        if not hasattr(self, "value"):
            # Removing any commas from the string and set the value using the
            # appropriate type
            self.value = py_type(raw_value.replace(",", ""))

    def __str__(self):
        val = self.value
        if isinstance(val, list):
            val_str = ":".join(str(x) for x in val)
            return f"{self.gfa_name}\t{val_str}"
        else:
            num_fmt = ",.2f" if isinstance(val, float) else ",d"
            return f"{self.gfa_name}\t{self.value:{num_fmt}}"

    def __repr__(self):
        out = StringIO()
        out.write(f"{__class__.__name__}(\n")
        for attr in self.__slots__:
            val = getattr(self, attr) if hasattr(self, attr) else None
            out.write(f"      {attr}={val!r},\n")
        out.write("    ),")
        return out.getvalue()


class GfaStats:
    def __init__(self, asm_file: Path):
        self.__asm_file = asm_file
        self.__stats: list[Stat] = []

    def __repr__(self):
        return (
            "GfaStats(\n"
            f"  asm_file={self.__asm_file!r},\n"
            "  stats=[\n" + "\n".join(f"    {stat!r}" for stat in self.stats()) + "\n])"
        )

    def run(self) -> "GfaStats":
        io = self.gfa_stats
        self.parse(io)
        return self

    def parse(self, io: StringIO) -> None:
        io.seek(0)
        for line in io:
            name, value = line.rstrip().split("\t")
            self.__stats.append(Stat(name, value))

    def stats(self) -> Iterable[Stat]:
        yield from self.__stats

    def short_stats(self) -> Iterable[Stat]:
        for stat in self.__stats:
            yield stat
            if stat.gfa_name == "# paths":
                break

    def stats_string(self) -> str:
        out = StringIO()
        for stat in self.stats():
            out.write(f"{stat}\n")
        return out.getvalue()

    def short_stats_string(self) -> str:
        out = StringIO()
        for stat in self.short_stats():
            out.write(f"{stat}\n")
        return out.getvalue()

    @cached_property
    def gfa_stats(self):
        cmd = [
            "gfastats",
            "--tabular",
            "--threads=2",
            "--locale=C",
            "--nstar-report",
            self.__asm_file,
        ]
        try:
            gfa = subprocess.run(cmd, check=True, capture_output=True, text=True)  # noqa: S603, S607
        except subprocess.CalledProcessError as cpe:
            cmd_str = " ".join(str(x) for x in cmd)
            msg = f"Error running: '{cmd_str}':\n  {cpe.stdout}"
            raise GfaStatsError(msg) from None

        return StringIO(gfa.stdout)


class VgpStats:
    """
    Accumulates and writes VGP before / after stats files such as:

    ```
    scaffolds
    total 926654326 926979538
    count 59 54
    N50 35473549 35473549
    L50 11 11
    N90 21253049 21253049
    L90 24 24

    contigs
    total 926649126 926972438
    count 111 115
    N50 21977639 20662820
    L50 15 15
    N90 7215027 7200072
    L90 41 42
    ```
    """

    def __init__(self):
        self.__stats = {
            "scaffold": VgpStatSet("scaffolds"),
            "contig": VgpStatSet("contigs"),
        }

    def __str__(self) -> str:
        return "\n".join(str(x) for x in self.__stats.values())

    def load_before_stats(self, gfa: GfaStats):
        for stat in gfa.stats():
            self.add_stat(stat, 0)

    def load_after_stats(self, gfa: GfaStats):
        for stat in gfa.stats():
            self.add_stat(stat, 1)

    def add_stat(self, stat: Stat, i: Literal[0, 1]):
        if stat_set := self.__stats.get(stat.category):
            attr = stat.measure
            if re.search(r"^.\d+$", attr):
                attr = attr.upper()
            if not hasattr(stat_set, attr):
                return
            if before_after := getattr(stat_set, attr):
                before_after[i] = stat.value


class VgpStatSet:
    """
    The stats for either "scaffolds" or "contigs"
    """

    __slots__ = "name", "total", "count", "N50", "L50", "N90", "L90"

    def __init__(self, name: str):
        self.name = name
        # Initialise each attribute apart from name
        for attr in self.__slots__[1:]:
            setattr(self, attr, [None, None])

    def __str__(self) -> str:
        out = StringIO()
        for attr in self.__slots__:
            value = getattr(self, attr)
            if attr == "name":
                out.write(value + "\n")
            else:
                out.write(" ".join(str(x) for x in (attr, *value)) + "\n")
        return out.getvalue()
