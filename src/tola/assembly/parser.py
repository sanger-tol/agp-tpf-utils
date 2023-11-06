import re
import string
from functools import cache

from tola.assembly.assembly import Assembly
from tola.assembly.fragment import Fragment
from tola.assembly.gap import Gap
from tola.assembly.scaffold import Scaffold


def parse_agp(file, name):
    asm = Assembly(name)
    scaffold = None
    scaffold_name = ""
    strand_dict = {"?": 0, "+": 1, "-": -1}
    for line in file:
        if re.match(r"\s*$", line):
            # Skip blank lines
            continue

        if line.startswith("##"):
            # We don't care about ## lines, do we?
            continue

        if line.startswith("#"):
            # This pattern match skips blank header lines
            if h := re.match(r"[#\s]+(.+)", line):
                asm.add_header_line(h.group(1))
            continue

        # AGP from pretext has stray tab on end of line, but rstrip() will
        # remove it.
        fields = line.rstrip().split("\t")

        if fields[0] != scaffold_name:
            scaffold_name = fields[0]
            scaffold = Scaffold(scaffold_name)
            asm.add_scaffold(scaffold)

        if fields[4] in ("U", "N"):
            scaffold.add_row(
                Gap(
                    length=fields[5],
                    gap_type=fields[6],
                ),
            )
        else:
            scaffold.add_row(
                Fragment(
                    name=fields[5],
                    start=fields[6],
                    end=fields[7],
                    strand=strand_dict[fields[8]],
                    # Tenth fields onwards added as tags metadata
                    tags=tuple(fields[9:]),
                ),
            )

    return asm


def parse_tpf(file, name):
    asm = Assembly(name)
    scaffold = None
    scaffold_name = ""
    strand_dict = {"PLUS": 1, "MINUS": -1}
    gap_type_dict = {
        "TYPE-2": "scaffold",
        "TYPE-3": "contig",
    }
    tr = lowercase_and_dash_to_underscore()
    for line in file:
        if re.match(r"\s*$", line):
            # Skip blank lines
            continue

        if line.startswith("#"):
            # This pattern match skips blank header lines
            if h := re.match(r"[#\s]+(.+)", line):
                asm.add_header_line(h.group(1))
            continue

        # Remove possible line endings and split on tabs
        fields = line.rstrip("\r\n").split("\t")

        if fields[0] == "GAP":
            if scaffold:
                scaffold.add_row(
                    Gap(
                        length=fields[2],
                        gap_type=gap_type_dict.get(
                            fields[1],
                            fields[1].translate(tr),
                        ),
                    ),
                )
            else:
                msg = f"Gap line before first sequence fragment: '{line}'"
                raise ValueError(msg)
        elif len(fields) == 4:
            if fields[2] != scaffold_name:
                scaffold_name = fields[2]
                scaffold = Scaffold(scaffold_name)
                asm.add_scaffold(scaffold)
            if m := re.match(r"(.+):(\d+)-(\d+)$", fields[1]):
                scaffold.add_row(
                    Fragment(
                        name=m.group(1),
                        start=m.group(2),
                        end=m.group(3),
                        strand=strand_dict[fields[3]],
                    ),
                )
            else:
                msg = f"Unexpected name format '{fields[1]}'"
                raise ValueError(msg)
        else:
            msg = f"Wrong field count {len(fields)}; 4 expected in line: '{line}'"
            raise ValueError(msg)

    return asm


@cache
def lowercase_and_dash_to_underscore():
    return str.maketrans(
        string.ascii_uppercase + "-",
        string.ascii_lowercase + "_",
    )
