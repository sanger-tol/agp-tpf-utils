import re

from tola.assembly.assembly import Assembly
from tola.assembly.fragment import Fragment
from tola.assembly.gap import Gap
from tola.assembly.scaffold import Scaffold


def parse_agp(file, name):
    asm = Assembly(name)
    scaffold = None
    scaffold_name = ""
    strand_dict = {"+": 1, "-": -1}
    for line in file:
        if re.match(r"\s*$", line):
            # Skip blank lines
            continue

        if line.startswith("##"):
            # We don't care about ## lines, do we?
            continue

        if line.startswith('#'):
            # This pattern match skips blank header lines
            if h := re.match(r"[#\s]+(.+)", line):
                asm.add_header_line(h.group(1))
            continue

        # Remove possible line endings and split on tabs
        fields = line.rstrip("\r\n").split("\t")

        if fields[0] != scaffold_name:
            scaffold_name = fields[0]
            scaffold = Scaffold(scaffold_name)
            asm.add_scaffold(scaffold)

        if fields[4] in ("D", "U"):
            scaffold.add_row(
                Gap(
                    length=fields[5],
                    gap_type=fields[6],
                )
            )
        else:
            scaffold.add_row(
                Fragment(
                    name=fields[5],
                    start=fields[6],
                    end=fields[7],
                    strand=strand_dict[fields[8]],
                    # Add the tenth field as meta if present
                    meta=(fields[9] if len(fields) > 9 else None),
                )
            )

    return asm


def parse_tpf(file, name):
    pass
