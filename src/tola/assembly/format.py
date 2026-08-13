import string
from functools import cache
from typing import TextIO

from tola.assembly.assembly import Assembly
from tola.assembly.gap import Gap


class FormatAssemblyError(Exception):
    """
    Error formatting an Assembly
    """


def format_agp(asm: Assembly, file: TextIO):
    STRAND_STR = "?", "+", "-"  # noqa: N806
    for line in asm.header:
        file.write(f"# {line}\n")

    seen = set()
    for scffld in asm.scaffolds:
        scffld_name = scffld.name

        # Check for duplicate Scaffold names
        if scffld_name in seen:
            msg = f"More than one Scaffold named {scffld_name!r}"
            raise FormatAssemblyError(msg)
        seen.add(scffld_name)

        p = 0
        for i, row in enumerate(scffld.rows):
            cols = [
                scffld_name,
                str(p + 1),
                str(p + row.length),
                str(i + 1),
            ]
            p += row.length
            if isinstance(row, Gap):
                cols.extend(
                    (
                        "U",
                        str(row.length),
                        str(row.gap_type),
                        "yes",
                        "proximity_ligation",
                    ),
                )
            else:
                cols.extend(
                    (
                        "W",
                        row.name,
                        str(row.start),
                        str(row.end),
                        STRAND_STR[row.strand],
                    ),
                )
                if m := row.tags:
                    cols.extend(m)
            file.write("\t".join(cols))
            file.write("\n")


def format_tpf(asm: Assembly, file: TextIO):
    STRAND_STR = "UNKNOWN", "PLUS", "MINUS"  # noqa: N806
    gap_type_dict = {
        "scaffold": "TYPE-2",
        "contig": "TYPE-3",
    }
    tr = uppercase_and_underscore_to_dash()
    for line in asm.header:
        file.write(f"## {line}\n")

    seen = set()
    for scffld in asm.scaffolds:
        scffld_name = scffld.name

        # Check for duplicate Scaffold names
        if scffld_name in seen:
            msg = f"More than one Scaffold named {scffld_name!r}"
            raise FormatAssemblyError(msg)
        seen.add(scffld_name)

        for row in scffld.rows:
            if isinstance(row, Gap):
                file.write(
                    "\t".join(
                        (
                            "GAP",
                            gap_type_dict.get(
                                row.gap_type,
                                row.gap_type.translate(tr),
                            ),
                            str(row.length),
                        ),
                    ),
                )
            else:
                file.write(
                    "\t".join(
                        (
                            "?",
                            f"{row.name}:{row.start}-{row.end}",
                            scffld_name,
                            STRAND_STR[row.strand],
                        ),
                    ),
                )
            file.write("\n")


@cache
def uppercase_and_underscore_to_dash() -> dict[int, int]:
    return str.maketrans(
        string.ascii_lowercase + "_",
        string.ascii_uppercase + "-",
    )
