import string

from functools import cache

from tola.assembly.gap import Gap


def format_agp(asm, file):
    STRAND_STR = "?", "+", "-"
    for line in asm.header:
        file.write(f"# {line}\n")
    for scffld in asm.scaffolds:
        scffld_name = scffld.name
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
                    )
                )
            else:
                cols.extend(
                    (
                        "W",
                        str(row.start),
                        str(row.end),
                        STRAND_STR[row.strand],
                    )
                )
                if m := row.meta:
                    cols.append(m)
            file.write("\t".join(cols))
            file.write("\n")


def format_tpf(asm, file):
    STRAND_STR = "UNKNOWN", "PLUS", "MINUS"
    gap_type_dict = {
        "scaffold": "TYPE-2",
        "contig": "TYPE-3",
    }
    tr = uppercase_and_underscore_to_dash()
    for line in asm.header:
        file.write(f"## {line}\n")
    for scffld in asm.scaffolds:
        scffld_name = scffld.name
        for row in scffld.rows:
            if isinstance(row, Gap):
                file.write(
                    "\t".join(
                        (
                            'GAP',
                            gap_type_dict.get(
                                row.gap_type,
                                row.gap_type.translate(tr),
                            ),
                            str(row.length),
                        )
                    )
                )
            else:
                file.write(
                    "\t".join(
                        (
                            "?",
                            f"{row.name}:{row.start}-{row.end}",
                            scffld_name,
                            STRAND_STR[row.strand],
                        )
                    )
                )
            file.write("\n")


@cache
def uppercase_and_underscore_to_dash():
    return str.maketrans(
        string.ascii_lowercase + '_',
        string.ascii_uppercase + '-',
    )
