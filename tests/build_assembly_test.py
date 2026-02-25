import math
import random
import sys

import pytest

from tola.assembly.assembly import Assembly
from tola.assembly.build_assembly import BuildAssembly, LongScaffoldCuttingError
from tola.assembly.fragment import Fragment
from tola.assembly.gap import Gap
from tola.assembly.indexed_assembly import IndexedAssembly
from tola.assembly.naming_utils import ChrGroup, ChrNamer, ScaffoldNamer, TaggingError
from tola.assembly.scaffold import Scaffold


def test_multi_chr_list():
    assert ChrGroup.multi_chr_list("3", 1) == ["3"]
    assert ChrGroup.multi_chr_list("3", 3) == ["3A", "3B", "3C"]


def list_chr_naming_tests():
    for test_data in [
        {
            # Changed to allow consecutive chromosomes in the first haplotype
            "input": [
                ("S1", "Hap1", 2_000_000, "S1"),
                ("S2", "Hap1", 1_000_000, "S2"),
                ("S3", "Hap3", 1_000_000, "S3"),
                ("S4", "Hap3", 2_000_000, "S4"),
            ],
            "expected": {
                "Hap1": [
                    ("SUPER_1A_Hap1", "Hap1", 2_000_000, "S1"),
                    ("SUPER_1B_Hap1", "Hap1", 1_000_000, "S2"),
                ],
                "Hap3": [
                    ("SUPER_1A_Hap3", "Hap3", 1_000_000, "S3"),
                    ("SUPER_1B_Hap3", "Hap3", 2_000_000, "S4"),
                ],
            },
        },
        {
            "input": [
                ("S1", "Hap1", 3_000_000, "S1", "Singleton"),
                ("S1_unloc_1", "Hap1", 12_000, "S1"),
                ("S2", "Hap1", 2_000_000, "S2"),
                ("S3", "Hap3", 2_000_000, "S3"),
                ("S4", "Hap1", 1_000_000, "S4", "Singleton"),
            ],
            "expected": {
                "Hap1": [
                    ("SUPER_1_Hap1", "Hap1", 3_000_000, "S1"),
                    ("SUPER_1_Hap1_unloc_1", "Hap1", 12_000, "S1"),
                    ("SUPER_2_Hap1", "Hap1", 2_000_000, "S2"),
                    ("SUPER_3_Hap1", "Hap1", 1_000_000, "S4"),
                ],
                "Hap3": [
                    ("SUPER_2_Hap3", "Hap3", 2_000_000, "S3"),
                ],
            },
        },
        {
            "input": [
                # SUPER_2
                ("P1", None, 10_000_000, "P1"),
                ("P1_unloc_2", None, 9_000, "P1"),
                ("P1_unloc_1", None, 10_000, "P1"),
                # SUPER_1
                ("P2", None, 20_000_000, "P2"),
            ],
            "expected": {
                None: [
                    ("SUPER_1", None, 20_000_000, "P2"),
                    ("SUPER_2", None, 10_000_000, "P1"),
                    ("SUPER_2_unloc_1", None, 10_000, "P1"),
                    ("SUPER_2_unloc_2", None, 9_000, "P1"),
                ]
            },
        },
        {
            "input": [
                # SUPER_2
                ("S4", "Hap1", 20_000_000, "S4"),
                ("S5", "Hap2", 10_000_000, "S5"),
                ("S5_unloc_2", "Hap2", 10_000, "S5"),
                ("S5_unloc_1", "Hap2", 19_000, "S5"),
                # SUPER_3
                ("S1", "Hap1", 10_000_000, "S1"),
                ("S2", "Hap2", 7_000_000, "S2"),
                ("S3", "Hap2", 6_000_000, "S3"),
                ("S3_unloc_1", "Hap2", 12_000, "S3"),
                # SUPER_1
                ("S6", "Hap1", 40_000_000, "S6"),
                ("S7", "Hap2", 2_000_000, "S7"),
            ],
            "expected": {
                "Hap1": [
                    ("SUPER_1_Hap1", "Hap1", 40_000_000, "S6"),
                    ("SUPER_2_Hap1", "Hap1", 20_000_000, "S4"),
                    ("SUPER_3_Hap1", "Hap1", 10_000_000, "S1"),
                ],
                "Hap2": [
                    ("SUPER_1_Hap2", "Hap2", 2_000_000, "S7"),
                    ("SUPER_2_Hap2", "Hap2", 10_000_000, "S5"),
                    ("SUPER_2_Hap2_unloc_1", "Hap2", 19_000, "S5"),
                    ("SUPER_2_Hap2_unloc_2", "Hap2", 10_000, "S5"),
                    ("SUPER_3A_Hap2", "Hap2", 7_000_000, "S2"),
                    ("SUPER_3B_Hap2", "Hap2", 6_000_000, "S3"),
                    ("SUPER_3B_Hap2_unloc_1", "Hap2", 12_000, "S3"),
                ],
            },
        },
    ]:
        yield (
            test_data["input"],
            test_data.get("expected", {}),
            test_data.get("exception"),
        )


@pytest.mark.parametrize("input_data,expected,exception", list_chr_naming_tests())
def test_chr_namer(input_data, expected, exception):
    if exception:
        exptn_type, exptn_pattern = exception
        with pytest.raises(exptn_type, match=exptn_pattern):
            run_chr_namer(input_data, expected)
    else:
        run_chr_namer(input_data, expected)


def run_chr_namer(input_data, expected):
    cn = ChrNamer()
    assemblies = {}
    for name, haplotype, length, orig, *tags in input_data:
        asm = assemblies.setdefault(haplotype, Assembly(haplotype))
        scffld = Scaffold(
            name,
            rows=[Fragment(name, 1, length, 1)],
            original_name=orig,
            original_tags=set(tags),
            rank=1,
        )
        print(scffld, file=sys.stderr)
        asm.add_scaffold(scffld)
        cn.add_scaffold(haplotype, scffld)
    cn.name_chromosomes()
    for haplotype, asm in assemblies.items():
        asm.smart_sort_scaffolds()
        fingerprint = [
            (x.name, haplotype, x.fragments_length, x.original_name)
            for x in asm.scaffolds
        ]
        assert fingerprint == expected[haplotype]


def test_unpainted_unloc():
    namer = ScaffoldNamer()
    with pytest.raises(
        TaggingError,
        match=r"Unloc in unpainted scaffold 'Scaffold_7': Bad_Unloc:1-100\(-\) Unloc",
    ):
        namer.label_scaffold(
            Scaffold("TEST"),
            Fragment("Bad_Unloc", 1, 100, -1, ("Unloc",)),
            set(),  # No "Painted" tag in scaffold tags
            "Scaffold_7",
        )


def test_no_coord_changes():
    ia1 = make_random_assembly(seed="Random assembly")

    p1 = Assembly("pretext")
    for out_n, scf in enumerate(ia1.scaffolds):
        ns = Scaffold(f"New_{out_n + 1}")
        ns.add_row(
            # Make a fragment which is the length of the input Scaffold to
            # test trivial remapping
            Fragment(scf.name, 1, scf.length, 1)
        )

    ba1 = BuildAssembly(ia1.name, bp_per_texel=100)
    ba1.remap_to_input_assembly(p1, ia1)
    a2 = ba1.assemblies_with_scaffolds_fused()[None]
    assert str(ia1).replace("IndexedAssembly", "Assembly") == str(a2)


def test_cut_long_scaffold():
    start = Fragment("start", 1, 1_000_000_000, 1)
    middle = Fragment("middle", 1, 1000, 1)
    end = Fragment("end", 1, 999_999_000, 1)
    gap = Gap(10, "scaffold")

    ba = BuildAssembly("Scissors", max_contig_length=2_000_000_000)

    big_scffld_1 = Scaffold(
        "too_big",
        rows=[start, gap, middle, gap, end],
    )
    parts_1 = ba.cut_scaffold_if_too_long(big_scffld_1)
    assert repr(parts_1) == repr(
        [
            Scaffold(
                name="too_big_1",
                rows=[start],
            ),
            Scaffold(
                name="too_big_2",
                rows=[middle, gap, end],
            ),
        ]
    )

    big_scffld_2 = Scaffold(
        "too_big",
        rows=[end, gap, middle, gap, start],
    )
    parts_2 = ba.cut_scaffold_if_too_long(big_scffld_2)
    assert repr(parts_2) == repr(
        [
            Scaffold(
                name="too_big_1",
                rows=[end, gap, middle],
            ),
            Scaffold(
                name="too_big_2",
                rows=[start],
            ),
        ]
    )

    big_scffld_3 = Scaffold(
        "too_big",
        rows=[start, gap, start],
    )
    parts_3 = ba.cut_scaffold_if_too_long(big_scffld_3)
    assert repr(parts_3) == repr(
        [
            Scaffold(
                name="too_big_1",
                rows=[start],
            ),
            Scaffold(
                name="too_big_2",
                rows=[start],
            ),
        ]
    )

    big_scffld_4 = Scaffold(
        "too_big",
        rows=[end, gap, end],
    )
    parts_4 = ba.cut_scaffold_if_too_long(big_scffld_4)
    assert repr(parts_4) == repr(
        [
            Scaffold(
                name="too_big",
                rows=[end, gap, end],
            ),
        ]
    )

    big_scffld_5 = Scaffold(
        "too_big",
        rows=[end, middle, gap, middle, gap, end],
    )
    parts_5 = ba.cut_scaffold_if_too_long(big_scffld_5)
    assert repr(parts_5) == repr(
        [
            Scaffold(
                name="too_big_1",
                rows=[end, middle],
            ),
            Scaffold(
                name="too_big_2",
                rows=[middle, gap, end],
            ),
        ]
    )

    huge = Fragment("contiguous", 1, 2_000_000_001, 1)
    with pytest.raises(
        LongScaffoldCuttingError,
        match="Failed to find a gap before or after",
    ):
        ba.cut_scaffold_if_too_long(Scaffold("no_gaps", [huge]))

    with pytest.raises(
        LongScaffoldCuttingError,
        match="longer than max contig length",
    ):
        ba.cut_scaffold_if_too_long(Scaffold("gap but too ", [huge, gap, middle]))


def test_shuffled_assembly(seed="Shuffled assembly"):
    ia1 = make_random_assembly(
        seed=seed,
        scaffolds=5,
        rows=100,
        fragment_length=100_000,
    )

    m1, p2, m2 = shuffle_and_remap_assembly(ia1, f"seed={seed!r}")
    assert str(m1) == str(m2)


def shuffle_and_remap_assembly(asm, name):
    p1 = shuffled_assembly(asm, name)
    # print(p1)
    ba1 = BuildAssembly(name, default_gap=Gap(200, "scaffold"))
    ba1.remap_to_input_assembly(p1, asm)
    m1 = ba1.assemblies_with_scaffolds_fused()[None]
    # print(m1)

    p2 = fuzz_coordinates(p1)
    ba2 = BuildAssembly(name, default_gap=Gap(200, "scaffold"))
    ba2.remap_to_input_assembly(p2, asm)
    m2 = ba2.assemblies_with_scaffolds_fused()[None]
    print(ba2)
    ba2.log_multi_scaffolds()
    # print(m2)

    ### Are short contigs never lost?

    return m1, p2, m2


def shuffled_assembly(asm, name):
    # Make a list of baits from random length chunks of rows
    baits = []
    max_chunk = 10
    for scffld in asm.scaffolds:
        i = 0
        p = 0
        while i < len(scffld.rows):
            j = i + random.randint(1, max_chunk)  # noqa: S311
            chunk = scffld.rows[i:j]
            start = p + 1
            end = p + sum(f.length for f in chunk)
            strand = 1 if random.random() < 0.8 else -1  # noqa: S311
            if any(isinstance(x, Fragment) for x in chunk):
                # Don't add any baits which are only gaps
                baits.append(Fragment(scffld.name, start, end, strand))
            i = j
            p = end
    random.shuffle(baits)

    ptxt = Assembly(name)
    ptxt.bp_per_texel = 1
    sn = 0
    while len(baits):
        sn += 1
        scffld = Scaffold(f"RL_{sn}")
        ptxt.add_scaffold(scffld)
        j = random.randint(1, 4 * max_chunk)  # noqa: S311
        scffld.rows = baits[0:j]
        del baits[0:j]
    return ptxt


def fuzz_coordinates(asm):
    assembly_length = sum(x.length for x in asm.scaffolds)
    bp_per_texel = assembly_length / 2**15
    print(
        f"Assembly '{asm.name}' = {assembly_length} bp == {bp_per_texel} bp per texel"
    )
    new = Assembly(asm.name)
    new.bp_per_texel = bp_per_texel
    for scffld in asm.scaffolds:
        new_scffld = Scaffold(scffld.name)
        new.add_scaffold(new_scffld)
        for row in scffld.rows:
            if isinstance(row, Gap):
                new_scffld.add_row(row)
            else:
                start = row.start
                end = row.end

                # Simulate edit of long fragments when Pretext is zoomed out
                # screen_res = row.length / 1000
                # if screen_res > bp_per_texel:
                #     start = 1 + round_down(start, screen_res)
                #     end = round_down(end, screen_res)

                start = 1 + round_down(start, bp_per_texel)
                end = round_down(end, bp_per_texel)
                new_scffld.add_row(Fragment(row.name, start, end, row.strand))

    # if pairs := new.find_overlapping_fragments():
    #     msg = (
    #         f"Fuzzing ({bp_per_texel:.2f}) caused overlapping fragments:\n"
    #         + "\n".join(
    #             (
    #                 f"\nOverlap: {p[0][0].overlap_length(p[1][0])}\n"
    #                 f"{p[0][1].name} {p[0][0]}\n{p[1][1].name} {p[1][0]}"
    #             )
    #             for p in pairs
    #         )
    #     )
    #     raise Exception(msg)
    return new


def round_down(x, float_res):
    res = int(float_res)
    return int(res * math.floor(x / res))


def make_random_assembly(
    seed=None,
    scaffolds=3,
    rows=20,
    fragment_length=100_000,
):
    if seed:
        random.seed(seed)
    a1 = IndexedAssembly(seed if seed else "rand_asmbly")
    g1 = Gap(200, "scaffold")
    for sn in range(1, scaffolds + 1):
        s = Scaffold(f"scaffold_{sn}")
        s.rank = 3
        p = 0
        for rn in range(random.randint(1, rows)):  # noqa: S311
            # Don't add a Gap if it's the first row
            if rn != 0:
                s.add_row(g1)
                p += g1.length

            end_offset = max(999, random.randint(1, fragment_length))  # noqa: S311

            # Add a Fragment
            f = Fragment(
                s.name,
                p + 1,
                p + end_offset,
                1 if random.random() < 0.8 else -1,  # noqa: S311
            )
            s.add_row(f)
            p += f.length
        a1.add_scaffold(s)
    return a1


if __name__ == "__main__":
    test_shuffled_assembly()

    if True:
        # rndm = 8808117548524440979
        rndm = 7850684667499316592
        asm = make_random_assembly(
            seed=rndm,
            scaffolds=1,
            fragment_length=100_000_000,
        )
        m1, p2, m2 = shuffle_and_remap_assembly(asm, f"seed={rndm!r}")
        if str(m1) != str(m2):
            print(
                m1,
                "Fuzzed:",
                p2,
                "From fuzzed:",
                m2,
                sep="\n",
            )

    else:
        while True:
            random.seed()
            rndm = random.randint(0, sys.maxsize)  # noqa: S311
            asm = make_random_assembly(
                seed=rndm,
                scaffolds=1,
                fragment_length=100_000_000,
            )
            m1, p2, m2 = shuffle_and_remap_assembly(asm, f"seed={rndm!r}")
            if str(m1) != str(m2):
                print(
                    m1,
                    "Fuzzed:",
                    p2,
                    "From fuzzed:",
                    m2,
                    sep="\n",
                )
                break
