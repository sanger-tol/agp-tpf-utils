import math
import random
import sys

from tola.assembly.assembly import Assembly
from tola.assembly.build_assembly import BuildAssembly
from tola.assembly.fragment import Fragment
from tola.assembly.gap import Gap
from tola.assembly.indexed_assembly import IndexedAssembly
from tola.assembly.scaffold import Scaffold


def test_no_coord_changes():
    ia1 = make_random_assembly(seed="Random assembly")

    p1 = Assembly("pretext")
    out_n = 0
    for scf in ia1.scaffolds:
        out_n += 1
        ns = Scaffold(f"New_{out_n}")
        ns.add_row(
            # Make a fragment which is the length of the input Scaffold to
            # test trivial remapping
            Fragment(scf.name, 1, scf.length, 1)
        )

    ba1 = BuildAssembly(ia1.name, bp_per_texel=100)
    ba1.remap_to_input_assembly(p1, ia1)
    a2 = ba1.assembly_with_scaffolds_fused()
    assert str(ia1).replace("IndexedAssembly", "Assembly") == str(a2)


def test_shuffled_assembly(seed="Shuffled assembly"):
    ia1 = make_random_assembly(
        seed=seed,
        scaffolds=5,
        rows=100,
        fragment_length=100_000,
    )

    m1, m2 = shuffle_and_remap_assembly(ia1, f"seed='{seed}'")
    assert str(m1) == str(m2)


def shuffle_and_remap_assembly(asm, name):
    p1 = shuffled_assembly(asm, name)
    # print(p1)
    ba1 = BuildAssembly(name, default_gap=Gap(200, "scaffold"))
    ba1.remap_to_input_assembly(p1, asm)
    m1 = ba1.assembly_with_scaffolds_fused()
    # print(m1)

    p2 = fuzz_coordinates(p1)
    ba2 = BuildAssembly(name, default_gap=Gap(200, "scaffold"))
    ba2.remap_to_input_assembly(p2, asm)
    m2 = ba2.assembly_with_scaffolds_fused()
    ba2.log_multi_scaffolds()
    # print(m2)

    ### Are short contigs never lost?

    return m1, m2


def shuffled_assembly(asm, name):
    # Make a list of baits from random length chunks of rows
    baits = []
    max_chunk = 10
    for scffld in asm.scaffolds:
        i = 0
        p = 0
        while i < len(scffld.rows):
            j = i + random.randint(1, max_chunk)
            chunk = scffld.rows[i:j]
            start = p + 1
            end = p + sum(f.length for f in chunk)
            strand = 1 if random.random() < 0.8 else -1
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
        j = random.randint(1, 4 * max_chunk)
        scffld.rows = baits[0:j]
        del baits[0:j]
    return ptxt


def fuzz_coordinates(asm):
    assembly_length = sum(x.length for x in asm.scaffolds)
    bp_per_texel = assembly_length / 2**15
    print(f"Assembly = {assembly_length} bp == {bp_per_texel} bp per texel")
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
                screen_res = row.length / 1000
                if screen_res > bp_per_texel:
                    start = 1 + round_down(start, screen_res)
                    end = round_down(end, screen_res)

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
        p = 0
        for rn in range(random.randint(1, rows)):
            # Don't add a Gap if it's the first row
            if rn != 0:
                s.add_row(g1)
                p += g1.length

            end_offset = max(999, random.randint(1, fragment_length))

            # Add a Fragment
            f = Fragment(
                s.name,
                p + 1,
                p + end_offset,
                1 if random.random() < 0.8 else -1,
            )
            s.add_row(f)
            p += f.length
        a1.add_scaffold(s)
    return a1


if __name__ == "__main__":
    # test_shuffled_assembly()
    while True:
        # Seed 8808117548524440979 fails
        random.seed()
        rndm = random.randint(0, sys.maxsize)
        asm = make_random_assembly(
            seed=rndm,
            scaffolds=1,
            fragment_length=100_000_000
        )
        m1, m2 = shuffle_and_remap_assembly(asm, f"seed='{rndm}'")
        if str(m1) != str(m2):
            print(m2)
            break
