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
    # print(ia1)

    # Make a list of baits from random length chunks of rows
    baits = []
    max_chunk = 10
    for scffld in ia1.scaffolds:
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

    p1 = Assembly(f"Edited {seed}")
    p1.bp_per_texel = 1
    sn = 0
    while len(baits):
        sn += 1
        scffld = Scaffold(f"RL_{sn}")
        p1.add_scaffold(scffld)
        j = random.randint(1, 4 * max_chunk)
        scffld.rows = baits[0:j]
        del baits[0:j]

    ba1 = BuildAssembly(f"Build {seed}", default_gap=Gap(200, "scaffold"))
    ba1.remap_to_input_assembly(p1, ia1)
    assert len(ba1.problem_scaffolds) == 0
    m1 = ba1.assembly_with_scaffolds_fused()
    # print(m1)

    p2 = fuzz_coordinates(p1, 1000)
    ba2 = BuildAssembly(ba1.name, default_gap=Gap(200, "scaffold"))
    ba2.remap_to_input_assembly(p2, ia1)
    m2 = ba2.assembly_with_scaffolds_fused()
    ba2.log_problem_scaffolds()
    if pairs := m2.find_overlapping_fragments():
        report_overlaps(pairs)
    # print(m2)

    assert str(m1) == str(m2)


def report_overlaps(pairs):
    for pr in pairs:
        f1, s1 = pr[0]
        f2, s2 = pr[1]
        print(f"\nDuplicated component:\n{s1.name} {f1.length:6d} {f1}\n{s2.name} {f2.length:6d} {f2}")


def fuzz_coordinates(asm, bp_per_texel):
    new = Assembly(asm.name)
    new.bp_per_texel = bp_per_texel
    fuzz = bp_per_texel
    for scffld in asm.scaffolds:
        new_scffld = Scaffold(scffld.name)
        new.add_scaffold(new_scffld)
        prev_fuzz = random.randint(-fuzz, fuzz)
        for row in scffld.rows:
            if isinstance(row, Gap):
                new_scffld.add_row(row)
            else:
                new_fuzz = random.randint(-fuzz, fuzz)
                start = row.start + prev_fuzz
                if start < 1:
                    start = 1
                end = row.end + new_fuzz
                if end < 1:
                    end = 1
                if start > end:
                    start, end = end, start
                new_scffld.add_row(Fragment(row.name, start, end, row.strand))
                prev_fuzz = new_fuzz
    return new


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

            # Add a Fragment
            f = Fragment(
                s.name,
                p + 1,
                p + random.randint(fragment_length // 50, fragment_length),
                1 if random.random() < 0.8 else -1,
            )
            s.add_row(f)
            p += f.length
        a1.add_scaffold(s)
    return a1


if __name__ == "__main__":
    # test_shuffled_assembly()
    random.seed()
    rndm = random.randint(0, sys.maxsize)
    test_shuffled_assembly(rndm)
