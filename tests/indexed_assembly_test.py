import random

import pytest
from tola.assembly.assembly import Assembly
from tola.assembly.fragment import Fragment
from tola.assembly.gap import Gap
from tola.assembly.indexed_assembly import IndexedAssembly
from tola.assembly.scaffold import Scaffold


def test_create():
    g1 = Gap(200, "scaffold")

    s1 = Scaffold(
        "scaffold_1",
        rows=[
            Fragment("cmpnt_1", 1, 100, 1),
            Fragment("cmpnt_1", 101, 200, 1),
            g1,
            Fragment("cmpnt_1", 201, 300, -1),
        ],
    )

    s2 = Scaffold(
        "scaffold_2",
        rows=[
            Fragment("cmpnt_2", 1, 100, 1),
            Fragment("cmpnt_2", 101, 200, 1),
        ],
    )

    i1 = IndexedAssembly("asm_1", scaffolds=[s1, s2])
    a1 = Assembly("asm_1", scaffolds=[s1, s2])
    i2 = IndexedAssembly.new_from_assembly(a1)

    assert str(i1) == str(i2)
    assert repr(i1) == repr(i2)

    with pytest.raises(ValueError, match=r"Already have Scaffold named 'scaffold_1'"):
        i1.add_scaffold(s1)

    assert i1.scaffold_by_name("scaffold_1") is s1
    with pytest.raises(ValueError, match=r"No such Scaffold"):
        x = i1.scaffold_by_name("nonesuch")

    i3 = IndexedAssembly("empty scaffold", scaffolds=[Scaffold("no data")])
    with pytest.raises(ValueError, match=r"Scaffold 'no data' is empty"):
        x = i3.find_overlaps(Fragment("no data", 1, 100, 1))

    # Do not set field directly in IndexedAssembly object
    i3._scaffold_dict["scaffold_1"] = s1
    with pytest.raises(ValueError, match=r"Scaffold 'scaffold_1' is not indexed"):
        x = i3.find_overlaps(Fragment("scaffold_1", 1, 100, 1))


def test_find_overlapping():
    """
    IndexedAssembly: Random assembly

      scaffold_1
                  1        9613  s1_fragment_0:1-9613(-)
               9614       15573  s1_fragment_1:1-5960(+)
              15574       15773  Gap:200 scaffold
              15774       19116  s1_fragment_2:1-3343(-)
              19117       20030  s1_fragment_3:1-914(-)
              20031       28894  s1_fragment_4:1-8864(+)
              28895       30532  s1_fragment_5:1-1638(+)
              30533       33334  s1_fragment_6:1-2802(-)

      scaffold_2
                  1         547  s2_fragment_0:1-547(+)
                548       10328  s2_fragment_1:1-9781(-)
              10329       10528  Gap:200 scaffold
              10529       20000  s2_fragment_2:1-9472(+)
              20001       20200  Gap:200 scaffold
              20201       27594  s2_fragment_3:1-7394(+)
              27595       27794  Gap:200 scaffold
              27795       33962  s2_fragment_4:1-6168(+)
              33963       37543  s2_fragment_5:1-3581(+)
              37544       42112  s2_fragment_6:1-4569(+)
              42113       50048  s2_fragment_7:1-7936(+)
              50049       54218  s2_fragment_8:1-4170(-)
              54219       54418  Gap:200 scaffold
              54419       57333  s2_fragment_9:1-2915(+)
              57334       65403  s2_fragment_10:1-8070(-)
              65404       68622  s2_fragment_11:1-3219(+)
              68623       76647  s2_fragment_12:1-8025(+)

      scaffold_3
                  1        9394  s3_fragment_0:1-9394(+)
               9395       12385  s3_fragment_1:1-2991(+)
              12386       16518  s3_fragment_2:1-4133(+)
              16519       25077  s3_fragment_3:1-8559(+)
              25078       25277  Gap:200 scaffold
              25278       31645  s3_fragment_4:1-6368(+)
              31646       31845  Gap:200 scaffold
              31846       34376  s3_fragment_5:1-2531(+)
    """
    asm = make_random_assembly(seed="Random assembly")
    bait_list = (
        Fragment("scaffold_1", 1, 1, 1),  # First base of Scaffold
        Fragment("scaffold_1", 33334, 33334, 1),  # Last base of Scaffold
        Fragment(
            "scaffold_1", 33335, 33335, 1
        ),  # One base beyond last base of Scaffold
        Fragment("scaffold_2", 27595, 54418, 1),  # Gaps on either end of overlap
        Fragment("scaffold_2", 10328, 65404, 1),
        Fragment("scaffold_3", 25078, 25277, 1),  # Only one gap
        Fragment("scaffold_3", 34376, 34376, 1),
    )

    for bait in bait_list:
        s = asm.scaffold_by_name(bait.name)
        overlaps = []
        offset = 0
        for row in s.rows:
            asm_row = Fragment(bait.name, offset + 1, offset + row.length, 1)
            if asm_row.overlaps(bait):
                overlaps.append(row)
            offset += row.length

        # Remove leading and trailing Gaps
        while overlaps and isinstance(overlaps[0], Gap):
            overlaps.pop(0)
        while overlaps and isinstance(overlaps[-1], Gap):
            overlaps.pop(-1)

        if not overlaps:
            assert asm.find_overlaps(bait) is None
        else:
            found = asm.find_overlaps(bait)
            assert found is not None
            assert found.bait is bait
            print(found)
            print(Scaffold("correct answer", overlaps))
            assert overlaps == found.rows

    asm2 = IndexedAssembly(
        "single entry scaffold",
        scaffolds=[
            Scaffold(
                "scaffold_1",
                [bait_list[0]],
            ),
        ],
    )
    print(asm2)
    assert asm2.find_overlaps(bait_list[0]).rows == [bait_list[0]]


def make_random_assembly(
    seed=None,
    scaffolds=3,
    rows=20,
    fragment_length=10_000,
    gap_threshold=0.4,
):
    if seed:
        random.seed(seed)
    a1 = IndexedAssembly(seed if seed else "rand_asmbly")
    g1 = Gap(200, "scaffold")
    for sn in range(1, scaffolds + 1):
        s = Scaffold(f"scaffold_{sn}")
        for rn in range(random.randint(1, rows)):
            # Maybe add a Gap, but not on the first row
            if rn != 0 and random.random() <= gap_threshold:
                s.add_row(g1)

            # Add a Fragment
            s.add_row(
                Fragment(
                    f"s{sn}_fragment_{rn}",
                    1,
                    random.randint(1, fragment_length),
                    1 if random.random() < 0.8 else -1,
                )
            )
        a1.add_scaffold(s)
    return a1


if __name__ == "__main__":
    asm1 = make_random_assembly(seed="Random assembly")
    print(asm1)
