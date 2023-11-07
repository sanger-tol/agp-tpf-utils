import pytest
import random

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
        x = i3.fetch_overlaps(Fragment("no data", 1, 100, 1))

    # Do not set field directly in IndexedAssembly object
    i3._scaffold_dict["scaffold_1"] = s1
    with pytest.raises(ValueError, match=r"Scaffold 'scaffold_1' is not indexed"):
        x = i3.fetch_overlaps(Fragment("scaffold_1", 1, 100, 1))


def test_fetch_overlapping():
    """
    IndexedAssembly: rand_asmbly

      scaffold_1
                  1        9613  fragment_0:1-9613(-)
               9614       15573  fragment_1:1-5960(+)
              15574       15773  Gap:200 scaffold
              15774       19116  fragment_2:1-3343(-)
              19117       20030  fragment_3:1-914(-)
              20031       28894  fragment_4:1-8864(+)
              28895       30532  fragment_5:1-1638(+)
              30533       33334  fragment_6:1-2802(-)

      scaffold_2
                  1         547  fragment_0:1-547(+)
                548       10328  fragment_1:1-9781(-)
              10329       10528  Gap:200 scaffold
              10529       20000  fragment_2:1-9472(+)
              20001       20200  Gap:200 scaffold
              20201       27594  fragment_3:1-7394(+)
              27595       27794  Gap:200 scaffold
              27795       33962  fragment_4:1-6168(+)
              33963       37543  fragment_5:1-3581(+)
              37544       42112  fragment_6:1-4569(+)
              42113       50048  fragment_7:1-7936(+)
              50049       54218  fragment_8:1-4170(-)
              54219       54418  Gap:200 scaffold
              54419       57333  fragment_9:1-2915(+)
              57334       65403  fragment_10:1-8070(-)
              65404       68622  fragment_11:1-3219(+)
              68623       76647  fragment_12:1-8025(+)

      scaffold_3
                  1        9394  fragment_0:1-9394(+)
               9395       12385  fragment_1:1-2991(+)
              12386       16518  fragment_2:1-4133(+)
              16519       25077  fragment_3:1-8559(+)
              25078       25277  Gap:200 scaffold
              25278       31645  fragment_4:1-6368(+)
              31646       31845  Gap:200 scaffold
              31846       34376  fragment_5:1-2531(+)
    """
    asm = example_assembly()
    bait_list = (
        # First base of Scaffold:
        Fragment("scaffold_1", 1, 1, 1),
        # Last base of Scaffold:
        Fragment("scaffold_1", 33334, 33334, 1),
        # One base beyond last base of Scaffold:
        Fragment("scaffold_1", 33335, 33335, 1),
        Fragment("scaffold_2", 10328, 65404, 1),
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
        if not overlaps:
            assert None == asm.fetch_overlaps(bait)
        else:
            found = asm.fetch_overlaps(bait)
            print(Scaffold("found", found))
            print(Scaffold("correct answer", overlaps))
            assert overlaps == asm.fetch_overlaps(bait)

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
    assert asm2.fetch_overlaps(bait_list[0]) == [bait_list[0]]


def example_assembly():
    """
    Made from make_random_assembly(seed="Random assembly")
    """
    return IndexedAssembly(
        name="rand_assmbly",
        scaffolds=[
            Scaffold(
                name="scaffold_1",
                rows=[
                    Fragment(name="fragment_0", start=1, end=9613, strand=-1),
                    Fragment(name="fragment_1", start=1, end=5960, strand=1),
                    Gap(length=200, gap_type="scaffold"),
                    Fragment(name="fragment_2", start=1, end=3343, strand=-1),
                    Fragment(name="fragment_3", start=1, end=914, strand=-1),
                    Fragment(name="fragment_4", start=1, end=8864, strand=1),
                    Fragment(name="fragment_5", start=1, end=1638, strand=1),
                    Fragment(name="fragment_6", start=1, end=2802, strand=-1),
                ],
            ),
            Scaffold(
                name="scaffold_2",
                rows=[
                    Fragment(name="fragment_0", start=1, end=547, strand=1),
                    Fragment(name="fragment_1", start=1, end=9781, strand=-1),
                    Gap(length=200, gap_type="scaffold"),
                    Fragment(name="fragment_2", start=1, end=9472, strand=1),
                    Gap(length=200, gap_type="scaffold"),
                    Fragment(name="fragment_3", start=1, end=7394, strand=1),
                    Gap(length=200, gap_type="scaffold"),
                    Fragment(name="fragment_4", start=1, end=6168, strand=1),
                    Fragment(name="fragment_5", start=1, end=3581, strand=1),
                    Fragment(name="fragment_6", start=1, end=4569, strand=1),
                    Fragment(name="fragment_7", start=1, end=7936, strand=1),
                    Fragment(name="fragment_8", start=1, end=4170, strand=-1),
                    Gap(length=200, gap_type="scaffold"),
                    Fragment(name="fragment_9", start=1, end=2915, strand=1),
                    Fragment(name="fragment_10", start=1, end=8070, strand=-1),
                    Fragment(name="fragment_11", start=1, end=3219, strand=1),
                    Fragment(name="fragment_12", start=1, end=8025, strand=1),
                ],
            ),
            Scaffold(
                name="scaffold_3",
                rows=[
                    Fragment(name="fragment_0", start=1, end=9394, strand=1),
                    Fragment(name="fragment_1", start=1, end=2991, strand=1),
                    Fragment(name="fragment_2", start=1, end=4133, strand=1),
                    Fragment(name="fragment_3", start=1, end=8559, strand=1),
                    Gap(length=200, gap_type="scaffold"),
                    Fragment(name="fragment_4", start=1, end=6368, strand=1),
                    Gap(length=200, gap_type="scaffold"),
                    Fragment(name="fragment_5", start=1, end=2531, strand=1),
                ],
            ),
        ],
    )


def make_random_assembly(seed=None, scaffolds=3):
    if seed:
        random.seed(seed)
    a1 = IndexedAssembly("rand_asmbly")
    g1 = Gap(200, "scaffold")
    for sn in range(1, scaffolds + 1):
        s = Scaffold(f"scaffold_{sn}")
        for rn in range(random.randint(1, 20)):
            # Maybe add a Gap, but not on the first row
            if rn != 0 and random.random() < 0.4:
                s.add_row(g1)

            # Add a Fragment
            s.add_row(
                Fragment(
                    f"fragment_{rn}",
                    1,
                    random.randint(1, 10_000),
                    1 if random.random() < 0.8 else -1,
                )
            )
        a1.add_scaffold(s)
    return a1


if __name__ == "__main__":
    asm1 = make_random_assembly(seed="Random assembly")
    print(asm1)
