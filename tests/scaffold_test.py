import sys

from tola.assembly.scaffold import Scaffold, Fragment, Gap


def test_create():
    s1 = Scaffold(
        name="scaffold_1",
        rows=[
            Fragment("scaffold_12", 1, 20_000, 1),
            Gap(200, "Type-2"),
            Fragment("scaffold_12", 23_200, 140_112, -1),
        ],
    )
    assert isinstance(s1, Scaffold)


def test_iterators():
    f1 = Fragment("scaffold_12", 1, 20_000, 1)
    f2 = Fragment("scaffold_12", 23_200, 140_112, -1)
    f3 = Fragment("scaffold_12", 140_113, 244_491, -1)
    f4 = Fragment("scaffold_3", 1, 244_232, -1)

    g1 = Gap(100, "Type-2")
    g2 = Gap(200, "Type-2")

    s1 = Scaffold(
        name="scaffold_1",
        rows=[
            f1,  # 0
            g1,  # 1
            f2,  # 2
            g2,  # 3
            f3,  # 4
            f4,  # 5
        ],
    )

    assert list(s1.fragments()) == [f1, f2, f3, f4]
    assert list(s1.idx_fragments()) == [
        (0, f1),
        (2, f2),
        (4, f3),
        (5, f4),
    ]
    assert list(s1.gaps()) == [g1, g2]
    assert list(s1.idx_gaps()) == [
        (1, g1),
        (3, g2),
    ]
