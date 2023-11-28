from tola.assembly.scaffold import Fragment, Gap, Scaffold


def test_create():
    s1 = Scaffold(
        name="scaffold_1",
        rows=[
            Fragment("scaffold_12", 1, 20_000, 1),
            Gap(200, "Type-2"),
            Fragment("scaffold_12", 23_201, 140_112, -1),
        ],
    )
    assert isinstance(s1, Scaffold)
    assert s1.length == 137_112
    assert s1.fragments_length == 136_912
    assert s1.gaps_length == 200


def test_iterators():
    f1 = Fragment("scaffold_12", 1, 20_000, 1)
    f2 = Fragment("scaffold_12", 23_200, 140_112, -1)
    f3 = Fragment("scaffold_12", 140_113, 244_491, 1)
    f4 = Fragment("scaffold_3", 1, 244_232, 1)

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

    assert s1.last_row_is_fragment == True
    s1.add_row(g1)
    assert s1.last_row_is_fragment == False


def test_reverse():
    f1 = Fragment("scaffold_12", 1, 20_000, 1)
    f2 = Fragment("scaffold_12", 23_200, 140_112, -1)
    f3 = Fragment("scaffold_12", 140_113, 244_491, 1)
    g1 = Gap(100, "Type-2")
    s1 = Scaffold(name="rev. test", rows=[f1, g1, f2, f3])
    s1r = s1.reverse()
    assert [x.strand for x in s1r.fragments()] == [-1, 1, -1]
    assert s1r.rows[2] is g1
