import pytest

from tola.assembly.fragment import Fragment


def test_create():
    x = Fragment("chr1", 1, 20_000, -1, "Painted")
    assert isinstance(x, Fragment)


def test_bad_attributes():
    with pytest.raises(ValueError):
        f1 = Fragment("chr1", "x", 20_000, 1)

    with pytest.raises(ValueError):
        f2 = Fragment("chr1", 1, "20000.2", 1)

    with pytest.raises(ValueError):
        f3 = Fragment("chr1", 1, 20_000, "+")

    with pytest.raises(ValueError):
        f4 = Fragment("chr1", 1, 20_000, 2)

    with pytest.raises(ValueError):
        f5 = Fragment("chr1", 101, 100, 1)


def test_length():
    f1 = Fragment("chr1", 10, 10, 1)
    f2 = Fragment("chr1", 101, 200, -1)

    assert f1.length == 1
    assert f2.length == 100


def test_equals():
    f1 = Fragment("chr1", 1, 20_000, 1)
    f2 = Fragment("chr1", "1", "20000", 1)
    f3 = Fragment("chr1", 1, 20_000, 1, "Painted")
    f4 = Fragment(
        strand=1,
        start=1,
        end=20_000,
        name="chr1",
    )
    f5 = Fragment("chr2", 1, 20_000, 1)

    assert f1 == f1
    assert f1 == f2
    assert f1 != f3
    assert f1 == f4
    assert f1 != f5


def test_overlaps():
    f1 = Fragment("chr1", 1, 100, 1)
    f2 = Fragment("chr1", 100, 120, 1)
    f3 = Fragment("chr1", 121, 140, 1)
    f4 = Fragment("chr1", 140, 160, -1)

    assert f1.overlaps(f2)
    assert f2.overlaps(f1)
    assert not f2.overlaps(f3)
    assert f3.overlaps(f4)


def test_overlap_length():
    f1 = Fragment("chr20", 1, 102, 1)
    f2 = Fragment("chr20", 100, 120, 1)
    f3 = Fragment("chr20", 121, 140, 1)

    assert f1.overlap_length(f2) == 3
    assert f2.overlap_length(f1) == 3
    assert f2.overlap_length(f3) is None


def test_gap_between():
    f1 = Fragment("chr20", 1, 102, 1)
    f2 = Fragment("chr20", 100, 120, 1)
    f3 = Fragment("chr20", 123, 140, 1)
    f4 = Fragment("chr20", 141, 150, -1)

    assert f1.gap_between(f2) is None
    assert f2.gap_between(f1) is None
    assert f2.gap_between(f3) == 2
    assert f3.gap_between(f4) == 0


def test_reverse():
    f1 = Fragment("chr20", 1, 102, 1, ("Painted", "Haplotig"))
    assert f1.reverse() == Fragment("chr20", 1, 102, -1, ("Painted", "Haplotig"))


def test_abuts():
    f1 = Fragment("chr1", 1, 100, 1)
    f2 = Fragment("chr1", 101, 120, 1)
    f3 = Fragment("chr1", 102, 120, 1)
    f4 = Fragment("chr1", 100, 120, 1)

    assert f1.abuts(f2)
    assert f2.abuts(f1)
    assert not f1.abuts(f3)
    assert not f1.abuts(f4)


def test_stringify():
    f1 = Fragment("chr1", 1, 20_000, 1)
    f2 = Fragment("chr1", 1, 20_000, -1)
    f3 = Fragment("chr1", 1, 20_000, 0)
    f4 = Fragment("chr1", 1, 20_000, 1, ("Painted",))
    f5 = Fragment("chr1", 1, 20_000, 1, ("Painted", "X"))

    assert str(f1) == "chr1:1-20000(+)"
    assert str(f2) == "chr1:1-20000(-)"
    assert str(f3) == "chr1:1-20000(.)"
    assert str(f4) == "chr1:1-20000(+) Painted"
    assert str(f5) == "chr1:1-20000(+) Painted X"


def test_tuples():
    f1 = Fragment("chr1", 1, 20_000, 1)
    f2 = Fragment("chr2", 2, 30_000, -1)
    f3 = Fragment("chr3", 3, 40_000, 1)
    f4 = Fragment("chr4", 4, 50_000, -1)

    assert f1.key_tuple == ("chr1", 1, 20_000)
    assert f2.key_tuple == ("chr2", 2, 30_000)
    assert f1.junction_tuple(f3) == ("chr1", 20_000, "chr3", 3)
    assert f1.junction_tuple(f2) == ("chr1", 20_000, 30_000, "chr2")
    assert f2.junction_tuple(f1) == (2, "chr2", "chr1", 1)
    assert f2.junction_tuple(f4) == ("chr4", 50_000, "chr2", 2)
