import pytest

from tola.assembly.fragment import Fragment
from tola.assembly.gap import Gap
from tola.assembly.overlap_result import OverlapResult

from .utils import strip_leading_spaces


def example_overlap_result():
    return OverlapResult(
        name="T1",
        bait=Fragment(name="chr_X", start=101_001, end=134_500, strand=-1),
        start=100_001,
        end=134_000,
        rows=[
            Fragment(name="frag_1", start=1, end=10_000, strand=1),
            Gap(length=200, gap_type="scaffold"),
            Fragment(name="frag_2", start=1, end=10_000, strand=1),
            Gap(length=200, gap_type="scaffold"),
            Fragment(name="frag_3", start=1, end=10_000, strand=1),
        ],
    )


def test_simple():
    o1 = OverlapResult(
        name="Simple example",
        bait=Fragment(name="small_bait", start=3, end=20, strand=-1),
        start=1,
        end=13,
        rows=[
            Fragment(name="tiny_a", start=1, end=5, strand=1),
            Gap(length=2, gap_type="scaffold"),
            Fragment(name="tiny_b", start=1, end=6, strand=1),
        ],
    )
    print(o1)
    assert o1.start == 1
    assert o1.end == 13
    assert o1.length == 13
    assert o1.start_overhang == 2
    assert o1.end_overhang == -7
    assert o1.start_row_bait_overlap == 3
    assert o1.end_row_bait_overlap == 6
    assert o1.length_error == -5
    assert o1.overhang_if_start_removed() == -5
    assert o1.overhang_if_end_removed() == -15
    assert o1.length_error_in_texels(10) == 5 / 10


def test_properties():
    o1 = example_overlap_result()
    print(o1)
    assert o1.start == 100_001
    assert o1.end == 134_000
    assert o1.length == 34_000
    assert o1.start_overhang == 1000
    assert o1.end_overhang == -500
    assert o1.start_row_bait_overlap == 9000
    assert o1.end_row_bait_overlap == 10_000
    assert o1.length_error == 500
    assert o1.overhang_if_start_removed() == -9200
    assert o1.overhang_if_end_removed() == -10700
    assert o1.length_error_in_texels(750) == 2 / 3


def test_manipulations():
    o1 = example_overlap_result()
    with pytest.raises(NotImplementedError):
        o1.reverse()
    with pytest.raises(NotImplementedError):
        o1.append_scaffold()

    o1.discard_start()
    assert len(o1.rows) == 3
    assert [x.name for x in o1.fragments()] == ["frag_2", "frag_3"]
    assert isinstance(o1.rows[1], Gap)
    assert o1.start_overhang == -9_200

    o2 = example_overlap_result()
    o2.discard_end()
    assert len(o2.rows) == 3
    assert [x.name for x in o2.fragments()] == ["frag_1", "frag_2"]
    assert isinstance(o2.rows[1], Gap)
    assert o2.end_overhang == -10_700


def test_trim_overhangs():
    o1 = example_overlap_result()

    # Change the start Fragment to one that overhangs by 71_000
    o1.rows[0] = Fragment(name="frag_x", start=1, end=80_000, strand=1)
    o1.start -= 70_000
    assert o1.start_overhang == 71_000
    assert o1.overhang_if_start_removed() == -9200
    assert o1.start_row_bait_overlap == 9000

    o1.trim_large_overhangs(20_000)
    assert len(o1.rows) == 3
    assert [x.name for x in o1.fragments()] == ["frag_2", "frag_3"]
    assert o1.start_overhang == -9_200

    # Check running again doesn't remove anything
    o1.trim_large_overhangs(20_000)
    assert len(o1.rows) == 3

    o2 = example_overlap_result()
    o2.discard_end()
    print(o2)
    assert o2.end_row_bait_overlap == 10_000


def test_to_scaffold():
    o1 = example_overlap_result()
    print(o1.to_scaffold())
    assert str(o1.to_scaffold()) == strip_leading_spaces(
        """
        T1
                  10_000  frag_3:1-10000(-)
                     200  Gap:200 scaffold
                  10_000  frag_2:1-10000(-)
                     200  Gap:200 scaffold
                  10_000  frag_1:1-10000(-)
        """
    )


def test_str_repr():
    o1 = OverlapResult(
        name="R4",
        bait=Fragment(name="scaffold_3", start=2463426, end=2540013, strand=-1),
        start=2463439,
        end=2540137,
        rows=[
            Fragment(name="scaffold_3", start=2463439, end=2490815, strand=1),
            Gap(length=200, gap_type="scaffold"),
            Fragment(name="scaffold_3", start=2491016, end=2528230, strand=1),
            Gap(length=200, gap_type="scaffold"),
            Fragment(name="scaffold_3", start=2528431, end=2540137, strand=1),
        ],
    )
    assert str(o1) == strip_leading_spaces(
        """
        R4
          length:         76_699
          bait:           76_588  scaffold_3:2463426-2540013(-)
          diff:              111
          overhang: -13
                 2_463_439      2_490_815      27_377  scaffold_3:2463439-2490815(+)
                 2_490_816      2_491_015         200  Gap:200 scaffold
                 2_491_016      2_528_230      37_215  scaffold_3:2491016-2528230(+)
                 2_528_231      2_528_430         200  Gap:200 scaffold
                 2_528_431      2_540_137      11_707  scaffold_3:2528431-2540137(+)
          overhang: 124"""
    )
    assert repr(o1) == strip_leading_spaces(
        """
        OverlapResult(
            name='R4',
            bait=Fragment(name='scaffold_3', start=2463426, end=2540013, strand=-1),
            start=2463439,
            end=2540137,
            rows=[
                Fragment(name='scaffold_3', start=2463439, end=2490815, strand=1),
                Gap(length=200, gap_type='scaffold'),
                Fragment(name='scaffold_3', start=2491016, end=2528230, strand=1),
                Gap(length=200, gap_type='scaffold'),
                Fragment(name='scaffold_3', start=2528431, end=2540137, strand=1),
            ],
        )"""
    )


if __name__ == "__main__":
    test_simple()
    test_properties()
    test_manipulations()
    test_trim_overhangs()
    test_to_scaffold()
    test_str_repr()
