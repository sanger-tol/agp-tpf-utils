from tola.assembly.assembly import Assembly
from tola.assembly.fragment import Fragment
from tola.assembly.gap import Gap
from tola.assembly.scaffold import Scaffold

from .utils import strip_leading_spaces


def test_natural_sort():
    s1 = Scaffold(name="chr12.34x")
    assert Assembly.name_natural_key(s1) == ("chr", 12, ".", 34, "x")

    a1 = Assembly(
        name="test",
        scaffolds=[
            Scaffold(name="chr2"),
            Scaffold(name="chr1"),
            Scaffold(name="chr22"),
            Scaffold(name="chr10"),
        ],
    )
    scfflds_sorted = a1.scaffolds_sorted_by_name()
    assert list(s.name for s in scfflds_sorted) == [
        "chr1",
        "chr2",
        "chr10",
        "chr22",
    ]


def test_str_and_repr():
    """
    Tests the implementation of __str__ and __repr__ in Assembly. Because it
    calls those methods on Scaffold, which calls them on Fragment and Gap, it
    tests __str__ and __repr__ implementations all the way down the tree.
    Changes in __str__ and __repr__ in Scaffold, Fragment or Gap will need
    changes here too.
    """
    f1 = Fragment("scaffold_12", 1, 20_000, 1)
    f2 = Fragment("scaffold_12", 23_200, 140_112, -1)
    f3 = Fragment("scaffold_12", 140_113, 244_491, 1, ("Painted",))
    f4 = Fragment("scaffold_3", 1, 244_232, 1)
    f5 = Fragment("scaffold_7", 1, 11_033_114_755, 1, ("Painted", "X"))
    f6 = Fragment("scaffold_7", 11_049_229_141, 11_049_229_150, 1, ("X", "Painted"))
    g1 = Gap(100, "Type-2")
    g2 = Gap(200, "Type-2")

    s1 = Scaffold(
        name="chr1",
        rows=[f1, g1, f2, g2, f3, f4],
    )

    s2 = Scaffold(
        name="chrX",
        rows=[f5, f6],
    )

    # Make an empty Assembly and test str() and repr()
    a1 = Assembly(name="hap1")
    assert str(a1) == "Assembly: hap1\n"
    assert repr(a1) == strip_leading_spaces(
        """
        Assembly(
            name='hap1',
        )""",
    )

    # Add header lines and test
    a1.add_header_line("DESCRIPTION: Generated by PretextView Version 0.2.5")
    a1.add_header_line("HiC MAP RESOLUTION: 8666.611572 bp/texel")
    assert str(a1) == strip_leading_spaces(
        """
        Assembly: hap1
          # DESCRIPTION: Generated by PretextView Version 0.2.5
          # HiC MAP RESOLUTION: 8666.611572 bp/texel
        """,
    )
    assert repr(a1) == strip_leading_spaces(
        """
        Assembly(
            name='hap1',
            header=[
                'DESCRIPTION: Generated by PretextView Version 0.2.5',
                'HiC MAP RESOLUTION: 8666.611572 bp/texel',
            ],
        )""",
    )
    assert a1.bp_per_texel == 8666.611572
    a1.bp_per_texel = 2000
    assert a1.bp_per_texel == 2000

    # Add Scaffolds and test str() and repr() again
    a1.add_scaffold(s1)
    a1.add_scaffold(s2)
    assert str(a1) == strip_leading_spaces(
        """
        Assembly: hap1
          # DESCRIPTION: Generated by PretextView Version 0.2.5
          # HiC MAP RESOLUTION: 8666.611572 bp/texel

          chr1
                      1       20000  scaffold_12:1-20000(+)
                  20001       20100  Gap:100 Type-2
                  20101      137013  scaffold_12:23200-140112(-)
                 137014      137213  Gap:200 Type-2
                 137214      241592  scaffold_12:140113-244491(+) Painted
                 241593      485824  scaffold_3:1-244232(+)

          chrX
                      1 11033114755  scaffold_7:1-11033114755(+) Painted X
            11033114756 11033114765  scaffold_7:11049229141-11049229150(+) X Painted
        """,
    )
    assert repr(a1) == strip_leading_spaces(
        """
        Assembly(
            name='hap1',
            header=[
                'DESCRIPTION: Generated by PretextView Version 0.2.5',
                'HiC MAP RESOLUTION: 8666.611572 bp/texel',
            ],
            scaffolds=[
                Scaffold(
                    name='chr1',
                    rows=[
                        Fragment(name='scaffold_12', start=1, end=20000, strand=1),
                        Gap(length=100, gap_type='Type-2'),
                        Fragment(name='scaffold_12', start=23200, end=140112, strand=-1),
                        Gap(length=200, gap_type='Type-2'),
                        Fragment(name='scaffold_12', start=140113, end=244491, strand=1, tags=('Painted',)),
                        Fragment(name='scaffold_3', start=1, end=244232, strand=1),
                    ],
                ),
                Scaffold(
                    name='chrX',
                    rows=[
                        Fragment(name='scaffold_7', start=1, end=11033114755, strand=1, tags=('Painted', 'X')),
                        Fragment(name='scaffold_7', start=11049229141, end=11049229150, strand=1, tags=('X', 'Painted')),
                    ],
                ),
            ],
        )""",
    )


if __name__ == "__main__":
    test_natural_sort()
    test_str()
