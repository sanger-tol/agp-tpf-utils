from tola.assembly.assembly import Assembly
from tola.assembly.fragment import Fragment
from tola.assembly.gap import Gap
from tola.assembly.scaffold import Scaffold


def test_find_overlaps():
    x1 = Fragment(name="F1", start=101, end=120, strand=1)
    x2 = Fragment(name="F1", start=121, end=130, strand=1)
    x3 = Fragment(name="F1", start=1, end=101, strand=1)
    x4 = Fragment(name="F1", start=141, end=150, strand=1)
    x5 = Fragment(name="F1", start=130, end=140, strand=-1)

    asm = Assembly(
        name="has_overlaps",
        scaffolds=[
            Scaffold(
                name="S1",
                rows=[x1, x2, x3],
            ),
            Scaffold(
                name="S2",
                rows=[x4, x5],
            ),
        ],
    )

    over_pairs = asm.find_overlapping_fragments()

    assert len(over_pairs) == 2

    m0, m1 = over_pairs[0]
    assert m0[0] is x1
    assert m1[0] is x3

    m2, m3 = over_pairs[1]
    assert m2[0] is x2
    assert m3[0] is x5


if __name__ == "__main__":
    test_find_overlaps()
