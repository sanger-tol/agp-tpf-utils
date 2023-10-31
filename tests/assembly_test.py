from tola.assembly.assembly import Assembly, Scaffold


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
    a1.sort_scaffolds_by_name()
    assert list(s.name for s in a1.scaffolds) == [
        "chr1",
        "chr2",
        "chr10",
        "chr22",
    ]
