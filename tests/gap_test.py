import pytest
from tola.assembly.gap import Gap


def test_create():
    x = Gap(200, "scaffold")
    assert isinstance(x, Gap)
    assert x.length == 200

    # Gap creation is @cache'd
    y = Gap(200, "scaffold")
    assert x is y

    z = Gap(200, "contig")
    assert x is not z

def test_bad_attributes():
    with pytest.raises(ValueError):
        g1 = Gap("x", "scaffold")


def test_stringify():
    g1 = Gap(100, "scaffold")
    g2 = Gap(200, "scaffold")

    assert str(g1) == "Gap:100 scaffold"
    assert str(g2) == "Gap:200 scaffold"

    assert repr(g1) == "Gap(length=100, gap_type='scaffold')"
    assert repr(g2) == "Gap(length=200, gap_type='scaffold')"
