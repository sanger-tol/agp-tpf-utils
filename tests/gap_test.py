import pytest
from tola.assembly.gap import Gap


def test_create():
    x = Gap(200, "Type-2")
    assert isinstance(x, Gap)
    assert x.length == 200


def test_bad_attributes():
    with pytest.raises(ValueError):
        g1 = Gap("x", "scaffold")


def test_stringify():
    g1 = Gap(100, "scaffold")
    g2 = Gap(200, "scaffold")

    assert str(g1) == "Gap:100 scaffold"
    assert str(g2) == "Gap:200 scaffold"
