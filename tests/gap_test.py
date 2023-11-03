import pytest
from tola.assembly.gap import Gap


def test_create():
    x = Gap(200, "Type-2")
    assert isinstance(x, Gap)
    assert x.length == 200


def test_bad_attributes():
    with pytest.raises(ValueError):
        g1 = Gap("x", "Type-2")

    with pytest.raises(ValueError):
        g1 = Gap("x", "Type-2")


def test_stringify():
    g1 = Gap(100, "Type-2")
    g2 = Gap(200, "Type-2", "proximity_ligation")

    assert str(g1) == "Gap:100 Type-2"
    assert str(g2) == "Gap:200 Type-2 proximity_ligation"
