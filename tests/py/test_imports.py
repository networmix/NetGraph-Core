"""The public surface of the 1a build: the graph on units, Eligibility and the conversions."""

import _netgraph_core

import netgraph_core as ngc


def test_public_surface():
    assert set(ngc.__all__) == {
        "__version__",
        "StrictMultiDiGraph",
        "Eligibility",
        "to_units",
        "as_float",
    }
    for name in ngc.__all__:
        assert hasattr(ngc, name)
    assert isinstance(ngc.__version__, str)


def test_extension_surface_is_exactly_the_1a_set():
    # The old library is deleted in 1a: the extension exports the graph on units, the
    # conversions and the private flow-core entry, nothing else.
    public = {n for n in dir(_netgraph_core) if not n.startswith("__")}
    assert public == {
        "StrictMultiDiGraph",
        "Eligibility",
        "to_units",
        "as_float",
        "_place_global_units",
    }
    assert not hasattr(ngc, "_place_global_units")
