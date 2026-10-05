"""The public surface of the 1b build: the graph on units, Eligibility, the conversions, the
route structures, the selection and placement axes and the ledger."""

import _netgraph_core

import netgraph_core as ngc

PUBLIC_1B = {
    "StrictMultiDiGraph",
    "Eligibility",
    "to_units",
    "as_float",
    "PredDAG",
    "Path",
    "EdgeTieBreak",
    "EdgePriority",
    "PathSelection",
    "FlowPlacement",
    "Admission",
    "FlowGraph",
    "FlowRecord",
}

# Private recorder entries: the flow core (1a), the single-bundle Pinned transaction, the
# internal search and the apportionment (1b). Never re-exported by netgraph_core.
PRIVATE_1B = {
    "_place_global_units",
    "_PinnedFlow",
    "_spf",
    "_round_robin",
    "_highest_averages",
}


def test_public_surface():
    assert set(ngc.__all__) == PUBLIC_1B | {"__version__"}
    for name in ngc.__all__:
        assert hasattr(ngc, name)
    assert isinstance(ngc.__version__, str)


def test_extension_surface_is_exactly_the_1b_set():
    names = {n for n in dir(_netgraph_core) if not n.startswith("__")}
    assert names == PUBLIC_1B | PRIVATE_1B
    for name in PRIVATE_1B:
        assert not hasattr(ngc, name)
