"""The public surface: the graph on units, Eligibility, the conversions, the route structures,
the selection and placement axes, the ledger, the policy engine, the DagCache, the evaluators
and the public search; the extension's private entry, which netgraph_core never re-exports;
the retired spellings and the retired search entry points, which exist nowhere, without
aliases; and the argument checks of the bindings, each with its exception class and whole
message."""

from types import SimpleNamespace

import _netgraph_core
import numpy as np
import pytest
from helpers import RES, config, exact, sel
from private_hooks import HOOKS

import netgraph_core as ngc

PUBLIC = {
    "Graph",
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
    "FlowLedger",
    "FlowRecord",
    "Absolute",
    "Stretch",
    "DiminishingReturns",
    "CostOnly",
    "CapacityAware",
    "Global",
    "Pinned",
    "FlowPolicyConfig",
    "OutputOptions",
    "Termination",
    "PlacementResult",
    "FlowPolicy",
    "MAX_FLOWS_PER_CALL",
    "DagCache",
    "DagCacheMiss",
    "DagCacheStats",
    "CapacityOptions",
    "SensitivityOptions",
    "SensitivityEntry",
    "SensitivityResult",
    "capacity",
    "batch_capacity",
    "sensitivity",
    "Direction",
    "SpfOptions",
    "SpfResult",
    "UNREACHABLE",
    "spf",
    "KspOptions",
    "ksp",
    "resolve_to_paths",
}

# The shipped private entry: the flow core's _place_global_units.
PRIVATE = {"_place_global_units"}

# The private hooks (NETGRAPH_CORE_PRIVATE_HOOKS: the development build only): module functions
# and class members.
HOOK_FUNCTIONS = {
    "_max_flows_per_call",
    "_fail_python_allocation",
    "_fail_native_allocation",
    "_round_robin",
    "_highest_averages",
}
HOOK_MEMBERS = {
    "FlowLedger": {"_state"},
    "FlowPolicy": {"_place_traced", "_state", "_owner"},
}


def test_public_surface():
    assert set(ngc.__all__) == PUBLIC | {"__version__"}
    # Every name the package exposes, not only the ones __all__ lists: the public names and the
    # private metadata reader behind __version__ (dunder module attributes aside).
    assert {n for n in dir(ngc) if not n.startswith("__")} == PUBLIC | {"_version"}
    for name in ngc.__all__:
        assert hasattr(ngc, name)
    # Every public name is the compiled module's own object, re-exported.
    for name in PUBLIC:
        assert getattr(ngc, name) is getattr(_netgraph_core, name), name
    assert isinstance(ngc.__version__, str)
    assert ngc.MAX_FLOWS_PER_CALL == 10000


def test_extension_surface_is_exactly_the_public_and_private_sets():
    names = {n for n in dir(_netgraph_core) if not n.startswith("__")}
    assert names == PUBLIC | PRIVATE | (HOOK_FUNCTIONS if HOOKS else set())
    for name in PRIVATE | HOOK_FUNCTIONS:
        assert not hasattr(ngc, name)
    # Every private class member is a hook (beside pybind11's own conduit); the hooks are all
    # present or all absent.
    for cls in PUBLIC:
        obj = getattr(_netgraph_core, cls)
        private = {
            n
            for n in dir(obj)
            if n.startswith("_")
            and not n.startswith("__")
            and n != "_pybind11_conduit_v1_"
        }
        want = HOOK_MEMBERS.get(cls, set()) if HOOKS else set()
        assert private == want, cls


def test_signatures_name_registered_types_only():
    # pybind11 writes a signature's types when the binding is made: a type registered after a
    # binding that takes, returns or defaults to it appears under its C++ name (the
    # registration files run in an order where every such type comes first).
    docs = []
    for name, obj in vars(_netgraph_core).items():
        docs.append((name, getattr(obj, "__doc__", None)))
        if isinstance(obj, type):
            for member, raw in vars(obj).items():
                target = raw.__func__ if isinstance(raw, staticmethod) else raw
                docs.append((f"{name}.{member}", getattr(target, "__doc__", None)))
                if isinstance(raw, property) and raw.fget is not None:
                    docs.append((f"{name}.{member}.fget", raw.fget.__doc__))
    signed = [(n, d) for n, d in docs if isinstance(d, str) and "(" in d]
    assert len(signed) > 200  # not vacuous: every function, method and property
    for name, doc in signed:
        assert "::" not in doc, (name, doc.splitlines()[0])


# Retired spellings of this API's names: no alias.
RETIRED = {
    "StrictMultiDiGraph",
    "FlowGraph",
    "Outputs",
    "Diminishing",
    "kMaxFlowsPerCall",
}

# The retired search entry points, their option types and the private name _spf: never bound.
RETIRED_SEARCH = {
    "shortest_paths",
    "shortest_paths_to",
    "spf_dag",
    "spf_to",
    "spf_reverse",
    "k_shortest_paths",
    "_spf",
    "EdgeSelection",
    "PathAlg",
}


def test_retired_spellings_are_absent():
    for name in RETIRED | RETIRED_SEARCH:
        assert not hasattr(_netgraph_core, name), name
        assert not hasattr(ngc, name), name
    assert ngc.Graph.__name__ == "Graph"
    assert ngc.FlowLedger.__name__ == "FlowLedger"
    assert ngc.OutputOptions.__name__ == "OutputOptions"
    assert ngc.DiminishingReturns.__name__ == "DiminishingReturns"
    # PredDAG and Path carry no resolver or alternative constructors.
    for name in ("resolve_to_paths", "from_edges", "from_hops"):
        assert not hasattr(ngc.PredDAG, name), name
    for name in ("from_edges", "from_hops"):
        assert not hasattr(ngc.Path, name), name
    # The search queue is C++ only.
    assert not hasattr(ngc, "SpfQueue")
    assert ngc.UNREACHABLE == 2**63 - 1


def test_policy_takes_its_ledger_first_or_by_the_name_ledger():
    cfg = config(ngc.CostOnly(), selection=sel(prio=RES))
    g = ngc.Graph(2, [0], [1], [3], [1], 1.0)
    ledger = ngc.FlowLedger(g)
    p = ngc.FlowPolicy(ledger=ledger, sources=0, destinations=1, config=cfg)
    assert p.place(2).placed == 2
    assert ledger.records()[0].delivered == 2
    # Positionally the ledger comes first.
    assert ngc.FlowPolicy(ledger, 0, 1, cfg).place(1).placed == 1
    with pytest.raises(TypeError):
        ngc.FlowPolicy(flow_graph=ledger, sources=0, destinations=1, config=cfg)
    with pytest.raises(ValueError, match=exact("ledger must be a FlowLedger")):
        ngc.FlowPolicy(None, 0, 1, cfg)
    assert "ledger: _netgraph_core.FlowLedger" in ngc.FlowPolicy.__init__.__doc__


def _refusal_inputs():
    """What the refused calls below share: the chain 0 -> 1 -> 2 and valid options."""
    g = ngc.Graph(3, [0, 1], [1, 2], [5, 5], [1, 1], 1.0)
    cfg = config(ngc.CostOnly())
    return SimpleNamespace(
        g=g,
        ledger=ngc.FlowLedger(g),
        cfg=cfg,
        capacity=ngc.CapacityOptions(cfg),
        sensitivity=ngc.SensitivityOptions(cfg),
        spf=ngc.SpfOptions(ngc.Direction.FORWARD, sel()),
        ksp=ngc.KspOptions(1, sel()),
        dag=ngc.Path(g, [[0], [1]]).to_dag(),
        path=ngc.Path(g, [[0]]),
    )


GRAPH = "graph must be a Graph"
INTERSECT = "Eligibility.intersect takes two Eligibility values"
PAIR = "each pair must be (sources, destinations)"
ROUTING = "routing must be CostOnly, CapacityAware, Global or Pinned (got {})"
SEQUENCE = "Pinned bundles must be a sequence of PredDAG"
BUNDLE = "Pinned bundles must be PredDAG values (a Path converts with to_dag())"

# Each wrong argument, the exception class and the whole message it is refused with. A None
# where a Graph, an Eligibility or a PredDAG is expected reaches the binding as a null value
# (pybind11 refuses every other wrong type first), so these checks also keep the call from
# reading it.
REFUSALS = {
    "spf-graph-None": (lambda s: ngc.spf(None, 0, s.spf), ValueError, GRAPH),
    "ksp-graph-None": (lambda s: ngc.ksp(None, 0, 2, s.ksp), ValueError, GRAPH),
    "resolve-graph-None": (
        lambda s: ngc.resolve_to_paths(None, s.dag, 0, 2),
        ValueError,
        GRAPH,
    ),
    "resolve-dag-None": (
        lambda s: ngc.resolve_to_paths(s.g, None, 0, 2),
        ValueError,
        "dag must be a PredDAG",
    ),
    "DagCache-graph-None": (lambda s: ngc.DagCache(None), ValueError, GRAPH),
    "FlowLedger-graph-None": (lambda s: ngc.FlowLedger(None), ValueError, GRAPH),
    "FlowLedger-graph-None-with-residual": (
        lambda s: ngc.FlowLedger(None, [5, 5]),
        ValueError,
        GRAPH,
    ),
    "Path-graph-None": (lambda s: ngc.Path(None, [[0]]), ValueError, GRAPH),
    "from_csr-graph-None": (
        lambda s: ngc.PredDAG.from_csr(None, [0, 0, 0, 0], [], []),
        ValueError,
        GRAPH,
    ),
    "capacity-graph-None": (
        lambda s: ngc.capacity(None, 0, 2, s.capacity),
        ValueError,
        GRAPH,
    ),
    "batch_capacity-graph-None": (
        lambda s: ngc.batch_capacity(None, [(0, 2)], s.capacity),
        ValueError,
        GRAPH,
    ),
    "sensitivity-graph-None": (
        lambda s: ngc.sensitivity(None, 0, 2, s.sensitivity),
        ValueError,
        GRAPH,
    ),
    "intersect-second-None": (
        lambda s: ngc.Eligibility.intersect(s.g.eligibility(), None),
        ValueError,
        INTERSECT,
    ),
    "intersect-first-None": (
        lambda s: ngc.Eligibility.intersect(None, s.g.eligibility()),
        ValueError,
        INTERSECT,
    ),
    # Shapes and types of sequence and field arguments.
    "batch_capacity-pairs-not-a-sequence": (
        lambda s: ngc.batch_capacity(s.g, 5, s.capacity),
        ValueError,
        "pairs must be a sequence of (sources, destinations)",
    ),
    "batch_capacity-pair-of-one": (
        lambda s: ngc.batch_capacity(s.g, [(0,)], s.capacity),
        ValueError,
        PAIR,
    ),
    "batch_capacity-pair-not-a-sequence": (
        lambda s: ngc.batch_capacity(s.g, [5], s.capacity),
        ValueError,
        PAIR,
    ),
    "batch_capacity-per_pair_elig-not-a-sequence": (
        lambda s: ngc.batch_capacity(s.g, [(0, 2)], s.capacity, 5),
        ValueError,
        "per_pair_elig must be a sequence of Eligibility or None",
    ),
    "Path-hops-not-a-sequence": (
        lambda s: ngc.Path(s.g, 5),
        ValueError,
        "Path: hops must be a sequence of hops (sequences of EdgeIds)",
    ),
    "FlowPolicy-sources-2-D": (
        lambda s: ngc.FlowPolicy(s.ledger, np.array([[0]]), 2, s.cfg),
        ValueError,
        "sources must be a 1-D array (got 2 dimensions)",
    ),
    "FlowLedger-record-unknown-ledger": (
        lambda s: s.ledger.record((-1, 0)),
        ValueError,
        "FlowId (-1, 0) is retired or unknown",
    ),
    "config-routing-int": (
        lambda s: config(5),
        TypeError,
        ROUTING.format("int"),
    ),
    "config-routing-None": (
        lambda s: config(None),
        TypeError,
        ROUTING.format("NoneType"),
    ),
    "CapacityAware-flow_count-None": (
        lambda s: ngc.CapacityAware(flow_count=None),
        ValueError,
        "flow_count must be an integer (got None)",
    ),
    "Absolute-cost-None": (
        lambda s: ngc.Absolute(None),
        ValueError,
        "Absolute cost must be an integer (got None)",
    ),
    "DiminishingReturns-window-None": (
        lambda s: ngc.DiminishingReturns(None, 0.5),
        ValueError,
        "diminishing window must be an integer (got None)",
    ),
    "Stretch-factor-str": (
        lambda s: ngc.Stretch("1.5"),
        ValueError,
        "Stretch factor must be a number (got str)",
    ),
    "Stretch-factor-bool": (
        lambda s: ngc.Stretch(True),
        ValueError,
        "Stretch factor must be a number (got bool)",
    ),
    "DiminishingReturns-fraction-str": (
        lambda s: ngc.DiminishingReturns(1, "0.5"),
        ValueError,
        "diminishing fraction must be a number (got str)",
    ),
    "Pinned-bundles-int": (lambda s: ngc.Pinned(5), TypeError, SEQUENCE),
    "Pinned-bundles-None": (lambda s: ngc.Pinned(None), TypeError, SEQUENCE),
    "Pinned-bundle-int": (lambda s: ngc.Pinned([5]), TypeError, BUNDLE),
    "Pinned-bundle-None": (lambda s: ngc.Pinned([None]), TypeError, BUNDLE),
    "Pinned-bundle-Path": (lambda s: ngc.Pinned([s.path]), TypeError, BUNDLE),
}


@pytest.mark.parametrize(
    ("call", "error", "message"), REFUSALS.values(), ids=REFUSALS.keys()
)
def test_a_wrong_argument_is_refused_with_its_message(call, error, message):
    s = _refusal_inputs()
    with pytest.raises(error, match=exact(message)) as info:
        call(s)
    assert type(info.value) is error
