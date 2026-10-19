"""The public surface: the graph on units, Eligibility, the conversions, the route structures,
the selection and placement axes, the ledger, the policy engine, the DagCache, the evaluators
and the public search; the extension's private entry, which netgraph_core never re-exports;
the retired spellings and the retired search entry points, which exist nowhere, without
aliases; and the argument checks of the bindings, each with its exception class and whole
message."""

import re
from pathlib import Path
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
    "Origin",
    "FlowLedger",
    "FlowRecord",
    "Absolute",
    "Stretch",
    "Tunnels",
    "CostOnly",
    "CapacityAware",
    "Global",
    "Pinned",
    "WeightBasis",
    "TunnelStatus",
    "Tunnel",
    "Route",
    "ForwardingProgram",
    "TunnelFlow",
    "TunnelPlan",
    "FlowPolicyConfig",
    "OutputOptions",
    "Termination",
    "PlacementResult",
    "SourceFlow",
    "FlowPolicy",
    "DagCache",
    "DagCacheMiss",
    "DagCacheStats",
    "QueryOptions",
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
    # 52 entries with __version__ (DiminishingReturns and MAX_FLOWS_PER_CALL are gone, Tunnels
    # and TunnelPlan are new).
    assert len(ngc.__all__) == 52


def test_enum_values_are_pinned():
    # The integer values of the public enumerations (contract 14.1, 14.2).
    assert {m.name: int(m) for m in ngc.TunnelStatus.__members__.values()} == {
        "UP": 0,
        "NO_BANDWIDTH": 1,
        "COST_CEILING": 2,
        "FAILED": 3,
        "BLACKHOLED": 4,
    }
    assert {m.name: int(m) for m in ngc.Termination.__members__.values()} == {
        "COMPLETE": 0,
        "NO_ROUTE": 1,
        "EXHAUSTED": 2,
        "NO_ACTIVE_WEIGHT": 3,
    }


def test_path_selection_defaults():
    # multipath None is the model's natural shape (contract 3.3).
    s = ngc.PathSelection()
    assert (s.tie_break, s.priority, s.multipath, s.multi_edge) == (
        ngc.EdgeTieBreak.DETERMINISTIC,
        ngc.EdgePriority.CAPACITY,
        None,
        True,
    )
    assert repr(s) == (
        "PathSelection(tie_break=DETERMINISTIC, priority=CAPACITY, multipath=None, multi_edge=True)"
    )
    assert repr(ngc.PathSelection(multipath=False)).count("multipath=False") == 1
    assert hash(ngc.PathSelection()) == hash(ngc.PathSelection(multipath=None))
    assert ngc.PathSelection(multipath=None) != ngc.PathSelection(multipath=True)


def test_tunnel_values_read_only():
    # The T1 configuration values have read-only fields and no T0 keyword.
    t = ngc.Tunnels(2)
    ca = ngc.CapacityAware(t, max_path_cost=ngc.Absolute(3))
    with pytest.raises(AttributeError):
        t.count = 3
    with pytest.raises(AttributeError):
        t.up_to = True
    with pytest.raises(AttributeError):
        ca.tunnels = ngc.Tunnels(3)
    with pytest.raises(AttributeError):
        ca.min_unreserved_bandwidth = 1
    with pytest.raises(TypeError):
        ngc.CapacityAware(flow_count=2)
    assert (
        repr(t) == "Tunnels(2)"
        and repr(ngc.Tunnels(4, up_to=True)) == "Tunnels(4, up_to=True)"
    )
    assert (
        repr(ca)
        == "CapacityAware(Tunnels(2), max_path_cost=Absolute(3), min_unreserved_bandwidth=0)"
    )
    assert "max_path_cost=None," in repr(ngc.CapacityAware(t))
    assert repr(
        ngc.CapacityAware(
            ngc.Tunnels(3), max_path_cost=ngc.Stretch(2.0), min_unreserved_bandwidth=5
        )
    ) == (
        "CapacityAware(Tunnels(3), max_path_cost=Stretch(2.0), min_unreserved_bandwidth=5)"
    )
    assert (
        t == ngc.Tunnels(2)
        and hash(t) == hash(ngc.Tunnels(2))
        and t != ngc.Tunnels(2, up_to=True)
    )


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
    "DiminishingReturns",
    "MAX_FLOWS_PER_CALL",
    "kMaxFlowsPerCall",
    "_max_flows_per_call",
    "CapacityOptions",
    "SensitivityOptions",
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


# Retired members: PlacementResult's per-source termination list (the result reports one
# termination and per-source rows in by_source), FlowRecord's kernel (only ever written),
# FlowPolicy's usable-bundle count and Pinned's bundle list (Pinned forwards a ForwardingProgram;
# FlowPolicy.program and PlacementResult.by_tunnel report what is forwarded), the growth
# computation's CapacityAware fields and terminations (T1 sets up its tunnels once per call).
RETIRED_MEMBERS = {
    "PlacementResult": {"local_outcomes"},
    "FlowRecord": {"placement"},
    "FlowPolicy": {"usable_bundles"},
    "Pinned": {"bundles"},
    "CapacityAware": {
        "flow_count",
        "max_flow_count",
        "min_flow_bandwidth",
        "reoptimize",
        "diminishing",
    },
    "Termination": {"LIMIT_FLOWS", "LIMIT_DIMINISHING"},
}


# The files that tell a hooks build from a shipped one by probing one private hook.
HOOK_PROBES = (
    ".github/workflows/tests.yml",
    ".github/workflows/release.yml",
    "dev/run-checks.sh",
    "dev/stubs/generate.py",
    "tests/py/private_hooks.py",
)


def test_hook_probes_name_a_bound_hook():
    # A probe of a deleted hook makes a hooks build look like a shipped one: CI then refuses a
    # correct editable build, and the release check that a wheel has no hooks can never fail.
    root = Path(__file__).resolve().parents[2]
    if not all((root / rel).is_file() for rel in HOOK_PROBES):
        pytest.skip("needs a git checkout: the sdist leaves out .github/ and dev/")
    names = {}
    for rel in HOOK_PROBES:
        found = re.findall(r"hasattr\(\w+, ['\"](\w+)['\"]\)", (root / rel).read_text())
        assert found, rel
        for name in found:
            names.setdefault(name, []).append(rel)
    assert set(names) == {"_fail_native_allocation"}, names
    assert not HOOKS or all(hasattr(_netgraph_core, name) for name in names)


def test_retired_spellings_are_absent():
    for name in RETIRED | RETIRED_SEARCH:
        assert not hasattr(_netgraph_core, name), name
        assert not hasattr(ngc, name), name
    for cls, members in RETIRED_MEMBERS.items():
        for member in members:
            assert not hasattr(getattr(ngc, cls), member), (cls, member)
    # FlowPolicyConfig takes keywords only: no positional spelling.
    with pytest.raises(TypeError):
        ngc.FlowPolicyConfig(ngc.CostOnly(), ngc.FlowPlacement.MAX_FLOW)
    assert ngc.Graph.__name__ == "Graph"
    assert ngc.FlowLedger.__name__ == "FlowLedger"
    assert ngc.OutputOptions.__name__ == "OutputOptions"
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
        query=ngc.QueryOptions(cfg),
        spf=ngc.SpfOptions(ngc.Direction.FORWARD, sel()),
        ksp=ngc.KspOptions(1, sel()),
        dag=ngc.Path(g, [[0], [1]]).to_dag(),
        route=ngc.Route(0, ngc.Path(g, [[0], [1]])),
    )


GRAPH = "graph must be a Graph"
INTERSECT = "Eligibility.intersect takes two Eligibility values"
PAIR = "each pair must be (sources, destinations)"
ROUTING = "routing must be CostOnly, CapacityAware, Global or Pinned (got {})"
ROUTES = "routes must be a sequence of Route"
BUNDLE = "Route bundle must be a Path or a PredDAG (got {})"

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
        lambda s: ngc.capacity(None, 0, 2, s.query),
        ValueError,
        GRAPH,
    ),
    "batch_capacity-graph-None": (
        lambda s: ngc.batch_capacity(None, [(0, 2)], s.query),
        ValueError,
        GRAPH,
    ),
    "sensitivity-graph-None": (
        lambda s: ngc.sensitivity(None, 0, 2, s.query),
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
        lambda s: ngc.batch_capacity(s.g, 5, s.query),
        ValueError,
        "pairs must be a sequence of (sources, destinations)",
    ),
    "batch_capacity-pair-of-one": (
        lambda s: ngc.batch_capacity(s.g, [(0,)], s.query),
        ValueError,
        PAIR,
    ),
    "batch_capacity-pair-not-a-sequence": (
        lambda s: ngc.batch_capacity(s.g, [5], s.query),
        ValueError,
        PAIR,
    ),
    "batch_capacity-per_pair_elig-not-a-sequence": (
        lambda s: ngc.batch_capacity(s.g, [(0, 2)], s.query, 5),
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
    "Tunnels-count-None": (
        lambda s: ngc.Tunnels(None),
        ValueError,
        "tunnel count must be an integer (got None)",
    ),
    "Tunnels-count-bool": (
        lambda s: ngc.Tunnels(True),
        ValueError,
        "tunnel count must be an integer (got bool)",
    ),
    "CapacityAware-min_unreserved_bandwidth-None": (
        lambda s: ngc.CapacityAware(ngc.Tunnels(1), min_unreserved_bandwidth=None),
        ValueError,
        "min_unreserved_bandwidth must be an integer (got None)",
    ),
    "Absolute-cost-None": (
        lambda s: ngc.Absolute(None),
        ValueError,
        "Absolute cost must be an integer (got None)",
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
    "Pinned-program-None": (
        lambda s: ngc.Pinned(None),
        TypeError,
        "Pinned program must be a ForwardingProgram",
    ),
    "Route-bundle-int": (lambda s: ngc.Route(0, 5), TypeError, BUNDLE.format("int")),
    "Route-bundle-None": (
        lambda s: ngc.Route(0, None),
        TypeError,
        BUNDLE.format("NoneType"),
    ),
    "from_routes-graph-None": (
        lambda s: ngc.ForwardingProgram.from_routes(None, [0], [2], [s.route]),
        ValueError,
        GRAPH,
    ),
    "from_routes-routes-int": (
        lambda s: ngc.ForwardingProgram.from_routes(s.g, [0], [2], 5),
        TypeError,
        ROUTES,
    ),
    "from_routes-route-int": (
        lambda s: ngc.ForwardingProgram.from_routes(s.g, [0], [2], [5]),
        TypeError,
        ROUTES,
    ),
    "from_routes-route-None": (
        lambda s: ngc.ForwardingProgram.from_routes(s.g, [0], [2], [None]),
        TypeError,
        ROUTES,
    ),
}


@pytest.mark.parametrize(
    ("call", "error", "message"), REFUSALS.values(), ids=REFUSALS.keys()
)
def test_a_wrong_argument_is_refused_with_its_message(call, error, message):
    s = _refusal_inputs()
    with pytest.raises(error, match=exact(message)) as info:
        call(s)
    assert type(info.value) is error


def test_t2_members():
    # Global.tunnels, Global.prefer_headroom, Tunnel.planned and the T2 fields of TunnelPlan,
    # read-only.
    gl = ngc.Global(tunnels=ngc.Tunnels(16), prefer_headroom=True)
    assert (gl.tunnels, gl.prefer_headroom) == (ngc.Tunnels(16), True)
    assert (ngc.Global().tunnels, ngc.Global().prefer_headroom) == (None, False)
    with pytest.raises(AttributeError):
        gl.tunnels = ngc.Tunnels(2)
    with pytest.raises(AttributeError):
        gl.prefer_headroom = False
    for name in (
        "live_heads",
        "insufficient_resolution",
        "fluid_share",
        "common_share",
        "planned_per_tunnel",
        "planned_total",
        "route_cost",
        "headroom_penalty",
        "forwarding_bound",
    ):
        assert isinstance(getattr(ngc.TunnelPlan, name), property), name
        assert getattr(ngc.TunnelPlan, name).fset is None, name
    assert isinstance(ngc.Tunnel.planned, property) and ngc.Tunnel.planned.fset is None
