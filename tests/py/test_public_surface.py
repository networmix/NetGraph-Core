"""The public surface: the graph on units, Eligibility, the conversions, the route structures,
the route choice and placement axes, the ledger, the policy engine, the ShortestPathCache, the evaluators
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
from helpers import HIGHER_AVAIL, assert_read_only, config, exact
from private_hooks import HOOKS

import netgraph_core as ngc

PUBLIC = {
    "Graph",
    "Eligibility",
    "to_units",
    "as_float",
    "PredecessorDAG",
    "Path",
    "RouteTieBreak",
    "RouteAllocation",
    "Admission",
    "SourceAllocation",
    "FlowLedger",
    "FlowRecord",
    "AbsoluteCost",
    "CostStretch",
    "TunnelCount",
    "ShortestPaths",
    "DistributedTE",
    "CentralizedTE",
    "Pinned",
    "WeightBasis",
    "TunnelStatus",
    "Tunnel",
    "RouteDefinition",
    "ForwardingProgram",
    "TunnelFlow",
    "TunnelPlan",
    "FlowPolicyConfig",
    "OutputOptions",
    "Termination",
    "PlacementResult",
    "SourceFlow",
    "FlowPolicy",
    "ShortestPathCache",
    "ShortestPathCacheMiss",
    "ShortestPathCacheStats",
    "EvaluationOptions",
    "SensitivityEntry",
    "SensitivityResult",
    "evaluate_demand",
    "evaluate_pairs",
    "evaluate_edge_failures",
    "Direction",
    "SpfOptions",
    "SpfResult",
    "UNREACHABLE",
    "spf",
    "KspOptions",
    "ksp",
    "resolve_to_paths",
}

# The shipped private entry: the flow core's _solve_min_cost_flow.
PRIVATE = {"_solve_min_cost_flow"}

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
    # 50 entries with __version__ (RouteTieBreak replaces the three earlier route-choice types).
    assert len(ngc.__all__) == 50


def test_enum_values_are_pinned():
    # The promised integer values of the public enumerations (docs/semantics.md, "Forwarding
    # programs" and "Termination").
    assert {m.name: int(m) for m in ngc.TunnelStatus.__members__.values()} == {
        "UP": 0,
        "INSUFFICIENT_BANDWIDTH": 1,
        "COST_CEILING_EXCEEDED": 2,
        "FAILED": 3,
        "BLACKHOLED": 4,
    }
    assert {m.name: int(m) for m in ngc.Termination.__members__.values()} == {
        "COMPLETE": 0,
        "NO_ROUTE": 1,
        "EXHAUSTED": 2,
        "NO_ACTIVE_WEIGHT": 3,
    }


def test_enum_declaration_order():
    # Every enumeration lists its members in the promised order (the interface listing,
    # contract_api.py's source); no other integer value is promised.
    order = {
        ngc.RouteTieBreak: [
            "FEWER_HOPS",
            "HIGHER_CAPACITY",
            "HIGHER_AVAILABLE_CAPACITY",
            "HIGHER_HEADROOM",
        ],
        ngc.RouteAllocation: ["EQUAL", "WEIGHTED", "OPTIMIZED"],
        ngc.SourceAllocation: ["EQUAL", "OPTIMIZED"],
        ngc.WeightBasis: ["CAPACITY", "EXPLICIT", "RESIDUAL_SNAPSHOT"],
        ngc.Admission: ["LOSSLESS", "LOSSY"],
        ngc.Direction: ["FORWARD", "REVERSE"],
        ngc.TunnelStatus: [
            "UP",
            "INSUFFICIENT_BANDWIDTH",
            "COST_CEILING_EXCEEDED",
            "FAILED",
            "BLACKHOLED",
        ],
        ngc.Termination: ["COMPLETE", "NO_ROUTE", "EXHAUSTED", "NO_ACTIVE_WEIGHT"],
    }
    for enum, names in order.items():
        assert list(enum.__members__) == names, enum.__name__
        values = [int(enum.__members__[n]) for n in names]
        assert values == sorted(values), enum.__name__


def test_route_choice_defaults_and_strict_bools():
    # ShortestPaths() is ECMP: FEWER_HOPS, every equal-cost route, every parallel edge.
    # The route choice is keywords only; a bool field takes bool or numpy.bool_ alone, so None
    # (no longer a multipath value) and 1 are TypeErrors.
    sp = ngc.ShortestPaths()
    assert (sp.tie_break, sp.multipath, sp.parallel_edges) == (
        ngc.RouteTieBreak.FEWER_HOPS,
        True,
        True,
    )
    assert repr(sp) == (
        "ShortestPaths(allocation=RouteAllocation.EQUAL, source_allocation=None,"
        " weight_basis=None, edge_weights=None, tie_break=RouteTieBreak.FEWER_HOPS,"
        " multipath=True, parallel_edges=True)"
    )
    assert_read_only(
        sp,
        "allocation",
        "source_allocation",
        "weight_basis",
        "edge_weights",
        "tie_break",
        "multipath",
        "parallel_edges",
    )
    one = ngc.ShortestPaths(multipath=np.bool_(False))
    assert one.multipath is False
    ca = ngc.DistributedTE(ngc.TunnelCount(1))
    assert (ca.tie_break, ca.parallel_edges) == (ngc.RouteTieBreak.FEWER_HOPS, True)
    assert_read_only(ca, "tie_break", "parallel_edges")
    with pytest.raises(TypeError):
        ngc.ShortestPaths(ngc.RouteTieBreak.HIGHER_CAPACITY)
    with pytest.raises(TypeError):
        ngc.DistributedTE(ngc.TunnelCount(1), multipath=False)  # T1 is single-path
    for bad in (None, 1, "no"):
        with pytest.raises(TypeError):
            ngc.ShortestPaths(multipath=bad)
        with pytest.raises(TypeError):
            ngc.ShortestPaths(parallel_edges=bad)
        with pytest.raises(TypeError):
            ngc.DistributedTE(ngc.TunnelCount(1), parallel_edges=bad)
        with pytest.raises(TypeError):
            ngc.OutputOptions(edge_flows=bad)
        with pytest.raises(TypeError):
            ngc.TunnelCount(2, up_to=bad)
        with pytest.raises(TypeError):
            ngc.CentralizedTE(
                source_allocation=ngc.SourceAllocation.EQUAL, prefer_headroom=bad
            )


def test_tunnel_values_read_only():
    # The T1 configuration values have read-only fields and no T0 keyword.
    t = ngc.TunnelCount(2)
    ca = ngc.DistributedTE(t, cost_ceiling=ngc.AbsoluteCost(3))
    assert_read_only(t, "per_head", "up_to")
    assert_read_only(ca, "tunnel_count", "cost_ceiling", "min_unreserved_bandwidth")
    with pytest.raises(TypeError):
        ngc.DistributedTE(flow_count=2)
    assert (
        repr(t) == "TunnelCount(2)"
        and repr(ngc.TunnelCount(4, up_to=True)) == "TunnelCount(4, up_to=True)"
    )
    assert repr(ca) == (
        "DistributedTE(TunnelCount(2), cost_ceiling=AbsoluteCost(3), min_unreserved_bandwidth=0,"
        " tie_break=RouteTieBreak.FEWER_HOPS, parallel_edges=True)"
    )
    assert "cost_ceiling=None," in repr(ngc.DistributedTE(t))
    assert repr(
        ngc.DistributedTE(
            ngc.TunnelCount(3),
            cost_ceiling=ngc.CostStretch(2.0),
            min_unreserved_bandwidth=5,
            tie_break=ngc.RouteTieBreak.HIGHER_HEADROOM,
            parallel_edges=False,
        )
    ) == (
        "DistributedTE(TunnelCount(3), cost_ceiling=CostStretch(Fraction(2, 1)), min_unreserved_bandwidth=5,"
        " tie_break=RouteTieBreak.HIGHER_HEADROOM, parallel_edges=False)"
    )
    assert (
        t == ngc.TunnelCount(2)
        and hash(t) == hash(ngc.TunnelCount(2))
        and t != ngc.TunnelCount(2, up_to=True)
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
    "PathSelection",
    "EdgeTieBreak",
    "EdgePriority",
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
# computation's DistributedTE fields and terminations (T1 sets up its tunnels once per call).
RETIRED_MEMBERS = {
    "PlacementResult": {"local_outcomes"},
    "FlowRecord": {"placement"},
    "FlowPolicy": {"usable_bundles"},
    "Pinned": {"bundles"},
    "DistributedTE": {
        "flow_count",
        "max_flow_count",
        "min_flow_bandwidth",
        "reoptimize",
        "diminishing",
        "selection",
    },
    "FlowPolicyConfig": {"selection", "placement", "origin", "weights", "edge_weights"},
    "SpfOptions": {"selection"},
    "KspOptions": {"selection", "multipath"},
    "Termination": {"LIMIT_FLOWS", "LIMIT_DIMINISHING"},
}


# The files that tell a hooks build from a shipped one by probing one private hook.
HOOK_PROBES = (
    "dev/pytest_guard.py",
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
        ngc.FlowPolicyConfig(ngc.ShortestPaths(), ngc.Admission.LOSSLESS)
    assert ngc.Graph.__name__ == "Graph"
    assert ngc.FlowLedger.__name__ == "FlowLedger"
    assert ngc.OutputOptions.__name__ == "OutputOptions"
    # PredecessorDAG and Path carry no resolver or alternative constructors.
    for name in ("resolve_to_paths", "from_edges", "from_hops"):
        assert not hasattr(ngc.PredecessorDAG, name), name
    for name in ("from_edges", "from_hops"):
        assert not hasattr(ngc.Path, name), name
    # The search queue is C++ only.
    assert not hasattr(ngc, "SpfQueue")
    assert ngc.UNREACHABLE == 2**63 - 1


def test_policy_takes_its_ledger_first_or_by_the_name_ledger():
    cfg = config(ngc.ShortestPaths(tie_break=HIGHER_AVAIL, parallel_edges=False))
    g = ngc.Graph(2, [0], [1], [3], [1], 1.0)
    ledger = ngc.FlowLedger(g)
    p = ngc.FlowPolicy(ledger=ledger, sources=0, destinations=1, config=cfg)
    assert p.place(2).call_delivered == 2
    assert ledger.records()[0].delivered == 2
    # Positionally the ledger comes first.
    assert ngc.FlowPolicy(ledger, 0, 1, cfg).place(1).call_delivered == 1
    with pytest.raises(TypeError):
        ngc.FlowPolicy(flow_graph=ledger, sources=0, destinations=1, config=cfg)
    with pytest.raises(TypeError, match=exact("ledger must be a FlowLedger")) as info:
        ngc.FlowPolicy(None, 0, 1, cfg)
    assert type(info.value) is TypeError
    assert "ledger: _netgraph_core.FlowLedger" in ngc.FlowPolicy.__init__.__doc__


def _refusal_inputs():
    """What the refused calls below share: the chain 0 -> 1 -> 2 and valid options."""
    g = ngc.Graph(3, [0, 1], [1, 2], [5, 5], [1, 1], 1.0)
    cfg = config(ngc.ShortestPaths())
    return SimpleNamespace(
        g=g,
        ledger=ngc.FlowLedger(g),
        cfg=cfg,
        query=ngc.EvaluationOptions(cfg),
        spf=ngc.SpfOptions(ngc.Direction.FORWARD),
        ksp=ngc.KspOptions(1),
        dag=ngc.Path(g, [[0], [1]]).to_dag(),
        route=ngc.RouteDefinition(0, ngc.Path(g, [[0], [1]])),
    )


GRAPH = "graph must be a Graph"
INTERSECT = "Eligibility.intersect takes two Eligibility values"
PAIR = "each pair must be (sources, destinations)"
ROUTING = (
    "routing must be ShortestPaths, DistributedTE, CentralizedTE or Pinned (got {})"
)
ROUTES = "routes must be a sequence of RouteDefinition"
ROUTE_TYPE = "RouteDefinition route must be a Path or a PredecessorDAG (got {})"

# Each wrong argument, the exception class and the whole message it is refused with: a wrong type
# (None included) TypeError, a wrong value (a pair of another length, an unknown FlowId)
# ValueError. A None where a Graph, an Eligibility or a PredecessorDAG is expected reaches the
# binding as a null value (pybind11 refuses every other wrong type first), so these checks also keep
# the call from reading it.
REFUSALS = {
    "spf-graph-None": (lambda s: ngc.spf(None, 0, s.spf), TypeError, GRAPH),
    "ksp-graph-None": (lambda s: ngc.ksp(None, 0, 2, s.ksp), TypeError, GRAPH),
    "resolve-graph-None": (
        lambda s: ngc.resolve_to_paths(None, s.dag, 0, 2),
        TypeError,
        GRAPH,
    ),
    "resolve-dag-None": (
        lambda s: ngc.resolve_to_paths(s.g, None, 0, 2),
        TypeError,
        "dag must be a PredecessorDAG",
    ),
    "ShortestPathCache-graph-None": (
        lambda s: ngc.ShortestPathCache(None),
        TypeError,
        GRAPH,
    ),
    "FlowLedger-graph-None": (lambda s: ngc.FlowLedger(None), TypeError, GRAPH),
    "FlowLedger-graph-None-with-residual": (
        lambda s: ngc.FlowLedger(None, [5, 5]),
        TypeError,
        GRAPH,
    ),
    "Path-graph-None": (lambda s: ngc.Path(None, [[0]]), TypeError, GRAPH),
    "from_csr-graph-None": (
        lambda s: ngc.PredecessorDAG.from_csr(None, [0, 0, 0, 0], [], []),
        TypeError,
        GRAPH,
    ),
    "capacity-graph-None": (
        lambda s: ngc.evaluate_demand(None, 0, 2, s.query),
        TypeError,
        GRAPH,
    ),
    "evaluate_pairs-graph-None": (
        lambda s: ngc.evaluate_pairs(None, [(0, 2)], s.query),
        TypeError,
        GRAPH,
    ),
    "sensitivity-graph-None": (
        lambda s: ngc.evaluate_edge_failures(None, 0, 2, s.query),
        TypeError,
        GRAPH,
    ),
    "intersect-second-None": (
        lambda s: ngc.Eligibility.intersect(s.g.eligibility(), None),
        TypeError,
        INTERSECT,
    ),
    "intersect-first-None": (
        lambda s: ngc.Eligibility.intersect(None, s.g.eligibility()),
        TypeError,
        INTERSECT,
    ),
    # Shapes and types of sequence and field arguments.
    "evaluate_pairs-pairs-not-a-sequence": (
        lambda s: ngc.evaluate_pairs(s.g, 5, s.query),
        TypeError,
        "pairs must be a sequence of (sources, destinations)",
    ),
    "evaluate_pairs-pair-of-one": (
        lambda s: ngc.evaluate_pairs(s.g, [(0,)], s.query),
        ValueError,
        PAIR,
    ),
    "evaluate_pairs-pair-not-a-sequence": (
        lambda s: ngc.evaluate_pairs(s.g, [5], s.query),
        TypeError,
        PAIR,
    ),
    "evaluate_pairs-per_pair_eligibility-not-a-sequence": (
        lambda s: ngc.evaluate_pairs(s.g, [(0, 2)], s.query, 5),
        TypeError,
        "per_pair_eligibility must be a sequence of Eligibility or None",
    ),
    "Path-hops-not-a-sequence": (
        lambda s: ngc.Path(s.g, 5),
        TypeError,
        "Path: hops must be a sequence of hops (sequences of EdgeIds)",
    ),
    "FlowPolicy-sources-2-D": (
        lambda s: ngc.FlowPolicy(s.ledger, np.array([[0]]), 2, s.cfg),
        TypeError,
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
    "TunnelCount-per_head-None": (
        lambda s: ngc.TunnelCount(None),
        TypeError,
        "TunnelCount per_head must be an integer (got NoneType)",
    ),
    "TunnelCount-per_head-bool": (
        lambda s: ngc.TunnelCount(True),
        TypeError,
        "TunnelCount per_head must be an integer (got bool)",
    ),
    "DistributedTE-min_unreserved_bandwidth-None": (
        lambda s: ngc.DistributedTE(ngc.TunnelCount(1), min_unreserved_bandwidth=None),
        TypeError,
        "min_unreserved_bandwidth must be an integer (got NoneType)",
    ),
    "AbsoluteCost-cost-None": (
        lambda s: ngc.AbsoluteCost(None),
        TypeError,
        "AbsoluteCost cost must be an integer (got NoneType)",
    ),
    "Stretch-factor-str": (
        lambda s: ngc.CostStretch("1.5"),
        TypeError,
        "CostStretch factor must be a Fraction, a float or an integer (got str)",
    ),
    "Stretch-factor-bool": (
        lambda s: ngc.CostStretch(True),
        TypeError,
        "CostStretch factor must be a Fraction, a float or an integer (got bool)",
    ),
    "Pinned-program-None": (
        lambda s: ngc.Pinned(None),
        TypeError,
        "Pinned program must be a ForwardingProgram",
    ),
    "RouteDefinition-route-int": (
        lambda s: ngc.RouteDefinition(0, 5),
        TypeError,
        ROUTE_TYPE.format("int"),
    ),
    "RouteDefinition-route-None": (
        lambda s: ngc.RouteDefinition(0, None),
        TypeError,
        ROUTE_TYPE.format("NoneType"),
    ),
    "from_routes-graph-None": (
        lambda s: ngc.ForwardingProgram.from_routes(None, [0], [2], [s.route]),
        TypeError,
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
    # CentralizedTE.tunnel_count, CentralizedTE.prefer_headroom, Tunnel.planned_bandwidth and the T2 fields of TunnelPlan,
    # read-only.
    eq = ngc.SourceAllocation.EQUAL
    gl = ngc.CentralizedTE(
        source_allocation=eq, tunnel_count=ngc.TunnelCount(16), prefer_headroom=True
    )
    assert (gl.source_allocation, gl.tunnel_count, gl.prefer_headroom) == (
        eq,
        ngc.TunnelCount(16),
        True,
    )
    plain = ngc.CentralizedTE(source_allocation=eq)
    assert (plain.tunnel_count, plain.prefer_headroom) == (None, False)
    assert_read_only(
        gl, "source_allocation", "support", "tunnel_count", "prefer_headroom"
    )
    for name in (
        "heads_with_capacity",
        "insufficient_resolution",
        "max_fractional_bandwidth_per_tunnel",
        "max_bandwidth_per_tunnel",
        "planned_bandwidth_per_tunnel",
        "total_planned_bandwidth",
        "total_tunnel_cost",
        "headroom_penalty",
        "max_lossless_offer",
    ):
        assert isinstance(getattr(ngc.TunnelPlan, name), property), name
        assert getattr(ngc.TunnelPlan, name).fset is None, name
    assert (
        isinstance(ngc.Tunnel.planned_bandwidth, property)
        and ngc.Tunnel.planned_bandwidth.fset is None
    )


def _top_level_labels(text):
    # The `name=` labels of a repr at its own level, not those of a nested repr.
    labels, depth, word = [], 0, ""
    for ch in text:
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
        if depth == 1 and ch == "=" and word:
            labels.append(word)
        word = word + ch if (ch.isalnum() or ch == "_") else ""
    return labels


def test_repr_labels_name_current_attributes():
    # Every labelled __repr__ starts with its class name, and each of its own labels is an
    # attribute of that class, except the three that print a count (no such attribute).
    counts = {
        ("PredecessorDAG", "entries"),
        ("SpfResult", "dag_entries"),
        ("Pinned", "tunnels"),
    }
    g = ngc.Graph(4, [0, 0, 1, 2], [1, 2, 3, 3], [8, 8, 8, 8], [1, 1, 1, 1], 1.0)
    cfg = ngc.FlowPolicyConfig(
        routing=ngc.DistributedTE(ngc.TunnelCount(2)), admission=ngc.Admission.LOSSY
    )
    fg = ngc.FlowLedger(g)
    p = ngc.FlowPolicy(fg, 0, 3, cfg)
    r = p.replace(4, design_demand=4)
    search = ngc.spf(g, 0, ngc.SpfOptions(ngc.Direction.FORWARD))
    path = ngc.ksp(g, 0, 3, ngc.KspOptions(1))[0][1]
    objects = [
        r,
        r.plan,
        r.program.tunnels[0],
        r.by_tunnel[0],
        r.by_source[0],
        p.records()[0],
        ngc.RouteDefinition(0, path, weight=2),
        path,
        search,
        search.dag,
        cfg,
        cfg.routing,
        ngc.DistributedTE(ngc.TunnelCount(2), cost_ceiling=ngc.AbsoluteCost(3)),
        ngc.CentralizedTE(
            source_allocation=ngc.SourceAllocation.EQUAL,
            support=search.dag,
            tunnel_count=ngc.TunnelCount(2, up_to=True),
        ),
        ngc.CentralizedTE(
            source_allocation=ngc.SourceAllocation.EQUAL,
            tunnel_count=ngc.TunnelCount(2),
            prefer_headroom=True,
        ),
        ngc.ShortestPaths(
            allocation=ngc.RouteAllocation.WEIGHTED,
            weight_basis=ngc.WeightBasis.EXPLICIT,
            edge_weights=[1, 2, 3, 4],
        ),
        ngc.TunnelCount(2, up_to=True),
        ngc.Pinned(r.program),
        ngc.ShortestPaths(tie_break=HIGHER_AVAIL, multipath=False),
    ]
    for obj in objects:
        name = type(obj).__name__
        text = repr(obj)
        assert text.startswith(name + "("), text
        for label in _top_level_labels(text):
            assert (name, label) in counts or hasattr(type(obj), label), (text, label)
    # The labels that once named retired fields, exactly.
    assert repr(r).startswith("PlacementResult(call_delivered=4, termination=")
    assert repr(r.plan) == "TunnelPlan(tunnels_per_head=2, total_reserved_bandwidth=4)"
    assert repr(r.program.tunnels[0]).startswith("Tunnel(entry_id=0, ")
    assert repr(r.by_tunnel[0]).startswith("TunnelFlow(entry_id=0, ")
    assert repr(p.records()[0]).startswith("FlowRecord(flow_id=(")
    assert ", policy_id=" in repr(p.records()[0])
    assert repr(path).startswith("Path(source=0, destination=3, ")
    assert (
        repr(ngc.RouteDefinition(0, path, weight=2))
        == "RouteDefinition(head=0, weight=2)"
    )


def _stub_members():
    """The installed stub's public names: per class its members (a property, field or enum member is
    "attribute"; a callable, every overload's parameters as contract_api writes them)."""
    import ast
    import pathlib

    tree = ast.parse((pathlib.Path(ngc.__file__).parent / "__init__.pyi").read_text())

    def params(fn):
        a = fn.args
        pos = a.posonlyargs + a.args
        defaults = [None] * (len(pos) - len(a.defaults)) + list(a.defaults)
        out = [
            (arg.arg, "pos", None if d is None else ast.unparse(d))
            for arg, d in zip(pos, defaults, strict=True)
        ]
        if a.vararg:
            out.append(("*" + a.vararg.arg, "var", None))
        out += [
            (arg.arg, "kw", None if d is None else ast.unparse(d))
            for arg, d in zip(a.kwonlyargs, a.kw_defaults, strict=True)
        ]
        if a.kwarg:
            out.append(("**" + a.kwarg.arg, "varkw", None))
        return [p for p in out if p[0] != "self"]

    names = {}
    for node in tree.body:
        if isinstance(node, ast.FunctionDef):
            names.setdefault(node.name, {}).setdefault("__call__", []).append(
                params(node)
            )
        elif isinstance(node, ast.ClassDef):
            members = {}
            for item in node.body:
                if isinstance(item, ast.FunctionDef):
                    decorators = [ast.unparse(d) for d in item.decorator_list]
                    if "property" in decorators:
                        members[item.name] = "attribute"
                    else:
                        members.setdefault(item.name, []).append(params(item))
                elif isinstance(item, ast.AnnAssign):
                    members[item.target.id] = "attribute"
            names[node.name] = members
    return names


def _runtime_default(cls, method, param):
    # The stub prints an enum default as `...`; pybind11's runtime signature names the value.
    doc = getattr(getattr(ngc, cls), method).__doc__ or ""
    found = re.search(rf"\b{param}: [\w.]+ = ([\w.]+)", doc)
    return found.group(1) if found else None


def _runtime_signatures(owner, member):
    """Every overload's parameters (name, kind) as pybind11's runtime docstring states them."""
    target = (
        getattr(ngc, owner)
        if member == "__call__"
        else getattr(getattr(ngc, owner), member)
    )
    name = owner if member == "__call__" else member
    out = []
    for line in (target.__doc__ or "").splitlines():
        found = re.match(rf"^(?:\d+\. )?{re.escape(name)}\((.*)\) -> ", line)
        if not found:
            continue
        params, depth, part, kind = [], 0, "", "pos"
        for ch in found.group(1) + ",":
            if ch in "[(":
                depth += 1
            elif ch in "])":
                depth -= 1
            if ch == "," and depth == 0:
                part = part.strip()
                if part == "*":
                    kind = "kw"
                elif part and part != "/":
                    pname = part.split(":")[0].split("=")[0].strip()
                    if pname != "self":
                        params.append((pname, kind))
                part = ""
            else:
                part += ch
        out.append(params)
    return out


def test_runtime_signatures_match_the_contract():
    # The extension's own signatures (pybind11's), parameter names and kinds, against the listing:
    # a binding changed without its stub fails here, not only in the stub's drift check.
    from contract_api import CONTRACT_API

    differences = []
    for name, members in CONTRACT_API.items():
        for member, want in members.items():
            if want in ("attribute", "declared") or not hasattr(ngc, name):
                continue
            have = _runtime_signatures(name, member)
            for overload in want:
                if overload == [("_no_public_constructor", "pos", None)]:
                    continue
                if [(n, k) for n, k, _ in overload] not in have:
                    differences.append((name, member, overload, have))
    assert differences == []


def test_api_matches_the_contract():
    # Every public name of the interface listing (contract_api.py; the batch operations, not yet in
    # the library, excepted), every member, and every parameter's name, kind and default; no other public member.
    from contract_api import CONTRACT_API

    stub = _stub_members()
    differences = []
    for name, members in CONTRACT_API.items():
        if name not in stub:
            differences.append(("missing", name))
            continue
        have = stub[name]
        for member, want in members.items():
            if member not in have:
                differences.append(("missing", name, member))
            elif want in ("attribute", "declared"):
                if want == "attribute" and have[member] != "attribute":
                    differences.append(("kind", name, member))
            elif have[member] == "attribute":
                differences.append(("kind", name, member))
            else:
                for overload in want:
                    match = [
                        got
                        for got in have[member]
                        if [(n, k) for n, k, _ in got]
                        == [(n, k) for n, k, _ in overload]
                    ]
                    if not match:
                        differences.append(("parameters", name, member, overload))
                        continue
                    for (param, _, default), (_, _, got) in zip(
                        overload, match[0], strict=True
                    ):
                        if got != default and not (
                            got == "..."
                            and _runtime_default(name, member, param) == default
                        ):
                            differences.append(("default", name, member, param, got))
        extra = sorted(
            m
            for m in have
            if not m.startswith("_") and m not in members and m not in ("name", "value")
        )
        if extra:
            differences.append(("extra", name, extra))
    assert differences == []


def _public_docstrings():
    """(where, text) of every public docstring: the package's, each public class's (its
    constructor's text without the signature lines), each public method's and property's, and each
    function's."""
    import inspect

    out = [("netgraph_core", ngc.__doc__ or "")]

    def text(obj):
        doc = getattr(obj, "__doc__", None) or ""
        # pybind11 prefixes each overload's text with its signature: names, not prose.
        return "\n".join(
            line
            for line in doc.splitlines()
            if not re.match(r"^\s*(\d+\. )?\w+\(.*\) -> ", line)
        )

    for name in ngc.__all__:
        obj = getattr(ngc, name)
        if isinstance(obj, str):
            continue
        out.append((name, text(obj)))
        if inspect.isclass(obj) and not issubclass(obj, BaseException):
            for member, value in vars(obj).items():
                if member.startswith("_") and member != "__init__":
                    continue
                if isinstance(value, (staticmethod, classmethod)):
                    value = (
                        value.__func__
                    )  # the function's own text, not the descriptor's help
                if isinstance(value, property) or callable(value):
                    out.append((f"{name}.{member}", text(value)))
    return out


def _stub_names():
    import ast

    tree = ast.parse((Path(ngc.__file__).parent / "__init__.pyi").read_text())
    names = set(ngc.__all__)
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.ClassDef)):
            names.add(node.name)
        elif isinstance(node, ast.arg):
            names.add(node.arg)
        elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
            names.add(node.target.id)
    return names


def test_docstring_names_exist():
    # Every name a public docstring writes as code (a snake_case name with an underscore, or a
    # call `name()`), after a dot too (`ledger.carried()`, `options.edge_weights`), is a public
    # name, member or parameter: a renamed field cannot stay behind in the prose.
    known = _stub_names()
    missing = set()
    for where, doc in _public_docstrings():
        for m in re.finditer(
            r"(?<!\w)([a-z][a-z0-9]*(?:_[a-z0-9]+)+|[a-z]\w*(?=\(\)))", doc
        ):
            if m.group(1) not in known:
                missing.add((where, m.group(1)))
    assert not missing, sorted(missing)


def test_class_docstrings_name_every_field():
    # A class's text (its docstring and its constructor's) names each of its fields, so the text
    # that lists them follows a rename.
    import inspect

    missing = []
    for name in ngc.__all__:
        cls = getattr(ngc, name)
        if (
            not inspect.isclass(cls)
            or issubclass(cls, BaseException)
            or hasattr(cls, "__members__")
        ):
            continue
        fields = [
            m
            for m, v in vars(cls).items()
            if not m.startswith("_") and isinstance(v, property)
        ]
        init = getattr(cls.__init__, "__doc__", None) or ""
        init = "\n".join(
            line for line in init.splitlines() if not line.startswith("__init__(")
        )
        words = set(re.findall(r"\w+", (cls.__doc__ or "") + init))
        missing += [f"{name}.{f}" for f in fields if f not in words]
    assert not missing, missing
