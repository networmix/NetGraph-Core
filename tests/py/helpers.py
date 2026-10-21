"""Builders and constants shared by the Python tests (not collected: no test_ prefix): graphs,
EdgeId mappings, DAGs, the two policy configurations most tests use, results
and ledgers as comparable data, the sweep of a call's Python object allocations, and the
pytest.raises pattern of an exact message."""

from __future__ import annotations

import math
import re

import _netgraph_core
import pytest

import netgraph_core as ngc

INF = math.inf
L62 = 1 << 62
FEWER = ngc.RouteTieBreak.FEWER_HOPS
HIGHER_CAP = ngc.RouteTieBreak.HIGHER_CAPACITY
HIGHER_AVAIL = ngc.RouteTieBreak.HIGHER_AVAILABLE_CAPACITY
HIGHER_HEADROOM = ngc.RouteTieBreak.HIGHER_HEADROOM
# Every optional output of a placement or an evaluator.
ALL_OUT = ngc.OutputOptions(
    edge_flows=True, residual=True, saturated_boundary_edges=True
)
# A chain 0 -> 1 -> 2 -> 3 (cost 1 per edge) with the bypasses 0 -> 2 and 1 -> 3 (cost 3), every
# capacity 1: (u, v, capacity, cost) in input order.
CHAIN = [(0, 1, 1, 1), (1, 2, 1, 1), (2, 3, 1, 1), (0, 2, 1, 3), (1, 3, 1, 3)]


def exact(text):
    """The pytest.raises pattern of exactly this message: the whole text, so a message with
    anything more, a trailing newline included, does not match."""
    return rf"\A{re.escape(text)}\Z"


def assert_read_only(obj, *names):
    """Each name is a read-only property of obj: the class defines it without a setter, it reads
    (a misspelt or removed name fails here, not as the AttributeError below), and assigning the
    value it holds raises AttributeError and leaves it unchanged."""
    for name in names:
        prop = getattr(type(obj), name)
        assert isinstance(prop, property) and prop.fset is None, name
        value = getattr(obj, name)
        with pytest.raises(AttributeError):
            setattr(obj, name, value)
        assert getattr(obj, name) == value, name


def build(n, edges, q=1.0):
    """A Graph from (u, v, capacity, cost) rows in input order, at resolution q."""
    return ngc.Graph(
        n,
        [e[0] for e in edges],
        [e[1] for e in edges],
        [e[2] for e in edges],
        [e[3] for e in edges],
        q,
    )


def to_ids(g, values_in):
    """A per-edge list in input order, permuted to EdgeId order."""
    return [values_in[int(i)] for i in g.input_positions()]


def by_input(g, values):
    """A per-edge list in EdgeId order, permuted to input order (the inverse of to_ids)."""
    out = [0] * len(values)
    for e, x in enumerate(values):
        out[int(g.input_positions()[e])] = x
    return out


def stats_of(cache):
    """The five statistics of a ShortestPathCache as one comparable tuple."""
    s = cache.stats()
    return (s.hits, s.misses, s.entries, s.bytes, s.live_values)


def ids(g, positions):
    """The EdgeIds of the edges at the given input positions, ascending."""
    wanted = {int(p) for p in positions}
    return sorted(int(e) for e, p in enumerate(g.input_positions()) if int(p) in wanted)


def assert_ints(values, optional=False):
    """Every value is exactly an int: not a bool, not a float equal to it (None too when optional).
    Equality alone passes 3.0 for 3."""
    for v in values:
        assert type(v) is int or (optional and v is None), (v, type(v).__name__)


def dag_of(g, edges):
    """The PredecessorDAG of a set of EdgeIds, through the public PredecessorDAG.from_csr."""
    per = [[] for _ in range(g.num_nodes())]
    for e in sorted(edges):
        per[int(g.edge_destinations()[e])].append(e)
    offsets, parents, via = [0], [], []
    for v in range(g.num_nodes()):
        for e in per[v]:
            parents.append(int(g.edge_sources()[e]))
            via.append(e)
        offsets.append(len(via))
    return ngc.PredecessorDAG.from_csr(g, offsets, parents, via)


def shortest(model=None, **changes):
    """A ShortestPaths model: `model` (default ShortestPaths()) with the named fields replaced."""
    model = model or ngc.ShortestPaths()
    fields = {
        name: getattr(model, name)
        for name in (
            "allocation",
            "source_allocation",
            "weight_basis",
            "edge_weights",
            "tie_break",
            "multipath",
            "parallel_edges",
        )
    }
    fields.update(changes)
    return ngc.ShortestPaths(**fields)


def config(routing, allocation=None, admission=None):
    """A FlowPolicyConfig, Lossless by default; `allocation` (ShortestPaths only) replaces the
    model's own."""
    if allocation is not None:
        routing = shortest(routing, allocation=allocation)
    return ngc.FlowPolicyConfig(
        routing=routing, admission=admission or ngc.Admission.LOSSLESS
    )


def centralized(source_allocation=None, **fields):
    """CentralizedTE with its required source allocation: OPTIMIZED (the exact maximum flow) by
    default, EQUAL with a tunnel count (T2)."""
    if source_allocation is None:
        source_allocation = (
            ngc.SourceAllocation.EQUAL
            if fields.get("tunnel_count") is not None
            else ngc.SourceAllocation.OPTIMIZED
        )
    return ngc.CentralizedTE(source_allocation=source_allocation, **fields)


def t1_config(count, admission=None, **routing):
    """A T1 FlowPolicyConfig: DistributedTE(TunnelCount(count), **routing), Lossless by default."""
    return ngc.FlowPolicyConfig(
        routing=ngc.DistributedTE(ngc.TunnelCount(count), **routing),
        admission=admission or ngc.Admission.LOSSLESS,
    )


# The exact maximum flow (CentralizedTE routing, the optimized source allocation) and ShortestPaths Equal, both
# Lossless (ShortestPaths() is ECMP: FEWER_HOPS, every route and parallel edge). Configurations
# are immutable, so tests share these values.
GLOBAL = config(centralized())
COSTONLY = config(ngc.ShortestPaths())


# The result classes whose public attributes are all fields (no methods).
_RESULTS = (
    ngc.PlacementResult,
    ngc.SourceFlow,
    ngc.TunnelFlow,
    ngc.ForwardingProgram,
    ngc.Tunnel,
    ngc.TunnelPlan,
    ngc.SensitivityResult,
    ngc.SensitivityEntry,
    ngc.SpfResult,
)


def plain(value):
    """A result as comparable data, with no object identities: a result object (a program and its
    tunnels included) by every public field, a PredecessorDAG by its three arrays, a Path by its
    endpoints, hops and cost, an Eligibility by its two masks, containers element by element,
    arrays as lists."""
    if isinstance(value, getattr(_netgraph_core, "_PlacementTrace", ())):
        return plain(value.as_dict())
    if isinstance(value, _RESULTS):
        fields = [
            f
            for f in dir(value)
            if not f.startswith("_") and not callable(getattr(value, f))
        ]
        return {f: plain(getattr(value, f)) for f in fields}
    if isinstance(value, ngc.PredecessorDAG):
        return [
            plain(a)
            for a in (value.parent_offsets(), value.parents(), value.via_edges())
        ]
    if isinstance(value, ngc.Eligibility):
        return [plain(value.node_allowed()), plain(value.edge_allowed())]
    if isinstance(value, ngc.Path):
        return [value.source(), value.destination(), plain(value.hops()), value.cost()]
    if isinstance(value, dict):
        return {k: plain(x) for k, x in value.items()}
    if isinstance(value, (list, tuple)):
        return [plain(x) for x in value]
    return value.tolist() if hasattr(value, "tolist") else value


def ledger_state(ledger):
    """A FlowLedger's public state as comparable data, with no object identity (its uid, which
    is also the first half of every FlowId, and its graph's uid, which every route carries):
    both projections (the residual and the allocations, the holds and the unreserved capacity)
    and every public field of every record, the route as its three arrays."""
    return [
        ledger.residual().tolist(),
        ledger.carried().tolist(),
        ledger.reserved().tolist(),
        ledger.unreserved().tolist(),
        [
            (
                r.flow_id[1],
                r.policy_id,
                r.sources,
                r.destinations,
                r.admission,
                r.entry_id,
                [
                    r.route.parent_offsets().tolist(),
                    r.route.parents().tolist(),
                    r.route.via_edges().tolist(),
                ],
                r.min_route_cost,
                r.offered,
                r.delivered,
                r.edge_flows,
                r.edge_drops,
                r.edge_reservations,
            )
            for r in ledger.records()
        ],
    ]


def full_ledger_state(ledger):
    """The ledger's observable state, its next FlowId and bound and orphaned owners included
    (``_state()``: private hooks needed)."""
    return [ledger._state(), ledger_state(ledger)]


def sweep_python_allocation(make, limit=400, collect=True, allocation_call=None):
    """The all-or-nothing sweep of a call's Python object allocations (private hooks needed).
    make() gives fresh (fn, args, snapshot); attempt k = 0, 1, ... makes the k-th Python object
    allocation of fn(*args) fail, until a call completes. Every reached failure raises
    MemoryError (exactly that class) with no value and leaves snapshot() as it was, and the same call
    then returns the reference (the plain() value of a call on fresh state); the call that
    completes returns the reference after each of its allocations failed exactly once. Returns
    the number of failures, for the caller's check that the sweep is not vacuous. collect (the
    default): a full collection (gc.collect()) before each attempt empties CPython's free lists, so
    every int, float, list and tuple the call makes is a new allocation that the attempt can fail."""
    import gc

    import _netgraph_core as ext

    # Install native probes outside the checked allocation window, when needed.
    if allocation_call is None:
        allocation_call = ext._fail_python_allocation

    fn, args, _ = make()
    want = plain(fn(*args))
    failures = 0
    for k in range(limit):
        fn, args, snapshot = make()
        before = snapshot()
        if collect:
            gc.collect()
        hit, seen, value, error = allocation_call(fn, args, k)
        where = (getattr(fn, "__name__", fn), k)
        if not hit:
            assert error is None, (where, repr(error))
            assert plain(value) == want, where
            assert seen == failures, (where, seen)  # every allocation failed once
            return failures
        failures += 1
        # A reached failure propagates: a call that absorbed it and returned a value fails here.
        assert value is None, (where, value)
        assert type(error) is MemoryError, (where, repr(error))
        assert snapshot() == before, where
        assert plain(fn(*args)) == want, where
    pytest.fail(f"{getattr(fn, '__name__', fn)} never completed")
