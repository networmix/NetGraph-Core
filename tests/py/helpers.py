"""Builders and constants shared by the Python tests (not collected: no test_ prefix): graphs,
EdgeId mappings, DAGs, path selections, the two policy configurations most tests use, results
and ledgers as comparable data, the sweep of a call's Python object allocations, and the
pytest.raises pattern of an exact message."""

from __future__ import annotations

import math
import re

import pytest

import netgraph_core as ngc

INF = math.inf
L62 = 1 << 62
DET = ngc.EdgeTieBreak.DETERMINISTIC
PREFER = ngc.EdgeTieBreak.PREFER_HIGHER_PRIORITY
CAP = ngc.EdgePriority.CAPACITY
RES = ngc.EdgePriority.RESIDUAL
# Every optional output of a placement or an evaluator.
ALL_OUT = ngc.OutputOptions(edge_flows=True, residual=True, bottleneck=True)
# A chain 0 -> 1 -> 2 -> 3 (cost 1 per edge) with the bypasses 0 -> 2 and 1 -> 3 (cost 3), every
# capacity 1: (u, v, capacity, cost) in input order.
CHAIN = [(0, 1, 1, 1), (1, 2, 1, 1), (2, 3, 1, 1), (0, 2, 1, 3), (1, 3, 1, 3)]


def exact(text):
    """The pytest.raises pattern of exactly this message: the whole text, so a message with
    anything more, a trailing newline included, does not match."""
    return rf"\A{re.escape(text)}\Z"


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
    return [values_in[int(i)] for i in g.input_index()]


def by_input(g, values):
    """A per-edge list in EdgeId order, permuted to input order (the inverse of to_ids)."""
    out = [0] * len(values)
    for e, x in enumerate(values):
        out[int(g.input_index()[e])] = x
    return out


def stats_of(cache):
    """The five statistics of a DagCache as one comparable tuple."""
    s = cache.stats()
    return (s.hits, s.misses, s.entries, s.bytes, s.live_values)


def ids(g, positions):
    """The EdgeIds of the edges at the given input positions, ascending."""
    wanted = {int(p) for p in positions}
    return sorted(int(e) for e, p in enumerate(g.input_index()) if int(p) in wanted)


def dag_of(g, edges):
    """The PredDAG of a set of EdgeIds, through the public PredDAG.from_csr."""
    per = [[] for _ in range(g.num_nodes())]
    for e in sorted(edges):
        per[int(g.edge_dst()[e])].append(e)
    offsets, parents, via = [0], [], []
    for v in range(g.num_nodes()):
        for e in per[v]:
            parents.append(int(g.edge_src()[e]))
            via.append(e)
        offsets.append(len(via))
    return ngc.PredDAG.from_csr(g, offsets, parents, via)


def sel(tie=DET, prio=CAP, mp=True, me=True):
    """A PathSelection: deterministic ties, capacity priority, multipath, multi-edge by default."""
    return ngc.PathSelection(tie, prio, mp, me)


def config(
    routing, placement=ngc.FlowPlacement.EQUAL_BALANCED, admission=None, selection=None
):
    """A FlowPolicyConfig: EqualBalanced, Lossless and sel() by default; Global takes neither a
    placement nor a selection, so they are left out for it."""
    if isinstance(routing, ngc.Global):
        return ngc.FlowPolicyConfig(
            routing=routing, admission=admission or ngc.Admission.LOSSLESS
        )
    return ngc.FlowPolicyConfig(
        routing=routing,
        placement=placement,
        admission=admission or ngc.Admission.LOSSLESS,
        selection=selection or sel(),
    )


# The exact maximum flow (Global routing, the optimized origin) and CostOnly EqualBalanced, both
# Lossless (CostOnly with sel()). Configurations are immutable, so tests share these values.
GLOBAL = config(ngc.Global())
COSTONLY = config(ngc.CostOnly())


# The result classes whose public attributes are all fields (no methods).
_RESULTS = (
    ngc.PlacementResult,
    ngc.SourceFlow,
    ngc.SensitivityResult,
    ngc.SensitivityEntry,
    ngc.SpfResult,
)


def plain(value):
    """A result as comparable data, with no object identities: a result object by every public
    field, a PredDAG by its three arrays, a Path by its endpoints, hops and cost, containers
    element by element, arrays as lists."""
    if isinstance(value, _RESULTS):
        fields = [f for f in dir(value) if not f.startswith("_")]
        return {f: plain(getattr(value, f)) for f in fields}
    if isinstance(value, ngc.PredDAG):
        return [
            plain(a)
            for a in (value.parent_offsets(), value.parents(), value.via_edges())
        ]
    if isinstance(value, ngc.Path):
        return [value.src(), value.dst(), plain(value.hops()), value.cost()]
    if isinstance(value, dict):
        return {k: plain(x) for k, x in value.items()}
    if isinstance(value, (list, tuple)):
        return [plain(x) for x in value]
    return value.tolist() if hasattr(value, "tolist") else value


def ledger_state(ledger):
    """A FlowLedger's public state as comparable data, with no object identity (its uid, which
    is also the first half of every FlowId, and its graph's uid, which every bundle carries):
    the residual, the allocations and every public field of every record, the bundle as its
    three arrays."""
    return [
        ledger.residual().tolist(),
        ledger.allocated().tolist(),
        [
            (
                r.flow_id[1],
                r.owner,
                r.sources,
                r.destinations,
                r.admission,
                [
                    r.bundle.parent_offsets().tolist(),
                    r.bundle.parents().tolist(),
                    r.bundle.via_edges().tolist(),
                ],
                r.route_cost,
                r.offered,
                r.delivered,
                r.edges,
                r.drops,
            )
            for r in ledger.records()
        ],
    ]


def full_ledger_state(ledger):
    """The ledger's observable state, its next FlowId and bound and orphaned owners included
    (``_state()``: private hooks needed)."""
    return [ledger._state(), ledger_state(ledger)]


def sweep_python_allocation(make, limit=400):
    """The all-or-nothing sweep of a call's Python object allocations (private hooks needed).
    make() gives fresh (fn, args, snapshot); attempt k = 0, 1, ... makes the k-th Python object
    allocation of fn(*args) fail, until a call completes. Every reached failure raises
    MemoryError or RuntimeError with no value and leaves snapshot() as it was, and the same call
    then returns the reference (the plain() value of a call on fresh state); the call that
    completes returns the reference after each of its allocations failed exactly once. Returns
    the number of failures, for the caller's check that the sweep is not vacuous."""
    import _netgraph_core as ext

    fn, args, _ = make()
    want = plain(fn(*args))
    failures = 0
    for k in range(limit):
        fn, args, snapshot = make()
        before = snapshot()
        hit, seen, value, error = ext._fail_python_allocation(fn, args, k)
        where = (getattr(fn, "__name__", fn), k)
        if not hit:
            assert error is None, (where, repr(error))
            assert plain(value) == want, where
            assert seen == failures, (where, seen)  # every allocation failed once
            return failures
        failures += 1
        # A reached failure propagates: a call that absorbed it and returned a value fails here.
        assert value is None, (where, value)
        assert isinstance(error, (MemoryError, RuntimeError)), (where, repr(error))
        assert snapshot() == before, where
        assert plain(fn(*args)) == want, where
    pytest.fail(f"{getattr(fn, '__name__', fn)} never completed")
