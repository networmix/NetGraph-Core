"""The origin and the result of one demand (docs/semantics.md, sections 7 and 10): the
configuration's defaults and table with its error texts and their order, the pooled single-path
rule, the controllers (Global and CostOnly MAX_FLOW) under both origins, min_cut as its own
output, zero calls before and after preparation, the result's accounting and per-source rows,
the terminations and the evaluators. The values are derived by hand (dev/design step-2 LEDGER
rows C, Q, M, G, S, T, I5, I6 and V); the random tests compare the controllers with the
independent equal-origin and optimized references of ``origin_oracle.py``."""

from __future__ import annotations

import random

import pytest
from helpers import build, dag_of, exact, sel
from origin_oracle import (
    INF,
    PooledReference,
    distances_to,
    equal_origin,
    optimized_origin,
    rr,
)

import netgraph_core as ngc

T = ngc.Termination
P = ngc.FlowPlacement
A = ngc.Admission
ORIGIN = ngc.Origin
ALL = ngc.OutputOptions(edge_flows=True, residual=True, bottleneck=True, min_cut=True)
CUT = ngc.OutputOptions(min_cut=True)

POOL = [(0, 2, 10, 1), (1, 2, 10, 1), (2, 3, 1, 1), (2, 3, 1, 1)]  # {0, 1} -> 3
B1 = [(0, 3, 2, 1), (1, 2, 10, 1), (2, 3, 10, 2)]  # {0, 1} -> 3
B2 = [(0, 1, 4, 1), (0, 2, 2, 1), (1, 2, 10, 1)]  # {0, 1} -> 2
B2S = [(0, 2, 2, 1), (1, 2, 10, 1)]  # {0, 1} -> 2
B3Z = [(0, 2, 0, 1), (1, 2, 10, 1)]  # {0, 1} -> 2
B3P = [(0, 2, 10, 1), (1, 2, 0, 1)]  # {0, 1} -> 2
ES = [(0, 1, 10, 1), (0, 2, 10, 1), (1, 3, 2, 1), (2, 3, 10, 2)]  # 0 -> 3
R4 = [(0, 2, 2, 9), (1, 2, 10, 1)]  # {0, 1} -> 2
ONE = [(0, 1, 10, 1)]
PAR = [(0, 1, 5, 1), (0, 1, 8, 1)]
ZERO = [(0, 1, 0, 1)]
SPF = [(0, 1, 0, 1), (0, 2, 5, 1), (2, 1, 5, 1)]  # 0 -> 1
ROWS = [(0, 3, 5, 1), (1, 2, 5, 1), (4, 3, 5, 1)]  # {0, 1, 4} -> 3, node 4 masked
WIT = [(0, 1, 1, 1), (0, 2, 1, 1), (1, 0, 1, 1)]  # 0 -> 2
PRIV = [(0, 2, 4, 1), (1, 2, 6, 1)]  # {0, 1} -> 2

Q12 = "CostOnly with multiple input sources rejects multipath=false with PREFER_HIGHER_PRIORITY"
PREFER_PATH = sel(tie=ngc.EdgeTieBreak.PREFER_HIGHER_PRIORITY, mp=False)


def pinned_route(g, src, dst, edges, placement=None, admission=None, **fields):
    """Pinned over a program of one split route (the edges' DAG) from src to dst."""
    program = ngc.ForwardingProgram.from_routes(
        g, src, dst, [ngc.Route(src, dag_of(g, edges))], placement=placement
    )
    return ngc.FlowPolicyConfig(
        routing=ngc.Pinned(program),
        admission=admission or ngc.Admission.LOSSLESS,
        **fields,
    )


def graph(edges):
    return build(1 + max(max(u, v) for u, v, *_ in edges), edges)


def glob(origin=None):
    return ngc.FlowPolicyConfig(routing=ngc.Global(), origin=origin)


def cost_only(
    placement=P.EQUAL_BALANCED, admission=A.LOSSLESS, selection=None, origin=None
):
    return ngc.FlowPolicyConfig(
        routing=ngc.CostOnly(),
        placement=placement,
        admission=admission,
        selection=selection,
        origin=origin,
    )


def rows(r):
    return [
        (x.source, x.offered, x.admitted, x.delivered, x.source_loss, x.dropped)
        for x in r.by_source
    ]


def fresh(edges, sources, dests, config, elig=None, residual=None):
    g = graph(edges)
    ledger = ngc.FlowLedger(g, residual)
    e = elig(g) if elig is not None else None
    return g, ledger, ngc.FlowPolicy(ledger, sources, dests, config, e)


# -- configuration ------------------------------------------------------------------------


def test_defaults_resolve():
    # C1, C2: the defaults are None on the value and resolve at the policy.
    c = ngc.FlowPolicyConfig(routing=ngc.CostOnly(), placement=P.EQUAL_BALANCED)
    assert (c.origin, c.admission, c.selection) == (None, A.LOSSLESS, None)
    g = graph(B1)
    p = ngc.FlowPolicy(ngc.FlowLedger(g), [0, 1], 3, c)
    assert (p.config.origin, p.config.selection) == (ORIGIN.EQUAL, ngc.PathSelection())
    for cfg in (cost_only(P.MAX_FLOW), glob()):
        assert (
            ngc.FlowPolicy(ngc.FlowLedger(g), [0, 1], 3, cfg).config.origin
            == ORIGIN.OPTIMIZED
        )
    for cfg in (cost_only(P.MAX_FLOW, origin=ORIGIN.EQUAL), glob(ORIGIN.EQUAL)):
        assert (
            ngc.FlowPolicy(ngc.FlowLedger(g), [0, 1], 3, cfg).config.origin
            == ORIGIN.EQUAL
        )
    assert ngc.FlowPolicy(ngc.FlowLedger(g), [0, 1], 3, glob()).config.selection is None
    pinned = pinned_route(g, 0, 3, [0])
    assert ngc.FlowPolicy(ngc.FlowLedger(g), 0, 3, pinned).config.origin is None


def test_conversion_errors():
    # C3: conversion first: positional arguments, an unknown keyword, a string origin.
    with pytest.raises(TypeError):
        ngc.FlowPolicyConfig(ngc.CostOnly(), P.EQUAL_BALANCED)
    with pytest.raises(TypeError):
        ngc.FlowPolicyConfig(routing=ngc.CostOnly(), bundles=None)
    with pytest.raises(TypeError):
        ngc.FlowPolicyConfig(routing=ngc.Global(), origin="EQUAL")


def test_path_selection_defaults():
    # C17.
    assert ngc.PathSelection() == ngc.PathSelection(
        ngc.EdgeTieBreak.DETERMINISTIC, ngc.EdgePriority.CAPACITY, True, True
    )


@pytest.mark.parametrize(
    ("config", "message"),
    [
        # C4 to C12.
        (ngc.FlowPolicyConfig(routing=ngc.CostOnly()), "CostOnly requires placement"),
        (
            ngc.FlowPolicyConfig(routing=ngc.Global(), placement=P.MAX_FLOW),
            "Global does not accept placement",
        ),
        (
            ngc.FlowPolicyConfig(routing=ngc.Global(), selection=ngc.PathSelection()),
            "Global does not accept selection",
        ),
        (
            cost_only(P.EQUAL_BALANCED, origin=ORIGIN.OPTIMIZED),
            "CostOnly with EQUAL_BALANCED or WEIGHTED requires origin=EQUAL",
        ),
        (
            cost_only(P.WEIGHTED, origin=ORIGIN.OPTIMIZED),
            "CostOnly with EQUAL_BALANCED or WEIGHTED requires origin=EQUAL",
        ),
        (cost_only(P.MAX_FLOW, A.LOSSY), "MAX_FLOW requires admission=LOSSLESS"),
        (
            ngc.FlowPolicyConfig(routing=ngc.Global(), admission=A.LOSSY),
            "MAX_FLOW requires admission=LOSSLESS",
        ),
        (
            ngc.FlowPolicyConfig(
                routing=ngc.CapacityAware(),
                placement=P.MAX_FLOW,
                origin=ORIGIN.OPTIMIZED,
            ),
            "T1 requires origin=EQUAL",
        ),
        (
            ngc.FlowPolicyConfig(routing=ngc.CapacityAware()),
            "CapacityAware requires placement",
        ),
        (
            ngc.FlowPolicyConfig(
                routing=ngc.CapacityAware(),
                placement=P.EQUAL_BALANCED,
                admission=A.LOSSY,
            ),
            "CapacityAware requires admission=LOSSLESS",
        ),
        # C13, C14, Q6: the row, then the columns placement, origin, admission, selection.
        (
            ngc.FlowPolicyConfig(
                routing=ngc.Global(),
                placement=P.MAX_FLOW,
                admission=A.LOSSY,
                selection=ngc.PathSelection(),
            ),
            "Global does not accept placement",
        ),
        (
            ngc.FlowPolicyConfig(
                routing=ngc.Global(),
                origin=ORIGIN.EQUAL,
                admission=A.LOSSY,
                selection=ngc.PathSelection(),
            ),
            "MAX_FLOW requires admission=LOSSLESS",
        ),
        (
            cost_only(P.EQUAL_BALANCED, selection=PREFER_PATH, origin=ORIGIN.OPTIMIZED),
            "CostOnly with EQUAL_BALANCED or WEIGHTED requires origin=EQUAL",
        ),
    ],
)
def test_matrix_texts(config, message):
    g = graph(B1)
    with pytest.raises(ValueError, match=exact(message)):
        ngc.FlowPolicy(ngc.FlowLedger(g), [0, 1], 3, config)
    with pytest.raises(ValueError, match=exact(message)):
        ngc.capacity(g, [0, 1], 3, ngc.QueryOptions(config))


def test_pinned_rows():
    # C8, C10, C12 for Pinned: no placement and no origin (its program has the injection).
    g = graph(B1)
    for config, message in (
        (
            pinned_route(g, 0, 3, [0], origin=ORIGIN.EQUAL),
            "Pinned does not accept origin",
        ),
        (
            ngc.FlowPolicyConfig(
                routing=pinned_route(g, 0, 3, [0]).routing, placement=P.MAX_FLOW
            ),
            "Pinned does not accept placement",
        ),
    ):
        with pytest.raises(ValueError, match=exact(message)):
            ngc.FlowPolicy(ngc.FlowLedger(g), 0, 3, config)


def test_cache_rule():
    # C11, C16: only CostOnly reads a DagCache; the table comes first.
    g = graph(B1)
    cache = ngc.DagCache(g)
    pinned = pinned_route(g, 0, 3, [0])
    capacity_aware = ngc.FlowPolicyConfig(
        routing=ngc.CapacityAware(), placement=P.MAX_FLOW
    )
    for config, sources, message in (
        (glob(), [0, 1], "Global does not use DagCache"),
        (capacity_aware, [0, 1], "T1 does not use DagCache"),
        (pinned, 0, "Pinned does not use DagCache"),
        (
            ngc.FlowPolicyConfig(routing=ngc.Global(), placement=P.MAX_FLOW),
            [0, 1],
            "Global does not accept placement",
        ),
    ):
        with pytest.raises(ValueError, match=exact(message)):
            ngc.FlowPolicy(ngc.FlowLedger(g), sources, 3, config, None, cache)
        with pytest.raises(ValueError, match=exact(message)):
            ngc.capacity(g, sources, 3, ngc.QueryOptions(config, cache=cache))
        with pytest.raises(ValueError, match=exact(message)):
            ngc.sensitivity(g, sources, 3, ngc.QueryOptions(config, cache=cache))
    assert cache.stats().misses == 0 and cache.stats().hits == 0


def test_endpoints_before_ownership():
    # C15.
    g = graph(B1)
    foreign = graph(B1).eligibility(excluded_nodes=[0])
    with pytest.raises(ValueError, match="sources") as e:
        ngc.FlowPolicy(ngc.FlowLedger(g), [], 3, glob(), foreign)
    assert "foreign" not in str(e.value)


def test_amount_texts():
    # C18, C19.
    _, _, p = fresh(ONE, 0, 1, cost_only())
    for bad in (-1, 2**62, float("nan"), 1.5, -INF, True):
        with pytest.raises(
            ValueError,
            match=exact("offer must be an integer in [0, 2^62) or permitted +inf"),
        ):
            p.place(bad)
    p.place(INF)
    for finite in (1, 0):
        with pytest.raises(
            ValueError,
            match=exact(
                "finite placement after infinite placement requires a new epoch"
            ),
        ):
            p.place(finite)
    p.remove()
    assert p.place(1).placed == 1


# -- the pooled single-path rule (Q12) ------------------------------------------------------


def test_q12():
    g = graph(POOL)
    for placement in (P.EQUAL_BALANCED, P.WEIGHTED, P.MAX_FLOW):
        config = cost_only(placement, selection=PREFER_PATH)
        # Q1, Q3: two input sources, also when the mask leaves one.
        for elig in (None, g.eligibility(excluded_nodes=[1])):
            with pytest.raises(ValueError, match=exact(Q12)):
                ngc.FlowPolicy(ngc.FlowLedger(g), [0, 1], 3, config, elig)
            with pytest.raises(ValueError, match=exact(Q12)):
                ngc.capacity(g, [0, 1], 3, ngc.QueryOptions(config, elig=elig))
        # Q2: one source; Q4: multipath.
        assert ngc.FlowPolicy(ngc.FlowLedger(g), [0, 0], 3, config).place(1).placed == 1
        multipath = cost_only(
            placement, selection=sel(tie=ngc.EdgeTieBreak.PREFER_HIGHER_PRIORITY)
        )
        assert (
            ngc.FlowPolicy(ngc.FlowLedger(g), [0, 1], 3, multipath).place(2).placed == 2
        )
    # Q7: batch_capacity checks each pair; pair 1 is the lowest failing.
    opts = ngc.QueryOptions(cost_only(selection=PREFER_PATH))
    with pytest.raises(ValueError, match=exact(Q12)):
        ngc.batch_capacity(g, [([0], 3), ([0, 1], 3)], opts)
    assert [r.placed for r in ngc.batch_capacity(g, [([0], 3), ([1], 3)], opts)] == [
        1,
        1,
    ]


def test_deterministic_single_paths_pool():
    # Q5: canonical paths 0: [e0, e2], 1: [e1, e2] (equal cost and hops, lower EdgeId e2).
    det = sel(mp=False)
    _, _, p = fresh(POOL, [0, 1], 3, cost_only(selection=det))
    r = p.place(INF)
    assert (r.placed, [x.admitted for x in r.by_source]) == (1, [1, 0])
    # Lossy 2: e2's arrivals (1, 1), carried H(1, [1, 1]) = [1, 0].
    _, _, p = fresh(POOL, [0, 1], 3, cost_only(admission=A.LOSSY, selection=det))
    r = p.place(2)
    assert (r.placed, r.call_drops) == (1, [(2, 1)])
    assert rows(r) == [(0, 1, 1, 1, 0, 0), (1, 1, 1, 0, 0, 1)]
    assert [list(rec.bundle.via_edges()) for rec in p.records()] == [[0, 2], [1, 2]]


# -- the controllers --------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("edges", "sources", "dest", "origin", "placed", "cost", "flows"),
    [
        (B1, [0, 1], 3, ORIGIN.EQUAL, 4, 8, [2, 2, 2]),  # G1 B1-equal
        (B1, [0, 1], 3, ORIGIN.OPTIMIZED, 12, 32, [2, 10, 10]),  # G2 B1-optimized
        (B2, [0, 1], 2, ORIGIN.EQUAL, 12, 16, [4, 2, 10]),  # G3 B2-equal-global
        (B2, [0, 1], 2, ORIGIN.OPTIMIZED, 12, 12, [0, 2, 10]),  # G4 B2-optimized
        (B2S, [0, 1], 2, ORIGIN.EQUAL, 4, 4, [2, 2]),  # G5 B2-shortest
        (B2S, [0, 1], 2, ORIGIN.OPTIMIZED, 12, 12, [2, 10]),  # G5 T2-OPT-global
        (B3Z, [0, 1], 2, ORIGIN.EQUAL, 0, 0, [0, 0]),  # G6 B3-zero-first
        (B3P, [0, 1], 2, ORIGIN.EQUAL, 1, 1, [1, 0]),  # G6 B3-positive-first
        (B3Z, [0, 1], 2, ORIGIN.OPTIMIZED, 10, 10, [0, 10]),  # G6
        (ES, 0, 3, ORIGIN.OPTIMIZED, 12, 34, [2, 10, 2, 10]),  # G7 E-support
    ],
)
def test_global_examples(edges, sources, dest, origin, placed, cost, flows):
    request = 10 if edges in (B3Z, B3P) else 12
    _, _, p = fresh(edges, sources, dest, glob(origin))
    r = p.place(request, ngc.OutputOptions(edge_flows=True))
    term = T.COMPLETE if placed == request else T.EXHAUSTED
    assert (r.placed, r.remaining, r.termination) == (placed, request - placed, term)
    assert (r.total_cost, r.edge_flows.tolist()) == (cost, flows)
    if origin == ORIGIN.EQUAL:
        assert r.cost_distribution is None
        assert [x.admitted for x in r.by_source] == rr(placed, len(r.by_source))
    else:
        assert sum(u for _, u in r.cost_distribution) == placed


def test_global_records_and_masks():
    # G1 rows; G3 records: source 0's 6 is 2 direct and 4 through node 1.
    _, _, p = fresh(B1, [0, 1], 3, glob(ORIGIN.EQUAL))
    assert rows(p.place(12)) == [(0, None, 2, 2, 0, 0), (1, None, 2, 2, 0, 0)]
    _, _, p = fresh(B2, [0, 1], 2, glob(ORIGIN.EQUAL))
    p.place(12)
    assert [(rec.edges, rec.delivered) for rec in p.records()] == [
        ([(0, 4), (1, 2), (2, 4)], 6),
        ([(2, 6)], 6),
    ]
    # G6 B3-delete-zero: the masked source leaves the roster.
    _, _, p = fresh(
        B3Z, [0, 1], 2, glob(ORIGIN.EQUAL), lambda g: g.eligibility(excluded_nodes=[0])
    )
    assert p.place(10).placed == 10
    # G11 R3-residual-optimal.
    _, _, p = fresh(WIT, 0, 2, glob())
    r = p.place(1, ngc.OutputOptions(edge_flows=True))
    assert (r.edge_flows.tolist(), r.total_cost) == ([0, 1, 0], 1)


def test_global_increments():
    # G8 R4-lossless-equal: 3 (RR(3) = (2, 1) fits), then the target 6 admits 4 (RR(5) puts
    # 3 on e0 of 2): allocation (2, 2), cost 2 * 9 + 2 * 1.
    _, _, p = fresh(R4, [0, 1], 2, glob(ORIGIN.EQUAL))
    a = p.place(3)
    b = p.place(3, ngc.OutputOptions(edge_flows=True))
    assert (a.placed, a.termination, b.placed, b.termination) == (
        3,
        T.COMPLETE,
        1,
        T.EXHAUSTED,
    )
    assert (b.edge_flows.tolist(), b.total_cost, b.requested, b.unadmitted) == (
        [2, 2],
        20,
        6,
        2,
    )
    # G9 R4-refusal: refused units are not offered again.
    g = graph(ONE)
    ledger = ngc.FlowLedger(g)
    other = ngc.FlowPolicy(ledger, 0, 1, cost_only())
    other.place(8)
    p = ngc.FlowPolicy(ledger, 0, 1, glob())
    assert p.place(5).placed == 2
    other.remove()
    other.place(5)
    assert p.place(1).placed == 1
    # G10: equal increments end where one call does (the final total solved from zero).
    _, _, one = fresh(B2, [0, 1], 2, glob(ORIGIN.EQUAL))
    whole = one.place(12, ngc.OutputOptions(edge_flows=True))
    _, _, two = fresh(B2, [0, 1], 2, glob(ORIGIN.EQUAL))
    two.place(6)
    split = two.place(6, ngc.OutputOptions(edge_flows=True))
    assert split.edge_flows.tolist() == whole.edge_flows.tolist()
    assert split.total_cost == whole.total_cost
    assert [(r.edges, r.delivered) for r in two.records()] == [
        (r.edges, r.delivered) for r in one.records()
    ]


@pytest.mark.parametrize(
    ("edges", "sources", "dest", "origin", "placed", "cost", "records"),
    [
        # M1, M2: one controller on the union of the sources' routes.
        (B1, [0, 1], 3, ORIGIN.OPTIMIZED, 12, 32, [[(0, 2)], [(1, 10), (2, 10)]]),
        (B1, [0, 1], 3, ORIGIN.EQUAL, 4, 8, [[(0, 2)], [(1, 2), (2, 2)]]),
        # M3: the support leaves out e0 (0 -> 1 -> 2 costs 2).
        (B2, [0, 1], 2, ORIGIN.EQUAL, 4, 4, [[(1, 2)], [(2, 2)]]),
        (B2, [0, 1], 2, ORIGIN.OPTIMIZED, 12, 12, [[(1, 2)], [(2, 10)]]),
        # M4 E-support-shortest; M5: the support is the zero edge alone, never bypassed.
        (ES, 0, 3, ORIGIN.OPTIMIZED, 2, 4, [[(0, 2), (2, 2)]]),
        (SPF, 0, 1, ORIGIN.OPTIMIZED, 0, 0, [[]]),
    ],
)
def test_cost_only_max_flow(edges, sources, dest, origin, placed, cost, records):
    request = 5 if edges is SPF else 12
    _, _, p = fresh(edges, sources, dest, cost_only(P.MAX_FLOW, origin=origin))
    r = p.place(request, CUT)
    term = T.COMPLETE if placed == request else T.EXHAUSTED
    assert (r.placed, r.termination, r.total_cost) == (placed, term, cost)
    assert [rec.edges for rec in p.records()] == records
    # M7: no curve and no cut.
    assert (r.cost_distribution, r.min_cut) == (None, None)


def test_cost_only_max_flow_increments():
    # M6.
    _, _, p = fresh(R4, [0, 1], 2, cost_only(P.MAX_FLOW, origin=ORIGIN.EQUAL))
    assert [p.place(3).placed, p.place(3).placed] == [3, 1]
    assert p.place(0).total_cost == 20


def test_min_cut_flag():
    # G13 MINCUT-B1 and MINCUT-zero-edge.
    for edges, dest in ((B1, 3), (B3Z, 2)):
        _, _, p = fresh(edges, [0, 1], dest, glob())
        assert p.place(INF, CUT).min_cut == [0, 1]
    # G12: a zero maximum over existing routes is EXHAUSTED; its cut is the zero edge.
    _, _, p = fresh(ZERO, 0, 1, glob())
    r = p.place(5, CUT)
    assert (r.placed, r.termination, r.min_cut) == (0, T.EXHAUSTED, [0])
    # G14: bottleneck alone, the equal origin, a request below the maximum, a zero call.
    _, _, p = fresh(B1, [0, 1], 3, glob())
    assert p.place(INF, ngc.OutputOptions(bottleneck=True)).min_cut is None
    _, _, p = fresh(B1, [0, 1], 3, glob(ORIGIN.EQUAL))
    assert p.place(INF, CUT).min_cut is None
    _, _, p = fresh(B1, [0, 1], 3, glob())
    assert p.place(5, CUT).min_cut is None
    assert p.place(0, CUT).min_cut is None
    # G15: an empty roster.
    _, _, p = fresh([(1, 0, 5, 1)], 0, 1, glob())
    r = p.place(5, CUT)
    assert (r.termination, r.min_cut, r.cost_distribution, r.total_cost) == (
        T.NO_ROUTE,
        [],
        [],
        0,
    )


P40, P39, P41 = 2**40, 2**39, 2**41
LO32, HI32 = 2**32 - 1, 2**32 + 1


@pytest.mark.parametrize(
    ("edges", "sources", "dest", "optimized", "equal"),
    [
        # One link of capacity = cost = 2^40.
        ([(0, 1, P40, P40)], 0, 1, 2**80, 2**80),
        # Optimized fills both links: 2^80 + 2^80; equal stops at T = 2^40 + 1 (shares
        # 2^39 + 1 and 2^39): (2^39 + 1) * 2^40 + 2^39 * 2^41.
        (
            [(0, 2, P40, P40), (1, 2, P39, P41)],
            [0, 1],
            2,
            2**81,
            3 * 2**79 + 2**40,
        ),
        # Two products of 2^64 - 1: the sum carries into the high half.
        ([(0, 1, LO32, HI32), (1, 2, LO32, HI32)], 0, 2, 2**65 - 2, 2**65 - 2),
    ],
    ids=["product", "two-sources", "carry"],
)
def test_total_cost_beyond_64_bits(edges, sources, dest, optimized, equal):
    # The controllers' total_cost is exact beyond 64 bits, under both origins, through place,
    # capacity and batch_capacity with workers.
    g = graph(edges)
    for origin, want in ((ORIGIN.OPTIMIZED, optimized), (ORIGIN.EQUAL, equal)):
        assert want >= 2**64
        for config in (glob(origin), cost_only(P.MAX_FLOW, origin=origin)):
            p = ngc.FlowPolicy(ngc.FlowLedger(g), sources, dest, config)
            assert p.place(INF).total_cost == want
            assert (
                ngc.capacity(g, sources, dest, ngc.QueryOptions(config)).total_cost
                == want
            )
            batch = ngc.batch_capacity(
                g, [(sources, dest)] * 2, ngc.QueryOptions(config), workers=2
            )
            assert [r.total_cost for r in batch] == [want, want]


# -- zero calls, the result, terminations ----------------------------------------------------


@pytest.mark.parametrize("routing", ["costonly", "global"])
def test_zero_unprepared(routing):
    # I5 RC02-zero-unprepared; I6: after a positive call and after remove().
    config = cost_only() if routing == "costonly" else glob()
    g, ledger, p = fresh(PAR, 0, 1, config, residual=[4, 6])

    def unprepared(r):
        assert (r.placed, r.remaining, r.delivery_delta, r.termination) == (
            0,
            0,
            0,
            T.COMPLETE,
        )
        assert (
            r.call_drops,
            r.call_source_loss,
            r.routes_changed,
            r.weights_changed,
        ) == (
            [],
            0,
            False,
            False,
        )
        assert (r.edge_flows.tolist(), r.residual.tolist()) == ([0, 0], [4, 6])
        assert (r.min_cut, r.cost_distribution, r.requested) == (None, None, 0)
        assert rows(r) == [(0, None, 0, 0, 0, 0)]
        assert r.total_cost == (0 if routing == "global" else None)
        assert p.records() == []

    unprepared(p.place(0, ALL))
    p.place(3)
    r = p.place(0, ALL)
    assert (r.placed, r.routes_changed, r.delivered, r.requested) == (0, False, 3, 3)
    if routing == "global":
        assert r.cost_distribution == [(1, 3)] and r.min_cut is None
    p.remove()
    unprepared(p.place(0, ALL))


def test_result_accounting():
    # S2: G1 Lossless.
    _, _, p = fresh(B1, [0, 1], 3, glob(ORIGIN.EQUAL))
    r = p.place(12)
    assert (r.requested, r.offered, r.admitted, r.delivered, r.unadmitted) == (
        12,
        None,
        4,
        4,
        8,
    )
    assert (
        r.source_loss,
        r.drops,
        r.call_drops,
        r.call_source_loss,
        r.delivery_delta,
    ) == (
        0,
        [],
        [],
        0,
        4,
    )
    # S5: a row for every input source: 4 is masked, 1 has no route.
    _, _, p = fresh(
        ROWS, [0, 1, 4], 3, cost_only(), lambda g: g.eligibility(excluded_nodes=[4])
    )
    assert rows(p.place(3)) == [
        (0, None, 3, 3, 0, 0),
        (1, None, 0, 0, 0, 0),
        (4, None, 0, 0, 0, 0),
    ]
    # S6: total_cost only for the controllers.
    g = graph(B1)
    for config in (
        cost_only(),
        cost_only(P.WEIGHTED),
        ngc.FlowPolicyConfig(routing=ngc.CapacityAware(), placement=P.EQUAL_BALANCED),
    ):
        assert (
            ngc.FlowPolicy(ngc.FlowLedger(g), [0, 1], 3, config).place(1).total_cost
            is None
        )
    pinned = pinned_route(g, 0, 3, [0])
    assert ngc.FlowPolicy(ngc.FlowLedger(g), 0, 3, pinned).place(1).total_cost is None
    support = ngc.FlowPolicyConfig(routing=ngc.Global(support=dag_of(g, [0])))
    assert ngc.FlowPolicy(ngc.FlowLedger(g), 0, 3, support).place(1).total_cost == 1


def test_capacity_aware_and_pinned_in_the_new_shape():
    # S7: CapacityAware keeps its per-source computation: offers 6 and 6.
    config = ngc.FlowPolicyConfig(
        routing=ngc.CapacityAware(), placement=P.EQUAL_BALANCED
    )
    _, _, p = fresh(B1, [0, 1], 3, config)
    r = p.place(12)
    assert (r.placed, r.termination, r.requested, r.unadmitted) == (
        8,
        T.EXHAUSTED,
        12,
        4,
    )
    assert rows(r) == [(0, None, 2, 2, 0, 0), (1, None, 6, 6, 0, 0)]
    # S8: Weighted over the zero edge has no member: its offer is lost at the source.
    g = graph(ZERO)
    for placement, want in (
        (P.WEIGHTED, (0, T.NO_ACTIVE_WEIGHT, 2, [], (0, 2, 2, 0, 2, 0))),
        (P.EQUAL_BALANCED, (0, T.COMPLETE, 0, [(0, 2)], (0, 2, 2, 0, 0, 2))),  # S9
    ):
        config = pinned_route(g, 0, 1, [0], placement=placement, admission=A.LOSSY)
        r = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 1, config).place(2)
        assert (
            r.placed,
            r.termination,
            r.call_source_loss,
            r.call_drops,
            rows(r)[0],
        ) == want
        assert (r.offered, r.source_loss) == (2, want[2])


def test_terminations():
    # T1, T2: an empty roster (no edge reaches 1 from 0).
    for config in (glob(), cost_only(), cost_only(admission=A.LOSSY)):
        _, _, p = fresh([(1, 0, 5, 1)], 0, 1, config)
        assert p.place(0).termination == T.COMPLETE
        r = p.place(5)
        assert (r.placed, r.termination) == (0, T.NO_ROUTE)
        if config.admission == A.LOSSY:
            # O2: nothing enters: the offer is all source loss and the row stays zero.
            assert (r.offered, r.source_loss, r.call_source_loss) == (5, 5, 5)
            assert rows(r) == [(0, 0, 0, 0, 0, 0)]
        else:
            assert (r.offered, r.unadmitted) == (None, 5)
    # T3: weighted splitting over the zero edge: no positive weight.
    _, _, p = fresh(ZERO, 0, 1, cost_only(P.WEIGHTED, A.LOSSY))
    r = p.place(2)
    assert (r.termination, r.call_source_loss, r.offered) == (T.NO_ACTIVE_WEIGHT, 2, 2)
    _, _, p = fresh(ZERO, 0, 1, cost_only(P.WEIGHTED))
    r = p.place(2)
    assert (r.termination, r.unadmitted) == (T.NO_ACTIVE_WEIGHT, 2)
    # T6: a controller at zero over existing routes.
    for config in (glob(), cost_only(P.MAX_FLOW)):
        _, _, p = fresh(
            SPF if config.routing.__class__ is ngc.CostOnly else ZERO, 0, 1, config
        )
        assert p.place(5).termination == T.EXHAUSTED


def test_result_objects():
    # S11.
    _, _, p = fresh(B1, [0, 1], 3, glob(ORIGIN.EQUAL))
    r = p.place(INF)
    assert isinstance(r.by_source, tuple) and all(
        isinstance(x, ngc.SourceFlow) for x in r.by_source
    )
    assert (r.requested, r.remaining) == (INF, INF)
    assert isinstance(r.requested, float)
    with pytest.raises(AttributeError):
        r.by_source[0].admitted = 1
    with pytest.raises(TypeError):
        ngc.SourceFlow()
    assert "SourceFlow(source=0" in repr(r.by_source[0])


# -- evaluators ---------------------------------------------------------------------------------


def test_evaluators():
    g = graph(B1)
    # V1, V2.
    for config, placed in (
        (glob(ORIGIN.EQUAL), 4),
        (glob(), 12),
        (cost_only(P.MAX_FLOW), 12),
        (cost_only(P.MAX_FLOW, origin=ORIGIN.EQUAL), 4),
        (cost_only(), 4),
    ):
        assert ngc.capacity(g, [0, 1], 3, ngc.QueryOptions(config)).placed == placed
    # V3, V4.
    assert (
        ngc.capacity(graph(POOL), [0, 1], 3, ngc.QueryOptions(cost_only())).placed == 2
    )
    r = ngc.capacity(
        graph(PRIV), [0, 1], 2, ngc.QueryOptions(cost_only(admission=A.LOSSY))
    )
    assert (r.delivered, r.offered, r.requested, r.termination) == (
        10,
        10,
        INF,
        T.EXHAUSTED,
    )
    # V5.
    pairs = [([0, 1], 3), (0, 3), (1, 3)]
    got = ngc.batch_capacity(g, pairs, ngc.QueryOptions(glob(ORIGIN.EQUAL)))
    assert [x.placed for x in got] == [4, 2, 10]
    # V6 B3-delete-zero: every allowed edge is a candidate under the equal origin.
    for workers in (1, 2):
        s = ngc.sensitivity(
            graph(B3Z), [0, 1], 2, ngc.QueryOptions(glob(ORIGIN.EQUAL)), workers
        )
        assert (s.baseline.placed, [(e.edge, e.loss) for e in s.entries]) == (
            0,
            [(0, -10), (1, 0)],
        )
    # V7: pooled CostOnly: deleting e0 takes source 0 out of the roster.
    s = ngc.sensitivity(g, [0, 1], 3, ngc.QueryOptions(cost_only()))
    assert (s.baseline.placed, [(e.edge, e.loss) for e in s.entries]) == (
        4,
        [(0, -6), (1, 2), (2, 2)],
    )


# -- random controllers against the references ------------------------------------------------


def random_controller_case(rng):
    n = rng.randint(3, 6)
    nodes = list(range(n))
    rng.shuffle(nodes)
    k = rng.randint(2, min(3, n - 1)) if rng.random() < 0.8 else 1
    sources = nodes[:k]
    dests = nodes[k : k + rng.randint(1, min(2, n - k))]
    edges = []
    for _ in range(rng.randint(n, 2 * n + 2)):
        u = rng.randrange(n)
        v = rng.choice(dests) if rng.random() < 0.3 else rng.randrange(n)
        if u != v:
            edges.append((u, v, rng.randint(0, 4), rng.randint(1, 4)))
    if not edges:
        edges.append((sources[0], dests[0], 1, 1))
    return n, edges, sources, dests


@pytest.mark.parametrize("seed", range(300))
def test_controllers_random(seed):
    # G16, M8: Global and CostOnly MAX_FLOW under both origins on random graphs with masks and
    # nonuniform costs: the delivery, the minimum cost and (equal) the supplies RR(T) equal the
    # references; increments re-solve the admitted total plus the request.
    rng = random.Random(seed)
    n, edges, sources, dests = random_controller_case(rng)
    g = build(n, edges)
    m = g.num_edges()
    xn = [
        v
        for v in range(n)
        if v not in sources and v not in dests and rng.random() < 0.15
    ]
    xe = [e for e in range(m) if rng.random() < 0.1]
    elig = g.eligibility(excluded_nodes=xn, excluded_edges=xe)
    src = [int(x) for x in g.edge_src()]
    dst = [int(x) for x in g.edge_dst()]
    cost = [int(x) for x in g.cost()]
    cap = [int(x) for x in g.capacity_units()]
    node_ok = [v not in xn for v in range(n)]
    edge_ok = [e not in xe and node_ok[src[e]] and node_ok[dst[e]] for e in range(m)]
    dist = distances_to(n, src, dst, cost, set(dests), node_ok, edge_ok)
    roster = sorted(s for s in set(sources) if dist[s] is not None)
    ref = PooledReference(
        g, sources, dests, "EB", False, cap, excluded_nodes=xn, excluded_edges=xe
    )
    support = [False] * m
    if ref.roster:
        ref._prepare()
        for s in ref.roster:
            for e in ref.route[s]:
                support[e] = True
    for routing, allowed in (("global", edge_ok), ("costonly", support)):
        for origin in (ORIGIN.EQUAL, ORIGIN.OPTIMIZED):
            config = (
                glob(origin)
                if routing == "global"
                else cost_only(P.MAX_FLOW, origin=origin)
            )
            oracle = equal_origin if origin == ORIGIN.EQUAL else optimized_origin
            best, best_cost = oracle(g, roster, dests, cap, allowed)
            r = ngc.capacity(g, sources, dests, ngc.QueryOptions(config, elig=elig))
            where = (seed, routing, origin.name)
            assert (r.placed, r.total_cost) == (best, best_cost), where
            if origin == ORIGIN.EQUAL:
                admitted = {x.source: x.admitted for x in r.by_source}
                assert [admitted[s] for s in roster] == rr(best, len(roster)), where
            # Increments: a then b admit min(max, a) and then min(max, that + b).
            a, b = rng.randint(0, 6), rng.randint(0, 6)
            p = ngc.FlowPolicy(ngc.FlowLedger(g), sources, dests, config, elig)
            first = p.place(a)
            second = p.place(b)
            want_a = min(best, a)
            want_b = min(best, want_a + b)
            assert (first.placed, second.placed) == (want_a, want_b - want_a), where
            _, cost_b = oracle(g, roster, dests, cap, allowed, want_b)
            assert second.total_cost == cost_b, where
