"""T2 (CentralizedTE(tunnel_count=TunnelCount(N)), TunnelCount(L, up_to=True)) through the binding: the worked T2
examples and calls (docs/semantics.md, "T2: the controller's equal tunnels"), the configuration and call checks, sizing, objectives,
count selection, operations, programs and replays, records, flags and terminations, the
evaluators, the numeric ranges of T2's sums, zero offers and everything masked, a property test
against the independent tunnel_oracle.py, and the Python-allocation sweeps of T2's results.
Expected values are derived by hand from those rules, never taken from a build."""

from __future__ import annotations

import random
from fractions import Fraction

import pytest
import tunnel_oracle as orc
from helpers import (
    COSTONLY,
    GLOBAL,
    INF,
    L62,
    build,
    dag_of,
    exact,
    full_ledger_state,
    sweep_python_allocation,
)
from private_hooks import needs_hooks

import netgraph_core as ngc

T = ngc.Termination
S = ngc.TunnelStatus
UP, FAILED, BLACKHOLED = S.UP, S.FAILED, S.BLACKHOLED
LOSSY = ngc.Admission.LOSSY
LOSSLESS = ngc.Admission.LOSSLESS


def graph(rows, n=None):
    """A graph whose rows are in EdgeId order (asserted): (u, v, capacity, cost)."""
    n = n if n is not None else 1 + max(max(u, v) for u, v, _, _ in rows)
    g = build(n, rows)
    assert [int(i) for i in g.input_positions()] == list(range(len(rows)))
    return g


def t2_config(
    count, admission=None, headroom=False, up_to=False, source_allocation=None
):
    """A T2 FlowPolicyConfig: CentralizedTE(source_allocation=source_allocation (EQUAL by default),
    tunnel_count=TunnelCount(count, up_to=up_to), prefer_headroom=headroom), Lossless by
    default."""
    return ngc.FlowPolicyConfig(
        routing=ngc.CentralizedTE(
            source_allocation=source_allocation or ngc.SourceAllocation.EQUAL,
            tunnel_count=ngc.TunnelCount(count, up_to=up_to),
            prefer_headroom=headroom,
        ),
        admission=admission or LOSSLESS,
    )


def policy(
    g,
    s,
    d,
    n,
    admission=None,
    ledger=None,
    eligibility=None,
    headroom=False,
    up_to=False,
    source_allocation=None,
):
    ledger = ledger if ledger is not None else ngc.FlowLedger(g)
    return ngc.FlowPolicy(
        ledger,
        s,
        d,
        t2_config(n, admission, headroom, up_to, source_allocation),
        eligibility,
    )


def routes(r):
    return [list(t.edges) for t in r.program.tunnels]


def counts(r, m):
    x = [0] * m
    for t in r.program.tunnels:
        for e in t.edges:
            x[e] += 1
    return x


def rows(r):
    return [
        (t.entry_id, t.head, t.status, t.assigned, t.delivered, t.total_edge_drops)
        for t in r.by_tunnel
    ]


def source_rows(r):
    return [
        (
            s.source,
            s.offered,
            s.admitted,
            s.delivered,
            s.source_drops,
            s.total_edge_drops,
        )
        for s in r.by_source
    ]


def plan_tuple(p):
    return (
        p.tunnels_per_head,
        p.design_demand,
        p.design_by_head,
        p.total_reserved_bandwidth,
        p.heads_with_capacity,
        p.insufficient_resolution,
        p.max_fractional_bandwidth_per_tunnel,
        p.max_bandwidth_per_tunnel,
        p.planned_bandwidth_per_tunnel,
        p.total_planned_bandwidth,
        p.total_tunnel_cost,
        p.headroom_penalty,
        p.max_lossless_offer,
    )


G1 = [(0, 1, 8192, 1), (0, 2, 40960, 1), (1, 3, 8192, 1), (2, 3, 40960, 1)]
G4 = [(0, i, 4096, 1) for i in range(1, 16)] + [(i, 16, 4096, 1) for i in range(1, 16)]
G10 = [
    (0, 1, 40960, 1), (0, 5, 40960, 1), (0, 3, 24576, 2), (1, 0, 40960, 1), (1, 2, 40960, 1), (1, 4, 16384, 2),
    (2, 1, 40960, 1), (2, 3, 40960, 1), (2, 5, 20480, 2), (3, 2, 40960, 1), (3, 4, 40960, 1), (3, 0, 24576, 2),
    (4, 3, 40960, 1), (4, 5, 40960, 1), (4, 1, 16384, 2), (5, 0, 40960, 1), (5, 4, 40960, 1), (5, 2, 20480, 2),
]  # fmt: skip
G10_CAPS = {
    "G10 head 0": [40960, 40960, 24576, 40960, 40960, 16384, 40960, 40960, 20480,
                40960, 40960, 24576, 40960, 40960, 16384, 40960, 40960, 20480],
    "G10 head 1": [28672, 28672, 0, 40960, 28672, 16384, 40960, 28672, 20480,
                40960, 40960, 24576, 28672, 40960, 16384, 40960, 28672, 20480],
    "G10 head 2": [28672, 10240, 0, 22528, 13312, 1024, 40960, 13312, 20480,
                40960, 25600, 24576, 28672, 40960, 16384, 40960, 10240, 20480],
}  # fmt: skip
G1_PATHS = [[0, 2]] * 2 + [[1, 3]] * 14


def twoh(a, b):
    return graph([(0, 2, a, 1), (1, 2, b, 1)])


def one10():
    return graph([(0, 1, 10, 1)])


# ---------------------------------------------------------------- the worked examples (X)

# id: (rows, heads, dst, N, D, free capacity or None, prefer_headroom, expected plan fields, routes, lossy offer, delivery,
#      source loss). Plan fields: (U, u, fluid, F, planned total, route cost, Phi or None, live heads, indices or None).
EXAMPLES = {
    "G1 design 49152": (G1, [0], 3, 16, 49152, None, True,
             (2925, 2925, Fraction(20480, 7), 46812, 46800, 32, 73610035200, (0,), None), G1_PATHS, 49152, 47104, 0),
    "G1 design 8192": (G1, [0], 3, 16, 8192, None, True,
                    (2925, 512, Fraction(20480, 7), 46812, 8192, 32, 12884901876, (0,), None), G1_PATHS, None, None, None),
    "G4 design 61440": (G4, [0], 16, 16, 61440, None, True,
             (2048, 2048, Fraction(2048), 32768, 32768, 32, 73014444032, (0,), None),
             [[0, 15], [0, 15]] + [[i, 15 + i] for i in range(1, 15)], 61440, 57856, 0),
    "heads behind 3 and 8": ([(0, 2, 3, 1), (1, 2, 8, 1)], [0, 1], 2, 2, 6, None, True,
                (1, 1, Fraction(3, 2), 6, 4, 4, 5905580031, (0, 1), [0, 1, 2, 3]), [[0], [1], [0], [1]], None, None, None),
    "a cheap link of 2": ([(0, 1, 2, 1), (0, 1, 10, 2)], [0], 1, 2, 2, None, True,
                      (5, 1, Fraction(5), 2, 2, 2, 6442450944, (0,), None), [[0], [0]], None, None, None),
    "G10 head 0": (G10, [0], 3, 16, 49152, G10_CAPS["G10 head 0"], True,
                (6144, 3072, Fraction(6144), 49160, 49152, 40, 38654705652, (0,), None),
                [[0, 4, 7]] * 4 + [[1, 16, 12]] * 4 + [[2]] * 8, None, None, None),
    "G10 head 1": (G10, [1], 4, 16, 49152, G10_CAPS["G10 head 1"], True,
                (4096, 3072, Fraction(4096), 52431, 49152, 43, None, (1,), None),
                [[3, 1, 16]] * 6 + [[4, 7, 10]] * 5 + [[5]] * 5, None, None, None),
    "G10 head 2": (G10, [2], 5, 16, 49152, G10_CAPS["G10 head 2"], True,
                (2560, 2560, Fraction(2560), 40960, 40960, 40, None, (2,), None),
                [[6, 3, 1]] * 4 + [[7, 10, 13]] * 4 + [[8]] * 8, None, None, None),
    "a dead head": ([(0, 2, 0, 1), (1, 2, 10, 1)], [0, 1], 2, 2, 10, None, True,
                     (5, 3, Fraction(5), 0, 6, 2, 3865470565, (1,), [1, 3]), [[1], [1]], 10, 5, 5),
    "heads behind 2 and 10": ([(0, 2, 2, 1), (1, 2, 10, 1)], [0, 1], 2, 2, None, None, True,
                           (1, 1, Fraction(1), 4, 4, 4, 7730941132, (0, 1), [0, 1, 2, 3]), [[0], [1], [0], [1]],
                           None, None, None),
    "an idle tunnel": ([(0, 1, 10, 1)], [0], 1, 4, 3, None, True,
                (2, 1, Fraction(5, 2), 10, 4, 4, 4294967294, (0,), None), [[0]] * 4, 3, 3, 0),
    "a zero design": ([(0, 1, 10, 1)], [0], 1, 2, 0, None, True,
                       (5, 1, Fraction(5), 10, 2, 2, None, (0,), None), [[0]] * 2, 4, 4, 0),
    "the remainder on the first link": ([(0, 1, 3, 1), (0, 1, 2, 1)], [0], 1, 2, 5, None, True,
                           (2, 2, Fraction(2), 5, 4, 2, None, (0,), None), [[0], [1]], None, None, None),
    "the remainder on the second link": ([(0, 1, 2, 1), (0, 1, 3, 1)], [0], 1, 2, 5, None, True,
                            (2, 2, Fraction(2), 4, 4, 2, None, (0,), None), [[0], [1]], None, None, None),
    "a tiny head": ([(0, 2, 1, 1), (1, 2, 100, 1)], [0, 1], 2, 16, None, None, True,
                     (0, 0, Fraction(1, 16), 0, 0, 0, 0, (0, 1), None), [], None, None, None),
}  # fmt: skip


@pytest.mark.parametrize("eid", list(EXAMPLES))
def test_contract_examples(eid):
    rows_, heads, dst, N, D, free, hr, want, rts, offer, delivery, loss = EXAMPLES[eid]
    g = graph(rows_)
    fg = ngc.FlowLedger(g, free) if free is not None else ngc.FlowLedger(g)
    r = policy(g, heads, dst, N, ledger=fg, headroom=hr).replace(
        INF, design_demand=INF if D is None else D
    )
    p = r.plan
    U, u, fluid, F, total, cost, phi, live, indices = want
    assert (
        p.max_bandwidth_per_tunnel,
        p.planned_bandwidth_per_tunnel,
        p.max_fractional_bandwidth_per_tunnel,
        p.max_lossless_offer,
    ) == (U, u, fluid, F)
    assert (p.total_planned_bandwidth, p.total_tunnel_cost, p.heads_with_capacity) == (
        total,
        cost,
        live,
    )
    if phi is not None:
        assert p.headroom_penalty == phi
    assert routes(r) == rts
    if indices is not None:
        assert [t.entry_id for t in r.program.tunnels] == indices
    assert all(
        (t.planned_bandwidth, t.weight, t.status) == (u, 1, UP)
        for t in r.program.tunnels
    )
    assert (r.call_delivered, r.termination) == (
        F,
        T.EXHAUSTED,
    )  # the unlimited query admits F
    assert p.insufficient_resolution == (U == 0)
    if offer is None:
        return
    ly = policy(
        g,
        heads,
        dst,
        N,
        LOSSY,
        ledger=ngc.FlowLedger(g, free) if free else None,
        headroom=hr,
    )
    s = ly.replace(offer, design_demand=D)
    assert (s.call_delivered, s.call_source_drops, s.termination) == (
        delivery,
        loss,
        T.COMPLETE,
    )
    assert routes(s) == rts  # the offer never changes the layout


UP_TO = {
    # id: (rows, heads, dst, L, D, chosen count, routes, lossless delivery of D (None: unlimited, the bound),
    #      the EXACT(k) bounds F'(k) of the all-live examples)
    "G1 up to 6": (G1, [0], 3, 6, 49152, 6, [[0, 2]] + [[1, 3]] * 5, 49152, [40960] * 5 + [49152]),
    "up to, a common size": ([(0, 3, 6, 1), (1, 3, 3, 1), (1, 3, 3, 1)], [0, 1], 3, 2, 12, 2, [[0], [1], [0], [2]], 12,
                        [7, 12]),
    "up to, a dead head": ([(0, 2, 0, 1), (1, 2, 10, 1)], [0, 1], 2, 2, 10, 1, [[1]], 0, None),
    "up to, a tiny head": ([(0, 2, 1, 1), (1, 2, 100, 1)], [0, 1], 2, 3, None, 1, [[0], [1]], 2, [2, 0, 0]),
}  # fmt: skip


@pytest.mark.parametrize("eid", list(UP_TO))
def test_up_to_examples(eid):
    rows_, heads, dst, L, D, k, rts, delivered, bounds = UP_TO[eid]
    g = graph(rows_)
    r = policy(g, heads, dst, L, up_to=True, headroom=True).replace(
        INF if D is None else D
    )
    assert (r.plan.tunnels_per_head, routes(r), r.call_delivered) == (k, rts, delivered)
    if bounds is None:
        return
    # With every head live, F'(k) is EXACT(k)'s bound at U(k) (an unlimited design).
    for c, b in enumerate(bounds, 1):
        assert (
            policy(g, heads, dst, c, headroom=True).replace(INF).call_delivered == b
        ), c


def test_g1_call():
    # The G1 call, then its program's fields and the ledger.
    q = 2.0**-12
    g = ngc.Graph(
        4,
        edge_sources=[0, 0, 1, 2],
        edge_destinations=[1, 2, 3, 3],
        edge_capacity_units=ngc.to_units([2, 10, 2, 10], q),
        edge_costs=[1, 1, 1, 1],
        resolution=q,
    )
    te16 = t2_config(16, ngc.Admission.LOSSY, headroom=True)
    fg = ngc.FlowLedger(g)
    p = ngc.FlowPolicy(fg, 0, 3, te16)
    r = p.replace(49152)
    assert (
        r.plan.max_bandwidth_per_tunnel,
        r.plan.planned_bandwidth_per_tunnel,
        r.plan.max_lossless_offer,
    ) == (2925, 2925, 46812)
    assert r.call_delivered == 47104
    assert plan_tuple(r.plan) == (
        16,
        49152,
        ((0, 49152),),
        0,
        (0,),
        False,
        Fraction(20480, 7),
        2925,
        2925,
        46800,
        32,
        73610035200,
        46812,
    )
    assert [list(t.edges) for t in r.program.tunnels] == G1_PATHS
    assert [
        (
            t.requested_bandwidth,
            t.reserved_bandwidth,
            t.planned_bandwidth,
            t.status,
            t.weight,
        )
        for t in r.program.tunnels
    ] == [(0, 0, 2925, UP, 1)] * 16
    assert [t.delivered for t in r.by_tunnel] == [3072] * 2 + [2926] * 10 + [2925] * 4
    assert r.call_edge_drops == [(1, 2048)]
    assert r.reservations_changed is False and (
        r.routes_changed,
        r.weights_changed,
    ) == (True, True)
    # A plan books nothing.
    assert fg.reserved().tolist() == [0, 0, 0, 0]
    assert all(rec.edge_reservations == [] for rec in p.records())
    assert ngc.CentralizedTE(
        source_allocation=ngc.SourceAllocation.EQUAL,
        tunnel_count=ngc.TunnelCount(256, up_to=True),
    ).tunnel_count == ngc.TunnelCount(256, up_to=True)


def test_masked_replay_call():
    # The replay under the mask, then the same program stale (black holes).
    g = graph(G1)
    p = policy(g, 0, 3, 16, LOSSY, headroom=True)
    p.replace(49152)
    base = p.program
    fail = g.eligibility(excluded_edges=[2])
    replay = ngc.FlowPolicy(
        ngc.FlowLedger(g),
        0,
        3,
        ngc.FlowPolicyConfig(routing=ngc.Pinned(base), admission=LOSSY),
        eligibility=fail,
    )
    r = replay.replace(49152)
    assert [t.status for t in r.by_tunnel][:3] == [FAILED, FAILED, UP]
    assert r.call_delivered == 40960 and r.plan is None
    stale = ngc.FlowPolicy(
        ngc.FlowLedger(g),
        0,
        3,
        ngc.FlowPolicyConfig(
            routing=ngc.Pinned(base, black_holes=fail), admission=LOSSY
        ),
    )
    s = stale.replace(49152)
    assert [t.status for t in s.by_tunnel] == [BLACKHOLED] * 2 + [UP] * 14
    assert (s.call_edge_drops, s.call_delivered) == ([(1, 2048), (2, 6144)], 40960)
    assert [t.status for t in base.tunnels] == [
        UP
    ] * 16  # the program's own statuses stay


# ---------------------------------------------------------------- validation (V)


SEVERAL = "CentralizedTE with tunnel_count and several sources requires source_allocation=EQUAL"


def cte(source_allocation=ngc.SourceAllocation.EQUAL, **fields):
    return ngc.CentralizedTE(source_allocation=source_allocation, **fields)


def test_matrix():
    # T2's row of the matrix, in order (the evaluators run the same checks). CentralizedTE has
    # no route choice and no weights: those keywords are no keywords of it.
    g = twoh(10, 10)
    fg = ngc.FlowLedger(g)
    t2 = cte(tunnel_count=ngc.TunnelCount(2))
    opt = cte(ngc.SourceAllocation.OPTIMIZED, tunnel_count=ngc.TunnelCount(2))
    dag = dag_of(g, [0, 1])
    both = cte(support=dag, tunnel_count=ngc.TunnelCount(2), prefer_headroom=True)
    for kw in (
        {"allocation": ngc.RouteAllocation.EQUAL},
        {"weight_basis": ngc.WeightBasis.CAPACITY},
        {"edge_weights": [1, 1]},
        {"tie_break": ngc.RouteTieBreak.FEWER_HOPS},
        {"multipath": False},
        {"parallel_edges": False},
    ):
        with pytest.raises(TypeError):
            cte(tunnel_count=ngc.TunnelCount(2), **kw)
    for kw in (
        {"placement": ngc.RouteAllocation.EQUAL},
        {"origin": ngc.SourceAllocation.EQUAL},
    ):
        with pytest.raises(TypeError):
            ngc.FlowPolicyConfig(routing=t2, **kw)
    cases = [
        # The several-sources rule first, then prefer_headroom without tunnel_count.
        (ngc.FlowPolicyConfig(routing=opt), SEVERAL),
        (ngc.FlowPolicyConfig(routing=cte(ngc.SourceAllocation.OPTIMIZED, support=dag,
                                          tunnel_count=ngc.TunnelCount(2), prefer_headroom=True),
                              admission=LOSSY), SEVERAL),
        (ngc.FlowPolicyConfig(routing=cte(prefer_headroom=True)), "prefer_headroom requires tunnel_count"),
        (ngc.FlowPolicyConfig(routing=cte(support=dag, prefer_headroom=True)),
         "prefer_headroom requires tunnel_count"),
    ]  # fmt: skip
    for cfg, text in cases:
        with pytest.raises(ValueError, match=exact(text)):
            ngc.FlowPolicy(fg, [0, 1], 2, cfg)
        with pytest.raises(ValueError, match=exact(text)):
            ngc.evaluate_demand(g, [0, 1], 2, ngc.EvaluationOptions(cfg))
    with pytest.raises(
        ValueError, match=exact("CentralizedTE does not use ShortestPathCache")
    ):
        ngc.FlowPolicy(
            fg,
            [0, 1],
            2,
            ngc.FlowPolicyConfig(routing=t2),
            cache=ngc.ShortestPathCache(g),
        )
    # The several-sources rule counts the normalized input sources, before the mask; one source
    # with OPTIMIZED is accepted (the same calculation).
    with pytest.raises(ValueError, match=exact(SEVERAL)):
        ngc.FlowPolicy(
            fg,
            [0, 1],
            2,
            ngc.FlowPolicyConfig(routing=opt),
            g.eligibility(excluded_nodes=[1]),
        )
    one = ngc.FlowPolicy(fg, [0, 0], 2, ngc.FlowPolicyConfig(routing=opt))
    assert one.config.routing.source_allocation == ngc.SourceAllocation.OPTIMIZED
    # Accepted: both admissions, with and without prefer_headroom, with and without a support.
    for adm in (LOSSLESS, LOSSY):
        for hr in (False, True):
            ngc.FlowPolicy(fg, [0, 1], 2, t2_config(2, adm, hr))
        ngc.FlowPolicy(fg, [0, 1], 2, ngc.FlowPolicyConfig(routing=both, admission=adm))
    got = ngc.FlowPolicy(fg, [0, 1], 2, ngc.FlowPolicyConfig(routing=t2)).config
    assert (got.routing.source_allocation, got.routing.tunnel_count) == (
        ngc.SourceAllocation.EQUAL,
        ngc.TunnelCount(2),
    )
    # The arguments of CentralizedTE itself: source_allocation required; prefer_headroom False
    # by default.
    for bad in ({}, {"source_allocation": None}):
        with pytest.raises(TypeError):
            ngc.CentralizedTE(**bad)
    with pytest.raises(TypeError):
        cte(tunnel_count=2)
    with pytest.raises(TypeError):
        cte(tunnel_count=ngc.TunnelCount(2), prefer_headroom="yes")
    plain = cte()
    assert (plain.support, plain.tunnel_count, plain.prefer_headroom) == (
        None,
        None,
        False,
    )
    assert cte(prefer_headroom=True).prefer_headroom is True
    eq = "source_allocation=SourceAllocation.EQUAL"
    assert (
        repr(cte(tunnel_count=ngc.TunnelCount(256, up_to=True)))
        == f"CentralizedTE({eq}, tunnel_count=TunnelCount(256, up_to=True))"
    )
    assert (
        repr(cte(tunnel_count=ngc.TunnelCount(16), prefer_headroom=True))
        == f"CentralizedTE({eq}, tunnel_count=TunnelCount(16), prefer_headroom=True)"
    )
    assert (
        repr(cte(prefer_headroom=True)) == f"CentralizedTE({eq}, prefer_headroom=True)"
    )
    assert repr(cte()) == f"CentralizedTE({eq})"
    assert (
        repr(cte(ngc.SourceAllocation.OPTIMIZED))
        == "CentralizedTE(source_allocation=SourceAllocation.OPTIMIZED)"
    )
    # Every set field, the invalid support-and-tunnels combination included.
    dag = dag_of(twoh(10, 10), [0, 1])
    assert repr(cte(support=dag)) == f"CentralizedTE({eq}, support=...)"
    assert (
        repr(cte(support=dag, tunnel_count=ngc.TunnelCount(2)))
        == f"CentralizedTE({eq}, support=..., tunnel_count=TunnelCount(2))"
    )
    assert (
        repr(cte(support=dag, prefer_headroom=True))
        == f"CentralizedTE({eq}, support=..., prefer_headroom=True)"
    )


def test_global_prefer_headroom_texts():
    # CentralizedTE's own cases in the model matrix's order, at a policy and at an evaluator:
    # prefer_headroom without tunnel_count is refused, with or without a support; support and
    # tunnels together are T2 on the support.
    g = twoh(10, 10)
    fg = ngc.FlowLedger(g)
    dag = dag_of(g, [0, 1])
    t2 = ngc.TunnelCount(2)
    cases = [
        (dict(prefer_headroom=True), "prefer_headroom requires tunnel_count"),
        (dict(prefer_headroom=True, support=dag), "prefer_headroom requires tunnel_count"),
        (dict(prefer_headroom=True, support=dag, tunnel_count=t2), None),
        (dict(prefer_headroom=True, tunnel_count=t2), None),
        ({}, None),
    ]  # fmt: skip
    for kwargs, text in cases:
        cfg = ngc.FlowPolicyConfig(routing=cte(**kwargs))
        if text is None:
            ngc.FlowPolicy(fg, 0, 2, cfg)
            ngc.evaluate_demand(g, 0, 2, ngc.EvaluationOptions(cfg))
            continue
        with pytest.raises(ValueError, match=exact(text)):
            ngc.FlowPolicy(fg, 0, 2, cfg)
        with pytest.raises(ValueError, match=exact(text)):
            ngc.evaluate_demand(g, 0, 2, ngc.EvaluationOptions(cfg))
    # No old spelling of the objective is accepted.
    with pytest.raises(TypeError):
        cte(tunnel_count=t2, priority=ngc.RouteTieBreak.HIGHER_HEADROOM)


# ---------------------------------------------------------------- a supplied support


SUPPORTED = [(0, 1, 6, 1), (0, 1, 6, 1), (0, 1, 100, 2)]


@pytest.mark.parametrize(
    ("excluded", "want"),
    [
        # T2-support: the support holds links 0 and 1; U = 6, one tunnel on each, the third link
        # never used, the lossless maximum 12.
        (None, (6, [1, 1, 0], [[0], [1]], [0, 1], 12)),
        # T2-support-masked: the mask also excludes link 1: U = 3, both tunnels on link 0.
        ([1], (3, [2, 0, 0], [[0], [0]], [0, 1], 6)),
        # T2-support-mask-outside: a mask outside the support changes nothing.
        ([2], (6, [1, 1, 0], [[0], [1]], [0, 1], 12)),
    ],
)
def test_support_examples(excluded, want):
    g = graph(SUPPORTED, 2)
    cfg = ngc.FlowPolicyConfig(
        routing=cte(support=dag_of(g, [0, 1]), tunnel_count=ngc.TunnelCount(2))
    )
    mask = None if excluded is None else g.eligibility(excluded_edges=excluded)
    r = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 1, cfg, mask).replace(INF)
    got = (
        r.plan.max_bandwidth_per_tunnel,
        counts(r, 3),
        routes(r),
        [t.entry_id for t in r.program.tunnels],
        r.plan.max_lossless_offer,
    )
    assert got == want
    assert r.call_delivered == want[4]


def on_sup(g, support=(0, 1), admission=LOSSLESS, headroom=False):
    return ngc.FlowPolicyConfig(
        routing=cte(
            support=dag_of(g, list(support)),
            tunnel_count=ngc.TunnelCount(2),
            prefer_headroom=headroom,
        ),
        admission=admission,
    )


def test_support_equals_mask():
    # T2 on {e0, e1} against T2 with a mask excluding e2: the same plan, tunnels, rows, deliveries,
    # drops and flags, Lossless and Lossy at offers 6, 12 and 13.
    g = graph(SUPPORTED, 2)
    outside = g.eligibility(excluded_edges=[2])
    for admission in (LOSSLESS, LOSSY):
        for offer in (6, 12, 13):
            got = []
            for cfg, mask in (
                (on_sup(g, admission=admission), None),
                (t2_config(2, admission), outside),
            ):
                p = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 1, cfg, mask)
                r = p.replace(offer)
                got.append(
                    (
                        plan_tuple(r.plan),
                        routes(r),
                        [t.planned_bandwidth for t in r.program.tunnels],
                        rows(r),
                        source_rows(r),
                        r.call_delivered,
                        r.call_edge_drops,
                        r.termination,
                        (r.routes_changed, r.weights_changed, r.reservations_changed),
                    )
                )
            assert got[0] == got[1], (admission, offer)


def test_support_edge_failures():
    # Lossless, no offer: the baseline 12; deleting e0 or e1 halves it, deleting e2 (outside the
    # support) changes nothing.
    g = graph(SUPPORTED, 2)
    for workers in (1, 2):
        s = ngc.evaluate_edge_failures(
            g, 0, 1, ngc.EvaluationOptions(on_sup(g)), workers
        )
        assert s.baseline.call_delivered == 12
        assert [(x.edge, x.delivery_reduction) for x in s.entries] == [
            (0, 6),
            (1, 6),
            (2, 0),
        ]


def test_support_empty_and_masked():
    # An empty support: no head reaches a destination (NO_ROUTE, the Lossy offer lost at the
    # source; a zero call COMPLETE); a mask excluding every edge does the same.
    g = graph(SUPPORTED, 2)
    empty = on_sup(g, support=(), admission=LOSSY)
    r = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 1, empty).replace(2)
    assert (r.termination, r.call_source_drops, r.call_delivered) == (T.NO_ROUTE, 2, 0)
    r = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 1, empty).replace(0)
    assert (r.termination, r.call_delivered) == (T.COMPLETE, 0)
    every = g.eligibility(excluded_edges=[0, 1, 2])
    r = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 1, on_sup(g), every).replace(2)
    assert (r.termination, r.call_delivered) == (T.NO_ROUTE, 0)


def test_support_headroom():
    # prefer_headroom on the support, N = 2, design 2: U = 6, u = 1; the second tunnel takes e1 (a
    # second slot on e0 costs more); the penalty is 2 floor(2^32 1 / 6).
    g = graph(SUPPORTED, 2)
    r = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 1, on_sup(g, headroom=True)).replace(
        INF, design_demand=2
    )
    assert (r.plan.max_bandwidth_per_tunnel, r.plan.planned_bandwidth_per_tunnel) == (
        6,
        1,
    )
    assert counts(r, 3) == [1, 1, 0]
    assert r.plan.headroom_penalty == 2 * (2**32 // 6) == 1431655764


def _composed(g, support_edges, mask_nodes, mask_edges):
    """The policy's mask with every edge outside the support excluded too."""
    outside = [e for e in range(g.num_edges()) if e not in set(support_edges)]
    return g.eligibility(
        excluded_nodes=mask_nodes, excluded_edges=sorted(set(mask_edges) | set(outside))
    )


@pytest.mark.parametrize("seed", range(80))
def test_support_is_the_composed_mask(seed):
    # T2 on a support gives what the same T2 configuration gives with the support's edges composed
    # into the mask: the plan, the routes, the rows, the deliveries, at a policy (replace, then a
    # forward) and at the three evaluators.
    rng = random.Random(seed)
    n = rng.randint(3, 6)
    # Forward links only: every support is a DAG.
    rows_ = []
    for _ in range(rng.randint(n, 3 * n)):
        u, v = sorted(rng.sample(range(n), 2))
        rows_.append((u, v, rng.choice([0, 1, 3, 6, 10]), rng.randint(1, 4)))
    g = build(n, rows_)
    m = g.num_edges()
    support_edges = sorted(rng.sample(range(m), rng.randint(1, m)))
    mask_nodes = [v for v in range(1, n - 1) if rng.random() < 0.15]
    mask_edges = [e for e in range(m) if rng.random() < 0.15]
    mask = g.eligibility(excluded_nodes=mask_nodes, excluded_edges=mask_edges)
    composed = _composed(g, support_edges, mask_nodes, mask_edges)
    sources = sorted(rng.sample(range(n - 1), rng.randint(1, min(2, n - 1))))
    count = ngc.TunnelCount(rng.randint(1, 3), up_to=rng.random() < 0.3)
    headroom = rng.random() < 0.5
    admission = rng.choice([LOSSLESS, LOSSY])
    on_support = ngc.FlowPolicyConfig(
        routing=cte(
            support=dag_of(g, support_edges),
            tunnel_count=count,
            prefer_headroom=headroom,
        ),
        admission=admission,
    )
    plain_t2 = ngc.FlowPolicyConfig(
        routing=cte(tunnel_count=count, prefer_headroom=headroom), admission=admission
    )
    offer = rng.choice([0, 3, 11, 40])
    design = rng.choice([None, 0, 5, 30])

    def call(cfg, eligibility):
        p = ngc.FlowPolicy(ngc.FlowLedger(g), sources, n - 1, cfg, eligibility)
        r = p.replace(offer, design_demand=design)
        f = p.forward(offer + 1)
        return (
            plan_tuple(r.plan),
            routes(r),
            rows(r),
            source_rows(r),
            r.call_delivered,
            r.call_edge_drops,
            r.call_source_drops,
            rows(f),
            f.call_delivered,
        )

    assert call(on_support, mask) == call(plain_t2, composed), seed
    options = dict(offer=offer if admission == LOSSY else None, design_demand=design)
    a = ngc.EvaluationOptions(on_support, eligibility=mask, **options)
    b = ngc.EvaluationOptions(plain_t2, eligibility=composed, **options)
    ra, rb = (
        ngc.evaluate_demand(g, sources, n - 1, a),
        ngc.evaluate_demand(g, sources, n - 1, b),
    )
    assert (plan_tuple(ra.plan), routes(ra), ra.call_delivered) == (
        plan_tuple(rb.plan),
        routes(rb),
        rb.call_delivered,
    ), seed
    (pa,) = ngc.evaluate_pairs(g, [(sources, n - 1)], a)
    assert pa.call_delivered == ra.call_delivered
    # The deletions of the support's allowed edges match; a deletion outside the support changes
    # nothing.
    sa = ngc.evaluate_edge_failures(g, sources, n - 1, a, 1)
    sb = ngc.evaluate_edge_failures(g, sources, n - 1, b, 1)
    reduction_b = {x.edge: x.delivery_reduction for x in sb.entries}
    for x in sa.entries:
        assert x.delivery_reduction == reduction_b.get(x.edge, 0), (seed, x.edge)


@pytest.mark.parametrize("on", [False, True])
def test_prefer_headroom_example(on):
    # prefer_headroom off and on: two equal-cost parallel links of 10, one head, N = 2, D = 2:
    # U = 10, u = 1. Off, the canonical tie puts both tunnels on link 0 (bound 10; Lossy 12
    # delivers 10). On, the second tunnel pays floor(2^32 / 10) on link 1 against
    # floor(2^32 2 / 10) for a second slot on link 0: Phi = 2 floor(2^32 / 10) = 858993458, bound
    # 20, Lossy 12 delivers 12. Both plan 2.
    g = graph([(0, 1, 10, 1), (0, 1, 10, 1)])
    r = policy(g, 0, 1, 2, LOSSY, headroom=on).replace(12, design_demand=2)
    p = r.plan
    assert (
        p.max_bandwidth_per_tunnel,
        p.planned_bandwidth_per_tunnel,
        p.total_planned_bandwidth,
        p.total_tunnel_cost,
    ) == (
        10,
        1,
        2,
        2,
    )
    assert routes(r) == ([[0], [1]] if on else [[0], [0]])
    assert p.headroom_penalty == (858993458 if on else None)
    assert (p.max_lossless_offer, r.call_delivered) == ((20, 12) if on else (10, 10))


def test_count_range_and_tunnel_limit():
    # The count range (the count and the largest count of "up to"); the tunnel limit uses L for
    # "up to".
    fg = ngc.FlowLedger(one10())
    for n in (0, 10001, -1):
        for up in (False, True):
            with pytest.raises(
                ValueError, match=exact("tunnel count must be in [1, 10000]")
            ):
                policy(None, 0, 1, n, ledger=fg, up_to=up)
    policy(None, 0, 1, 10000, ledger=fg)
    policy(None, 0, 1, 10000, ledger=fg, up_to=True)
    # 1000 heads x 10000 = 10^7 tunnels: within the limit (2^30 - 1).
    rows_ = [(h, 1000, 1, 1) for h in range(1000)]
    policy(graph(rows_), list(range(1000)), 1000, 10000, up_to=True)


@needs_hooks
def test_traced_place_checks_the_operation_first():
    p = policy(one10(), 0, 1, 2)
    for v in (5, 0, -1, L62, 2**100):
        with pytest.raises(
            ValueError,
            match=exact(
                "CentralizedTE with tunnel_count does not support incremental placement"
            ),
        ):
            p._place_traced(v)


def test_operations_and_design_validation():
    # The operations of T2 and the design's domain, with their texts.
    g = one10()
    p = policy(g, 0, 1, 2)
    # The operation comes before the offer: also for an offer outside its domain.
    for v in (5, 0, -1, L62, 2**100):
        with pytest.raises(
            ValueError,
            match=exact(
                "CentralizedTE with tunnel_count does not support incremental placement"
            ),
        ):
            p.place(v)
    with pytest.raises(
        ValueError, match=exact("forward requires an installed program")
    ):
        p.forward(5)
    for d in (-1, L62):
        with pytest.raises(
            ValueError,
            match=exact(
                "design_demand must be an integer in [0, 2^62) or permitted +inf"
            ),
        ):
            p.replace(5, design_demand=d)
    with pytest.raises(
        ValueError,
        match=exact("offer must be an integer in [0, 2^62) or permitted +inf"),
    ):
        p.replace(-1, design_demand=-1)
    ly = policy(g, 0, 1, 2, LOSSY)
    with pytest.raises(
        ValueError, match=exact("LOSSY forwarding requires a finite offer")
    ):
        ly.replace(INF)
    assert p.replace(INF).plan.design_demand == INF  # an unlimited design
    for cfg in (COSTONLY, GLOBAL):
        with pytest.raises(
            ValueError,
            match=exact(
                f"{'ShortestPaths' if cfg is COSTONLY else 'CentralizedTE'} does not accept design_demand"
            ),
        ):
            ngc.FlowPolicy(ngc.FlowLedger(g), 0, 1, cfg).replace(5, design_demand=5)


# ---------------------------------------------------------------- sizing, objectives and count (U)


def test_headroom_off_and_finite_cost_tie():
    # G1 without prefer_headroom (the cost tie goes to the smallest keys); PAR412 with and without.
    g = graph(G1)
    r = policy(g, 0, 3, 16, LOSSY).replace(49152)
    assert (
        counts(r, 4),
        r.plan.max_lossless_offer,
        r.call_delivered,
        r.plan.headroom_penalty,
    ) == ([2, 14, 2, 14], 46812, 47104, None)
    r = policy(g, 0, 3, 16).replace(INF, design_demand=8192)
    assert (
        r.plan.planned_bandwidth_per_tunnel,
        counts(r, 4),
        r.plan.max_lossless_offer,
    ) == (
        512,
        [16, 0, 16, 0],
        8192,
    )
    r = policy(g, 0, 3, 16, headroom=True).replace(INF, design_demand=8192)
    assert (counts(r, 4), r.plan.max_lossless_offer) == ([2, 14, 2, 14], 46812)
    p412 = graph([(0, 1, 4, 1), (0, 1, 12, 1)])
    r = policy(p412, 0, 1, 2, headroom=True).replace(INF, design_demand=2)
    assert (
        r.plan.max_bandwidth_per_tunnel,
        r.plan.planned_bandwidth_per_tunnel,
        routes(r),
        r.plan.max_lossless_offer,
    ) == (6, 1, [[1], [1]], 12)
    assert r.plan.headroom_penalty == 1073741823
    r = policy(p412, 0, 1, 2).replace(INF, design_demand=2)
    assert (routes(r), r.plan.max_lossless_offer) == ([[0], [0]], 4)


# The graphs of T2's three canonical tie rules (u, v, capacity, cost), each derived in tests/cpp's twin.
ROOT_TIE_A = [
    (0, 2, 1, 1),
    (0, 4, 1, 1),
    (1, 2, 1, 1),
    (1, 5, 1, 1),
    (2, 3, 1, 1),
    (4, 3, 1, 1),
    (5, 3, 1, 1),
]
ROOT_TIE_B = [
    (0, 2, 1, 1),
    (0, 5, 1, 1),
    (1, 2, 1, 1),
    (1, 4, 1, 1),
    (2, 3, 1, 1),
    (4, 3, 1, 1),
    (5, 3, 1, 1),
]
FEWEST_A = [(0, 1, 1, 1), (0, 2, 1, 2), (1, 2, 1, 1)]
FEWEST_B = [(0, 1, 1, 1), (0, 3, 1, 3), (1, 2, 1, 1), (2, 3, 1, 1)]
CUT_A = [(0, 1, 2, 2), (1, 2, 2, 2), (1, 2, 3, 2)]
CUT_B = [
    (0, 2, 1, 1),
    (1, 2, 1, 1),
    (2, 3, 1, 1),
    (2, 4, 1, 1),
    (3, 5, 1, 1),
    (4, 5, 1, 1),
]


def head_routes(r):
    return [(t.head, list(t.edges)) for t in r.program.tunnels]


@pytest.mark.parametrize("headroom", [False, True])
def test_root_key_tie(headroom):
    # The heads' first walks tie (price 2, four arcs): the smaller head NodeId goes first; the
    # other head then takes its four-arc walk, not the six-arc one through a reverse arc.
    r = policy(graph(ROOT_TIE_A), [0, 1], 3, 1, headroom=headroom).replace(INF)
    assert head_routes(r) == [(0, [0, 4]), (1, [3, 6])]
    r = policy(graph(ROOT_TIE_B), [0, 1], 3, 1, headroom=headroom).replace(INF)
    assert head_routes(r) == [(0, [0, 4]), (1, [3, 5])]


def test_fewest_arcs_before_keys():
    # The direct link ties the longer walk on cost; fewer arcs wins before the smallest keys.
    assert routes(policy(graph(FEWEST_A), 0, 2, 1).replace(INF)) == [[1]]
    assert routes(policy(graph(FEWEST_B), 0, 3, 1).replace(INF)) == [[1]]


def test_cut_into_unit_walks():
    # Live heads ascending, each walk along the smallest positive EdgeId: head 0 cuts first.
    r = policy(graph(CUT_A), [0, 1], 2, 2).replace(INF)
    assert head_routes(r) == [(0, [0, 1]), (1, [2]), (0, [0, 1]), (1, [2])]
    r = policy(graph(CUT_B), [0, 1], 5, 1).replace(INF)
    assert head_routes(r) == [(0, [0, 2, 4]), (1, [1, 3, 5])]


def test_reverse_arcs_and_preloaded_prices():
    # A reverse arc: the second unit cancels the first unit's middle link (reverse arcs undo an
    # earlier choice); a preload: the slot prices read the static capacity, not the free one.
    g = graph([(0, 1, 1, 1), (0, 2, 1, 5), (1, 2, 1, 1), (1, 3, 1, 5), (2, 3, 1, 1)])
    r = policy(g, 0, 3, 2).replace(INF)
    assert (r.plan.max_bandwidth_per_tunnel, routes(r), r.plan.total_tunnel_cost) == (
        1,
        [[0, 3], [1, 4]],
        12,
    )
    r = policy(
        one10(), 0, 1, 2, ledger=ngc.FlowLedger(one10(), [6]), headroom=True
    ).replace(INF, design_demand=4)
    assert (
        r.plan.max_bandwidth_per_tunnel,
        r.plan.planned_bandwidth_per_tunnel,
        r.plan.headroom_penalty,
    ) == (3, 2, 6012954213)


def test_unlimited_design_and_up_to():
    # An unlimited design (u = U); G1 up to 6 delivers 49152.
    g = graph(G1)
    r = policy(g, 0, 3, 16, headroom=True).replace(INF)
    assert (
        r.plan.planned_bandwidth_per_tunnel,
        r.plan.design_by_head,
        r.plan.design_demand,
    ) == (
        2925,
        (),
        INF,
    )
    assert (r.call_delivered, r.termination) == (46812, T.EXHAUSTED)
    r = policy(g, 0, 3, 6, up_to=True, headroom=True).replace(INF)
    assert (r.plan.tunnels_per_head, r.call_delivered) == (6, 49152)


def test_up_to_reaches_the_design_exactly():
    # Three parallel links of 6, 5 and 4, D = 11, up to 3. k = 1: U = 6, F'(1) = 6. k = 2: U = 5
    # (one tunnel on the 6-link, one on the 5-link), the bounds 10 and 11 leave the score open,
    # F'(2) = 11 (z = 11 gives 6 and 5) = D: chosen. k = 3: U = 4, its lower bound min(11, 3 * 4)
    # already reaches D, so a budget F'(2) one unit low would choose 3.
    g = graph([(0, 1, 6, 1), (0, 1, 5, 1), (0, 1, 4, 1)])
    r = policy(g, 0, 1, 3, up_to=True).replace(11)
    assert (
        r.plan.tunnels_per_head,
        r.plan.max_bandwidth_per_tunnel,
        r.plan.planned_bandwidth_per_tunnel,
    ) == (2, 5, 5)
    assert (routes(r), r.call_delivered) == ([[0], [1]], 11)


def test_all_dead_and_insufficient_resolution():
    # A plan of zeros without a live head; every head dead, forwarding without tunnels.
    g = twoh(0, 0)
    for up, count in ((False, 2), (True, 1)):
        for adm in (LOSSY, LOSSLESS):
            p = policy(g, [0, 1], 2, 2 if not up else 3, adm, up_to=up)
            r = p.replace(10, design_demand=10)
            assert plan_tuple(r.plan) == (
                count,
                10,
                ((0, 5), (1, 5)),
                0,
                (),
                False,
                Fraction(0),
                0,
                0,
                0,
                0,
                None,
                0,
            )
            assert (r.program.tunnels, r.program.sources) == ((), (0, 1))
            assert r.call_delivered == 0 and p.records() == []
            if adm == LOSSY:
                assert (r.call_source_drops, r.termination) == (10, T.COMPLETE)
                assert source_rows(r) == [(0, 5, 5, 0, 5, 0), (1, 5, 5, 0, 5, 0)]
            else:
                assert (r.unadmitted, r.termination) == (10, T.EXHAUSTED)
    g = twoh(1, 100)
    for adm in (LOSSY, LOSSLESS):
        p = policy(g, [0, 1], 2, 16, adm)
        r = p.replace(200, design_demand=200)
        assert (
            r.plan.insufficient_resolution,
            r.plan.max_fractional_bandwidth_per_tunnel,
            r.plan.max_bandwidth_per_tunnel,
        ) == (True, Fraction(1, 16), 0)
        assert (r.plan.heads_with_capacity, r.program.tunnels, r.call_delivered) == (
            (0, 1),
            (),
            0,
        )
        assert p.program is not None and p.program.sources == (0, 1)
        f = p.forward(10)
        if adm == LOSSY:
            assert source_rows(r) == [
                (0, 100, 100, 0, 100, 0),
                (1, 100, 100, 0, 100, 0),
            ]
            assert (f.call_source_drops, f.termination) == (10, T.COMPLETE)
            assert source_rows(f) == [(0, 5, 5, 0, 5, 0), (1, 5, 5, 0, 5, 0)]
        else:
            assert (r.termination, f.call_delivered, f.termination) == (
                T.EXHAUSTED,
                0,
                T.EXHAUSTED,
            )
        assert f.plan == r.plan or plan_tuple(f.plan) == plan_tuple(r.plan)


def test_dead_head_rows_and_records():
    # A dead head's indices, rows and records; the round order with heads behind 3 and 8.
    g = twoh(0, 10)
    p = policy(g, [0, 1], 2, 2, LOSSY)
    r = p.replace(10, design_demand=10)
    assert [t.entry_id for t in r.program.tunnels] == [1, 3]
    assert source_rows(r) == [(0, 5, 5, 0, 5, 0), (1, 5, 5, 5, 0, 0)]
    assert sorted(rec.entry_id for rec in p.records()) == [1, 3]
    assert [(t.entry_id, t.status, t.assigned, t.delivered) for t in r.by_tunnel] == [
        (1, UP, 3, 3),
        (3, UP, 2, 2),
    ]
    r = policy(graph([(0, 2, 3, 1), (1, 2, 8, 1)]), [0, 1], 2, 2).replace(
        INF, design_demand=6
    )
    assert [
        (t.entry_id, t.head, t.local_index, list(t.edges)) for t in r.program.tunnels
    ] == [
        (0, 0, 0, [0]),
        (1, 1, 0, [1]),
        (2, 0, 1, [0]),
        (3, 1, 1, [1]),
    ]


def test_optimized_one_source():
    # OPTIMIZED with one source: every field equal to the EQUAL call's.
    g = graph(G1)
    a = policy(g, 0, 3, 16, LOSSY, headroom=True).replace(49152)
    p = policy(
        g,
        0,
        3,
        16,
        LOSSY,
        headroom=True,
        source_allocation=ngc.SourceAllocation.OPTIMIZED,
    )
    b = p.replace(49152)
    assert plan_tuple(a.plan) == plan_tuple(b.plan)
    assert (routes(a), rows(a), a.call_delivered, a.call_edge_drops) == (
        routes(b),
        rows(b),
        b.call_delivered,
        b.call_edge_drops,
    )
    assert (p.config.routing.source_allocation, b.program.source_weights) == (
        ngc.SourceAllocation.OPTIMIZED,
        None,
    )


def test_planned_only_on_t2():
    # The planned bandwidth is T2's alone: 0 on T1's tunnels, on a controller's paths (CentralizedTE,
    # ShortestPaths OPTIMIZED) and on a supplied program's tunnels; Pinned over a T2 program keeps u.
    g = twoh(2, 10)

    def planned(cfg, call):
        p = ngc.FlowPolicy(ngc.FlowLedger(g), [0, 1], 2, cfg)
        return [t.planned_bandwidth for t in call(p).program.tunnels]

    t1 = ngc.FlowPolicyConfig(
        routing=ngc.DistributedTE(ngc.TunnelCount(2)), admission=LOSSY
    )
    assert planned(t1, lambda p: p.replace(12, design_demand=12)) == [0, 0, 0, 0]
    assert planned(GLOBAL, lambda p: p.place(12)) == [0, 0]
    max_flow = ngc.FlowPolicyConfig(
        routing=ngc.ShortestPaths(allocation=ngc.RouteAllocation.OPTIMIZED)
    )
    assert planned(max_flow, lambda p: p.place(12)) == [0, 0]
    supplied = ngc.ForwardingProgram.from_routes(
        g,
        [0, 1],
        [2],
        [
            ngc.RouteDefinition(0, ngc.Path(g, [[0]])),
            ngc.RouteDefinition(1, ngc.Path(g, [[1]])),
        ],
    )
    pinned = ngc.FlowPolicyConfig(routing=ngc.Pinned(supplied), admission=LOSSY)
    assert planned(pinned, lambda p: p.replace(4)) == [0, 0]
    laid = policy(g, [0, 1], 2, 2).replace(INF).program
    assert [t.planned_bandwidth for t in laid.tunnels] == [1, 1, 1, 1]
    replay = ngc.FlowPolicyConfig(routing=ngc.Pinned(laid), admission=LOSSY)
    assert planned(replay, lambda p: p.replace(4)) == [1, 1, 1, 1]


def test_two_heads_equal_against_global():
    # Heads behind 2 and 10 against the global optimum, N = 2: the common size 1, lossless
    # bound 4; optimized CentralizedTE delivers 12 on the same network.
    g = twoh(2, 10)
    r = policy(g, [0, 1], 2, 2).replace(INF)
    assert (
        r.plan.max_bandwidth_per_tunnel,
        r.plan.planned_bandwidth_per_tunnel,
        r.call_delivered,
    ) == (1, 1, 4)
    assert (
        ngc.evaluate_demand(g, [0, 1], 2, ngc.EvaluationOptions(GLOBAL)).call_delivered
        == 12
    )


def test_empty_roster_zero_design_forward_and_flags():
    # An empty roster, a zero design, forward, and the change flags.
    empty = graph([(1, 2, 5, 1)])
    for adm in (LOSSY, LOSSLESS):
        r = policy(empty, 0, 2, 2, adm, headroom=True).replace(5, design_demand=5)
        assert (
            r.termination,
            r.program.tunnels,
            r.program.sources,
            r.program.destinations,
        ) == (T.NO_ROUTE, (), (0,), (2,))
        assert plan_tuple(r.plan) == (
            2,
            5,
            (),
            0,
            (),
            False,
            Fraction(0),
            0,
            0,
            0,
            0,
            0,
            0,
        )
    r = policy(one10(), 0, 1, 2).replace(0)
    assert (
        r.plan.planned_bandwidth_per_tunnel,
        [t.planned_bandwidth for t in r.program.tunnels],
        r.call_delivered,
        r.termination,
    ) == (
        1,
        [1, 1],
        0,
        T.COMPLETE,
    )
    g = graph(G1)
    p = policy(g, 0, 3, 16, LOSSY, headroom=True)
    first = p.replace(49152)
    f = p.forward(32768)
    assert (
        f.call_delivered,
        f.routes_changed,
        f.weights_changed,
        f.reservations_changed,
    ) == (
        32768,
        False,
        False,
        False,
    )
    assert f.program is first.program and plan_tuple(f.plan) == plan_tuple(first.plan)
    q = policy(g, 0, 3, 16, headroom=True)
    q.replace(49152)
    assert q.forward(INF).call_delivered == 46812
    second = p.replace(49152, design_demand=8192)
    assert (
        second.routes_changed,
        second.weights_changed,
        second.reservations_changed,
    ) == (False, False, False)
    assert (
        first.plan.planned_bandwidth_per_tunnel,
        second.plan.planned_bandwidth_per_tunnel,
    ) == (
        2925,
        512,
    )


def test_no_certificates_and_plan_only_for_tunnels():
    # No certificates; the plan only for the tunnel models.
    g = graph(G1)
    r = policy(g, 0, 3, 16).replace(INF, outputs=ngc.OutputOptions(min_cut=True))
    assert (r.total_cost, r.min_cut, r.marginal_cost_curve) == (None, None, None)
    pin = ngc.FlowPolicy(
        ngc.FlowLedger(g), 0, 3, ngc.FlowPolicyConfig(routing=ngc.Pinned(r.program))
    )
    assert pin.replace(10).plan is None
    # An optimized-source-allocation T2 policy whose roster is empty (source 0 reaches nothing): no
    # certificate either, never an empty cut or curve.
    empty = ngc.FlowPolicy(
        ngc.FlowLedger(graph([(1, 2, 10, 1)], 3)),
        0,
        2,
        t2_config(2, source_allocation=ngc.SourceAllocation.OPTIMIZED),
    )
    for offer, termination in ((5, T.NO_ROUTE), (0, T.COMPLETE)):
        e = empty.replace(offer, outputs=ngc.OutputOptions(min_cut=True))
        assert (e.termination, e.call_delivered, e.plan is not None) == (
            termination,
            0,
            True,
        )
        assert (e.total_cost, e.min_cut, e.marginal_cost_curve) == (None, None, None)


def test_t2_replay_with_every_tunnel_failed():
    # The G1 call's program with e0 and e1 masked: no head reaches 3; sixteen zero records.
    g = graph(G1)
    base = policy(g, 0, 3, 16, LOSSY, headroom=True)
    base.replace(49152)
    pin = ngc.FlowPolicy(
        ngc.FlowLedger(g),
        0,
        3,
        ngc.FlowPolicyConfig(routing=ngc.Pinned(base.program), admission=LOSSY),
        g.eligibility(excluded_edges=[0, 1]),
    )
    r = pin.replace(49152)
    assert (r.termination, len(pin.records()), [t.status for t in r.by_tunnel]) == (
        T.NO_ROUTE,
        16,
        [FAILED] * 16,
    )
    assert all(rec.delivered == 0 for rec in pin.records())
    assert pin.replace(0).termination == T.COMPLETE


# ---------------------------------------------------------------- zero offers (E)


def test_zero_and_masked():
    # Zero offers and an Eligibility that excludes every edge.
    g = twoh(0, 0)
    for adm in (LOSSY, LOSSLESS):
        p = policy(g, [0, 1], 2, 2, adm)
        assert p.replace(0, design_demand=10).termination == T.COMPLETE
        assert p.evaluate(0, design_demand=10).termination == T.COMPLETE
        assert p.forward(0).termination == T.COMPLETE
        r = p.replace(10, design_demand=10)
        assert r.termination == (T.COMPLETE if adm == LOSSY else T.EXHAUSTED)
        assert (r.program.tunnels, p.records()) == ((), [])
    p = policy(twoh(1, 100), [0, 1], 2, 16, LOSSY)
    assert (p.replace(0).termination, p.forward(0).termination) == (
        T.COMPLETE,
        T.COMPLETE,
    )
    assert source_rows(p.forward(10)) == [(0, 5, 5, 0, 5, 0), (1, 5, 5, 0, 5, 0)]
    for adm in (LOSSY, LOSSLESS):
        p = policy(one10(), 0, 1, 2, adm)
        r = p.replace(0)
        assert (r.termination, [t.planned_bandwidth for t in r.program.tunnels]) == (
            T.COMPLETE,
            [1, 1],
        )
        assert [(rec.entry_id, rec.delivered) for rec in p.records()] == [
            (0, 0),
            (1, 0),
        ]
        ids = [rec.flow_id for rec in p.records()]
        e = p.evaluate(0)
        assert (
            e.termination == T.COMPLETE and [rec.flow_id for rec in p.records()] == ids
        )
        f = p.forward(0)
        assert (f.termination, f.routes_changed, f.weights_changed) == (
            T.COMPLETE,
            False,
            False,
        )
        assert len(p.records()) == 2 and [rec.flow_id for rec in p.records()] != ids
    g = twoh(10, 10)
    none = g.eligibility(excluded_edges=[0, 1])
    q = ngc.EvaluationOptions(t2_config(2), eligibility=none)
    r = ngc.evaluate_demand(g, [0, 1], 2, q)
    assert (r.termination, r.call_delivered) == (T.NO_ROUTE, 0)
    assert [
        x.termination for x in ngc.evaluate_pairs(g, [([0, 1], [2]), ([0], [2])], q)
    ] == [T.NO_ROUTE] * 2


# ---------------------------------------------------------------- evaluators (Q)


def test_capacity():
    # The design rule of capacity, nothing shared between pairs, a zero offered.
    g = graph(G1)
    hr16 = t2_config(16, headroom=True)
    assert (
        ngc.evaluate_demand(g, 0, 3, ngc.EvaluationOptions(hr16)).call_delivered
        == 46812
    )
    assert (
        ngc.evaluate_demand(
            g, 0, 3, ngc.EvaluationOptions(hr16, design_demand=8192)
        ).call_delivered
        == 46812
    )
    ly = t2_config(16, LOSSY, headroom=True)
    assert (
        ngc.evaluate_demand(
            g, 0, 3, ngc.EvaluationOptions(ly, design_demand=49152)
        ).call_delivered
        == 47104
    )
    assert (
        ngc.evaluate_demand(
            g, 0, 3, ngc.EvaluationOptions(ly, offer=49152)
        ).call_delivered
        == 47104
    )
    with pytest.raises(
        ValueError,
        match=exact("CentralizedTE LOSSY evaluation requires a finite offer"),
    ):
        ngc.evaluate_demand(g, 0, 3, ngc.EvaluationOptions(ly))
    with pytest.raises(
        ValueError,
        match=exact("CentralizedTE LOSSY evaluation requires a finite offer"),
    ):
        ngc.evaluate_demand(g, 0, 3, ngc.EvaluationOptions(ly, design_demand=INF))
    # Lossless with a finite offer and no design: the design is the offer, delivered up to it.
    r5 = ngc.evaluate_demand(g, 0, 3, ngc.EvaluationOptions(hr16, offer=5))
    assert (r5.call_delivered, r5.plan.design_demand) == (5, 5)
    with pytest.raises(
        ValueError,
        match=exact("design_demand must be an integer in [0, 2^62) or permitted +inf"),
    ):
        ngc.evaluate_demand(g, 0, 3, ngc.EvaluationOptions(hr16, design_demand=-1))
    r = ngc.evaluate_demand(g, 0, 3, ngc.EvaluationOptions(ly, offer=0))
    assert (r.termination, r.call_delivered) == (T.COMPLETE, 0)
    one = one10()
    got = ngc.evaluate_pairs(
        one, [([0], [1]), ([0], [1])], ngc.EvaluationOptions(t2_config(2))
    )
    assert [x.call_delivered for x in got] == [10, 10]


@pytest.mark.parametrize("offered", [None, INF, 0, 7])
@pytest.mark.parametrize("design", [None, 0, 4, INF])
def test_lossy_offer_rule(design, offered):
    # T2 LOSSY's offer rule in all three evaluators: a finite offer is used as given, an
    # absent one is the design if finite, an explicit inf is an error whatever the design (the
    # offer's domain, checked first). One link of 10, N = 2: both tunnels share the link,
    # delivery min(offer, 10).
    g = one10()
    q = ngc.EvaluationOptions(t2_config(2, LOSSY), design_demand=design, offer=offered)
    if offered is None:
        want = design if design is not None and design != INF else None
    else:
        want = None if offered == INF else offered
    if want is None:
        text = (
            "offer must be an integer in [0, 2^62)"
            if offered == INF
            else "CentralizedTE LOSSY evaluation requires a finite offer"
        )
        for call in (
            lambda: ngc.evaluate_demand(g, 0, 1, q),
            lambda: ngc.evaluate_pairs(g, [([0], [1]), ([0], [1])], q, workers=2),
            lambda: ngc.evaluate_edge_failures(g, 0, 1, q, workers=2),
        ):
            with pytest.raises(ValueError, match=exact(text)):
                call()
        return
    assert ngc.evaluate_demand(g, 0, 1, q).call_delivered == want
    got = ngc.evaluate_pairs(g, [([0], [1]), ([0], [1])], q, workers=2)
    assert [r.call_delivered for r in got] == [want, want]
    assert (
        ngc.evaluate_edge_failures(g, 0, 1, q, workers=2).baseline.call_delivered
        == want
    )


def test_edge_failures():
    # Every allowed link is a candidate; deleting head 0's zero link
    # removes it from the roster (a negative loss).
    g = twoh(0, 10)
    s = ngc.evaluate_edge_failures(
        g, [0, 1], 2, ngc.EvaluationOptions(t2_config(2), design_demand=10)
    )
    assert (
        s.baseline.call_delivered,
        [(e.edge, e.delivery_reduction) for e in s.entries],
    ) == (
        0,
        [(0, -10), (1, 0)],
    )
    g = graph([(0, 1, 10, 1), (2, 1, 5, 1)])  # edge 1 lies on no route of source 0
    s = ngc.evaluate_edge_failures(g, 0, 1, ngc.EvaluationOptions(t2_config(2)))
    assert [(e.edge, e.delivery_reduction) for e in s.entries] == [(0, 10), (1, 0)]
    s = ngc.evaluate_edge_failures(
        g, 0, 1, ngc.EvaluationOptions(t2_config(2, LOSSY), offer=0)
    )
    assert [e.delivery_reduction for e in s.entries] == [0, 0]


# ---------------------------------------------------------------- numeric ranges (W)


def test_ranges():
    # The numeric ranges: a route cost above 2^64, the largest sizes and budgets.
    wide = graph([(0, 1, 16, 2**61 + 1)])
    assert (
        policy(wide, 0, 1, 16).replace(INF, design_demand=16).plan.total_tunnel_cost
        == 36893488147419103248
    )
    onemax = graph([(0, 1, L62 - 1, 1)])
    r = policy(onemax, 0, 1, 1, headroom=True).replace(INF)
    assert (
        r.plan.max_bandwidth_per_tunnel,
        r.plan.planned_bandwidth_per_tunnel,
        r.plan.total_planned_bandwidth,
    ) == (
        L62 - 1,
        L62 - 1,
        L62 - 1,
    )
    assert (r.plan.headroom_penalty, r.plan.max_lossless_offer, r.call_delivered) == (
        2**32,
        L62 - 1,
        L62 - 1,
    )
    fan = graph([(h, 4, 2**60 - 1, 1) for h in range(4)])
    r = policy(fan, [0, 1, 2, 3], 4, 16).replace(INF)
    assert (
        r.plan.max_bandwidth_per_tunnel,
        r.plan.total_planned_bandwidth,
        r.plan.max_lossless_offer,
    ) == (
        2**56 - 1,
        L62 - 64,
        L62 - 4,
    )
    dead = graph([(h, 101, 0, 1) for h in range(100)] + [(100, 101, 2**61, 1)])
    r = policy(dead, list(range(101)), 101, 2, up_to=True).replace(INF)
    assert (
        r.plan.tunnels_per_head,
        r.plan.max_lossless_offer,
        r.plan.heads_with_capacity,
    ) == (1, 0, (100,))
    r = policy(dead, list(range(101)), 101, 1, LOSSY).replace(101)
    assert (r.call_delivered, r.call_source_drops) == (1, 100)
    costmax = graph([(0, 1, 10000, L62 - 1)])
    r = policy(costmax, 0, 1, 10000, headroom=True).replace(INF, design_demand=10000)
    assert (r.plan.max_bandwidth_per_tunnel, r.plan.planned_bandwidth_per_tunnel) == (
        1,
        1,
    )
    assert r.plan.total_tunnel_cost == 10000 * (L62 - 1) == 46116860184273879030000
    assert (
        r.plan.headroom_penalty
        == sum(2**32 * j // 10000 for j in range(1, 10001))
        == 21476983958656
    )


# ---------------------------------------------------------------- the property test


def _random_case(rng):
    """A small random graph with a route from 0 to the last node: rows (u, v, capacity, cost)."""
    n = rng.randint(3, 5)
    raw = [
        (i, i + 1, rng.choice([0, 2, 4, 6, 8]), rng.randint(1, 4)) for i in range(n - 1)
    ]
    for _ in range(rng.randint(1, 5)):
        u, v = rng.randrange(n - 1), rng.randrange(1, n)
        if u != v:
            raw.append((u, v, rng.choice([0, 1, 2, 3, 4, 6, 8]), rng.randint(1, 4)))
    rng.shuffle(raw)
    return n, raw


def test_property():
    rng = random.Random(20261019)
    laid = optimal = bounds = lossy_cmp = upto = live_all = 0
    for trial in range(300):
        n, rows_ = _random_case(rng)
        g = build(n, rows_)
        edges = [rows_[int(i)] for i in g.input_positions()]
        sources = sorted(rng.sample(range(n - 1), rng.randint(1, min(2, n - 1))))
        dests = [n - 1]
        mask = eligibility = None
        if rng.random() < 0.3:
            xe = [rng.randrange(len(rows_))]
            mask, eligibility = ([], xe), g.eligibility(excluded_edges=xe)
        caps = [c for _, _, c, _ in edges]
        free = [rng.randint(0, c) for c in caps] if rng.random() < 0.4 else list(caps)
        hr = rng.random() < 0.5
        N = rng.randint(1, 3)
        up = rng.random() < 0.3
        D = rng.choice([None, rng.randint(0, 14)])
        fg = ngc.FlowLedger(g, free)
        p = ngc.FlowPolicy(
            fg, sources, dests, t2_config(N, headroom=hr, up_to=up), eligibility
        )
        r = p.replace(INF, design_demand=INF if D is None else D)
        want = orc.t2(
            edges, sources, dests, N, D, up_to=up, V=free, mask=mask, headroom=hr
        )
        where = (trial, rows_, sources, free, N, up, D, hr, mask)
        pl = r.plan
        assert (
            pl.tunnels_per_head,
            list(pl.heads_with_capacity),
            pl.max_bandwidth_per_tunnel,
            pl.max_fractional_bandwidth_per_tunnel,
        ) == (
            want.count,
            want.live,
            want.U,
            want.fluid,
        ), where
        assert (pl.planned_bandwidth_per_tunnel, pl.insufficient_resolution) == (
            want.u,
            want.insufficient,
        ), where
        got = [(t.entry_id, t.head, list(t.edges)) for t in r.program.tunnels]
        assert got == [(t.index, t.head, t.path) for t in want.tunnels], where
        assert (pl.total_tunnel_cost, pl.headroom_penalty if hr else 0) == (
            want.J,
            want.phi,
        ), where
        assert pl.max_lossless_offer == want.F == r.call_delivered, where
        laid += 1
        if want.U:
            # U maximal: U + 1 does not fit; J then Phi minimal over every count vector at u.
            ok = orc.allowed(edges, mask)
            k = want.count
            assert not orc.t2_fits(
                edges, ok, want.roster, want.live, dests, k, Fraction(want.U + 1), free
            ), where
            best = orc.t2_best_objectives(
                edges, ok, want.live, dests, k, want.u, free, hr
            )
            assert (want.J, want.phi) == best, where
            optimal += 1
            K = len(want.live) * k
            if want.u == want.U and want.live == want.roster:
                # The remainder bound 0 <= F - K U < K.
                assert 0 <= want.F - K * want.U < K, where
                bounds += 1
        if up:
            upto += 1
        if want.live and want.live == want.roster:
            live_all += 1
        # Lossy: the oracle's loss rule, conservation, the layout independent of the offer.
        offer = rng.randint(0, 14)
        lp = ngc.FlowPolicy(
            ngc.FlowLedger(g, free),
            sources,
            dests,
            t2_config(N, LOSSY, hr, up),
            eligibility,
        )
        lr = lp.replace(offer, design_demand=INF if D is None else D)
        assert [(t.entry_id, list(t.edges)) for t in lr.program.tunnels] == [
            (a, c) for a, _, c in got
        ], where
        assert (
            lr.offered
            == lr.source_drops + sum(x for _, x in lr.edge_drops) + lr.delivered
        ), where
        if want.roster:
            try:
                delivered, lost, per, _ = orc.lossy(
                    edges, want.roster, want.tunnels, offer, free
                )
            except AssertionError:
                continue
            assert (lr.call_delivered, lr.call_source_drops) == (delivered, lost), where
            assert [
                (t.assigned, t.delivered, t.total_edge_drops) for t in lr.by_tunnel
            ] == per, where
            lossy_cmp += 1
    assert (
        laid > 250
        and optimal > 120
        and bounds > 40
        and lossy_cmp > 150
        and upto > 50
        and live_all > 100
    ), (
        laid,
        optimal,
        bounds,
        lossy_cmp,
        upto,
        live_all,
    )


# ---------------------------------------------------------------- Python allocation failures


@needs_hooks
@pytest.mark.parametrize("method", ["replace", "evaluate", "forward"])
def test_t2_results_survive_python_allocation_failures(method):
    # Every Python object of a T2 result (the plan, its fields, the tunnels) exists before anything
    # is published. Positional calls, as for T1 (the module docstring's limit 1).
    g = graph(G1)

    def make():
        fg = ngc.FlowLedger(g)
        p = ngc.FlowPolicy(fg, 0, 3, t2_config(4, LOSSY, True))
        p.replace(4096)

        def snapshot():
            return [full_ledger_state(fg), repr(p._state())]

        return getattr(p, method), (8192,), snapshot

    assert sweep_python_allocation(make) >= 3


GETTER_SWEEP = """
import _netgraph_core as ext, netgraph_core as ngc
g = ngc.Graph(4, [0, 0, 1, 2], [1, 2, 3, 3], [8192, 40960, 8192, 40960], [1, 1, 1, 1], 1.0)
cfg = ngc.FlowPolicyConfig(routing=ngc.CentralizedTE(source_allocation=ngc.SourceAllocation.EQUAL,
                                                     tunnel_count=ngc.TunnelCount(4), prefer_headroom=True),
                           admission=ngc.Admission.LOSSY)
p = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 3, cfg)
r = p.replace(4096)
gl = ngc.CentralizedTE(source_allocation=ngc.SourceAllocation.EQUAL, tunnel_count=ngc.TunnelCount(4))
fn = lambda: {getter}
expected = repr(fn())  # the getter reads its value first: a missing name fails here
for k in range(400):
    hit, seen, value, error = ext._fail_python_allocation(fn, (), k)
    if not hit:
        assert error is None and repr(value) == expected, (k, error, value)
        break
    assert value is None and type(error) is MemoryError, (k, error)
else:
    raise SystemExit("never completed")
print("failures", k)
"""


@needs_hooks
@pytest.mark.parametrize(
    "getter",
    [
        "r.plan",
        "r.plan.max_fractional_bandwidth_per_tunnel",
        "r.plan.heads_with_capacity",
        "r.plan.total_tunnel_cost",
        "r.plan.headroom_penalty",
        "r.program.tunnels",
        "gl.tunnel_count",
        "repr(r)",
    ],
)
def test_t2_getters_survive_python_allocation_failures(getter):
    import subprocess
    import sys

    done = subprocess.run(
        [sys.executable, "-c", GETTER_SWEEP.format(getter=getter)],
        capture_output=True,
        text=True,
        timeout=300,
    )
    assert done.returncode == 0, (getter, done.returncode, done.stderr[-400:])
    assert int(done.stdout.split()[-1]) >= 1
