"""T2 (Global(tunnels=Tunnels(N)), Tunnels(L, up_to=True)) through the binding: the contract's T2
examples and calls (7.1 to 7.9, 8.5), the configuration and call checks, sizing, objectives,
count selection, operations, programs and replays, records, flags and terminations, the
evaluators, the numeric ranges of T2's sums, zero offers and everything masked, a property test
against the independent tunnel_oracle.py, and the Python-allocation sweeps of T2's results.
Expected values are the contract's or derived from it by hand, never taken from a build."""

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
HEAD = ngc.EdgePriority.HEADROOM


def graph(rows, n=None):
    """A graph whose rows are in EdgeId order (asserted): (u, v, capacity, cost)."""
    n = n if n is not None else 1 + max(max(u, v) for u, v, _, _ in rows)
    g = build(n, rows)
    assert [int(i) for i in g.input_index()] == list(range(len(rows)))
    return g


def t2_config(count, admission=None, headroom=False, up_to=False, origin=None):
    """A T2 FlowPolicyConfig: Global(tunnels=Tunnels(count, up_to=up_to),
    prefer_headroom=headroom), Lossless by default."""
    return ngc.FlowPolicyConfig(
        routing=ngc.Global(
            tunnels=ngc.Tunnels(count, up_to=up_to), prefer_headroom=headroom
        ),
        admission=admission or LOSSLESS,
        origin=origin,
    )


def policy(
    g,
    s,
    d,
    n,
    admission=None,
    ledger=None,
    elig=None,
    headroom=False,
    up_to=False,
    origin=None,
):
    ledger = ledger if ledger is not None else ngc.FlowLedger(g)
    return ngc.FlowPolicy(
        ledger, s, d, t2_config(n, admission, headroom, up_to, origin), elig
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
        (t.index, t.head, t.status, t.assigned, t.delivered, t.dropped)
        for t in r.by_tunnel
    ]


def source_rows(r):
    return [
        (s.source, s.offered, s.admitted, s.delivered, s.source_loss, s.dropped)
        for s in r.by_source
    ]


def plan_tuple(p):
    return (
        p.count,
        p.design_demand,
        p.head_design,
        p.reserved_total,
        p.live_heads,
        p.insufficient_resolution,
        p.fluid_share,
        p.common_share,
        p.planned_per_tunnel,
        p.planned_total,
        p.route_cost,
        p.headroom_penalty,
        p.forwarding_bound,
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


# ---------------------------------------------------------------- the contract's examples (X)

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
        p.common_share,
        p.planned_per_tunnel,
        p.fluid_share,
        p.forwarding_bound,
    ) == (U, u, fluid, F)
    assert (p.planned_total, p.route_cost, p.live_heads) == (total, cost, live)
    if phi is not None:
        assert p.headroom_penalty == phi
    assert routes(r) == rts
    if indices is not None:
        assert [t.index for t in r.program.tunnels] == indices
    assert all((t.planned, t.weight, t.status) == (u, 1, UP) for t in r.program.tunnels)
    assert (r.placed, r.termination) == (F, T.EXHAUSTED)  # the unlimited query admits F
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
    assert (s.placed, s.call_source_loss, s.termination) == (delivery, loss, T.COMPLETE)
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
    assert (r.plan.count, routes(r), r.placed) == (k, rts, delivered)
    if bounds is None:
        return
    # With every head live, F'(k) is EXACT(k)'s bound at U(k) (an unlimited design).
    for c, b in enumerate(bounds, 1):
        assert policy(g, heads, dst, c, headroom=True).replace(INF).placed == b, c


def test_call_of_section_7_9():
    # The contract's runnable call, then its program's fields and the ledger.
    q = 2.0**-12
    g = ngc.Graph(
        4,
        src=[0, 0, 1, 2],
        dst=[1, 2, 3, 3],
        capacity_units=ngc.to_units([2, 10, 2, 10], q),
        cost=[1, 1, 1, 1],
        resolution=q,
    )
    te16 = ngc.FlowPolicyConfig(
        routing=ngc.Global(tunnels=ngc.Tunnels(16), prefer_headroom=True),
        admission=ngc.Admission.LOSSY,
    )
    fg = ngc.FlowLedger(g)
    p = ngc.FlowPolicy(fg, 0, 3, te16)
    r = p.replace(49152)
    assert (
        r.plan.common_share,
        r.plan.planned_per_tunnel,
        r.plan.forwarding_bound,
    ) == (2925, 2925, 46812)
    assert r.placed == 47104
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
        (t.requested, t.reserved, t.planned, t.status, t.weight)
        for t in r.program.tunnels
    ] == [(0, 0, 2925, UP, 1)] * 16
    assert [t.delivered for t in r.by_tunnel] == [3072] * 2 + [2926] * 10 + [2925] * 4
    assert r.call_drops == [(1, 2048)]
    assert r.reservations_changed is False and (
        r.routes_changed,
        r.weights_changed,
    ) == (True, True)
    # A plan books nothing.
    assert fg.reserved().tolist() == [0, 0, 0, 0]
    assert all(rec.reserved_edges == [] for rec in p.records())
    assert ngc.Global(tunnels=ngc.Tunnels(256, up_to=True)).tunnels == ngc.Tunnels(
        256, up_to=True
    )


def test_failure_replay_of_section_8_5():
    # The replay of 8.5 under the mask, then the same program stale (black holes).
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
        elig=fail,
    )
    r = replay.replace(49152)
    assert [t.status for t in r.by_tunnel][:3] == [FAILED, FAILED, UP]
    assert r.placed == 40960 and r.plan is None
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
    assert (s.call_drops, s.placed) == ([(1, 2048), (2, 6144)], 40960)
    assert [t.status for t in base.tunnels] == [
        UP
    ] * 16  # the program's own statuses stay


# ---------------------------------------------------------------- validation (V)


def test_matrix():
    # T2's row of the matrix and the Global selection rule, in order (the evaluators run the same
    # checks).
    g = twoh(10, 10)
    fg = ngc.FlowLedger(g)
    t2 = ngc.Global(tunnels=ngc.Tunnels(2))
    sel = ngc.PathSelection
    dag = dag_of(g, [0, 1])
    both = ngc.Global(support=dag, tunnels=ngc.Tunnels(2), prefer_headroom=True)
    no_sel = "Global does not accept selection"
    cases = [
        (ngc.FlowPolicyConfig(routing=t2, placement=ngc.FlowPlacement.EQUAL_BALANCED), "T2 does not accept placement"),
        (ngc.FlowPolicyConfig(routing=t2, weights=ngc.WeightBasis.CAPACITY), "T2 does not accept weights"),
        (ngc.FlowPolicyConfig(routing=t2, edge_weights=[1, 1]), "T2 does not accept edge_weights"),
        # Every Global rejects any selection (the default and HEADROOM included) before support and
        # tunnels, and so before placement; prefer_headroom without tunnels right after them.
        (ngc.FlowPolicyConfig(routing=t2, selection=sel()), no_sel),
        (ngc.FlowPolicyConfig(routing=t2, selection=sel(priority=HEAD)), no_sel),
        (ngc.FlowPolicyConfig(routing=ngc.Global(), selection=sel()), no_sel),
        (ngc.FlowPolicyConfig(routing=t2, placement=ngc.FlowPlacement.EQUAL_BALANCED, selection=sel()), no_sel),
        (ngc.FlowPolicyConfig(routing=both, selection=sel()), no_sel),
        (ngc.FlowPolicyConfig(routing=both), "Global support and tunnels are mutually exclusive"),
        (ngc.FlowPolicyConfig(routing=both, placement=ngc.FlowPlacement.EQUAL_BALANCED),
         "Global support and tunnels are mutually exclusive"),
        (ngc.FlowPolicyConfig(routing=ngc.Global(prefer_headroom=True)), "prefer_headroom requires tunnels"),
        (ngc.FlowPolicyConfig(routing=ngc.Global(support=dag, prefer_headroom=True),
                              placement=ngc.FlowPlacement.EQUAL_BALANCED),
         "prefer_headroom requires tunnels"),
        (ngc.FlowPolicyConfig(routing=t2, origin=ngc.Origin.OPTIMIZED), "T2 with several sources requires origin=EQUAL"),
    ]  # fmt: skip
    for cfg, text in cases:
        with pytest.raises(ValueError, match=exact(text)):
            ngc.FlowPolicy(fg, [0, 1], 2, cfg)
        with pytest.raises(ValueError, match=exact(text)):
            ngc.capacity(g, [0, 1], 2, ngc.QueryOptions(cfg))
    with pytest.raises(ValueError, match=exact("T2 does not use DagCache")):
        ngc.FlowPolicy(
            fg, [0, 1], 2, ngc.FlowPolicyConfig(routing=t2), cache=ngc.DagCache(g)
        )
    # The several-sources rule counts the normalized input sources, before the mask; one source
    # with OPTIMIZED is accepted (the same calculation).
    opt = ngc.FlowPolicyConfig(routing=t2, origin=ngc.Origin.OPTIMIZED)
    with pytest.raises(
        ValueError, match=exact("T2 with several sources requires origin=EQUAL")
    ):
        ngc.FlowPolicy(fg, [0, 1], 2, opt, g.eligibility(excluded_nodes=[1]))
    assert ngc.FlowPolicy(fg, [0, 0], 2, opt).config.origin == ngc.Origin.OPTIMIZED
    # Accepted: both admissions, with and without prefer_headroom; the snapshot keeps no selection.
    for adm in (LOSSLESS, LOSSY):
        for hr in (False, True):
            ngc.FlowPolicy(fg, [0, 1], 2, t2_config(2, adm, hr))
    got = ngc.FlowPolicy(fg, [0, 1], 2, ngc.FlowPolicyConfig(routing=t2)).config
    assert (got.selection, got.origin, got.placement, got.routing.tunnels) == (
        None,
        ngc.Origin.EQUAL,
        None,
        ngc.Tunnels(2),
    )
    # The arguments of Global itself; prefer_headroom is False by default.
    with pytest.raises(TypeError):
        ngc.Global(tunnels=2)
    with pytest.raises(TypeError):
        ngc.Global(tunnels=ngc.Tunnels(2), prefer_headroom="yes")
    assert (ngc.Global().tunnels, ngc.Global().prefer_headroom) == (None, False)
    assert ngc.Global(prefer_headroom=True).prefer_headroom is True
    assert (
        repr(ngc.Global(tunnels=ngc.Tunnels(256, up_to=True)))
        == "Global(tunnels=Tunnels(256, up_to=True))"
    )
    assert (
        repr(ngc.Global(tunnels=ngc.Tunnels(16), prefer_headroom=True))
        == "Global(tunnels=Tunnels(16), prefer_headroom=True)"
    )
    assert repr(ngc.Global(prefer_headroom=True)) == "Global(prefer_headroom=True)"
    assert repr(ngc.Global()) == "Global()"
    # Every set field, the invalid support-and-tunnels combination included.
    dag = dag_of(twoh(10, 10), [0, 1])
    assert repr(ngc.Global(support=dag)) == "Global(support=...)"
    assert (
        repr(ngc.Global(support=dag, tunnels=ngc.Tunnels(2)))
        == "Global(support=..., tunnels=Tunnels(2))"
    )
    assert (
        repr(ngc.Global(support=dag, prefer_headroom=True))
        == "Global(support=..., prefer_headroom=True)"
    )


def test_global_selection_and_prefer_headroom_texts():
    # A Global with a selection: the contract's seven cases in its order, at a policy and at an
    # evaluator. Global with any selection first; then support and tunnels; then prefer_headroom
    # without tunnels.
    g = twoh(10, 10)
    fg = ngc.FlowLedger(g)
    dag = dag_of(g, [0, 1])
    t2 = ngc.Tunnels(2)
    cases = [
        (dict(tunnels=t2), ngc.PathSelection(), "Global does not accept selection"),
        (dict(tunnels=t2), ngc.PathSelection(priority=HEAD), "Global does not accept selection"),
        ({}, ngc.PathSelection(), "Global does not accept selection"),
        (dict(prefer_headroom=True), None, "prefer_headroom requires tunnels"),
        (dict(prefer_headroom=True, support=dag, tunnels=t2), None,
         "Global support and tunnels are mutually exclusive"),
        (dict(prefer_headroom=True, tunnels=t2), None, None),
        ({}, None, None),
    ]  # fmt: skip
    for kwargs, sel, text in cases:
        cfg = ngc.FlowPolicyConfig(routing=ngc.Global(**kwargs), selection=sel)
        if text is None:
            ngc.FlowPolicy(fg, 0, 2, cfg)
            ngc.capacity(g, 0, 2, ngc.QueryOptions(cfg))
            continue
        with pytest.raises(ValueError, match=exact(text)):
            ngc.FlowPolicy(fg, 0, 2, cfg)
        with pytest.raises(ValueError, match=exact(text)):
            ngc.capacity(g, 0, 2, ngc.QueryOptions(cfg))
    # No old spelling of the objective is accepted.
    with pytest.raises(TypeError):
        ngc.Global(tunnels=t2, priority=HEAD)


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
    assert (p.common_share, p.planned_per_tunnel, p.planned_total, p.route_cost) == (
        10,
        1,
        2,
        2,
    )
    assert routes(r) == ([[0], [1]] if on else [[0], [0]])
    assert p.headroom_penalty == (858993458 if on else None)
    assert (p.forwarding_bound, r.placed) == ((20, 12) if on else (10, 10))


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


def test_operations_and_design_validation():
    # The operations of T2 and the design's domain, with their texts.
    g = one10()
    p = policy(g, 0, 1, 2)
    for v in (5, 0):
        with pytest.raises(
            ValueError, match=exact("T2 does not support incremental placement")
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
                f"{'CostOnly' if cfg is COSTONLY else 'Global'} does not accept design_demand"
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
        r.plan.forwarding_bound,
        r.placed,
        r.plan.headroom_penalty,
    ) == ([2, 14, 2, 14], 46812, 47104, None)
    r = policy(g, 0, 3, 16).replace(INF, design_demand=8192)
    assert (r.plan.planned_per_tunnel, counts(r, 4), r.plan.forwarding_bound) == (
        512,
        [16, 0, 16, 0],
        8192,
    )
    r = policy(g, 0, 3, 16, headroom=True).replace(INF, design_demand=8192)
    assert (counts(r, 4), r.plan.forwarding_bound) == ([2, 14, 2, 14], 46812)
    p412 = graph([(0, 1, 4, 1), (0, 1, 12, 1)])
    r = policy(p412, 0, 1, 2, headroom=True).replace(INF, design_demand=2)
    assert (
        r.plan.common_share,
        r.plan.planned_per_tunnel,
        routes(r),
        r.plan.forwarding_bound,
    ) == (6, 1, [[1], [1]], 12)
    assert r.plan.headroom_penalty == 1073741823
    r = policy(p412, 0, 1, 2).replace(INF, design_demand=2)
    assert (routes(r), r.plan.forwarding_bound) == ([[0], [0]], 4)


# The graphs of the three tie rules of 7.3 (u, v, capacity, cost), each derived in tests/cpp's twin.
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
    # A reverse arc: the second unit cancels the first unit's middle link (7.3: reverse arcs undo an
    # earlier choice); a preload: the slot prices read the static capacity, not the free one.
    g = graph([(0, 1, 1, 1), (0, 2, 1, 5), (1, 2, 1, 1), (1, 3, 1, 5), (2, 3, 1, 1)])
    r = policy(g, 0, 3, 2).replace(INF)
    assert (r.plan.common_share, routes(r), r.plan.route_cost) == (
        1,
        [[0, 3], [1, 4]],
        12,
    )
    r = policy(
        one10(), 0, 1, 2, ledger=ngc.FlowLedger(one10(), [6]), headroom=True
    ).replace(INF, design_demand=4)
    assert (
        r.plan.common_share,
        r.plan.planned_per_tunnel,
        r.plan.headroom_penalty,
    ) == (3, 2, 6012954213)


def test_unlimited_design_and_up_to():
    # An unlimited design (u = U); G1 up to 6 delivers 49152.
    g = graph(G1)
    r = policy(g, 0, 3, 16, headroom=True).replace(INF)
    assert (r.plan.planned_per_tunnel, r.plan.head_design, r.plan.design_demand) == (
        2925,
        (),
        INF,
    )
    assert (r.placed, r.termination) == (46812, T.EXHAUSTED)
    r = policy(g, 0, 3, 6, up_to=True, headroom=True).replace(INF)
    assert (r.plan.count, r.placed) == (6, 49152)


def test_up_to_reaches_the_design_exactly():
    # Three parallel links of 6, 5 and 4, D = 11, up to 3. k = 1: U = 6, F'(1) = 6. k = 2: U = 5
    # (one tunnel on the 6-link, one on the 5-link), the bounds 10 and 11 leave the score open,
    # F'(2) = 11 (z = 11 gives 6 and 5) = D: chosen. k = 3: U = 4, its lower bound min(11, 3 * 4)
    # already reaches D, so a budget F'(2) one unit low would choose 3.
    g = graph([(0, 1, 6, 1), (0, 1, 5, 1), (0, 1, 4, 1)])
    r = policy(g, 0, 1, 3, up_to=True).replace(11)
    assert (r.plan.count, r.plan.common_share, r.plan.planned_per_tunnel) == (2, 5, 5)
    assert (routes(r), r.placed) == ([[0], [1]], 11)


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
            assert r.placed == 0 and p.records() == []
            if adm == LOSSY:
                assert (r.call_source_loss, r.termination) == (10, T.COMPLETE)
                assert source_rows(r) == [(0, 5, 5, 0, 5, 0), (1, 5, 5, 0, 5, 0)]
            else:
                assert (r.unadmitted, r.termination) == (10, T.EXHAUSTED)
    g = twoh(1, 100)
    for adm in (LOSSY, LOSSLESS):
        p = policy(g, [0, 1], 2, 16, adm)
        r = p.replace(200, design_demand=200)
        assert (
            r.plan.insufficient_resolution,
            r.plan.fluid_share,
            r.plan.common_share,
        ) == (True, Fraction(1, 16), 0)
        assert (r.plan.live_heads, r.program.tunnels, r.placed) == ((0, 1), (), 0)
        assert p.program is not None and p.program.sources == (0, 1)
        f = p.forward(10)
        if adm == LOSSY:
            assert source_rows(r) == [
                (0, 100, 100, 0, 100, 0),
                (1, 100, 100, 0, 100, 0),
            ]
            assert (f.call_source_loss, f.termination) == (10, T.COMPLETE)
            assert source_rows(f) == [(0, 5, 5, 0, 5, 0), (1, 5, 5, 0, 5, 0)]
        else:
            assert (r.termination, f.placed, f.termination) == (
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
    assert [t.index for t in r.program.tunnels] == [1, 3]
    assert source_rows(r) == [(0, 5, 5, 0, 5, 0), (1, 5, 5, 5, 0, 0)]
    assert sorted(rec.entry for rec in p.records()) == [1, 3]
    assert [(t.index, t.status, t.assigned, t.delivered) for t in r.by_tunnel] == [
        (1, UP, 3, 3),
        (3, UP, 2, 2),
    ]
    r = policy(graph([(0, 2, 3, 1), (1, 2, 8, 1)]), [0, 1], 2, 2).replace(
        INF, design_demand=6
    )
    assert [
        (t.index, t.head, t.local_index, list(t.edges)) for t in r.program.tunnels
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
    p = policy(g, 0, 3, 16, LOSSY, headroom=True, origin=ngc.Origin.OPTIMIZED)
    b = p.replace(49152)
    assert plan_tuple(a.plan) == plan_tuple(b.plan)
    assert (routes(a), rows(a), a.placed, a.call_drops) == (
        routes(b),
        rows(b),
        b.placed,
        b.call_drops,
    )
    assert (p.config.origin, b.program.source_weights) == (ngc.Origin.OPTIMIZED, None)


def test_planned_only_on_t2():
    # The planned bandwidth is T2's alone: 0 on T1's tunnels, on a controller's paths (Global,
    # CostOnly MAX_FLOW) and on a supplied program's tunnels; Pinned over a T2 program keeps u.
    g = twoh(2, 10)

    def planned(cfg, call):
        p = ngc.FlowPolicy(ngc.FlowLedger(g), [0, 1], 2, cfg)
        return [t.planned for t in call(p).program.tunnels]

    t1 = ngc.FlowPolicyConfig(
        routing=ngc.CapacityAware(ngc.Tunnels(2)), admission=LOSSY
    )
    assert planned(t1, lambda p: p.replace(12, design_demand=12)) == [0, 0, 0, 0]
    assert planned(GLOBAL, lambda p: p.place(12)) == [0, 0]
    max_flow = ngc.FlowPolicyConfig(
        routing=ngc.CostOnly(), placement=ngc.FlowPlacement.MAX_FLOW
    )
    assert planned(max_flow, lambda p: p.place(12)) == [0, 0]
    supplied = ngc.ForwardingProgram.from_routes(
        g,
        [0, 1],
        [2],
        [ngc.Route(0, ngc.Path(g, [[0]])), ngc.Route(1, ngc.Path(g, [[1]]))],
    )
    pinned = ngc.FlowPolicyConfig(routing=ngc.Pinned(supplied), admission=LOSSY)
    assert planned(pinned, lambda p: p.replace(4)) == [0, 0]
    laid = policy(g, [0, 1], 2, 2).replace(INF).program
    assert [t.planned for t in laid.tunnels] == [1, 1, 1, 1]
    replay = ngc.FlowPolicyConfig(routing=ngc.Pinned(laid), admission=LOSSY)
    assert planned(replay, lambda p: p.replace(4)) == [1, 1, 1, 1]


def test_two_heads_equal_against_global():
    # Heads behind 2 and 10 against the global optimum, N = 2: the common size 1, lossless
    # bound 4; optimized Global delivers 12 on the same network.
    g = twoh(2, 10)
    r = policy(g, [0, 1], 2, 2).replace(INF)
    assert (r.plan.common_share, r.plan.planned_per_tunnel, r.placed) == (1, 1, 4)
    assert ngc.capacity(g, [0, 1], 2, ngc.QueryOptions(GLOBAL)).placed == 12


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
        r.plan.planned_per_tunnel,
        [t.planned for t in r.program.tunnels],
        r.placed,
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
    assert (f.placed, f.routes_changed, f.weights_changed, f.reservations_changed) == (
        32768,
        False,
        False,
        False,
    )
    assert f.program is first.program and plan_tuple(f.plan) == plan_tuple(first.plan)
    q = policy(g, 0, 3, 16, headroom=True)
    q.replace(49152)
    assert q.forward(INF).placed == 46812
    second = p.replace(49152, design_demand=8192)
    assert (
        second.routes_changed,
        second.weights_changed,
        second.reservations_changed,
    ) == (False, False, False)
    assert (first.plan.planned_per_tunnel, second.plan.planned_per_tunnel) == (
        2925,
        512,
    )


def test_no_certificates_and_plan_only_for_tunnels():
    # No certificates; the plan only for the tunnel models.
    g = graph(G1)
    r = policy(g, 0, 3, 16).replace(INF, outputs=ngc.OutputOptions(min_cut=True))
    assert (r.total_cost, r.min_cut, r.cost_distribution) == (None, None, None)
    pin = ngc.FlowPolicy(
        ngc.FlowLedger(g), 0, 3, ngc.FlowPolicyConfig(routing=ngc.Pinned(r.program))
    )
    assert pin.replace(10).plan is None
    # An optimized-origin T2 policy whose roster is empty (source 0 reaches nothing): no
    # certificate either, never an empty cut or curve.
    empty = ngc.FlowPolicy(
        ngc.FlowLedger(graph([(1, 2, 10, 1)], 3)),
        0,
        2,
        t2_config(2, origin=ngc.Origin.OPTIMIZED),
    )
    for offer, termination in ((5, T.NO_ROUTE), (0, T.COMPLETE)):
        e = empty.replace(offer, outputs=ngc.OutputOptions(min_cut=True))
        assert (e.termination, e.placed, e.plan is not None) == (termination, 0, True)
        assert (e.total_cost, e.min_cut, e.cost_distribution) == (None, None, None)


def test_t2_replay_with_every_tunnel_failed():
    # The 7.9 call's program with e0 and e1 masked: no head reaches 3; sixteen zero records.
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
        assert (r.termination, [t.planned for t in r.program.tunnels]) == (
            T.COMPLETE,
            [1, 1],
        )
        assert [(rec.entry, rec.delivered) for rec in p.records()] == [(0, 0), (1, 0)]
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
    q = ngc.QueryOptions(t2_config(2), elig=none)
    r = ngc.capacity(g, [0, 1], 2, q)
    assert (r.termination, r.placed) == (T.NO_ROUTE, 0)
    assert [
        x.termination for x in ngc.batch_capacity(g, [([0, 1], [2]), ([0], [2])], q)
    ] == [T.NO_ROUTE] * 2


# ---------------------------------------------------------------- evaluators (Q)


def test_capacity():
    # The design rule of capacity, nothing shared between pairs, a zero offered.
    g = graph(G1)
    hr16 = t2_config(16, headroom=True)
    assert ngc.capacity(g, 0, 3, ngc.QueryOptions(hr16)).placed == 46812
    assert (
        ngc.capacity(g, 0, 3, ngc.QueryOptions(hr16, design_demand=8192)).placed
        == 46812
    )
    ly = t2_config(16, LOSSY, headroom=True)
    assert (
        ngc.capacity(g, 0, 3, ngc.QueryOptions(ly, design_demand=49152)).placed == 47104
    )
    assert ngc.capacity(g, 0, 3, ngc.QueryOptions(ly, offered=49152)).placed == 47104
    with pytest.raises(
        ValueError, match=exact("T2 LOSSY capacity requires finite offered")
    ):
        ngc.capacity(g, 0, 3, ngc.QueryOptions(ly))
    with pytest.raises(
        ValueError, match=exact("T2 LOSSY capacity requires finite offered")
    ):
        ngc.capacity(g, 0, 3, ngc.QueryOptions(ly, design_demand=INF))
    with pytest.raises(
        ValueError, match=exact("T2 LOSSLESS capacity does not accept offered")
    ):
        ngc.capacity(g, 0, 3, ngc.QueryOptions(hr16, offered=5))
    with pytest.raises(
        ValueError,
        match=exact("design_demand must be an integer in [0, 2^62) or permitted +inf"),
    ):
        ngc.capacity(g, 0, 3, ngc.QueryOptions(hr16, design_demand=-1))
    r = ngc.capacity(g, 0, 3, ngc.QueryOptions(ly, offered=0))
    assert (r.termination, r.placed) == (T.COMPLETE, 0)
    one = one10()
    got = ngc.batch_capacity(
        one, [([0], [1]), ([0], [1])], ngc.QueryOptions(t2_config(2))
    )
    assert [x.placed for x in got] == [10, 10]


@pytest.mark.parametrize("offered", [None, INF, 0, 7])
@pytest.mark.parametrize("design", [None, 0, 4, INF])
def test_lossy_offer_rule(design, offered):
    # T2 LOSSY's offer rule in all three evaluators: a finite offer is used as given, an
    # absent one is the design if finite, an explicit inf is an error whatever the design.
    # One link of 10, N = 2: both tunnels share the link, delivery min(offer, 10).
    g = one10()
    q = ngc.QueryOptions(t2_config(2, LOSSY), design_demand=design, offered=offered)
    if offered is None:
        want = design if design is not None and design != INF else None
    else:
        want = None if offered == INF else offered
    if want is None:
        for call in (
            lambda: ngc.capacity(g, 0, 1, q),
            lambda: ngc.batch_capacity(g, [([0], [1]), ([0], [1])], q, workers=2),
            lambda: ngc.sensitivity(g, 0, 1, q, workers=2),
        ):
            with pytest.raises(
                ValueError, match=exact("T2 LOSSY capacity requires finite offered")
            ):
                call()
        return
    assert ngc.capacity(g, 0, 1, q).placed == want
    got = ngc.batch_capacity(g, [([0], [1]), ([0], [1])], q, workers=2)
    assert [r.placed for r in got] == [want, want]
    assert ngc.sensitivity(g, 0, 1, q, workers=2).baseline.placed == want


def test_sensitivity():
    # Every allowed link is a candidate; deleting head 0's zero link
    # removes it from the roster (a negative loss).
    g = twoh(0, 10)
    s = ngc.sensitivity(g, [0, 1], 2, ngc.QueryOptions(t2_config(2), design_demand=10))
    assert (s.baseline.placed, [(e.edge, e.loss) for e in s.entries]) == (
        0,
        [(0, -10), (1, 0)],
    )
    g = graph([(0, 1, 10, 1), (2, 1, 5, 1)])  # edge 1 lies on no route of source 0
    s = ngc.sensitivity(g, 0, 1, ngc.QueryOptions(t2_config(2)))
    assert [(e.edge, e.loss) for e in s.entries] == [(0, 10), (1, 0)]
    s = ngc.sensitivity(g, 0, 1, ngc.QueryOptions(t2_config(2, LOSSY), offered=0))
    assert [e.loss for e in s.entries] == [0, 0]


# ---------------------------------------------------------------- numeric ranges (W)


def test_ranges():
    # The numeric ranges: a route cost above 2^64, the largest sizes and budgets.
    wide = graph([(0, 1, 16, 2**61 + 1)])
    assert (
        policy(wide, 0, 1, 16).replace(INF, design_demand=16).plan.route_cost
        == 36893488147419103248
    )
    onemax = graph([(0, 1, L62 - 1, 1)])
    r = policy(onemax, 0, 1, 1, headroom=True).replace(INF)
    assert (r.plan.common_share, r.plan.planned_per_tunnel, r.plan.planned_total) == (
        L62 - 1,
        L62 - 1,
        L62 - 1,
    )
    assert (r.plan.headroom_penalty, r.plan.forwarding_bound, r.placed) == (
        2**32,
        L62 - 1,
        L62 - 1,
    )
    fan = graph([(h, 4, 2**60 - 1, 1) for h in range(4)])
    r = policy(fan, [0, 1, 2, 3], 4, 16).replace(INF)
    assert (r.plan.common_share, r.plan.planned_total, r.plan.forwarding_bound) == (
        2**56 - 1,
        L62 - 64,
        L62 - 4,
    )
    dead = graph([(h, 101, 0, 1) for h in range(100)] + [(100, 101, 2**61, 1)])
    r = policy(dead, list(range(101)), 101, 2, up_to=True).replace(INF)
    assert (r.plan.count, r.plan.forwarding_bound, r.plan.live_heads) == (1, 0, (100,))
    r = policy(dead, list(range(101)), 101, 1, LOSSY).replace(101)
    assert (r.placed, r.call_source_loss) == (1, 100)
    costmax = graph([(0, 1, 10000, L62 - 1)])
    r = policy(costmax, 0, 1, 10000, headroom=True).replace(INF, design_demand=10000)
    assert (r.plan.common_share, r.plan.planned_per_tunnel) == (1, 1)
    assert r.plan.route_cost == 10000 * (L62 - 1) == 46116860184273879030000
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
        edges = [rows_[int(i)] for i in g.input_index()]
        sources = sorted(rng.sample(range(n - 1), rng.randint(1, min(2, n - 1))))
        dests = [n - 1]
        mask = elig = None
        if rng.random() < 0.3:
            xe = [rng.randrange(len(rows_))]
            mask, elig = ([], xe), g.eligibility(excluded_edges=xe)
        caps = [c for _, _, c, _ in edges]
        free = [rng.randint(0, c) for c in caps] if rng.random() < 0.4 else list(caps)
        hr = rng.random() < 0.5
        N = rng.randint(1, 3)
        up = rng.random() < 0.3
        D = rng.choice([None, rng.randint(0, 14)])
        fg = ngc.FlowLedger(g, free)
        p = ngc.FlowPolicy(
            fg, sources, dests, t2_config(N, headroom=hr, up_to=up), elig
        )
        r = p.replace(INF, design_demand=INF if D is None else D)
        want = orc.t2(
            edges, sources, dests, N, D, up_to=up, V=free, mask=mask, headroom=hr
        )
        where = (trial, rows_, sources, free, N, up, D, hr, mask)
        pl = r.plan
        assert (pl.count, list(pl.live_heads), pl.common_share, pl.fluid_share) == (
            want.count,
            want.live,
            want.U,
            want.fluid,
        ), where
        assert (pl.planned_per_tunnel, pl.insufficient_resolution) == (
            want.u,
            want.insufficient,
        ), where
        got = [(t.index, t.head, list(t.edges)) for t in r.program.tunnels]
        assert got == [(t.index, t.head, t.path) for t in want.tunnels], where
        assert (pl.route_cost, pl.headroom_penalty if hr else 0) == (
            want.J,
            want.phi,
        ), where
        assert pl.forwarding_bound == want.F == r.placed, where
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
                assert 0 <= want.F - K * want.U < K, where  # 7.6's remainder bound
                bounds += 1
        if up:
            upto += 1
        if want.live and want.live == want.roster:
            live_all += 1
        # Lossy: the oracle's loss rule, conservation, the layout independent of the offer.
        offer = rng.randint(0, 14)
        lp = ngc.FlowPolicy(
            ngc.FlowLedger(g, free), sources, dests, t2_config(N, LOSSY, hr, up), elig
        )
        lr = lp.replace(offer, design_demand=INF if D is None else D)
        assert [(t.index, list(t.edges)) for t in lr.program.tunnels] == [
            (a, c) for a, _, c in got
        ], where
        assert (
            lr.offered == lr.source_loss + sum(x for _, x in lr.drops) + lr.delivered
        ), where
        if want.roster:
            try:
                delivered, lost, per, _ = orc.lossy(
                    edges, want.roster, want.tunnels, offer, free
                )
            except AssertionError:
                continue
            assert (lr.placed, lr.call_source_loss) == (delivered, lost), where
            assert [
                (t.assigned, t.delivered, t.dropped) for t in lr.by_tunnel
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
cfg = ngc.FlowPolicyConfig(routing=ngc.Global(tunnels=ngc.Tunnels(4), prefer_headroom=True),
                           admission=ngc.Admission.LOSSY)
p = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 3, cfg)
r = p.replace(4096)
gl = ngc.Global(tunnels=ngc.Tunnels(4))
fn = lambda: {getter}
for k in range(400):
    hit, seen, value, error = ext._fail_python_allocation(fn, (), k)
    if not hit:
        break
    assert value is None and isinstance(error, (MemoryError, RuntimeError)), (k, error)
else:
    raise SystemExit("never completed")
print("failures", k)
"""


@needs_hooks
@pytest.mark.parametrize(
    "getter",
    [
        "r.plan",
        "r.plan.fluid_share",
        "r.plan.live_heads",
        "r.plan.route_cost",
        "r.plan.headroom_penalty",
        "r.program.tunnels",
        "gl.tunnels",
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
