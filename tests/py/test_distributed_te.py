"""T1 (DistributedTE(TunnelCount(N))) through the binding: the worked T1 examples, the configuration and call checks, setup, forwarding, operations, programs and replays,
records, flags and terminations, the evaluators, the numeric ranges, zero offers and everything
masked, a property test against the independent tunnel_oracle.py, and the Python-allocation
sweeps of T1's results. The ledger's two projections are test_ledger_projections.py's."""

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
    exact,
    full_ledger_state,
    ledger_state,
    sweep_python_allocation,
    t1_config,
)
from private_hooks import HOOKS, needs_hooks

import netgraph_core as ngc

T = ngc.Termination
S = ngc.TunnelStatus
UP, NB, CC, FAILED, BLACKHOLED = (
    S.UP,
    S.INSUFFICIENT_BANDWIDTH,
    S.COST_CEILING_EXCEEDED,
    S.FAILED,
    S.BLACKHOLED,
)
LOSSY = ngc.Admission.LOSSY
LOSSLESS = ngc.Admission.LOSSLESS
TIE = ngc.RouteTieBreak
HOLDS = ngc.OutputOptions(edge_reservations=True)


def ledger_now(fg):
    """The ledger's public state; with the private hooks also its next FlowId and owners."""
    return full_ledger_state(fg) if HOOKS else ledger_state(fg)


def graph(rows, n=None):
    """A graph whose rows are in EdgeId order (asserted): (u, v, capacity, cost)."""
    n = n if n is not None else 1 + max(max(u, v) for u, v, _, _ in rows)
    g = build(n, rows)
    assert [int(i) for i in g.input_positions()] == list(range(len(rows)))
    return g


def one10():
    return graph([(0, 1, 10, 1)])


def twoh(a, b):
    return graph([(0, 2, a, 1), (1, 2, b, 1)])


def par612():
    return graph([(0, 1, 6, 1), (0, 1, 6, 2)])


DIAM5 = [(0, 1, 6, 1), (0, 2, 6, 3), (1, 2, 4, 1), (1, 3, 6, 3), (2, 3, 6, 1)]


def diam5():
    return graph(DIAM5)


def x8():
    return graph([(0, 2, 0, 1), (1, 2, 8, 3)])


def policy(g, s, d, n, admission=None, ledger=None, eligibility=None, **routing):
    ledger = ledger if ledger is not None else ngc.FlowLedger(g)
    return ngc.FlowPolicy(ledger, s, d, t1_config(n, admission, **routing), eligibility)


def statuses(r):
    return [t.status for t in r.program.tunnels]


def routes(r):
    return [list(t.edges) if t.edges is not None else [] for t in r.program.tunnels]


def requests(r):
    return [t.requested_bandwidth for t in r.program.tunnels]


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


def err(fn):
    with pytest.raises(ValueError) as e:
        fn()
    return str(e.value)


# ---------------------------------------------------------------- the worked examples (X)

G1 = [(0, 1, 8192, 1), (0, 2, 40960, 1), (1, 3, 8192, 1), (2, 3, 40960, 1)]
BRANCH_ROUTES = [[0, 2]] * 2 + [[1, 3]] * 13 + [[]]
BRANCH_SHARES = [3277] * 12 + [3276] * 3 + [0]
EXAMPLES = {
    # id: (rows, heads, dst, D, N, routing keywords, requests, routes, statuses, reserved,
    #      reserved per edge, lossless bound, offer, Lossy shares per entry, delivery, source loss)
    "four tunnels on two parallel links": (
        [(0, 1, 6, 1), (0, 1, 6, 2)], [0], 1, 12, 4, {}, [3] * 4, [[0], [0], [1], [1]],
        [UP] * 4, 12, [6, 6], 12, 13, [4, 3, 3, 3], 12, 0),
    "uneven requests 3 and 2": (
        [(0, 1, 3, 1), (0, 1, 2, 2)], [0], 1, 5, 2, {}, [3, 2], [[0], [1]], [UP, UP], 5,
        [3, 2], 5, 5, [3, 2], 5, 0),
    "a share below the tunnel count": (
        [(0, 1, 6, 1), (0, 1, 6, 2)], [0], 1, 1, 2, {}, [1, 0], [[0], [0]],
        [UP, UP], 1, [1, 0], 6, 1, [1, 0], 1, 0),
    "carried above the holds": (
        [(0, 1, 5, 1), (0, 1, 5, 1)], [0], 1, 8, 2, {}, [4, 4], [[0], [1]],
        [UP, UP], 8, [4, 4], 10, 10, [5, 5], 10, 0),
    "a zero request on a zero-capacity link": (
        [(0, 1, 0, 1)], [0], 1, 0, 1, {}, [0], [[0]], [UP], 0, [0], 0, 1, [1],
        0, 0),
    "the last tunnel without bandwidth": (
        [(0, 1, 2, 1), (0, 1, 10, 2)], [0], 1, 12, 4, {}, [3] * 4, [[1], [1], [1], []],
        [UP, UP, UP, NB], 9, [0, 9], 10, 12, [4, 4, 4, 0], 10, 0),
    "two heads sharing a link": (
        [(0, 2, 6, 1), (1, 2, 6, 1), (2, 3, 6, 1)], [0, 1], 3, 12, 2, {}, [3] * 4,
        [[0, 2], [1, 2], [], []], [UP, UP, NB, NB], 6, [3, 3, 6], 6, 12, [6, 6, 0, 0], 6, 0),
    "three heads, one with a zero request": (
        [(0, 3, 20, 1), (1, 3, 20, 1), (2, 3, 20, 1)], [0, 1, 2], 3, 2, 1, {}, [1, 1, 0],
        [[0], [1], [2]], [UP] * 3, 2, [1, 1, 0], 60, 3, [1, 1, 1], 3, 0),
    "two equal links, one tunnel each": (
        [(0, 1, 5, 1), (0, 1, 5, 1)], [0], 1, 8, 2, {}, [4, 4], [[0], [1]], [UP, UP], 8,
        [4, 4], 10, 8, [4, 4], 8, 0),
    "a head behind a zero-capacity link": (
        [(0, 2, 0, 1), (1, 2, 8, 3)], [0, 1], 2, 4, 1, {}, [2, 2], [[], [1]], [NB, UP], 2,
        [0, 2], 0, 4, [0, 2], 2, 2),
    "a head behind a zero-capacity link, AbsoluteCost(3)": (
        [(0, 2, 0, 1), (1, 2, 8, 3)], [0, 1], 2, 4, 1,
        {"cost_ceiling": ngc.AbsoluteCost(3)}, [2, 2], [[], [1]], [NB, UP], 2,
        [0, 2], 0, None, None, None, None),
    "a head behind a zero-capacity link, CostStretch(2)": (
        [(0, 2, 0, 1), (1, 2, 8, 3)], [0, 1], 2, 4, 1,
        {"cost_ceiling": ngc.CostStretch(2)}, [2, 2], [[], []], [NB, CC], 0, [0, 0],
        0, None, None, None, None),
    # A-X8-Stretch-exact: 1.2 is exactly 6/5, so c0 = 5 admits the cost-6 path at equality;
    # A-X8-Stretch-cost7: a cost-7 path exceeds the ceiling.
    "CostStretch(1.2) admits the next path at equality": (
        [(0, 1, 0, 5), (0, 1, 1, 6)], [0], 1, 1, 1,
        {"cost_ceiling": ngc.CostStretch(1.2)}, [1], [[1]], [UP], 1, [0, 1], 1,
        1, [1], 1, 0),
    "CostStretch(6/5) below a path of cost 7": (
        [(0, 1, 0, 5), (0, 1, 1, 7)], [0], 1, 1, 1,
        {"cost_ceiling": ngc.CostStretch(Fraction(6, 5))}, [1], [[]], [CC], 0, [0, 0], 0,
        1, [0], 0, 1),
    "the diamond, design 8": (
        DIAM5, [0], 3, 8, 2, {}, [4, 4], [[0, 2, 4], []], [UP, NB], 4, [4, 0, 4, 0, 4], 4,
        None, None, None, None),
    "the diamond, design 10": (
        DIAM5, [0], 3, 10, 2, {}, [5, 5], [[0, 3], [1, 4]], [UP, UP], 10,
        [5, 5, 0, 5, 5], 12, None, None, None, None),
    "the diamond, design 14": (
        DIAM5, [0], 3, 14, 2, {}, [7, 7], [[], []], [NB, NB], 0, [0] * 5, 0, 14, [0, 0],
        0, 14),
    "two heads, the second with bandwidth": (
        [(0, 2, 2, 1), (1, 2, 10, 1)], [0, 1], 2, 12, 2, {}, [3] * 4,
        [[], [1], [], [1]], [NB, UP, NB, UP], 6, [0, 6], 0, 12, [0, 3, 0, 3], 6, 6),
    "more tunnels than units": (
        [(0, 1, 10, 1)], [0], 1, 3, 4, {}, [1, 1, 1, 0], [[0]] * 4, [UP] * 4, 3, [3], 10,
        3, [1, 1, 1, 0], 3, 0),
    "an idle tunnel on a zero-capacity link": (
        [(0, 1, 0, 1), (0, 1, 10, 2)], [0], 1, 1, 2, {}, [1, 0], [[1], [0]],
        [UP, UP], 1, [0, 1], 1, 2, [1, 1], 1, 0),
    "a zero design": (
        [(0, 1, 10, 1)], [0], 1, 0, 2, {}, [0, 0], [[0], [0]], [UP, UP], 0, [0],
        10, 4, [2, 2], 4, 0),
    # The sixteen tunnels' Lossy delivery is not in the example: by hand 6554 within 8192 plus
    # 40960 of 42598.
    "sixteen tunnels on two branches": (
        G1, [0], 3, 49152, 16, {}, [3072] * 16, BRANCH_ROUTES, [UP] * 15 + [NB], 46080,
        [6144, 39936, 6144, 39936], 47262, 49152, BRANCH_SHARES, 47514, 0),
    "sixteen tunnels, a design above the network": (
        G1, [0], 3, 409600, 16, {}, [25600] * 16, [[1, 3]] + [[]] * 15,
        [UP] + [NB] * 15, 25600, [0, 25600, 0, 25600], 40960, None, None, None,
        None),
}  # fmt: skip


@pytest.mark.parametrize("eid", list(EXAMPLES))
def test_contract_examples(eid):
    (rows_, heads, dst, D, N, kw, req, rts, st, reserved, per_edge, bound) = EXAMPLES[
        eid
    ][:12]
    offer, shares, delivery, loss = EXAMPLES[eid][12:]
    g = graph(rows_)
    fg = ngc.FlowLedger(g)
    p = policy(g, heads, dst, N, ledger=fg, **kw)
    r = p.replace(INF, design_demand=D, outputs=HOLDS)
    assert (requests(r), routes(r), statuses(r)) == (req, rts, st)
    assert r.plan.total_reserved_bandwidth == reserved
    assert r.edge_reservations.tolist() == per_edge == fg.reserved().tolist()
    assert (r.call_delivered, r.termination) == (bound, T.EXHAUSTED)
    if offer is None:
        return
    ly = policy(g, heads, dst, N, LOSSY, **kw)
    s = ly.replace(offer, design_demand=D)
    assert [t.assigned for t in s.by_tunnel] == shares
    assert (s.call_delivered, s.call_source_drops, s.termination) == (
        delivery,
        loss,
        T.COMPLETE,
    )
    assert (routes(s), statuses(s)) == (rts, st)  # the offer never changes the setup


def test_twoh_call():
    g = ngc.Graph(
        3,
        edge_sources=[0, 1],
        edge_destinations=[2, 2],
        edge_capacity_units=[2, 10],
        edge_costs=[1, 1],
        resolution=1.0,
    )
    cfg = ngc.FlowPolicyConfig(
        routing=ngc.DistributedTE(ngc.TunnelCount(2)), admission=LOSSY
    )
    fg = ngc.FlowLedger(g)
    p = ngc.FlowPolicy(fg, [0, 1], 2, cfg)
    r = p.replace(12, design_demand=12)
    assert [t.status for t in r.program.tunnels] == [NB, UP, NB, UP]
    assert (r.call_delivered, r.call_source_drops) == (6, 6)
    assert r.plan.total_reserved_bandwidth == 6
    assert (r.plan.tunnels_per_head, r.plan.design_demand, r.plan.design_by_head) == (
        2,
        12,
        ((0, 6), (1, 6)),
    )
    f = p.forward(8)
    assert (f.call_delivered, f.call_source_drops, f.delivery_delta) == (4, 4, -2)
    assert [t.assigned for t in f.by_tunnel] == [0, 2, 0, 2]
    assert fg.reserved().tolist() == [0, 6]


# ---------------------------------------------------------------- validation (V)


def test_matrix():
    # T1's row of the matrix, in order.
    g = one10()
    fg = ngc.FlowLedger(g)
    ca = ngc.DistributedTE(ngc.TunnelCount(2))
    # The fields of the other models are no keywords of DistributedTE or FlowPolicyConfig.
    for kw in (
        {"allocation": ngc.RouteAllocation.EQUAL},
        {"source_allocation": ngc.SourceAllocation.EQUAL},
        {"weight_basis": ngc.WeightBasis.CAPACITY},
        {"edge_weights": [1]},
        {"multipath": False},
    ):
        with pytest.raises(TypeError):
            ngc.DistributedTE(ngc.TunnelCount(2), **kw)
    for kw in (
        {"placement": ngc.RouteAllocation.EQUAL},
        {"origin": ngc.SourceAllocation.EQUAL},
    ):
        with pytest.raises(TypeError):
            ngc.FlowPolicyConfig(routing=ca, **kw)
    upto = ngc.FlowPolicyConfig(
        routing=ngc.DistributedTE(ngc.TunnelCount(2, up_to=True))
    )
    text = "DistributedTE requires an exact tunnel count"
    with pytest.raises(ValueError, match=exact(text)):
        ngc.FlowPolicy(fg, 0, 1, upto)
    with pytest.raises(ValueError, match=exact(text)):
        ngc.evaluate_demand(g, 0, 1, ngc.EvaluationOptions(upto, design_demand=1))
    with pytest.raises(
        ValueError, match=exact("DistributedTE does not use ShortestPathCache")
    ):
        ngc.FlowPolicy(
            fg, 0, 1, ngc.FlowPolicyConfig(routing=ca), cache=ngc.ShortestPathCache(g)
        )
    # Both admissions with every tie_break and parallel_edges (T1 always orders single paths:
    # no tie rule is ineffective); T1 has no multipath field.
    for tie in TIE.__members__.values():
        for parallel_edges in (True, False):
            for adm in (LOSSLESS, LOSSY):
                routing = ngc.DistributedTE(
                    ngc.TunnelCount(2), tie_break=tie, parallel_edges=parallel_edges
                )
                ngc.FlowPolicy(
                    fg, 0, 1, ngc.FlowPolicyConfig(routing=routing, admission=adm)
                )
    got = ngc.FlowPolicy(fg, 0, 1, ngc.FlowPolicyConfig(routing=ca)).config
    assert repr(got.routing) == repr(ca)
    assert (got.routing.tie_break, got.routing.parallel_edges) == (TIE.FEWER_HOPS, True)


def test_validation_order():
    # After the matrix: the count's range, then "up to", the ceiling, the search minimum.
    fg = ngc.FlowLedger(one10())
    for routing, text in (
        (
            ngc.DistributedTE(ngc.TunnelCount(0, up_to=True)),
            "tunnel count must be in [1, 10000]",
        ),
        (
            ngc.DistributedTE(
                ngc.TunnelCount(2, up_to=True), cost_ceiling=ngc.AbsoluteCost(-1)
            ),
            "DistributedTE requires an exact tunnel count",
        ),
        (
            ngc.DistributedTE(
                ngc.TunnelCount(2),
                cost_ceiling=ngc.AbsoluteCost(-1),
                min_unreserved_bandwidth=L62,
            ),
            "cost_ceiling AbsoluteCost cost must be >= 0 (got -1)",
        ),
        (
            ngc.DistributedTE(ngc.TunnelCount(2), min_unreserved_bandwidth=L62),
            "min_unreserved_bandwidth must be an integer in [0, 2^62)",
        ),
    ):
        with pytest.raises(ValueError, match=exact(text)):
            ngc.FlowPolicy(fg, 0, 1, ngc.FlowPolicyConfig(routing=routing))


def test_count_range():
    # [1, 10000] at the policy; a bool or a float at TunnelCount itself (a configuration value's
    # wrong type: TypeError).
    fg = ngc.FlowLedger(one10())
    for n in (0, 10001, -1):
        with pytest.raises(
            ValueError, match=exact("tunnel count must be in [1, 10000]")
        ):
            policy(None, 0, 1, n, ledger=fg)
    policy(None, 0, 1, 10000, ledger=fg)
    for bad, kind in ((True, "bool"), (1.0, "float")):
        with pytest.raises(
            TypeError,
            match=exact(f"TunnelCount per_head must be an integer (got {kind})"),
        ):
            ngc.TunnelCount(bad)


def test_parameters():
    # The search minimum's domain, the unreserved texts and the ceilings' kept texts.
    fg = ngc.FlowLedger(one10())
    for bad in (-1, L62):
        with pytest.raises(
            ValueError,
            match=exact("min_unreserved_bandwidth must be an integer in [0, 2^62)"),
        ):
            policy(None, 0, 1, 1, ledger=fg, min_unreserved_bandwidth=bad)
    policy(None, 0, 1, 1, ledger=fg, min_unreserved_bandwidth=L62 - 1)
    with pytest.raises(
        ValueError, match=exact("cost_ceiling AbsoluteCost cost must be >= 0 (got -1)")
    ):
        policy(None, 0, 1, 1, ledger=fg, cost_ceiling=ngc.AbsoluteCost(-1))
    for factor, shown in ((0.5, "1/2"), (Fraction(-(2**63)), "-9223372036854775808/1")):
        with pytest.raises(
            ValueError,
            match=exact(f"cost_ceiling CostStretch factor must be >= 1 (got {shown})"),
        ):
            policy(None, 0, 1, 1, ledger=fg, cost_ceiling=ngc.CostStretch(factor))


def test_tunnel_limit():
    # Heads after the mask times the count; FAN(107374) with N = 10000 just fits.
    def fan(h):
        return build(h + 1, [(i, h, 1, 1) for i in range(h)])

    big = fan(107375)
    with pytest.raises(
        ValueError,
        match=exact("tunnel limit exceeded: 107375 heads x 10000 tunnels > 1073741823"),
    ):
        policy(big, list(range(107375)), 107375, 10000)
    policy(
        big,
        list(range(107375)),
        107375,
        10000,
        eligibility=big.eligibility(excluded_nodes=[0]),
    )
    policy(fan(107374), list(range(107374)), 107374, 10000)


def test_operations_and_design_validation():
    # Place rejected (also 0), forward needs a program, design_demand only
    # for T1, in its domain, finite; the offer's checks first.
    p = policy(one10(), 0, 1, 1)
    for v in (5, 0):
        with pytest.raises(
            ValueError,
            match=exact("DistributedTE does not support incremental placement"),
        ):
            p.place(v)
    with pytest.raises(
        ValueError, match=exact("forward requires an installed program")
    ):
        p.forward(5)
    dom = "design_demand must be an integer in [0, 2^62) or permitted +inf"
    # 2**63 and -(2**64) overflow int64 in the binding: still the domain's error.
    for bad in (-1, L62, 2**63, -(2**64), True, 1.5):
        with pytest.raises(ValueError, match=exact(dom)):
            p.replace(5, design_demand=bad)
    with pytest.raises(
        ValueError,
        match=exact("offer must be an integer in [0, 2^62) or permitted +inf"),
    ):
        p.replace(-1, design_demand=-1)
    for call in (
        lambda: p.replace(5, design_demand=INF),
        lambda: p.replace(INF),
        lambda: p.evaluate(INF),
    ):
        with pytest.raises(
            ValueError, match=exact("DistributedTE requires a finite design demand")
        ):
            call()
    ly = policy(one10(), 0, 1, 1, LOSSY)
    with pytest.raises(
        ValueError, match=exact("LOSSY forwarding requires a finite offer")
    ):
        ly.replace(INF, design_demand=5)
    g = one10()
    pin = ngc.ForwardingProgram.from_routes(
        g, 0, 1, [ngc.RouteDefinition(0, ngc.Path(g, [[0]]))]
    )
    for cfg, model in (
        (COSTONLY, "ShortestPaths"),
        (GLOBAL, "CentralizedTE"),
        (ngc.FlowPolicyConfig(routing=ngc.Pinned(pin)), "Pinned"),
    ):
        q = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 1, cfg)
        for call in (q.replace, q.evaluate):
            with pytest.raises(
                ValueError, match=exact(f"{model} does not accept design_demand")
            ):
                call(5, design_demand=5)
    assert p.records() == []


# ---------------------------------------------------------------- setup (T)


def test_tie_breaks():
    # PARX, N = 2, D = 8.
    g = graph([(0, 1, 4, 1), (0, 1, 10, 1)])
    for s, want, held in (
        (TIE.FEWER_HOPS, [[0], [1]], [4, 4]),
        (TIE.HIGHER_CAPACITY, [[1], [1]], [0, 8]),
        (TIE.HIGHER_AVAILABLE_CAPACITY, [[1], [1]], [0, 8]),
        (TIE.HIGHER_HEADROOM, [[0], [1]], [4, 4]),
    ):
        fg = ngc.FlowLedger(g)
        r = ngc.FlowPolicy(fg, 0, 1, t1_config(2, tie_break=s)).replace(
            8, design_demand=8
        )
        assert (routes(r), fg.reserved().tolist()) == (want, held)


def test_tie_break_reads_unreserved():
    # PVIEW with unreserved [5, 8], N = 1, D = 4.
    g = graph([(0, 1, 10, 1), (0, 1, 8, 1)])
    for s, want in (
        (TIE.FEWER_HOPS, [0]),
        (TIE.HIGHER_CAPACITY, [0]),
        (TIE.HIGHER_AVAILABLE_CAPACITY, [1]),
        (TIE.HIGHER_HEADROOM, [1]),
    ):
        fg = ngc.FlowLedger(g, reservation_baseline=[5, 8])
        r = ngc.FlowPolicy(fg, 0, 1, t1_config(1, tie_break=s)).replace(
            4, design_demand=4
        )
        assert routes(r) == [want]


def test_tie_scores_read_unreserved():
    # PAR(10, 10) with the reservation baseline [3, 8], N = 1, D = 1: equal capacities, so
    # HIGHER_CAPACITY ties and falls to hops, then the lower EdgeId (e0); the availability is the
    # unreserved bandwidth: HIGHER_AVAILABLE_CAPACITY 3 < 8 and HIGHER_HEADROOM 3/10 < 8/10 take
    # e1. The service view (the full capacities) would tie both.
    g = graph([(0, 1, 10, 1), (0, 1, 10, 1)])
    for tie, want in (
        (TIE.FEWER_HOPS, [0]),
        (TIE.HIGHER_CAPACITY, [0]),
        (TIE.HIGHER_AVAILABLE_CAPACITY, [1]),
        (TIE.HIGHER_HEADROOM, [1]),
    ):
        fg = ngc.FlowLedger(g, reservation_baseline=[3, 8])
        r = ngc.FlowPolicy(fg, 0, 1, t1_config(1, tie_break=tie)).replace(
            1, design_demand=1
        )
        assert routes(r) == [want], tie


def test_filter_before_collapse_and_min_unreserved():
    # The filter before the collapse.
    g = graph([(0, 1, 2, 1), (0, 1, 10, 1)])
    for parallel_edges in (False, True):
        r = ngc.FlowPolicy(
            ngc.FlowLedger(g), 0, 1, t1_config(1, parallel_edges=parallel_edges)
        ).replace(6, design_demand=6)
        assert (routes(r), statuses(r)) == ([[1]], [UP])
    # PARC, N = 2, D = 2, minimum 5 and 0.
    parc = graph([(0, 1, 4, 1), (0, 1, 10, 2)])
    for minimum, path, held, bound in ((5, [1], [0, 2], 10), (0, [0], [2, 0], 4)):
        fg = ngc.FlowLedger(parc)
        r = policy(parc, 0, 1, 2, ledger=fg, min_unreserved_bandwidth=minimum).replace(
            INF, design_demand=2
        )
        assert (routes(r), fg.reserved().tolist(), r.call_delivered) == (
            [path, path],
            held,
            bound,
        )
    # ZL with minimum 1.
    zl = graph([(0, 1, 0, 1), (0, 1, 10, 2)])
    fg = ngc.FlowLedger(zl)
    r = policy(zl, 0, 1, 2, LOSSY, ledger=fg, min_unreserved_bandwidth=1).replace(
        2, design_demand=1
    )
    assert (requests(r), routes(r), fg.reserved().tolist(), r.call_delivered) == (
        [1, 0],
        [[1], [1]],
        [0, 1],
        2,
    )


def test_other_owners_holds():
    g = one10()
    fg = ngc.FlowLedger(g)
    a = policy(g, 0, 1, 1, ledger=fg)
    b = policy(g, 0, 1, 1, LOSSY, ledger=fg)
    assert statuses(a.replace(6, design_demand=6)) == [UP]
    assert (fg.reserved().tolist(), fg.unreserved().tolist()) == ([6], [4])
    down = b.replace(6, design_demand=6)
    assert (statuses(down), down.call_source_drops) == ([NB], 6)
    a.remove()
    assert statuses(b.replace(6, design_demand=6)) == [UP]
    assert fg.reserved().tolist() == [6]


def test_holds_do_not_queue():
    # Setup reads unreserved 10, forwarding the residual 0.
    g = one10()
    fg = ngc.FlowLedger(g)
    c = ngc.FlowPolicy(
        fg,
        0,
        1,
        ngc.FlowPolicyConfig(routing=ngc.ShortestPaths(), admission=LOSSY),
    )
    assert c.replace(10).call_delivered == 10
    r = policy(g, 0, 1, 1, LOSSY, ledger=fg).replace(4, design_demand=4)
    assert (
        statuses(r),
        fg.reserved().tolist(),
        r.call_edge_drops,
        r.call_delivered,
    ) == (
        [UP],
        [4],
        [(0, 4)],
        0,
    )


def test_replace_evaluate_and_zero():
    fg = ngc.FlowLedger(one10())
    p = policy(None, 0, 1, 1, ledger=fg)
    p.replace(6, design_demand=6)
    first = p.records()[0].flow_id
    r = p.replace(6, design_demand=6)
    assert (statuses(r), fg.reserved().tolist()) == ([UP], [6])
    assert not (r.routes_changed or r.weights_changed or r.reservations_changed)
    assert p.records()[0].flow_id != first
    state = ledger_now(fg)
    e = p.evaluate(3, design_demand=8, outputs=HOLDS)
    assert (
        statuses(e),
        requests(e),
        e.program.tunnels[0].reserved_bandwidth,
        e.plan.total_reserved_bandwidth,
    ) == ([UP], [8], 8, 8)
    assert (e.edge_reservations.tolist(), e.call_delivered) == ([8], 3)
    assert ledger_now(fg) == state
    fg2 = ngc.FlowLedger(one10())
    q = policy(None, 0, 1, 2, ledger=fg2)
    z = q.replace(0)
    assert (
        statuses(z),
        requests(z),
        fg2.reserved().tolist(),
        z.termination,
        z.call_delivered,
    ) == (
        [UP, UP],
        [0, 0],
        [0],
        T.COMPLETE,
        0,
    )
    z12 = q.replace(0, design_demand=12)
    assert (requests(z12), statuses(z12), fg2.reserved().tolist()) == (
        [6, 6],
        [UP, NB],
        [6],
    )


def test_infinite_offer_masked_source_and_empty_roster():
    r = policy(par612(), 0, 1, 4).replace(INF, design_demand=12)
    assert (
        r.call_delivered,
        r.call_undelivered,
        r.requested,
        r.unadmitted,
        r.termination,
    ) == (
        12,
        INF,
        INF,
        None,
        T.EXHAUSTED,
    )
    g = twoh(10, 10)
    r = policy(
        g, [0, 1], 2, 2, LOSSY, eligibility=g.eligibility(excluded_nodes=[0])
    ).replace(10, design_demand=10)
    assert (requests(r), r.plan.design_by_head) == ([5, 5], ((1, 10),))
    assert source_rows(r) == [(0, 0, 0, 0, 0, 0), (1, 10, 10, 10, 0, 0)]
    # EMPTY.
    empty = build(3, [(1, 2, 5, 1)])
    p = policy(empty, 0, 2, 2, LOSSY)
    r = p.replace(5, design_demand=5)
    assert (
        r.termination,
        r.program.sources,
        r.program.destinations,
        r.program.tunnels,
    ) == (
        T.NO_ROUTE,
        (0,),
        (2,),
        (),
    )
    assert (
        r.plan.tunnels_per_head,
        r.plan.design_demand,
        r.plan.design_by_head,
        r.plan.total_reserved_bandwidth,
    ) == (2, 5, (), 0)
    assert (r.call_source_drops, source_rows(r)) == (5, [(0, 0, 0, 0, 0, 0)])
    for z in (
        p.replace(0, design_demand=5),
        p.evaluate(0, design_demand=5),
        p.forward(0),
    ):
        assert (z.termination, z.program.tunnels) == (T.COMPLETE, ())
    assert p.records() == []


def test_all_down_and_round_order():
    # Every tunnel down (design 14 asks 7 per tunnel, no path has 7); the round order.
    for offer in (5, 14):
        r = policy(diam5(), 0, 3, 2, LOSSY).replace(offer, design_demand=14)
        assert (statuses(r), r.call_source_drops, r.call_delivered, r.termination) == (
            [NB, NB],
            offer,
            0,
            T.COMPLETE,
        )
        s = policy(diam5(), 0, 3, 2).replace(offer, design_demand=14)
        assert (s.call_delivered, s.termination) == (0, T.EXHAUSTED)
    r = policy(twoh(10, 10), [0, 1], 2, 2).replace(10, design_demand=10)
    assert requests(r) == [3, 3, 2, 2]
    assert [(t.entry_id, t.head, t.local_index) for t in r.program.tunnels] == [
        (0, 0, 0),
        (1, 1, 0),
        (2, 0, 1),
        (3, 1, 1),
    ]


def test_by_tunnel_and_stretch():
    r = policy(par612(), 0, 1, 4, LOSSY).replace(13, design_demand=12)
    assert rows(r) == [
        (0, 0, UP, 4, 4, 0),
        (1, 0, UP, 3, 2, 1),
        (2, 0, UP, 3, 3, 0),
        (3, 0, UP, 3, 3, 0),
    ]
    assert (r.call_edge_drops, r.call_delivered) == ([(0, 1)], 12)
    # The pool baseline; the largest numerator compares exactly.
    for f, second in ((2, CC), (2**63 - 1, UP)):
        r = policy(x8(), [0, 1], 2, 1, LOSSY, cost_ceiling=ngc.CostStretch(f)).replace(
            4, design_demand=4
        )
        assert statuses(r) == [NB, second]


def test_flags():
    p = policy(diam5(), 0, 3, 2)
    a = p.replace(8, design_demand=8)
    b = p.replace(10, design_demand=10)
    f = p.forward(10)
    assert (a.routes_changed, a.weights_changed, a.reservations_changed) == (
        True,
        True,
        True,
    )
    assert (b.routes_changed, b.weights_changed, b.reservations_changed) == (
        True,
        True,
        True,
    )
    assert (f.routes_changed, f.weights_changed, f.reservations_changed) == (False,) * 3
    d = policy(diam5(), 0, 3, 2, LOSSY).replace(5, design_demand=14)
    assert (d.routes_changed, d.weights_changed, d.reservations_changed) == (
        True,
        False,
        True,
    )
    # Equal requests, changed holds: two links of 4, N = 2, D = 8. Another owner's hold of 1 on
    # e0 leaves tunnel 1 down (reserved 4, 0); once it leaves, both come up (reserved 4, 4).
    g = graph([(0, 1, 4, 1), (0, 1, 4, 1)])
    fg = ngc.FlowLedger(g)
    other = policy(None, 0, 1, 1, ledger=fg)
    other.replace(1, design_demand=1)
    mine = policy(None, 0, 1, 2, ledger=fg)
    held = mine.replace(8, design_demand=8)
    other.remove()
    freed = mine.replace(8, design_demand=8)
    assert (statuses(held), statuses(freed)) == ([UP, NB], [UP, UP])
    assert [t.reserved_bandwidth for t in held.program.tunnels] == [4, 0]
    assert requests(held) == requests(freed) == [4, 4]
    assert freed.reservations_changed


def test_pinned_over_a_t1_program_reserves_nothing():
    # Pinned forwards an imported T1 program and books no holds: its first evaluate, replace
    # and place report no reservation change, the program keeps the imported requests and holds
    # as data, and the ledger and the records hold nothing. Also for an all-down program whose
    # requests are positive (a request of 20 on one link of 10).
    g = one10()
    up = policy(g, 0, 1, 1).replace(0, design_demand=6).program
    down = policy(g, 0, 1, 1).replace(0, design_demand=20).program
    assert [
        (t.status, t.requested_bandwidth, t.reserved_bandwidth) for t in up.tunnels
    ] == [(UP, 6, 6)]
    assert [
        (t.status, t.requested_bandwidth, t.reserved_bandwidth) for t in down.tunnels
    ] == [(NB, 20, 0)]
    for prog in (up, down):
        for op in ("evaluate", "replace", "place"):
            fg = ngc.FlowLedger(g)
            p = pinned(g, prog, 0, 1, LOSSLESS, ledger=fg)
            r = getattr(p, op)(1)
            assert not r.reservations_changed, (prog.tunnels[0].status, op)
            assert r.plan is None
            assert [
                (t.requested_bandwidth, t.reserved_bandwidth) for t in r.program.tunnels
            ] == [(t.requested_bandwidth, t.reserved_bandwidth) for t in prog.tunnels]
            assert fg.reserved().tolist() == [0]
            assert all(rec.edge_reservations == [] for rec in p.records())


# ---------------------------------------------------------------- records, rows, terminations (F)


def test_rows():
    # The TWOH call's policy's rows and records; forward recreates the records with their holds.
    p = policy(twoh(2, 10), [0, 1], 2, 2, LOSSY)
    r = p.replace(12, design_demand=12)
    assert rows(r) == [
        (0, 0, NB, 0, 0, 0),
        (1, 1, UP, 3, 3, 0),
        (2, 0, NB, 0, 0, 0),
        (3, 1, UP, 3, 3, 0),
    ]
    assert source_rows(r) == [(0, 6, 6, 0, 6, 0), (1, 6, 6, 6, 0, 0)]
    recs = p.records()
    assert [rec.entry_id for rec in recs] == [0, 1, 2, 3]
    assert [rec.min_route_cost for rec in recs] == [None, 1, None, 1]
    assert [len(rec.route.via_edges()) for rec in recs] == [0, 1, 0, 1]
    ids = [rec.flow_id for rec in recs]
    f = p.forward(8)
    assert not (f.routes_changed or f.weights_changed or f.reservations_changed)
    after = p.records()
    assert [rec.edge_reservations for rec in after] == [
        rec.edge_reservations for rec in recs
    ]
    assert all(a.flow_id != b for a, b in zip(after, ids, strict=True))


def test_terminations():
    # Zero offer, Lossy finite, Lossless full and partial, unlimited; a zero-weight route.
    p = policy(one10(), 0, 1, 1)
    assert p.replace(0, design_demand=5).termination == T.COMPLETE
    assert p.replace(10, design_demand=10).termination == T.COMPLETE
    assert p.replace(11, design_demand=10).termination == T.EXHAUSTED
    assert p.replace(INF, design_demand=10).termination == T.EXHAUSTED
    assert (
        policy(one10(), 0, 1, 1, LOSSY).replace(11, design_demand=10).termination
        == T.COMPLETE
    )
    g = one10()
    prog = ngc.ForwardingProgram.from_routes(
        g, 0, 1, [ngc.RouteDefinition(0, ngc.Path(g, [[0]]), weight=0)]
    )
    for adm in (LOSSLESS, LOSSY):
        z = ngc.FlowPolicy(
            ngc.FlowLedger(g),
            0,
            1,
            ngc.FlowPolicyConfig(routing=ngc.Pinned(prog), admission=adm),
        )
        assert z.replace(4).termination == T.NO_ACTIVE_WEIGHT


def test_no_certificates_and_plan_only_for_t1():
    g = twoh(2, 10)
    r = policy(g, [0, 1], 2, 2).replace(
        INF, design_demand=12, outputs=ngc.OutputOptions(min_cut=True)
    )
    assert (r.total_cost, r.min_cut, r.marginal_cost_curve) == (None, None, None)
    assert isinstance(r.plan, ngc.TunnelPlan)
    pin = ngc.FlowPolicy(
        ngc.FlowLedger(g),
        [0, 1],
        2,
        ngc.FlowPolicyConfig(routing=ngc.Pinned(r.program)),
    )
    assert pin.replace(12).plan is None
    for cfg in (COSTONLY, GLOBAL):
        q = ngc.FlowPolicy(ngc.FlowLedger(g), [0, 1], 2, cfg)
        assert q.place(5).plan is None and q.replace(5).plan is None
    with pytest.raises(TypeError):
        ngc.TunnelPlan()


# ---------------------------------------------------------------- programs (P)


def test_program_fields():
    r = policy(twoh(2, 10), [0, 1], 2, 2, LOSSY).replace(12, design_demand=12)
    got = [
        (
            t.entry_id,
            t.head,
            t.local_index,
            t.status,
            t.requested_bandwidth,
            t.reserved_bandwidth,
            t.weight,
            t.destination,
            t.edges,
            t.path_cost,
        )
        for t in r.program.tunnels
    ]
    assert got == [
        (0, 0, 0, NB, 3, 0, 0, None, None, None),
        (1, 1, 0, UP, 3, 3, 1, 2, (1,), 1),
        (2, 0, 1, NB, 3, 0, 0, None, None, None),
        (3, 1, 1, UP, 3, 3, 1, 2, (1,), 1),
    ]
    assert (r.program.sources, r.program.destinations, r.program.source_weights) == (
        (0, 1),
        (2,),
        None,
    )
    # Two-hop tunnels end at the destination: DIAM5, N = 2, D = 10 (requests 5, 5; e2 of 4 is
    # filtered out): tunnel 0 takes e0 e3, tunnel 1 e1 e4, both cost 4, both to node 3.
    two = policy(diam5(), 0, 3, 2).replace(10, design_demand=10)
    assert [(t.edges, t.destination, t.path_cost) for t in two.program.tunnels] == [
        ((0, 3), 3, 4),
        ((1, 4), 3, 4),
    ]
    g = twoh(2, 10)
    for t in (
        ngc.FlowPolicy(ngc.FlowLedger(g), [0, 1], 2, GLOBAL).place(12).program.tunnels
    ):
        assert (t.status, t.requested_bandwidth, t.reserved_bandwidth) == (UP, 0, 0)
        assert t.destination == 2 and t.edges is not None and t.path_cost is not None
    supplied = ngc.ForwardingProgram.from_routes(
        g, [0, 1], 2, [ngc.RouteDefinition(0, ngc.Path(g, [[0]]))]
    )
    assert [
        (t.status, t.requested_bandwidth, t.reserved_bandwidth, t.destination)
        for t in supplied.tunnels
    ] == [(UP, 0, 0, 2)]


def pinned(
    g, prog, s, d, admission=LOSSY, eligibility=None, black_holes=None, ledger=None
):
    cfg = ngc.FlowPolicyConfig(
        routing=ngc.Pinned(prog, black_holes=black_holes), admission=admission
    )
    ledger = ledger if ledger is not None else ngc.FlowLedger(g)
    return ngc.FlowPolicy(ledger, s, d, cfg, eligibility)


def test_t1_replay_stale_and_down_tunnels():
    g = par612()
    prog = policy(g, 0, 1, 4, LOSSY).replace(13, design_demand=12).program
    mask = g.eligibility(excluded_edges=[0])
    r = pinned(g, prog, 0, 1, eligibility=mask).replace(13)
    assert rows(r) == [
        (0, 0, FAILED, 0, 0, 0),
        (1, 0, FAILED, 0, 0, 0),
        (2, 0, UP, 7, 3, 4),
        (3, 0, UP, 6, 3, 3),
    ]
    assert (r.call_delivered, r.call_edge_drops) == (6, [(1, 7)])
    assert (
        pinned(g, prog, 0, 1, LOSSLESS, eligibility=mask).replace(13).call_delivered
        == 6
    )
    # Node 2 a black hole; the live forward the same, holds kept.
    g = diam5()
    fg = ngc.FlowLedger(g)
    p = policy(g, 0, 3, 2, LOSSY, ledger=fg)
    prog = p.replace(10, design_demand=10).program
    dark = g.eligibility(excluded_nodes=[2])
    r = pinned(g, prog, 0, 3, black_holes=dark).replace(10)
    assert rows(r) == [(0, 0, UP, 5, 5, 0), (1, 0, BLACKHOLED, 5, 0, 5)]
    assert (r.call_edge_drops, r.call_delivered) == ([(1, 5)], 5)
    assert (
        pinned(g, prog, 0, 3, LOSSLESS, black_holes=dark).replace(10).call_delivered
        == 1
    )
    held = fg.reserved().tolist()
    f = p.forward(10, black_holes=dark)
    assert (rows(f), f.call_delivered, fg.reserved().tolist()) == (rows(r), 5, held)
    assert [t.status for t in p.program.tunnels] == [UP, UP]
    g = twoh(2, 10)
    live = policy(g, [0, 1], 2, 2, LOSSY).replace(12, design_demand=12)
    pin = pinned(g, live.program, [0, 1], 2)
    r = pin.replace(12)
    assert (r.call_delivered, r.call_source_drops, rows(r), r.plan, r.termination) == (
        6,
        6,
        rows(live),
        None,
        T.COMPLETE,
    )
    # The replay roster comes from the topology basis: head 0 keeps its share and loses it.
    rg = graph([(0, 2, 10, 1), (0, 3, 10, 5), (1, 2, 10, 1), (3, 2, 10, 5)])
    prog = policy(rg, [0, 1], 2, 1, LOSSY).replace(10, design_demand=10).program
    x = pinned(
        rg, prog, [0, 1], 2, eligibility=rg.eligibility(excluded_edges=[0])
    ).replace(10)
    assert (x.call_delivered, x.call_source_drops, x.termination) == (5, 5, T.COMPLETE)
    assert source_rows(x) == [(0, 5, 5, 0, 5, 0), (1, 5, 5, 5, 0, 0)]


# ---------------------------------------------------------------- evaluators (Q)


def q(cfg, **kw):
    return ngc.EvaluationOptions(cfg, **kw)


def test_capacity():
    g = par612()
    assert (
        ngc.evaluate_demand(g, 0, 1, q(t1_config(4), design_demand=12)).call_delivered
        == 12
    )
    r = ngc.evaluate_demand(g, 0, 1, q(t1_config(4, LOSSY), design_demand=12, offer=13))
    assert (r.call_delivered, r.call_edge_drops, r.plan.total_reserved_bandwidth) == (
        12,
        [(0, 1)],
        12,
    )
    for opts, text in (
        (q(t1_config(4)), "DistributedTE evaluation requires design_demand"),
        (q(t1_config(4, LOSSY)), "DistributedTE evaluation requires design_demand"),
        (
            q(t1_config(4), design_demand=INF),
            "DistributedTE requires a finite design demand",
        ),
        (
            q(t1_config(4, LOSSY), design_demand=12),
            "DistributedTE LOSSY evaluation requires a finite offer",
        ),
        (
            q(t1_config(4, LOSSY), design_demand=12, offer=INF),
            "offer must be an integer in [0, 2^62)",
        ),
        (q(COSTONLY, design_demand=5), "ShortestPaths does not accept design_demand"),
        (q(GLOBAL, design_demand=5), "CentralizedTE does not accept design_demand"),
    ):
        with pytest.raises(ValueError, match=exact(text)):
            ngc.evaluate_demand(g, 0, 1, opts)
    # Lossless with a finite offer: the tunnels' delivery up to it.
    r5 = ngc.evaluate_demand(g, 0, 1, q(t1_config(4), design_demand=12, offer=5))
    assert r5.call_delivered == 5
    # The design is kept as given and checked by the evaluator.
    assert q(t1_config(4), design_demand=-1).design_demand == -1
    with pytest.raises(
        ValueError,
        match=exact("design_demand must be an integer in [0, 2^62) or permitted +inf"),
    ):
        ngc.evaluate_demand(g, 0, 1, q(t1_config(4), design_demand=-1))
    assert q(t1_config(4), design_demand=INF).design_demand == INF
    assert (
        q(t1_config(4), design_demand=7).design_demand == 7
        and q(COSTONLY).design_demand is None
    )
    prog = ngc.ForwardingProgram.from_routes(
        g, 0, 1, [ngc.RouteDefinition(0, ngc.Path(g, [[0]]))]
    )
    with pytest.raises(ValueError, match=exact("Pinned does not accept design_demand")):
        ngc.evaluate_demand(
            g, 0, 1, q(ngc.FlowPolicyConfig(routing=ngc.Pinned(prog)), design_demand=5)
        )
    rs = ngc.evaluate_pairs(one10(), [(0, 1), (0, 1)], q(t1_config(1), design_demand=6))
    assert [r.call_delivered for r in rs] == [10, 10]


def test_edge_failures():
    # X8, design 4: the bound is 0 (head 0's tunnel is down and its head keeps a share);
    # deleting the zero-capacity e0 drops head 0, so the bound rises by 8; deleting e1 loses 0.
    s = ngc.evaluate_edge_failures(x8(), [0, 1], 2, q(t1_config(1), design_demand=4))
    assert (
        s.baseline.call_delivered,
        [(e.edge, e.delivery_reduction) for e in s.entries],
    ) == (
        0,
        [(0, -8), (1, 0)],
    )
    g = graph([(0, 1, 6, 1), (0, 1, 6, 2), (2, 0, 5, 1)])
    s = ngc.evaluate_edge_failures(g, 0, 1, q(t1_config(4), design_demand=12))
    assert [(e.edge, e.delivery_reduction) for e in s.entries] == [
        (0, 6),
        (1, 6),
        (2, 0),
    ]
    zero = q(t1_config(1, LOSSY), design_demand=4, offer=0)
    r = ngc.evaluate_demand(x8(), [0, 1], 2, zero)
    assert (r.termination, r.call_delivered) == (T.COMPLETE, 0)
    assert all(
        e.delivery_reduction == 0
        for e in ngc.evaluate_edge_failures(x8(), [0, 1], 2, zero).entries
    )
    # Everything masked.
    gm = x8()
    masked = q(
        t1_config(1), design_demand=4, eligibility=gm.eligibility(excluded_edges=[0, 1])
    )
    r = ngc.evaluate_demand(gm, [0, 1], 2, masked)
    assert (r.termination, r.call_delivered) == (T.NO_ROUTE, 0)
    assert [x.termination for x in ngc.evaluate_pairs(gm, [([0, 1], 2)], masked)] == [
        T.NO_ROUTE
    ]


# ---------------------------------------------------------------- multipath is a bool (M)


def test_multipath_is_a_bool_and_the_pooled_rule():
    # None is no multipath value (it was the model's natural shape before); the pooled
    # single-path rule fires on multipath=False with a HIGHER_* tie_break only.
    with pytest.raises(TypeError):
        ngc.SpfOptions(ngc.Direction.FORWARD, multipath=None)
    with pytest.raises(TypeError):
        ngc.ShortestPaths(multipath=None)
    g = twoh(10, 10)

    def pooled(routing):
        cfg = ngc.FlowPolicyConfig(routing=routing)
        return ngc.FlowPolicy(ngc.FlowLedger(g), [0, 1], 2, cfg)

    pooled(ngc.ShortestPaths(tie_break=TIE.HIGHER_CAPACITY, parallel_edges=False))
    pooled(ngc.ShortestPaths(multipath=False))
    with pytest.raises(
        ValueError,
        match=exact(
            "ShortestPaths with multiple input sources requires tie_break=FEWER_HOPS"
            " when multipath=False"
        ),
    ):
        pooled(ngc.ShortestPaths(tie_break=TIE.HIGHER_CAPACITY, multipath=False))


# ---------------------------------------------------------------- numeric ranges (W)


def test_ranges():
    # ONEMAX.
    m = L62 - 1
    g = graph([(0, 1, m, 1)])
    fg = ngc.FlowLedger(g)
    r = policy(g, 0, 1, 1, LOSSY, ledger=fg).replace(m, design_demand=m)
    assert (
        statuses(r),
        r.plan.total_reserved_bandwidth,
        fg.reserved().tolist(),
        fg.unreserved().tolist(),
        r.call_delivered,
    ) == (
        [UP],
        m,
        [m],
        [0],
        m,
    )
    assert policy(g, 0, 1, 1).replace(INF, design_demand=m).call_delivered == m
    fg = ngc.FlowLedger(g)
    half = 1 << 61
    assert statuses(
        policy(g, 0, 1, 1, ledger=fg).replace(half, design_demand=half)
    ) == [UP]
    assert statuses(
        policy(g, 0, 1, 1, ledger=fg).replace(half, design_demand=half)
    ) == [NB]
    assert fg.unreserved().tolist() == [half - 1]
    # FAN4MAX with N = 10000: 40000 tunnels, all up, the design carried whole.
    cap = (1 << 60) - 1
    g = graph([(i, 4, cap, 1) for i in range(4)])
    fg = ngc.FlowLedger(g)
    D = L62 - 4
    r = policy(g, [0, 1, 2, 3], 4, 10000, ledger=fg).replace(D, design_demand=D)
    reqs = [t.requested_bandwidth for t in r.program.tunnels]
    assert len(reqs) == 40000 and all(t.status == UP for t in r.program.tunnels)
    assert (reqs.count(115292150460685), reqs.count(115292150460684)) == (
        4 * 6975,
        4 * 3025,
    )
    assert (
        sum(reqs),
        r.plan.total_reserved_bandwidth,
        fg.reserved().tolist(),
        r.call_delivered,
    ) == (
        D,
        D,
        [cap] * 4,
        D,
    )


# ---------------------------------------------------------------- zero offers and masks (E)


def test_zero_and_masked():
    for adm in (LOSSLESS, LOSSY):
        fg = ngc.FlowLedger(twoh(2, 10))
        p = policy(None, [0, 1], 2, 2, adm, ledger=fg)
        r = p.replace(0, design_demand=12)
        assert (
            r.termination,
            r.call_delivered,
            statuses(r),
            fg.reserved().tolist(),
        ) == (
            T.COMPLETE,
            0,
            [NB, UP, NB, UP],
            [0, 6],
        )
        assert [(rec.entry_id, rec.delivered) for rec in p.records()] == [
            (0, 0),
            (1, 0),
            (2, 0),
            (3, 0),
        ]
        state = ledger_now(fg)
        assert (
            p.evaluate(0, design_demand=12).termination == T.COMPLETE
            and ledger_now(fg) == state
        )
        f = p.forward(0)
        assert (
            f.termination,
            f.program is r.program or f.program.tunnels == r.program.tunnels,
        ) == (T.COMPLETE, True)
        assert not (f.routes_changed or f.weights_changed or f.reservations_changed)
        assert fg.reserved().tolist() == [0, 6]
        d = policy(diam5(), 0, 3, 2, adm)
        assert d.replace(0, design_demand=14).termination == T.COMPLETE
        assert d.replace(5, design_demand=14).termination == (
            T.COMPLETE if adm == LOSSY else T.EXHAUSTED
        )
        assert [rec.entry_id for rec in d.records()] == [0, 1]
        assert d.forward(0).termination == T.COMPLETE
        g = twoh(10, 10)
        e = policy(
            g, [0, 1], 2, 2, adm, eligibility=g.eligibility(excluded_edges=[0, 1])
        )
        assert e.replace(5, design_demand=5).termination == T.NO_ROUTE
        assert (
            e.replace(0, design_demand=5).termination == T.COMPLETE
            and e.records() == []
        )
        # Pinned over the TWOH call's program.
        g = twoh(2, 10)
        prog = policy(g, [0, 1], 2, 2, LOSSY).replace(12, design_demand=12).program
        x = pinned(
            g, prog, [0, 1], 2, adm, eligibility=g.eligibility(excluded_edges=[1])
        )
        r = x.replace(12)
        assert r.termination == T.NO_ACTIVE_WEIGHT
        assert [t.status for t in r.by_tunnel] == [NB, FAILED, NB, FAILED]
        if adm == LOSSY:
            assert source_rows(r) == [(0, 12, 12, 0, 12, 0), (1, 0, 0, 0, 0, 0)]
        assert x.replace(0).termination == T.COMPLETE
        assert [rec.entry_id for rec in x.records()] == [0, 1, 2, 3]
        y = pinned(
            g, prog, [0, 1], 2, adm, eligibility=g.eligibility(excluded_nodes=[0, 1])
        )
        assert (
            y.replace(12).termination,
            y.replace(0).termination,
            len(y.records()),
        ) == (T.NO_ROUTE, T.COMPLETE, 4)


# ---------------------------------------------------------------- the property against the oracle


def _random_case(rng):
    """A small random graph with a route from 0 to the last node: rows (u, v, capacity, cost)."""
    n = rng.randint(3, 6)
    raw = [
        (i, i + 1, rng.choice([0, 2, 4, 6, 8]), rng.randint(1, 4)) for i in range(n - 1)
    ]
    for _ in range(rng.randint(1, 7)):
        u, v = rng.randrange(n - 1), rng.randrange(1, n)
        if u != v:
            raw.append((u, v, rng.choice([0, 1, 2, 3, 4, 6, 8]), rng.randint(1, 4)))
    rng.shuffle(raw)
    return n, raw


def test_property():
    rng = random.Random(20261018)
    setups = lossless = lossy_cmp = all_up = 0
    for trial in range(400):
        n, rows_ = _random_case(rng)
        if not rows_:
            continue
        g = build(n, rows_)
        # The oracle's edges in EdgeId order (the graph's own order of the input rows).
        edges = [rows_[int(i)] for i in g.input_positions()]
        sources = sorted(rng.sample(range(n - 1), rng.randint(1, min(2, n - 1))))
        dests = [n - 1]
        mask = None
        eligibility = None
        if rng.random() < 0.3:
            xe = [rng.randrange(len(rows_))]
            mask = ([], xe)
            eligibility = g.eligibility(excluded_edges=xe)
        caps = [c for _, _, c, _ in edges]
        unres = [rng.randint(0, c) for c in caps] if rng.random() < 0.4 else list(caps)
        resid = [rng.randint(0, c) for c in caps] if rng.random() < 0.4 else list(caps)
        fg = ngc.FlowLedger(g, resid, reservation_baseline=unres)
        held = [0] * len(rows_)
        if rng.random() < 0.3 and sources:
            other = ngc.FlowPolicy(fg, sources[:1], dests, t1_config(1))
            other.replace(0, design_demand=rng.randint(0, 4))
            held = fg.reserved().tolist()
        php = rng.random() < 0.5
        prio = rng.choice(["CAPACITY", "AVAILABLE_CAPACITY", "HEADROOM"])
        me = rng.random() < 0.7
        tie = f"HIGHER_{prio}" if php else "FEWER_HOPS"
        s = {"tie_break": getattr(TIE, tie), "parallel_edges": me}
        minimum = rng.choice([0, 0, 1, 2])
        ceiling = rng.choice(
            [
                None,
                None,
                ("abs", rng.randint(1, 6)),
                ("str", rng.choice([1.0, 1.2, 1.5, 2.0])),
            ]
        )
        kw = {"min_unreserved_bandwidth": minimum}
        okw = {"minimum": minimum}
        if ceiling and ceiling[0] == "abs":
            kw["cost_ceiling"], okw["absolute"] = (
                ngc.AbsoluteCost(ceiling[1]),
                ceiling[1],
            )
        elif ceiling:
            kw["cost_ceiling"], okw["stretch"] = ngc.CostStretch(ceiling[1]), ceiling[1]
        N = rng.randint(1, 3)
        D = rng.randint(0, 12)
        want = orc.setup(
            edges,
            sources,
            dests,
            D,
            N,
            sel=orc.Selection(tie, me),
            mask=mask,
            unreserved=fg.unreserved().tolist(),
            held=None,
            service=fg.residual().tolist(),
            **okw,
        )
        before_holds = fg.reserved().tolist()
        V = fg.residual().tolist()
        p = ngc.FlowPolicy(fg, sources, dests, t1_config(N, **s, **kw), eligibility)
        r = p.replace(INF, design_demand=D)
        where = (trial, rows_, sources, D, N, s, kw)
        got = [
            (t.status.name, list(t.edges or []), t.requested_bandwidth)
            for t in r.program.tunnels
        ]
        assert got == [(t.status, t.path or [], t.requested) for t in want.tunnels], (
            where
        )
        setups += 1
        # Every up path passed its filter; holds = the up requests; within R0.
        assert [
            a - b for a, b in zip(fg.reserved().tolist(), before_holds, strict=True)
        ] == want.holds, where
        assert all(
            x <= u for x, u in zip(fg.reserved().tolist(), unres, strict=True)
        ), where
        # The program's lossless bound on the service view.
        if want.roster:
            assert r.call_delivered == orc.lossless_bound(
                edges, want.roster, want.tunnels, V
            ), where
            lossless += 1
        # Every tunnel up and U <= V: the design is carried losslessly.
        U = [u - h for u, h in zip(unres, held, strict=True)]
        if (
            want.roster
            and all(t.up for t in want.tunnels)
            and all(a <= b for a, b in zip(U, V, strict=True))
        ):
            assert p.replace(D, design_demand=D).call_delivered == D, where
            all_up += 1
        # Lossy: the oracle's loss rule, conservation and the setup independent of the offer.
        offer = rng.randint(0, 14)
        # The same setup view: the other owner's holds taken out of the reservation baseline.
        seen = [u - h for u, h in zip(unres, held, strict=True)]
        lp = ngc.FlowPolicy(
            ngc.FlowLedger(g, resid, reservation_baseline=seen),
            sources,
            dests,
            t1_config(N, LOSSY, **s, **kw),
            eligibility,
        )
        lr = lp.replace(offer, design_demand=D)
        assert [(t.status.name, list(t.edges or [])) for t in lr.program.tunnels] == [
            (a, b) for a, b, _ in got
        ], where
        assert (
            lr.offered
            == lr.source_drops + sum(x for _, x in lr.edge_drops) + lr.delivered
        ), where
        assert (
            lr.offered == sum(row.offered for row in lr.by_source) or not want.roster
        ), where
        if want.roster:
            try:
                delivered, lost, per, _ = orc.lossy(
                    edges, want.roster, want.tunnels, offer, V
                )
            except AssertionError:
                continue
            assert (lr.call_delivered, lr.call_source_drops) == (delivered, lost), where
            assert [
                (t.assigned, t.delivered, t.total_edge_drops) for t in lr.by_tunnel
            ] == per, where
            lossy_cmp += 1
    # Not vacuous.
    assert setups > 300 and lossless > 250 and lossy_cmp > 200 and all_up > 20, (
        setups,
        lossless,
        lossy_cmp,
        all_up,
    )


# ---------------------------------------------------------------- Python allocation failures


@needs_hooks
@pytest.mark.parametrize("method", ["replace", "evaluate", "forward"])
def test_t1_results_survive_python_allocation_failures(method):
    # Every Python object of a T1 result (the plan, tunnels with None fields, edge_reservations)
    # exists before anything is published. The calls are positional (the design is the offer):
    # pybind11's parsing of a call that passes keywords does not check its allocations, so a
    # keyword call may crash there (limit 1 of the module docstring).
    g = twoh(2, 10)

    def make():
        fg = ngc.FlowLedger(g)
        p = ngc.FlowPolicy(fg, [0, 1], 2, t1_config(2, LOSSY))
        p.replace(4, design_demand=12)

        def snapshot():
            return [full_ledger_state(fg), repr(p._state())]

        return getattr(p, method), (8,), snapshot

    assert sweep_python_allocation(make) >= 3


GETTER_SWEEP = """
import _netgraph_core as ext, netgraph_core as ngc
g = ngc.Graph(3, [0, 1], [2, 2], [2, 10], [1, 1], 1.0)
fg = ngc.FlowLedger(g)
cfg = ngc.FlowPolicyConfig(routing=ngc.DistributedTE(ngc.TunnelCount(2)), admission=ngc.Admission.LOSSY)
p = ngc.FlowPolicy(fg, [0, 1], 2, cfg)
r = p.replace(12, design_demand=12)
fid = p.records()[0].flow_id
ca = cfg.routing
t0 = r.program.tunnels[0]
abs3 = ngc.DistributedTE(ngc.TunnelCount(2), cost_ceiling=ngc.AbsoluteCost(3))
st2 = ngc.DistributedTE(ngc.TunnelCount(2), cost_ceiling=ngc.CostStretch(2.0))
sp = ngc.ShortestPaths(tie_break=ngc.RouteTieBreak.HIGHER_HEADROOM, multipath=False)
fn = lambda: {getter}
expected = repr(fn())  # the getter reads its value first: a missing name fails here
for k in range(400):
    hit, seen, value, error = ext._fail_python_allocation(fn, (), k)
    if not hit:
        assert error is None and repr(value) == expected, (k, error, value)
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
        "r.program.tunnels",
        "fg.reserved()",
        "fg.unreserved()",
        "r.plan",
        "r.plan.design_by_head",
        "ca.tunnel_count",
        "t0.status",
        "repr(abs3)",
        "repr(st2)",
        "st2.cost_ceiling.factor",
        "st2.cost_ceiling",
        "cfg.routing",
        "repr(r)",
        "repr(sp)",
        "p.records()",
        "p.record(fid)",
        "fg.records()",
        "fg.record(fid)",
    ],
)
def test_t1_getters_survive_python_allocation_failures(getter):
    # Each getter's sweep runs in its own process: a crash fails the test instead of the run.
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


# ---------------------------------------------------------------- the exact stretch (S)


def test_stretch_exact():
    # A-X8-Stretch-exact / -cost7: links of costs 5 (capacity 0) and 6 or 7, D = 1, N = 1. c0 = 5
    # counts the zero-capacity link; 1.2 is exactly 6/5, so cost 6 is admitted at equality (up on
    # e1, 1 delivered, Lossless bound 1) and cost 7 is not (the share lost at the source).
    for cost, status, route, delivered in ((6, UP, [1], 1), (7, CC, [], 0)):
        g = graph([(0, 1, 0, 5), (0, 1, 1, cost)])
        factor = 1.2 if cost == 6 else Fraction(6, 5)
        lossy = policy(g, 0, 1, 1, LOSSY, cost_ceiling=ngc.CostStretch(factor))
        r = lossy.replace(1, design_demand=1)
        assert (statuses(r), routes(r)) == ([status], [route]), cost
        assert (r.call_delivered, r.call_source_drops) == (delivered, 1 - delivered), (
            cost
        )
        lossless = policy(g, 0, 1, 1, cost_ceiling=ngc.CostStretch(factor))
        assert lossless.replace(INF, design_demand=1).call_delivered == delivered, cost


def test_stretch_below_one():
    # T1-stretch-below-one: 1/2 is refused where the policy checks the ceiling (after the tunnel
    # count); the lowest numerator too.
    g = graph([(0, 1, 10, 5)])
    for factor, shown in ((Fraction(1, 2), "1/2"), (-(2**63), f"{-(2**63)}/1")):
        text = f"cost_ceiling CostStretch factor must be >= 1 (got {shown})"
        with pytest.raises(ValueError, match=exact(text)):
            policy(g, 0, 1, 1, cost_ceiling=ngc.CostStretch(factor))
    with pytest.raises(ValueError, match=exact("tunnel count must be in [1, 10000]")):
        policy(g, 0, 1, 0, cost_ceiling=ngc.CostStretch(0.5))


def test_stretch_wide_products():
    # BIG(2^60): links 0 -> 1 of capacity 1 and costs 2^60 and 2^60 + 1, D = 2, N = 2: tunnel 0
    # (request 1) is up on e0; tunnel 1 finds e0 full and e1 at cost 2^60 + 1: up with
    # (2^60 + 1) / 2^60 (equality), COST_CEILING_EXCEEDED with (2^61 + 1) / 2^61.
    c0 = 2**60
    g = graph([(0, 1, 1, c0), (0, 1, 1, c0 + 1)])
    for factor, second in (
        (Fraction(2**60 + 1, 2**60), UP),
        (Fraction(2**61 + 1, 2**61), CC),
    ):
        r = policy(g, 0, 1, 2, cost_ceiling=ngc.CostStretch(factor)).replace(
            INF, design_demand=2
        )
        assert statuses(r) == [UP, second], factor
