"""Edges of capacity 0 stay in routing, through the binding: the worked examples of zero
capacity (SPF and KSP over the zero edge and their residual forms, the search graph S under every
model, a source behind a zero link, a zero first or last link of a Pinned split route, a link of
service 0, the zero-weight WCMP cases, T1's roster and its CostStretch pool baseline) and the minimum
cut. The graphs list their rows in EdgeId order (asserted), so an
EdgeId is the row's position. S: 0: 0 -> 1 cap 0, 1: 0 -> 2 cap 5, 2: 2 -> 1 cap 5 (costs 1)."""

from __future__ import annotations

import math
from fractions import Fraction

import pytest
from helpers import (
    GLOBAL,
    HIGHER_CAP,
    HIGHER_HEADROOM,
    SHORTEST_PATHS,
    config,
    dag_of,
    t1_config,
)

import netgraph_core as ngc

EB = ngc.RouteAllocation.EQUAL
W = ngc.RouteAllocation.WEIGHTED
MF = ngc.RouteAllocation.OPTIMIZED
LOSSY = ngc.Admission.LOSSY
FWD = ngc.Direction.FORWARD
S = [(0, 1, 0, 1), (0, 2, 5, 1), (2, 1, 5, 1)]


def graph(n, rows):
    g = ngc.Graph(n, *map(list, zip(*rows, strict=True)), 1.0)
    assert [int(i) for i in g.input_positions()] == list(range(len(rows)))
    return g


def edges(dag):
    return sorted(int(e) for e in dag.via_edges())


def walks(result):
    return [(c, [hop[0] for hop in p.hops()]) for c, p in result]


def place(g, s, d, cfg, offer, eligibility=None):
    policy = ngc.FlowPolicy(ngc.FlowLedger(g), s, d, cfg, eligibility)
    return policy.place(offer), policy


def test_spf_and_ksp_search_the_zero_edge():
    g = graph(3, S)
    r = ngc.spf(g, 0, ngc.SpfOptions(FWD, targets=1))
    assert (r.target_cost, edges(r.dag)) == (1, [0])  # SPF over the zero edge
    gated = ngc.spf(g, 0, ngc.SpfOptions(FWD, targets=1, residual=[0, 5, 5]))
    assert (gated.target_cost, edges(gated.dag)) == (2, [1, 2])
    full = ngc.spf(g, 0, ngc.SpfOptions(FWD))
    assert (full.costs.tolist(), edges(full.dag)) == ([0, 1, 1], [0, 1])
    assert walks(ngc.ksp(g, 0, 1, ngc.KspOptions(2))) == [(1, [0]), (2, [1, 2])]
    opts = ngc.KspOptions(2, residual=[0, 5, 5])
    assert walks(ngc.ksp(g, 0, 1, opts)) == [(2, [1, 2])]
    # CostStretch(1.2), exactly 6/5: the static baseline 5 counts the zero-capacity edge, and
    # 6 * 5 <= 6 * 5 admits the cost-6 walk at equality; 11/10 does not.
    h = graph(2, [(0, 1, 0, 5), (0, 1, 1, 6)])
    ceiled = ngc.KspOptions(2, cost_ceiling=ngc.CostStretch(1.2))
    assert walks(ngc.ksp(h, 0, 1, ceiled)) == [(5, [0]), (6, [1])]
    below = ngc.KspOptions(2, cost_ceiling=ngc.CostStretch(Fraction(11, 10)))
    assert walks(ngc.ksp(h, 0, 1, below)) == [(5, [0])]


def test_headroom_without_a_residual_reads_one_or_zero():
    d = graph(4, [(0, 1, 1, 1), (0, 2, 5, 1), (1, 3, 1, 1), (2, 3, 5, 1)])
    z = graph(4, [(0, 1, 0, 1), (0, 2, 5, 1), (1, 3, 5, 1), (2, 3, 5, 1)])

    def path(g, tie):
        o = ngc.SpfOptions(FWD, tie_break=tie, multipath=False, targets=3)
        return edges(ngc.spf(g, 0, o).dag)

    assert path(d, HIGHER_HEADROOM) == [
        0,
        2,
    ]  # every positive edge reads 1: the EdgeId order decides
    assert path(d, HIGHER_CAP) == [1, 3]
    assert path(z, HIGHER_HEADROOM) == [1, 3]  # 0 on the zero-capacity edge


def test_g8_for_every_model():
    g = graph(3, S)
    r, p = place(g, 0, 1, SHORTEST_PATHS, 5)
    assert (r.call_delivered, r.call_undelivered, r.termination) == (
        0,
        5,
        ngc.Termination.EXHAUSTED,
    )
    [rec] = p.records()
    assert (rec.delivered, edges(rec.route), rec.min_route_cost, rec.edge_flows) == (
        0,
        [0],
        1,
        [],
    )
    r, p = place(g, 0, 1, config(ngc.ShortestPaths(), EB, LOSSY), 5)
    assert (r.call_delivered, r.edge_drops, r.termination) == (
        0,
        [(0, 5)],
        ngc.Termination.COMPLETE,
    )
    for admission in (ngc.Admission.LOSSLESS, LOSSY):
        r, p = place(g, 0, 1, config(ngc.ShortestPaths(), W, admission), 5)
        assert (r.call_delivered, r.edge_drops) == (
            0,
            [],
        )  # weight 0: no member, never a bypass
    assert place(g, 0, 1, config(ngc.ShortestPaths(), MF), 5)[0].call_delivered == 0
    assert place(g, 0, 1, GLOBAL, 5)[0].call_delivered == 5
    # T1: the request of 5 filters the zero link out; the tunnel comes up over 0 -> 2 -> 1.
    t1 = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 1, t1_config(1)).replace(
        5, design_demand=5
    )
    assert (t1.call_delivered, t1.program.tunnels[0].edges) == (5, (1, 2))


def test_equal_and_weighted_splits_over_a_zero_member():
    first = graph(2, [(0, 1, 0, 1), (0, 1, 10, 1)])
    last = graph(2, [(0, 1, 10, 1), (0, 1, 0, 1)])
    assert place(first, 0, 1, SHORTEST_PATHS, 10)[0].call_delivered == 0
    assert place(last, 0, 1, SHORTEST_PATHS, 10)[0].call_delivered == 1
    r, _ = place(first, 0, 1, config(ngc.ShortestPaths(), EB, LOSSY), 10)
    assert (r.call_delivered, r.edge_drops) == (5, [(0, 5)])
    assert (
        place(first, 0, 1, config(ngc.ShortestPaths(), W), 10)[0].call_delivered == 10
    )
    # A WCMP member of capacity 0, and all weights zero.
    wcmp = graph(2, [(0, 1, 0, 1), (0, 1, 40, 1)])
    assert (
        place(wcmp, 0, 1, config(ngc.ShortestPaths(), W, LOSSY), 40)[0].call_delivered
        == 40
    )
    one = graph(2, [(0, 1, 0, 1)])
    r, p = place(one, 0, 1, config(ngc.ShortestPaths(), W, LOSSY), 2)
    [rec] = p.records()
    assert (r.call_delivered, r.edge_drops, rec.offered, rec.delivered) == (0, [], 2, 0)
    # A link of service 0: both units dropped on the zero-capacity link.
    r, _ = place(one, 0, 1, config(ngc.ShortestPaths(), EB, LOSSY), 2)
    assert (r.call_delivered, r.edge_drops) == (0, [(0, 2)])


@pytest.mark.parametrize(
    ("caps", "lossless", "drops"),
    [((0, 10), 0, [(0, 5)]), ((10, 0), 1, [(1, 5)])],
    ids=["zero-first", "positive-first"],
)
def test_two_sources_one_behind_a_zero_link(caps, lossless, drops):
    # The pooled equal source allocation: a LOSSLESS total a is admitted only when RR(a) fits, so a source
    # behind the zero link caps it at its round-robin position: 0 when source 0 is behind it
    # (RR(1) gives source 0 the unit), 1 when source 1 is (RR(2) gives source 1 a unit). LOSSY
    # offers RR(10) = 5 and 5 and drops the share behind the zero link.
    g = graph(3, [(0, 2, caps[0], 1), (1, 2, caps[1], 1)])
    r, p = place(g, [0, 1], 2, SHORTEST_PATHS, 10)
    assert (r.call_delivered, r.termination) == (lossless, ngc.Termination.EXHAUSTED)
    assert [row.source for row in r.by_source] == [0, 1]
    assert len(p.records()) == 2
    r, _ = place(g, [0, 1], 2, config(ngc.ShortestPaths(), EB, LOSSY), 10)
    assert (r.call_delivered, r.edge_drops) == (5, drops)


def test_b3_delete_zero_and_edge_failures():
    g = graph(3, [(0, 2, 0, 1), (1, 2, 10, 1)])
    for eligibility in (
        g.eligibility(excluded_nodes=[0]),
        g.eligibility(excluded_edges=[0]),
    ):
        assert (
            place(g, [0, 1], 2, SHORTEST_PATHS, 10, eligibility)[0].call_delivered == 10
        )
        assert place(g, [0, 1], 2, GLOBAL, 10, eligibility)[0].call_delivered == 10
    # The zero-capacity source first, under the pooled equal source allocation, admits 0; removing the zero link takes source 0
    # out of the program, so the deletion gains 10 (a negative loss); removing 1 -> 2 leaves 0.
    s = ngc.evaluate_edge_failures(g, [0, 1], 2, ngc.EvaluationOptions(SHORTEST_PATHS))
    assert (
        s.baseline.call_delivered,
        [(e.edge, e.delivery_reduction) for e in s.entries],
    ) == (
        0,
        [(0, -10), (1, 0)],
    )
    # A negative loss for the zero-capacity link.
    s8 = ngc.evaluate_edge_failures(
        graph(3, S), 0, 1, ngc.EvaluationOptions(SHORTEST_PATHS)
    )
    assert [(e.edge, e.delivery_reduction) for e in s8.entries] == [
        (0, -5),
        (1, 0),
        (2, 0),
    ]


def test_t1_roster_and_stretch_baseline():
    # A-X8: head 0's link has no unreserved capacity: its tunnel is down and it loses its 2.
    x8 = graph(3, [(0, 2, 0, 1), (1, 2, 8, 3)])
    p = ngc.FlowPolicy(ngc.FlowLedger(x8), [0, 1], 2, t1_config(1, LOSSY))
    r = p.replace(4, design_demand=4)
    assert (r.call_delivered, r.call_source_drops, r.termination) == (
        2,
        2,
        ngc.Termination.COMPLETE,
    )
    assert [t.status for t in r.by_tunnel] == [
        ngc.TunnelStatus.INSUFFICIENT_BANDWIDTH,
        ngc.TunnelStatus.UP,
    ]
    # CostStretch(1.2), exactly 6/5: c0 = 5 (the zero-capacity link counts), and the cost-6 path is
    # admitted at equality; at 11/10 it exceeds the ceiling.
    h = graph(2, [(0, 1, 0, 5), (0, 1, 1, 6)])
    for factor, delivered, status in (
        (1.2, 1, ngc.TunnelStatus.UP),
        (Fraction(11, 10), 0, ngc.TunnelStatus.COST_CEILING_EXCEEDED),
    ):
        cfg = t1_config(1, LOSSY, cost_ceiling=ngc.CostStretch(factor))
        p = ngc.FlowPolicy(ngc.FlowLedger(h), 0, 1, cfg)
        r = p.replace(1, design_demand=1)
        assert (r.call_delivered, r.by_tunnel[0].status) == (delivered, status), factor
    assert [rec.entry_id for rec in p.records()] == [0]


@pytest.mark.parametrize(
    ("caps", "bound"), [((0, 4), 0), ((4, 0), 1)], ids=["zero-first", "zero-last"]
)
def test_pinned_routes_keep_their_zero_links(caps, bound):
    # Two split routes, one over the zero link: the head's round robin gives the first unit to
    # route 0, so a zero first route blocks every total and a zero last route admits 1.
    g = graph(2, [(0, 1, caps[0], 1), (0, 1, caps[1], 1)])
    program = ngc.ForwardingProgram.from_routes(
        g,
        0,
        1,
        [
            ngc.RouteDefinition(0, dag_of(g, [0])),
            ngc.RouteDefinition(0, dag_of(g, [1])),
        ],
    )
    r, p = place(g, 0, 1, ngc.FlowPolicyConfig(routing=ngc.Pinned(program)), math.inf)
    assert (r.call_delivered, len(p.records())) == (bound, 2)


def test_global_min_cut_lists_zero_capacity_links_leaving_the_closure():
    want = ngc.OutputOptions(min_cut=True)
    b = graph(3, [(0, 2, 0, 1), (1, 2, 10, 1)])
    r = ngc.evaluate_demand(b, [0, 1], 2, ngc.EvaluationOptions(GLOBAL, outputs=want))
    assert (r.call_delivered, r.min_cut) == (10, [0, 1])
    s = graph(3, S)
    r = ngc.evaluate_demand(s, 0, 1, ngc.EvaluationOptions(GLOBAL, outputs=want))
    assert (r.call_delivered, r.min_cut) == (5, [0, 1])
    masked = ngc.EvaluationOptions(
        GLOBAL, eligibility=s.eligibility(excluded_edges=[0]), outputs=want
    )
    assert ngc.evaluate_demand(s, 0, 1, masked).min_cut == [1]
