"""Edges of capacity 0 stay in routing, through the binding: the contract's step-1 examples
(SPF-zero-edge, KSP-zero-edge and their residual forms, the owner page's G8, B3, E-zero-first and
E-zero-last, D-zero-service, the zero-weight WCMP cases, A-X8's delivery and A-X8-Stretch-exact's
baseline) and the minimum cut. The graphs list their rows in EdgeId order (asserted), so an
EdgeId is the row's position. S: 0: 0 -> 1 cap 0, 1: 0 -> 2 cap 5, 2: 2 -> 1 cap 5 (costs 1)."""

from __future__ import annotations

import math

import pytest
from helpers import COSTONLY, GLOBAL, PREFER, config, dag_of, sel

import netgraph_core as ngc

EB = ngc.FlowPlacement.EQUAL_BALANCED
W = ngc.FlowPlacement.WEIGHTED
MF = ngc.FlowPlacement.MAX_FLOW
LOSSY = ngc.Admission.LOSSY
HEAD = ngc.EdgePriority.HEADROOM
FWD = ngc.Direction.FORWARD
S = [(0, 1, 0, 1), (0, 2, 5, 1), (2, 1, 5, 1)]


def graph(n, rows):
    g = ngc.Graph(n, *map(list, zip(*rows, strict=True)), 1.0)
    assert [int(i) for i in g.input_index()] == list(range(len(rows)))
    return g


def edges(dag):
    return sorted(int(e) for e in dag.via_edges())


def walks(result):
    return [(c, [hop[0] for hop in p.hops()]) for c, p in result]


def place(g, s, d, cfg, volume, elig=None):
    policy = ngc.FlowPolicy(ngc.FlowLedger(g), s, d, cfg, elig)
    return policy.place(volume), policy


def test_spf_and_ksp_search_the_zero_edge():
    g = graph(3, S)
    r = ngc.spf(g, 0, ngc.SpfOptions(FWD, sel(), targets=1))
    assert (r.target_cost, edges(r.dag)) == (1, [0])  # SPF-zero-edge
    gated = ngc.spf(g, 0, ngc.SpfOptions(FWD, sel(), targets=1, residual=[0, 5, 5]))
    assert (gated.target_cost, edges(gated.dag)) == (2, [1, 2])
    full = ngc.spf(g, 0, ngc.SpfOptions(FWD, sel()))
    assert (full.costs.tolist(), edges(full.dag)) == ([0, 1, 1], [0, 1])
    assert walks(ngc.ksp(g, 0, 1, ngc.KspOptions(2, sel()))) == [(1, [0]), (2, [1, 2])]
    opts = ngc.KspOptions(2, sel(), residual=[0, 5, 5])
    assert walks(ngc.ksp(g, 0, 1, opts)) == [(2, [1, 2])]
    # A-X8-Stretch-exact: the static baseline 5 counts the zero-capacity edge; 1.2 * 5 < 6.
    h = graph(2, [(0, 1, 0, 5), (0, 1, 1, 6)])
    ceiled = ngc.KspOptions(2, sel(), max_path_cost=ngc.Stretch(1.2))
    assert walks(ngc.ksp(h, 0, 1, ceiled)) == [(5, [0])]


def test_headroom_without_a_residual_reads_one_or_zero():
    d = graph(4, [(0, 1, 1, 1), (0, 2, 5, 1), (1, 3, 1, 1), (2, 3, 5, 1)])
    z = graph(4, [(0, 1, 0, 1), (0, 2, 5, 1), (1, 3, 5, 1), (2, 3, 5, 1)])

    def path(g, prio):
        o = ngc.SpfOptions(FWD, sel(PREFER, prio, mp=False), targets=3)
        return edges(ngc.spf(g, 0, o).dag)

    assert path(d, HEAD) == [
        0,
        2,
    ]  # every positive edge reads 1: the EdgeId order decides
    assert path(d, ngc.EdgePriority.CAPACITY) == [1, 3]
    assert path(z, HEAD) == [1, 3]  # 0 on the zero-capacity edge


def test_g8_for_every_model():
    g = graph(3, S)
    r, p = place(g, 0, 1, COSTONLY, 5)
    assert (r.placed, r.remaining, r.termination) == (0, 5, ngc.Termination.EXHAUSTED)
    [rec] = p.records()
    assert (rec.delivered, edges(rec.bundle), rec.route_cost, rec.edges) == (
        0,
        [0],
        1,
        [],
    )
    r, p = place(g, 0, 1, config(ngc.CostOnly(), EB, LOSSY), 5)
    assert (r.placed, r.drops, r.termination) == (0, [(0, 5)], ngc.Termination.COMPLETE)
    for admission in (ngc.Admission.LOSSLESS, LOSSY):
        r, p = place(g, 0, 1, config(ngc.CostOnly(), W, admission), 5)
        assert (r.placed, r.drops) == (0, [])  # weight 0: no member, never a bypass
    assert place(g, 0, 1, config(ngc.CostOnly(), MF), 5)[0].placed == 0
    assert place(g, 0, 1, GLOBAL, 5)[0].placed == 5
    assert place(g, 0, 1, config(ngc.CapacityAware()), 5)[0].placed == 5


def test_equal_and_weighted_splits_over_a_zero_member():
    first = graph(2, [(0, 1, 0, 1), (0, 1, 10, 1)])
    last = graph(2, [(0, 1, 10, 1), (0, 1, 0, 1)])
    assert place(first, 0, 1, COSTONLY, 10)[0].placed == 0
    assert place(last, 0, 1, COSTONLY, 10)[0].placed == 1
    r, _ = place(first, 0, 1, config(ngc.CostOnly(), EB, LOSSY), 10)
    assert (r.placed, r.drops) == (5, [(0, 5)])
    assert place(first, 0, 1, config(ngc.CostOnly(), W), 10)[0].placed == 10
    # STALE-wcmp-zero-capacity and STALE-all-zero-weight.
    wcmp = graph(2, [(0, 1, 0, 1), (0, 1, 40, 1)])
    assert place(wcmp, 0, 1, config(ngc.CostOnly(), W, LOSSY), 40)[0].placed == 40
    one = graph(2, [(0, 1, 0, 1)])
    r, p = place(one, 0, 1, config(ngc.CostOnly(), W, LOSSY), 2)
    [rec] = p.records()
    assert (r.placed, r.drops, rec.offered, rec.delivered) == (0, [], 2, 0)
    # D-zero-service: both units dropped on the zero-capacity link.
    r, _ = place(one, 0, 1, config(ngc.CostOnly(), EB, LOSSY), 2)
    assert (r.placed, r.drops) == (0, [(0, 2)])


@pytest.mark.parametrize(
    ("caps", "lossless", "drops"),
    [((0, 10), 0, [(0, 5)]), ((10, 0), 1, [(1, 5)])],
    ids=["B3-zero-first", "B3-positive-first"],
)
def test_two_sources_one_behind_a_zero_link(caps, lossless, drops):
    # The pooled equal origin: a LOSSLESS total a is admitted only when RR(a) fits, so a source
    # behind the zero link caps it at its round-robin position: 0 when source 0 is behind it
    # (RR(1) gives source 0 the unit), 1 when source 1 is (RR(2) gives source 1 a unit). LOSSY
    # offers RR(10) = 5 and 5 and drops the share behind the zero link.
    g = graph(3, [(0, 2, caps[0], 1), (1, 2, caps[1], 1)])
    r, p = place(g, [0, 1], 2, COSTONLY, 10)
    assert (r.placed, r.termination) == (lossless, ngc.Termination.EXHAUSTED)
    assert [row.source for row in r.by_source] == [0, 1]
    assert len(p.records()) == 2
    r, _ = place(g, [0, 1], 2, config(ngc.CostOnly(), EB, LOSSY), 10)
    assert (r.placed, r.drops) == (5, drops)


def test_b3_delete_zero_and_sensitivity():
    g = graph(3, [(0, 2, 0, 1), (1, 2, 10, 1)])
    for elig in (g.eligibility(excluded_nodes=[0]), g.eligibility(excluded_edges=[0])):
        assert place(g, [0, 1], 2, COSTONLY, 10, elig)[0].placed == 10
        assert place(g, [0, 1], 2, GLOBAL, 10, elig)[0].placed == 10
    # B3-zero-first under the pooled equal origin admits 0; removing the zero link takes source 0
    # out of the program, so the deletion gains 10 (a negative loss); removing 1 -> 2 leaves 0.
    s = ngc.sensitivity(g, [0, 1], 2, ngc.QueryOptions(COSTONLY))
    assert (s.baseline.placed, [(e.edge, e.loss) for e in s.entries]) == (
        0,
        [(0, -10), (1, 0)],
    )
    # G8: a negative loss for the zero-capacity link.
    s8 = ngc.sensitivity(graph(3, S), 0, 1, ngc.QueryOptions(COSTONLY))
    assert [(e.edge, e.loss) for e in s8.entries] == [(0, -5), (1, 0), (2, 0)]


def test_capacity_aware_roster_and_stretch_baseline():
    x8 = graph(3, [(0, 2, 0, 1), (1, 2, 8, 3)])
    r, _ = place(x8, [0, 1], 2, config(ngc.CapacityAware()), 4)
    assert (r.placed, r.termination) == (
        2,
        ngc.Termination.NO_ROUTE,
    )  # A-X8: head 0 loses 2
    h = graph(2, [(0, 1, 0, 5), (0, 1, 1, 6)])
    r, p = place(h, 0, 1, config(ngc.CapacityAware(max_path_cost=ngc.Stretch(1.2))), 1)
    assert (r.placed, p.records()) == (0, [])


@pytest.mark.parametrize(
    ("caps", "bound"), [((0, 4), 0), ((4, 0), 1)], ids=["E-zero-first", "E-zero-last"]
)
def test_pinned_bundles_keep_their_zero_links(caps, bound):
    g = graph(2, [(0, 1, caps[0], 1), (0, 1, caps[1], 1)])
    cfg = config(ngc.Pinned([dag_of(g, [0]), dag_of(g, [1])]))
    r, p = place(g, 0, 1, cfg, math.inf)
    assert (r.placed, p.usable_bundles()) == (bound, 2)


def test_global_min_cut_lists_zero_capacity_links_leaving_the_closure():
    want = ngc.OutputOptions(min_cut=True)
    b = graph(3, [(0, 2, 0, 1), (1, 2, 10, 1)])
    r = ngc.capacity(b, [0, 1], 2, ngc.QueryOptions(GLOBAL, outputs=want))
    assert (r.placed, r.min_cut) == (10, [0, 1])
    s = graph(3, S)
    r = ngc.capacity(s, 0, 1, ngc.QueryOptions(GLOBAL, outputs=want))
    assert (r.placed, r.min_cut) == (5, [0, 1])
    masked = ngc.QueryOptions(
        GLOBAL, elig=s.eligibility(excluded_edges=[0]), outputs=want
    )
    assert ngc.capacity(s, 0, 1, masked).min_cut == [1]
