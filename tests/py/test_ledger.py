"""The integer ledger FlowGraph through the bindings (APPROACH 8): baselines on unnarrowed
values, snapshots and records, FlowIds, bottleneck, and the lifetime table driven by the
private single-bundle Pinned entry of sub-PR 1b (``_netgraph_core._PinnedFlow``)."""

from __future__ import annotations

import copy
import gc
import math

import _netgraph_core as ext
import numpy as np
import pytest

import netgraph_core as ngc

L62 = 1 << 62
INF = math.inf


def build(n, edges, q=1.0):
    return ngc.StrictMultiDiGraph(
        n,
        [e[0] for e in edges],
        [e[1] for e in edges],
        [e[2] for e in edges],
        [e[3] for e in edges],
        q,
    )


def dag_of(g, edges):
    """The PredDAG of an edge set (ascending via per node)."""
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


def pinned(fg, src, dst, edges, placement, admission=None, priority=None, elig=None):
    g = fg.graph()
    return ext._PinnedFlow(
        fg,
        src,
        dst,
        dag_of(g, edges),
        placement,
        admission if admission is not None else ngc.Admission.LOSSLESS,
        priority if priority is not None else ngc.EdgePriority.RESIDUAL,
        elig,
    )


EB = ngc.FlowPlacement.EQUAL_BALANCED
MF = ngc.FlowPlacement.MAX_FLOW
W = ngc.FlowPlacement.WEIGHTED
LOSSY = ngc.Admission.LOSSY


def test_baseline_validation_order_and_snapshots():
    g = build(2, [(0, 1, 5, 1), (1, 0, 3, 1)])
    fg = ngc.FlowGraph(g)
    assert list(fg.residual()) == [5, 3]
    assert fg.residual().dtype == np.int64
    assert list(fg.allocated()) == [0, 0]
    assert list(fg.capacity_units()) == [5, 3]
    assert fg.records() == []
    assert fg.graph() is g
    for bad, msg in (
        ([6, 0], r"residual\[0\] = 6 is outside \[0, capacity_units\[0\] = 5\]"),
        ([1], "residual has length 1, expected 2"),
        ([-1, 0], r"residual\[0\] = -1"),
        ([L62 * 4, 0], r"residual\[0\] = 18446744073709551616"),
        ([2.5, 0], "integer"),
        ([True, 0], "integer"),
        (np.zeros((1, 2), dtype=np.int64), "1-D"),
    ):
        with pytest.raises(ValueError, match=msg):
            ngc.FlowGraph(g, bad)
    fg2 = ngc.FlowGraph(graph=g, residual=np.array([1, 2], dtype=np.uint64))
    assert list(fg2.residual()) == [1, 2]
    assert fg2.uid() != fg.uid()


def test_lifetime_bound_rows_and_epochs():
    g = build(2, [(0, 1, 5, 1)])
    fg = ngc.FlowGraph(g)
    p = pinned(fg, 0, 1, [0], EB)
    r = p.place(0)
    assert (r["termination"], r["routes_changed"], fg.records()) == (
        "complete",
        False,
        [],
    )
    assert fg._state()["next_flow_id"] == 0
    with pytest.raises(ValueError, match="bound"):
        fg.reset()
    with pytest.raises(ValueError, match="bound"):
        fg.copy()
    with pytest.raises(ValueError, match="bound"):
        copy.copy(fg)
    r = p.place(2)
    assert (r["placed"], r["routes_changed"], p.flow_id()[1]) == (2, True, 0)
    r = p.place(INF)
    assert (r["placed"], r["termination"], r["remaining"]) == (3, "exhausted", INF)
    with pytest.raises(ValueError, match=r"after \+inf"):
        p.place(1)
    with pytest.raises(ValueError, match=r"after \+inf"):
        p.place(0)
    old = p.flow_id()
    p.remove()
    assert list(fg.allocated()) == [0]
    with pytest.raises(ValueError, match="retired or unknown"):
        fg.record(old)
    r = p.place(1)
    assert (r["placed"], p.flow_id()[1]) == (1, 1)
    p.remove()
    p.remove()
    assert (p.requested(), p.infinite()) == (0, False)
    assert fg._state() == {"next_flow_id": 2, "bound": [0], "orphans": []}


def test_orphans_from_destruction_reset_and_copy():
    g = build(2, [(0, 1, 5, 1)])
    fg = ngc.FlowGraph(g)
    q = pinned(fg, 0, 1, [0], EB)
    q.place(2)
    del q
    gc.collect()  # destroyed without remove(): an orphan
    assert fg._state()["orphans"] == [0]
    assert fg.records()[0].delivered == 2
    p = pinned(fg, 0, 1, [0], EB)
    assert p.place(5)["placed"] == 3
    p.release()
    with pytest.raises(Exception, match="released"):
        p.place(1)
    c = fg.copy()
    assert c.uid() != fg.uid()
    assert [r.flow_id for r in c.records()] == [(c.uid(), 0), (c.uid(), 1)]
    assert list(c.allocated()) == [5]
    d = copy.deepcopy(fg)
    assert list(d.allocated()) == [5]
    fg.reset([4])
    assert (fg.records(), list(fg.residual())) == ([], [4])
    assert fg._state() == {"next_flow_id": 2, "bound": [], "orphans": []}
    assert len(c.records()) == 2  # the copy is independent


def test_policy_keeps_the_ledger_alive():
    g = build(2, [(0, 1, 5, 1)])
    fg = ngc.FlowGraph(g)
    p = pinned(fg, 0, 1, [0], MF)
    del fg
    gc.collect()
    assert p.place(2)["placed"] == 2
    assert p.placed_demand() == 2


def test_records_fields_and_counters_beyond_64_bits():
    g = build(4, [(0, 1, 10, 1), (0, 2, 10, 1), (1, 3, 1, 1), (2, 3, 10, 1)])
    fg = ngc.FlowGraph(g)
    p = pinned(fg, 0, 3, range(4), EB, LOSSY)
    r = p.place(10)
    assert (r["placed"], r["termination"], r["drops"]) == (6, "complete", [(2, 4)])
    rec = fg.records()[0]
    assert rec.flow_id == (fg.uid(), 0)
    assert (rec.owner, rec.sources, rec.destinations) == (0, [0], [3])
    assert (rec.placement, rec.admission) == (EB, LOSSY)
    assert list(rec.bundle.via_edges()) == [0, 1, 2, 3]
    assert (rec.route_cost, rec.offered, rec.delivered) == (2, 10, 6)
    assert rec.edges == [(0, 5), (1, 5), (2, 1), (3, 5)]
    assert rec.drops == [(2, 4)]
    assert fg.record(rec.flow_id).delivered == 6
    # Lossless records: no offered; a route without a walk at creation: route_cost None.
    g2 = build(2, [(0, 1, 1, 1)])
    fg2 = ngc.FlowGraph(g2)
    pinned(fg2, 0, 1, [0], MF).place(1)
    assert fg2.records()[0].offered is None
    # Every unit dropped on a saturated link, nine times: the counters pass 2^64.
    sat = pinned(fg2, 0, 1, [0], EB, LOSSY)
    for _ in range(9):
        sat.place(L62 - 1)
    lossy = fg2.records()[1]
    assert lossy.offered == 9 * (L62 - 1)
    assert lossy.drops == [(0, 9 * (L62 - 1))]
    assert sat.requested() == 9 * (L62 - 1)
    assert lossy.offered > 1 << 64


def test_empty_split_graph_keeps_the_route_cost():
    # Weighted Residual with every route edge at residual 0 under another owner's flow: the
    # record exists at zero delivery with an empty split graph; route_cost is the bundle's
    # minimum walk cost (APPROACH 8), not a split-graph value.
    g = build(2, [(0, 1, 2, 3)])
    fg = ngc.FlowGraph(g)
    pinned(fg, 0, 1, [0], MF).place(2)
    w = pinned(fg, 0, 1, [0], W)
    r = w.place(1)
    assert (r["placed"], r["termination"], r["routes_changed"]) == (
        0,
        "exhausted",
        True,
    )
    rec = fg.record(w.flow_id())
    assert (rec.route_cost, rec.delivered, rec.edges) == (3, 0, [])
    assert w.plan()["members"] == []


def test_flow_ids_of_other_ledgers_and_bad_shapes():
    g = build(2, [(0, 1, 5, 1)])
    a, b = ngc.FlowGraph(g), ngc.FlowGraph(g)
    p = pinned(a, 0, 1, [0], MF)
    p.place(1)
    with pytest.raises(ValueError, match="another ledger"):
        b.record(p.flow_id())
    with pytest.raises(ValueError, match="tuple"):
        a.record(0)
    with pytest.raises(ValueError, match="retired or unknown"):
        a.record((a.uid(), L62 * 8))


def test_requests_are_integer_units_or_inf():
    g = build(2, [(0, 1, 5, 1)])
    fg = ngc.FlowGraph(g)
    p = pinned(fg, 0, 1, [0], EB)
    for bad, msg in (
        (1.5, "integer units or"),
        (-1, "out of range"),
        (L62, "out of range"),
        (True, "integer"),
        (-INF, "integer units or"),
        (float("nan"), "integer units or"),
    ):
        with pytest.raises(ValueError, match=msg):
            p.place(bad)
    assert fg.records() == []
    assert p.place(np.int64(2))["placed"] == 2


def test_bottleneck_versus_min_cut():
    g = build(4, [(0, 1, 1, 1), (0, 2, 2, 3), (1, 2, 1, 1), (2, 3, 2, 1)])
    fg = ngc.FlowGraph(g)
    p = pinned(fg, 0, 3, range(4), MF)
    assert p.place(INF)["placed"] == 2
    inv = [0] * 4
    for e, pos in enumerate(g.input_index()):
        inv[pos] = e
    assert fg.bottleneck(0) == sorted([inv[0], inv[3]])
    assert fg.bottleneck([0, 0]) == sorted([inv[0], inv[3]])
    assert fg.bottleneck([0], g.eligibility(excluded_nodes=[0])) == []
    with pytest.raises(ValueError, match="empty"):
        fg.bottleneck([])
    with pytest.raises(ValueError, match="out of node range"):
        fg.bottleneck([L62 * 8])
    with pytest.raises(ValueError, match="foreign"):
        fg.bottleneck([0], build(4, [(0, 1, 1, 1)]).eligibility(excluded_edges=[0]))


def test_pinned_construction_errors():
    g = build(3, [(0, 1, 5, 1), (2, 1, 5, 1)])
    fg = ngc.FlowGraph(g)
    with pytest.raises(ValueError, match="malformed bundle"):
        pinned(fg, 0, 2, [0, 1], MF)
    with pytest.raises(ValueError, match="Lossy"):
        pinned(fg, 0, 1, [0], MF, LOSSY)
    with pytest.raises(ValueError, match="Headroom"):
        pinned(fg, 0, 1, [0], W, priority=ngc.EdgePriority.HEADROOM)
    with pytest.raises(ValueError, match="both endpoint sets"):
        pinned(fg, 1, 1, [0], MF)
    other = build(3, [(0, 1, 5, 1), (2, 1, 5, 1)])
    with pytest.raises(ValueError, match="foreign PredDAG"):
        ext._PinnedFlow(
            fg,
            0,
            1,
            dag_of(other, [0]),
            MF,
            ngc.Admission.LOSSLESS,
            ngc.EdgePriority.RESIDUAL,
            None,
        )
    assert fg._state()["bound"] == []
