"""PredDAG, Path, PathSelection, the enums and Eligibility.intersect through the bindings
(APPROACH 4, 5, 6): validation in the documented order on unnarrowed values, read-only views,
structural equality."""

from __future__ import annotations

import gc

import numpy as np
import pytest

import netgraph_core as ngc

L64 = 1 << 64


def build(n, edges, q=1.0):
    return ngc.StrictMultiDiGraph(
        n,
        [e[0] for e in edges],
        [e[1] for e in edges],
        [e[2] for e in edges],
        [e[3] for e in edges],
        q,
    )


# ids by (src, cost, dst, input): 0: 0->1, 1: 0->1, 2: 0->2 (cost 3), 3: 1->2, 4: 2->0
G_EDGES = [(0, 1, 5, 1), (0, 1, 5, 1), (1, 2, 5, 1), (0, 2, 5, 3), (2, 0, 5, 1)]


def test_from_csr_round_trip_and_views():
    g = build(3, G_EDGES)
    dag = ngc.PredDAG.from_csr(g, [0, 0, 2, 4], [0, 0, 0, 1], [0, 1, 2, 3])
    assert dag.graph_uid() == g.uid()
    assert list(dag.parent_offsets()) == [0, 0, 2, 4]
    assert list(dag.parents()) == [0, 0, 0, 1]
    assert list(dag.via_edges()) == [0, 1, 2, 3]
    view = dag.via_edges()
    assert view.dtype == np.int32
    with pytest.raises(ValueError):
        view[0] = 9
    del dag
    gc.collect()
    assert list(view) == [0, 1, 2, 3]  # the view keeps the DAG alive
    same = ngc.PredDAG.from_csr(
        g, np.array([0, 0, 2, 4]), np.array([0, 0, 0, 1]), np.array([0, 1, 2, 3])
    )
    again = ngc.PredDAG.from_csr(g, [0, 0, 2, 4], [0, 0, 0, 1], [0, 1, 2, 3])
    assert same == again and hash(same) == hash(again)
    assert same != ngc.PredDAG.from_csr(g, [0, 0, 0, 0], [], [])


@pytest.mark.parametrize(
    ("offsets", "parents", "via", "message"),
    [
        ([0, 0, 0], [], [], "offsets has length 3"),
        ([1, 0, 2, 1], [], [], r"offsets\[0\] = 1"),
        ([L64, 0, 0, 0], [], [], r"offsets\[0\] = 18446744073709551616"),
        ([0, 2, 1, 2], [0, 0], [0, 1], "not monotone at 2"),
        ([0, L64, 0, 0], [], [], "not monotone at 2"),
        ([0, 0, 2, 2], [0], [0, 1], "parents has 1"),
        ([0, 0, 2, 3], [9, 0], [0, 1], r"offsets\[N\] = 3"),
        ([0, 0, 2, 2], [3, 0], [0, 1], r"parents\[0\] = 3"),
        ([0, 0, 2, 2], [-L64, 0], [0, 1], r"parents\[0\] = -18446744073709551616"),
        ([0, 0, 2, 2], [0, 0], [0, L64], r"via\[1\] = 18446744073709551616"),
        ([0, 0, 2, 2], [1, 0], [0, 1], r"joins \(0, 1\), not \(1, 1\)"),
        ([0, 0, 2, 2], [0, 0], [1, 1], "duplicate entry"),
        ([0, 1, 2, 3], [2, 0, 1], [4, 0, 3], "cycle"),
        ([0, 0, 2, 2], [0, 0], [1, 0], "not in ascending via EdgeId"),
    ],
)
def test_from_csr_validation_order(offsets, parents, via, message):
    g = build(3, G_EDGES)
    with pytest.raises(ValueError, match=message):
        ngc.PredDAG.from_csr(g, offsets, parents, via)


@pytest.mark.parametrize(
    "offsets",
    [
        [0, 2**100, 2**99],
        (0, 2**100, 2**99),
        [0, 2**64 - 1, 2**64 - 2],
        np.array([0, 2**64 - 1, 2**64 - 2], dtype=np.uint64),
        [0, 100, 99],
    ],
)
def test_from_csr_compares_offsets_beyond_int64_exactly(offsets):
    # review-1b D2: two offsets above int64 compare by value, so a decreasing pair is the
    # monotonicity error at its index, before the offsets[N] count check.
    g = build(2, [(0, 1, 1, 1)])
    with pytest.raises(ValueError, match="not monotone at 2"):
        ngc.PredDAG.from_csr(g, offsets, [], [])
    increasing = [0, 2**99, 2**100]
    with pytest.raises(
        ValueError, match=r"offsets\[N\] = 1267650600228229401496703205376"
    ):
        ngc.PredDAG.from_csr(g, increasing, [], [])


def test_from_csr_type_and_shape_defects():
    g = build(3, G_EDGES)
    with pytest.raises(ValueError, match="1-D"):
        ngc.PredDAG.from_csr(g, np.zeros((2, 2), dtype=np.int64), [], [])
    with pytest.raises(ValueError, match="integer"):
        ngc.PredDAG.from_csr(g, [0, 0, 0, 0.0], [], [])
    with pytest.raises(ValueError, match="integer"):
        ngc.PredDAG.from_csr(g, [0, 0, 0, True], [], [])


def test_path_cost_hops_and_to_dag():
    g = build(3, G_EDGES)
    p = ngc.Path(g, [[1, 0], [3]])
    assert (p.src(), p.dst(), p.cost()) == (0, 2, 2)
    assert p.hops() == [[1, 0], [3]]  # the caller's member order is kept
    assert p.nodes() == [0, 1, 2]
    dag = p.to_dag()
    assert list(dag.parent_offsets()) == [0, 0, 2, 3]
    assert list(dag.via_edges()) == [0, 1, 3]
    assert ngc.Path(g, [[1, 0], [3]]) == p
    assert ngc.Path(g, [[0, 1], [3]]) != p
    assert hash(ngc.Path(g, [[1, 0], [3]])) == hash(p)


@pytest.mark.parametrize(
    ("hops", "message"),
    [
        ([[], [7]], "member 7 out of edge range"),
        ([[L64]], "member 18446744073709551616 out of edge range"),
        ([[0, 0], []], "duplicate member 0"),
        ([], "hop list is empty"),
        ([[0], []], "hop 1 is empty"),
        ([[0, 2]], r"do not share \(u, v\)"),
        ([[0], [2]], "not contiguous"),
        ([[2], [4]], "node 0 is repeated"),
    ],
)
def test_path_validation_order(hops, message):
    g = build(3, G_EDGES)
    with pytest.raises(ValueError, match=message):
        ngc.Path(g, hops)


def test_path_selection_and_enums():
    sel = ngc.PathSelection(
        ngc.EdgeTieBreak.PREFER_HIGHER_PRIORITY, ngc.EdgePriority.HEADROOM
    )
    assert (sel.multipath, sel.multi_edge) == (True, True)
    assert sel.tie_break == ngc.EdgeTieBreak.PREFER_HIGHER_PRIORITY
    assert sel == ngc.PathSelection(
        ngc.EdgeTieBreak.PREFER_HIGHER_PRIORITY, ngc.EdgePriority.HEADROOM, True, True
    )
    assert sel != ngc.PathSelection(
        ngc.EdgeTieBreak.PREFER_HIGHER_PRIORITY,
        ngc.EdgePriority.HEADROOM,
        multipath=False,
    )
    assert "PREFER_HIGHER_PRIORITY" in repr(sel)
    with pytest.raises(TypeError):
        ngc.PathSelection()  # type: ignore[call-arg]  # tie_break and priority are explicit
    assert {m.name for m in ngc.FlowPlacement.__members__.values()} == {
        "MAX_FLOW",
        "WEIGHTED",
        "EQUAL_BALANCED",
    }
    assert set(ngc.Admission.__members__) == {"LOSSLESS", "LOSSY"}
    assert set(ngc.EdgePriority.__members__) == {"RESIDUAL", "CAPACITY", "HEADROOM"}
    assert set(ngc.EdgeTieBreak.__members__) == {
        "DETERMINISTIC",
        "PREFER_HIGHER_PRIORITY",
    }


def test_eligibility_intersect_returns_the_operand_object():
    g = build(3, G_EDGES)
    everything = g.eligibility()
    a = g.eligibility(excluded_edges=[0])
    b = g.eligibility(excluded_edges=[1])
    ab = g.eligibility(excluded_edges=[0, 1])
    assert ngc.Eligibility.intersect(everything, a) is a
    assert ngc.Eligibility.intersect(a, everything) is a
    assert ngc.Eligibility.intersect(a, a) is a
    assert ngc.Eligibility.intersect(a, ab) is ab
    x = ngc.Eligibility.intersect(a, b)
    y = ngc.Eligibility.intersect(a, b)
    assert x.uid() not in (0, a.uid(), b.uid(), y.uid())
    assert list(x.edge_allowed()) == list(ab.edge_allowed())
    other = build(3, G_EDGES)
    with pytest.raises(ValueError, match="different graphs"):
        ngc.Eligibility.intersect(a, other.eligibility())
