"""PredDAG.from_hops: explicit routes whose hops may use several parallel edges."""

from __future__ import annotations

import numpy as np
import pytest

import netgraph_core as ngc


def _graph():
    # Four parallel 10-unit links 0 -> 1 (ids 0-3), 1 -> 0 back (id 4) and
    # 1 -> 2 (id 5): edge ids follow the graph's (cost, src, dst) ordering.
    return ngc.StrictMultiDiGraph.from_arrays(
        3,
        src=np.array([0, 0, 0, 0, 1, 1], dtype=np.int32),
        dst=np.array([1, 1, 1, 1, 2, 0], dtype=np.int32),
        capacity=np.array([10.0, 10.0, 10.0, 10.0, 100.0, 1.0], dtype=np.float64),
        cost=np.array([1, 1, 1, 1, 1, 1], dtype=np.int64),
        ext_edge_ids=np.arange(6, dtype=np.int64),
    )


def test_from_hops_layout_and_order():
    g = _graph()
    dag = ngc.PredDAG.from_hops(g, [[3, 1], [5]])
    assert list(dag.parent_offsets) == [0, 0, 2, 3]
    assert list(dag.parents) == [0, 0, 1]
    assert list(dag.via_edges) == [3, 1, 5]
    single = ngc.PredDAG.from_edges(g, [2, 5])
    assert list(single.via_edges) == [2, 5]


def test_from_hops_bundle_carries_every_parallel_edge():
    g = _graph()
    fg = ngc.FlowGraph(g)
    one = ngc.PredDAG.from_edges(g, [0, 5])
    bundle = ngc.PredDAG.from_hops(g, [[0, 1, 2, 3], [5]])
    idx_a = ngc.FlowIndex(0, 2, 0, 1)
    idx_b = ngc.FlowIndex(0, 2, 0, 2)
    assert fg.place(
        idx_a, 0, 2, one, 100.0, ngc.FlowPlacement.PROPORTIONAL
    ) == pytest.approx(10.0)
    fg.remove(idx_a)
    assert fg.place(
        idx_b, 0, 2, bundle, 100.0, ngc.FlowPlacement.PROPORTIONAL
    ) == pytest.approx(40.0)


@pytest.mark.parametrize(
    "hops, message",
    [
        ([], "non-empty"),
        ([[0], []], "empty"),
        ([[0, 9]], "out of range"),
        ([[0, 0]], "repeated"),
        ([[4, 5]], "not parallel"),  # 1 -> 0 and 1 -> 2
        ([[0], [5], [4]], "contiguous"),  # 0 -> 1 -> 2, then 1 -> 0
        ([[0], [5]], None),  # valid: 0 -> 1 -> 2
    ],
)
def test_from_hops_validation(hops, message):
    g = _graph()
    if message is None:
        ngc.PredDAG.from_hops(g, hops)
        return
    with pytest.raises(ValueError, match=message):
        ngc.PredDAG.from_hops(g, hops)


def test_from_hops_rejects_revisited_node():
    g = _graph()
    with pytest.raises(ValueError, match="simple path"):
        ngc.PredDAG.from_hops(g, [[0], [4]])  # 0 -> 1 -> 0
