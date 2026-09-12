"""spf_to: reverse SPF toward one destination, with forced fan-out edges."""

from __future__ import annotations

import numpy as np
import pytest

import netgraph_core as ngc


def _graph(src, dst, cap, cost, n):
    return ngc.StrictMultiDiGraph.from_arrays(
        num_nodes=n,
        src=np.array(src, dtype=np.int32),
        dst=np.array(dst, dtype=np.int32),
        capacity=np.array(cap, dtype=np.float64),
        cost=np.array(cost, dtype=np.int64),
        ext_edge_ids=np.arange(len(src), dtype=np.int64),
    )


SEL = ngc.EdgeSelection(
    multi_edge=True, require_capacity=False, tie_break=ngc.EdgeTieBreak.DETERMINISTIC
)


@pytest.fixture
def fanout_graph():
    """P(0) -> S1(1) cost 0, P -> S2(2) cost 0; S1 -> T(4) cap 100; S2 -> M(3) cap 20; M -> T cap 60."""
    g = _graph(
        [0, 0, 1, 2, 3],
        [1, 2, 4, 3, 4],
        [1e6, 1e6, 100.0, 20.0, 60.0],
        [0, 0, 1, 1, 1],
        5,
    )
    algs = ngc.Algorithms(ngc.Backend.cpu())
    return g, algs, algs.build_graph(g)


def test_distances_to_dst_match_forward_spf(fanout_graph):
    g, algs, h = fanout_graph
    dist_to, dag = algs.spf_to(h, 4, selection=SEL)
    for s in range(g.num_nodes()):
        d, _ = algs.spf(h, s, 4, selection=SEL)
        assert dist_to[s] == d[4]
    assert dag.parent_offsets.shape == (g.num_nodes() + 1,)


def test_fanout_forces_all_sources_and_admission_scales_globally(fanout_graph):
    g, algs, h = fanout_graph
    dist, dag = algs.spf_to(h, 4, selection=SEL, fanout_edges=[0, 1])
    assert dist[1] == 1.0 and dist[2] == 2.0
    # S2 has P as parent although P->S2 is not on a shortest P->T path.
    assert list(dag.via_edges[dag.parent_offsets[2] : dag.parent_offsets[3]]) == [1]
    fg = ngc.FlowGraph(g)
    placed = fg.place(
        ngc.FlowIndex(0, 4, 0, 0),
        0,
        4,
        dag,
        100.0,
        ngc.FlowPlacement.EQUAL_BALANCED,
    )
    # Even 50/50 split; S2's 20-unit link admits 20 of 50, so the demand scales to 0.4.
    assert placed == pytest.approx(40.0)
    assert fg.edge_flow_view()[2] == pytest.approx(20.0)


def test_fanout_lossy_delivers_what_each_branch_carries(fanout_graph):
    g, algs, h = fanout_graph
    _, dag = algs.spf_to(h, 4, selection=SEL, fanout_edges=[0, 1])
    fg = ngc.FlowGraph(g)
    placed, drops = fg.place_with_drops(
        ngc.FlowIndex(0, 4, 0, 0),
        0,
        4,
        dag,
        100.0,
        ngc.FlowPlacement.EQUAL_BALANCED_LOSSY,
    )
    assert placed == pytest.approx(70.0)  # 50 via S1 + 20 via S2
    assert drops == [(3, pytest.approx(30.0))]


def test_without_fanout_the_pseudo_source_reaches_only_the_nearest_source(fanout_graph):
    g, algs, h = fanout_graph
    _, dag = algs.spf(h, 0, 4, selection=SEL)
    # T is reached only via S1 -> T (edge 2): S2's branch costs 2 and never
    # carries flow placed from P.
    assert list(dag.via_edges[dag.parent_offsets[4] : dag.parent_offsets[5]]) == [2]


def test_fanout_respects_masks_and_reachability(fanout_graph):
    g, algs, h = fanout_graph
    edge_mask = np.ones(g.num_edges(), dtype=bool)
    edge_mask[3] = False  # S2 -> M down
    dist, dag = algs.spf_to(
        h, 4, selection=SEL, edge_mask=edge_mask, fanout_edges=[0, 1]
    )
    assert np.isinf(dist[2])
    assert dag.parent_offsets[3] - dag.parent_offsets[2] == 0


def test_fanout_validation(fanout_graph):
    g, algs, h = fanout_graph
    with pytest.raises(ValueError):
        algs.spf_to(
            h, 4, selection=SEL, fanout_edges=[4]
        )  # M has an incoming DAG entry
    with pytest.raises(ValueError):
        algs.spf_to(h, 4, selection=SEL, fanout_edges=[99])
    with pytest.raises(ValueError):
        algs.spf_to(h, 42, selection=SEL)
