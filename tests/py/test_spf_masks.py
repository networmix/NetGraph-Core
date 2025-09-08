"""SPF masks: behavior of node_mask and edge_mask in shortest-path routines."""

import numpy as np

import netgraph_core as ngc


def test_node_mask_blocks_path(build_graph, algs, to_handle):
    # A=0 -> B=1 -> C=2, all cost=1
    edges = [
        (0, 1, 1, 1, 0),
        (1, 2, 1, 1, 1),
    ]
    g = build_graph(3, edges)
    node_mask = np.array([True, False, True], dtype=bool)
    sel = ngc.EdgeSelection(
        multi_edge=True,
        require_capacity=False,
        tie_break=ngc.EdgeTieBreak.DETERMINISTIC,
    )
    dist, dag = algs.spf(to_handle(g), 0, 2, selection=sel, node_mask=node_mask)
    # Destination unreachable because B is masked out
    assert np.isinf(dist[2])
    # No parents recorded for 2
    off = np.asarray(dag.parent_offsets)
    assert off[2] == off[3]


def test_edge_mask_filters_parallel_edges(build_graph, algs, to_handle):
    # A=0 -> B=1 with two parallel edges of different costs; only cheap edge should be allowed
    edges = [
        (0, 1, 1, 1, 10),
        (0, 1, 2, 1, 20),
    ]
    g = build_graph(2, edges)
    edge_mask = np.array([True, False], dtype=bool)
    sel = ngc.EdgeSelection(
        multi_edge=True,
        require_capacity=False,
        tie_break=ngc.EdgeTieBreak.DETERMINISTIC,
    )
    dist, dag = algs.spf(to_handle(g), 0, 1, selection=sel, edge_mask=edge_mask)
    assert np.isclose(dist[1], 1.0)
    # Only one predecessor edge should appear
    offsets = np.asarray(dag.parent_offsets)
    via = np.asarray(dag.via_edges)
    start, end = int(offsets[1]), int(offsets[2])
    assert end - start == 1
    # Map via edge id -> ensure it's a valid edge id
    e = int(via[start])
    assert 0 <= e < g.num_edges()
