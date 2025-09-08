from __future__ import annotations

import netgraph_core as ngc


def test_resolve_no_split_parallel_edges(
    build_graph, assert_paths_concrete, algs, to_handle
):
    # Graph: 0->1 has two equal-cost edges; 1->2 single edge
    edges = [
        (0, 1, 1.0, 1.0),
        (0, 1, 1.0, 1.0),
        (1, 2, 1.0, 1.0),
    ]
    g = build_graph(3, edges)
    sel = ngc.EdgeSelection(
        multi_edge=True,
        require_capacity=False,
        tie_break=ngc.EdgeTieBreak.DETERMINISTIC,
    )
    dist, dag = algs.spf(to_handle(g), 0, 2, selection=sel)
    paths = dag.resolve_to_paths(0, 2, split_parallel_edges=False)
    # Expect one path grouping both parallel edges at hop 0->1
    assert len(paths) == 1
    assert_paths_concrete(paths, 0, 2, expect_split=False)


def test_resolve_with_split_parallel_edges_cartesian(
    build_graph, assert_paths_concrete, algs, to_handle
):
    # Graph: 0->1 two edges, 1->2 two edges (all equal cost)
    edges = [
        (0, 1, 1.0, 1.0),
        (0, 1, 1.0, 1.0),
        (1, 2, 1.0, 1.0),
        (1, 2, 1.0, 1.0),
    ]
    g = build_graph(3, edges)
    sel = ngc.EdgeSelection(
        multi_edge=True,
        require_capacity=False,
        tie_break=ngc.EdgeTieBreak.DETERMINISTIC,
    )
    dist, dag = algs.spf(to_handle(g), 0, 2, selection=sel)
    paths = dag.resolve_to_paths(0, 2, split_parallel_edges=True)
    # Expect 2*2 = 4 concrete paths
    assert len(paths) == 4
    assert_paths_concrete(paths, 0, 2, expect_split=True)


def test_resolve_max_paths_cap(build_graph, assert_paths_concrete, algs, to_handle):
    edges = [
        (0, 1, 1.0, 1.0),
        (0, 1, 1.0, 1.0),
        (1, 2, 1.0, 1.0),
        (1, 2, 1.0, 1.0),
    ]
    g = build_graph(3, edges)
    sel = ngc.EdgeSelection(
        multi_edge=True,
        require_capacity=False,
        tie_break=ngc.EdgeTieBreak.DETERMINISTIC,
    )
    dist, dag = algs.spf(to_handle(g), 0, 2, selection=sel)
    paths = dag.resolve_to_paths(0, 2, split_parallel_edges=True, max_paths=3)
    assert len(paths) == 3
    assert_paths_concrete(paths, 0, 2, expect_split=True)
