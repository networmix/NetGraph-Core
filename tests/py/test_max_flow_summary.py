"""Max-flow FlowSummary semantics: shapes and MinCut validity."""

from __future__ import annotations

import numpy as np

import netgraph_core as ngc


def _two_path_graph(build_graph):
    # S=0, A=1, B=2, T=3
    # Path1 cost=2 cap=5; Path2 cost=4 cap=3
    edges = [
        (0, 1, 1, 5, 0),
        (1, 3, 1, 5, 1),
        (0, 2, 2, 3, 2),
        (2, 3, 2, 3, 3),
    ]
    return build_graph(4, edges)


def test_flow_summary_shapes_and_min_cut_valid(
    build_graph, assert_edge_flows_shape, assert_valid_min_cut, algs, to_handle
):
    g = _two_path_graph(build_graph)
    total, summary = algs.max_flow(
        to_handle(g),
        0,
        3,
        flow_placement=ngc.FlowPlacement.PROPORTIONAL,
        shortest_path=False,
        with_edge_flows=True,
    )
    assert np.isclose(total, 8.0)
    assert_edge_flows_shape(g, summary, expected_present=True)
    assert_valid_min_cut(g, summary.min_cut)
    # Cost distribution checked in test_max_flow_cost_distribution.py


def test_flow_summary_without_edge_flows(
    build_graph, assert_edge_flows_shape, algs, to_handle
):
    g = _two_path_graph(build_graph)
    total, summary = algs.max_flow(
        to_handle(g),
        0,
        3,
        flow_placement=ngc.FlowPlacement.PROPORTIONAL,
        shortest_path=False,
        with_edge_flows=False,
    )
    assert np.isclose(total, 8.0)
    assert_edge_flows_shape(g, summary, expected_present=False)
