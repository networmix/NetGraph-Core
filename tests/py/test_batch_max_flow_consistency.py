"""Batch max-flow consistency with per-pair max_flow results."""

import numpy as np

import netgraph_core as ngc


def test_batch_max_flow_matches_individuals_proportional(
    square1_graph, algs, to_handle
):
    g = square1_graph
    pairs = np.array([[0, 2], [0, 3]], dtype=np.int32)
    out = algs.batch_max_flow(
        to_handle(g),
        pairs,
        flow_placement=ngc.FlowPlacement.PROPORTIONAL,
        shortest_path=False,
        with_edge_flows=False,
    )
    # Individual calculations
    v1, _ = algs.max_flow(
        to_handle(g),
        0,
        2,
        flow_placement=ngc.FlowPlacement.PROPORTIONAL,
        shortest_path=False,
        with_edge_flows=False,
    )
    v2, _ = algs.max_flow(
        to_handle(g),
        0,
        3,
        flow_placement=ngc.FlowPlacement.PROPORTIONAL,
        shortest_path=False,
        with_edge_flows=False,
    )
    assert len(out) == 2
    assert np.isclose(out[0].total_flow, v1)
    assert np.isclose(out[1].total_flow, v2)


def test_batch_max_flow_matches_individuals_equal_balanced(
    line1_graph, algs, to_handle
):
    g = line1_graph
    pairs = np.array([[0, 2], [1, 2]], dtype=np.int32)
    out = algs.batch_max_flow(
        to_handle(g),
        pairs,
        flow_placement=ngc.FlowPlacement.EQUAL_BALANCED,
        shortest_path=False,
        with_edge_flows=False,
    )
    v1, _ = algs.max_flow(
        to_handle(g),
        0,
        2,
        flow_placement=ngc.FlowPlacement.EQUAL_BALANCED,
        shortest_path=False,
        with_edge_flows=False,
    )
    v2, _ = algs.max_flow(
        to_handle(g),
        1,
        2,
        flow_placement=ngc.FlowPlacement.EQUAL_BALANCED,
        shortest_path=False,
        with_edge_flows=False,
    )
    assert len(out) == 2
    assert np.isclose(out[0].total_flow, v1)
    assert np.isclose(out[1].total_flow, v2)
