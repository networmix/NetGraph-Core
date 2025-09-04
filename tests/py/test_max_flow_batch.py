import numpy as np

import netgraph_core as ngc


def test_batch_max_flow_matches_individuals_proportional(square1_graph):
    g = square1_graph
    pairs = np.array([[0, 2], [0, 3]], dtype=np.int32)
    out = ngc.batch_max_flow(
        g,
        pairs,
        flow_placement=ngc.FlowPlacement.PROPORTIONAL,
        shortest_path=False,
        with_edge_flows=False,
    )
    # Individual calculations
    v1, _ = ngc.max_flow(
        g,
        0,
        2,
        flow_placement=ngc.FlowPlacement.PROPORTIONAL,
        shortest_path=False,
        with_edge_flows=False,
    )
    v2, _ = ngc.max_flow(
        g,
        0,
        3,
        flow_placement=ngc.FlowPlacement.PROPORTIONAL,
        shortest_path=False,
        with_edge_flows=False,
    )
    assert len(out) == 2
    assert np.isclose(out[0].total_flow, v1)
    assert np.isclose(out[1].total_flow, v2)


def test_batch_max_flow_matches_individuals_equal_balanced(line1_graph):
    g = line1_graph
    pairs = np.array([[0, 2], [1, 2]], dtype=np.int32)
    out = ngc.batch_max_flow(
        g,
        pairs,
        flow_placement=ngc.FlowPlacement.EQUAL_BALANCED,
        shortest_path=False,
        with_edge_flows=False,
    )
    v1, _ = ngc.max_flow(
        g,
        0,
        2,
        flow_placement=ngc.FlowPlacement.EQUAL_BALANCED,
        shortest_path=False,
        with_edge_flows=False,
    )
    v2, _ = ngc.max_flow(
        g,
        1,
        2,
        flow_placement=ngc.FlowPlacement.EQUAL_BALANCED,
        shortest_path=False,
        with_edge_flows=False,
    )
    assert len(out) == 2
    assert np.isclose(out[0].total_flow, v1)
    assert np.isclose(out[1].total_flow, v2)
