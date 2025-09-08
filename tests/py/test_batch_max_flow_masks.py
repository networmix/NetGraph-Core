"""Batch max-flow with node/edge masks and per-edge flow outputs."""

from __future__ import annotations

import numpy as np

import netgraph_core as ngc


def test_batch_max_flow_with_masks_and_edge_flows(square1_graph, algs, to_handle):
    g = square1_graph
    pairs = np.array([[0, 2], [0, 3]], dtype=np.int32)
    node_masks = [
        np.array([True, True, True, True], dtype=bool),
        np.array([True, True, True, False], dtype=bool),  # block node 3 for second pair
    ]
    # Supply a proper mask for each pair; the first pair allows all edges
    edge_masks = [
        np.ones(g.num_edges(), dtype=bool),
        np.array([True, True, False, False], dtype=bool),  # block the alt path edges
    ]
    out = algs.batch_max_flow(
        to_handle(g),
        pairs,
        flow_placement=ngc.FlowPlacement.PROPORTIONAL,
        shortest_path=False,
        with_edge_flows=True,
        node_masks=node_masks,
        edge_masks=edge_masks,  # type: ignore[arg-type]
    )
    assert len(out) == 2
    # Second pair should be zero due to masks blocking all routes
    assert out[1].total_flow == 0.0
    # Edge flows length matches num_edges for with_edge_flows=True
    assert np.asarray(out[0].edge_flows).shape == (g.num_edges(),)
