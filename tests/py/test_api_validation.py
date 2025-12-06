"""API validation: argument shapes, bounds, and error handling for frontdoor APIs."""

from __future__ import annotations

import numpy as np
import pytest

import netgraph_core as ngc


def _make_graph(num_nodes, src, dst, cap, cost):
    """Helper to build graph with auto-generated ext_edge_ids."""
    ext_edge_ids = np.arange(len(src), dtype=np.int64)
    return ngc.StrictMultiDiGraph.from_arrays(
        num_nodes, src, dst, cap, cost, ext_edge_ids
    )


def _small_graph(build_graph):
    # 0 -> 1 with cost=1, cap=1
    return build_graph(2, [(0, 1, 1.0, 1.0, 0)])


def test_spf_invalid_node_mask_shape_raises(build_graph, algs, to_handle):
    g = _small_graph(build_graph)
    sel = ngc.EdgeSelection(
        multi_edge=True,
        require_capacity=False,
        tie_break=ngc.EdgeTieBreak.DETERMINISTIC,
    )
    bad_mask = np.array([True], dtype=bool)  # wrong length
    with pytest.raises(TypeError):
        algs.spf(to_handle(g), 0, 1, selection=sel, node_mask=bad_mask)


def test_spf_invalid_edge_mask_shape_raises(build_graph, algs, to_handle):
    g = _small_graph(build_graph)
    sel = ngc.EdgeSelection(
        multi_edge=True,
        require_capacity=False,
        tie_break=ngc.EdgeTieBreak.DETERMINISTIC,
    )
    bad_edge_mask = np.array([True, False, True], dtype=bool)  # wrong length
    with pytest.raises(TypeError):
        algs.spf(to_handle(g), 0, 1, selection=sel, edge_mask=bad_edge_mask)


def test_spf_invalid_source_node_raises(build_graph, algs, to_handle):
    g = _small_graph(build_graph)
    with pytest.raises(ValueError):
        algs.spf(to_handle(g), 10)  # out of range src


def test_max_flow_invalid_masks_shape_raises(build_graph, algs, to_handle):
    g = _small_graph(build_graph)
    bad_edge_mask = np.array([True, False, True], dtype=bool)
    with pytest.raises(TypeError):
        algs.max_flow(to_handle(g), 0, 1, edge_mask=bad_edge_mask)

    bad_node_mask = np.array([True], dtype=bool)
    with pytest.raises(TypeError):
        algs.max_flow(to_handle(g), 0, 1, node_mask=bad_node_mask)


def test_ksp_invalid_k_and_dst_raises(build_graph, algs, to_handle):
    g = _small_graph(build_graph)
    with pytest.raises(ValueError):
        algs.ksp(to_handle(g), 0, 1, k=0, max_cost_factor=None, unique=True)
    with pytest.raises(ValueError):
        algs.ksp(to_handle(g), 0, 100, k=1, max_cost_factor=None, unique=True)


def test_batch_max_flow_invalid_pairs_shape_raises(build_graph, algs, to_handle):
    g = _small_graph(build_graph)
    pairs = np.array([[0, 1, 2]], dtype=np.int32)  # wrong shape
    with pytest.raises(TypeError):
        algs.batch_max_flow(to_handle(g), pairs)


def test_batch_max_flow_masks_length_mismatch_raises(build_graph, algs, to_handle):
    g = _small_graph(build_graph)
    pairs = np.array([[0, 1], [0, 1]], dtype=np.int32)
    node_masks = [np.array([True, True], dtype=bool)]  # len 1, need 2
    with pytest.raises(TypeError):
        algs.batch_max_flow(to_handle(g), pairs, node_masks=node_masks)


def test_flowstate_views_are_readonly():
    """Verify that FlowState views return read-only arrays."""
    g = _make_graph(
        3,
        np.array([0, 1], dtype=np.int32),
        np.array([1, 2], dtype=np.int32),
        np.array([10.0, 10.0]),
        np.ones(2, dtype=np.int64),
    )

    fs = ngc.FlowState(g)

    # Test capacity_view is read-only
    capacity_view = fs.capacity_view()
    assert not capacity_view.flags.writeable, "capacity_view should be read-only"
    with pytest.raises((ValueError, RuntimeError)):
        capacity_view[0] = 5.0

    # Test residual_view is read-only
    residual_view = fs.residual_view()
    assert not residual_view.flags.writeable, "residual_view should be read-only"
    with pytest.raises((ValueError, RuntimeError)):
        residual_view[0] = 5.0

    # Test edge_flow_view is read-only
    edge_flow_view = fs.edge_flow_view()
    assert not edge_flow_view.flags.writeable, "edge_flow_view should be read-only"
    with pytest.raises((ValueError, RuntimeError)):
        edge_flow_view[0] = 5.0


def test_flowgraph_views_are_readonly():
    """Verify that FlowGraph views return read-only arrays."""
    g = _make_graph(
        3,
        np.array([0, 1], dtype=np.int32),
        np.array([1, 2], dtype=np.int32),
        np.array([10.0, 10.0]),
        np.ones(2, dtype=np.int64),
    )

    fg = ngc.FlowGraph(g)

    # Test capacity_view is read-only
    capacity_view = fg.capacity_view()
    assert not capacity_view.flags.writeable, "capacity_view should be read-only"
    with pytest.raises((ValueError, RuntimeError)):
        capacity_view[0] = 999.0

    # Test residual_view is read-only
    residual_view = fg.residual_view()
    assert not residual_view.flags.writeable, "residual_view should be read-only"
    with pytest.raises((ValueError, RuntimeError)):
        residual_view[0] = 999.0

    # Test edge_flow_view is read-only
    edge_flow_view = fg.edge_flow_view()
    assert not edge_flow_view.flags.writeable, "edge_flow_view should be read-only"
    with pytest.raises((ValueError, RuntimeError)):
        edge_flow_view[0] = 999.0
