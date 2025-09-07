"""API validation: argument shapes, bounds, and error handling for frontdoor APIs."""

from __future__ import annotations

import numpy as np
import pytest

import netgraph_core as ngc


def _small_graph(build_graph):
    # 0 -> 1 with cost=1, cap=1
    return build_graph(2, [(0, 1, 1.0, 1.0, 0)])


def test_spf_invalid_node_mask_shape_raises(build_graph):
    g = _small_graph(build_graph)
    sel = ngc.EdgeSelection(
        multipath=True, require_capacity=False, tie_break=ngc.EdgeTieBreak.DETERMINISTIC
    )
    bad_mask = np.array([True], dtype=bool)  # wrong length
    with pytest.raises(TypeError):
        ngc.spf(g, 0, 1, selection=sel, node_mask=bad_mask)


def test_spf_invalid_edge_mask_shape_raises(build_graph):
    g = _small_graph(build_graph)
    sel = ngc.EdgeSelection(
        multipath=True, require_capacity=False, tie_break=ngc.EdgeTieBreak.DETERMINISTIC
    )
    bad_edge_mask = np.array([True, False, True], dtype=bool)  # wrong length
    with pytest.raises(TypeError):
        ngc.spf(g, 0, 1, selection=sel, edge_mask=bad_edge_mask)


def test_spf_invalid_source_node_raises(build_graph):
    g = _small_graph(build_graph)
    with pytest.raises(ValueError):
        ngc.spf(g, 10)  # out of range src


def test_max_flow_invalid_masks_shape_raises(build_graph):
    g = _small_graph(build_graph)
    bad_edge_mask = np.array([True, False, True], dtype=bool)
    with pytest.raises(TypeError):
        ngc.max_flow(g, 0, 1, edge_mask=bad_edge_mask)

    bad_node_mask = np.array([True], dtype=bool)
    with pytest.raises(TypeError):
        ngc.max_flow(g, 0, 1, node_mask=bad_node_mask)


def test_ksp_invalid_k_and_dst_raises(build_graph):
    g = _small_graph(build_graph)
    with pytest.raises(ValueError):
        ngc.ksp(g, 0, 1, k=0, max_cost_factor=None, unique=True)
    with pytest.raises(ValueError):
        ngc.ksp(g, 0, 100, k=1, max_cost_factor=None, unique=True)


def test_batch_max_flow_invalid_pairs_shape_raises(build_graph):
    g = _small_graph(build_graph)
    pairs = np.array([[0, 1, 2]], dtype=np.int32)  # wrong shape
    with pytest.raises(TypeError):
        ngc.batch_max_flow(g, pairs)


def test_batch_max_flow_masks_length_mismatch_raises(build_graph):
    g = _small_graph(build_graph)
    pairs = np.array([[0, 1], [0, 1]], dtype=np.int32)
    node_masks = [np.array([True, True], dtype=bool)]  # len 1, need 2
    with pytest.raises(TypeError):
        ngc.batch_max_flow(g, pairs, node_masks=node_masks)
