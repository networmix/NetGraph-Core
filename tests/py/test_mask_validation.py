"""Mask validation tests.

Validates that API entry points enforce mask length requirements
and reject malformed masks with clear error messages.
"""

import numpy as np
import pytest

import netgraph_core as ngc


def _make_graph(num_nodes, src, dst, cap, cost):
    """Helper to build graph with auto-generated ext_edge_ids."""
    ext_edge_ids = np.arange(len(src), dtype=np.int64)
    return ngc.StrictMultiDiGraph.from_arrays(
        num_nodes, src, dst, cap, cost, ext_edge_ids
    )


@pytest.fixture
def simple_graph():
    """Create a simple 3-node graph: 0->1->2."""
    return _make_graph(
        3,
        np.array([0, 1], np.int32),
        np.array([1, 2], np.int32),
        np.array([10.0, 10.0], np.float64),
        np.array([1, 1], np.int64),
    )


def test_spf_rejects_wrong_node_mask_length(simple_graph):
    """shortest_paths rejects node_mask with wrong length."""
    algs = ngc.Algorithms(ngc.Backend.cpu())
    gh = algs.build_graph(simple_graph)

    with pytest.raises(TypeError, match="node_mask length must equal"):
        algs.spf(gh, 0, 2, node_mask=np.array([True, True], dtype=bool))


def test_spf_rejects_wrong_edge_mask_length(simple_graph):
    """shortest_paths rejects edge_mask with wrong length."""
    algs = ngc.Algorithms(ngc.Backend.cpu())
    gh = algs.build_graph(simple_graph)

    with pytest.raises(TypeError, match="edge_mask length must equal"):
        algs.spf(gh, 0, 2, edge_mask=np.array([True, True, True], dtype=bool))


def test_ksp_rejects_wrong_node_mask_length(simple_graph):
    """k_shortest_paths rejects node_mask with wrong length."""
    algs = ngc.Algorithms(ngc.Backend.cpu())
    gh = algs.build_graph(simple_graph)

    with pytest.raises(TypeError, match="node_mask length must equal"):
        algs.ksp(gh, 0, 2, k=2, node_mask=np.array([True], dtype=bool))


def test_ksp_rejects_wrong_edge_mask_length(simple_graph):
    """k_shortest_paths rejects edge_mask with wrong length."""
    algs = ngc.Algorithms(ngc.Backend.cpu())
    gh = algs.build_graph(simple_graph)

    with pytest.raises(TypeError, match="edge_mask length must equal"):
        algs.ksp(gh, 0, 2, k=2, edge_mask=np.array([True], dtype=bool))


def test_max_flow_rejects_wrong_node_mask_length(simple_graph):
    """calc_max_flow rejects node_mask with wrong length."""
    algs = ngc.Algorithms(ngc.Backend.cpu())
    gh = algs.build_graph(simple_graph)

    with pytest.raises(TypeError, match="node_mask length must equal"):
        algs.max_flow(
            gh, 0, 2, node_mask=np.array([True, True, True, True], dtype=bool)
        )


def test_max_flow_rejects_wrong_edge_mask_length(simple_graph):
    """calc_max_flow rejects edge_mask with wrong length."""
    algs = ngc.Algorithms(ngc.Backend.cpu())
    gh = algs.build_graph(simple_graph)

    with pytest.raises(TypeError, match="edge_mask length must equal"):
        algs.max_flow(gh, 0, 2, edge_mask=np.array([True], dtype=bool))


def test_batch_max_flow_rejects_wrong_mask_length(simple_graph):
    """batch_max_flow validates all masks upfront."""
    algs = ngc.Algorithms(ngc.Backend.cpu())
    gh = algs.build_graph(simple_graph)

    # Wrong node mask length in batch (batch_max_flow expects numpy array of pairs)
    pairs = np.array([[0, 2], [1, 2]], dtype=np.int32)
    with pytest.raises(TypeError, match="node_masks length must equal"):
        algs.batch_max_flow(
            gh,
            pairs,
            node_masks=[
                np.array([True, True, True], dtype=bool),
                np.array([True], dtype=bool),
            ],
        )

    # Wrong edge mask length in batch
    pairs = np.array([[0, 2]], dtype=np.int32)
    with pytest.raises(TypeError, match="edge_masks length must equal"):
        algs.batch_max_flow(
            gh, pairs, edge_masks=[np.array([True, True, True], dtype=bool)]
        )


def test_flow_state_compute_min_cut_rejects_wrong_mask_length(simple_graph):
    """FlowState.compute_min_cut rejects wrong mask lengths."""
    fs = ngc.FlowState(simple_graph)

    with pytest.raises(TypeError, match="node_mask length must equal"):
        fs.compute_min_cut(0, node_mask=np.array([True], dtype=bool))

    with pytest.raises(TypeError, match="edge_mask length must equal"):
        fs.compute_min_cut(0, edge_mask=np.array([True, True, True], dtype=bool))


def test_flow_state_place_max_flow_rejects_wrong_mask_length(simple_graph):
    """FlowState.place_max_flow rejects wrong mask lengths."""
    fs = ngc.FlowState(simple_graph)

    with pytest.raises(TypeError, match="node_mask length must equal"):
        fs.place_max_flow(
            0,
            2,
            ngc.FlowPlacement.PROPORTIONAL,
            node_mask=np.array([True, True], dtype=bool),
        )

    with pytest.raises(TypeError, match="edge_mask length must equal"):
        fs.place_max_flow(
            0,
            2,
            ngc.FlowPlacement.PROPORTIONAL,
            edge_mask=np.array([True, True, True], dtype=bool),
        )


def test_flow_policy_rejects_wrong_mask_length_at_construction(simple_graph):
    """FlowPolicy validates mask lengths at construction time."""
    algs = ngc.Algorithms(ngc.Backend.cpu())
    gh = algs.build_graph(simple_graph)
    cfg = ngc.FlowPolicyConfig()

    with pytest.raises(TypeError, match="node_mask length must equal"):
        ngc.FlowPolicy(algs, gh, cfg, node_mask=np.array([True, True], dtype=bool))

    with pytest.raises(TypeError, match="edge_mask length must equal"):
        ngc.FlowPolicy(algs, gh, cfg, edge_mask=np.array([True], dtype=bool))


def test_flow_state_place_max_flow_honors_node_mask(simple_graph):
    """FlowState.place_max_flow respects node_mask."""
    fs = ngc.FlowState(simple_graph)

    # Mask out source
    node_mask = np.array([False, True, True], dtype=bool)
    total = fs.place_max_flow(0, 2, ngc.FlowPlacement.PROPORTIONAL, node_mask=node_mask)
    assert total == 0.0


def test_flow_state_place_max_flow_honors_edge_mask(simple_graph):
    """FlowState.place_max_flow respects edge_mask."""
    fs = ngc.FlowState(simple_graph)

    # Mask out first edge (0->1)
    edge_mask = np.array([False, True], dtype=bool)
    total = fs.place_max_flow(0, 2, ngc.FlowPlacement.PROPORTIONAL, edge_mask=edge_mask)
    assert total == 0.0


def test_flow_state_place_max_flow_with_valid_masks(simple_graph):
    """FlowState.place_max_flow works correctly with valid masks."""
    fs = ngc.FlowState(simple_graph)

    # Allow all nodes and edges
    node_mask = np.array([True, True, True], dtype=bool)
    edge_mask = np.array([True, True], dtype=bool)
    total = fs.place_max_flow(
        0, 2, ngc.FlowPlacement.PROPORTIONAL, node_mask=node_mask, edge_mask=edge_mask
    )
    assert total == pytest.approx(10.0, abs=1e-9)


def test_all_apis_accept_empty_masks(simple_graph):
    """All APIs accept empty masks (no masking)."""
    algs = ngc.Algorithms(ngc.Backend.cpu())
    gh = algs.build_graph(simple_graph)

    # SPF with no masks
    dist, dag = algs.spf(gh, 0, 2)
    assert not np.isinf(dist[2])

    # KSP with no masks
    paths = algs.ksp(gh, 0, 2, k=1)
    assert len(paths) == 1

    # Max flow with no masks
    total, _ = algs.max_flow(gh, 0, 2)
    assert total > 0

    # Batch max flow with no masks (expects numpy array)
    pairs = np.array([[0, 2]], dtype=np.int32)
    summaries = algs.batch_max_flow(gh, pairs)
    assert len(summaries) == 1

    # FlowState with no masks
    fs = ngc.FlowState(simple_graph)
    cut = fs.compute_min_cut(0)
    assert isinstance(cut, ngc.MinCut)

    total = fs.place_max_flow(0, 2, ngc.FlowPlacement.PROPORTIONAL)
    assert total > 0


def test_mask_validation_happens_before_expensive_operations(simple_graph):
    """Mask validation fails fast before expensive graph operations."""
    algs = ngc.Algorithms(ngc.Backend.cpu())
    gh = algs.build_graph(simple_graph)

    # Invalid mask should be caught before Dijkstra runs
    with pytest.raises(TypeError, match="node_mask length must equal"):
        algs.spf(gh, 0, 2, node_mask=np.array([True], dtype=bool))
