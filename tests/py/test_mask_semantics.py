"""Quick CI probes for critical masking semantics.

Validates that masks behave correctly end-to-end across all APIs,
as suggested in the final masking review.
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
def tiny_graph():
    """Tiny graph: 0 -> 1 -> 2"""
    return _make_graph(
        3,
        np.array([0, 1], np.int32),
        np.array([1, 2], np.int32),
        np.array([1.0, 1.0], np.float64),
        np.array([1, 1], np.int64),
    )


def test_spf_masked_source_returns_empty_dag_inf_distances(tiny_graph):
    """SPF: masked source -> empty DAG + inf distances."""
    algs = ngc.Algorithms(ngc.Backend.cpu())
    G = algs.build_graph(tiny_graph)

    dist, dag = algs.spf(G, 0, dst=2, node_mask=np.array([False, True, True], bool))

    assert np.isinf(dist[2]), "Distance to dst should be infinity"
    assert dag.parent_offsets[-1] == 0, "DAG should be empty"


def test_spf_edge_mask_excludes_path(tiny_graph):
    """SPF: edge mask excludes (1->2) -> unreachable."""
    algs = ngc.Algorithms(ngc.Backend.cpu())
    G = algs.build_graph(tiny_graph)

    # Mask out edge 1 (1->2)
    dist, dag = algs.spf(G, 0, dst=2, edge_mask=np.array([True, False], bool))

    assert np.isinf(dist[2]), "Node 2 should be unreachable"


def test_max_flow_destination_masked_returns_zero(tiny_graph):
    """MaxFlow: destination masked -> zero flow, empty cut."""
    algs = ngc.Algorithms(ngc.Backend.cpu())
    G = algs.build_graph(tiny_graph)

    # Mask out destination (node 2)
    total, summary = algs.max_flow(
        G, 0, 2, node_mask=np.array([True, True, False], bool)
    )

    assert total == 0.0, "Flow should be zero"
    assert len(summary.min_cut.edges) == 0, "Min-cut should be empty"


def test_flow_state_place_max_flow_accepts_masks(tiny_graph):
    """FlowState.place_max_flow now accepts masks (strict validation)."""
    fs = ngc.FlowState(tiny_graph)

    # Mask out source
    placed = fs.place_max_flow(
        0,
        2,
        ngc.FlowPlacement.PROPORTIONAL,
        node_mask=np.array([False, True, True], bool),
    )

    assert placed == 0.0, "No flow should be placed with masked source"


def test_flow_state_place_max_flow_with_valid_masks(tiny_graph):
    """FlowState.place_max_flow works correctly with valid masks."""
    fs = ngc.FlowState(tiny_graph)

    # Allow all nodes
    placed = fs.place_max_flow(
        0,
        2,
        ngc.FlowPlacement.PROPORTIONAL,
        node_mask=np.array([True, True, True], bool),
        edge_mask=np.array([True, True], bool),
    )

    assert placed > 0, "Flow should be placed with all nodes/edges allowed"


def test_batch_max_flow_validates_mask_lengths(tiny_graph):
    """Batch: shape mismatches throw (node_masks length != num_nodes)."""
    algs = ngc.Algorithms(ngc.Backend.cpu())
    G = algs.build_graph(tiny_graph)

    with pytest.raises(TypeError, match="length"):
        # Wrong node mask length (2 instead of 3)
        algs.batch_max_flow(
            G,
            np.array([[0, 2]], np.int32),
            node_masks=[np.ones(2, bool)],  # Wrong length!
        )


def test_flow_policy_with_masks_end_to_end(tiny_graph):
    """FlowPolicy with masks: end-to-end place_demand."""
    algs = ngc.Algorithms(ngc.Backend.cpu())
    G = algs.build_graph(tiny_graph)

    # Create policy with destination masked out
    cfg = ngc.FlowPolicyConfig()
    policy = ngc.FlowPolicy(algs, G, cfg, node_mask=np.array([True, True, False], bool))

    fg = ngc.FlowGraph(tiny_graph)
    placed, leftover = policy.place_demand(fg, 0, 2, flowClass=0, volume=1.0)

    assert placed == 0.0, "No flow placed with masked destination"
    assert leftover == pytest.approx(1.0, abs=1e-9), "All volume should remain"


def test_ksp_honors_node_mask(tiny_graph):
    """KSP: honors node masks (masked intermediate node)."""
    algs = ngc.Algorithms(ngc.Backend.cpu())
    G = algs.build_graph(tiny_graph)

    # Mask out intermediate node 1
    paths = algs.ksp(G, 0, 2, k=2, node_mask=np.array([True, False, True], bool))

    assert len(paths) == 0, "No paths should exist with node 1 masked"


def test_min_cut_with_masked_source(tiny_graph):
    """FlowState.compute_min_cut with masked source returns empty cut."""
    fs = ngc.FlowState(tiny_graph)

    # Place some flow first
    fs.place_max_flow(0, 2, ngc.FlowPlacement.PROPORTIONAL)

    # Compute min-cut with source masked
    cut = fs.compute_min_cut(0, node_mask=np.array([False, True, True], bool))

    assert len(cut.edges) == 0, "Min-cut should be empty with masked source"


def test_all_apis_accept_none_masks(tiny_graph):
    """All APIs accept None masks (no masking)."""
    algs = ngc.Algorithms(ngc.Backend.cpu())
    G = algs.build_graph(tiny_graph)

    # SPF
    dist, dag = algs.spf(G, 0, 2, node_mask=None, edge_mask=None)
    assert not np.isinf(dist[2])

    # KSP
    paths = algs.ksp(G, 0, 2, k=1, node_mask=None, edge_mask=None)
    assert len(paths) == 1

    # Max flow
    total, _ = algs.max_flow(G, 0, 2, node_mask=None, edge_mask=None)
    assert total > 0

    # FlowState.place_max_flow
    fs = ngc.FlowState(tiny_graph)
    placed = fs.place_max_flow(
        0, 2, ngc.FlowPlacement.PROPORTIONAL, node_mask=None, edge_mask=None
    )
    assert placed > 0

    # FlowState.compute_min_cut
    cut = fs.compute_min_cut(0, node_mask=None, edge_mask=None)
    assert isinstance(cut, ngc.MinCut)
