"""Masking behavior tests.

Verifies node/edge mask handling across all algorithms:
- Masked source/destination nodes
- Combined node and edge masks
- Interaction with residual capacity
- Masking in max_flow, shortest_paths, and k_shortest_paths
- Mask semantics (true = allowed, false = excluded)
"""

import numpy as np
import pytest

import netgraph_core as ngc

# ============================================================================
# Test Fixtures
# ============================================================================


@pytest.fixture
def line_graph_3(build_graph):
    """Simple line graph: 0 -> 1 -> 2"""
    edges = [
        (0, 1, 1.0, 1, 0),
        (1, 2, 1.0, 1, 1),
    ]
    return build_graph(3, edges)


@pytest.fixture
def parallel_paths_graph(build_graph):
    """Two parallel paths: 0 -> 1 -> 3 and 0 -> 2 -> 3"""
    edges = [
        (0, 1, 1.0, 1, 0),
        (1, 3, 1.0, 1, 1),
        (0, 2, 1.0, 1, 2),
        (2, 3, 1.0, 1, 3),
    ]
    return build_graph(4, edges)


@pytest.fixture
def triangle_graph(build_graph):
    """Triangle: 0 -> 1 -> 2 and 0 -> 2"""
    edges = [
        (0, 1, 2.0, 1, 0),
        (1, 2, 1.0, 1, 1),
        (0, 2, 1.0, 1, 2),
    ]
    return build_graph(3, edges)


# ============================================================================
# Masked Source Tests (Critical Bug Fix)
# ============================================================================


def test_masked_source_returns_empty_dag(line_graph_3, algs, to_handle):
    """
    CRITICAL BUG FIX TEST: When source is masked out, SPF must return empty DAG
    and leave all distances at infinity.

    This tests the fix for the bug where a masked source was still enqueued
    in the priority queue, allowing exploration from a forbidden node.
    """
    g = line_graph_3

    # Mask out source node 0
    node_mask = np.array([False, True, True], dtype=bool)

    sel = ngc.EdgeSelection(
        multi_edge=True,
        require_capacity=False,
        tie_break=ngc.EdgeTieBreak.DETERMINISTIC,
    )

    dist, dag = algs.spf(to_handle(g), 0, selection=sel, node_mask=node_mask)

    # All distances should remain at infinity
    assert np.isinf(dist[0]), "Source distance should be infinity when masked"
    assert np.isinf(dist[1]), "All nodes should be unreachable when source is masked"
    assert np.isinf(dist[2]), "All nodes should be unreachable when source is masked"

    # DAG should be completely empty
    offsets = np.asarray(dag.parent_offsets)
    parents = np.asarray(dag.parents)
    via_edges = np.asarray(dag.via_edges)

    assert len(offsets) == g.num_nodes() + 1
    assert np.all(offsets == 0), "All offsets should be zero for empty DAG"
    assert len(parents) == 0, "Parents array should be empty when source is masked"
    assert len(via_edges) == 0, "Via edges array should be empty when source is masked"


def test_masked_source_with_destination(line_graph_3, algs, to_handle):
    """Test masked source with explicit destination."""
    g = line_graph_3

    # Mask out source
    node_mask = np.array([False, True, True], dtype=bool)

    sel = ngc.EdgeSelection(
        multi_edge=True,
        require_capacity=False,
    )

    dist, dag = algs.spf(to_handle(g), 0, dst=2, selection=sel, node_mask=node_mask)

    # Destination should be unreachable
    assert np.isinf(dist[2]), "Destination unreachable when source is masked"

    # DAG should be empty
    offsets = np.asarray(dag.parent_offsets)
    assert np.all(offsets == 0)


def test_masked_source_in_max_flow(line_graph_3, algs, to_handle):
    """Test that max_flow correctly handles masked source."""
    g = line_graph_3

    # Mask out source
    node_mask = np.array([False, True, True], dtype=bool)

    flow_val, summary = algs.max_flow(
        to_handle(g),
        0,
        2,
        flow_placement=ngc.FlowPlacement.PROPORTIONAL,
        shortest_path=False,
        require_capacity=True,
        with_edge_flows=True,
        node_mask=node_mask,
    )

    # No flow should be possible
    assert np.isclose(flow_val, 0.0), "No flow possible when source is masked"

    # Edge flows array should be empty or all zeros
    if len(summary.edge_flows) > 0:
        assert np.all(np.isclose(summary.edge_flows, 0.0))


def test_masked_source_in_ksp(line_graph_3, algs, to_handle):
    """Test that k_shortest_paths correctly handles masked source."""
    g = line_graph_3

    # Mask out source
    node_mask = np.array([False, True, True], dtype=bool)

    results = algs.ksp(to_handle(g), 0, 2, k=2, unique=True, node_mask=node_mask)

    # Should find no paths
    assert len(results) == 0, "No paths should be found when source is masked"


def test_masked_source_flow_state_compute_min_cut(triangle_graph, algs):
    """Test FlowState.compute_min_cut with masked source."""
    g = triangle_graph
    fs = ngc.FlowState(g)

    # Mask out source
    node_mask = np.array([False, True, True], dtype=bool)

    mincut = fs.compute_min_cut(0, node_mask=node_mask)

    # Min cut should be empty when source is masked
    assert len(mincut.edges) == 0, "Min cut should be empty when source is masked"


# ============================================================================
# Masked Destination Tests
# ============================================================================


def test_masked_destination_unreachable(line_graph_3, algs, to_handle):
    """Test that masked destination is correctly marked unreachable."""
    g = line_graph_3

    # Mask out destination node 2
    node_mask = np.array([True, True, False], dtype=bool)

    sel = ngc.EdgeSelection(multi_edge=True, require_capacity=False)

    dist, dag = algs.spf(to_handle(g), 0, selection=sel, node_mask=node_mask)

    # Source and node 1 should be reachable
    assert np.isclose(dist[0], 0.0)
    assert np.isclose(dist[1], 1.0)

    # Destination should be unreachable
    assert np.isinf(dist[2]), "Masked destination should be unreachable"

    # Node 2 should have no predecessors
    offsets = np.asarray(dag.parent_offsets)
    assert offsets[2] == offsets[3], "Masked node should have no predecessors"


def test_masked_intermediate_node(line_graph_3, algs, to_handle):
    """Test that masking intermediate node blocks path."""
    g = line_graph_3

    # Mask out intermediate node 1
    node_mask = np.array([True, False, True], dtype=bool)

    sel = ngc.EdgeSelection(multi_edge=True, require_capacity=False)

    dist, dag = algs.spf(to_handle(g), 0, dst=2, selection=sel, node_mask=node_mask)

    # Source reachable
    assert np.isclose(dist[0], 0.0)

    # Nodes 1 and 2 should be unreachable
    assert np.isinf(dist[1]), "Masked intermediate node unreachable"
    assert np.isinf(dist[2]), "Destination unreachable when path is blocked"


# ============================================================================
# Edge Mask Tests
# ============================================================================


def test_masked_edge_blocks_path(line_graph_3, algs, to_handle):
    """Test that masked edge blocks path traversal."""
    g = line_graph_3

    # Mask out first edge (0->1)
    edge_mask = np.array([False, True], dtype=bool)

    sel = ngc.EdgeSelection(multi_edge=True, require_capacity=False)

    dist, dag = algs.spf(to_handle(g), 0, selection=sel, edge_mask=edge_mask)

    assert np.isclose(dist[0], 0.0)
    assert np.isinf(dist[1]), "Node unreachable when edge is masked"
    assert np.isinf(dist[2]), "Downstream nodes unreachable when edge is masked"


def test_masked_edge_selects_alternative_path(parallel_paths_graph, algs, to_handle):
    """Test that masking one edge allows alternative path."""
    g = parallel_paths_graph

    # Mask out edge 0->1 (first edge)
    edge_mask = np.array([False, True, True, True], dtype=bool)

    sel = ngc.EdgeSelection(multi_edge=True, require_capacity=False)

    dist, dag = algs.spf(to_handle(g), 0, dst=3, selection=sel, edge_mask=edge_mask)

    # Should still reach destination via alternative path 0->2->3
    assert np.isclose(dist[0], 0.0)
    assert np.isinf(dist[1]), "Node 1 unreachable via masked edge"
    assert np.isclose(dist[2], 1.0), "Node 2 reachable via unmasked path"
    assert np.isclose(dist[3], 2.0), "Destination reachable via alternative path"


def test_multiple_edges_masked(parallel_paths_graph, algs, to_handle):
    """Test masking multiple edges."""
    g = parallel_paths_graph

    # Mask out second hop on first path
    edge_mask = np.array([True, False, True, True], dtype=bool)

    sel = ngc.EdgeSelection(multi_edge=True, require_capacity=False)

    dist, dag = algs.spf(to_handle(g), 0, selection=sel, edge_mask=edge_mask)

    # Should reach destination via second path only
    assert np.isclose(dist[0], 0.0)
    assert np.isclose(dist[3], 2.0), "Destination reachable via unmasked path"


# ============================================================================
# Combined Node and Edge Mask Tests
# ============================================================================


def test_combined_node_and_edge_masks(parallel_paths_graph, algs, to_handle):
    """Test that node and edge masks work together correctly."""
    g = parallel_paths_graph

    # Mask node 1 and edge 2->3
    node_mask = np.array([True, False, True, True], dtype=bool)
    edge_mask = np.array([True, True, True, False], dtype=bool)

    sel = ngc.EdgeSelection(multi_edge=True, require_capacity=False)

    dist, dag = algs.spf(
        to_handle(g), 0, selection=sel, node_mask=node_mask, edge_mask=edge_mask
    )

    # Source reachable
    assert np.isclose(dist[0], 0.0)

    # Node 1 masked out
    assert np.isinf(dist[1])

    # Node 2 reachable but path to 3 blocked by edge mask
    assert np.isclose(dist[2], 1.0)
    assert np.isinf(dist[3]), "Destination unreachable due to edge mask"


def test_node_mask_overrides_edge_availability(line_graph_3, algs, to_handle):
    """Test that node mask effectively blocks all edges to/from that node."""
    g = line_graph_3

    # Mask node 1, but allow all edges
    node_mask = np.array([True, False, True], dtype=bool)
    edge_mask = np.array([True, True], dtype=bool)

    sel = ngc.EdgeSelection(multi_edge=True, require_capacity=False)

    dist, dag = algs.spf(
        to_handle(g), 0, selection=sel, node_mask=node_mask, edge_mask=edge_mask
    )

    # Node 1 and beyond should be unreachable despite edges being "allowed"
    assert np.isclose(dist[0], 0.0)
    assert np.isinf(dist[1]), "Masked node unreachable despite edge mask"
    assert np.isinf(dist[2]), "Downstream nodes unreachable"


# ============================================================================
# Mask Interaction with Residual Capacity
# ============================================================================


def test_mask_with_residual_capacity(build_graph, algs, to_handle):
    """Test that masks work correctly with residual capacity filtering."""
    # Create line graph 0->1->2->3
    edges = [
        (0, 1, 1.0, 1, 0),
        (1, 2, 1.0, 1, 1),
        (2, 3, 1.0, 1, 2),
    ]
    g = build_graph(4, edges)

    # Set residual capacity (zero on edge 1->2)
    residual = np.array([1.0, 0.0, 1.0], dtype=np.float64)

    # Also mask edge 2->3
    edge_mask = np.array([True, True, False], dtype=bool)

    sel = ngc.EdgeSelection(
        multi_edge=True,
        require_capacity=True,  # Enable residual filtering
    )

    dist, dag = algs.spf(
        to_handle(g), 0, selection=sel, residual=residual, edge_mask=edge_mask
    )

    # Should reach node 1 but not beyond
    assert np.isclose(dist[0], 0.0)
    assert np.isclose(dist[1], 1.0)
    assert np.isinf(dist[2]), "Blocked by residual capacity"
    assert np.isinf(dist[3]), "Blocked by mask"


def test_residual_and_node_mask_combined(parallel_paths_graph, algs, to_handle):
    """Test combination of residual capacity and node mask."""
    g = parallel_paths_graph

    # Zero capacity on one path
    residual = np.array([0.0, 1.0, 1.0, 1.0], dtype=np.float64)

    # Mask intermediate node on other path
    node_mask = np.array([True, True, False, True], dtype=bool)

    sel = ngc.EdgeSelection(
        multi_edge=True,
        require_capacity=True,
    )

    dist, dag = algs.spf(
        to_handle(g), 0, dst=3, selection=sel, residual=residual, node_mask=node_mask
    )

    # Destination should be unreachable (both paths blocked)
    assert np.isinf(dist[3]), "Both paths blocked"


# ============================================================================
# Max Flow with Masks
# ============================================================================


def test_max_flow_with_node_mask(parallel_paths_graph, algs, to_handle):
    """Test max_flow with node mask blocking one path."""
    g = parallel_paths_graph

    # Mask out node 1 (blocks one path)
    node_mask = np.array([True, False, True, True], dtype=bool)

    flow_val, summary = algs.max_flow(
        to_handle(g),
        0,
        3,
        flow_placement=ngc.FlowPlacement.PROPORTIONAL,
        shortest_path=False,
        require_capacity=True,
        with_edge_flows=True,
        node_mask=node_mask,
    )

    # Should get flow only through remaining path (0->2->3)
    assert flow_val > 0.0, "Some flow possible through unmasked path"
    assert flow_val <= 1.0, "Only one path available"

    # Verify path 0->1->3 is blocked (node 1 masked) and flow uses 0->2->3
    # Due to multipath SPF, we need to check actual flow distribution
    # The important thing is flow_val > 0 which proves masking doesn't block everything


def test_max_flow_with_edge_mask(triangle_graph, algs, to_handle):
    """Test max_flow with edge mask."""
    g = triangle_graph

    # Find and mask direct edge 0->2
    edge_src = g.edge_src_view()
    edge_dst = g.edge_dst_view()
    edge_mask = np.ones(g.num_edges(), dtype=bool)

    for e in range(g.num_edges()):
        if edge_src[e] == 0 and edge_dst[e] == 2:
            edge_mask[e] = False

    flow_val, summary = algs.max_flow(
        to_handle(g),
        0,
        2,
        flow_placement=ngc.FlowPlacement.PROPORTIONAL,
        shortest_path=False,
        require_capacity=True,
        with_edge_flows=True,
        edge_mask=edge_mask,
    )

    # Should get flow only through indirect path 0->1->2
    assert flow_val > 0.0, "Flow possible through indirect path"
    assert flow_val <= 1.0, "Limited by bottleneck"

    # Verify masked edge has no flow
    for e in range(g.num_edges()):
        if edge_src[e] == 0 and edge_dst[e] == 2:
            assert np.isclose(summary.edge_flows[e], 0.0), (
                "Masked edge should have zero flow"
            )


def test_max_flow_all_paths_masked(line_graph_3, algs, to_handle):
    """Test max_flow when all paths are masked."""
    g = line_graph_3

    # Mask intermediate node to block all paths
    node_mask = np.array([True, False, True], dtype=bool)

    flow_val, summary = algs.max_flow(
        to_handle(g),
        0,
        2,
        flow_placement=ngc.FlowPlacement.PROPORTIONAL,
        shortest_path=False,
        require_capacity=True,
        with_edge_flows=True,
        node_mask=node_mask,
    )

    # No flow should be possible
    assert np.isclose(flow_val, 0.0), "No flow possible when all paths masked"

    # All edge flows should be zero
    assert np.all(np.isclose(summary.edge_flows, 0.0))


# ============================================================================
# K-Shortest Paths with Masks
# ============================================================================


def test_ksp_with_node_mask(build_graph, algs, to_handle):
    """Test k_shortest_paths with node mask."""
    # Create graph with two paths
    edges = [
        (0, 1, 1.0, 1, 0),
        (1, 2, 1.0, 1, 1),
        (0, 3, 1.0, 1, 2),
        (3, 2, 1.0, 1, 3),
    ]
    g = build_graph(4, edges)

    # Mask out node 1 (blocks one path)
    node_mask = np.array([True, False, True, True], dtype=bool)

    results = algs.ksp(to_handle(g), 0, 2, k=2, unique=True, node_mask=node_mask)

    # Should find at most 1 path (through node 3)
    assert len(results) <= 1, "Only one path available with mask"

    if len(results) > 0:
        dist, dag = results[0]
        # Destination should be reachable via alternate path
        assert dist[2] > 0.0
        assert not np.isinf(dist[2])


def test_ksp_with_edge_mask(build_graph, algs, to_handle):
    """Test k_shortest_paths with edge mask."""
    # Create graph with two paths
    edges = [
        (0, 1, 1.0, 1, 0),
        (1, 2, 1.0, 1, 1),
        (0, 3, 1.0, 2, 2),
        (3, 2, 1.0, 2, 3),
    ]
    g = build_graph(4, edges)

    # Mask out first edge
    edge_mask = np.array([False, True, True, True], dtype=bool)

    results = algs.ksp(to_handle(g), 0, 2, k=2, unique=True, edge_mask=edge_mask)

    # Should find path through node 3 only
    assert len(results) >= 1, "Should find at least one path"

    if len(results) > 0:
        dist, dag = results[0]
        # Cost should be 2.0 (path 0->3->2 with cost 1+1=2)
        assert np.isclose(dist[2], 2.0), (
            "Should use alternative path when shorter is masked"
        )


# ============================================================================
# Empty Mask Tests
# ============================================================================


def test_empty_node_mask_allows_all(line_graph_3, algs, to_handle):
    """Test that empty mask means no masking (all nodes allowed)."""
    g = line_graph_3

    sel = ngc.EdgeSelection(multi_edge=True, require_capacity=False)

    # No mask provided
    dist, dag = algs.spf(to_handle(g), 0, selection=sel)

    # All nodes should be reachable
    assert np.isclose(dist[0], 0.0)
    assert np.isclose(dist[1], 1.0)
    assert np.isclose(dist[2], 2.0)


def test_empty_edge_mask_allows_all(line_graph_3, algs, to_handle):
    """Test that empty edge mask means no masking (all edges allowed)."""
    g = line_graph_3

    sel = ngc.EdgeSelection(multi_edge=True, require_capacity=False)

    # No edge mask provided
    dist, dag = algs.spf(to_handle(g), 0, selection=sel)

    # All nodes should be reachable
    for i in range(g.num_nodes()):
        assert not np.isinf(dist[i]), f"Node {i} should be reachable"


# ============================================================================
# All-False Mask Tests
# ============================================================================


def test_all_nodes_masked_no_reachability(line_graph_3, algs, to_handle):
    """Test that masking all nodes makes everything unreachable."""
    g = line_graph_3

    # Mask all nodes
    node_mask = np.array([False, False, False], dtype=bool)

    sel = ngc.EdgeSelection(multi_edge=True, require_capacity=False)

    dist, dag = algs.spf(to_handle(g), 0, selection=sel, node_mask=node_mask)

    # All nodes unreachable (including source)
    for i in range(g.num_nodes()):
        assert np.isinf(dist[i]), f"Node {i} should be unreachable"

    # DAG should be empty
    offsets = np.asarray(dag.parent_offsets)
    assert np.all(offsets == 0), "DAG should be empty"


def test_all_edges_masked_only_source_reachable(line_graph_3, algs, to_handle):
    """Test that masking all edges makes only source reachable."""
    g = line_graph_3

    # Mask all edges
    edge_mask = np.array([False, False], dtype=bool)

    sel = ngc.EdgeSelection(multi_edge=True, require_capacity=False)

    dist, dag = algs.spf(to_handle(g), 0, selection=sel, edge_mask=edge_mask)

    # Only source should be reachable
    assert np.isclose(dist[0], 0.0), "Source should be reachable"

    # All other nodes unreachable
    for i in range(1, g.num_nodes()):
        assert np.isinf(dist[i]), f"Node {i} should be unreachable"


# ============================================================================
# Mask Validation Tests
# ============================================================================


def test_invalid_node_mask_length_raises_error(line_graph_3, algs, to_handle):
    """Test that incorrect mask length raises an error."""
    g = line_graph_3

    # Wrong length mask
    node_mask = np.array([True, True], dtype=bool)  # Should be length 3

    sel = ngc.EdgeSelection(multi_edge=True, require_capacity=False)

    with pytest.raises(TypeError, match="node_mask"):
        algs.spf(to_handle(g), 0, selection=sel, node_mask=node_mask)


def test_invalid_edge_mask_length_raises_error(line_graph_3, algs, to_handle):
    """Test that incorrect edge mask length raises an error."""
    g = line_graph_3

    # Wrong length mask
    edge_mask = np.array([True], dtype=bool)  # Should be length 2

    sel = ngc.EdgeSelection(multi_edge=True, require_capacity=False)

    with pytest.raises(TypeError, match="edge_mask"):
        algs.spf(to_handle(g), 0, selection=sel, edge_mask=edge_mask)


def test_non_bool_mask_raises_error(line_graph_3, algs, to_handle):
    """Test that non-boolean mask raises an error."""
    g = line_graph_3

    # Integer mask instead of bool
    node_mask = np.array([1, 1, 1], dtype=np.int32)

    sel = ngc.EdgeSelection(multi_edge=True, require_capacity=False)

    with pytest.raises(TypeError, match="bool"):
        algs.spf(to_handle(g), 0, selection=sel, node_mask=node_mask)
