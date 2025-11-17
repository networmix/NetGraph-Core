"""FlowPolicy path distribution tests.

Verifies that LSPs distribute correctly across available equal-cost paths:
- All available paths are utilized
- Flow distribution is balanced (not clustered)
- Per-path volumes match expected values
- ECMP equal-balancing constraint is satisfied
"""

from __future__ import annotations

import pytest

import netgraph_core as ngc

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================


def analyze_path_usage(
    fg: ngc.FlowGraph, graph: ngc.StrictMultiDiGraph, source: int, target: int
) -> dict[int, float]:
    """Analyze which parallel paths are used and their volumes.

    Returns dict mapping middle_node_id -> total_flow_volume
    """
    edge_flows = fg.edge_flow_view()
    edge_src = graph.edge_src_view()
    edge_dst = graph.edge_dst_view()

    path_usage = {}

    # Find all edges from source
    for edge_id in range(graph.num_edges()):
        if edge_src[edge_id] == source:
            middle_node = edge_dst[edge_id]
            flow = edge_flows[edge_id]
            if flow > 0.001:
                path_usage[middle_node] = flow

    return path_usage


def count_lsps_per_path(
    policy: ngc.FlowPolicy,
    fg: ngc.FlowGraph,
    graph: ngc.StrictMultiDiGraph,
    source: int,
    target: int,
) -> dict[int, list[int]]:
    """Count how many LSPs use each path.

    Returns dict mapping middle_node_id -> list of flow_ids
    """
    edge_src = graph.edge_src_view()
    edge_dst = graph.edge_dst_view()

    path_lsps = {}

    # For each flow, determine which path it uses
    for flow_id in range(policy.flow_count()):
        try:
            flow_index = ngc.FlowIndex(source, target, 0, flow_id)
            flow_edges = fg.get_flow_edges(flow_index)

            # Find the middle node this flow uses
            # flow_edges is list of (edge_id, volume) tuples
            for edge_id, _ in flow_edges:
                src = edge_src[edge_id]
                dst = edge_dst[edge_id]
                if src == source:
                    middle_node = dst
                    if middle_node not in path_lsps:
                        path_lsps[middle_node] = []
                    path_lsps[middle_node].append(flow_id)
                    break
        except (KeyError, RuntimeError):
            pass  # Flow might not exist or have no edges

    return path_lsps


# ============================================================================
# TEST FIXTURES: Parallel Path Topologies
# ============================================================================


@pytest.fixture
def parallel_3_paths(build_graph):
    """Three parallel 2-hop paths, each with capacity 100."""
    edges = [
        (0, 1, 1, 100.0),  # S -> M1
        (0, 2, 1, 100.0),  # S -> M2
        (0, 3, 1, 100.0),  # S -> M3
        (1, 4, 1, 100.0),  # M1 -> T
        (2, 4, 1, 100.0),  # M2 -> T
        (3, 4, 1, 100.0),  # M3 -> T
    ]
    return build_graph(5, edges)


@pytest.fixture
def parallel_5_paths(build_graph):
    """Five parallel 2-hop paths."""
    edges = [
        (0, 1, 1, 100.0),
        (0, 2, 1, 100.0),
        (0, 3, 1, 100.0),
        (0, 4, 1, 100.0),
        (0, 5, 1, 100.0),
        (1, 6, 1, 100.0),
        (2, 6, 1, 100.0),
        (3, 6, 1, 100.0),
        (4, 6, 1, 100.0),
        (5, 6, 1, 100.0),
    ]
    return build_graph(7, edges)


@pytest.fixture
def parallel_2_paths(build_graph):
    """Two parallel 2-hop paths."""
    edges = [
        (0, 1, 1, 100.0),
        (0, 2, 1, 100.0),
        (1, 3, 1, 100.0),
        (2, 3, 1, 100.0),
    ]
    return build_graph(4, edges)


# ============================================================================
# PATH UTILIZATION TESTS
# ============================================================================


def test_all_paths_used_with_multipath_false_3_paths_6_lsps(
    parallel_3_paths, algs, to_handle
):
    """CRITICAL: When multipath=false, ALL available paths must be used.

    This is the test that would have caught the bug!

    With 3 paths and 6 LSPs:
    - Expected: 2 LSPs per path (optimal distribution)
    - Bug behavior: Only 2 paths used, 3 LSPs each, path-3 idle
    """
    g = parallel_3_paths
    fg = ngc.FlowGraph(g)

    config = ngc.FlowPolicyConfig()
    config.path_alg = ngc.PathAlg.SPF
    config.flow_placement = ngc.FlowPlacement.EQUAL_BALANCED
    config.selection = ngc.EdgeSelection(
        multi_edge=False,
        require_capacity=True,
        tie_break=ngc.EdgeTieBreak.PREFER_HIGHER_RESIDUAL,
    )
    config.multipath = False  # Each LSP uses ONE path
    config.min_flow_count = 6
    config.max_flow_count = 6

    policy = ngc.FlowPolicy(algs, to_handle(g), config)

    demand = 250.0
    placed, remaining = policy.place_demand(fg, 0, 4, flowClass=0, volume=demand)

    # Check path usage
    path_usage = analyze_path_usage(fg, g, 0, 4)

    # CRITICAL ASSERTION: All 3 paths must be used!
    assert len(path_usage) == 3, (
        f"Expected all 3 paths to be used, but only {len(path_usage)} were used. "
        f"Paths with flow: {list(path_usage.keys())}"
    )

    # All paths should have equal flow (100 each with optimal distribution)
    # With 6 LSPs and 3 paths: 2 LSPs per path = 100 capacity / 2 = 50 per LSP
    # Total per path: 2 * 50 = 100 units (if using reoptimization)
    # OR if balanced: each LSP gets demand/6, so 250/6 = 41.67 per LSP
    # With 2 LSPs per path: 2 * 41.67 = 83.33 per path

    for path_id, flow in path_usage.items():
        assert flow > 0.0, f"Path {path_id} has zero flow"

    # Check LSP distribution
    lsps_per_path = count_lsps_per_path(policy, fg, g, 0, 4)

    # With 6 LSPs and 3 paths, each path should have 2 LSPs
    for path_id in [1, 2, 3]:  # Middle nodes
        lsp_count = len(lsps_per_path.get(path_id, []))
        assert lsp_count == 2, (
            f"Path via node {path_id} should have 2 LSPs, but has {lsp_count}. "
            f"Distribution: {lsps_per_path}"
        )

    # Expected capacity: 6 LSPs * 50 units = 300 units (if perfectly balanced)
    # With demand=250, we should place it all
    assert placed >= 240.0, (  # Allow some tolerance
        f"Expected to place ~250 units with optimal distribution, got {placed}"
    )


def test_all_paths_used_with_multipath_false_5_paths_3_lsps(
    parallel_5_paths, algs, to_handle
):
    """When paths > LSPs, only use as many paths as needed, but distribute optimally."""
    g = parallel_5_paths
    fg = ngc.FlowGraph(g)

    config = ngc.FlowPolicyConfig()
    config.path_alg = ngc.PathAlg.SPF
    config.flow_placement = ngc.FlowPlacement.EQUAL_BALANCED
    config.selection = ngc.EdgeSelection(
        multi_edge=False,
        require_capacity=True,
        tie_break=ngc.EdgeTieBreak.PREFER_HIGHER_RESIDUAL,
    )
    config.multipath = False
    config.min_flow_count = 3
    config.max_flow_count = 3

    policy = ngc.FlowPolicy(algs, to_handle(g), config)

    demand = 150.0
    placed, remaining = policy.place_demand(fg, 0, 6, flowClass=0, volume=demand)

    # Check path usage
    path_usage = analyze_path_usage(fg, g, 0, 6)

    # With 3 LSPs and 5 paths: should use 3 paths (1 LSP each)
    assert len(path_usage) == 3, (
        f"Expected 3 paths to be used (1 LSP each), but {len(path_usage)} were used"
    )

    # Check LSP distribution
    lsps_per_path = count_lsps_per_path(policy, fg, g, 0, 6)

    # Each used path should have exactly 1 LSP
    for path_id, lsp_list in lsps_per_path.items():
        assert len(lsp_list) == 1, (
            f"Path {path_id} should have 1 LSP, has {len(lsp_list)}"
        )

    # With 1 LSP per path and 100 capacity: 3 * 100 = 300 available
    # Demand is 150, so should place it all
    assert placed == pytest.approx(150.0, abs=1.0), (
        f"Expected to place all 150 units, got {placed}"
    )


def test_all_paths_used_with_multipath_false_2_paths_6_lsps(
    parallel_2_paths, algs, to_handle
):
    """When LSPs > paths, distribute evenly across available paths."""
    g = parallel_2_paths
    fg = ngc.FlowGraph(g)

    config = ngc.FlowPolicyConfig()
    config.path_alg = ngc.PathAlg.SPF
    config.flow_placement = ngc.FlowPlacement.EQUAL_BALANCED
    config.selection = ngc.EdgeSelection(
        multi_edge=False,
        require_capacity=True,
        tie_break=ngc.EdgeTieBreak.PREFER_HIGHER_RESIDUAL,
    )
    config.multipath = False
    config.min_flow_count = 6
    config.max_flow_count = 6

    policy = ngc.FlowPolicy(algs, to_handle(g), config)

    demand = 200.0
    placed, remaining = policy.place_demand(fg, 0, 3, flowClass=0, volume=demand)

    # Check path usage
    path_usage = analyze_path_usage(fg, g, 0, 3)

    # Both paths must be used
    assert len(path_usage) == 2, (
        f"Expected both paths to be used, but only {len(path_usage)} were used"
    )

    # Check LSP distribution
    lsps_per_path = count_lsps_per_path(policy, fg, g, 0, 3)

    # With 6 LSPs and 2 paths: each path should have 3 LSPs
    for path_id in [1, 2]:  # Middle nodes
        lsp_count = len(lsps_per_path.get(path_id, []))
        assert lsp_count == 3, (
            f"Path via node {path_id} should have 3 LSPs, but has {lsp_count}"
        )

    # With 3 LSPs per path and 100 capacity: 100/3 = 33.33 per LSP
    # Total: 6 * 33.33 = 200 units
    assert placed == pytest.approx(200.0, abs=1.0), (
        f"Expected to place all 200 units, got {placed}"
    )


# ============================================================================
# TEST: Path Distribution Balance
# ============================================================================


def test_balanced_path_distribution_multipath_false(parallel_3_paths, algs, to_handle):
    """Verify LSPs are evenly distributed across paths, not clustered."""
    g = parallel_3_paths
    fg = ngc.FlowGraph(g)

    config = ngc.FlowPolicyConfig()
    config.path_alg = ngc.PathAlg.SPF
    config.flow_placement = ngc.FlowPlacement.EQUAL_BALANCED
    config.selection = ngc.EdgeSelection(
        multi_edge=False,
        require_capacity=True,
        tie_break=ngc.EdgeTieBreak.DETERMINISTIC,
    )
    config.multipath = False
    config.min_flow_count = 9  # 9 LSPs across 3 paths = 3 per path
    config.max_flow_count = 9

    policy = ngc.FlowPolicy(algs, to_handle(g), config)

    demand = 300.0
    placed, remaining = policy.place_demand(fg, 0, 4, flowClass=0, volume=demand)

    # Check LSP distribution
    lsps_per_path = count_lsps_per_path(policy, fg, g, 0, 4)

    # All 3 paths must be used
    assert len(lsps_per_path) == 3, "All 3 paths must be used"

    # Each path should have exactly 3 LSPs (9 / 3 = 3)
    for path_id in [1, 2, 3]:
        lsp_count = len(lsps_per_path.get(path_id, []))
        assert lsp_count == 3, (
            f"Path {path_id} should have 3 LSPs for balanced distribution, "
            f"but has {lsp_count}"
        )

    # With 3 LSPs per path and 100 capacity: 100/3 = 33.33 per LSP
    # Total: 9 * 33.33 = 300 units
    assert placed == pytest.approx(300.0, abs=1.0)


# ============================================================================
# TEST: ECMP Equal-Balancing Constraint
# ============================================================================


def test_ecmp_all_lsps_equal_volume(parallel_3_paths, algs, to_handle):
    """Verify ECMP constraint: all LSPs must carry equal volume."""
    g = parallel_3_paths
    fg = ngc.FlowGraph(g)

    config = ngc.FlowPolicyConfig()
    config.path_alg = ngc.PathAlg.SPF
    config.flow_placement = ngc.FlowPlacement.EQUAL_BALANCED
    config.selection = ngc.EdgeSelection(
        multi_edge=False,
        require_capacity=True,
        tie_break=ngc.EdgeTieBreak.PREFER_HIGHER_RESIDUAL,
    )
    config.multipath = False
    config.min_flow_count = 6
    config.max_flow_count = 6

    policy = ngc.FlowPolicy(algs, to_handle(g), config)

    demand = 240.0
    placed, remaining = policy.place_demand(fg, 0, 4, flowClass=0, volume=demand)

    # Extract flow volumes
    flow_volumes = []
    for flow_id in range(policy.flow_count()):
        flow_key = (0, 4, 0, flow_id)
        if flow_key in policy.flows:
            flow_info = policy.flows[flow_key]
            volume = flow_info[3]
            flow_volumes.append(volume)

    # All flows should have equal volume (ECMP constraint)
    if len(flow_volumes) > 0:
        expected_volume = flow_volumes[0]
        for vol in flow_volumes:
            assert vol == pytest.approx(expected_volume, abs=0.01), (
                f"ECMP violation: LSPs have different volumes: {flow_volumes}"
            )


# ============================================================================
# TEST: Multipath=True Should Use All Paths Differently
# ============================================================================


def test_multipath_true_uses_all_paths(parallel_3_paths, algs, to_handle):
    """With multipath=true, each flow should split across all paths."""
    g = parallel_3_paths
    fg = ngc.FlowGraph(g)

    config = ngc.FlowPolicyConfig()
    config.path_alg = ngc.PathAlg.SPF
    config.flow_placement = ngc.FlowPlacement.EQUAL_BALANCED
    config.selection = ngc.EdgeSelection(
        multi_edge=False,
        require_capacity=True,
        tie_break=ngc.EdgeTieBreak.PREFER_HIGHER_RESIDUAL,
    )
    config.multipath = True  # Each flow SPLITS across all paths
    config.min_flow_count = 6
    config.max_flow_count = 6

    policy = ngc.FlowPolicy(algs, to_handle(g), config)

    demand = 250.0
    placed, remaining = policy.place_demand(fg, 0, 4, flowClass=0, volume=demand)

    # Check path usage
    path_usage = analyze_path_usage(fg, g, 0, 4)

    # All 3 paths must be used
    assert len(path_usage) == 3, (
        f"With multipath=true, all paths should be used, "
        f"but only {len(path_usage)} were used"
    )

    # With multipath=true, flows are evenly distributed
    # 6 flows, each taking 1/3 of each path = 2 "effective" flows per path
    # Each flow: 250/6 = 41.67 units, split 1/3 per path = 13.89 per path per flow
    # Total per path: 6 * 13.89 = 83.33 units

    for flow in path_usage.values():
        assert flow == pytest.approx(83.33, abs=1.0), (
            "With multipath=true, paths should be equally loaded"
        )

    # Should place all demand
    assert placed == pytest.approx(250.0, abs=1.0)


# ============================================================================
# TEST: Progressive Flow Addition
# ============================================================================


def test_incremental_lsp_placement(parallel_3_paths, algs, to_handle):
    """Test placing LSPs incrementally to verify path selection logic."""
    g = parallel_3_paths
    fg = ngc.FlowGraph(g)

    config = ngc.FlowPolicyConfig()
    config.path_alg = ngc.PathAlg.SPF
    config.flow_placement = ngc.FlowPlacement.EQUAL_BALANCED
    config.selection = ngc.EdgeSelection(
        multi_edge=False,
        require_capacity=True,
        tie_break=ngc.EdgeTieBreak.PREFER_HIGHER_RESIDUAL,
    )
    config.multipath = False
    config.min_flow_count = 1
    config.max_flow_count = 1
    config.reoptimize_flows_on_each_placement = False

    # Place LSPs one at a time using DIFFERENT flowClasses (each represents a separate LSP)
    lsp_volume = 30.0
    paths_used_over_time = []

    for i in range(6):
        policy = ngc.FlowPolicy(algs, to_handle(g), config)  # New policy for each "LSP"
        placed, remaining = policy.place_demand(
            fg,
            0,
            4,
            flowClass=i,
            volume=lsp_volume,  # Different flowClass = different LSP
        )

        path_usage = analyze_path_usage(fg, g, 0, 4)
        paths_used_over_time.append(set(path_usage.keys()))

        # After 3 LSPs, all 3 paths should be in use
        if i >= 2:
            assert len(path_usage) == 3, (
                f"After placing {i + 1} LSPs, all 3 paths should be used, "
                f"but only {len(path_usage)} are used: {list(path_usage.keys())}"
            )


# ============================================================================
# TEST: Capacity Constraint Enforcement
# ============================================================================


def test_capacity_limits_with_multipath_false(parallel_3_paths, algs, to_handle):
    """Verify capacity constraints are properly enforced per path."""
    g = parallel_3_paths
    fg = ngc.FlowGraph(g)

    config = ngc.FlowPolicyConfig()
    config.path_alg = ngc.PathAlg.SPF
    config.flow_placement = ngc.FlowPlacement.EQUAL_BALANCED
    config.selection = ngc.EdgeSelection(
        multi_edge=False,
        require_capacity=True,
        tie_break=ngc.EdgeTieBreak.PREFER_HIGHER_RESIDUAL,
    )
    config.multipath = False
    config.min_flow_count = 6
    config.max_flow_count = 6

    policy = ngc.FlowPolicy(algs, to_handle(g), config)

    # Demand exceeds total capacity (3 * 100 = 300)
    demand = 400.0
    placed, remaining = policy.place_demand(fg, 0, 4, flowClass=0, volume=demand)

    # Should place at most 300 (total capacity)
    assert placed <= 300.0, f"Placed {placed} exceeds total capacity of 300"
    assert remaining >= 100.0, f"Should have at least 100 remaining, got {remaining}"

    # Check no path is over-utilized
    edge_flows = fg.edge_flow_view()
    edge_cap = g.capacity_view()

    for edge_id in range(g.num_edges()):
        flow = edge_flows[edge_id]
        cap = edge_cap[edge_id]
        assert flow <= cap + 0.1, (  # Small tolerance for floating point
            f"Edge {edge_id} flow ({flow}) exceeds capacity ({cap})"
        )
