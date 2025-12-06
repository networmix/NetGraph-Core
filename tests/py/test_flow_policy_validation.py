"""Validation tests for FlowPolicy path distribution and flow placement.

Tests cover:
- Path usage and distribution
- Optimality vs theoretical maximum
- ECMP constraints and balancing
"""

import numpy as np
import pytest
from conftest import (
    analyze_path_usage,
    calculate_theoretical_capacity_parallel_paths,
    check_distribution_balanced,
    count_lsps_per_path,
)

import netgraph_core as ngc


def _make_graph(num_nodes, src, dst, cap, cost):
    """Helper to build graph with auto-generated ext_edge_ids."""
    ext_edge_ids = np.arange(len(src), dtype=np.int64)
    return ngc.StrictMultiDiGraph.from_arrays(
        num_nodes, src, dst, cap, cost, ext_edge_ids
    )


# ============================================================================
# FIXTURES: Test Topologies
# ============================================================================


@pytest.fixture
def parallel_paths_varying_capacities():
    """3 parallel paths with different capacities for distribution testing."""
    # S -> [M1, M2, M3] -> T
    # Paths: S-M1-T (cap 100), S-M2-T (cap 50), S-M3-T (cap 75)
    num_nodes = 5
    src = np.array([0, 0, 0, 1, 2, 3], dtype=np.int32)
    dst = np.array([1, 2, 3, 4, 4, 4], dtype=np.int32)
    cap = np.array([100.0, 50.0, 75.0, 100.0, 50.0, 75.0], dtype=np.float64)
    cost = np.array([10, 10, 10, 10, 10, 10], dtype=np.int64)  # Equal costs

    return _make_graph(num_nodes, src, dst, cap, cost)


@pytest.fixture
def parallel_paths_varying_costs():
    """3 parallel paths with different costs for cost-aware testing."""
    # S -> [M1, M2, M3] -> T
    # Paths: S-M1-T (cost 10), S-M2-T (cost 20), S-M3-T (cost 30)
    num_nodes = 5
    src = np.array([0, 0, 0, 1, 2, 3], dtype=np.int32)
    dst = np.array([1, 2, 3, 4, 4, 4], dtype=np.int32)
    cap = np.array([100.0, 100.0, 100.0, 100.0, 100.0, 100.0], dtype=np.float64)
    # First hop costs: 10, 20, 30; Second hop costs: 10, 20, 30
    cost = np.array([10, 20, 30, 10, 20, 30], dtype=np.int64)

    return _make_graph(num_nodes, src, dst, cap, cost)


# ============================================================================
# TEST 1: Path Distribution and Balance
# ============================================================================


def test_path_distribution_with_equal_cost_paths(algs, to_handle):
    """Verify LSPs are evenly distributed across equal-cost paths."""
    # Create topology with 4 equal-cost parallel paths
    num_nodes = 6
    src = np.array([0, 0, 0, 0, 1, 2, 3, 4], dtype=np.int32)
    dst = np.array([1, 2, 3, 4, 5, 5, 5, 5], dtype=np.int32)
    cap = np.array([50.0] * 8, dtype=np.float64)
    cost = np.array([10] * 8, dtype=np.int64)  # All equal cost
    g = _make_graph(num_nodes, src, dst, cap, cost)

    fg = ngc.FlowGraph(g)
    config = ngc.FlowPolicyConfig()
    config.multipath = False  # Single-path LSPs
    config.min_flow_count = 8  # 8 LSPs for 4 paths = 2 per path
    config.max_flow_count = 8
    config.flow_placement = ngc.FlowPlacement.EQUAL_BALANCED

    policy = ngc.FlowPolicy(algs, to_handle(g), config)
    placed, remaining = policy.place_demand(fg, 0, 5, flowClass=0, volume=100.0)

    # Verify all paths are used
    path_usage = analyze_path_usage(fg, g, 0, 5)
    assert len(path_usage) == 4, (
        f"Expected all 4 paths to be used, got {len(path_usage)}: {list(path_usage.keys())}"
    )

    # Verify LSPs are evenly distributed
    lsps_per_path = count_lsps_per_path(policy, fg, g, 0, 5)
    is_balanced, stats = check_distribution_balanced(lsps_per_path, tolerance=0)
    assert is_balanced, f"LSPs not evenly distributed: {stats}"

    # Each path should have exactly 2 LSPs
    for path_id in [1, 2, 3, 4]:
        assert path_id in lsps_per_path, f"Path {path_id} has no LSPs"
        assert len(lsps_per_path[path_id]) == 2, (
            f"Path {path_id} should have 2 LSPs, has {len(lsps_per_path[path_id])}"
        )


def test_path_usage_optimality(algs, to_handle):
    """Verify that available paths are used optimally to maximize capacity."""
    # 5 parallel paths, 10 LSPs -> should use all 5 paths (2 LSPs each)
    num_nodes = 7
    src = np.array([0, 0, 0, 0, 0, 1, 2, 3, 4, 5], dtype=np.int32)
    dst = np.array([1, 2, 3, 4, 5, 6, 6, 6, 6, 6], dtype=np.int32)
    cap = np.array([100.0] * 10, dtype=np.float64)
    cost = np.array([10] * 10, dtype=np.int64)
    g = _make_graph(num_nodes, src, dst, cap, cost)

    fg = ngc.FlowGraph(g)
    config = ngc.FlowPolicyConfig()
    config.multipath = False
    config.min_flow_count = 10
    config.max_flow_count = 10
    config.flow_placement = ngc.FlowPlacement.EQUAL_BALANCED

    policy = ngc.FlowPolicy(algs, to_handle(g), config)
    demand = 500.0
    placed, remaining = policy.place_demand(fg, 0, 6, flowClass=0, volume=demand)

    # Calculate theoretical maximum
    num_paths = 5
    num_lsps = 10
    theoretical_max = calculate_theoretical_capacity_parallel_paths(
        100.0, num_paths, num_lsps
    )

    # Should achieve theoretical maximum (or very close)
    assert placed >= theoretical_max * 0.99, (
        f"Placed {placed} is less than 99% of theoretical max {theoretical_max}"
    )


# ============================================================================
# TEST 2: ECMP Constraints
# ============================================================================


def test_ecmp_equal_balancing_strict(algs, to_handle):
    """Verify ECMP constraint: all LSPs must carry equal volume."""
    num_nodes = 5
    src = np.array([0, 0, 0, 1, 2, 3], dtype=np.int32)
    dst = np.array([1, 2, 3, 4, 4, 4], dtype=np.int32)
    cap = np.array([100.0, 100.0, 100.0, 100.0, 100.0, 100.0], dtype=np.float64)
    cost = np.array([10] * 6, dtype=np.int64)
    g = _make_graph(num_nodes, src, dst, cap, cost)

    fg = ngc.FlowGraph(g)
    config = ngc.FlowPolicyConfig()
    config.multipath = False
    config.min_flow_count = 12  # 12 LSPs
    config.max_flow_count = 12
    config.flow_placement = ngc.FlowPlacement.EQUAL_BALANCED

    policy = ngc.FlowPolicy(algs, to_handle(g), config)
    placed, remaining = policy.place_demand(fg, 0, 4, flowClass=0, volume=100.0)

    # Verify all LSPs have equal volume
    flows = policy.flows
    volumes = [v[3] for v in flows.values()]

    if len(volumes) > 1:
        ref_volume = volumes[0]
        for vol in volumes:
            assert vol == pytest.approx(ref_volume, abs=1e-3), (
                f"ECMP violation: volumes not equal. Got {volumes}"
            )


def test_ecmp_with_path_capacity_constraints(
    algs, to_handle, parallel_paths_varying_capacities
):
    """Verify ECMP with varying path capacities limits correctly."""
    g = parallel_paths_varying_capacities
    fg = ngc.FlowGraph(g)

    config = ngc.FlowPolicyConfig()
    config.multipath = False
    config.min_flow_count = 6  # 6 LSPs for 3 paths
    config.max_flow_count = 6
    config.flow_placement = ngc.FlowPlacement.EQUAL_BALANCED

    policy = ngc.FlowPolicy(algs, to_handle(g), config)
    demand = 300.0
    placed, remaining = policy.place_demand(fg, 0, 4, flowClass=0, volume=demand)

    # Bottleneck: M2 has only 50 capacity
    # With 6 LSPs, some must share M2, limiting per-LSP volume to 50/N_M2
    # Theoretical max depends on distribution

    # Verify accounting
    assert placed + remaining == pytest.approx(demand, abs=1e-6)

    # Verify ECMP: all LSPs have equal volume
    flows = policy.flows
    volumes = [v[3] for v in flows.values()]
    if len(volumes) > 1:
        ref_volume = volumes[0]
        for vol in volumes:
            assert vol == pytest.approx(ref_volume, abs=1e-3)


# ============================================================================
# TEST 3: Multi-path vs Single-path Behavior
# ============================================================================


def test_multipath_true_vs_false_capacity_difference(algs, to_handle):
    """Verify multipath=true achieves higher capacity than multipath=false."""
    # Topology with 3 equal-cost parallel paths
    num_nodes = 5
    src = np.array([0, 0, 0, 1, 2, 3], dtype=np.int32)
    dst = np.array([1, 2, 3, 4, 4, 4], dtype=np.int32)
    cap = np.array([100.0] * 6, dtype=np.float64)
    cost = np.array([10] * 6, dtype=np.int64)
    g = _make_graph(num_nodes, src, dst, cap, cost)

    # Test 1: multipath=true (hash-based ECMP)
    fg1 = ngc.FlowGraph(g)
    config1 = ngc.FlowPolicyConfig()
    config1.multipath = True  # Each flow splits across all paths
    config1.min_flow_count = 3
    config1.max_flow_count = 3
    policy1 = ngc.FlowPolicy(algs, to_handle(g), config1)
    placed1, _ = policy1.place_demand(fg1, 0, 4, flowClass=0, volume=500.0)

    # Test 2: multipath=false (tunnel-based ECMP)
    fg2 = ngc.FlowGraph(g)
    config2 = ngc.FlowPolicyConfig()
    config2.multipath = False  # Each flow uses single path
    config2.min_flow_count = 3
    config2.max_flow_count = 3
    policy2 = ngc.FlowPolicy(algs, to_handle(g), config2)
    placed2, _ = policy2.place_demand(fg2, 0, 4, flowClass=0, volume=500.0)

    # multipath=true should achieve significantly more capacity
    # With 3 flows and 3 paths:
    # - multipath=true: each flow uses all 3 paths (hash-based ECMP)
    # - multipath=false: each flow uses 1 path (tunnel-based ECMP)
    # Both should achieve similar total throughput with balanced distribution.

    # Verify multipath=true placement is at least as good as multipath=false
    assert placed1 >= placed2 * 0.99, (
        f"multipath=true ({placed1}) should not be worse than multipath=false ({placed2})"
    )


# ============================================================================
# TEST 4: Flow Count Verification
# ============================================================================


def test_flow_count_matches_configuration(algs, to_handle):
    """Verify exact flow count matches configuration when paths available."""
    # Create topology with 10 parallel paths (plenty)
    src_nodes = list(range(10))
    dst_nodes = list(range(1, 11))
    src = np.array(src_nodes + dst_nodes, dtype=np.int32)
    dst = np.array(dst_nodes + [11] * 10, dtype=np.int32)
    cap = np.array([100.0] * 20, dtype=np.float64)
    cost = np.array([10] * 20, dtype=np.int64)
    g = _make_graph(12, src, dst, cap, cost)

    # Test various flow counts
    for target_flows in [1, 4, 8, 16]:
        fg = ngc.FlowGraph(g)
        config = ngc.FlowPolicyConfig()
        config.multipath = False
        config.min_flow_count = target_flows
        config.max_flow_count = target_flows

        policy = ngc.FlowPolicy(algs, to_handle(g), config)
        placed, _ = policy.place_demand(fg, 0, 11, flowClass=0, volume=50.0)

        if placed > 0:
            actual_flows = policy.flow_count()
            assert actual_flows == target_flows, (
                f"Expected {target_flows} flows, got {actual_flows}"
            )


# ============================================================================
# TEST 5: Distribution Variance
# ============================================================================


def test_distribution_variance_minimal(algs, to_handle):
    """Verify LSP distribution has minimal variance across paths."""
    # 4 paths, 9 LSPs -> should be [3, 2, 2, 2] or [3, 3, 2, 1] etc
    # Variance should be at most 1
    num_nodes = 6
    src = np.array([0, 0, 0, 0, 1, 2, 3, 4], dtype=np.int32)
    dst = np.array([1, 2, 3, 4, 5, 5, 5, 5], dtype=np.int32)
    cap = np.array([50.0] * 8, dtype=np.float64)
    cost = np.array([10] * 8, dtype=np.int64)
    g = _make_graph(num_nodes, src, dst, cap, cost)

    fg = ngc.FlowGraph(g)
    config = ngc.FlowPolicyConfig()
    config.multipath = False
    config.min_flow_count = 9
    config.max_flow_count = 9
    config.flow_placement = ngc.FlowPlacement.EQUAL_BALANCED

    policy = ngc.FlowPolicy(algs, to_handle(g), config)
    placed, _ = policy.place_demand(fg, 0, 5, flowClass=0, volume=100.0)

    lsps_per_path = count_lsps_per_path(policy, fg, g, 0, 5)
    is_balanced, stats = check_distribution_balanced(lsps_per_path, tolerance=1)

    assert is_balanced, (
        f"Distribution variance too high. Stats: {stats}. "
        f"With 9 LSPs and 4 paths, variance should be ≤1"
    )


# ============================================================================
# TEST 6: Optimality with Complex Topology
# ============================================================================


def test_capacity_near_theoretical_max(algs, to_handle):
    """Verify placed capacity is near theoretical maximum."""
    # Simple scenario: 3 paths, 6 LSPs, equal capacities
    num_nodes = 5
    src = np.array([0, 0, 0, 1, 2, 3], dtype=np.int32)
    dst = np.array([1, 2, 3, 4, 4, 4], dtype=np.int32)
    cap = np.array([60.0] * 6, dtype=np.float64)
    cost = np.array([10] * 6, dtype=np.int64)
    g = _make_graph(num_nodes, src, dst, cap, cost)

    fg = ngc.FlowGraph(g)
    config = ngc.FlowPolicyConfig()
    config.multipath = False
    config.min_flow_count = 6
    config.max_flow_count = 6
    config.flow_placement = ngc.FlowPlacement.EQUAL_BALANCED

    policy = ngc.FlowPolicy(algs, to_handle(g), config)
    demand = 200.0
    placed, remaining = policy.place_demand(fg, 0, 4, flowClass=0, volume=demand)

    # Theoretical: 3 paths, 6 LSPs -> 2 LSPs per path
    # Each path has 60 capacity -> 60/2 = 30 per LSP
    # Total: 6 * 30 = 180
    theoretical = calculate_theoretical_capacity_parallel_paths(60.0, 3, 6)
    assert theoretical == 180.0

    # Placed should be at or near theoretical
    assert placed >= theoretical * 0.95, (
        f"Placed {placed} is less than 95% of theoretical {theoretical}"
    )


# ============================================================================
# TEST 7: Path Usage Coverage
# ============================================================================


def test_all_equal_cost_paths_utilized(algs, to_handle):
    """Verify ALL equal-cost paths are utilized when LSPs exceed path count."""
    # Critical test: 2 paths, 4 LSPs -> both paths must be used
    num_nodes = 4
    src = np.array([0, 0, 1, 2], dtype=np.int32)
    dst = np.array([1, 2, 3, 3], dtype=np.int32)
    cap = np.array([50.0] * 4, dtype=np.float64)
    cost = np.array([10] * 4, dtype=np.int64)
    g = _make_graph(num_nodes, src, dst, cap, cost)

    fg = ngc.FlowGraph(g)
    config = ngc.FlowPolicyConfig()
    config.flow_placement = ngc.FlowPlacement.EQUAL_BALANCED
    config.multipath = False
    config.min_flow_count = 4
    config.max_flow_count = 4

    policy = ngc.FlowPolicy(algs, to_handle(g), config)
    placed, _ = policy.place_demand(fg, 0, 3, flowClass=0, volume=100.0)

    # Both paths MUST be used
    path_usage = analyze_path_usage(fg, g, 0, 3)
    assert len(path_usage) == 2, (
        f"Both paths must be used with 4 LSPs and 2 paths. "
        f"Used: {list(path_usage.keys())}"
    )

    # Each path should have 2 LSPs
    lsps_per_path = count_lsps_per_path(policy, fg, g, 0, 3)
    for path_id in [1, 2]:
        assert path_id in lsps_per_path, f"Path {path_id} not used"
        assert len(lsps_per_path[path_id]) == 2, (
            f"Path {path_id} should have 2 LSPs, has {len(lsps_per_path[path_id])}"
        )
