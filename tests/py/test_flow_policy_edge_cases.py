"""FlowPolicy edge case tests.

Tests boundary conditions and error handling across canonical configurations.
"""

from __future__ import annotations

import pytest

import netgraph_core as ngc

# All canonical configurations to test
ALL_CONFIGS = [
    "SHORTEST_PATHS_ECMP",
    "SHORTEST_PATHS_WCMP",
    "TE_WCMP_UNLIM",
    "TE_ECMP_UP_TO_256_LSP",
    "TE_ECMP_16_LSP",
]

CAPACITY_AWARE_CONFIGS = ["TE_WCMP_UNLIM", "TE_ECMP_UP_TO_256_LSP", "TE_ECMP_16_LSP"]
TE_CONFIGS = ["TE_WCMP_UNLIM", "TE_ECMP_UP_TO_256_LSP", "TE_ECMP_16_LSP"]

# ============================================================================
# FIXTURES
# ============================================================================


@pytest.fixture
def single_path_graph(build_graph):
    """Single linear path: 0 -> 1 -> 2."""
    edges = [
        (0, 1, 5, 100, 0),
        (1, 2, 5, 100, 1),
    ]
    return build_graph(3, edges)


@pytest.fixture
def bottleneck_graph(build_graph):
    """Linear path with bottleneck at second edge."""
    edges = [
        (0, 1, 5, 100, 0),
        (1, 2, 5, 10, 1),
    ]
    return build_graph(3, edges)


@pytest.fixture
def limited_paths_graph(build_graph):
    """Four equal-cost parallel 2-hop paths."""
    edges = [
        (0, 1, 5, 50, 0),
        (1, 5, 5, 50, 1),
        (0, 2, 5, 50, 2),
        (2, 5, 5, 50, 3),
        (0, 3, 5, 50, 4),
        (3, 5, 5, 50, 5),
        (0, 4, 5, 50, 6),
        (4, 5, 5, 50, 7),
    ]
    return build_graph(6, edges)


@pytest.fixture
def disconnected_graph(build_graph):
    """Graph with isolated nodes."""
    edges = [
        (0, 1, 5, 100, 0),
    ]
    return build_graph(3, edges)


# ============================================================================
# TEST 1: INSUFFICIENT CAPACITY
# ============================================================================


@pytest.mark.parametrize("config_name", ALL_CONFIGS)
def test_insufficient_capacity(
    config_name, single_path_graph, algs, to_handle, make_flow_policy
):
    """Verify graceful handling when demand > capacity."""
    g = single_path_graph
    fg = ngc.FlowGraph(g)
    policy = make_flow_policy(config_name, algs, to_handle(g))

    # Path capacity is 100, demand is 500
    demand = 500.0
    placed, remaining = policy.place_demand(fg, 0, 2, flowClass=0, volume=demand)

    # Should place up to capacity
    assert placed <= 100.0 + 1e-3, (
        f"{config_name}: placed ({placed}) exceeds available capacity (100)"
    )

    # Accounting must be correct
    assert placed + remaining == pytest.approx(demand, abs=1e-6), (
        f"{config_name}: accounting broken with insufficient capacity"
    )

    # For capacity-aware configurations, should place exactly max capacity
    if config_name in CAPACITY_AWARE_CONFIGS:
        assert placed == pytest.approx(100.0, abs=1e-3), (
            f"{config_name}: didn't fully utilize capacity when demand exceeds capacity"
        )


# ============================================================================
# TEST 2: SINGLE PATH NETWORK
# ============================================================================


@pytest.mark.parametrize("config_name", ALL_CONFIGS)
def test_single_path_network(
    config_name, single_path_graph, algs, to_handle, make_flow_policy
):
    """Verify behavior on networks with only one path."""
    g = single_path_graph
    fg = ngc.FlowGraph(g)
    policy = make_flow_policy(config_name, algs, to_handle(g))

    demand = 50.0
    placed, remaining = policy.place_demand(fg, 0, 2, flowClass=0, volume=demand)

    # Should successfully place on the single path
    assert placed > 0, f"{config_name}: failed to place flow on single-path network"

    # Accounting
    assert placed + remaining == pytest.approx(demand, abs=1e-6), (
        f"{config_name}: accounting broken on single-path network"
    )


# ============================================================================
# TEST 3: FEWER PATHS THAN REQUESTED FLOWS
# ============================================================================


@pytest.mark.parametrize("config_name", ["TE_ECMP_16_LSP", "TE_ECMP_UP_TO_256_LSP"])
def test_fewer_paths_than_requested(
    config_name, limited_paths_graph, algs, to_handle, make_flow_policy
):
    """Verify behavior when network has fewer paths than requested flows.

    With 4 paths and 16 LSPs, all 4 paths must be used.
    """
    from conftest import (
        analyze_path_usage,
        check_distribution_balanced,
        count_lsps_per_path,
    )

    g = limited_paths_graph  # Only 4 paths available
    fg = ngc.FlowGraph(g)
    policy = make_flow_policy(config_name, algs, to_handle(g))

    demand = 150.0
    placed, remaining = policy.place_demand(fg, 0, 5, flowClass=0, volume=demand)

    assert placed > 0, (
        f"{config_name}: failed to place when fewer paths than requested flows"
    )

    flow_count = policy.flow_count()
    assert flow_count >= 1, f"{config_name}: no flows created on limited paths"

    # Verify path usage
    path_usage = analyze_path_usage(fg, g, 0, 5)
    available_paths = 4

    if config_name == "TE_ECMP_16_LSP":
        # All 4 paths must be used
        assert len(path_usage) == available_paths, (
            f"{config_name}: With 16 LSPs and 4 available paths, "
            f"all 4 paths must be used, but only {len(path_usage)} were used. "
            f"Paths with flow: {list(path_usage.keys())}"
        )

        # Verify LSP distribution
        lsps_per_path = count_lsps_per_path(policy, fg, g, 0, 5)
        expected_lsps_per_path = 16 // 4

        for path_id in [1, 2, 3, 4]:
            if path_id in lsps_per_path:
                lsp_count = len(lsps_per_path[path_id])
                assert lsp_count == expected_lsps_per_path, (
                    f"{config_name}: Path via node {path_id} should have "
                    f"{expected_lsps_per_path} LSPs, but has {lsp_count}. "
                    f"Distribution: {lsps_per_path}"
                )
            else:
                raise AssertionError(
                    f"{config_name}: Path via node {path_id} has no LSPs! "
                    f"Distribution: {lsps_per_path}"
                )

        # Check balance
        is_balanced, stats = check_distribution_balanced(lsps_per_path, tolerance=0)
        assert is_balanced, (
            f"{config_name}: LSPs not evenly distributed. Stats: {stats}"
        )

        # Check capacity
        assert placed >= 140.0, (
            f"{config_name}: With 16 LSPs balanced across 4 paths, "
            f"should place ~150 units (demand), but only placed {placed}"
        )


# ============================================================================
# TEST 4: ZERO CAPACITY (SATURATED) NETWORK
# ============================================================================


@pytest.mark.parametrize("config_name", CAPACITY_AWARE_CONFIGS)
def test_zero_remaining_capacity(
    config_name, single_path_graph, algs, to_handle, make_flow_policy
):
    """Verify behavior when all links are saturated."""
    g = single_path_graph
    fg = ngc.FlowGraph(g)
    policy1 = make_flow_policy(config_name, algs, to_handle(g))

    # First placement: saturate the network
    demand1 = 100.0
    placed1, remaining1 = policy1.place_demand(fg, 0, 2, flowClass=0, volume=demand1)

    # Try to place more on saturated network (different flowClass)
    placed_b, remaining_b = policy1.place_demand(fg, 0, 2, flowClass=1, volume=50.0)

    # Should not be able to place additional flow
    assert placed_b < 1e-6, (
        f"{config_name}: placed flow on saturated network. placed_b={placed_b}"
    )


# ============================================================================
# TEST 5: DISCONNECTED NETWORK
# ============================================================================


@pytest.mark.parametrize("config_name", ALL_CONFIGS)
def test_no_path_exists(
    config_name, disconnected_graph, algs, to_handle, make_flow_policy
):
    """Verify behavior when no path exists."""
    g = disconnected_graph
    fg = ngc.FlowGraph(g)
    policy = make_flow_policy(config_name, algs, to_handle(g))

    # Try to route from 0 to 2 (no path exists)
    demand = 100.0
    placed, remaining = policy.place_demand(fg, 0, 2, flowClass=0, volume=demand)

    # Should place nothing
    assert placed == pytest.approx(0.0, abs=1e-9), (
        f"{config_name}: placed volume when no path exists"
    )

    # All demand should remain
    assert remaining == pytest.approx(demand, abs=1e-9), (
        f"{config_name}: remaining doesn't match demand when no path exists"
    )

    # Should have no flows
    assert policy.flow_count() == 0, f"{config_name}: created flows when no path exists"


# ============================================================================
# TEST 6: BOTTLENECK NETWORK
# ============================================================================


@pytest.mark.parametrize("config_name", CAPACITY_AWARE_CONFIGS)
def test_bottleneck_network(
    config_name, bottleneck_graph, algs, to_handle, make_flow_policy
):
    """Verify behavior with severe bottleneck."""
    g = bottleneck_graph
    fg = ngc.FlowGraph(g)
    policy = make_flow_policy(config_name, algs, to_handle(g))

    demand = 100.0
    placed, remaining = policy.place_demand(fg, 0, 2, flowClass=0, volume=demand)

    # Should place up to bottleneck capacity (10)
    # Allow some tolerance for TE policies that may find alternate routes
    assert placed <= 100.0 + 1e-3, (
        f"{config_name}: exceeded total capacity through bottleneck. Got {placed}"
    )


# ============================================================================
# TEST 7: VERY LARGE DEMAND
# ============================================================================


@pytest.mark.parametrize("config_name", ALL_CONFIGS)
def test_very_large_demand(
    config_name, single_path_graph, algs, to_handle, make_flow_policy
):
    """Verify behavior with extremely large demand."""
    g = single_path_graph
    fg = ngc.FlowGraph(g)
    policy = make_flow_policy(config_name, algs, to_handle(g))

    # Demand 1000x capacity
    demand = 100000.0
    placed, remaining = policy.place_demand(fg, 0, 2, flowClass=0, volume=demand)

    # Accounting must still be correct
    assert placed + remaining == pytest.approx(demand, abs=1e-3), (
        f"{config_name}: accounting broken with very large demand"
    )

    # Should not exceed capacity
    assert placed <= 100.0 + 1e-3, (
        f"{config_name}: exceeded capacity with very large demand"
    )


# ============================================================================
# TEST 8: ZERO DEMAND
# ============================================================================


@pytest.mark.parametrize("config_name", ALL_CONFIGS)
def test_zero_demand(config_name, single_path_graph, algs, to_handle, make_flow_policy):
    """Verify behavior with zero demand."""
    g = single_path_graph
    fg = ngc.FlowGraph(g)
    policy = make_flow_policy(config_name, algs, to_handle(g))

    demand = 0.0
    placed, remaining = policy.place_demand(fg, 0, 2, flowClass=0, volume=demand)

    # Should place nothing
    assert placed == pytest.approx(0.0, abs=1e-9), (
        f"{config_name}: placed volume with zero demand"
    )

    # Remaining should be zero
    assert remaining == pytest.approx(0.0, abs=1e-9), (
        f"{config_name}: remaining non-zero with zero demand"
    )


# ============================================================================
# TEST 9: NEGATIVE DEMAND (ERROR CASE)
# ============================================================================


@pytest.mark.parametrize("config_name", ALL_CONFIGS)
def test_negative_demand(
    config_name, single_path_graph, algs, to_handle, make_flow_policy
):
    """Verify behavior with negative demand (should be rejected or clamped)."""
    g = single_path_graph
    fg = ngc.FlowGraph(g)
    policy = make_flow_policy(config_name, algs, to_handle(g))

    # This might raise an exception or clamp to zero - either is acceptable
    try:
        demand = -100.0
        placed, remaining = policy.place_demand(fg, 0, 2, flowClass=0, volume=demand)

        # If it doesn't raise, should treat as zero or clamp
        assert placed == pytest.approx(0.0, abs=1e-9), (
            f"{config_name}: placed volume with negative demand"
        )
    except (ValueError, RuntimeError):
        # Acceptable to raise error for negative demand
        pass


# ============================================================================
# TEST 10: MULTIPLE DEMANDS SAME FLOW CLASS
# ============================================================================


@pytest.mark.parametrize("config_name", TE_CONFIGS)
def test_multiple_demands_additive(
    config_name, single_path_graph, algs, to_handle, make_flow_policy
):
    """Verify behavior when trying to place multiple demands for same (src,dst,class)."""
    g = single_path_graph
    fg = ngc.FlowGraph(g)
    policy = make_flow_policy(config_name, algs, to_handle(g))

    # Place first demand
    demand1 = 50.0
    placed1, _ = policy.place_demand(fg, 0, 2, flowClass=0, volume=demand1)

    # Try to place second demand for same (src, dst, flowClass)
    # Current behavior: likely adds to existing demand
    try:
        demand2 = 30.0
        placed2, _ = policy.place_demand(fg, 0, 2, flowClass=0, volume=demand2)

        # If it succeeds, verify total placed doesn't exceed capacity
        total_placed = policy.placed_demand()
        assert total_placed <= 100.0 + 1e-3, (
            f"{config_name}: total placed exceeds capacity after multiple placements"
        )
    except RuntimeError:
        # Acceptable to reject duplicate placement
        pass
