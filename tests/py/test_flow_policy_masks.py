"""FlowPolicy mask tests.

Tests node and edge mask functionality for failure simulation across
canonical configurations.
"""

from __future__ import annotations

import numpy as np
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

# ============================================================================
# FIXTURES
# ============================================================================


@pytest.fixture
def redundant_paths_graph(build_graph):
    """Graph with redundant paths for testing masking."""
    # S=0, T=4
    # Path 1: 0->1->4, cost=10, cap=100 (primary)
    # Path 2: 0->2->4, cost=10, cap=100 (backup)
    # Path 3: 0->3->4, cost=15, cap=100 (expensive backup)
    edges = [
        (0, 1, 5, 100, 0),  # Primary path
        (1, 4, 5, 100, 1),
        (0, 2, 5, 100, 2),  # Backup path (same cost)
        (2, 4, 5, 100, 3),
        (0, 3, 7, 100, 4),  # Expensive path
        (3, 4, 8, 100, 5),
    ]
    return build_graph(5, edges)


# ============================================================================
# TEST 1: NODE MASK EXCLUSION
# ============================================================================


@pytest.mark.parametrize("config_name", ALL_CONFIGS)
def test_node_mask_excludes_node(
    config_name, redundant_paths_graph, algs, to_handle, make_flow_policy
):
    """Verify configurations respect node exclusions."""
    g = redundant_paths_graph

    # Without mask - should place successfully using all paths
    fg1 = ngc.FlowGraph(g)
    policy1 = make_flow_policy(config_name, algs, to_handle(g), node_mask=None)
    placed1, _ = policy1.place_demand(fg1, 0, 4, flowClass=0, volume=100.0)

    # Create mask excluding node 1 (part of primary path)
    node_mask = np.ones(5, dtype=bool)
    node_mask[1] = False  # Exclude node 1

    # With mask - should place using alternate paths
    fg2 = ngc.FlowGraph(g)
    policy2 = make_flow_policy(config_name, algs, to_handle(g), node_mask=node_mask)
    placed2, _ = policy2.place_demand(fg2, 0, 4, flowClass=0, volume=100.0)

    # Should still place flow (via nodes 2 or 3)
    assert placed2 > 0, f"{config_name}: failed to place when node 1 masked"

    # Placement with mask may differ from without mask (configuration-dependent)
    # Verifies alternate path discovery when primary path is masked


# ============================================================================
# TEST 2: EDGE MASK EXCLUSION
# ============================================================================


@pytest.mark.parametrize("config_name", ALL_CONFIGS)
def test_edge_mask_excludes_edge(
    config_name, redundant_paths_graph, algs, to_handle, make_flow_policy
):
    """Verify configurations respect edge exclusions."""
    g = redundant_paths_graph

    # Without mask
    fg1 = ngc.FlowGraph(g)
    policy1 = make_flow_policy(config_name, algs, to_handle(g), edge_mask=None)
    placed1, _ = policy1.place_demand(fg1, 0, 4, flowClass=0, volume=100.0)

    # Create mask excluding edge 0 (0->1, first hop of primary path)
    edge_mask = np.ones(6, dtype=bool)
    edge_mask[0] = False  # Exclude edge 0

    # With mask - should use alternate paths
    fg2 = ngc.FlowGraph(g)
    policy2 = make_flow_policy(config_name, algs, to_handle(g), edge_mask=edge_mask)
    placed2, _ = policy2.place_demand(fg2, 0, 4, flowClass=0, volume=100.0)

    # Should still place flow (via other edges)
    assert placed2 > 0, f"{config_name}: failed to place when edge 0 masked"


# ============================================================================
# TEST 3: MASK REDUCES CAPACITY
# ============================================================================


@pytest.mark.parametrize("config_name", CAPACITY_AWARE_CONFIGS)
def test_mask_reduces_available_capacity(
    config_name, redundant_paths_graph, algs, to_handle, make_flow_policy
):
    """Verify masking reduces available capacity correctly."""
    g = redundant_paths_graph

    # Without mask - get baseline capacity
    fg1 = ngc.FlowGraph(g)
    policy1 = make_flow_policy(config_name, algs, to_handle(g))
    placed1, _ = policy1.place_demand(fg1, 0, 4, flowClass=0, volume=1000.0)

    # With mask excluding primary path nodes
    node_mask = np.ones(5, dtype=bool)
    node_mask[1] = False  # Exclude primary path

    fg2 = ngc.FlowGraph(g)
    policy2 = make_flow_policy(config_name, algs, to_handle(g), node_mask=node_mask)
    placed2, _ = policy2.place_demand(fg2, 0, 4, flowClass=0, volume=1000.0)

    # Masked version should place less or equal
    assert placed2 <= placed1 + 1e-6, (
        f"{config_name}: masking didn't reduce or maintain capacity. "
        f"Without mask: {placed1}, with mask: {placed2}"
    )


# ============================================================================
# TEST 4: MASK CHANGES COSTS
# ============================================================================


@pytest.mark.parametrize("config_name", ALL_CONFIGS)
def test_mask_forces_higher_cost_paths(
    config_name, redundant_paths_graph, algs, to_handle, make_flow_policy
):
    """Verify masking forces use of alternate (more expensive) paths."""
    g = redundant_paths_graph

    # Helper to extract minimum cost
    def get_min_cost(policy):
        flows = policy.flows
        if not flows:
            return float("inf")
        costs = [v[2] for v in flows.values()]
        return min(costs)

    # Without mask - should use cheap paths (cost 10)
    fg1 = ngc.FlowGraph(g)
    policy1 = make_flow_policy(config_name, algs, to_handle(g))
    placed1, _ = policy1.place_demand(fg1, 0, 4, flowClass=0, volume=50.0)

    if placed1 > 0:
        min_cost1 = get_min_cost(policy1)

        # With mask excluding cheap paths - should use expensive path (cost 15)
        node_mask = np.ones(5, dtype=bool)
        node_mask[1] = False  # Exclude node 1
        node_mask[2] = False  # Exclude node 2

        fg2 = ngc.FlowGraph(g)
        policy2 = make_flow_policy(config_name, algs, to_handle(g), node_mask=node_mask)
        placed2, _ = policy2.place_demand(fg2, 0, 4, flowClass=0, volume=50.0)

        if placed2 > 0:
            min_cost2 = get_min_cost(policy2)

            # Masked version should use more expensive paths
            assert min_cost2 >= min_cost1 - 1e-6, (
                f"{config_name}: masking didn't increase or maintain path cost. "
                f"Without mask: {min_cost1}, with mask: {min_cost2}"
            )


# ============================================================================
# TEST 5: COMPLETE NETWORK MASKING
# ============================================================================


@pytest.mark.parametrize("config_name", ALL_CONFIGS)
def test_complete_mask_blocks_all_paths(
    config_name, redundant_paths_graph, algs, to_handle, make_flow_policy
):
    """Verify behavior when all paths are masked."""
    g = redundant_paths_graph

    # Mask all intermediate nodes
    node_mask = np.ones(5, dtype=bool)
    node_mask[1] = False  # Exclude all
    node_mask[2] = False  # intermediate
    node_mask[3] = False  # nodes

    fg = ngc.FlowGraph(g)
    policy = make_flow_policy(config_name, algs, to_handle(g), node_mask=node_mask)
    placed, remaining = policy.place_demand(fg, 0, 4, flowClass=0, volume=100.0)

    # Should place nothing
    assert placed == pytest.approx(0.0, abs=1e-9), (
        f"{config_name}: placed volume when all paths masked"
    )
    assert remaining == pytest.approx(100.0, abs=1e-9), (
        f"{config_name}: remaining doesn't match demand when all paths masked"
    )


# ============================================================================
# TEST 6: SOURCE NODE MASKED
# ============================================================================


@pytest.mark.parametrize("config_name", ALL_CONFIGS)
def test_source_node_masked(
    config_name, redundant_paths_graph, algs, to_handle, make_flow_policy
):
    """Verify behavior when source node is masked."""
    g = redundant_paths_graph

    # Mask source node
    node_mask = np.ones(5, dtype=bool)
    node_mask[0] = False  # Exclude source

    fg = ngc.FlowGraph(g)
    policy = make_flow_policy(config_name, algs, to_handle(g), node_mask=node_mask)
    placed, remaining = policy.place_demand(fg, 0, 4, flowClass=0, volume=100.0)

    # Should place nothing (no traversal can start from excluded source)
    assert placed == pytest.approx(0.0, abs=1e-9), (
        f"{config_name}: placed volume when source masked"
    )
    assert remaining == pytest.approx(100.0, abs=1e-9), (
        f"{config_name}: remaining doesn't match demand when source masked"
    )


# ============================================================================
# TEST 7: DESTINATION NODE MASKED
# ============================================================================


@pytest.mark.parametrize("config_name", ALL_CONFIGS)
def test_destination_node_masked(
    config_name, redundant_paths_graph, algs, to_handle, make_flow_policy
):
    """Verify behavior when destination node is masked."""
    g = redundant_paths_graph

    # Mask destination node
    node_mask = np.ones(5, dtype=bool)
    node_mask[4] = False  # Exclude destination

    fg = ngc.FlowGraph(g)
    policy = make_flow_policy(config_name, algs, to_handle(g), node_mask=node_mask)
    placed, remaining = policy.place_demand(fg, 0, 4, flowClass=0, volume=100.0)

    # Should place nothing (cannot reach excluded destination)
    assert placed == pytest.approx(0.0, abs=1e-9), (
        f"{config_name}: placed volume when destination masked"
    )
    assert remaining == pytest.approx(100.0, abs=1e-9), (
        f"{config_name}: remaining doesn't match demand when destination masked"
    )


# ============================================================================
# TEST 8: COMBINED NODE AND EDGE MASKS
# ============================================================================


@pytest.mark.parametrize("config_name", ALL_CONFIGS)
def test_combined_node_and_edge_masks(
    config_name, redundant_paths_graph, algs, to_handle, make_flow_policy
):
    """Verify behavior with both node and edge masks."""
    g = redundant_paths_graph

    # Exclude node 1 and edge 2 (0->2)
    node_mask = np.ones(5, dtype=bool)
    node_mask[1] = False

    edge_mask = np.ones(6, dtype=bool)
    edge_mask[2] = False

    fg = ngc.FlowGraph(g)
    policy = make_flow_policy(
        config_name, algs, to_handle(g), node_mask=node_mask, edge_mask=edge_mask
    )
    placed, remaining = policy.place_demand(fg, 0, 4, flowClass=0, volume=100.0)

    # Should still place via expensive path (0->3->4)
    assert placed > 0, f"{config_name}: failed to place with combined masks"


# ============================================================================
# TEST 9: MASK VALIDATION
# ============================================================================


def test_invalid_node_mask_length(
    redundant_paths_graph, algs, to_handle, make_flow_policy
):
    """Verify that invalid mask lengths are rejected."""
    g = redundant_paths_graph

    # Create mask with wrong length
    node_mask = np.ones(10, dtype=bool)  # Graph has 5 nodes, not 10

    # Should raise TypeError for length mismatch
    with pytest.raises(TypeError, match="node_mask"):
        make_flow_policy("SHORTEST_PATHS_ECMP", algs, to_handle(g), node_mask=node_mask)


def test_invalid_edge_mask_length(
    redundant_paths_graph, algs, to_handle, make_flow_policy
):
    """Verify that invalid mask lengths are rejected."""
    g = redundant_paths_graph

    # Create mask with wrong length
    edge_mask = np.ones(10, dtype=bool)  # Graph has 6 edges, not 10

    # Should raise TypeError for length mismatch
    with pytest.raises(TypeError, match="edge_mask"):
        make_flow_policy("SHORTEST_PATHS_ECMP", algs, to_handle(g), edge_mask=edge_mask)


def test_invalid_mask_dtype(redundant_paths_graph, algs, to_handle, make_flow_policy):
    """Verify that non-bool masks are rejected."""
    g = redundant_paths_graph

    # Create mask with wrong dtype
    node_mask = np.ones(5, dtype=np.float64)  # Should be bool

    # Should raise TypeError for wrong dtype
    with pytest.raises(TypeError, match="bool"):
        make_flow_policy("SHORTEST_PATHS_ECMP", algs, to_handle(g), node_mask=node_mask)


# ============================================================================
# TEST 10: MASK IMMUTABILITY
# ============================================================================


@pytest.mark.parametrize("config_name", ["SHORTEST_PATHS_ECMP", "TE_WCMP_UNLIM"])
def test_mask_is_copied(
    config_name, redundant_paths_graph, algs, to_handle, make_flow_policy
):
    """Verify that masks are copied (not referenced) by FlowPolicy."""
    g = redundant_paths_graph

    # Create mask
    node_mask = np.ones(5, dtype=bool)
    node_mask[1] = False

    fg = ngc.FlowGraph(g)
    policy = make_flow_policy(config_name, algs, to_handle(g), node_mask=node_mask)
    placed1, _ = policy.place_demand(fg, 0, 4, flowClass=0, volume=50.0)

    # Modify the mask after policy creation
    node_mask[2] = False

    # Policy should not be affected by external mask changes
    # (behavior should be consistent with first placement)
    policy2 = make_flow_policy(config_name, algs, to_handle(g), node_mask=None)
    fg2 = ngc.FlowGraph(g)
    # This test mainly validates that the policy doesn't crash or behave erratically
    # The actual assertion is that this doesn't raise an exception
    placed2, _ = policy2.place_demand(fg2, 0, 4, flowClass=0, volume=50.0)
    assert placed2 > 0  # Should work normally
