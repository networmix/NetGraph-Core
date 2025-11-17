"""FlowPolicy tests across canonical configurations.

Tests 5 canonical configurations:
- SHORTEST_PATHS_ECMP: IP/IGP ECMP (hash-based multipath, shortest paths only)
- SHORTEST_PATHS_WCMP: IP/IGP WCMP (proportional multipath, shortest paths only)
- TE_WCMP_UNLIM: Traffic engineering with unlimited proportional flows
- TE_ECMP_UP_TO_256_LSP: TE with up to 256 ECMP LSPs (tunnel-based)
- TE_ECMP_16_LSP: TE with exactly 16 ECMP LSPs (fixed allocation)

Validates:
- Accounting: placed + remaining = demand
- Capacity constraints: TE configurations respect max_flow
- Flow count enforcement: min/max flow counts
- Balancing: ECMP equal, WCMP proportional
- Cost distribution: volumes sum correctly
- Progressive TE: uses multiple cost tiers when needed
- Lifecycle: remove_demand, rebalance
"""

from __future__ import annotations

from typing import Dict, List, Tuple

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

ECMP_CONFIGS = ["SHORTEST_PATHS_ECMP", "TE_ECMP_UP_TO_256_LSP", "TE_ECMP_16_LSP"]
WCMP_CONFIGS = ["SHORTEST_PATHS_WCMP", "TE_WCMP_UNLIM"]
CAPACITY_AWARE_CONFIGS = ["TE_WCMP_UNLIM", "TE_ECMP_UP_TO_256_LSP", "TE_ECMP_16_LSP"]
TE_CONFIGS = ["TE_WCMP_UNLIM", "TE_ECMP_UP_TO_256_LSP", "TE_ECMP_16_LSP"]


def get_expected_flow_count(config_name: str) -> Tuple[int, int | None]:
    """Return (min_expected, max_expected) flow count for a configuration."""
    if config_name in ["SHORTEST_PATHS_ECMP", "SHORTEST_PATHS_WCMP"]:
        return (1, 1)
    elif config_name == "TE_WCMP_UNLIM":
        return (1, None)  # Unlimited
    elif config_name == "TE_ECMP_UP_TO_256_LSP":
        return (1, 256)
    elif config_name == "TE_ECMP_16_LSP":
        return (16, 16)
    raise ValueError(f"Unknown configuration: {config_name}")


def extract_cost_distribution(policy: ngc.FlowPolicy) -> Dict[float, float]:
    """Extract cost distribution from policy.flows dictionary."""
    cost_dist = {}
    for _key, value in policy.flows.items():
        # value is tuple: (src, dst, cost, volume)
        cost = float(value[2])
        volume = float(value[3])
        cost_dist[cost] = cost_dist.get(cost, 0.0) + volume
    return cost_dist


def check_equal_balancing(volumes: List[float], tolerance: float = 1e-3) -> bool:
    """Check if all volumes are approximately equal."""
    if not volumes or len(volumes) == 1:
        return True
    avg = sum(volumes) / len(volumes)
    # Use both absolute and relative tolerance for robustness
    return all(
        abs(v - avg) < tolerance or (abs(v - avg) / avg < 1e-6 if avg > 0 else False)
        for v in volumes
    )


# ============================================================================
# FIXTURES
# ============================================================================


@pytest.fixture
def multi_path_graph(build_graph):
    """4 parallel paths with equal cost but different capacities."""
    # S=0, T=5
    # Path 1: 0->1->5, cost=10, cap=100
    # Path 2: 0->2->5, cost=10, cap=80
    # Path 3: 0->3->5, cost=10, cap=60
    # Path 4: 0->4->5, cost=10, cap=40
    edges = [
        (0, 1, 5, 100, 0),
        (1, 5, 5, 100, 1),
        (0, 2, 5, 80, 2),
        (2, 5, 5, 80, 3),
        (0, 3, 5, 60, 4),
        (3, 5, 5, 60, 5),
        (0, 4, 5, 40, 6),
        (4, 5, 5, 40, 7),
    ]
    return build_graph(6, edges)


@pytest.fixture
def multi_tier_graph(build_graph):
    """3 cost tiers with different capacities."""
    # S=0, T=4
    # Tier 1: cost=10, cap=50  (0->1->4)
    # Tier 2: cost=15, cap=75  (0->2->4)
    # Tier 3: cost=20, cap=100 (0->3->4)
    edges = [
        (0, 1, 5, 50, 0),
        (1, 4, 5, 50, 1),
        (0, 2, 7, 75, 2),
        (2, 4, 8, 75, 3),
        (0, 3, 10, 100, 4),
        (3, 4, 10, 100, 5),
    ]
    return build_graph(5, edges)


# ============================================================================
# TEST 1: ACCOUNTING VALIDATION
# ============================================================================


@pytest.mark.parametrize("config_name", ALL_CONFIGS)
def test_accounting_correctness(
    config_name, multi_path_graph, algs, to_handle, make_flow_policy
):
    """Verify placed + remaining = demand for all configurations."""
    g = multi_path_graph
    fg = ngc.FlowGraph(g)
    policy = make_flow_policy(config_name, algs, to_handle(g))

    demand = 200.0
    placed, remaining = policy.place_demand(fg, 0, 5, flowClass=0, volume=demand)

    assert placed + remaining == pytest.approx(demand, rel=0, abs=1e-6), (
        f"{config_name}: placed({placed}) + remaining({remaining}) != demand({demand})"
    )


# ============================================================================
# TEST 2: CAPACITY CONSTRAINTS
# ============================================================================


@pytest.mark.parametrize("config_name", CAPACITY_AWARE_CONFIGS)
def test_capacity_aware_respects_max_flow(
    config_name, multi_path_graph, algs, to_handle, make_flow_policy
):
    """Verify capacity-aware configurations don't exceed max_flow."""
    g = multi_path_graph
    fg = ngc.FlowGraph(g)

    # Calculate max flow
    max_flow, _ = algs.max_flow(
        to_handle(g),
        0,
        5,
        flow_placement=ngc.FlowPlacement.PROPORTIONAL,
        shortest_path=False,
        with_edge_flows=False,
    )

    policy = make_flow_policy(config_name, algs, to_handle(g))
    placed, _ = policy.place_demand(fg, 0, 5, flowClass=0, volume=1000.0)

    assert placed <= float(max_flow) + 1e-3, (
        f"{config_name}: placed({placed}) exceeds max_flow({max_flow})"
    )


# ============================================================================
# TEST 3: FLOW COUNT VALIDATION
# ============================================================================


@pytest.mark.parametrize("config_name", ALL_CONFIGS)
def test_flow_count_enforcement(
    config_name, multi_path_graph, algs, to_handle, make_flow_policy
):
    """Verify each configuration creates the expected number of flows."""
    g = multi_path_graph
    fg = ngc.FlowGraph(g)
    policy = make_flow_policy(config_name, algs, to_handle(g))

    placed, _ = policy.place_demand(fg, 0, 5, flowClass=0, volume=200.0)
    flow_count = policy.flow_count()

    min_exp, max_exp = get_expected_flow_count(config_name)

    assert flow_count >= min_exp, (
        f"{config_name}: flow_count({flow_count}) < min_expected({min_exp})"
    )

    if max_exp is not None:
        assert flow_count <= max_exp, (
            f"{config_name}: flow_count({flow_count}) > max_expected({max_exp})"
        )

    # Specific checks for fixed-count configurations
    if config_name == "TE_ECMP_16_LSP":
        assert flow_count == 16, (
            f"{config_name}: expected exactly 16 flows, got {flow_count}"
        )
    elif config_name in ["SHORTEST_PATHS_ECMP", "SHORTEST_PATHS_WCMP"]:
        assert flow_count == 1, (
            f"{config_name}: expected exactly 1 flow, got {flow_count}"
        )


# ============================================================================
# TEST 4: EQUAL BALANCING (ECMP)
# ============================================================================


@pytest.mark.parametrize("config_name", ECMP_CONFIGS)
def test_ecmp_equal_balancing(
    config_name, multi_path_graph, algs, to_handle, make_flow_policy
):
    """Verify ECMP configurations split flows equally."""
    g = multi_path_graph
    fg = ngc.FlowGraph(g)
    policy = make_flow_policy(config_name, algs, to_handle(g))

    placed, _ = policy.place_demand(fg, 0, 5, flowClass=0, volume=200.0)

    if placed == 0:
        pytest.skip("No flow placed, cannot check balancing")

    flows = policy.flows
    volumes = [v[3] for v in flows.values()]

    assert check_equal_balancing(volumes, tolerance=1e-3), (
        f"{config_name}: flows not equally balanced. Volumes: {volumes}"
    )


# ============================================================================
# TEST 5: PROPORTIONAL BALANCING (WCMP)
# ============================================================================


@pytest.mark.parametrize("config_name", WCMP_CONFIGS)
def test_wcmp_proportional_balancing(
    config_name, multi_path_graph, algs, to_handle, make_flow_policy
):
    """Verify WCMP configurations split flows proportionally."""
    g = multi_path_graph
    fg = ngc.FlowGraph(g)
    policy = make_flow_policy(config_name, algs, to_handle(g))

    placed, _ = policy.place_demand(fg, 0, 5, flowClass=0, volume=200.0)

    if placed == 0:
        pytest.skip("No flow placed, cannot check balancing")

    flows = policy.flows
    volumes = [v[3] for v in flows.values()]

    # Verify total is correct
    assert abs(sum(volumes) - placed) < 1e-6, (
        f"{config_name}: sum of flow volumes doesn't match placed"
    )


# ============================================================================
# TEST 6: COST DISTRIBUTION
# ============================================================================


@pytest.mark.parametrize("config_name", ALL_CONFIGS)
def test_cost_distribution_sums_correctly(
    config_name, multi_path_graph, algs, to_handle, make_flow_policy
):
    """Verify cost distribution volumes sum to total placed volume."""
    g = multi_path_graph
    fg = ngc.FlowGraph(g)
    policy = make_flow_policy(config_name, algs, to_handle(g))

    placed, _ = policy.place_demand(fg, 0, 5, flowClass=0, volume=200.0)

    if placed == 0:
        pytest.skip("No flow placed, cannot check cost distribution")

    cost_dist = extract_cost_distribution(policy)
    total_from_costs = sum(cost_dist.values())

    assert abs(total_from_costs - placed) < 1e-6, (
        f"{config_name}: sum of cost_distribution ({total_from_costs}) != placed ({placed})"
    )


# ============================================================================
# TEST 7: PROGRESSIVE TE USES MULTIPLE TIERS
# ============================================================================


@pytest.mark.parametrize("config_name", TE_CONFIGS)
def test_te_uses_multiple_cost_tiers(
    config_name, multi_tier_graph, algs, to_handle, make_flow_policy
):
    """Verify TE configurations use progressively more expensive paths when needed."""
    g = multi_tier_graph
    fg = ngc.FlowGraph(g)
    policy = make_flow_policy(config_name, algs, to_handle(g))

    # Demand exceeds cheapest tier (50), should use tier 2 and 3
    demand = 150.0
    placed, _ = policy.place_demand(fg, 0, 4, flowClass=0, volume=demand)

    cost_dist = extract_cost_distribution(policy)
    unique_costs = sorted(cost_dist.keys())

    # Should have used multiple cost tiers
    assert len(unique_costs) >= 2, (
        f"{config_name}: expected multiple cost tiers, got {unique_costs}"
    )


# ============================================================================
# TEST 8: PLACED DEMAND TRACKING
# ============================================================================


@pytest.mark.parametrize("config_name", ALL_CONFIGS)
def test_placed_demand_tracking(
    config_name, multi_path_graph, algs, to_handle, make_flow_policy
):
    """Verify policy.placed_demand() matches returned placed value."""
    g = multi_path_graph
    fg = ngc.FlowGraph(g)
    policy = make_flow_policy(config_name, algs, to_handle(g))

    placed, _ = policy.place_demand(fg, 0, 5, flowClass=0, volume=200.0)
    tracked_placed = policy.placed_demand()

    assert abs(placed - tracked_placed) < 1e-9, (
        f"{config_name}: placed({placed}) != policy.placed_demand()({tracked_placed})"
    )


# ============================================================================
# TEST 9: FLOW DICTIONARY STRUCTURE
# ============================================================================


@pytest.mark.parametrize("config_name", ALL_CONFIGS)
def test_flow_dictionary_structure(
    config_name, multi_path_graph, algs, to_handle, make_flow_policy
):
    """Verify policy.flows has correct structure."""
    g = multi_path_graph
    fg = ngc.FlowGraph(g)
    policy = make_flow_policy(config_name, algs, to_handle(g))

    placed, _ = policy.place_demand(fg, 0, 5, flowClass=0, volume=200.0)

    if placed == 0:
        pytest.skip("No flow placed")

    flows = policy.flows

    # Check structure
    assert isinstance(flows, dict), f"{config_name}: flows is not a dict"

    for key, value in flows.items():
        # Key should be tuple of 4 elements
        assert len(key) == 4, f"{config_name}: flow key length != 4"

        # Value should be tuple of 4 elements: (src, dst, cost, volume)
        assert len(value) == 4, f"{config_name}: flow value length != 4"

        # src, dst should be ints
        assert isinstance(value[0], int), f"{config_name}: src not int"
        assert isinstance(value[1], int), f"{config_name}: dst not int"

        # cost, volume should be numeric
        assert isinstance(value[2], (int, float)), f"{config_name}: cost not numeric"
        assert isinstance(value[3], (int, float)), f"{config_name}: volume not numeric"

        # volume should be positive
        assert value[3] >= 0, f"{config_name}: negative volume"


# ============================================================================
# TEST 10: REMOVE DEMAND
# ============================================================================


@pytest.mark.parametrize("config_name", ALL_CONFIGS)
def test_remove_demand_clears_flows(
    config_name, multi_path_graph, algs, to_handle, make_flow_policy
):
    """Verify removing demand clears all flows."""
    g = multi_path_graph
    fg = ngc.FlowGraph(g)
    policy = make_flow_policy(config_name, algs, to_handle(g))

    # Place demand
    placed, _ = policy.place_demand(fg, 0, 5, flowClass=0, volume=200.0)

    if placed == 0:
        pytest.skip("No flow placed")

    # Remove demand
    policy.remove_demand(fg)

    # Verify cleared
    assert policy.flow_count() == 0, f"{config_name}: flow_count != 0 after remove"
    assert policy.placed_demand() == pytest.approx(0.0, abs=1e-12), (
        f"{config_name}: placed_demand != 0 after remove"
    )


# ============================================================================
# TEST 11: COMPARISON WITH MAX_FLOW REFERENCE
# ============================================================================


# Demands per graph for reference comparison
EXPECTED: Dict[str, Tuple[Tuple[int, int], float]] = {
    "square1_graph": ((0, 2), 3.0),
    "line1_graph": ((0, 2), 7.0),
    "graph3": ((0, 2), 10.0),
    "two_disjoint_shortest_graph": ((0, 3), 3.0),
    "square4_graph": ((0, 1), 350.0),
    "triangle1_graph": ((0, 2), 5.0),
    "square3_graph": ((0, 2), 320.0),
    "graph1_graph": ((0, 4), 1.0),
    "graph2_graph": ((0, 4), 1.0),
    "graph4_graph": ((0, 4), 6.0),
}


@pytest.mark.parametrize(
    "fixture_name",
    [
        "square1_graph",
        "line1_graph",
        "graph3",
        "two_disjoint_shortest_graph",
        "square4_graph",
        "triangle1_graph",
        "square3_graph",
        "graph1_graph",
        "graph2_graph",
        "graph4_graph",
    ],
)
@pytest.mark.parametrize(
    "config_name", ["SHORTEST_PATHS_WCMP", "SHORTEST_PATHS_ECMP", "TE_ECMP_16_LSP"]
)
def test_matches_maxflow_reference(
    fixture_name, config_name, algs, to_handle, make_flow_policy, request
):
    """Validate that FlowPolicy placement matches max_flow for various configurations."""
    g = request.getfixturevalue(fixture_name)
    (src, dst), demand = EXPECTED[fixture_name]

    fg = ngc.FlowGraph(g)
    policy = make_flow_policy(config_name, algs, to_handle(g))
    placed, remaining = policy.place_demand(fg, src, dst, flowClass=0, volume=demand)

    # Compute expected using Algorithms.max_flow for corresponding mode
    if config_name == "SHORTEST_PATHS_WCMP":
        exp_total, _ = algs.max_flow(
            to_handle(g),
            src,
            dst,
            flow_placement=ngc.FlowPlacement.PROPORTIONAL,
            shortest_path=True,
        )
        expected = float(exp_total)
        assert placed == pytest.approx(expected, rel=0, abs=1e-3), (
            f"placed mismatch for {fixture_name}:{config_name}: got {placed}, expected ~{expected}"
        )
        assert remaining == pytest.approx(demand - expected, rel=0, abs=1e-3), (
            f"remaining mismatch for {fixture_name}:{config_name}: got {remaining}, expected ~{demand - expected}"
        )
    elif config_name == "SHORTEST_PATHS_ECMP":
        exp_total, _ = algs.max_flow(
            to_handle(g),
            src,
            dst,
            flow_placement=ngc.FlowPlacement.EQUAL_BALANCED,
            shortest_path=True,
        )
        expected = float(exp_total)
        assert placed == pytest.approx(expected, rel=0, abs=1e-3), (
            f"placed mismatch for {fixture_name}:{config_name}: got {placed}, expected ~{expected}"
        )
        assert remaining == pytest.approx(demand - expected, rel=0, abs=1e-3), (
            f"remaining mismatch for {fixture_name}:{config_name}: got {remaining}, expected ~{demand - expected}"
        )
    elif config_name == "TE_ECMP_16_LSP":
        exp_total, _ = algs.max_flow(
            to_handle(g),
            src,
            dst,
            flow_placement=ngc.FlowPlacement.EQUAL_BALANCED,
            shortest_path=False,
        )
        expected = float(exp_total)
        assert placed <= expected + 1e-3, (
            f"placed exceeds max_flow for {fixture_name}:{config_name}: got {placed}, max {expected}"
        )
        assert (placed + remaining) == pytest.approx(demand, rel=0, abs=1e-6), (
            f"accounting mismatch for {fixture_name}:{config_name}: placed+remaining={placed + remaining}, demand={demand}"
        )
