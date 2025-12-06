"""Edge case and special behavior tests for NetGraph-Core.

This module tests:
- Unsafe API usage patterns that lack runtime guards
- Edge cases in EqualBalanced min-flow gating with parallel edges
- Exact integer distance calculations (int64) for large cost values

Note: Object lifetime tests are in test_lifetime_safety.py
"""

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


def test_eb_min_flow_gating_does_not_prune_valid_parallel_edges(algs):
    """EB path gating: do not prune paths where per-edge residual < min_flow but group supports it."""
    # Graph: 0->1 cap 1.2; 1->2 has two parallel edges cap 0.6 each; all costs=1
    src = np.array([0, 1, 1], dtype=np.int32)
    dst = np.array([1, 2, 2], dtype=np.int32)
    cap = np.array([1.2, 0.6, 0.6], dtype=np.float64)
    cost = np.array([1, 1, 1], dtype=np.int64)
    g = _make_graph(3, src, dst, cap, cost)
    fg = ngc.FlowGraph(g)
    sel = ngc.EdgeSelection(
        multi_edge=True, require_capacity=True, tie_break=ngc.EdgeTieBreak.DETERMINISTIC
    )
    # Use EqualBalanced with max_flow_count=1 so per_target == requested volume
    cfg = ngc.FlowPolicyConfig(
        path_alg=ngc.PathAlg.SPF,
        flow_placement=ngc.FlowPlacement.EQUAL_BALANCED,
        selection=sel,
        min_flow_count=1,
        max_flow_count=1,
        shortest_path=True,
    )
    policy = ngc.FlowPolicy(algs, algs.build_graph(g), cfg)
    placed, leftover = policy.place_demand(fg, 0, 2, 0, 1.0)
    # Should place 1.0 (group capacity min(0.6)*2 = 1.2 supports it)
    assert placed >= 1.0 - 1e-9


def test_spf_distance_dtype_int64_exact(algs):
    """Distances should support dtype='int64' to avoid float rounding for large sums."""
    # Build a chain 0->1->2 with very large costs that are not exactly representable in float64.
    big = np.int64(2**53 + 1)  # not exactly representable in float64
    src = np.array([0, 1], dtype=np.int32)
    dst = np.array([1, 2], dtype=np.int32)
    cap = np.array([1.0, 1.0], dtype=np.float64)
    cost = np.array([big, big], dtype=np.int64)
    g = _make_graph(3, src, dst, cap, cost)
    try:
        dist_i64, _ = algs.spf(algs.build_graph(g), 0, 2, dtype="int64")
    except TypeError:
        pytest.skip("Algorithms.spf does not support dtype= argument yet")
    else:
        dist_i64 = np.asarray(dist_i64)
        assert dist_i64.dtype == np.int64
        assert int(dist_i64[2]) == int(big + big)


def test_ksp_distance_dtype_int64_exact(algs):
    """KSP distances should support dtype='int64' to avoid float rounding."""
    big = np.int64(2**53 + 1)
    src = np.array([0, 1], dtype=np.int32)
    dst = np.array([1, 2], dtype=np.int32)
    cap = np.array([1.0, 1.0], dtype=np.float64)
    cost = np.array([big, big], dtype=np.int64)
    g = _make_graph(3, src, dst, cap, cost)
    try:
        items = algs.ksp(algs.build_graph(g), 0, 2, k=1, dtype="int64")
    except TypeError:
        pytest.skip("Algorithms.ksp does not support dtype= argument yet")
    else:
        assert len(items) >= 1
        dist_i64, _ = items[0]
        dist_i64 = np.asarray(dist_i64)
        assert dist_i64.dtype == np.int64
        assert int(dist_i64[2]) == int(big + big)


def test_max_flow_zero_capacity_on_shortest_path_equal_balanced(algs, to_handle):
    """EqualBalanced with require_capacity=False returns 0 when shortest path has no capacity.

    With require_capacity=False (true IP/IGP semantics), routing is cost-only.
    If the shortest path has a zero-capacity edge, no flow can be placed.

    Topology:
        Shortest path: 0->1->2 (cost 2, but 0->1 has cap=0)
        Longer path:   0->3->2 (cost 4, cap=100)
    """
    src = np.array([0, 1, 0, 3], dtype=np.int32)
    dst = np.array([1, 2, 3, 2], dtype=np.int32)
    cap = np.array(
        [0.0, 10.0, 100.0, 100.0], dtype=np.float64
    )  # 0->1 has zero capacity
    cost = np.array([1, 1, 2, 2], dtype=np.int64)
    g = _make_graph(4, src, dst, cap, cost)

    # Test with EqualBalanced + require_capacity=False + shortest_path=True
    total, summary = algs.max_flow(
        to_handle(g),
        0,
        2,
        flow_placement=ngc.FlowPlacement.EQUAL_BALANCED,
        shortest_path=True,
        require_capacity=False,
        with_edge_flows=True,
    )

    # Should return 0 because the shortest path (cost 2) has no capacity
    assert np.isclose(total, 0.0), (
        f"Expected 0 flow when shortest path has zero capacity, got {total}"
    )

    # Verify no edge flows were placed
    edge_flows = np.asarray(summary.edge_flows)
    assert np.allclose(edge_flows, 0.0), f"Expected no edge flows, got {edge_flows}"


def test_max_flow_zero_capacity_on_shortest_path_proportional(algs, to_handle):
    """Proportional with require_capacity=False returns 0 when shortest path has no capacity.

    With require_capacity=False, if the shortest path has a zero-capacity edge,
    no flow can be placed. Both placements behave consistently.
    """
    src = np.array([0, 1, 0, 3], dtype=np.int32)
    dst = np.array([1, 2, 3, 2], dtype=np.int32)
    cap = np.array([0.0, 10.0, 100.0, 100.0], dtype=np.float64)
    cost = np.array([1, 1, 2, 2], dtype=np.int64)
    g = _make_graph(4, src, dst, cap, cost)

    total, summary = algs.max_flow(
        to_handle(g),
        0,
        2,
        flow_placement=ngc.FlowPlacement.PROPORTIONAL,
        shortest_path=True,
        require_capacity=False,
        with_edge_flows=True,
    )

    assert np.isclose(total, 0.0), (
        f"Expected 0 flow when shortest path has zero capacity, got {total}"
    )


def test_max_flow_require_capacity_true_finds_alternative(algs, to_handle):
    """With require_capacity=True, max_flow should find an alternative path with capacity.

    This is the contrast test: when require_capacity=True, the algorithm should
    skip the zero-capacity edge and find the longer path that has capacity.
    """
    src = np.array([0, 1, 0, 3], dtype=np.int32)
    dst = np.array([1, 2, 3, 2], dtype=np.int32)
    cap = np.array([0.0, 10.0, 100.0, 100.0], dtype=np.float64)
    cost = np.array([1, 1, 2, 2], dtype=np.int64)
    g = _make_graph(4, src, dst, cap, cost)

    total, summary = algs.max_flow(
        to_handle(g),
        0,
        2,
        flow_placement=ngc.FlowPlacement.EQUAL_BALANCED,
        shortest_path=True,
        require_capacity=True,  # This should skip zero-capacity edges
        with_edge_flows=True,
    )

    # Should find the longer path 0->3->2 and place 100 flow
    assert np.isclose(total, 100.0), (
        f"Expected 100 flow via alternative path, got {total}"
    )
