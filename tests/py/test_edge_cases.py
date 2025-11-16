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


def test_eb_min_flow_gating_does_not_prune_valid_parallel_edges(algs):
    """EB path gating: do not prune paths where per-edge residual < min_flow but group supports it."""
    # Graph: 0->1 cap 1.2; 1->2 has two parallel edges cap 0.6 each; all costs=1
    src = np.array([0, 1, 1], dtype=np.int32)
    dst = np.array([1, 2, 2], dtype=np.int32)
    cap = np.array([1.2, 0.6, 0.6], dtype=np.float64)
    cost = np.array([1, 1, 1], dtype=np.int64)
    g = ngc.StrictMultiDiGraph.from_arrays(3, src, dst, cap, cost)
    fg = ngc.FlowGraph(g)
    sel = ngc.EdgeSelection(
        multi_edge=True, require_capacity=True, tie_break=ngc.EdgeTieBreak.DETERMINISTIC
    )
    # Use EqualBalanced with max_flow_count=1 so per_target == requested volume
    policy = ngc.FlowPolicy(
        algs,
        algs.build_graph(g),
        path_alg=ngc.PathAlg.SPF,
        flow_placement=ngc.FlowPlacement.EQUAL_BALANCED,
        selection=sel,
        min_flow_count=1,
        max_flow_count=1,
        shortest_path=True,
    )
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
    g = ngc.StrictMultiDiGraph.from_arrays(3, src, dst, cap, cost)
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
    g = ngc.StrictMultiDiGraph.from_arrays(3, src, dst, cap, cost)
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
