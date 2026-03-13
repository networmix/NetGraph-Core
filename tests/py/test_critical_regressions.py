from __future__ import annotations

import numpy as np

import netgraph_core as ngc


def _graph_overflow_with_alternative():
    kmax = np.iinfo(np.int64).max
    src = np.array([0, 1, 0], dtype=np.int32)
    dst = np.array([1, 2, 2], dtype=np.int32)
    cap = np.array([1.0, 1.0, 1.0], dtype=np.float64)
    cost = np.array([kmax - 5, 10, 100], dtype=np.int64)
    ext = np.arange(3, dtype=np.int64)
    return ngc.StrictMultiDiGraph.from_arrays(3, src, dst, cap, cost, ext)


def _graph_overflow_single_path():
    kmax = np.iinfo(np.int64).max
    src = np.array([0, 1], dtype=np.int32)
    dst = np.array([1, 2], dtype=np.int32)
    cap = np.array([1.0, 1.0], dtype=np.float64)
    cost = np.array([kmax - 5, 10], dtype=np.int64)
    ext = np.arange(2, dtype=np.int64)
    return ngc.StrictMultiDiGraph.from_arrays(3, src, dst, cap, cost, ext)


def _graph_zero_cost_cycle():
    src = np.array([0, 1, 2, 2, 1], dtype=np.int32)
    dst = np.array([1, 2, 1, 3, 3], dtype=np.int32)
    cap = np.array([1.0, 1.0, 1.0, 1.0, 1.0], dtype=np.float64)
    cost = np.array([0, 0, 0, 1, 1], dtype=np.int64)
    ext = np.arange(5, dtype=np.int64)
    return ngc.StrictMultiDiGraph.from_arrays(4, src, dst, cap, cost, ext)


def test_spf_overflow_does_not_wrap_negative(algs, to_handle):
    g = _graph_overflow_with_alternative()
    gh = to_handle(g)
    dist, _ = algs.spf(gh, 0, dtype="int64")
    assert int(dist[2]) == 100


def test_ksp_overflow_saturates_finite_cost(algs, to_handle):
    g = _graph_overflow_single_path()
    gh = to_handle(g)
    items = algs.ksp(gh, 0, 2, k=1, dtype="int64")
    assert len(items) == 1
    kmax = np.iinfo(np.int64).max
    assert int(items[0][0][2]) == int(kmax - 1)


def test_max_flow_summary_costs_no_negative_wrap(algs, to_handle):
    g = _graph_overflow_single_path()
    gh = to_handle(g)
    flow, summary = algs.max_flow(gh, 0, 2)
    assert flow == 1.0
    assert int(summary.costs[0]) == int(np.iinfo(np.int64).max - 1)


def test_flow_policy_max_cost_factor_handles_large_cost(algs, to_handle):
    g = _graph_overflow_single_path()
    gh = to_handle(g)
    fg = ngc.FlowGraph(g)
    cfg = ngc.FlowPolicyConfig(max_flow_count=1, max_path_cost_factor=2.0)
    policy = ngc.FlowPolicy(algs, gh, cfg)
    placed, left = policy.place_demand(fg, 0, 2, 0, 1.0)
    assert placed == 1.0
    assert left == 0.0


def test_resolve_to_paths_terminates_on_zero_cost_cycle(algs, to_handle):
    g = _graph_zero_cost_cycle()
    gh = to_handle(g)
    _, dag = algs.spf(gh, 0, dst=3, multipath=True, dtype="int64")
    paths = dag.resolve_to_paths(0, 3)
    assert len(paths) == 2
    for path in paths:
        nodes = [int(n) for n, _ in path]
        assert len(nodes) == len(set(nodes))


def test_equal_balanced_max_flow_zero_cost_cycle(algs, to_handle):
    g = _graph_zero_cost_cycle()
    gh = to_handle(g)
    flow_prop, _ = algs.max_flow(
        gh,
        0,
        3,
        flow_placement=ngc.FlowPlacement.PROPORTIONAL,
        shortest_path=False,
        require_capacity=True,
    )
    flow_eb, _ = algs.max_flow(
        gh,
        0,
        3,
        flow_placement=ngc.FlowPlacement.EQUAL_BALANCED,
        shortest_path=False,
        require_capacity=True,
    )
    assert flow_prop == 1.0
    assert flow_eb == 1.0
