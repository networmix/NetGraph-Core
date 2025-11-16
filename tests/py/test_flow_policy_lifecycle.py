"""FlowPolicy lifecycle: place, rebalance, remove, and flow count bounds."""

from __future__ import annotations

import pytest

import netgraph_core as ngc


def _make_sel():
    return ngc.EdgeSelection(
        multi_edge=True,
        require_capacity=True,
        tie_break=ngc.EdgeTieBreak.DETERMINISTIC,
    )


def _small_graph(build_graph):
    # 0 -> 1 (cap 3), 0 -> 2 -> 1 (cap 2)
    edges = [
        (0, 1, 1.0, 3.0, 0),
        (0, 2, 1.0, 2.0, 1),
        (2, 1, 1.0, 2.0, 2),
    ]
    return build_graph(3, edges)


def test_flow_policy_place_rebalance_remove(build_graph, algs, to_handle):
    g = _small_graph(build_graph)
    fg = ngc.FlowGraph(g)
    cfg = ngc.FlowPolicyConfig(
        path_alg=ngc.PathAlg.SPF,
        flow_placement=ngc.FlowPlacement.EQUAL_BALANCED,
        selection=_make_sel(),
        min_flow_count=2,
        max_flow_count=4,
    )
    policy = ngc.FlowPolicy(algs, to_handle(g), cfg)

    placed, remaining = policy.place_demand(fg, 0, 1, flowClass=0, volume=4.0)
    assert placed + remaining == pytest.approx(4.0, rel=0, abs=1e-9)
    assert policy.flow_count() >= 1
    assert policy.placed_demand() == pytest.approx(placed, rel=0, abs=1e-9)
    assert isinstance(policy.flows, dict)
    assert all(len(k) == 4 for k in policy.flows.keys())

    # Rebalance to a higher target
    placed2, remaining2 = policy.rebalance_demand(fg, 0, 1, flowClass=0, target=5.0)
    # Rebalance keeps total demand constant; target is per-flow guidance
    assert placed2 + remaining2 == pytest.approx(placed, rel=0, abs=1e-9)
    # Flows should be roughly equalized after rebalance
    vals = [v[3] for v in policy.flows.values()]
    if len(vals) >= 2:
        assert max(vals) - min(vals) <= 1e-6

    # Remove and ensure counters reset
    policy.remove_demand(fg)
    assert policy.flow_count() == 0
    assert policy.placed_demand() == pytest.approx(0.0, rel=0, abs=1e-12)


def test_flow_policy_flow_count_bounds(build_graph, algs, to_handle):
    g = _small_graph(build_graph)
    fg = ngc.FlowGraph(g)
    # Constrain to exactly 3 flows
    cfg = ngc.FlowPolicyConfig(
        path_alg=ngc.PathAlg.SPF,
        flow_placement=ngc.FlowPlacement.EQUAL_BALANCED,
        selection=_make_sel(),
        min_flow_count=3,
        max_flow_count=3,
    )
    policy = ngc.FlowPolicy(algs, to_handle(g), cfg)
    placed, _ = policy.place_demand(fg, 0, 1, flowClass=0, volume=3.0)
    assert placed >= 0.0
    assert policy.flow_count() == 3
