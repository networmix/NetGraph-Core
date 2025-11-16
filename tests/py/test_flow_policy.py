"""FlowPolicy behavior across canonical graphs.

Systematic checks for PROPORTIONAL vs EQUAL_BALANCED and multipath True/False
on the canonical graphs used elsewhere. Assertions compare total placed and
remaining volumes to expected capacities.
"""

from __future__ import annotations

import math
from typing import Dict, Tuple

import pytest

import netgraph_core as ngc

# Policy configs covered (subset):
# - SHORTEST_PATHS_WCMP  -> Proportional split over equal-cost SPF (multipath)
# - SHORTEST_PATHS_ECMP  -> Equal-balanced over equal-cost SPF (multipath)
# - TE_ECMP_16_LSP       -> Equal-balanced, capacity-aware single-path selection, 16 flows
POLICY_CONFIGS = (
    "SHORTEST_PATHS_WCMP",
    "SHORTEST_PATHS_ECMP",
    "TE_ECMP_16_LSP",
)


def make_policy(config: str, algs: ngc.Algorithms, gh) -> ngc.FlowPolicy:
    if config == "SHORTEST_PATHS_WCMP":
        sel = ngc.EdgeSelection(
            multi_edge=True,
            require_capacity=True,
            tie_break=ngc.EdgeTieBreak.DETERMINISTIC,
        )
        cfg = ngc.FlowPolicyConfig(
            path_alg=ngc.PathAlg.SPF,
            flow_placement=ngc.FlowPlacement.PROPORTIONAL,
            selection=sel,
            shortest_path=True,
        )
        return ngc.FlowPolicy(algs, gh, cfg)
    if config == "SHORTEST_PATHS_ECMP":
        sel = ngc.EdgeSelection(
            multi_edge=True,
            require_capacity=True,
            tie_break=ngc.EdgeTieBreak.DETERMINISTIC,
        )
        cfg = ngc.FlowPolicyConfig(
            path_alg=ngc.PathAlg.SPF,
            flow_placement=ngc.FlowPlacement.EQUAL_BALANCED,
            selection=sel,
            shortest_path=True,
            max_flow_count=1,
        )
        return ngc.FlowPolicy(algs, gh, cfg)
    if config == "TE_ECMP_16_LSP":
        sel = ngc.EdgeSelection(
            multi_edge=False,
            require_capacity=True,
            tie_break=ngc.EdgeTieBreak.PREFER_HIGHER_RESIDUAL,
        )
        cfg = ngc.FlowPolicyConfig(
            path_alg=ngc.PathAlg.SPF,
            flow_placement=ngc.FlowPlacement.EQUAL_BALANCED,
            selection=sel,
            min_flow_count=16,
            max_flow_count=16,
            reoptimize_flows_on_each_placement=True,
            # TE placement should not stop early due to tiny increments
            # to ensure full utilization
        )
        return ngc.FlowPolicy(algs, gh, cfg)
    raise AssertionError(f"Unknown config {config}")


# Demands per graph
EXPECTED: Dict[str, Tuple[Tuple[int, int], float]] = {
    # key: graph fixture name; value: ((src,dst), demand)
    "square1_graph": ((0, 2), 3.0),  # total cap 3
    "line1_graph": ((0, 2), 7.0),  # total cap limited by 0->1 = 5
    "graph3": ((0, 2), 10.0),  # total cap 10
    "two_disjoint_shortest_graph": ((0, 3), 3.0),  # total cap 3
    "square4_graph": ((0, 1), 350.0),  # total cap 350
    "triangle1_graph": ((0, 2), 5.0),  # A->C direct cap 5, cheaper detours huge
    "square3_graph": ((0, 2), 320.0),  # large demand to exercise tiers
    "graph1_graph": ((0, 4), 1.0),
    "graph2_graph": ((0, 4), 1.0),
    "graph4_graph": ((0, 4), 6.0),
}


def _almost_equal(a: float, b: float, tol: float = 1e-9) -> bool:
    return math.isclose(a, b, rel_tol=0, abs_tol=tol)


def _run_case(
    g: ngc.StrictMultiDiGraph,
    label: str,
    src: int,
    dst: int,
    demand: float,
    algs: ngc.Algorithms,
    to_handle,
):
    for key in POLICY_CONFIGS:
        fg = ngc.FlowGraph(g)
        policy = make_policy(key, algs, to_handle(g))
        placed, remaining = policy.place_demand(
            fg, src, dst, flowClass=0, volume=demand
        )
        # Compute expected using Algorithms.max_flow for corresponding mode
        if key == "SHORTEST_PATHS_WCMP":
            exp_total, _ = algs.max_flow(
                to_handle(g),
                src,
                dst,
                flow_placement=ngc.FlowPlacement.PROPORTIONAL,
                shortest_path=True,
            )
            expected = float(exp_total)
            assert placed == pytest.approx(expected, rel=0, abs=1e-3), (
                f"placed mismatch for {label}:{key}: got {placed}, expected ~{expected}"
            )
            assert remaining == pytest.approx(demand - expected, rel=0, abs=1e-3), (
                f"remaining mismatch for {label}:{key}: got {remaining}, expected ~{demand - expected}"
            )
        elif key == "SHORTEST_PATHS_ECMP":
            exp_total, _ = algs.max_flow(
                to_handle(g),
                src,
                dst,
                flow_placement=ngc.FlowPlacement.EQUAL_BALANCED,
                shortest_path=True,
            )
            expected = float(exp_total)
            assert placed == pytest.approx(expected, rel=0, abs=1e-3), (
                f"placed mismatch for {label}:{key}: got {placed}, expected ~{expected}"
            )
            assert remaining == pytest.approx(demand - expected, rel=0, abs=1e-3), (
                f"remaining mismatch for {label}:{key}: got {remaining}, expected ~{demand - expected}"
            )
        elif key == "TE_ECMP_16_LSP":
            # TE-like full mode: placement should not exceed max_flow, and accounting must hold
            exp_total, _ = algs.max_flow(
                to_handle(g),
                src,
                dst,
                flow_placement=ngc.FlowPlacement.EQUAL_BALANCED,
                shortest_path=False,
            )
            expected = float(exp_total)
            assert placed <= expected + 1e-3, (
                f"placed exceeds max_flow for {label}:{key}: got {placed}, max {expected}"
            )
            assert (placed + remaining) == pytest.approx(demand, rel=0, abs=1e-6), (
                f"accounting mismatch for {label}:{key}: placed+remaining={placed + remaining}, demand={demand}"
            )
        else:
            raise AssertionError(f"Unknown policy key {key}")


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
def test_policy_placement_modes_match_maxflow_reference(
    fixture_name, algs, to_handle, request
):
    """Validates that FlowPolicy placement matches max_flow for various configurations.

    Tests PROPORTIONAL vs EQUAL_BALANCED vs TE_ECMP across different topologies,
    ensuring placed+remaining accounting is correct and results match max_flow reference.
    """
    g = request.getfixturevalue(fixture_name)
    (src, dst), dem = EXPECTED[fixture_name]
    _run_case(g, fixture_name, src, dst, dem, algs, to_handle)
