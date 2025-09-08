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
        return ngc.FlowPolicy(
            algs,
            gh,
            path_alg=ngc.PathAlg.SPF,
            flow_placement=ngc.FlowPlacement.PROPORTIONAL,
            selection=sel,
            shortest_path=True,
        )
    if config == "SHORTEST_PATHS_ECMP":
        sel = ngc.EdgeSelection(
            multi_edge=True,
            require_capacity=True,
            tie_break=ngc.EdgeTieBreak.DETERMINISTIC,
        )
        return ngc.FlowPolicy(
            algs,
            gh,
            path_alg=ngc.PathAlg.SPF,
            flow_placement=ngc.FlowPlacement.EQUAL_BALANCED,
            selection=sel,
            shortest_path=True,
            max_flow_count=1,
        )
    if config == "TE_ECMP_16_LSP":
        sel = ngc.EdgeSelection(
            multi_edge=False,
            require_capacity=True,
            tie_break=ngc.EdgeTieBreak.PREFER_HIGHER_RESIDUAL,
        )
        return ngc.FlowPolicy(
            algs,
            gh,
            path_alg=ngc.PathAlg.SPF,
            flow_placement=ngc.FlowPlacement.EQUAL_BALANCED,
            selection=sel,
            min_flow_count=16,
            max_flow_count=16,
            reoptimize_flows_on_each_placement=True,
            # TE placement should not stop early due to tiny increments
            # to ensure full utilization
        )
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

# Expected placed/remaining per policy variant (recomputed via dev/recalc_expected_core.py)
EXPECTED_MATRIX = {
    "square1_graph": {
        "SHORTEST_PATHS_WCMP": (1.0, 2.0),
        "SHORTEST_PATHS_ECMP": (1.0, 2.0),
        "TE_ECMP_16_LSP": (2.666955787246, 0.333044212754),
    },
    "line1_graph": {
        "SHORTEST_PATHS_WCMP": (4.0, 3.0),
        "SHORTEST_PATHS_ECMP": (2.0, 5.0),
        "TE_ECMP_16_LSP": (4.0, 3.0),
    },
    "graph3": {
        "SHORTEST_PATHS_WCMP": (10.0, 0.0),
        "SHORTEST_PATHS_ECMP": (4.0, 6.0),
        "TE_ECMP_16_LSP": (9.142875552177, 0.857124447823),
    },
    "two_disjoint_shortest_graph": {
        "SHORTEST_PATHS_WCMP": (3.0, 0.0),
        "SHORTEST_PATHS_ECMP": (2.0, 1.0),
        "TE_ECMP_16_LSP": (2.90917557478, 0.0908244252205),
    },
    "square4_graph": {
        "SHORTEST_PATHS_WCMP": (100.0, 250.0),
        "SHORTEST_PATHS_ECMP": (100.0, 250.0),
        "TE_ECMP_16_LSP": (320.0002320045314, 29.999767995468574),
    },
    "triangle1_graph": {
        "SHORTEST_PATHS_WCMP": (5.0, 0.0),
        "SHORTEST_PATHS_ECMP": (5.0, 0.0),
        "TE_ECMP_16_LSP": (5.0, 0.0),
    },
    "square3_graph": {
        "SHORTEST_PATHS_WCMP": (150.0, 170.0),
        "SHORTEST_PATHS_ECMP": (100.0, 220.0),
        "TE_ECMP_16_LSP": (160.00001430511475, 159.99998569488525),
    },
    "graph1_graph": {
        "SHORTEST_PATHS_WCMP": (1.0, 0.0),
        "SHORTEST_PATHS_ECMP": (1.0, 0.0),
        "TE_ECMP_16_LSP": (1.0, 0.0),
    },
    "graph2_graph": {
        "SHORTEST_PATHS_WCMP": (1.0, 0.0),
        "SHORTEST_PATHS_ECMP": (1.0, 0.0),
        "TE_ECMP_16_LSP": (1.0, 0.0),
    },
    "graph4_graph": {
        "SHORTEST_PATHS_WCMP": (1.0, 5.0),
        "SHORTEST_PATHS_ECMP": (1.0, 5.0),
        "TE_ECMP_16_LSP": (5.333504606494898, 0.6664953935051017),
    },
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
        exp_p, exp_r = EXPECTED_MATRIX[label][key]
        # Use approx matching to avoid encoding tiny floating drift from TE placement
        assert placed == pytest.approx(exp_p, rel=0, abs=1e-3), (
            f"placed mismatch for {label}:{key}: got {placed}, expected ~{exp_p}"
        )
        assert remaining == pytest.approx(exp_r, rel=0, abs=1e-3), (
            f"remaining mismatch for {label}:{key}: got {remaining}, expected ~{exp_r}"
        )


def test_flow_policy_on_square1(square1_graph, algs, to_handle):
    g = square1_graph
    (src, dst), dem = EXPECTED["square1_graph"]
    _run_case(g, "square1_graph", src, dst, dem, algs, to_handle)


def test_flow_policy_on_line1(line1_graph, algs, to_handle):
    g = line1_graph
    (src, dst), dem = EXPECTED["line1_graph"]
    _run_case(g, "line1_graph", src, dst, dem, algs, to_handle)


def test_flow_policy_on_graph3(graph3, algs, to_handle):
    g = graph3
    (src, dst), dem = EXPECTED["graph3"]
    _run_case(g, "graph3", src, dst, dem, algs, to_handle)


def test_flow_policy_on_two_disjoint(two_disjoint_shortest_graph, algs, to_handle):
    g = two_disjoint_shortest_graph
    (src, dst), dem = EXPECTED["two_disjoint_shortest_graph"]
    _run_case(g, "two_disjoint_shortest_graph", src, dst, dem, algs, to_handle)


def test_flow_policy_on_square4(square4_graph, algs, to_handle):
    g = square4_graph
    (src, dst), dem = EXPECTED["square4_graph"]
    _run_case(g, "square4_graph", src, dst, dem, algs, to_handle)


def test_flow_policy_on_triangle1(triangle1_graph, algs, to_handle):
    g = triangle1_graph
    (src, dst), dem = EXPECTED["triangle1_graph"]
    _run_case(g, "triangle1_graph", src, dst, dem, algs, to_handle)


def test_flow_policy_on_square3(square3_graph, algs, to_handle):
    g = square3_graph
    (src, dst), dem = EXPECTED["square3_graph"]
    _run_case(g, "square3_graph", src, dst, dem, algs, to_handle)


def test_flow_policy_on_graph1(graph1_graph, algs, to_handle):
    g = graph1_graph
    (src, dst), dem = EXPECTED["graph1_graph"]
    _run_case(g, "graph1_graph", src, dst, dem, algs, to_handle)


def test_flow_policy_on_graph2(graph2_graph, algs, to_handle):
    g = graph2_graph
    (src, dst), dem = EXPECTED["graph2_graph"]
    _run_case(g, "graph2_graph", src, dst, dem, algs, to_handle)


def test_flow_policy_on_graph4(graph4_graph, algs, to_handle):
    g = graph4_graph
    (src, dst), dem = EXPECTED["graph4_graph"]
    _run_case(g, "graph4_graph", src, dst, dem, algs, to_handle)
