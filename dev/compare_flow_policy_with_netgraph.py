from __future__ import annotations

import importlib
import json
import sys
from typing import Dict, Tuple


def main() -> int:
    try:
        # Optional: import external Python reference implementation if available
        importlib.import_module("ngraph")
        from ngraph.algorithms.base import EdgeSelect, PathAlg  # type: ignore
        from ngraph.algorithms.flow_init import init_flow_graph  # type: ignore
        from ngraph.algorithms.placement import FlowPlacement  # type: ignore
        from ngraph.flows.policy import FlowPolicy  # type: ignore
        from ngraph.graph.strict_multidigraph import StrictMultiDiGraph  # type: ignore
    except Exception as exc:
        print(
            "INFO: External reference (ngraph) not available; skipping external comparison.\n"
            f"Details: {exc}"
        )
        return 0

    # Import expected matrix from this repository's tests
    try:
        tfp = importlib.import_module("tests.py.test_flow_policy")
    except Exception:
        # Fallback when not installed as package; try relative path import
        sys.path.insert(0, ".")
        tfp = importlib.import_module("tests.py.test_flow_policy")

    EXPECTED: Dict[str, Tuple[Tuple[int, int], float]] = tfp.EXPECTED  # type: ignore
    EXPECTED_MATRIX: Dict[str, Dict[str, Tuple[float, float]]] = tfp.EXPECTED_MATRIX  # type: ignore

    # Builders (mirroring tests/py/conftest.py)
    def build_triangle1() -> StrictMultiDiGraph:
        g = StrictMultiDiGraph()
        for n in (0, 1, 2):
            g.add_node(n)
        g.add_edge(0, 1, cost=1, capacity=15)
        g.add_edge(1, 0, cost=1, capacity=15)
        g.add_edge(1, 2, cost=1, capacity=15)
        g.add_edge(2, 1, cost=1, capacity=15)
        g.add_edge(0, 2, cost=1, capacity=5)
        g.add_edge(2, 0, cost=1, capacity=5)
        init_flow_graph(g)
        return g

    def build_square3() -> StrictMultiDiGraph:
        g = StrictMultiDiGraph()
        for n in (0, 1, 2, 3):
            g.add_node(n)
        g.add_edge(0, 1, cost=1, capacity=100)
        g.add_edge(1, 2, cost=1, capacity=125)
        g.add_edge(0, 3, cost=1, capacity=75)
        g.add_edge(3, 2, cost=1, capacity=50)
        g.add_edge(1, 3, cost=1, capacity=50)
        g.add_edge(3, 1, cost=1, capacity=50)
        init_flow_graph(g)
        return g

    def build_graph1() -> StrictMultiDiGraph:
        g = StrictMultiDiGraph()
        for n in range(5):
            g.add_node(n)
        g.add_edge(0, 1, cost=1, capacity=1)
        g.add_edge(0, 2, cost=1, capacity=1)
        g.add_edge(1, 3, cost=1, capacity=1)
        g.add_edge(2, 3, cost=1, capacity=1)
        g.add_edge(1, 2, cost=1, capacity=1)
        g.add_edge(2, 1, cost=1, capacity=1)
        g.add_edge(3, 4, cost=1, capacity=1)
        init_flow_graph(g)
        return g

    def build_graph2() -> StrictMultiDiGraph:
        g = StrictMultiDiGraph()
        for n in range(5):
            g.add_node(n)
        g.add_edge(0, 1, cost=1, capacity=1)
        g.add_edge(1, 2, cost=1, capacity=1)
        g.add_edge(1, 3, cost=1, capacity=1)
        g.add_edge(2, 3, cost=1, capacity=1)
        g.add_edge(3, 2, cost=1, capacity=1)
        g.add_edge(2, 4, cost=1, capacity=1)
        g.add_edge(3, 4, cost=1, capacity=1)
        init_flow_graph(g)
        return g

    def build_graph4() -> StrictMultiDiGraph:
        g = StrictMultiDiGraph()
        for n in range(5):
            g.add_node(n)
        g.add_edge(0, 1, cost=1, capacity=1)
        g.add_edge(1, 4, cost=1, capacity=1)
        g.add_edge(0, 2, cost=2, capacity=2)
        g.add_edge(2, 4, cost=2, capacity=2)
        g.add_edge(0, 3, cost=3, capacity=3)
        g.add_edge(3, 4, cost=3, capacity=3)
        init_flow_graph(g)
        return g

    def build_square1() -> StrictMultiDiGraph:
        g = StrictMultiDiGraph()
        for n in range(4):
            g.add_node(n)
        g.add_edge(0, 1, cost=1, capacity=1)
        g.add_edge(1, 2, cost=1, capacity=1)
        g.add_edge(0, 3, cost=2, capacity=2)
        g.add_edge(3, 2, cost=2, capacity=2)
        init_flow_graph(g)
        return g

    def build_line1() -> StrictMultiDiGraph:
        g = StrictMultiDiGraph()
        for n in range(3):
            g.add_node(n)
        g.add_edge(0, 1, cost=1, capacity=5)
        g.add_edge(1, 2, cost=1, capacity=1)
        g.add_edge(1, 2, cost=1, capacity=3)
        g.add_edge(1, 2, cost=2, capacity=7)
        init_flow_graph(g)
        return g

    def build_graph3() -> StrictMultiDiGraph:
        g = StrictMultiDiGraph()
        for n in range(6):
            g.add_node(n)
        g.add_edge(0, 1, cost=1, capacity=2)
        g.add_edge(0, 1, cost=1, capacity=4)
        g.add_edge(0, 1, cost=1, capacity=6)
        g.add_edge(1, 2, cost=1, capacity=1)
        g.add_edge(1, 2, cost=1, capacity=2)
        g.add_edge(1, 2, cost=1, capacity=3)
        g.add_edge(2, 3, cost=2, capacity=3)
        g.add_edge(0, 4, cost=1, capacity=5)
        g.add_edge(4, 2, cost=1, capacity=4)
        g.add_edge(0, 3, cost=4, capacity=2)
        g.add_edge(2, 5, cost=1, capacity=1)
        g.add_edge(5, 3, cost=1, capacity=2)
        init_flow_graph(g)
        return g

    def build_two_disjoint() -> StrictMultiDiGraph:
        g = StrictMultiDiGraph()
        for n in range(4):
            g.add_node(n)
        g.add_edge(0, 1, cost=1, capacity=3)
        g.add_edge(1, 3, cost=1, capacity=2)
        g.add_edge(0, 2, cost=1, capacity=4)
        g.add_edge(2, 3, cost=1, capacity=1)
        init_flow_graph(g)
        return g

    def build_square4() -> StrictMultiDiGraph:
        g = StrictMultiDiGraph()
        for n in range(4):
            g.add_node(n)
        g.add_edge(0, 1, cost=1, capacity=100)
        g.add_edge(1, 2, cost=1, capacity=125)
        g.add_edge(0, 3, cost=1, capacity=75)
        g.add_edge(3, 2, cost=1, capacity=50)
        g.add_edge(1, 3, cost=1, capacity=50)
        g.add_edge(3, 1, cost=1, capacity=50)
        g.add_edge(0, 1, cost=2, capacity=200)
        g.add_edge(1, 3, cost=2, capacity=200)
        g.add_edge(3, 2, cost=2, capacity=200)
        init_flow_graph(g)
        return g

    BUILDERS = {
        "square1_graph": build_square1,
        "line1_graph": build_line1,
        "graph3": build_graph3,
        "two_disjoint_shortest_graph": build_two_disjoint,
        "square4_graph": build_square4,
        "triangle1_graph": build_triangle1,
        "square3_graph": build_square3,
        "graph1_graph": build_graph1,
        "graph2_graph": build_graph2,
        "graph4_graph": build_graph4,
    }

    POLICIES = [
        (
            "PROPORTIONAL_ALL_MIN_COST_mpTrue",
            FlowPlacement.PROPORTIONAL,
            EdgeSelect.ALL_MIN_COST,
            True,
            None,
        ),
        (
            "PROPORTIONAL_SINGLE_MIN_COST_WITH_CAP_mpFalse",
            FlowPlacement.PROPORTIONAL,
            EdgeSelect.SINGLE_MIN_COST_WITH_CAP_REMAINING,
            False,
            None,
        ),
        (
            "EQUAL_BALANCED_ALL_MIN_COST_mpTrue",
            FlowPlacement.EQUAL_BALANCED,
            EdgeSelect.ALL_MIN_COST,
            True,
            16,
        ),
        (
            "EQUAL_BALANCED_SINGLE_MIN_COST_WITH_CAP_mpFalse",
            FlowPlacement.EQUAL_BALANCED,
            EdgeSelect.SINGLE_MIN_COST_WITH_CAP_REMAINING,
            False,
            16,
        ),
    ]

    mismatches = []
    for gname, (sd, demand) in EXPECTED.items():
        src, dst = sd
        builder = BUILDERS[gname]
        g = builder()
        for key, placement, sel, mp, mfc in POLICIES:
            policy = FlowPolicy(
                path_alg=PathAlg.SPF,
                flow_placement=placement,
                edge_select=sel,
                multipath=mp,
                max_flow_count=mfc,
            )
            placed, remaining = policy.place_demand(g, src, dst, "cls", demand)
            exp_p, exp_r = EXPECTED_MATRIX[gname][key]
            if abs(placed - exp_p) > 1e-6 or abs(remaining - exp_r) > 1e-6:
                mismatches.append(
                    {
                        "graph": gname,
                        "policy": key,
                        "expected": (exp_p, exp_r),
                        "actual": (placed, remaining),
                    }
                )

    print(json.dumps({"mismatches": mismatches}, indent=2))
    return 0 if not mismatches else 1


if __name__ == "__main__":
    raise SystemExit(main())
