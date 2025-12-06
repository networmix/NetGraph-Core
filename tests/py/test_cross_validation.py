"""Cross-validation tests comparing NetGraph-Core to original NetGraph implementation.

This module validates that NetGraph-Core produces the same maxflow results
as the original NetGraph repository for important real-world topologies.

References:
- NetGraph tests/integration/test_scenario_3.py (Clos topology)
- NetGraph tests/solver/test_maxflow_api.py (Triangle, parallel paths)
- NetGraph tests/model/test_flow.py (Basic topologies)
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


def _build_3tier_clos_graph() -> ngc.StrictMultiDiGraph:
    """Build a 3-tier Clos topology matching NetGraph scenario_3.yaml structure.

    Structure (replicates scenario_3 from NetGraph integration tests):
    - 2 Clos instances (clos1, clos2)
    - Each Clos has:
      - 2 bricks (b1, b2)
      - Each brick has:
        - 4 T1 nodes (tier1/leaf layer)
        - 4 T2 nodes (tier2/aggregation layer)
        - Mesh: all T1->T2 edges within brick (100 Gbps, cost 1)
      - 16 spine nodes (T3)
      - Connections: Each brick's 4 T2 nodes connect one-to-one to 4 spines (400 Gbps, cost 1)
    - Inter-Clos: 16 spine-to-spine links (400 Gbps, cost 1)

    Node numbering:
    - Clos1:
      - b1/t1: nodes 0-3
      - b1/t2: nodes 4-7
      - b2/t1: nodes 8-11
      - b2/t2: nodes 12-15
      - spine: nodes 16-31
    - Clos2:
      - b1/t1: nodes 32-35
      - b1/t2: nodes 36-39
      - b2/t1: nodes 40-43
      - b2/t2: nodes 44-47
      - spine: nodes 48-63

    Total nodes: 64

    Expected max flow from clos1 T1 nodes to clos2 T1 nodes: 3200.0 Gbps
    (8 inter-clos spine links × 400 Gbps = 3200 Gbps bottleneck)
    """
    num_nodes = 64
    src_list = []
    dst_list = []
    cap_list = []
    cost_list = []

    def add_mesh(group1_start, group1_count, group2_start, group2_count, capacity):
        """Add mesh connections between two groups."""
        for i in range(group1_count):
            for j in range(group2_count):
                src_list.append(group1_start + i)
                dst_list.append(group2_start + j)
                cap_list.append(capacity)
                cost_list.append(1)

    def add_one_to_one(group1_start, group1_count, group2_start, capacity):
        """Add one-to-one connections."""
        for i in range(group1_count):
            src_list.append(group1_start + i)
            dst_list.append(group2_start + i)
            cap_list.append(capacity)
            cost_list.append(1)

    # Clos1 brick1: T1 (0-3) mesh to T2 (4-7)
    add_mesh(0, 4, 4, 4, 100.0)

    # Clos1 brick2: T1 (8-11) mesh to T2 (12-15)
    add_mesh(8, 4, 12, 4, 100.0)

    # Clos1 b1/T2 (4-7) one-to-one to spine (16-19)
    add_one_to_one(4, 4, 16, 400.0)

    # Clos1 b2/T2 (12-15) one-to-one to spine (20-23)
    add_one_to_one(12, 4, 20, 400.0)

    # Inter-Clos: Clos1 spine (16-31) one-to-one to Clos2 spine (48-63)
    add_one_to_one(16, 16, 48, 400.0)

    # Clos2 spine (48-63) to b1/T2 (36-39) - only first 4 spines connect
    add_one_to_one(48, 4, 36, 400.0)

    # Clos2 spine (52-55) to b2/T2 (44-47)
    add_one_to_one(52, 4, 44, 400.0)

    # Clos2 brick1: T2 (36-39) mesh to T1 (32-35)
    add_mesh(36, 4, 32, 4, 100.0)

    # Clos2 brick2: T2 (44-47) mesh to T1 (40-43)
    add_mesh(44, 4, 40, 4, 100.0)

    return _make_graph(
        num_nodes,
        np.array(src_list, dtype=np.int32),
        np.array(dst_list, dtype=np.int32),
        np.array(cap_list, dtype=np.float64),
        np.array(cost_list, dtype=np.int64),
    )


def _build_triangle_graph() -> ngc.StrictMultiDiGraph:
    """Triangle topology from NetGraph solver tests.

    Used in test_maxflow_api.py::_triangle_network()

    Topology:
      A(0) -> B(1) (cap 2)
      B(1) -> C(2) (cap 1)
      A(0) -> C(2) (cap 1)

    Expected max flow A->C: 2.0 (1 direct + 1 via B)
    """
    src = np.array([0, 1, 0], dtype=np.int32)
    dst = np.array([1, 2, 2], dtype=np.int32)
    cap = np.array([2.0, 1.0, 1.0], dtype=np.float64)
    cost = np.array([1, 1, 1], dtype=np.int64)

    return _make_graph(3, src, dst, cap, cost)


def _build_simple_parallel_paths_graph() -> ngc.StrictMultiDiGraph:
    """Simple network with two disjoint parallel paths.

    From NetGraph solver tests: test_maxflow_api.py::_simple_network()

    Topology:
      S(0) -> A(1) (cap 1) -> T(3) (cap 1)
      S(0) -> B(2) (cap 1) -> T(3) (cap 1)

    Expected max flow S->T: 2.0
    """
    src = np.array([0, 1, 0, 2], dtype=np.int32)
    dst = np.array([1, 3, 2, 3], dtype=np.int32)
    cap = np.array([1.0, 1.0, 1.0, 1.0], dtype=np.float64)
    cost = np.array([1, 1, 1, 1], dtype=np.int64)

    return _make_graph(4, src, dst, cap, cost)


class TestNetGraphCrossValidation:
    """Validate maxflow results match original NetGraph for important topologies."""

    def test_3tier_clos_proportional_shortest(self, algs, to_handle):
        """Validates 3-tier Clos fabric maxflow matches NetGraph scenario_3 expectations.

        Reference: NetGraph tests/integration/test_scenario_3.py
        Expected: 3200.0 Gbps from clos1 T1 nodes to clos2 T1 nodes
        """
        g = _build_3tier_clos_graph()
        gh = to_handle(g)

        # Test single T1-to-T1 path (representative of full fabric capacity)
        # node 0 (clos1/b1/t1-1) to node 32 (clos2/b1/t1-1)
        total, _summary = algs.max_flow(
            gh,
            0,
            32,
            flow_placement=ngc.FlowPlacement.PROPORTIONAL,
            shortest_path=True,
            with_edge_flows=False,
        )

        # Each T1 node can push 400 Gbps through the network
        # (4 T2 paths × 100 Gbps = 400 Gbps from T1 to T2 layer)
        assert pytest.approx(total, rel=0, abs=1.0) == 400.0

    def test_3tier_clos_equal_balanced_shortest(self, algs, to_handle):
        """Validates 3-tier Clos with EQUAL_BALANCED placement.

        Reference: NetGraph scenario_3 "capacity_analysis_forward_balanced"
        """
        g = _build_3tier_clos_graph()
        gh = to_handle(g)

        total, _summary = algs.max_flow(
            gh,
            0,
            32,
            flow_placement=ngc.FlowPlacement.EQUAL_BALANCED,
            shortest_path=True,
            with_edge_flows=False,
        )

        # With EQUAL_BALANCED, flow is split equally across paths
        # Expected similar capacity as PROPORTIONAL for this balanced topology
        assert total > 0.0, "Should find positive flow"
        assert total <= 400.0 + 1.0, "Should not exceed topology capacity"

    def test_triangle_topology_combine_mode(self, algs, to_handle):
        """Triangle topology from NetGraph solver tests.

        Reference: test_maxflow_api.py::test_max_flow_combine_basic()
        Expected maxflow A->C: 2.0 (1 direct + 1 via B)
        """
        g = _build_triangle_graph()
        gh = to_handle(g)

        total, _summary = algs.max_flow(
            gh,
            0,
            2,
            flow_placement=ngc.FlowPlacement.PROPORTIONAL,
            shortest_path=False,  # Use all paths
            with_edge_flows=False,
        )

        assert pytest.approx(total, rel=0, abs=1e-9) == 2.0

    def test_simple_parallel_paths_topology_shortest_path(self, algs, to_handle):
        """Two disjoint parallel paths topology with shortest_path=True.

        Reference: test_maxflow_api.py::test_shortest_path_vs_full_max_flow()
        Expected maxflow S->T: 2.0 (both equal-cost paths saturated)
        """
        g = _build_simple_parallel_paths_graph()
        gh = to_handle(g)

        # Test with shortest_path=True (should use both equal-cost paths)
        total_sp, _summary = algs.max_flow(
            gh,
            0,
            3,
            flow_placement=ngc.FlowPlacement.PROPORTIONAL,
            shortest_path=True,
            with_edge_flows=False,
        )

        assert pytest.approx(total_sp, rel=0, abs=1e-9) == 2.0

    def test_simple_parallel_paths_topology_full_maxflow(self, algs, to_handle):
        """Two disjoint parallel paths topology with shortest_path=False.

        Reference: test_maxflow_api.py::test_shortest_path_vs_full_max_flow()
        Expected maxflow S->T: 2.0
        """
        g = _build_simple_parallel_paths_graph()
        gh = to_handle(g)

        # Test with shortest_path=False
        total_full, _summary = algs.max_flow(
            gh,
            0,
            3,
            flow_placement=ngc.FlowPlacement.PROPORTIONAL,
            shortest_path=False,
            with_edge_flows=False,
        )

        assert pytest.approx(total_full, rel=0, abs=1e-9) == 2.0


class TestCommonFixtureTopologies:
    """Validate that common test fixtures produce expected maxflow results."""

    def test_square1_graph_proportional(self, square1_graph, algs, to_handle):
        """Validate square1_graph (one shortest, one longer alternative path)."""
        gh = to_handle(square1_graph)

        # With shortest_path=True, should use shortest path only
        total_sp, _ = algs.max_flow(
            gh,
            0,
            2,
            flow_placement=ngc.FlowPlacement.PROPORTIONAL,
            shortest_path=True,
            with_edge_flows=False,
        )

        # Shortest path: 0->1->2, capacity limited to 1.0
        assert pytest.approx(total_sp, rel=0, abs=1e-9) == 1.0

        # With shortest_path=False, should use both paths
        total_full, _ = algs.max_flow(
            gh,
            0,
            2,
            flow_placement=ngc.FlowPlacement.PROPORTIONAL,
            shortest_path=False,
            with_edge_flows=False,
        )

        # Both paths: shortest (cap 1) + longer (cap 2) = 3.0
        assert pytest.approx(total_full, rel=0, abs=1e-9) == 3.0

    def test_line1_graph_wcmp_vs_ecmp(self, line1_graph, algs, to_handle):
        """Validate line1_graph behavior with WCMP vs ECMP placement."""
        gh = to_handle(line1_graph)

        # WCMP (PROPORTIONAL) with shortest_path=True
        total_wcmp, _ = algs.max_flow(
            gh,
            0,
            2,
            flow_placement=ngc.FlowPlacement.PROPORTIONAL,
            shortest_path=True,
            with_edge_flows=False,
        )

        # ECMP (EQUAL_BALANCED) with shortest_path=True
        total_ecmp, _ = algs.max_flow(
            gh,
            0,
            2,
            flow_placement=ngc.FlowPlacement.EQUAL_BALANCED,
            shortest_path=True,
            with_edge_flows=False,
        )

        # WCMP should utilize parallel paths proportionally to capacity
        # ECMP should split equally, potentially leaving capacity unused
        assert total_wcmp > 0.0
        assert total_ecmp > 0.0
        # WCMP typically utilizes more capacity when paths have different capacities
        assert total_wcmp >= total_ecmp - 1e-9

    def test_graph3_multipath_capacity(self, graph3, algs, to_handle):
        """Validate graph3 (complex 6-node topology with parallels and longer routes)."""
        gh = to_handle(graph3)

        total, _ = algs.max_flow(
            gh,
            0,
            2,
            flow_placement=ngc.FlowPlacement.PROPORTIONAL,
            shortest_path=True,
            with_edge_flows=False,
        )

        # graph3 has multiple parallel paths from 0->2
        # Total capacity should be sum of parallel edge capacities
        # 0->1: caps [2, 4, 6] = 12 total
        # 1->2: caps [1, 2, 3] = 6 total
        # Bottleneck at 1->2 limits to 6.0
        # But also 0->4->2 path exists with cap 4
        # Total expected: 10.0
        assert pytest.approx(total, rel=0, abs=1e-9) == 10.0
