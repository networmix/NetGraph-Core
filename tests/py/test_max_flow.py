"""Max-flow end-to-end.

Validates proportional/equal-balanced placement, multi-tier augmentation,
shortest_path (single augmentation) vs full, per-edge flows, and min-cut.
"""

import numpy as np

import netgraph_core as ngc


def flows_by_eid(g: ngc.StrictMultiDiGraph, edge_flows: list[float]):
    out = {}
    for eid in range(g.num_edges()):
        out[int(eid)] = edge_flows[eid]
    return out


def test_max_flow_square1_proportional_with_edge_flows(square1_graph):
    g = square1_graph
    total, summary = ngc.max_flow(
        g,
        0,
        2,
        flow_placement=ngc.FlowPlacement.PROPORTIONAL,
        shortest_path=False,
        with_edge_flows=True,
    )
    # Multi-tier: cost 2 path carries 1, then cost 4 path carries 2 => total 3
    assert np.isclose(total, 3.0)
    fb = flows_by_eid(g, summary.edge_flows)
    # Check totals on 4 edges in this small graph (order is deterministic by compaction)
    assert len(fb) == g.num_edges()
    # Edge flows are per-edge amounts; summing across all edges counts each unit once per edge
    # along its path. For two-edge paths with total 3, the sum is 2 * 3 = 6.
    assert np.isclose(sum(summary.edge_flows), 6.0)
    # Cost distribution: {2.0: 1.0, 4.0: 2.0}
    d = {float(b.cost): float(b.share) for b in summary.cost_distribution.buckets}
    assert len(d) == 2
    assert np.isclose(d[2.0], 1.0)
    assert np.isclose(d[4.0], 2.0)


def test_square1_equal_balanced_min_cut_and_distribution(square1_graph):
    g = square1_graph
    total, summary = ngc.max_flow(
        g,
        0,
        2,
        flow_placement=ngc.FlowPlacement.EQUAL_BALANCED,
        shortest_path=False,
        with_edge_flows=True,
    )
    # Total flow should be 3 (1 along A->B->C and 2 along A->D->C)
    assert np.isclose(total, 3.0)
    # Min-cut may be edges out of A or into C; accept either set
    # Min-cut returns EdgeIds now; ensure it has size 2 and corresponds to cut around source or sink
    mc = set(map(int, summary.min_cut.edges))
    assert len(mc) == 2
    # Cost distribution by path cost: 2.0 -> 1, 4.0 -> 2
    d = {float(b.cost): float(b.share) for b in summary.cost_distribution.buckets}
    assert np.isclose(d[2.0], 1.0)
    assert np.isclose(d[4.0], 2.0)


def test_max_flow_line1_proportional_full_and_shortest(line1_graph):
    g = line1_graph
    total_full, _ = ngc.max_flow(
        g,
        0,
        2,
        flow_placement=ngc.FlowPlacement.PROPORTIONAL,
        shortest_path=False,
        with_edge_flows=False,
    )
    total_sp, _ = ngc.max_flow(
        g,
        0,
        2,
        flow_placement=ngc.FlowPlacement.PROPORTIONAL,
        shortest_path=True,
        with_edge_flows=False,
    )
    assert np.isclose(total_full, 5.0)
    assert np.isclose(total_sp, 4.0)


def test_max_flow_graph5_proportional_full_and_shortest(make_fully_connected_graph):
    # Fully connected graph with 5 nodes (A..E), capacity 1 per edge
    # Build: nodes 0..4 fully connected excluding self
    g = make_fully_connected_graph(5, cost=1.0, cap=1.0)
    total_full, _ = ngc.max_flow(
        g,
        0,
        1,
        flow_placement=ngc.FlowPlacement.PROPORTIONAL,
        shortest_path=False,
        with_edge_flows=False,
    )
    total_sp, _ = ngc.max_flow(
        g,
        0,
        1,
        flow_placement=ngc.FlowPlacement.PROPORTIONAL,
        shortest_path=True,
        with_edge_flows=False,
    )
    assert np.isclose(total_full, 4.0)
    assert np.isclose(total_sp, 1.0)


def test_max_flow_square1_shortest_path_single_augmentation(square1_graph):
    g = square1_graph
    total, summary = ngc.max_flow(
        g,
        0,
        2,
        flow_placement=ngc.FlowPlacement.PROPORTIONAL,
        shortest_path=True,
        with_edge_flows=True,
    )
    assert np.isclose(total, 1.0)
    fb = flows_by_eid(g, summary.edge_flows)
    # Should augment along A->B->C only once
    assert np.isclose(fb[0], 1.0)
    assert np.isclose(fb[1], 1.0)
    assert np.isclose(fb[2], 0.0)
    assert np.isclose(fb[3], 0.0)


def test_max_flow_line1_equal_balanced(line1_graph):
    g = line1_graph
    total, summary = ngc.max_flow(
        g,
        0,
        2,
        flow_placement=ngc.FlowPlacement.EQUAL_BALANCED,
        shortest_path=False,
        with_edge_flows=True,
    )
    # Equal-balanced across tiers: limited by A->B capacity => total 5
    assert np.isclose(total, 5.0)
    fb = flows_by_eid(g, summary.edge_flows)
    # Expect 2 across min-cost tier (1 + 1), then remaining 3 on min/higher edges by successive tiers
    assert np.isclose(fb[0], 5.0)  # A->B carries all
    # Edges are compacted deterministically; verify flows on B->C parallels
    assert np.isclose(fb[1], 1.0)
    assert np.isclose(fb[2], 3.0)
    assert np.isclose(fb[3], 1.0)


def test_max_flow_graph3_proportional_parallel_distribution(graph3):
    # A=0, B=1, C=2, D=3, E=4, F=5
    g = graph3
    total, summary = ngc.max_flow(
        g,
        0,
        2,
        flow_placement=ngc.FlowPlacement.PROPORTIONAL,
        shortest_path=False,
        with_edge_flows=True,
    )
    # Incoming to C via min-cost parents: from B (cap 1+2+3=6) and from E (cap 4) => total 10
    assert np.isclose(total, 10.0)
    fb = flows_by_eid(g, summary.edge_flows)
    # B->C parallels proportional to capacity: 1:2:3 over total 6 => 1,2,3
    assert np.isclose(fb[4], 1.0)
    assert np.isclose(fb[5], 2.0)
    assert np.isclose(fb[6], 3.0)
    # E->C should carry 4
    assert np.isclose(fb[8], 4.0)
    # A->B proportional to carry 6 (total to B): 2:4:6 over sum 12 => 1,2,3
    assert np.isclose(fb[0], 1.0)
    assert np.isclose(fb[1], 2.0)
    assert np.isclose(fb[2], 3.0)
    # A->E carries 4
    assert np.isclose(fb[8], 4.0)


def test_max_flow_two_disjoint_shortest_routes_proportional(
    two_disjoint_shortest_graph,
):
    # S=0, A=1, B=2, T=3
    # Two disjoint shortest paths S->A->T and S->B->T with equal total cost
    g = two_disjoint_shortest_graph
    total, summary = ngc.max_flow(
        g,
        0,
        3,
        flow_placement=ngc.FlowPlacement.PROPORTIONAL,
        shortest_path=False,
        with_edge_flows=True,
    )
    # Bottlenecks along S->A->T = min(3,2)=2 and S->B->T = min(4,1)=1 => total 3
    assert np.isclose(total, 3.0)
    fb = flows_by_eid(g, summary.edge_flows)
    assert np.isclose(fb[0], 2.0)  # S->A
    assert np.isclose(fb[2], 2.0)  # A->T
    assert np.isclose(fb[1], 1.0)  # S->B
    assert np.isclose(fb[3], 1.0)  # B->T
    # Cost distribution: all path costs equal (2.0) => single bucket with total
    d = {float(b.cost): float(b.share) for b in summary.cost_distribution.buckets}
    assert len(d) == 1
    assert np.isclose(d[2.0], 3.0)


def test_max_flow_square4_full_and_shortest_proportional(square4_graph):
    g = square4_graph
    total_full, _ = ngc.max_flow(
        g,
        0,
        1,
        flow_placement=ngc.FlowPlacement.PROPORTIONAL,
        shortest_path=False,
    )
    total_sp, _ = ngc.max_flow(
        g,
        0,
        1,
        flow_placement=ngc.FlowPlacement.PROPORTIONAL,
        shortest_path=True,
    )
    assert np.isclose(total_full, 350.0)
    assert np.isclose(total_sp, 100.0)
