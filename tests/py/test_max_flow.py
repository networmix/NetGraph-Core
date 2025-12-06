"""Max-flow end-to-end: placement modes, augmentation depth, and min-cut."""

import numpy as np

import netgraph_core as ngc


def _make_graph(num_nodes, src, dst, cap, cost):
    """Helper to build graph with auto-generated ext_edge_ids."""
    ext_edge_ids = np.arange(len(src), dtype=np.int64)
    return ngc.StrictMultiDiGraph.from_arrays(
        num_nodes, src, dst, cap, cost, ext_edge_ids
    )


def test_max_flow_square1_proportional_with_edge_flows(
    square1_graph,
    flows_by_eid,
    assert_edge_flows_shape,
    assert_valid_min_cut,
    algs,
    to_handle,
):
    g = square1_graph
    total, summary = algs.max_flow(
        to_handle(g),
        0,
        2,
        flow_placement=ngc.FlowPlacement.PROPORTIONAL,
        shortest_path=False,
        with_edge_flows=True,
    )
    # Multi-tier: cost 2 path carries 1, then cost 4 path carries 2 => total 3
    assert np.isclose(total, 3.0)
    assert_edge_flows_shape(g, summary, expected_present=True)
    assert_valid_min_cut(g, summary.min_cut)
    fb = flows_by_eid(g, summary.edge_flows)
    # Check totals on 4 edges in this small graph (order is deterministic by compaction)
    assert len(fb) == g.num_edges()
    # Edge flows are per-edge amounts; summing across all edges counts each unit once per edge
    # along its path. For two-edge paths with total 3, the sum is 2 * 3 = 6.
    assert np.isclose(np.asarray(summary.edge_flows).sum(), 6.0)
    # Cost distribution checks live in test_max_flow_cost_distribution.py


def test_square1_equal_balanced_min_cut_and_distribution(
    square1_graph, assert_edge_flows_shape, assert_valid_min_cut, algs, to_handle
):
    g = square1_graph
    total, summary = algs.max_flow(
        to_handle(g),
        0,
        2,
        flow_placement=ngc.FlowPlacement.EQUAL_BALANCED,
        shortest_path=False,
        with_edge_flows=True,
    )
    # Total flow should be 3 (1 along A->B->C and 2 along A->D->C)
    assert np.isclose(total, 3.0)
    # Min-cut may be edges out of A or into C; accept either set
    # Min-cut returns EdgeIds; ensure it has size 2 and corresponds to cut around source or sink
    mc = set(map(int, summary.min_cut.edges))
    assert len(mc) == 2
    assert_valid_min_cut(g, summary.min_cut)
    assert_edge_flows_shape(g, summary, expected_present=True)
    # Cost distribution checks live in test_max_flow_cost_distribution.py


def test_max_flow_line1_proportional_full_and_shortest(line1_graph, algs, to_handle):
    g = line1_graph
    total_full, _ = algs.max_flow(
        to_handle(g),
        0,
        2,
        flow_placement=ngc.FlowPlacement.PROPORTIONAL,
        shortest_path=False,
        with_edge_flows=False,
    )
    total_sp, _ = algs.max_flow(
        to_handle(g),
        0,
        2,
        flow_placement=ngc.FlowPlacement.PROPORTIONAL,
        shortest_path=True,
        with_edge_flows=False,
    )
    assert np.isclose(total_full, 5.0)
    assert np.isclose(total_sp, 4.0)


def test_max_flow_graph5_proportional_full_and_shortest(
    make_fully_connected_graph, algs, to_handle
):
    # Fully connected graph with 5 nodes (A..E), capacity 1 per edge
    # Build: nodes 0..4 fully connected excluding self
    g = make_fully_connected_graph(5, cost=1.0, cap=1.0)
    total_full, _ = algs.max_flow(
        to_handle(g),
        0,
        1,
        flow_placement=ngc.FlowPlacement.PROPORTIONAL,
        shortest_path=False,
        with_edge_flows=False,
    )
    total_sp, _ = algs.max_flow(
        to_handle(g),
        0,
        1,
        flow_placement=ngc.FlowPlacement.PROPORTIONAL,
        shortest_path=True,
        with_edge_flows=False,
    )
    assert np.isclose(total_full, 4.0)
    assert np.isclose(total_sp, 1.0)


def test_max_flow_square1_shortest_path_single_augmentation(
    square1_graph,
    flows_by_eid,
    assert_edge_flows_shape,
    assert_valid_min_cut,
    algs,
    to_handle,
):
    g = square1_graph
    total, summary = algs.max_flow(
        to_handle(g),
        0,
        2,
        flow_placement=ngc.FlowPlacement.PROPORTIONAL,
        shortest_path=True,
        with_edge_flows=True,
    )
    assert np.isclose(total, 1.0)
    assert_edge_flows_shape(g, summary, expected_present=True)
    assert_valid_min_cut(g, summary.min_cut)
    fb = flows_by_eid(g, summary.edge_flows)
    # Should augment along A->B->C only once
    assert np.isclose(fb[0], 1.0)
    assert np.isclose(fb[1], 1.0)
    assert np.isclose(fb[2], 0.0)
    assert np.isclose(fb[3], 0.0)


def test_max_flow_line1_equal_balanced(
    line1_graph,
    flows_by_eid,
    assert_edge_flows_shape,
    assert_valid_min_cut,
    algs,
    to_handle,
):
    g = line1_graph
    total, summary = algs.max_flow(
        to_handle(g),
        0,
        2,
        flow_placement=ngc.FlowPlacement.EQUAL_BALANCED,
        shortest_path=False,
        with_edge_flows=True,
    )
    # Equal-balanced across tiers: limited by A->B capacity => total 5
    assert np.isclose(total, 5.0)
    assert_edge_flows_shape(g, summary, expected_present=True)
    assert_valid_min_cut(g, summary.min_cut)
    fb = flows_by_eid(g, summary.edge_flows)
    # Expect 2 across min-cost tier (1 + 1), then remaining 3 on min/higher edges by successive tiers
    assert np.isclose(fb[0], 5.0)  # A->B carries all
    # Edges are compacted deterministically; verify flows on B->C parallels
    assert np.isclose(fb[1], 1.0)
    assert np.isclose(fb[2], 3.0)
    assert np.isclose(fb[3], 1.0)


def test_max_flow_ecmp3_equal_balanced(algs, to_handle):
    # 0->{1,2,3}->{4}, all caps=5, equal costs
    src = np.array([0, 0, 0, 1, 2, 3], dtype=np.int32)
    dst = np.array([1, 2, 3, 4, 4, 4], dtype=np.int32)
    cap = np.array([5.0, 5.0, 5.0, 5.0, 5.0, 5.0], dtype=np.float64)
    cost = np.array([1, 1, 1, 1, 1, 1], dtype=np.int64)
    g = _make_graph(5, src, dst, cap, cost)
    total, summary = algs.max_flow(
        to_handle(g),
        0,
        4,
        flow_placement=ngc.FlowPlacement.EQUAL_BALANCED,
        shortest_path=False,
        with_edge_flows=True,
    )
    # Total is limited by downstream 3*5 = 15
    assert np.isclose(total, 15.0)
    fb = np.asarray(summary.edge_flows)
    # First hop equalized: three edges out of 0 carry 5 each
    assert np.allclose(fb[:3], [5.0, 5.0, 5.0])
    # Second hop equalized: three edges into 4 carry 5 each
    assert np.allclose(fb[3:], [5.0, 5.0, 5.0])


def test_max_flow_downstream_equal_balanced(algs, to_handle):
    # Downstream split: 0->1 (10), then 1->{2,3}(5 each), then to 4 (5 each)
    src = np.array([0, 1, 1, 2, 3], dtype=np.int32)
    dst = np.array([1, 2, 3, 4, 4], dtype=np.int32)
    cap = np.array([10.0, 5.0, 5.0, 5.0, 5.0], dtype=np.float64)
    cost = np.array([1, 1, 1, 1, 1], dtype=np.int64)
    g = _make_graph(5, src, dst, cap, cost)
    total, summary = algs.max_flow(
        to_handle(g),
        0,
        4,
        flow_placement=ngc.FlowPlacement.EQUAL_BALANCED,
        shortest_path=False,
        with_edge_flows=True,
    )
    assert np.isclose(total, 10.0)
    fb = np.asarray(summary.edge_flows)
    # First hop carries all 10
    assert np.isclose(fb[0], 10.0)
    # Downstream equal split 5/5
    assert np.isclose(fb[1], 5.0)
    assert np.isclose(fb[2], 5.0)
    assert np.isclose(fb[3], 5.0)
    assert np.isclose(fb[4], 5.0)


def test_max_flow_graph3_proportional_parallel_distribution(
    graph3, flows_by_eid, assert_edge_flows_shape, assert_valid_min_cut, algs, to_handle
):
    # A=0, B=1, C=2, D=3, E=4, F=5
    g = graph3
    total, summary = algs.max_flow(
        to_handle(g),
        0,
        2,
        flow_placement=ngc.FlowPlacement.PROPORTIONAL,
        shortest_path=False,
        with_edge_flows=True,
    )
    # Incoming to C via min-cost parents: from B (cap 1+2+3=6) and from E (cap 4) => total 10
    assert np.isclose(total, 10.0)
    assert_edge_flows_shape(g, summary, expected_present=True)
    assert_valid_min_cut(g, summary.min_cut)
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
    assert np.isclose(fb[3], 4.0)


def test_max_flow_two_disjoint_shortest_routes_proportional(
    two_disjoint_shortest_graph,
    flows_by_eid,
    assert_edge_flows_shape,
    assert_valid_min_cut,
    algs,
    to_handle,
):
    # S=0, A=1, B=2, T=3
    # Two disjoint shortest paths S->A->T and S->B->T with equal total cost
    g = two_disjoint_shortest_graph
    total, summary = algs.max_flow(
        to_handle(g),
        0,
        3,
        flow_placement=ngc.FlowPlacement.PROPORTIONAL,
        shortest_path=False,
        with_edge_flows=True,
    )
    # Bottlenecks along S->A->T = min(3,2)=2 and S->B->T = min(4,1)=1 => total 3
    assert np.isclose(total, 3.0)
    assert_edge_flows_shape(g, summary, expected_present=True)
    assert_valid_min_cut(g, summary.min_cut)
    fb = flows_by_eid(g, summary.edge_flows)
    assert np.isclose(fb[0], 2.0)  # S->A
    assert np.isclose(fb[2], 2.0)  # A->T
    assert np.isclose(fb[1], 1.0)  # S->B
    assert np.isclose(fb[3], 1.0)  # B->T
    # Cost distribution checks live in test_max_flow_cost_distribution.py


def test_max_flow_square4_full_and_shortest_proportional(
    square4_graph, algs, to_handle
):
    g = square4_graph
    total_full, _ = algs.max_flow(
        to_handle(g),
        0,
        1,
        flow_placement=ngc.FlowPlacement.PROPORTIONAL,
        shortest_path=False,
    )
    total_sp, _ = algs.max_flow(
        to_handle(g),
        0,
        1,
        flow_placement=ngc.FlowPlacement.PROPORTIONAL,
        shortest_path=True,
    )
    assert np.isclose(total_full, 350.0)
    assert np.isclose(total_sp, 100.0)
