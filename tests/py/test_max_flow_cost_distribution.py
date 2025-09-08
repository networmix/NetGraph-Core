"""Max-flow cost distribution analytics across scenarios."""

import numpy as np

import netgraph_core as ngc


def test_cost_distribution_multiple_paths(build_graph, cost_dist_dict, algs, to_handle):
    # S=0, A=1, B=2, T=3
    # Path1: S->A->T cost=2, cap=5; Path2: S->B->T cost=4, cap=3
    edges = [
        (0, 1, 1, 5, 0),
        (1, 3, 1, 5, 1),
        (0, 2, 2, 3, 2),
        (2, 3, 2, 3, 3),
    ]
    g = build_graph(4, edges)
    total, summary = algs.max_flow(
        to_handle(g),
        0,
        3,
        flow_placement=ngc.FlowPlacement.PROPORTIONAL,
        shortest_path=False,
        with_edge_flows=False,
    )
    assert np.isclose(total, 8.0)
    d = cost_dist_dict(summary)
    assert np.isclose(d[2.0], 5.0)
    assert np.isclose(d[4.0], 3.0)


def test_cost_distribution_single_path(build_graph, cost_dist_dict, algs, to_handle):
    # A=0, B=1, C=2; cost=3+2=5; cap=10
    edges = [
        (0, 1, 3, 10, 0),
        (1, 2, 2, 10, 1),
    ]
    g = build_graph(3, edges)
    total, summary = algs.max_flow(
        to_handle(g),
        0,
        2,
        flow_placement=ngc.FlowPlacement.PROPORTIONAL,
        shortest_path=False,
        with_edge_flows=False,
    )
    assert np.isclose(total, 10.0)
    d = cost_dist_dict(summary)
    assert np.isclose(d[5.0], 10.0)


def test_cost_distribution_equal_cost_paths(
    build_graph, cost_dist_dict, algs, to_handle
):
    # S=0, A=1, B=2, T=3; two paths both cost=2, caps 4 and 6
    edges = [
        (0, 1, 1, 4, 0),
        (1, 3, 1, 4, 1),
        (0, 2, 1, 6, 2),
        (2, 3, 1, 6, 3),
    ]
    g = build_graph(4, edges)
    total, summary = algs.max_flow(
        to_handle(g),
        0,
        3,
        flow_placement=ngc.FlowPlacement.PROPORTIONAL,
        shortest_path=False,
        with_edge_flows=False,
    )
    assert np.isclose(total, 10.0)
    d = cost_dist_dict(summary)
    assert np.isclose(d[2.0], 10.0)


def test_cost_distribution_three_tiers(build_graph, cost_dist_dict, algs, to_handle):
    # S=0, A=1, B=2, C=3, T=4
    # Tier1: cost=1 cap=2 (0->1 cost1, 1->4 cost0)
    # Tier2: cost=3 cap=4 (0->2 cost2, 2->4 cost1)
    # Tier3: cost=6 cap=3 (0->3 cost3, 3->4 cost3)
    edges = [
        (0, 1, 1, 2, 0),
        (1, 4, 0, 2, 1),
        (0, 2, 2, 4, 2),
        (2, 4, 1, 4, 3),
        (0, 3, 3, 3, 4),
        (3, 4, 3, 3, 5),
    ]
    g = build_graph(5, edges)
    total, summary = algs.max_flow(
        to_handle(g),
        0,
        4,
        flow_placement=ngc.FlowPlacement.PROPORTIONAL,
        shortest_path=False,
        with_edge_flows=False,
    )
    assert np.isclose(total, 9.0)
    d = cost_dist_dict(summary)
    assert np.isclose(d[1.0], 2.0)
    assert np.isclose(d[3.0], 4.0)
    assert np.isclose(d[6.0], 3.0)


def test_cost_distribution_no_flow(build_graph, cost_dist_dict, algs, to_handle):
    g = build_graph(2, [])  # No edges
    total, summary = algs.max_flow(
        to_handle(g),
        0,
        1,
        flow_placement=ngc.FlowPlacement.PROPORTIONAL,
        shortest_path=False,
        with_edge_flows=False,
    )
    assert np.isclose(total, 0.0)
    assert np.asarray(summary.costs).size == 0 and np.asarray(summary.flows).size == 0


def test_cost_distribution_shortest_path_mode(
    build_graph, cost_dist_dict, algs, to_handle
):
    # Two paths: cost=2 cap=5; cost=4 cap=3; shortest_path=True -> only first tier
    edges = [
        (0, 1, 1, 5, 0),
        (1, 3, 1, 5, 1),
        (0, 2, 2, 3, 2),
        (2, 3, 2, 3, 3),
    ]
    g = build_graph(4, edges)
    total, summary = algs.max_flow(
        to_handle(g),
        0,
        3,
        flow_placement=ngc.FlowPlacement.PROPORTIONAL,
        shortest_path=True,
        with_edge_flows=False,
    )
    assert np.isclose(total, 5.0)
    d = cost_dist_dict(summary)
    assert len(d) == 1
    assert np.isclose(d[2.0], 5.0)


def test_cost_distribution_capacity_bottleneck(
    build_graph, cost_dist_dict, algs, to_handle
):
    # Cheap path bottlenecked at 2 (cost 1); expensive path cost 3 cap 5
    edges = [
        (0, 1, 1, 10, 0),
        (1, 3, 0, 2, 1),  # bottleneck at 2
        (0, 2, 2, 5, 2),
        (2, 3, 1, 5, 3),
    ]
    g = build_graph(4, edges)
    total, summary = algs.max_flow(
        to_handle(g),
        0,
        3,
        flow_placement=ngc.FlowPlacement.PROPORTIONAL,
        shortest_path=False,
        with_edge_flows=False,
    )
    assert np.isclose(total, 7.0)
    d = cost_dist_dict(summary)
    assert np.isclose(d[1.0], 2.0)
    assert np.isclose(d[3.0], 5.0)
