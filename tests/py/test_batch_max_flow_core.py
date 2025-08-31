import numpy as np

import netgraph_core as ngc


def build_graph(num_nodes, edges):
    import numpy as np

    src = np.array([e[0] for e in edges], dtype=np.int32)
    dst = np.array([e[1] for e in edges], dtype=np.int32)
    cost = np.array([float(e[2]) for e in edges], dtype=np.float64)
    cap = np.array([float(e[3]) for e in edges], dtype=np.float64)
    link_ids = np.array([int(e[4]) for e in edges], dtype=np.int64)
    return ngc.StrictMultiDiGraph.from_arrays(
        num_nodes, src, dst, cap, cost, link_ids, add_reverse=False
    )


def test_batch_max_flow_matches_individuals_proportional():
    # Square1 graph: A=0, B=1, C=2, D=3
    edges = [
        (0, 1, 1, 1, 0),  # A->B
        (1, 2, 1, 1, 1),  # B->C
        (0, 3, 2, 2, 2),  # A->D
        (3, 2, 2, 2, 3),  # D->C
    ]
    g = build_graph(4, edges)
    pairs = np.array([[0, 2], [0, 3]], dtype=np.int32)
    out = ngc.batch_max_flow(
        g,
        pairs,
        flow_placement=ngc.FlowPlacement.PROPORTIONAL,
        shortest_path=False,
        eps=1e-12,
        with_edge_flows=False,
    )
    # Individual calculations
    v1, _ = ngc.calc_max_flow(
        g,
        0,
        2,
        flow_placement=ngc.FlowPlacement.PROPORTIONAL,
        shortest_path=False,
        eps=1e-12,
        with_edge_flows=False,
    )
    v2, _ = ngc.calc_max_flow(
        g,
        0,
        3,
        flow_placement=ngc.FlowPlacement.PROPORTIONAL,
        shortest_path=False,
        eps=1e-12,
        with_edge_flows=False,
    )
    assert len(out) == 2
    assert np.isclose(out[0].total_flow, v1)
    assert np.isclose(out[1].total_flow, v2)


def test_batch_max_flow_matches_individuals_equal_balanced():
    # line1 graph: A=0, B=1, C=2
    edges = [
        (0, 1, 1, 5, 0),  # A->B cap 5
        (1, 2, 1, 1, 2),  # B->C cap 1 (min-cost)
        (1, 2, 1, 3, 4),  # B->C cap 3 (min-cost)
        (1, 2, 2, 7, 6),  # B->C cap 7 (higher cost)
    ]
    g = build_graph(3, edges)
    pairs = np.array([[0, 2], [1, 2]], dtype=np.int32)
    out = ngc.batch_max_flow(
        g,
        pairs,
        flow_placement=ngc.FlowPlacement.EQUAL_BALANCED,
        shortest_path=False,
        eps=1e-12,
        with_edge_flows=False,
    )
    v1, _ = ngc.calc_max_flow(
        g,
        0,
        2,
        flow_placement=ngc.FlowPlacement.EQUAL_BALANCED,
        shortest_path=False,
        eps=1e-12,
        with_edge_flows=False,
    )
    v2, _ = ngc.calc_max_flow(
        g,
        1,
        2,
        flow_placement=ngc.FlowPlacement.EQUAL_BALANCED,
        shortest_path=False,
        eps=1e-12,
        with_edge_flows=False,
    )
    assert len(out) == 2
    assert np.isclose(out[0].total_flow, v1)
    assert np.isclose(out[1].total_flow, v2)
