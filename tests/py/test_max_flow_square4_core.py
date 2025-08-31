import numpy as np

import netgraph_core as ngc


def build_square4():
    # A=0, B=1, C=2, D=3
    edges = [
        (0, 1, 1, 100, 0),  # A->B
        (1, 2, 1, 125, 1),  # B->C
        (0, 3, 1, 75, 2),  # A->D
        (3, 2, 1, 50, 3),  # D->C
        (1, 3, 1, 50, 4),  # B->D
        (3, 1, 1, 50, 5),  # D->B
        (0, 1, 2, 200, 6),  # A->B (higher cost)
        (1, 3, 2, 200, 7),  # B->D (higher cost)
        (3, 2, 2, 200, 8),  # D->C (higher cost)
    ]
    src = np.array([e[0] for e in edges], dtype=np.int32)
    dst = np.array([e[1] for e in edges], dtype=np.int32)
    cost = np.array([float(e[2]) for e in edges], dtype=np.float64)
    cap = np.array([float(e[3]) for e in edges], dtype=np.float64)
    link_ids = np.array([int(e[4]) for e in edges], dtype=np.int64)
    return ngc.StrictMultiDiGraph.from_arrays(
        4, src, dst, cap, cost, link_ids, add_reverse=False
    )


def test_square4_full_and_shortest_proportional():
    g = build_square4()
    total_full, summary_full = ngc.calc_max_flow(
        g,
        0,
        1,
        flow_placement=ngc.FlowPlacement.PROPORTIONAL,
        shortest_path=False,
        eps=1e-12,
        with_edge_flows=False,
    )
    total_sp, summary_sp = ngc.calc_max_flow(
        g,
        0,
        1,
        flow_placement=ngc.FlowPlacement.PROPORTIONAL,
        shortest_path=True,
        eps=1e-12,
        with_edge_flows=False,
    )
    assert np.isclose(total_full, 350.0)
    assert np.isclose(total_sp, 100.0)
