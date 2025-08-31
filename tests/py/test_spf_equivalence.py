import numpy as np

import netgraph_core as ngc


def build_graph(num_nodes, edges):
    """Build StrictMultiDiGraph from an iterable of (src, dst, cost, cap, link_id).

    - src, dst: int32
    - cost, cap: float64
    - link_id: int64 (used to match NetGraph expectations)
    """
    src = np.array([e[0] for e in edges], dtype=np.int32)
    dst = np.array([e[1] for e in edges], dtype=np.int32)
    cost = np.array([float(e[2]) for e in edges], dtype=np.float64)
    cap = np.array([float(e[3]) for e in edges], dtype=np.float64)
    link_ids = np.array([int(e[4]) for e in edges], dtype=np.int64)
    return ngc.StrictMultiDiGraph.from_arrays(
        num_nodes, src, dst, cap, cost, link_ids, add_reverse=False
    )


def dag_to_pred_map(g: ngc.StrictMultiDiGraph, dag: ngc.PredDAG):
    """Convert PredDAG to mapping: node -> {parent: [link_ids...]}, using link_id_of."""
    pred = {}
    offsets = np.asarray(dag.parent_offsets)
    parents = np.asarray(dag.parents)
    via = np.asarray(dag.via_edges)
    n = g.num_nodes()
    for v in range(n):
        start = int(offsets[v])
        end = int(offsets[v + 1])
        if start == end:
            continue
        group = {}
        for i in range(start, end):
            p = int(parents[i])
            e = int(via[i])
            lid = int(g.link_id_of(e))
            group.setdefault(p, []).append(lid)
        pred[v] = group
    # Ensure source has an empty dict if present in arrays but no parents
    if 0 not in pred and n > 0:
        pred[0] = {}
    return pred


def run_spf(
    g, src, dst=None, edge_select=ngc.EdgeSelect.ALL_MIN_COST, multipath=True, eps=1e-10
):
    dist, dag = ngc.spf(
        g, src, dst, edge_select=edge_select, multipath=multipath, eps=eps
    )
    return dist, dag_to_pred_map(g, dag)


def test_spf_line1_all_min_cost():
    # Nodes: A=0, B=1, C=2
    edges = [
        (0, 1, 1, 5, 0),  # A->B key 0
        (1, 2, 1, 1, 2),  # B->C key 2
        (1, 2, 1, 3, 4),  # B->C key 4
        (1, 2, 2, 7, 6),  # B->C key 6
    ]
    g = build_graph(3, edges)
    dist, pred = run_spf(g, 0)
    assert np.allclose(dist, np.array([0.0, 1.0, 2.0]))
    assert pred == {0: {}, 1: {0: [0]}, 2: {1: [2, 4]}}


def test_spf_square1_all_min_cost():
    # A=0, B=1, C=2, D=3
    edges = [
        (0, 1, 1, 1, 0),  # A->B key 0
        (1, 2, 1, 1, 1),  # B->C key 1
        (0, 3, 2, 2, 2),  # A->D key 2
        (3, 2, 2, 2, 3),  # D->C key 3
    ]
    g = build_graph(4, edges)
    dist, pred = run_spf(g, 0)
    expected_costs = np.array([0.0, 1.0, 2.0, 2.0])
    assert np.allclose(dist, expected_costs)
    assert pred == {0: {}, 1: {0: [0]}, 3: {0: [2]}, 2: {1: [1]}}


def test_spf_square2_all_min_cost():
    # A=0, B=1, C=2, D=3
    edges = [
        (0, 1, 1, 1, 0),  # A->B key 0
        (1, 2, 1, 1, 1),  # B->C key 1
        (0, 3, 1, 2, 2),  # A->D key 2
        (3, 2, 1, 2, 3),  # D->C key 3
    ]
    g = build_graph(4, edges)
    dist, pred = run_spf(g, 0)
    assert np.isclose(dist[0], 0.0)
    assert np.isclose(dist[1], 1.0)
    assert np.isclose(dist[3], 1.0)
    assert np.isclose(dist[2], 2.0)
    assert pred == {0: {}, 1: {0: [0]}, 3: {0: [2]}, 2: {1: [1], 3: [3]}}


def test_spf_graph3_all_min_cost_and_single():
    # A=0, B=1, C=2, D=3, E=4, F=5
    edges = [
        (0, 1, 1, 2, 0),  # A->B cost 1 (parallels)
        (0, 1, 1, 4, 1),
        (0, 1, 1, 6, 2),
        (1, 2, 1, 1, 3),  # B->C parallels
        (1, 2, 1, 2, 4),
        (1, 2, 1, 3, 5),
        (2, 3, 2, 3, 6),  # C->D cost 2
        (0, 4, 1, 5, 7),  # A->E cost 1
        (4, 2, 1, 4, 8),  # E->C cost 1
        (0, 3, 4, 2, 9),  # A->D cost 4
        (2, 5, 1, 1, 10),  # C->F cost 1
        (5, 3, 1, 2, 11),  # F->D cost 1
    ]
    g = build_graph(6, edges)
    # ALL_MIN_COST, multipath=True
    dist, pred = run_spf(g, 0, edge_select=ngc.EdgeSelect.ALL_MIN_COST, multipath=True)
    assert np.isclose(dist[0], 0.0)
    assert np.isclose(dist[1], 1.0)
    assert np.isclose(dist[4], 1.0)
    assert np.isclose(dist[2], 2.0)
    assert np.isclose(dist[5], 3.0)
    assert np.isclose(dist[3], 4.0)
    assert pred == {
        0: {},
        1: {0: [0, 1, 2]},
        4: {0: [7]},
        2: {1: [3, 4, 5], 4: [8]},
        5: {2: [10]},
        3: {0: [9], 2: [6], 5: [11]},
    }
    # SINGLE_MIN_COST, multipath=False
    dist2, pred2 = run_spf(
        g, 0, edge_select=ngc.EdgeSelect.SINGLE_MIN_COST, multipath=False
    )
    assert np.allclose(dist2, dist)
    assert pred2 == {
        0: {},
        1: {0: [0]},
        4: {0: [7]},
        2: {1: [3]},
        5: {2: [10]},
        3: {0: [9]},
    }


def test_spf_capacity_filter_all_min_cost_with_cap_remaining():
    # A=0, B=1
    edges = [
        (0, 1, 1, 0.0, 0),  # insufficient capacity
        (0, 1, 1, 1.0, 1),  # sufficient
    ]
    g = build_graph(2, edges)
    dist, pred = run_spf(
        g, 0, edge_select=ngc.EdgeSelect.ALL_MIN_COST_WITH_CAP_REMAINING
    )
    assert np.allclose(dist, np.array([0.0, 1.0]))
    assert pred == {0: {}, 1: {0: [1]}}


def test_spf_with_dst_early_exit_equivalence():
    # Use square2 and check dst early-exit produces same result
    edges = [
        (0, 1, 1, 1, 0),  # A->B key 0
        (1, 2, 1, 1, 1),  # B->C key 1
        (0, 3, 1, 2, 2),  # A->D key 2
        (3, 2, 1, 2, 3),  # D->C key 3
    ]
    g = build_graph(4, edges)
    dist1, pred1 = run_spf(g, 0)
    dist2, pred2 = run_spf(g, 0, dst=2)
    assert np.allclose(dist1, dist2)
    assert pred1 == pred2
