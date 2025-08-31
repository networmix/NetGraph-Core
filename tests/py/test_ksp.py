"""K-shortest paths (KSP) tests.

KSP must return SPF-compatible outputs: for each path, (dist: float64[N],
PredDAG). These are validated against canonical topologies by converting the
PredDAG into a predecessor mapping using link_id_of, similar to SPF tests.
"""

import numpy as np

import netgraph_core as ngc


def dag_to_pred_map(g: ngc.StrictMultiDiGraph, dag: ngc.PredDAG):
    pred = {}
    off = np.asarray(dag.parent_offsets)
    parents = np.asarray(dag.parents)
    via = np.asarray(dag.via_edges)
    for v in range(g.num_nodes()):
        s, e = int(off[v]), int(off[v + 1])
        if s == e:
            continue
        group = {}
        for i in range(s, e):
            p = int(parents[i])
            lid = int(g.link_id_of(int(via[i])))
            group.setdefault(p, []).append(lid)
        pred[v] = group
    if 0 not in pred and g.num_nodes() > 0:
        pred[0] = {}
    return pred


def t_cost(dist, t):
    return float(dist[int(t)])


def test_ksp_line1_two_paths(line1_graph):
    g = line1_graph
    items = ngc.ksp(g, 0, 2, k=5, max_cost_factor=None, unique=True)
    assert len(items) >= 2
    costs = [t_cost(dist, 2) for dist, _ in items[:2]]
    assert all(np.isclose(c, 2.0) for c in costs)


def test_ksp_square1_two_paths(square1_graph):
    g = square1_graph
    items = ngc.ksp(g, 0, 2, k=5, max_cost_factor=None, unique=True)
    assert len(items) >= 2
    costs = [t_cost(dist, 2) for dist, _ in items[:2]]
    # Expect costs 2 and 4 (A->B->C and A->D->C)
    assert set(round(c, 6) for c in costs) == {2.0, 4.0}


def test_ksp_fully_connected_costs(fully_connected_graph):
    g = fully_connected_graph(5, cost=1.0, cap=1.0)
    items = ngc.ksp(g, 0, 1, k=2, max_cost_factor=None, unique=True)
    assert len(items) == 2
    assert np.isclose(t_cost(items[0][0], 1), 1.0)
    assert np.isclose(t_cost(items[1][0], 1), 2.0)


def test_ksp_max_cost_factor_limit(fully_connected_graph):
    g = fully_connected_graph(5, cost=1.0, cap=1.0)
    # Restrict to paths with cost <= best_cost * 1.0 => only the direct path
    items = ngc.ksp(g, 0, 1, k=5, max_cost_factor=1.0, unique=True)
    assert len(items) == 1
    assert np.isclose(t_cost(items[0][0], 1), 1.0)


def test_ksp_graph5_thresholds(graph5):
    g = graph5
    # Best path cost = 1.0 (direct); with factor 2.0, paths up to cost 2 allowed.
    items = ngc.ksp(g, 0, 1, k=10, max_cost_factor=2.0, unique=True)
    # Ensure at least the direct and one via-node path exist
    assert len(items) >= 2
    costs = sorted(t_cost(dist, 1) for dist, _ in items[:2])
    assert np.isclose(costs[0], 1.0)
    assert np.isclose(costs[1], 2.0)
    # With factor 1.0, only direct route
    items2 = ngc.ksp(g, 0, 1, k=10, max_cost_factor=1.0, unique=True)
    assert len(items2) == 1
    assert np.isclose(t_cost(items2[0][0], 1), 1.0)


def test_ksp_square5_routes(square5_graph):
    g = square5_graph
    # Multiple routes from A(0) -> D(3): direct two-hop and via B<->C detours
    items = ngc.ksp(g, 0, 3, k=5, max_cost_factor=None, unique=True)
    assert len(items) >= 2
    top_costs = sorted(t_cost(dist, 3) for dist, _ in items[:2])
    # Expect two shortest (cost=2)
    assert np.isclose(top_costs[0], 2.0)
    assert np.isclose(top_costs[1], 2.0)
    # No route from A(0) to E(4)
    empty = ngc.ksp(g, 0, 4, k=5, max_cost_factor=None, unique=True)
    assert empty == []
