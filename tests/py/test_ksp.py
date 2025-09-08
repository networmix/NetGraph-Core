"""K-shortest paths (KSP) on canonical graphs; SPF-compatible outputs."""

import numpy as np


def test_ksp_line1_two_paths(
    line1_graph, t_cost, assert_pred_dag_integrity, algs, to_handle
):
    g = line1_graph
    items = algs.ksp(to_handle(g), 0, 2, k=5, max_cost_factor=None, unique=True)
    assert len(items) >= 2
    costs = [t_cost(dist, 2) for dist, _ in items[:2]]
    assert all(np.isclose(c, 2.0) for c in costs)
    for _dist, dag in items[:2]:
        assert_pred_dag_integrity(g, dag)


def test_ksp_square1_two_paths(
    square1_graph, t_cost, assert_pred_dag_integrity, algs, to_handle
):
    g = square1_graph
    items = algs.ksp(to_handle(g), 0, 2, k=5, max_cost_factor=None, unique=True)
    assert len(items) >= 2
    costs = [t_cost(dist, 2) for dist, _ in items[:2]]
    # Expect costs 2 and 4 (A->B->C and A->D->C)
    assert set(round(c, 6) for c in costs) == {2.0, 4.0}
    for _dist, dag in items[:2]:
        assert_pred_dag_integrity(g, dag)


def test_ksp_fully_connected_costs(
    make_fully_connected_graph, t_cost, assert_pred_dag_integrity, algs, to_handle
):
    g = make_fully_connected_graph(5, cost=1.0, cap=1.0)
    items = algs.ksp(to_handle(g), 0, 1, k=2, max_cost_factor=None, unique=True)
    assert len(items) == 2
    assert np.isclose(t_cost(items[0][0], 1), 1.0)
    assert np.isclose(t_cost(items[1][0], 1), 2.0)
    for _dist, dag in items:
        assert_pred_dag_integrity(g, dag)


def test_ksp_max_cost_factor_limit(
    make_fully_connected_graph, t_cost, assert_pred_dag_integrity, algs, to_handle
):
    g = make_fully_connected_graph(5, cost=1.0, cap=1.0)
    # Restrict to paths with cost <= best_cost * 1.0 => only the direct path
    items = algs.ksp(to_handle(g), 0, 1, k=5, max_cost_factor=1.0, unique=True)
    assert len(items) == 1
    assert np.isclose(t_cost(items[0][0], 1), 1.0)
    assert_pred_dag_integrity(g, items[0][1])


def test_ksp_graph5_thresholds(
    graph5, t_cost, assert_pred_dag_integrity, algs, to_handle
):
    g = graph5
    # Best path cost = 1.0 (direct); with factor 2.0, paths up to cost 2 allowed.
    items = algs.ksp(to_handle(g), 0, 1, k=10, max_cost_factor=2.0, unique=True)
    # Ensure at least the direct and one via-node path exist
    assert len(items) >= 2
    costs = sorted(t_cost(dist, 1) for dist, _ in items[:2])
    assert np.isclose(costs[0], 1.0)
    assert np.isclose(costs[1], 2.0)
    # With factor 1.0, only direct route
    items2 = algs.ksp(to_handle(g), 0, 1, k=10, max_cost_factor=1.0, unique=True)
    assert len(items2) == 1
    assert np.isclose(t_cost(items2[0][0], 1), 1.0)
    for _dist, dag in items[:2]:
        assert_pred_dag_integrity(g, dag)
    assert_pred_dag_integrity(g, items2[0][1])


def test_ksp_square5_routes(
    square5_graph, t_cost, assert_pred_dag_integrity, algs, to_handle
):
    g = square5_graph
    # Multiple routes from A(0) -> D(3): direct two-hop and via B<->C detours
    items = algs.ksp(to_handle(g), 0, 3, k=5, max_cost_factor=None, unique=True)
    assert len(items) >= 2
    top_costs = sorted(t_cost(dist, 3) for dist, _ in items[:2])
    # Expect two shortest (cost=2)
    assert np.isclose(top_costs[0], 2.0)
    assert np.isclose(top_costs[1], 2.0)
    # No route from A(0) to E(4)
    empty = algs.ksp(to_handle(g), 0, 4, k=5, max_cost_factor=None, unique=True)
    assert empty == []
    for _dist, dag in items[:2]:
        assert_pred_dag_integrity(g, dag)
