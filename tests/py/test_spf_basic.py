"""SPF basics.

Core shortest-path semantics: distances, predecessor DAG equivalence, output
shapes/offsets, and destination early-exit behavior.
"""

import numpy as np

import netgraph_core as ngc


def run_spf(
    g,
    src,
    dst=None,
    selection=None,
    *,
    multipath=True,
    _conv=None,
    algs=None,
    to_handle=None,
):
    if selection is None:
        selection = ngc.EdgeSelection(
            multi_edge=True,
            require_capacity=False,
            tie_break=ngc.EdgeTieBreak.DETERMINISTIC,
        )
    assert algs is not None and to_handle is not None
    gh = to_handle(g)
    dist, dag = algs.spf(gh, src, dst, selection=selection, multipath=multipath)
    return dist, (_conv or (lambda gg, dd: dd))(g, dag)


def test_spf_line1_all_min_cost(
    line1_graph, dag_to_pred_map, assert_pred_dag_integrity, algs, to_handle
):
    g = line1_graph
    dist, pred = run_spf(g, 0, _conv=dag_to_pred_map, algs=algs, to_handle=to_handle)
    assert np.allclose(dist, np.array([0.0, 1.0, 2.0]))
    assert pred == {0: {}, 1: {0: [0]}, 2: {1: [1, 2]}}
    # Validate raw DAG integrity from a direct call
    dist2, dag2 = algs.spf(to_handle(g), 0)
    assert_pred_dag_integrity(g, dag2)


def test_spf_square1_all_min_cost(
    square1_graph, dag_to_pred_map, assert_pred_dag_integrity, algs, to_handle
):
    g = square1_graph
    dist, pred = run_spf(g, 0, _conv=dag_to_pred_map, algs=algs, to_handle=to_handle)
    expected_costs = np.array([0.0, 1.0, 2.0, 2.0])
    assert np.allclose(dist, expected_costs)
    assert pred == {0: {}, 1: {0: [0]}, 3: {0: [2]}, 2: {1: [1]}}
    dist2, dag2 = algs.spf(to_handle(g), 0)
    assert_pred_dag_integrity(g, dag2)


def test_spf_square2_all_min_cost(
    square2_graph, dag_to_pred_map, assert_pred_dag_integrity, algs, to_handle
):
    g = square2_graph
    dist, pred = run_spf(g, 0, _conv=dag_to_pred_map, algs=algs, to_handle=to_handle)
    assert np.isclose(dist[0], 0.0)
    assert np.isclose(dist[1], 1.0)
    assert np.isclose(dist[3], 1.0)
    assert np.isclose(dist[2], 2.0)
    assert pred == {0: {}, 1: {0: [0]}, 3: {0: [1]}, 2: {1: [2], 3: [3]}}
    dist2, dag2 = algs.spf(to_handle(g), 0)
    assert_pred_dag_integrity(g, dag2)


def test_spf_graph3_all_min_cost_and_single(
    graph3, dag_to_pred_map, assert_pred_dag_integrity, algs, to_handle
):
    g = graph3
    # ALL_MIN_COST, multipath=True
    sel = ngc.EdgeSelection(
        multi_edge=True,
        require_capacity=False,
        tie_break=ngc.EdgeTieBreak.DETERMINISTIC,
    )
    dist, pred = run_spf(
        g, 0, selection=sel, _conv=dag_to_pred_map, algs=algs, to_handle=to_handle
    )
    assert np.isclose(dist[0], 0.0)
    assert np.isclose(dist[1], 1.0)
    assert np.isclose(dist[4], 1.0)
    assert np.isclose(dist[2], 2.0)
    assert np.isclose(dist[5], 3.0)
    assert np.isclose(dist[3], 4.0)
    assert pred == {
        0: {},
        1: {0: [0, 1, 2]},
        4: {0: [3]},
        2: {1: [4, 5, 6], 4: [8]},
        5: {2: [7]},
        3: {0: [11], 2: [10], 5: [9]},
    }
    # SINGLE_MIN_COST per adjacency with ECMP across parents
    sel2 = ngc.EdgeSelection(
        multi_edge=False,
        require_capacity=False,
        tie_break=ngc.EdgeTieBreak.DETERMINISTIC,
    )
    dist2, pred2 = run_spf(
        g,
        0,
        selection=sel2,
        multipath=True,
        _conv=dag_to_pred_map,
        algs=algs,
        to_handle=to_handle,
    )
    assert np.allclose(dist2, dist)
    assert pred2 == {
        0: {},
        1: {0: [0]},
        4: {0: [3]},
        2: {1: [4], 4: [8]},
        5: {2: [7]},
        3: {0: [11], 2: [10], 5: [9]},
    }
    # Validate integrity of returned DAGs
    dist_all, dag_all = algs.spf(to_handle(g), 0, selection=sel)
    assert_pred_dag_integrity(g, dag_all)
    dist_single, dag_single = algs.spf(to_handle(g), 0, selection=sel2)
    assert_pred_dag_integrity(g, dag_single)


def test_spf_capacity_filter_all_min_cost_with_cap_remaining(
    build_graph, dag_to_pred_map, assert_pred_dag_integrity, algs, to_handle
):
    g = build_graph(2, [(0, 1, 1, 0.0, 0), (0, 1, 1, 1.0, 1)])
    sel = ngc.EdgeSelection(
        multi_edge=True, require_capacity=True, tie_break=ngc.EdgeTieBreak.DETERMINISTIC
    )
    dist, pred = run_spf(
        g, 0, selection=sel, _conv=dag_to_pred_map, algs=algs, to_handle=to_handle
    )
    assert np.allclose(dist, np.array([0.0, 1.0]))
    assert pred == {0: {}, 1: {0: [1]}}
    dist2, dag2 = algs.spf(to_handle(g), 0, selection=sel)
    assert_pred_dag_integrity(g, dag2)


def test_spf_with_dst_early_exit_equivalence(
    square2_graph, dag_to_pred_map, assert_pred_dag_integrity, algs, to_handle
):
    g = square2_graph
    dist1, pred1 = run_spf(g, 0, _conv=dag_to_pred_map, algs=algs, to_handle=to_handle)
    dist2, pred2 = run_spf(
        g, 0, dst=2, _conv=dag_to_pred_map, algs=algs, to_handle=to_handle
    )
    assert np.allclose(dist1, dist2)
    assert pred1 == pred2
    d1, dag1 = algs.spf(to_handle(g), 0)
    d2, dag2 = algs.spf(to_handle(g), 0, dst=2)
    assert_pred_dag_integrity(g, dag1)
    assert_pred_dag_integrity(g, dag2)


def test_spf_pred_dag_shapes_and_monotonic_offsets(square2_graph, algs, to_handle):
    g = square2_graph
    sel = ngc.EdgeSelection(
        multi_edge=True,
        require_capacity=False,
        tie_break=ngc.EdgeTieBreak.DETERMINISTIC,
    )
    dist, dag = algs.spf(to_handle(g), 0, 2, selection=sel)
    # Distances shape
    assert isinstance(dist, np.ndarray)
    assert dist.shape == (g.num_nodes(),)
    # PredDAG arrays
    off = np.asarray(dag.parent_offsets)
    par = np.asarray(dag.parents)
    via = np.asarray(dag.via_edges)
    # Offsets length and monotonicity
    assert off.shape == (g.num_nodes() + 1,)
    assert np.all(off[:-1] <= off[1:])
    # Parents/via sizes match total entries
    total = int(off[-1])
    assert par.shape == (total,)
    assert via.shape == (total,)
