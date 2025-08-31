"""SPF basics.

Core shortest-path semantics: distances, predecessor DAG equivalence, output
shapes/offsets, and destination early-exit behavior.
"""

import numpy as np

import netgraph_core as ngc


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


def test_spf_line1_all_min_cost(line1_graph):
    g = line1_graph
    dist, pred = run_spf(g, 0)
    assert np.allclose(dist, np.array([0.0, 1.0, 2.0]))
    assert pred == {0: {}, 1: {0: [0]}, 2: {1: [2, 4]}}


def test_spf_square1_all_min_cost(square1_graph):
    g = square1_graph
    dist, pred = run_spf(g, 0)
    expected_costs = np.array([0.0, 1.0, 2.0, 2.0])
    assert np.allclose(dist, expected_costs)
    assert pred == {0: {}, 1: {0: [0]}, 3: {0: [2]}, 2: {1: [1]}}


def test_spf_square2_all_min_cost(square2_graph):
    g = square2_graph
    dist, pred = run_spf(g, 0)
    assert np.isclose(dist[0], 0.0)
    assert np.isclose(dist[1], 1.0)
    assert np.isclose(dist[3], 1.0)
    assert np.isclose(dist[2], 2.0)
    assert pred == {0: {}, 1: {0: [0]}, 3: {0: [2]}, 2: {1: [1], 3: [3]}}


def test_spf_graph3_all_min_cost_and_single(graph3):
    g = graph3
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


def test_spf_capacity_filter_all_min_cost_with_cap_remaining(build_graph):
    g = build_graph(2, [(0, 1, 1, 0.0, 0), (0, 1, 1, 1.0, 1)])
    dist, pred = run_spf(
        g, 0, edge_select=ngc.EdgeSelect.ALL_MIN_COST_WITH_CAP_REMAINING
    )
    assert np.allclose(dist, np.array([0.0, 1.0]))
    assert pred == {0: {}, 1: {0: [1]}}


def test_spf_with_dst_early_exit_equivalence(square2_graph):
    g = square2_graph
    dist1, pred1 = run_spf(g, 0)
    dist2, pred2 = run_spf(g, 0, dst=2)
    assert np.allclose(dist1, dist2)
    assert pred1 == pred2


def test_spf_pred_dag_shapes_and_monotonic_offsets(square2_graph):
    g = square2_graph
    dist, dag = ngc.spf(
        g, 0, 2, edge_select=ngc.EdgeSelect.ALL_MIN_COST, multipath=True, eps=1e-12
    )
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
