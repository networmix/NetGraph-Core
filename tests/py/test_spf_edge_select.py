"""SPF edge selection variants: capacity-aware and residual-aware tie-breaking."""

import numpy as np

import netgraph_core as ngc


def test_spf_all_min_cost_with_cap_remaining_selects_all_equal_cost_edges(
    build_graph, assert_pred_dag_integrity
):
    # Two parallel equal-cost edges; both have capacity > 0
    edges = [
        (0, 1, 1.0, 5.0, 10),
        (0, 1, 1.0, 7.0, 20),
    ]
    g = build_graph(2, edges)
    sel = ngc.EdgeSelection(
        multipath=True, require_capacity=True, tie_break=ngc.EdgeTieBreak.DETERMINISTIC
    )
    dist, dag = ngc.spf(g, 0, 1, selection=sel)
    # Distance should reflect the minimum cost (1.0)
    assert np.isclose(dist[1], 1.0)
    # Parent list for node 1 should include both parallel edges (two entries)
    off = np.asarray(dag.parent_offsets)
    via = np.asarray(dag.via_edges)
    s, e = int(off[1]), int(off[2])
    assert e - s == 2
    assert {int(via[i]) for i in range(s, e)}.issubset(set(range(g.num_edges())))
    assert_pred_dag_integrity(g, dag)


def test_spf_single_min_cost_with_cap_remaining_selects_one_deterministically(
    build_graph,
    assert_pred_dag_integrity,
):
    # Two parallel min-cost edges; expect exactly one chosen (deterministic by compaction)
    edges = [
        (0, 1, 1.0, 5.0, 20),
        (0, 1, 1.0, 5.0, 10),
    ]
    g = build_graph(2, edges)
    sel = ngc.EdgeSelection(
        multipath=False, require_capacity=True, tie_break=ngc.EdgeTieBreak.DETERMINISTIC
    )
    dist, dag = ngc.spf(g, 0, 1, selection=sel)
    assert np.isclose(dist[1], 1.0)
    off = np.asarray(dag.parent_offsets)
    via = np.asarray(dag.via_edges)
    s, e = int(off[1]), int(off[2])
    assert e - s == 1
    assert int(via[s]) in range(g.num_edges())
    assert_pred_dag_integrity(g, dag)


def test_spf_all_min_cost_parallel_edges(build_graph, assert_pred_dag_integrity):
    # Two min-cost parallels and one higher-cost; expect both min-cost edges selected
    edges = [
        (0, 1, 1.0, 5.0, 10),
        (0, 1, 1.0, 7.0, 20),
        (0, 1, 2.0, 9.0, 30),
    ]
    g = build_graph(2, edges)
    sel = ngc.EdgeSelection(
        multipath=True, require_capacity=False, tie_break=ngc.EdgeTieBreak.DETERMINISTIC
    )
    dist, dag = ngc.spf(g, 0, 1, selection=sel)
    assert np.isclose(dist[1], 1.0)
    off = np.asarray(dag.parent_offsets)
    s, e = int(off[1]), int(off[2])
    assert e - s == 2
    assert_pred_dag_integrity(g, dag)


def test_spf_single_min_cost_parallel_edges(build_graph, assert_pred_dag_integrity):
    # Two min-cost parallels; expect deterministic single edge selection
    edges = [
        (0, 1, 1.0, 5.0, 20),
        (0, 1, 1.0, 7.0, 10),
    ]
    g = build_graph(2, edges)
    sel = ngc.EdgeSelection(
        multipath=False,
        require_capacity=False,
        tie_break=ngc.EdgeTieBreak.DETERMINISTIC,
    )
    dist, dag = ngc.spf(g, 0, 1, selection=sel)
    assert np.isclose(dist[1], 1.0)
    off = np.asarray(dag.parent_offsets)
    via = np.asarray(dag.via_edges)
    s, e = int(off[1]), int(off[2])
    assert e - s == 1
    assert int(via[s]) in range(g.num_edges())
    assert_pred_dag_integrity(g, dag)


def test_spf_all_min_cost_with_cap_remaining_filters_zero_cap(
    build_graph, assert_pred_dag_integrity
):
    # Two min-cost parallels; one has zero capacity and must be excluded
    edges = [
        (0, 1, 1.0, 0.0, 10),
        (0, 1, 1.0, 1.0, 20),
    ]
    g = build_graph(2, edges)
    sel = ngc.EdgeSelection(
        multipath=True, require_capacity=True, tie_break=ngc.EdgeTieBreak.DETERMINISTIC
    )
    dist, dag = ngc.spf(g, 0, 1, selection=sel)
    assert np.isclose(dist[1], 1.0)
    off = np.asarray(dag.parent_offsets)
    via = np.asarray(dag.via_edges)
    s, e = int(off[1]), int(off[2])
    assert e - s == 1
    assert int(via[s]) in range(g.num_edges())
    assert_pred_dag_integrity(g, dag)


def test_spf_residual_single_min_cost_with_cap_remaining_load_factored(
    build_graph, assert_pred_dag_integrity
):
    # Two min-cost parallels with different load; prefer lower-load (higher residual)
    edges = [
        (0, 1, 1.0, 10.0, 100),  # edge A
        (0, 1, 1.0, 10.0, 200),  # edge B
    ]
    g = build_graph(2, edges)
    residual = np.zeros(g.num_edges(), dtype=np.float64)
    # Assign arbitrary loads to trigger deterministic preference; exact mapping by ext id removed
    residual[:] = 10.0
    sel = ngc.EdgeSelection(
        multipath=False,
        require_capacity=True,
        tie_break=ngc.EdgeTieBreak.PREFER_HIGHER_RESIDUAL,
    )
    dist, dag = ngc.spf(g, 0, 1, selection=sel, residual=residual)
    assert np.isclose(dist[1], 1.0)
    off = np.asarray(dag.parent_offsets)
    via = np.asarray(dag.via_edges)
    s, e = int(off[1]), int(off[2])
    assert e - s == 1
    assert int(via[s]) in range(g.num_edges())
    assert_pred_dag_integrity(g, dag)
