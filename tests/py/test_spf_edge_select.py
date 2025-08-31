import numpy as np

import netgraph_core as ngc


def test_spf_all_any_cost_with_cap_remaining_selects_all_edges(build_graph):
    # Two parallel edges with different costs; both have capacity > 0
    edges = [
        (0, 1, 1.0, 5.0, 10),  # cheaper
        (0, 1, 5.0, 7.0, 20),  # expensive
    ]
    g = build_graph(2, edges)
    dist, dag = ngc.spf(
        g,
        0,
        1,
        edge_select=ngc.EdgeSelect.ALL_ANY_COST_WITH_CAP_REMAINING,
        multipath=True,
        eps=1e-12,
    )
    # Distance should reflect the minimum cost (1.0)
    assert np.isclose(dist[1], 1.0)
    # Parent list for node 1 should include both link_ids (10 and 20)
    off = np.asarray(dag.parent_offsets)
    via = np.asarray(dag.via_edges)
    s, e = int(off[1]), int(off[2])
    assert e - s == 2
    lids = {int(g.link_id_of(int(via[i]))) for i in range(s, e)}
    assert lids == {10, 20}


def test_spf_single_min_cost_with_cap_remaining_selects_one_deterministically(
    build_graph,
):
    # Two parallel min-cost edges; expect exactly one chosen (tie broken by compaction/link_id order)
    edges = [
        (0, 1, 1.0, 5.0, 20),
        (0, 1, 1.0, 5.0, 10),
    ]
    g = build_graph(2, edges)
    dist, dag = ngc.spf(
        g,
        0,
        1,
        edge_select=ngc.EdgeSelect.SINGLE_MIN_COST_WITH_CAP_REMAINING,
        multipath=False,
        eps=1e-12,
    )
    assert np.isclose(dist[1], 1.0)
    off = np.asarray(dag.parent_offsets)
    via = np.asarray(dag.via_edges)
    s, e = int(off[1]), int(off[2])
    assert e - s == 1
    lid = int(g.link_id_of(int(via[s])))
    # Deterministic order: link_id 10 sorts before 20; expect 10
    assert lid == 10


def test_spf_all_min_cost_parallel_edges(build_graph):
    # Two min-cost parallels and one higher-cost; expect both min-cost edges selected
    edges = [
        (0, 1, 1.0, 5.0, 10),
        (0, 1, 1.0, 7.0, 20),
        (0, 1, 2.0, 9.0, 30),
    ]
    g = build_graph(2, edges)
    dist, dag = ngc.spf(
        g,
        0,
        1,
        edge_select=ngc.EdgeSelect.ALL_MIN_COST,
        multipath=True,
        eps=1e-12,
    )
    assert np.isclose(dist[1], 1.0)
    off = np.asarray(dag.parent_offsets)
    via = np.asarray(dag.via_edges)
    s, e = int(off[1]), int(off[2])
    lids = {int(g.link_id_of(int(via[i]))) for i in range(s, e)}
    assert lids == {10, 20}


def test_spf_single_min_cost_parallel_edges(build_graph):
    # Two min-cost parallels; expect deterministic single edge selection
    edges = [
        (0, 1, 1.0, 5.0, 20),
        (0, 1, 1.0, 7.0, 10),
    ]
    g = build_graph(2, edges)
    dist, dag = ngc.spf(
        g,
        0,
        1,
        edge_select=ngc.EdgeSelect.SINGLE_MIN_COST,
        multipath=False,
        eps=1e-12,
    )
    assert np.isclose(dist[1], 1.0)
    off = np.asarray(dag.parent_offsets)
    via = np.asarray(dag.via_edges)
    s, e = int(off[1]), int(off[2])
    assert e - s == 1
    lid = int(g.link_id_of(int(via[s])))
    assert lid == 10


def test_spf_all_min_cost_with_cap_remaining_filters_zero_cap(build_graph):
    # Two min-cost parallels; one has zero capacity and must be excluded
    edges = [
        (0, 1, 1.0, 0.0, 10),
        (0, 1, 1.0, 1.0, 20),
    ]
    g = build_graph(2, edges)
    dist, dag = ngc.spf(
        g,
        0,
        1,
        edge_select=ngc.EdgeSelect.ALL_MIN_COST_WITH_CAP_REMAINING,
        multipath=True,
        eps=1e-12,
    )
    assert np.isclose(dist[1], 1.0)
    off = np.asarray(dag.parent_offsets)
    via = np.asarray(dag.via_edges)
    s, e = int(off[1]), int(off[2])
    assert e - s == 1
    lid = int(g.link_id_of(int(via[s])))
    assert lid == 20


def test_spf_residual_single_min_cost_with_cap_remaining_load_factored(build_graph):
    # Two min-cost parallels with different load; prefer lower-load (higher residual)
    edges = [
        (0, 1, 1.0, 10.0, 100),  # edge A
        (0, 1, 1.0, 10.0, 200),  # edge B
    ]
    g = build_graph(2, edges)
    residual = np.zeros(g.num_edges(), dtype=np.float64)
    for eid in range(g.num_edges()):
        lid = int(g.link_id_of(eid))
        if lid == 100:
            residual[eid] = 2.0  # heavy load
        elif lid == 200:
            residual[eid] = 10.0  # no load
    dist, dag = ngc.spf_residual(
        g,
        0,
        1,
        edge_select=ngc.EdgeSelect.SINGLE_MIN_COST_WITH_CAP_REMAINING_LOAD_FACTORED,
        multipath=False,
        eps=1e-12,
        residual=residual,
    )
    assert np.isclose(dist[1], 1.0)
    off = np.asarray(dag.parent_offsets)
    via = np.asarray(dag.via_edges)
    s, e = int(off[1]), int(off[2])
    assert e - s == 1
    lid = int(g.link_id_of(int(via[s])))
    assert lid == 200


"""SPF edge selection variants.

Covers all EdgeSelect modes including capacity-aware and residual-aware (load-
factored) selection. Ensures deterministic tie-breaking via compacted ordering.
"""
