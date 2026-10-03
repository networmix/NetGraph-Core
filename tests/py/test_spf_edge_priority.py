"""EdgeSelection.priority: which per-edge value ranks equal-cost alternatives."""

from __future__ import annotations

import numpy as np

import netgraph_core as ngc


def _graph():
    # Two parallel equal-cost edges 0 -> 1: capacities 10 and 5, residuals 2 and 4.
    g = ngc.StrictMultiDiGraph.from_arrays(
        2,
        src=np.array([0, 0], dtype=np.int32),
        dst=np.array([1, 1], dtype=np.int32),
        capacity=np.array([10.0, 5.0], dtype=np.float64),
        cost=np.array([1, 1], dtype=np.int64),
        ext_edge_ids=np.array([100, 101], dtype=np.int64),
    )
    algs = ngc.Algorithms(ngc.Backend.cpu())
    return g, algs, algs.build_graph(g)


def _chosen_edge(algs, gh, priority, residual):
    sel = ngc.EdgeSelection(
        multi_edge=False,
        tie_break=ngc.EdgeTieBreak.PREFER_HIGHER_RESIDUAL,
        priority=priority,
    )
    _, dag = algs.spf(gh, 0, selection=sel, residual=residual)
    via = np.asarray(dag.via_edges)
    assert len(via) == 1
    return int(via[0])


def test_priority_default_and_attribute():
    sel = ngc.EdgeSelection()
    assert sel.priority == ngc.EdgePriority.RESIDUAL
    sel.priority = ngc.EdgePriority.HEADROOM
    assert sel.priority == ngc.EdgePriority.HEADROOM
    assert set(ngc.EdgePriority.__members__) == {"RESIDUAL", "CAPACITY", "HEADROOM"}


def test_priority_selects_the_ranking_value():
    _, algs, gh = _graph()
    residual = np.array([2.0, 4.0], dtype=np.float64)
    assert _chosen_edge(algs, gh, ngc.EdgePriority.RESIDUAL, residual) == 1
    assert _chosen_edge(algs, gh, ngc.EdgePriority.CAPACITY, residual) == 0
    assert _chosen_edge(algs, gh, ngc.EdgePriority.HEADROOM, residual) == 1  # 0.8 > 0.2


def test_priority_without_residual_falls_back_to_capacity():
    _, algs, gh = _graph()
    for pri in (
        ngc.EdgePriority.RESIDUAL,
        ngc.EdgePriority.CAPACITY,
        ngc.EdgePriority.HEADROOM,
    ):
        assert _chosen_edge(algs, gh, pri, None) == 0


def test_priority_does_not_change_eligibility():
    # A saturated edge stays excluded under capacity gating whatever the priority.
    _, algs, gh = _graph()
    residual = np.array([0.0, 4.0], dtype=np.float64)
    for pri in (
        ngc.EdgePriority.RESIDUAL,
        ngc.EdgePriority.CAPACITY,
        ngc.EdgePriority.HEADROOM,
    ):
        assert _chosen_edge(algs, gh, pri, residual) == 1
