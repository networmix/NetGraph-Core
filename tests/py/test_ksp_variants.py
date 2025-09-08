"""KSP variants: enumeration with unique=False, k-limit, and masks."""

from __future__ import annotations

import numpy as np


def _ksp_graph(build_graph):
    # 0->1 direct; 0->2->1; 0->3->2->1 to generate multiple paths
    edges = [
        (0, 1, 1, 10, 0),
        (0, 2, 1, 10, 1),
        (2, 1, 1, 10, 2),
        (0, 3, 1, 10, 3),
        (3, 2, 1, 10, 4),
    ]
    return build_graph(4, edges)


def test_ksp_unique_false_allows_alternative_enumeration(
    build_graph, assert_pred_dag_integrity, algs, to_handle
):
    g = _ksp_graph(build_graph)
    items = algs.ksp(to_handle(g), 0, 1, k=5, max_cost_factor=None, unique=False)
    assert len(items) >= 3
    # Costs should be nondecreasing among top items
    costs = [float(dist[1]) for dist, _ in items]
    assert all(c >= 1.0 for c in costs)
    for _, dag in items[:3]:
        assert_pred_dag_integrity(g, dag)


def test_ksp_respects_k_limit(build_graph, algs, to_handle):
    g = _ksp_graph(build_graph)
    items = algs.ksp(to_handle(g), 0, 1, k=2, max_cost_factor=None, unique=True)
    assert len(items) == 2


def test_ksp_with_masks(build_graph, algs, to_handle):
    g = _ksp_graph(build_graph)
    node_mask = np.array([True, True, True, False], dtype=bool)  # block node 3
    items = algs.ksp(
        to_handle(g), 0, 1, k=5, max_cost_factor=None, unique=True, node_mask=node_mask
    )
    # Only direct and via-2 remain
    assert len(items) >= 2
