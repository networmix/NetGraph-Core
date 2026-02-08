"""KSP variants: enumeration with unique=False, k-limit, masks, and out-of-order nodes."""

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
    for d, dag in items[:3]:
        assert_pred_dag_integrity(g, dag, dist=d)


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


def test_ksp_out_of_order_node_preddag_semantics(
    build_graph, assert_pred_dag_integrity, algs, to_handle
):
    """Verify PredDAG correctness when paths visit nodes out of numerical order.

    Topology: 0->1->2 (cost 2) vs 0->3->2 (cost 4).
    Path [0,3,2] visits node 3 before node 2, which previously triggered
    a CSR fill-order bug where parent entries were swapped between nodes.
    """
    edges = [
        (0, 1, 1, 1.0, 0),  # 0->1 cost 1
        (1, 2, 1, 1.0, 1),  # 1->2 cost 1
        (0, 3, 2, 2.0, 2),  # 0->3 cost 2
        (3, 2, 2, 2.0, 3),  # 3->2 cost 2
    ]
    g = build_graph(4, edges)
    items = algs.ksp(to_handle(g), 0, 2, k=5, max_cost_factor=None, unique=True)
    assert len(items) >= 2

    # Verify costs
    costs = sorted(float(dist[2]) for dist, _ in items[:2])
    assert np.isclose(costs[0], 2.0), f"Expected shortest cost 2.0, got {costs[0]}"
    assert np.isclose(costs[1], 4.0), f"Expected second cost 4.0, got {costs[1]}"

    for d, dag in items:
        # Full semantic validation including distance consistency
        assert_pred_dag_integrity(g, dag, dist=d)

        # Reconstruct concrete paths via resolve_to_paths
        concrete = dag.resolve_to_paths(0, 2, split_parallel_edges=True)
        assert len(concrete) >= 1, "Should reconstruct at least one path"
        for path in concrete:
            assert len(path) >= 2
            assert path[0][0] == 0, "Path should start at source"
            assert path[-1][0] == 2, "Path should end at destination"
            # Verify edge connectivity
            edge_src = np.asarray(g.edge_src_view())
            edge_dst = np.asarray(g.edge_dst_view())
            for i in range(len(path) - 1):
                from_node = path[i][0]
                to_node = path[i + 1][0]
                for eid in path[i][1]:
                    assert int(edge_src[eid]) == from_node
                    assert int(edge_dst[eid]) == to_node


def test_ksp_larger_out_of_order_topology(
    build_graph, assert_pred_dag_integrity, algs, to_handle
):
    """Larger topology where k-shortest paths traverse nodes in non-numerical order.

    Path 1: 0->5->3->1 (cost 3, severely out-of-order)
    Path 2: 0->2->4->1 (cost 6, also out-of-order)
    Path 3: 0->1       (cost 10, direct)
    """
    edges = [
        (0, 5, 1, 1.0, 0),  # 0->5 cost 1
        (5, 3, 1, 1.0, 1),  # 5->3 cost 1
        (3, 1, 1, 1.0, 2),  # 3->1 cost 1
        (0, 2, 2, 1.0, 3),  # 0->2 cost 2
        (2, 4, 2, 1.0, 4),  # 2->4 cost 2
        (4, 1, 2, 1.0, 5),  # 4->1 cost 2
        (0, 1, 10, 1.0, 6),  # 0->1 cost 10
    ]
    g = build_graph(6, edges)
    items = algs.ksp(to_handle(g), 0, 1, k=5, max_cost_factor=None, unique=True)
    assert len(items) >= 2

    # Check costs non-decreasing
    costs = [float(dist[1]) for dist, _ in items]
    for i in range(len(costs) - 1):
        assert costs[i] <= costs[i + 1] + 1e-9

    # First path should be 0->5->3->1 with cost 3
    assert np.isclose(costs[0], 3.0), f"Expected shortest cost 3.0, got {costs[0]}"

    for d, dag in items:
        assert_pred_dag_integrity(g, dag, dist=d)

        concrete = dag.resolve_to_paths(0, 1, split_parallel_edges=True)
        assert len(concrete) >= 1
        for path in concrete:
            # Verify no repeated nodes (simple path)
            nodes = [step[0] for step in path]
            assert len(nodes) == len(set(nodes)), (
                f"Path contains repeated nodes: {nodes}"
            )
