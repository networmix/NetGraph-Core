import numpy as np

import netgraph_core as ngc


def _make_graph(num_nodes, src, dst, capacity, cost):
    """Helper to build graph with auto-generated ext_edge_ids."""
    ext_edge_ids = np.arange(len(src), dtype=np.int64)
    return ngc.StrictMultiDiGraph.from_arrays(
        num_nodes, src, dst, capacity, cost, ext_edge_ids
    )


def test_sensitivity_simple():
    # S->T (cap 10)
    g = _make_graph(
        num_nodes=2,
        src=np.array([0], dtype=np.int32),
        dst=np.array([1], dtype=np.int32),
        capacity=np.array([10.0], dtype=np.float64),
        cost=np.array([1], dtype=np.int64),
    )
    alg = ngc.Algorithms(ngc.Backend.cpu())
    gh = alg.build_graph(g)

    res = alg.sensitivity_analysis(gh, 0, 1)
    assert len(res) == 1
    assert res[0][0] == 0  # Edge ID
    assert res[0][1] == 10.0


def test_sensitivity_parallel():
    # S->T (cap 10), S->T (cap 10)
    g = _make_graph(
        num_nodes=2,
        src=np.array([0, 0], dtype=np.int32),
        dst=np.array([1, 1], dtype=np.int32),
        capacity=np.array([10.0, 10.0], dtype=np.float64),
        cost=np.array([1, 1], dtype=np.int64),
    )
    alg = ngc.Algorithms(ngc.Backend.cpu())
    gh = alg.build_graph(g)

    res = alg.sensitivity_analysis(gh, 0, 1)
    assert len(res) == 2
    for _eid, delta in res:
        assert delta == 10.0


def test_sensitivity_partial():
    # S->A (10), S->B (5), A->T (5), B->T (10)
    # 0->1 (0), 0->2 (1), 1->3 (2), 2->3 (3)
    g = _make_graph(
        num_nodes=4,
        src=np.array([0, 0, 1, 2], dtype=np.int32),
        dst=np.array([1, 2, 3, 3], dtype=np.int32),
        capacity=np.array([10.0, 5.0, 5.0, 10.0], dtype=np.float64),
        cost=np.array([1, 1, 1, 1], dtype=np.int64),
    )
    alg = ngc.Algorithms(ngc.Backend.cpu())
    gh = alg.build_graph(g)

    # With all equal-cost edges, shortest_path=True/False give the same result.

    res = alg.sensitivity_analysis(gh, 0, 3)

    # Saturated: S->B (1), A->T (2)
    saturated = {1, 2}
    assert len(res) == 2
    for eid, delta in res:
        assert eid in saturated
        assert delta == 5.0


def test_sensitivity_masked():
    # Two parallel paths, cap 10. One masked out via input mask.
    g = _make_graph(
        num_nodes=2,
        src=np.array([0, 0], dtype=np.int32),
        dst=np.array([1, 1], dtype=np.int32),
        capacity=np.array([10.0, 10.0], dtype=np.float64),
        cost=np.array([1, 1], dtype=np.int64),
    )
    alg = ngc.Algorithms(ngc.Backend.cpu())
    gh = alg.build_graph(g)

    # Mask edge 1
    edge_mask = np.array([True, False], dtype=bool)

    res = alg.sensitivity_analysis(gh, 0, 1, edge_mask=edge_mask)

    # Should only see edge 0. Max flow is 10. Sensitivity of edge 0 is 10.
    assert len(res) == 1
    assert res[0][0] == 0
    assert res[0][1] == 10.0


def test_sensitivity_shortest_path_vs_max_flow():
    """Test that shortest_path mode produces different sensitivity results.

    Topology: S(0) -> A(1) -> T(2) [cost 1+1=2, cap 10 each]
              S(0) -> B(3) -> T(2) [cost 2+2=4, cap 5 each]

    With shortest_path=False (full max-flow):
      - Uses both paths: S->A->T (10) + S->B->T (5) = 15 total
      - All 4 edges are saturated and critical

    With shortest_path=True (single-pass, IP/IGP mode):
      - Only uses cheapest path: S->A->T (10)
      - Edges 2,3 (S->B->T path) are NOT used, so NOT critical
    """
    g = _make_graph(
        num_nodes=4,
        src=np.array([0, 1, 0, 3], dtype=np.int32),
        dst=np.array([1, 2, 3, 2], dtype=np.int32),
        capacity=np.array([10.0, 10.0, 5.0, 5.0], dtype=np.float64),
        cost=np.array([1, 1, 2, 2], dtype=np.int64),  # S->A->T costs 2, S->B->T costs 4
    )
    alg = ngc.Algorithms(ngc.Backend.cpu())
    gh = alg.build_graph(g)

    # Step 1: Verify baseline flow values with max_flow
    # This validates the routing semantics before testing sensitivity
    flow_full, _ = alg.max_flow(gh, 0, 2, shortest_path=False)
    flow_sp, _ = alg.max_flow(gh, 0, 2, shortest_path=True)

    assert abs(flow_full - 15.0) < 1e-9, f"Full max-flow should be 15, got {flow_full}"
    assert abs(flow_sp - 10.0) < 1e-9, f"Shortest-path flow should be 10, got {flow_sp}"

    # Step 2: Test sensitivity with shortest_path=False (full max-flow)
    res_full = alg.sensitivity_analysis(gh, 0, 2, shortest_path=False)

    # All 4 edges should be saturated and critical
    assert len(res_full) == 4, "Full max-flow should report all 4 edges as critical"
    edge_ids_full = {eid for eid, _ in res_full}
    assert edge_ids_full == {0, 1, 2, 3}

    # Step 3: Test sensitivity with shortest_path=True (IP/IGP mode)
    res_sp = alg.sensitivity_analysis(gh, 0, 2, shortest_path=True)

    # Only edges 0,1 (S->A->T path) should be critical
    assert len(res_sp) == 2, "Shortest-path mode should only report 2 edges as critical"
    edge_ids_sp = {eid for eid, _ in res_sp}
    assert edge_ids_sp == {0, 1}, f"Expected edges 0,1 but got {edge_ids_sp}"

    # Delta values: baseline=10, with edge removed traffic uses S->B->T (cap 5)
    for eid, delta in res_sp:
        assert abs(delta - 5.0) < 1e-9, f"Edge {eid} should have delta 5.0, got {delta}"
