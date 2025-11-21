import numpy as np

import netgraph_core as ngc


def test_sensitivity_simple():
    # S->T (cap 10)
    g = ngc.StrictMultiDiGraph.from_arrays(
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
    g = ngc.StrictMultiDiGraph.from_arrays(
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
    g = ngc.StrictMultiDiGraph.from_arrays(
        num_nodes=4,
        src=np.array([0, 0, 1, 2], dtype=np.int32),
        dst=np.array([1, 2, 3, 3], dtype=np.int32),
        capacity=np.array([10.0, 5.0, 5.0, 10.0], dtype=np.float64),
        cost=np.array([1, 1, 1, 1], dtype=np.int64),
    )
    alg = ngc.Algorithms(ngc.Backend.cpu())
    gh = alg.build_graph(g)

    # Note: shortest_path shouldn't matter for sensitivity_analysis default (which runs max flow)
    # but sensitivity_analysis API doesn't take shortest_path arg currently?
    # Let's check signature in module.cpp. It does NOT take shortest_path.
    # It hardcodes shortest_path=false in C++ implementation:
    # calc_max_flow(..., /*shortest_path=*/false, ...)

    res = alg.sensitivity_analysis(gh, 0, 3)

    # Saturated: S->B (1), A->T (2)
    saturated = {1, 2}
    assert len(res) == 2
    for eid, delta in res:
        assert eid in saturated
        assert delta == 5.0


def test_sensitivity_masked():
    # Two parallel paths, cap 10. One masked out via input mask.
    g = ngc.StrictMultiDiGraph.from_arrays(
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
