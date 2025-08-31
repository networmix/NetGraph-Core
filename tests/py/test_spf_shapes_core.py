import numpy as np

import netgraph_core as ngc


def build_graph(num_nodes, edges):
    src = np.array([e[0] for e in edges], dtype=np.int32)
    dst = np.array([e[1] for e in edges], dtype=np.int32)
    cost = np.array([float(e[2]) for e in edges], dtype=np.float64)
    cap = np.array([float(e[3]) for e in edges], dtype=np.float64)
    return ngc.StrictMultiDiGraph.from_arrays(
        num_nodes, src, dst, cap, cost, link_ids=None, add_reverse=False
    )


def test_spf_pred_dag_shapes_and_monotonic_offsets():
    # square2
    edges = [
        (0, 1, 1, 1),  # A->B
        (1, 2, 1, 1),  # B->C
        (0, 3, 1, 2),  # A->D
        (3, 2, 1, 2),  # D->C
    ]
    g = build_graph(4, edges)
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
