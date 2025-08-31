import numpy as np

import netgraph_core as ngc


def build_graph(num_nodes, edges):
    src = np.array([e[0] for e in edges], dtype=np.int32)
    dst = np.array([e[1] for e in edges], dtype=np.int32)
    cost = np.array([float(e[2]) for e in edges], dtype=np.float64)
    cap = np.array([float(e[3]) for e in edges], dtype=np.float64)
    link_ids = np.array([int(e[4]) for e in edges], dtype=np.int64)
    return ngc.StrictMultiDiGraph.from_arrays(
        num_nodes, src, dst, cap, cost, link_ids, add_reverse=False
    )


def test_node_mask_blocks_path():
    # A=0 -> B=1 -> C=2, all cost=1
    edges = [
        (0, 1, 1, 1, 0),
        (1, 2, 1, 1, 1),
    ]
    g = build_graph(3, edges)
    node_mask = np.array([True, False, True], dtype=bool)
    dist, dag = ngc.spf(
        g,
        0,
        2,
        edge_select=ngc.EdgeSelect.ALL_MIN_COST,
        multipath=True,
        node_mask=node_mask,
        eps=1e-12,
    )
    # Destination unreachable because B is masked out
    assert np.isinf(dist[2])
    # No parents recorded for 2
    off = np.asarray(dag.parent_offsets)
    assert off[2] == off[3]


def test_edge_mask_filters_parallel_edges():
    # A=0 -> B=1 with two parallel edges of different costs; only cheap edge should be allowed
    edges = [
        (0, 1, 1, 1, 10),
        (0, 1, 2, 1, 20),
    ]
    g = build_graph(2, edges)
    edge_mask = np.array([True, False], dtype=bool)
    dist, dag = ngc.spf(
        g,
        0,
        1,
        edge_select=ngc.EdgeSelect.ALL_MIN_COST,
        multipath=True,
        edge_mask=edge_mask,
        eps=1e-12,
    )
    assert np.isclose(dist[1], 1.0)
    # Only link_id 10 should appear in preds
    offsets = np.asarray(dag.parent_offsets)
    via = np.asarray(dag.via_edges)
    start, end = int(offsets[1]), int(offsets[2])
    assert end - start == 1
    # Map via edge id -> link id and check it's the first edge (link_id 10)
    e = int(via[start])
    assert int(g.link_id_of(e)) == 10
