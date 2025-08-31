import numpy as np

import netgraph_core as ngc


def test_compaction_sorts_and_groups_parallel_edges():
    # Nodes: 0->1 with mixed order edges; 0->2; expect CSR row for 0 sorted by col and grouped
    src = np.array([0, 0, 0, 1], dtype=np.int32)
    dst = np.array([2, 1, 1, 2], dtype=np.int32)
    cap = np.array([1.0, 1.0, 2.0, 1.0], dtype=np.float64)
    cost = np.array([2.0, 1.0, 1.5, 1.0], dtype=np.float64)
    g = ngc.StrictMultiDiGraph.from_arrays(
        3, src, dst, cap, cost, link_ids=None, add_reverse=False
    )
    row = np.asarray(g.row_offsets_view())
    col = np.asarray(g.col_indices_view())
    # Row for node 0 spans [row[0], row[1]) and should be non-decreasing by col
    s, e = int(row[0]), int(row[1])
    assert np.all(col[s:e][:-1] <= col[s:e][1:])
    # Parallel edges to the same neighbor should be consecutive (grouped)
    # Count transitions; for two neighbors (1 and 2), expect at most two groups
    groups = 1 if e - s <= 1 else int(np.sum(col[s:e][1:] != col[s:e][:-1]) + 1)
    assert groups <= 2
