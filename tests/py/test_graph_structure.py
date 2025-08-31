"""Graph structure tests.

Validate StrictMultiDiGraph compaction (deterministic ordering and grouping)
and CSR/reverse CSR integrity at a basic level.
"""

import numpy as np


def test_compaction_sorts_and_groups_parallel_edges(build_graph):
    # Nodes: 0->1 with mixed order edges; 0->2; expect CSR row for 0 sorted by col and grouped
    edges = [
        (0, 2, 2.0, 1.0),
        (0, 1, 1.0, 1.0),
        (0, 1, 1.5, 2.0),
        (1, 2, 1.0, 1.0),
    ]
    g = build_graph(3, edges)
    row = np.asarray(g.row_offsets_view())
    col = np.asarray(g.col_indices_view())
    # Row for node 0 spans [row[0], row[1]) and should be non-decreasing by col
    s, e = int(row[0]), int(row[1])
    assert np.all(col[s:e][:-1] <= col[s:e][1:])
    # Parallel edges to the same neighbor should be consecutive (grouped)
    # Count transitions; for two neighbors (1 and 2), expect at most two groups
    groups = 1 if e - s <= 1 else int(np.sum(col[s:e][1:] != col[s:e][:-1]) + 1)
    assert groups <= 2
