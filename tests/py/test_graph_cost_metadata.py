"""StrictMultiDiGraph.max_cost / cost_gcd: the inputs of the SPF queue dispatch rule."""

from __future__ import annotations

import numpy as np

import netgraph_core as ngc


def _graph(costs):
    n = len(costs) + 1
    return ngc.StrictMultiDiGraph.from_arrays(
        n,
        src=np.arange(len(costs), dtype=np.int32),
        dst=np.arange(1, n, dtype=np.int32),
        capacity=np.ones(len(costs), dtype=np.float64),
        cost=np.array(costs, dtype=np.int64),
        ext_edge_ids=np.arange(len(costs), dtype=np.int64),
    )


def test_cost_metadata():
    assert _graph([6, 9, 0]).max_cost() == 9
    assert _graph([6, 9, 0]).cost_gcd() == 3
    assert _graph([0, 0]).max_cost() == 0
    assert _graph([0, 0]).cost_gcd() == 0
    big = (1 << 40) + 7
    assert _graph([big]).max_cost() == big
    assert _graph([big]).cost_gcd() == big
