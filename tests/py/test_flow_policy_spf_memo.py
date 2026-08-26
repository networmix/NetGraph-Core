"""FlowPolicy's SPF memo must be invisible: identical results, fewer SPF runs.

EqualBalanced placement and rebalance rounds re-request path bundles against
residual state that repeats (94% of SPF calls in a measured place/rebalance
cycle were exact input repeats), so get_path_bundle memoizes raw SPF results,
keyed by (src, dst, residual content, residual-awareness). These tests pin the
one way that can go wrong: a stale hit -- serving a cached DAG computed for
residuals that have since changed. Each scenario mutates residual state between
placements and asserts routing reflects the *current* residuals.
"""

from __future__ import annotations

import numpy as np
import pytest

import netgraph_core as ngc


def _graph(num_nodes, edges):
    """edges: (src, dst, cost, cap)."""
    s = np.array([e[0] for e in edges], np.int32)
    d = np.array([e[1] for e in edges], np.int32)
    k = np.array([int(e[2]) for e in edges], np.int64)
    c = np.array([float(e[3]) for e in edges], np.float64)
    return ngc.StrictMultiDiGraph.from_arrays(
        num_nodes, s, d, c, k, np.arange(len(edges), dtype=np.int64)
    )


def _eb_policy(algs, gh, **kw):
    sel = ngc.EdgeSelection(
        multi_edge=True,
        require_capacity=True,
        tie_break=ngc.EdgeTieBreak.DETERMINISTIC,
    )
    cfg = ngc.FlowPolicyConfig(
        path_alg=ngc.PathAlg.SPF,
        flow_placement=ngc.FlowPlacement.EQUAL_BALANCED,
        selection=sel,
        **kw,
    )
    return ngc.FlowPolicy(algs, gh, cfg)


def _two_path_graph():
    """0->1->3 (cost 2) and 0->2->3 (cost 4); cheap path capacity 2."""
    return _graph(
        4,
        [
            (0, 1, 1, 2.0),
            (1, 3, 1, 2.0),
            (0, 2, 2, 10.0),
            (2, 3, 2, 10.0),
        ],
    )


def _edges_used(fg, g, idx):
    flows = fg.get_flow_edges(ngc.FlowIndex(*idx))
    src = np.asarray(g.edge_src_view())
    dst = np.asarray(g.edge_dst_view())
    return {(int(src[e]), int(dst[e])) for e, amt in flows if amt > 1e-12}


class TestMemoInvalidation:
    def test_replacement_sees_residuals_changed_by_another_policy(self, algs):
        """A cached bundle must not survive residual changes made between calls.

        Policy A places on the cheap path and is removed (restoring residuals,
        which reproduces the memo's key). Policy B then saturates the cheap
        path through the same FlowGraph. When A places again, its SPF runs
        against the *new* residuals and must route via the expensive path; a
        stale memo hit would re-serve the cheap-path DAG.
        """
        g = _two_path_graph()
        gh = algs.build_graph(g)
        fg = ngc.FlowGraph(g)

        pa = _eb_policy(algs, gh, max_flow_count=1, min_flow_count=1)
        placed, _ = pa.place_demand(fg, 0, 3, 0, 1.0)
        assert placed == pytest.approx(1.0)
        (idx,) = pa.flows.keys()
        assert _edges_used(fg, g, idx) == {(0, 1), (1, 3)}
        pa.remove_demand(fg)  # residual content is now back to the initial bytes

        pb = _eb_policy(algs, gh, max_flow_count=1, min_flow_count=1)
        placed_b, _ = pb.place_demand(fg, 0, 3, 1, 2.0)
        assert placed_b == pytest.approx(2.0)  # saturates 0->1->3

        placed2, _ = pa.place_demand(fg, 0, 3, 0, 1.0)
        assert placed2 == pytest.approx(1.0)
        (idx2,) = pa.flows.keys()
        assert _edges_used(fg, g, idx2) == {(0, 2), (2, 3)}, (
            "placement after residual change reused a stale cached path"
        )

    def test_alternating_flow_graphs_stay_isolated(self, algs):
        """One policy asked about two FlowGraphs must answer for the right one.

        fg_busy has the cheap path saturated; fg_idle does not. Alternating
        get-path queries (via place/remove cycles) between the two must route
        differently every time, regardless of which answer the memo last held.
        """
        g = _two_path_graph()
        gh = algs.build_graph(g)
        fg_idle = ngc.FlowGraph(g)
        fg_busy = ngc.FlowGraph(g)

        blocker = _eb_policy(algs, gh, max_flow_count=1, min_flow_count=1)
        placed, _ = blocker.place_demand(fg_busy, 0, 3, 9, 2.0)
        assert placed == pytest.approx(2.0)

        prober = _eb_policy(algs, gh, max_flow_count=1, min_flow_count=1)
        for _ in range(3):
            placed, _ = prober.place_demand(fg_idle, 0, 3, 0, 1.0)
            assert placed == pytest.approx(1.0)
            (idx,) = prober.flows.keys()
            assert _edges_used(fg_idle, g, idx) == {(0, 1), (1, 3)}
            prober.remove_demand(fg_idle)

            placed, _ = prober.place_demand(fg_busy, 0, 3, 0, 1.0)
            assert placed == pytest.approx(1.0)
            (idx,) = prober.flows.keys()
            assert _edges_used(fg_busy, g, idx) == {(0, 2), (2, 3)}
            prober.remove_demand(fg_busy)

    def test_eb_rebalance_churn_matches_fresh_policy(self, algs):
        """Heavy place/rebalance churn (the memoized path) must land exactly
        where a fresh policy with no memo history lands."""
        rng = np.random.default_rng(3)
        n = 24
        edges = []
        for _ in range(70):
            u, v = rng.integers(0, n, size=2)
            if u != v:
                edges.append(
                    (
                        int(u),
                        int(v),
                        int(rng.integers(1, 21)),
                        float(rng.integers(1, 6)),
                    )
                )
        g = _graph(n, edges)
        gh = algs.build_graph(g)

        def run(cycles):
            p = _eb_policy(algs, gh, max_flow_count=4, min_flow_count=4)
            fg = ngc.FlowGraph(g)
            for _ in range(cycles):
                p.place_demand(fg, 0, n - 1, 0, 8.0)
                p.rebalance_demand(fg, 0, n - 1, 0, 2.0)
            return np.asarray(fg.edge_flow_view()).copy()

        churned = run(cycles=4)  # memo-heavy history
        fresh = run(cycles=4)  # brand-new policy, identical inputs
        np.testing.assert_array_equal(churned, fresh)
