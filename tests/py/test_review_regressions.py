"""Regression tests for the 2026-08 code review findings.

Each test locks in the fixed behavior for a confirmed defect:

1. max_flow under-reported when the SSP tier loop needed flow cancellation;
   its own min_cut then contradicted total_flow (duality violation).
2. Zero-cost edges produced a cyclic PredDAG: EqualBalanced placement returned
   0 flow and resolve_to_paths never terminated.
3. k_shortest_paths enumerated every equal-cost spur path (exponential in ECMP
   width; a 22-stage ladder needed ~19 s and ~3.8 GB for k=3).
4. batch_max_flow results must be identical to per-pair calls (now computed in
   parallel).
5. FlowPolicy.place_demand silently routed a second (src, dst) pair over the
   first pair's paths.
6. StrictMultiDiGraph accepted cost totals that overflow int64 path arithmetic.
7. FlowState.compute_min_cut treated pre-existing usage from a custom
   residual_init as cancellable flow.
8. batch_max_flow rejected a genuine int32 `pairs` array on Windows, and read
   non-contiguous input as if it were packed.
"""

from __future__ import annotations

import numpy as np
import pytest

import netgraph_core as ngc


def _graph(num_nodes, edges):
    """edges: list of (src, dst, capacity, cost)."""
    s, d, c, k = zip(*edges, strict=False)
    return ngc.StrictMultiDiGraph.from_arrays(
        num_nodes,
        np.array(s, dtype=np.int32),
        np.array(d, dtype=np.int32),
        np.array(c, dtype=np.float64),
        np.array(k, dtype=np.int64),
        np.arange(len(s), dtype=np.int64),
    )


def _bidir(links):
    out = []
    for u, v, cap, cost in links:
        out.append((u, v, cap, cost))
        out.append((v, u, cap, cost))
    return out


@pytest.fixture
def algs():
    return ngc.Algorithms(ngc.Backend.cpu())


class TestMaxFlowCompletion:
    """Finding 1: forward-only SSP stranded capacity that needs cancellation."""

    def _cancellation_graph(self):
        # s=0, a=1, b=2, t=3; all links full duplex. True max flow is 3.0:
        # 0->1->3 (1), 0->2->1->3 (1), 0->2->3 (1); cut {0,2}|{1,3} has cap 3.
        return _graph(
            4,
            _bidir(
                [
                    (1, 3, 3.0, 17),
                    (0, 2, 8.0, 17),
                    (1, 2, 1.0, 1),
                    (1, 0, 1.0, 1),
                    (3, 2, 1.0, 1),
                ]
            ),
        )

    def test_exact_max_flow(self, algs):
        g = self._cancellation_graph()
        pg = algs.build_graph(g)
        total, _ = algs.max_flow(pg, 0, 3)
        assert total == pytest.approx(3.0)

    def test_min_cut_duality(self, algs):
        g = self._cancellation_graph()
        pg = algs.build_graph(g)
        total, summary = algs.max_flow(pg, 0, 3)
        cap = np.asarray(g.capacity_view())
        cut = np.asarray(summary.min_cut.edges)
        assert cut.size > 0
        assert float(cap[cut].sum()) == pytest.approx(total)

    def test_cost_distribution_accounts_for_total(self, algs):
        g = self._cancellation_graph()
        pg = algs.build_graph(g)
        total, summary = algs.max_flow(pg, 0, 3)
        assert float(np.asarray(summary.flows).sum()) == pytest.approx(total)

    def test_completion_respects_masks(self, algs):
        g = self._cancellation_graph()
        pg = algs.build_graph(g)
        mask = np.ones(4, dtype=bool)
        mask[2] = False  # only 0->1->3 remains; bottleneck is cap(0->1) = 1
        total, _ = algs.max_flow(pg, 0, 3, node_mask=mask)
        assert total == pytest.approx(1.0)


class TestZeroCostEdges:
    """Finding 2: zero-cost edges made the PredDAG cyclic."""

    def _zero_cost_graph(self):
        # 0 ->(1) 1 <->(0) 2 ->(1) 3
        return _graph(
            4,
            [(0, 1, 10.0, 1), (1, 2, 10.0, 0), (2, 1, 10.0, 0), (2, 3, 10.0, 1)],
        )

    def test_pred_dag_acyclic(self, algs):
        g = self._zero_cost_graph()
        pg = algs.build_graph(g)
        _, dag = algs.spf(pg, 0, None, multipath=True)
        po = np.asarray(dag.parent_offsets)
        pa = np.asarray(dag.parents)
        for v in range(4):
            for i in range(po[v], po[v + 1]):
                u = pa[i]
                assert v not in pa[po[u] : po[u + 1]], f"cycle {u}<->{v}"

    def test_equal_balanced_places_flow(self, algs):
        g = self._zero_cost_graph()
        pg = algs.build_graph(g)
        total, _ = algs.max_flow(
            pg, 0, 3, flow_placement=ngc.FlowPlacement.EQUAL_BALANCED
        )
        assert total == pytest.approx(10.0)

    def test_resolve_to_paths_terminates(self, algs):
        g = self._zero_cost_graph()
        pg = algs.build_graph(g)
        _, dag = algs.spf(pg, 0, None, multipath=True)
        paths = dag.resolve_to_paths(0, 3)  # previously never returned
        assert [n for n, _ in paths[0]] == [0, 1, 2, 3]


class TestKspEcmpLadder:
    """Finding 3: spur enumeration was exponential in equal-cost path count."""

    def test_ladder_is_tractable(self, algs):
        # 30 stages of parallel 2-hop diamonds: 2^30 equal-cost paths.
        # Pre-fix this was unreachable (22 stages already took ~19 s / 3.8 GB).
        stages = 30
        edges = []
        nid = 0

        def new():
            nonlocal nid
            nid += 1
            return nid - 1

        main = [new()]
        for _ in range(stages):
            nxt, a, b = new(), new(), new()
            for u, v in ((main[-1], a), (a, nxt), (main[-1], b), (b, nxt)):
                edges.append((u, v, 10.0, 1))
            main.append(nxt)
        g = _graph(nid, edges)
        pg = algs.build_graph(g)
        out = algs.ksp(pg, main[0], main[-1], k=4)
        assert len(out) == 4
        # All returned paths are shortest (equal cost) and distinct.
        costs = [dist[main[-1]] for dist, _ in out]
        assert costs == [2 * stages] * 4
        sigs = {tuple(dag.via_edges) for _, dag in out}
        assert len(sigs) == 4


class TestBatchMaxFlowParity:
    """Finding 4: batch results must equal per-pair results (now parallel)."""

    def test_batch_equals_serial(self, algs):
        rng = np.random.default_rng(7)
        n = 12
        edges = []
        for _ in range(40):
            u, v = rng.integers(0, n, size=2)
            if u != v:
                edges.append((int(u), int(v), float(rng.integers(1, 8)), 1))
        g = _graph(n, edges)
        pg = algs.build_graph(g)
        pairs = np.array([[i, (i + 5) % n] for i in range(8)], dtype=np.int32)
        batch = algs.batch_max_flow(pg, pairs, with_edge_flows=True)
        for i, (a, b) in enumerate(pairs):
            total, summary = algs.max_flow(pg, int(a), int(b), with_edge_flows=True)
            assert batch[i].total_flow == pytest.approx(total)
            np.testing.assert_array_equal(
                np.asarray(batch[i].edge_flows), np.asarray(summary.edge_flows)
            )


class TestFlowPolicySrcDstGuard:
    """Finding 5: reusing a policy for a different pair corrupted placement."""

    def test_mismatched_pair_raises(self, algs):
        g = _graph(4, [(0, 1, 10.0, 1), (2, 3, 10.0, 1)])
        pg = algs.build_graph(g)
        fg = ngc.FlowGraph(g)
        policy = ngc.FlowPolicy(algs, pg, ngc.FlowPolicyConfig())
        placed, _ = policy.place_demand(fg, 0, 1, 0, 5.0)
        assert placed == pytest.approx(5.0)
        with pytest.raises(ValueError, match="different \\(src, dst\\)"):
            policy.place_demand(fg, 2, 3, 0, 5.0)

    def test_same_pair_and_after_remove_ok(self, algs):
        g = _graph(4, [(0, 1, 10.0, 1), (2, 3, 10.0, 1)])
        pg = algs.build_graph(g)
        fg = ngc.FlowGraph(g)
        policy = ngc.FlowPolicy(algs, pg, ngc.FlowPolicyConfig())
        assert policy.place_demand(fg, 0, 1, 0, 5.0)[0] == pytest.approx(5.0)
        assert policy.place_demand(fg, 0, 1, 0, 2.0)[0] == pytest.approx(2.0)
        policy.remove_demand(fg)
        assert policy.place_demand(fg, 2, 3, 0, 5.0)[0] == pytest.approx(5.0)
        np.testing.assert_allclose(np.asarray(fg.edge_flow_view()), [0.0, 5.0])


class TestCostOverflowGuard:
    """Finding 6: cost totals near INT64_MAX silently wrapped in SPF."""

    def test_total_at_2_pow_62_rejected(self):
        big = 1 << 62
        with pytest.raises(ValueError, match="2\\^62"):
            _graph(3, [(0, 1, 10.0, big), (1, 2, 10.0, big)])

    def test_large_but_safe_costs_accepted(self, algs):
        big = 1 << 60
        g = _graph(3, [(0, 1, 10.0, big), (1, 2, 10.0, big)])
        pg = algs.build_graph(g)
        dist, _ = algs.spf(pg, 0, None, dtype="int64")
        assert dist[2] == 2 * big  # no wraparound


class TestFlowStateCustomResidual:
    """Finding 7: custom residual_init desynced min-cut's notion of flow."""

    def test_min_cut_uses_own_placed_flow(self):
        g = _graph(3, [(0, 1, 10.0, 1), (1, 2, 10.0, 1)])
        fs = ngc.FlowState(g, np.array([4.0, 4.0]))
        placed = fs.place_max_flow(0, 2, ngc.FlowPlacement.PROPORTIONAL)
        assert placed == pytest.approx(4.0)
        # Reverse arcs must be keyed on this state's own edge_flow (4.0), not
        # capacity - residual (10.0, which counts pre-existing usage).
        np.testing.assert_allclose(np.asarray(fs.edge_flow_view()), [4.0, 4.0])
        cut = np.asarray(fs.compute_min_cut(0).edges)
        assert cut.size > 0


class TestBatchPairsDtype:
    """Finding 8: `batch_max_flow` rejected a genuine int32 array on Windows.

    The dtype check compared buffer format strings. NumPy spells int32 as
    NPY_LONG ('l') on LLP64 (Windows) but NPY_INT ('i') on LP64, while
    pybind11's `format_descriptor<int32_t>` is always 'i', so the check threw on
    Windows for a correctly-typed array. Parametrizing over both 32-bit
    spellings exercises whichever one is the platform's alias for int32.
    """

    def _line_graph(self):
        return _graph(3, [(0, 1, 5.0, 1), (1, 2, 5.0, 1)])

    @pytest.mark.parametrize("dtype", [np.int32, np.intc])
    def test_accepts_every_int32_spelling(self, algs, dtype):
        g = self._line_graph()
        pg = algs.build_graph(g)
        pairs = np.array([[0, 2]], dtype=dtype)
        assert np.dtype(dtype).itemsize == 4 and np.dtype(dtype).kind == "i"
        out = algs.batch_max_flow(pg, pairs)
        assert out[0].total_flow == pytest.approx(5.0)

    def test_rejects_wrong_dtype(self, algs):
        pg = algs.build_graph(self._line_graph())
        with pytest.raises(TypeError, match="int32"):
            algs.batch_max_flow(pg, np.array([[0, 2]], dtype=np.int64))

    def test_rejects_non_contiguous(self, algs):
        # A strided view would otherwise be read as if it were packed, silently
        # turning [[0, 2]] into the pair (0, 99).
        pg = algs.build_graph(self._line_graph())
        strided = np.array([[0, 99, 2, 99]], dtype=np.int32)[:, ::2]
        assert not strided.flags["C_CONTIGUOUS"]
        with pytest.raises(TypeError, match="C-contiguous"):
            algs.batch_max_flow(pg, strided)
