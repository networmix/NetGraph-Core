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
9-10. FlowPolicy accepted a FlowGraph built from a different graph, and
   rebalance_demand bypassed the (src, dst) guard.
11-12. batch_max_flow and ksp validated inputs later (or not at all) compared
   with their single-pair counterparts.
13. _docs.py declared types that did not match runtime.
14. The same class of user error raised different exception types.
15. k_shortest_paths read paths.back() on an empty vector when
   max_cost_factor < 1.0 and k > 1.
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

    def test_batch_equals_serial(self, algs, certify_max_flow):
        rng = np.random.default_rng(7)
        n = 12
        edges = []
        for _ in range(40):
            u, v = rng.integers(0, n, size=2)
            if u != v:
                # Costs must be dispersed. With a single cost tier there is nothing
                # for shortest-path tier ordering to get wrong, so a uniform-cost
                # generator only ever samples the region where max-flow defects of
                # the Finding 1 kind cannot occur.
                cap = float(rng.integers(1, 8))
                cost = int(rng.integers(1, 21))
                edges.append((int(u), int(v), cap, cost))
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
            # Agreement between two code paths that share calc_max_flow only proves
            # consistency; both could be wrong together. Certify the result too.
            certify_max_flow(g, summary, total, int(a), int(b))


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


class TestFlowPolicyTargetGuards:
    """Findings 9-10: FlowPolicy accepted uses that silently produced wrong answers."""

    def _two_graphs_same_edge_count(self):
        # Same edge count, different topology: the residual-length check that caught
        # the mismatched-graph case by accident does not fire here.
        a = _graph(4, [(0, 1, 5.0, 1), (1, 3, 5.0, 1), (0, 2, 1.0, 1)])
        b = _graph(4, [(0, 2, 9.0, 1), (2, 3, 9.0, 1), (0, 1, 1.0, 1)])
        return a, b

    def test_flowgraph_from_a_different_graph_is_rejected(self, algs):
        ga, gb = self._two_graphs_same_edge_count()
        policy = ngc.FlowPolicy(algs, algs.build_graph(ga), ngc.FlowPolicyConfig())
        with pytest.raises(ValueError, match="different StrictMultiDiGraph"):
            policy.place_demand(ngc.FlowGraph(gb), 0, 3, 0, 100.0)

    def test_rebalance_demand_rejects_a_different_pair(self, algs):
        # remove_demand() empties flows_, so rebalance_demand must check the pair
        # itself rather than relying on the check inside place_demand.
        g = _graph(4, [(0, 1, 10.0, 1), (2, 3, 10.0, 1)])
        fg = ngc.FlowGraph(g)
        policy = ngc.FlowPolicy(algs, algs.build_graph(g), ngc.FlowPolicyConfig())
        policy.place_demand(fg, 0, 1, 0, 5.0)
        with pytest.raises(ValueError, match="different \\(src, dst\\)"):
            policy.rebalance_demand(fg, 2, 3, 0, 5.0)
        # The matching pair still rebalances, and retargeting after remove_demand works.
        policy.rebalance_demand(fg, 0, 1, 0, 2.0)
        policy.remove_demand(fg)
        assert policy.place_demand(fg, 2, 3, 0, 5.0)[0] == pytest.approx(5.0)


class TestBindingValidationConsistency:
    """Findings 11-12: validation that differed between analogous entry points."""

    def _line(self):
        return _graph(3, [(0, 1, 5.0, 1), (1, 2, 5.0, 1)])

    def test_batch_max_flow_rejects_out_of_range_like_max_flow(self, algs):
        pg = algs.build_graph(self._line())
        # Previously this returned total_flow=0.0 and the bad id vanished into the batch.
        with pytest.raises(ValueError, match="out of range"):
            algs.batch_max_flow(pg, np.array([[0, 99]], dtype=np.int32))
        with pytest.raises(ValueError, match="out of range"):
            algs.max_flow(pg, 0, 99)

    def test_ksp_validates_dtype_even_when_no_paths_exist(self, algs):
        pg = algs.build_graph(self._line())
        # 2 -> 0 is unreachable, so the results loop never runs; the dtype check used
        # to live inside that loop and silently accepted the bad value.
        with pytest.raises(ValueError, match="dtype must be"):
            algs.ksp(pg, 2, 0, k=1, dtype="float32")

    def test_algorithms_backend_argument_is_nameable(self):
        assert isinstance(ngc.Algorithms(backend=ngc.Backend.cpu()), ngc.Algorithms)

    def test_unreachable_flow_record_class_is_not_exported(self):
        import _netgraph_core

        # Nothing could construct or receive it, and the name collided with the C++
        # alias `using Flow = double`.
        assert not hasattr(_netgraph_core, "Flow")


class TestTypeStubAccuracy:
    """Finding 13: _docs.py declared types that did not match runtime."""

    def test_min_cut_edges_is_an_int32_array_not_a_list(self, algs):
        g = _graph(3, [(0, 1, 5.0, 1), (1, 2, 5.0, 1)])
        _, summary = algs.max_flow(algs.build_graph(g), 0, 2)
        edges = summary.min_cut.edges
        assert isinstance(edges, np.ndarray)
        assert edges.dtype == np.int32

    def test_stub_declares_only_names_that_exist(self):
        import _netgraph_core

        from netgraph_core import _docs

        declared = {
            n
            for n in dir(_docs)
            if not n.startswith("_") and isinstance(getattr(_docs, n), type)
        }
        missing = {n for n in declared if not hasattr(_netgraph_core, n)}
        assert not missing, f"_docs.py declares nonexistent runtime types: {missing}"


class TestErrorTypeConsistency:
    """Finding 14: the same class of user error raised different exception types."""

    @pytest.fixture
    def pg(self, algs):
        return algs.build_graph(_graph(3, [(0, 1, 5.0, 1), (1, 2, 5.0, 1)]))

    def test_spf_residual_length_raises_type_error_like_masks(self, algs, pg):
        # Previously the residual length check fell through to the C++ core, which
        # throws std::invalid_argument -> ValueError, while every mask length check
        # in the binding layer raises TypeError.
        with pytest.raises(TypeError, match="residual length must equal"):
            algs.spf(pg, 0, residual=np.zeros(5))
        with pytest.raises(TypeError, match="node_mask length must equal"):
            algs.spf(pg, 0, node_mask=np.ones(7, dtype=bool))

    def test_wrong_typed_arguments_raise_type_error_not_runtime_error(self, algs, pg):
        # These hand-rolled casts used to surface as
        # "RuntimeError: Unable to cast ... to C++ type '?'", which no `except
        # TypeError` catches and which names nothing actionable.
        with pytest.raises(TypeError, match="must be a StrictMultiDiGraph"):
            algs.build_graph(5)
        with pytest.raises(TypeError, match="must be an Algorithms"):
            ngc.FlowPolicy("not-algorithms", pg, ngc.FlowPolicyConfig())


class TestKspCostCeiling:
    """Finding 15: k>1 with max_cost_factor < 1.0 read paths.back() on an empty vector.

    The ceiling lands below the shortest path, so nothing is admitted. The k == 1
    branch returned early and masked the undefined behaviour; k >= 2 fell into the
    spur loop, which opens with `paths.back()`.
    """

    @pytest.mark.parametrize("k", [1, 2, 5])
    @pytest.mark.parametrize("factor", [0.5, 0.99])
    def test_sub_unit_factor_returns_empty_for_every_k(self, algs, k, factor):
        g = _graph(3, [(0, 1, 5.0, 1), (1, 2, 5.0, 1)])
        pg = algs.build_graph(g)
        assert algs.ksp(pg, 0, 2, k=k, max_cost_factor=factor) == []

    @pytest.mark.parametrize("k", [1, 2, 5])
    def test_unit_factor_still_admits_the_shortest_path(self, algs, k):
        g = _graph(3, [(0, 1, 5.0, 1), (1, 2, 5.0, 1)])
        pg = algs.build_graph(g)
        assert len(algs.ksp(pg, 0, 2, k=k, max_cost_factor=1.0)) == 1


class TestStubConstructorConsistency:
    """Finding 16: widened stubs declared no __init__ for classes that need one.

    `_docs.py` is wired to type checkers, so a stub whose constructor disagrees
    with the binding turns "unchecked" into "confidently wrong": pyright rejected
    every real `FlowIndex(src, dst, cls, id)` call in the downstream consumer.
    This repo's own pyright missed it because `tests/**` is excluded, so assert
    the correspondence at runtime instead.
    """

    def _runtime_classes(self):
        import _netgraph_core

        for name in sorted(n for n in dir(_netgraph_core) if not n.startswith("_")):
            obj = getattr(_netgraph_core, name)
            if isinstance(obj, type):
                yield name, obj

    def test_stub_declares_init_wherever_construction_needs_arguments(self):
        from netgraph_core import _docs

        missing = []
        for name, rt in self._runtime_classes():
            stub = getattr(_docs, name, None)
            if stub is None:
                continue  # not part of the typed surface
            try:
                rt()
                continue  # zero-arg construction works; a bare stub is fine
            except Exception:
                pass
            doc = rt.__init__.__doc__ or ""
            if "__init__" not in doc:
                continue  # no bound constructor at all (e.g. result-only types)
            if "__init__" not in vars(stub):
                missing.append(name)
        assert not missing, (
            f"stub classes need an __init__ matching the binding: {missing}"
        )

    def test_stub_attributes_exist_at_runtime(self):
        from netgraph_core import _docs

        wrong = []
        for name, rt in self._runtime_classes():
            stub = getattr(_docs, name, None)
            if stub is None:
                continue
            declared = {
                a
                for a in vars(stub)
                if not a.startswith("_") and not callable(vars(stub)[a])
            }
            declared |= set(getattr(stub, "__annotations__", {}))
            for attr in declared:
                if not hasattr(rt, attr):
                    wrong.append(f"{name}.{attr}")
        assert not wrong, f"_docs.py declares attributes absent at runtime: {wrong}"

    def test_flow_index_constructor_matches_the_binding(self):
        # The exact call shape the downstream consumer uses.
        idx = ngc.FlowIndex(1, 2, 3, 4)
        assert (idx.src, idx.dst, idx.flowClass, idx.flowId) == (1, 2, 3, 4)
        with pytest.raises(TypeError):
            ngc.FlowIndex(src=1, dst=2, flowClass=3, flowId=4)  # positional-only
        with pytest.raises(AttributeError):
            idx.src = 9  # read-only, as the stub's properties declare


class TestFlowPolicyMaxPathCostFactor:
    """`max_path_cost_factor` had no test coverage at all.

    The gate multiplies the best path cost by the factor and casts back to Cost.
    Before 0.8.0 it did that unconditionally, including while `best_path_cost_` was
    still the INT64_MAX "no best path yet" sentinel -- the product is then far outside
    the int64 range and the conversion is undefined behavior. The fix skips the
    relative bound until a best cost exists and drops it when the product would not
    fit. Neither branch was exercised by any test.

    Scope: these are *coverage* tests, not detectors of the original UB. They were run
    against 0.7.2 (pre-fix) and pass there too, because an out-of-range float->int
    conversion is undefined rather than reliably wrong -- in practice it saturates to
    something that still admits the path. What discriminates fixed from unfixed here is
    UBSan (`make sanitize-test`), and these tests are what gives it something to
    instrument on this path. Their standalone value is guarding the *functional*
    behaviour of the bound during future refactors.
    """

    @staticmethod
    def _policy(algs, gh, *, max_flow_count=1, **cfg_kwargs):
        sel = ngc.EdgeSelection(
            multi_edge=True,
            require_capacity=True,
            tie_break=ngc.EdgeTieBreak.DETERMINISTIC,
        )
        cfg = ngc.FlowPolicyConfig(
            path_alg=ngc.PathAlg.SPF,
            flow_placement=ngc.FlowPlacement.PROPORTIONAL,
            selection=sel,
            max_flow_count=max_flow_count,
            **cfg_kwargs,
        )
        return ngc.FlowPolicy(algs, gh, cfg)

    def test_factor_applies_before_any_best_cost_exists(self, algs):
        """First placement: the sentinel path must not be multiplied and cast."""
        g = _graph(3, [(0, 1, 1.0, 5), (1, 2, 1.0, 5)])
        gh = algs.build_graph(g)
        policy = self._policy(algs, gh, max_path_cost_factor=2.0)
        placed, left = policy.place_demand(ngc.FlowGraph(g), 0, 2, 0, 1.0)
        assert placed == pytest.approx(1.0)
        assert left == pytest.approx(0.0)

    def test_huge_factor_does_not_wrap_the_bound(self, algs):
        """best_cost * factor overflows int64; the bound must be dropped, not wrapped.

        A wrapped (negative) bound would reject the shortest path itself, so a
        successful placement is what distinguishes the guard from the bug.
        """
        big = 1 << 60  # legal: total stays below the 2**62 construction ceiling
        g = _graph(3, [(0, 1, 1.0, big), (1, 2, 1.0, big)])
        gh = algs.build_graph(g)
        policy = self._policy(algs, gh, max_path_cost_factor=1e9)
        placed, _ = policy.place_demand(ngc.FlowGraph(g), 0, 2, 0, 1.0)
        assert placed == pytest.approx(1.0)

    def test_factor_still_excludes_a_too_expensive_alternative(self, algs):
        """The bound must remain functional, not merely safe."""
        g = _graph(
            4,
            [
                (0, 1, 1.0, 1),
                (1, 3, 1.0, 1),  # cost 2 route
                (0, 2, 5.0, 50),
                (2, 3, 5.0, 50),  # cost 100 route, far beyond 2x
            ],
        )
        gh = algs.build_graph(g)
        policy = self._policy(
            algs, gh, max_flow_count=4, max_path_cost_factor=2.0, min_flow_count=1
        )
        placed, _ = policy.place_demand(ngc.FlowGraph(g), 0, 3, 0, 6.0)
        # Only the cost-2 route is within 2x; the cost-100 route must stay unused.
        assert placed == pytest.approx(1.0)
