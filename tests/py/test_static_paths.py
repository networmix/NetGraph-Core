"""FlowPolicy.set_static_paths: pinned path bundles (MPLS-style routing).

Semantics under test (design-panel reviewed):
- One flow per usable bundle, each permanently bound to its own bundle in
  supply order (the original Python semantics; the first C++ port bound every
  flow to bundle[0]).
- Bundles are validated against the graph and pruned against the policy's
  masks; a bundle with no surviving src->dst walk is DOWN and creates no flow.
- Flow cost is the min-cost src->dst walk of the PRUNED bundle.
- Pinned policies never grow their flow set and never reoptimize;
  max_path_cost/-factor, min_flow_count and reoptimize_flows_on_each_placement
  are inert.
- EqualBalanced spreads over the usable (up) bundles only, and the rebalance
  loop is iterative (many heterogeneous pinned bundles must not blow the stack).
- PredDAG.from_edges builds a single-path bundle from a contiguous edge list.
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


def _eid(g, u, v):
    """Edge id of u->v after the graph's (cost, src, dst) reordering."""
    es, ed = g.edge_src_view(), g.edge_dst_view()
    (idx,) = np.where((es == u) & (ed == v))
    assert idx.size == 1, f"expected exactly one {u}->{v} edge"
    return int(idx[0])


@pytest.fixture
def algs():
    return ngc.Algorithms(ngc.Backend.cpu())


@pytest.fixture
def square(algs):
    """0->1->2 (cost 1+1, cap 1 each) and 0->3->2 (cost 2+2, cap 2 each)."""
    g = _graph(4, [(0, 1, 1.0, 1), (1, 2, 1.0, 1), (0, 3, 2.0, 2), (3, 2, 2.0, 2)])
    pg = algs.build_graph(g)
    short = ngc.PredDAG.from_edges(g, [_eid(g, 0, 1), _eid(g, 1, 2)])
    long = ngc.PredDAG.from_edges(g, [_eid(g, 0, 3), _eid(g, 3, 2)])
    return g, pg, short, long


class TestFromEdges:
    def test_builds_usable_single_path_dag(self, algs, square):
        g, pg, short, _ = square
        fs = ngc.FlowState(g)
        placed = fs.place_on_dag(
            0, 2, short, float("inf"), ngc.FlowPlacement.PROPORTIONAL
        )
        assert placed == pytest.approx(1.0)  # bottleneck of the 0->1->2 path

    def test_rejects_empty_out_of_range_non_contiguous_and_revisits(self, algs):
        g = _graph(3, [(0, 1, 1.0, 1), (1, 2, 1.0, 1), (2, 0, 1.0, 1)])
        e01, e12, e20 = _eid(g, 0, 1), _eid(g, 1, 2), _eid(g, 2, 0)
        with pytest.raises(ValueError, match="non-empty"):
            ngc.PredDAG.from_edges(g, [])
        with pytest.raises(ValueError, match="out of range"):
            ngc.PredDAG.from_edges(g, [99])
        with pytest.raises(ValueError, match="not contiguous"):
            ngc.PredDAG.from_edges(g, [e01, e20])
        with pytest.raises(ValueError, match="simple path"):
            ngc.PredDAG.from_edges(g, [e01, e12, e20])  # returns to node 0


class TestPerFlowBinding:
    """The port bug this redesign fixes: every flow bound to bundle[0]."""

    def test_each_flow_bound_to_its_own_bundle_in_supply_order(self, algs, square):
        g, pg, short, long = square
        fg = ngc.FlowGraph(g)
        policy = ngc.FlowPolicy(algs, pg, ngc.FlowPolicyConfig())
        policy.set_static_paths(0, 2, [short, long])
        placed, left = policy.place_demand(fg, 0, 2, 0, 3.0)
        assert placed == pytest.approx(3.0)  # 1 on the short path + 2 on the long
        assert left == pytest.approx(0.0)
        # Ordinal binding: creation order == supply order; costs prove which
        # bundle each flow holds (short = 2, long = 4).
        by_creation = sorted(policy.flows.items(), key=lambda kv: kv[0][3])
        assert [v[2] for _, v in by_creation] == [2, 4]
        assert [v[3] for _, v in by_creation] == pytest.approx([1.0, 2.0])

    def test_supply_order_is_placement_precedence(self, algs, square):
        """Proportional pinned placement is greedy in bundle-supply order."""
        g, pg, short, long = square
        fg = ngc.FlowGraph(g)
        policy = ngc.FlowPolicy(algs, pg, ngc.FlowPolicyConfig())
        policy.set_static_paths(0, 2, [long, short])  # long first this time
        placed, _ = policy.place_demand(fg, 0, 2, 0, 2.0)
        assert placed == pytest.approx(2.0)
        by_creation = sorted(policy.flows.items(), key=lambda kv: kv[0][3])
        # First-supplied bundle (long, cost 4) absorbed everything.
        assert [v[2] for _, v in by_creation] == [4, 2]
        assert [v[3] for _, v in by_creation] == pytest.approx([2.0, 0.0])


class TestEqualBalanced:
    def test_eb_equalizes_to_bottleneck_bundle(self, algs, square):
        g, pg, short, long = square
        fg = ngc.FlowGraph(g)
        policy = ngc.FlowPolicy(
            algs,
            pg,
            ngc.FlowPolicyConfig(flow_placement=ngc.FlowPlacement.EQUAL_BALANCED),
        )
        policy.set_static_paths(0, 2, [short, long])
        placed, _ = policy.place_demand(fg, 0, 2, 0, 3.0)
        # ECMP over 2 LSPs; the short path caps each at ~1.0.
        assert placed == pytest.approx(2.0, abs=1e-3)
        vols = [v[3] for v in policy.flows.values()]
        assert max(vols) - min(vols) < 1e-3

    def test_eb_spreads_over_up_bundles_only(self, algs, square):
        g, pg, short, long = square
        em = np.ones(4, dtype=bool)
        em[_eid(g, 0, 1)] = False  # short path down
        fg = ngc.FlowGraph(g)
        policy = ngc.FlowPolicy(
            algs,
            pg,
            ngc.FlowPolicyConfig(flow_placement=ngc.FlowPlacement.EQUAL_BALANCED),
            edge_mask=em,
        )
        policy.set_static_paths(0, 2, [short, long])
        placed, _ = policy.place_demand(fg, 0, 2, 0, 3.0)
        # One up LSP; head-end ECMP over up tunnels only -> full 2.0 on the long
        # path, not volume/2 as a supplied-count divisor would give.
        assert policy.flow_count() == 1
        assert placed == pytest.approx(2.0, abs=1e-3)

    def test_many_heterogeneous_bundles_do_not_blow_the_stack(self, algs):
        """256 pinned parallel LSPs with one small bundle: the EB rebalance loop
        (previously recursion of depth ~U * ln(imbalance/kMinFlow)) must survive."""
        n_paths = 256
        edges = []
        # src=0, dst=1, path i via node 2+i with capacity 100 except one runt.
        for i in range(n_paths):
            mid = 2 + i
            cap = 0.5 if i == 0 else 100.0
            edges.append((0, mid, cap, 1))
            edges.append((mid, 1, cap, 1))
        g = _graph(2 + n_paths, edges)
        pg = algs.build_graph(g)
        bundles = [
            ngc.PredDAG.from_edges(g, [_eid(g, 0, 2 + i), _eid(g, 2 + i, 1)])
            for i in range(n_paths)
        ]
        fg = ngc.FlowGraph(g)
        policy = ngc.FlowPolicy(
            algs,
            pg,
            ngc.FlowPolicyConfig(flow_placement=ngc.FlowPlacement.EQUAL_BALANCED),
        )
        policy.set_static_paths(0, 1, bundles)
        placed, _ = policy.place_demand(fg, 0, 1, 0, 10_000.0)
        # ECMP: the runt bundle (cap 0.5) bounds every LSP.
        assert placed == pytest.approx(n_paths * 0.5, rel=1e-2)

    def test_second_placement_is_not_incremental(self, algs, square):
        """Documented: per-call volume defines the EB per-flow target, so a second
        smaller placement adds nothing (parity with dynamic EB + max_flow_count)."""
        g, pg, short, long = square
        fg = ngc.FlowGraph(g)
        policy = ngc.FlowPolicy(
            algs,
            pg,
            ngc.FlowPolicyConfig(flow_placement=ngc.FlowPlacement.EQUAL_BALANCED),
        )
        policy.set_static_paths(0, 2, [short, long])
        first, _ = policy.place_demand(fg, 0, 2, 0, 2.0)
        assert first == pytest.approx(2.0, abs=1e-3)
        # The smaller per-call volume yields a smaller per-flow target that the
        # flows already meet, so nothing new is placed and the new volume is
        # returned as leftover.
        second, second_left = policy.place_demand(fg, 0, 2, 0, 1.0)
        assert second == pytest.approx(0.0, abs=1e-3)
        assert second_left == pytest.approx(1.0, abs=1e-3)
        assert policy.placed_demand() == pytest.approx(first, abs=1e-3)


class TestDiminishingReturnsArming:
    def test_more_bundles_than_window_all_receive_volume(self, algs):
        """16 LSPs x cap 1 against volume 10000: the diminishing-returns window
        (8) must not fire before every pinned flow was visited once."""
        n_paths = 16
        edges = []
        for i in range(n_paths):
            mid = 2 + i
            edges.append((0, mid, 1.0, 1))
            edges.append((mid, 1, 1.0, 1))
        g = _graph(2 + n_paths, edges)
        pg = algs.build_graph(g)
        bundles = [
            ngc.PredDAG.from_edges(g, [_eid(g, 0, 2 + i), _eid(g, 2 + i, 1)])
            for i in range(n_paths)
        ]
        fg = ngc.FlowGraph(g)
        policy = ngc.FlowPolicy(algs, pg, ngc.FlowPolicyConfig())
        policy.set_static_paths(0, 1, bundles)
        placed, _ = policy.place_demand(fg, 0, 1, 0, 10_000.0)
        assert placed == pytest.approx(float(n_paths))  # all 16 units of capacity


class TestMaskSemantics:
    def test_down_bundle_creates_no_flow_and_cost_is_post_prune(self, algs):
        """A bundle whose cheapest walk is masked reports the surviving walk's
        cost; a bundle with no surviving walk is down."""
        # Two parallel two-hop walks 0->1: via node 2 (cost 1+1) and node 3 (cost 5+0).
        g = _graph(
            4,
            [(0, 2, 5.0, 1), (2, 1, 5.0, 1), (0, 3, 5.0, 5), (3, 1, 5.0, 0)],
        )
        pg = algs.build_graph(g)
        cheap = [_eid(g, 0, 2), _eid(g, 2, 1)]
        dear = [_eid(g, 0, 3), _eid(g, 3, 1)]
        # One multi-walk bundle containing both walks: merge the two single-path
        # DAGs by using an SPF-style DAG is not possible here (costs differ), so
        # build the union manually via two bundles for the down test, and use the
        # dear-only bundle for the cost check after masking the cheap walk.
        b_cheap = ngc.PredDAG.from_edges(g, cheap)
        b_dear = ngc.PredDAG.from_edges(g, dear)

        em = np.ones(4, dtype=bool)
        em[cheap[0]] = False  # kill the cheap walk entirely

        fg = ngc.FlowGraph(g)
        policy = ngc.FlowPolicy(algs, pg, ngc.FlowPolicyConfig(), edge_mask=em)
        policy.set_static_paths(0, 1, [b_cheap, b_dear])
        placed, _ = policy.place_demand(fg, 0, 1, 0, 100.0)
        assert policy.flow_count() == 1  # cheap bundle is DOWN
        ((_, flow),) = policy.flows.items()
        assert flow[2] == 5  # cost of the surviving (dear) walk
        assert placed == pytest.approx(5.0)

    def test_multiwalk_bundle_survives_partial_failure_with_surviving_cost(self, algs):
        """An SPF DAG used as one bundle renormalizes over surviving walks and
        reports the post-prune min cost (SR-TE-style semantics)."""
        # Equal-cost diamond: 0->1 via 2 or via 3, all edges cost 1 cap 5.
        g = _graph(4, [(0, 2, 5.0, 1), (2, 1, 5.0, 1), (0, 3, 5.0, 1), (3, 1, 5.0, 1)])
        pg = algs.build_graph(g)
        algs_local = algs
        _, dag = algs_local.spf(pg, 0, None, multipath=True)

        em = np.ones(4, dtype=bool)
        em[_eid(g, 0, 2)] = False  # one branch of the DAG fails

        fg = ngc.FlowGraph(g)
        policy = ngc.FlowPolicy(algs, pg, ngc.FlowPolicyConfig(), edge_mask=em)
        policy.set_static_paths(0, 1, [dag])
        placed, _ = policy.place_demand(fg, 0, 1, 0, 100.0)
        assert policy.flow_count() == 1
        assert placed == pytest.approx(5.0)  # surviving branch only
        ((_, flow),) = policy.flows.items()
        assert flow[2] == 2

    def test_all_bundles_down_places_nothing_and_is_idempotent(self, algs, square):
        g, pg, short, long = square
        nm = np.ones(4, dtype=bool)
        nm[1] = False
        nm[3] = False  # both intermediate nodes down
        fg = ngc.FlowGraph(g)
        policy = ngc.FlowPolicy(algs, pg, ngc.FlowPolicyConfig(), node_mask=nm)
        policy.set_static_paths(0, 2, [short, long])
        for _ in range(2):
            placed, left = policy.place_demand(fg, 0, 2, 0, 3.0)
            assert (placed, left) == (0.0, 3.0)
            assert policy.flow_count() == 0


class TestPinnedInertness:
    def test_reoptimize_on_each_placement_is_inert(self, algs, square):
        """With reoptimization enabled (as the TE presets set), pinned flows keep
        their bundles: deliberate divergence from the original Python, where
        reoptimization silently rerouted pinned flows onto SPF paths."""
        g, pg, short, long = square
        fg = ngc.FlowGraph(g)
        policy = ngc.FlowPolicy(
            algs, pg, ngc.FlowPolicyConfig(reoptimize_flows_on_each_placement=True)
        )
        policy.set_static_paths(0, 2, [short, long])
        placed, _ = policy.place_demand(fg, 0, 2, 0, 3.0)
        assert placed == pytest.approx(3.0)
        by_creation = sorted(policy.flows.items(), key=lambda kv: kv[0][3])
        assert [v[2] for _, v in by_creation] == [2, 4]  # bindings unchanged

    def test_max_path_cost_is_inert_for_pinned_bundles(self, algs, square):
        g, pg, short, long = square
        fg = ngc.FlowGraph(g)
        policy = ngc.FlowPolicy(algs, pg, ngc.FlowPolicyConfig(max_path_cost=1))
        policy.set_static_paths(0, 2, [short, long])  # costs 2 and 4, both > 1
        placed, _ = policy.place_demand(fg, 0, 2, 0, 3.0)
        assert placed == pytest.approx(3.0)  # pinned bundles bypass the cost gate

    def test_lifecycle_remove_and_rebalance(self, algs, square):
        g, pg, short, long = square
        fg = ngc.FlowGraph(g)
        policy = ngc.FlowPolicy(algs, pg, ngc.FlowPolicyConfig())
        policy.set_static_paths(0, 2, [short, long])
        policy.place_demand(fg, 0, 2, 0, 3.0)
        policy.remove_demand(fg)
        assert np.asarray(fg.edge_flow_view()).sum() == pytest.approx(0.0)
        placed, _ = policy.place_demand(fg, 0, 2, 0, 3.0)  # re-pins from bundles
        assert placed == pytest.approx(3.0)
        placed_r, _ = policy.rebalance_demand(fg, 0, 2, 0, 1.0)
        assert placed_r == pytest.approx(2.0)  # 1.0 per flow target


class TestValidation:
    def test_repin_before_placement_replaces(self, algs, square):
        g, pg, short, long = square
        policy = ngc.FlowPolicy(algs, pg, ngc.FlowPolicyConfig())
        policy.set_static_paths(0, 2, [short])
        policy.set_static_paths(0, 2, [short, long])  # no spurious throw
        fg = ngc.FlowGraph(g)
        placed, _ = policy.place_demand(fg, 0, 2, 0, 3.0)
        assert placed == pytest.approx(3.0)
        assert policy.flow_count() == 2

    def test_repin_with_masks_pruning_is_not_spurious(self, algs, square):
        """The design-panel blocker: derived max_flow_count must not make
        re-pinning the same bundles throw once masks prune some of them."""
        g, pg, short, long = square
        em = np.ones(4, dtype=bool)
        em[_eid(g, 0, 1)] = False  # short bundle will be down (U=1 < N=2)
        policy = ngc.FlowPolicy(algs, pg, ngc.FlowPolicyConfig(), edge_mask=em)
        policy.set_static_paths(0, 2, [short, long])
        policy.set_static_paths(0, 2, [short, long])  # must not throw
        assert policy.flow_count() == 0  # no flows until placement

    def test_rejects_bundle_from_a_different_graph(self, algs, square):
        g, pg, short, _ = square
        other = _graph(
            4, [(0, 3, 1.0, 7), (3, 2, 1.0, 7), (0, 1, 1.0, 9), (1, 2, 1.0, 9)]
        )
        foreign = ngc.PredDAG.from_edges(other, [_eid(other, 0, 1), _eid(other, 1, 2)])
        policy = ngc.FlowPolicy(algs, pg, ngc.FlowPolicyConfig())
        with pytest.raises(ValueError, match="different graph|does not connect"):
            policy.set_static_paths(0, 2, [foreign])

    def test_rejects_bundle_without_src_dst_walk(self, algs, square):
        g, pg, short, _ = square
        policy = ngc.FlowPolicy(algs, pg, ngc.FlowPolicyConfig())
        with pytest.raises(ValueError, match="no src->dst walk"):
            policy.set_static_paths(1, 3, [short])  # short connects 0->2, not 1->3

    def test_set_after_placement_and_wrong_demand_throw(self, algs, square):
        g, pg, short, long = square
        fg = ngc.FlowGraph(g)
        policy = ngc.FlowPolicy(algs, pg, ngc.FlowPolicyConfig())
        policy.set_static_paths(0, 2, [short])
        policy.place_demand(fg, 0, 2, 0, 1.0)
        with pytest.raises(ValueError, match="already holds flows"):
            policy.set_static_paths(0, 2, [long])
        policy.remove_demand(fg)
        policy.set_static_paths(0, 2, [long])  # allowed again after removal
        with pytest.raises(ValueError, match="do not match demand"):
            policy.place_demand(fg, 0, 1, 0, 1.0)

    def test_config_incompatibilities(self, algs, square):
        g, pg, short, long = square
        with pytest.raises(ValueError, match="max_flow_count"):
            ngc.FlowPolicy(
                algs, pg, ngc.FlowPolicyConfig(max_flow_count=3)
            ).set_static_paths(0, 2, [short, long])
        with pytest.raises(ValueError, match="shortest_path"):
            ngc.FlowPolicy(
                algs, pg, ngc.FlowPolicyConfig(shortest_path=True)
            ).set_static_paths(0, 2, [short])
        with pytest.raises(ValueError, match="non-empty"):
            ngc.FlowPolicy(algs, pg, ngc.FlowPolicyConfig()).set_static_paths(0, 2, [])
        with pytest.raises(ValueError, match="out of range or equal"):
            ngc.FlowPolicy(algs, pg, ngc.FlowPolicyConfig()).set_static_paths(
                2, 2, [short]
            )
