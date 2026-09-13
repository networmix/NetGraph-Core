"""EQUAL_BALANCED and EQUAL_BALANCED_LOSSY placement semantics.

Both model a hop-by-hop forwarding table that does not react to load: the
split set is every shortest-path edge with capacity. EQUAL_BALANCED admits
losslessly (a full member blocks admission on that DAG), EQUAL_BALANCED_LOSSY
forwards best-effort (each member carries what it can and the rest is dropped).
"""

from __future__ import annotations

import numpy as np
import pytest

import netgraph_core as ngc


def _graph(src, dst, cap, cost, n):
    return ngc.StrictMultiDiGraph.from_arrays(
        num_nodes=n,
        src=np.array(src, dtype=np.int32),
        dst=np.array(dst, dtype=np.int32),
        capacity=np.array(cap, dtype=np.float64),
        cost=np.array(cost, dtype=np.int64),
        ext_edge_ids=np.arange(len(src), dtype=np.int64),
    )


SEL = ngc.EdgeSelection(
    multi_edge=True, require_capacity=False, tie_break=ngc.EdgeTieBreak.DETERMINISTIC
)


@pytest.fixture
def pair():
    """0 -> 1 over two equal-cost parallel edges: edge 0 cap 10, edge 1 cap 100."""
    g = _graph([0, 0], [1, 1], [10.0, 100.0], [1, 1], 2)
    algs = ngc.Algorithms(ngc.Backend.cpu())
    handle = algs.build_graph(g)
    _, dag = algs.spf(
        handle, src=0, dst=None, selection=SEL, multipath=True, dtype="float64"
    )
    return g, algs, handle, dag


def test_enum_members():
    assert set(ngc.FlowPlacement.__members__) == {
        "PROPORTIONAL",
        "EQUAL_BALANCED",
        "EQUAL_BALANCED_LOSSY",
    }


def test_full_member_blocks_admission_on_the_same_dag(pair):
    g, _, _, dag = pair
    fg = ngc.FlowGraph(g)
    assert fg.place(
        ngc.FlowIndex(0, 1, 0, 0), 0, 1, dag, 20.0, ngc.FlowPlacement.EQUAL_BALANCED
    ) == pytest.approx(20.0)
    assert fg.place(
        ngc.FlowIndex(0, 1, 0, 1), 0, 1, dag, 10.0, ngc.FlowPlacement.EQUAL_BALANCED
    ) == pytest.approx(0.0)
    assert fg.edge_flow_view()[1] == pytest.approx(10.0)


def test_progress_comes_from_a_residual_aware_dag(pair):
    g, algs, handle, dag = pair
    fg = ngc.FlowGraph(g)
    assert fg.place(
        ngc.FlowIndex(0, 1, 0, 0), 0, 1, dag, 20.0, ngc.FlowPlacement.EQUAL_BALANCED
    ) == pytest.approx(20.0)
    residual = np.ascontiguousarray(fg.residual_view(), dtype=np.float64)
    _, fresh = algs.spf(
        handle,
        src=0,
        dst=None,
        selection=SEL,
        residual=residual,
        multipath=True,
        dtype="float64",
    )
    assert fg.place(
        ngc.FlowIndex(0, 1, 0, 1), 0, 1, fresh, 10.0, ngc.FlowPlacement.EQUAL_BALANCED
    ) == pytest.approx(10.0)
    assert fg.edge_flow_view()[1] == pytest.approx(20.0)


def test_cost_only_max_flow_is_single_pass(pair):
    """place_max_flow with require_capacity=False never re-places on the stale DAG."""
    g, _, _, _ = pair
    fs = ngc.FlowState(g)
    total = fs.place_max_flow(
        0,
        1,
        flow_placement=ngc.FlowPlacement.EQUAL_BALANCED,
        shortest_path=False,
        require_capacity=False,
    )
    assert total == pytest.approx(20.0)


def test_lossy_delivers_and_reports_drops(pair):
    g, _, _, dag = pair
    fg = ngc.FlowGraph(g)
    placed, drops = fg.place_with_drops(
        ngc.FlowIndex(0, 1, 0, 0),
        0,
        1,
        dag,
        100.0,
        ngc.FlowPlacement.EQUAL_BALANCED_LOSSY,
    )
    assert placed == pytest.approx(60.0)
    assert drops == [(0, pytest.approx(40.0))]
    assert fg.edge_flow_view().tolist() == pytest.approx([10.0, 50.0])
    placed2, drops2 = fg.place_with_drops(
        ngc.FlowIndex(0, 1, 0, 1),
        0,
        1,
        dag,
        10.0,
        ngc.FlowPlacement.EQUAL_BALANCED_LOSSY,
    )
    assert placed2 == pytest.approx(5.0)
    assert drops2 == [(0, pytest.approx(5.0))]


def test_place_with_drops_is_empty_for_other_modes(pair):
    g, _, _, dag = pair
    fg = ngc.FlowGraph(g)
    placed, drops = fg.place_with_drops(
        ngc.FlowIndex(0, 1, 0, 0), 0, 1, dag, 100.0, ngc.FlowPlacement.EQUAL_BALANCED
    )
    assert placed == pytest.approx(20.0)
    assert drops == []


def test_lossy_ledger_holds_carried_volume_only(pair):
    g, _, _, dag = pair
    fg = ngc.FlowGraph(g)
    idx = ngc.FlowIndex(0, 1, 0, 3)
    placed, _ = fg.place_with_drops(
        idx, 0, 1, dag, 100.0, ngc.FlowPlacement.EQUAL_BALANCED_LOSSY
    )
    assert sum(a for _, a in fg.get_flow_edges(idx)) == pytest.approx(placed)
    fg.remove(idx)
    assert fg.residual_view().tolist() == pytest.approx([10.0, 100.0])


def test_cost_only_flow_policy_does_not_reroute_around_saturation():
    """require_capacity=False must route on cost alone, even with an EB per-flow target.

    A -> B direct (cap 10, cost 1); A -> C -> B (cap 100, cost 5 + 5).
    """
    g = _graph([0, 0, 2], [1, 2, 1], [10.0, 100.0, 100.0], [1, 5, 5], 3)
    algs = ngc.Algorithms(ngc.Backend.cpu())
    handle = algs.build_graph(g)

    cfg = ngc.FlowPolicyConfig()
    cfg.path_alg = ngc.PathAlg.SPF
    cfg.flow_placement = ngc.FlowPlacement.EQUAL_BALANCED
    cfg.selection = SEL
    cfg.require_capacity = False
    cfg.shortest_path = True
    cfg.min_flow_count = 1
    cfg.max_flow_count = 1

    fg = ngc.FlowGraph(g)
    _, dag = algs.spf(
        handle, src=0, dst=None, selection=SEL, multipath=True, dtype="float64"
    )
    assert fg.place(
        ngc.FlowIndex(0, 1, 9, 0), 0, 1, dag, 10.0, ngc.FlowPlacement.EQUAL_BALANCED
    ) == pytest.approx(10.0)

    policy = ngc.FlowPolicy(algs, handle, cfg)
    placed, remaining = policy.place_demand(fg, 0, 1, 0, 50.0)
    assert placed == pytest.approx(0.0), (
        "cost-only routing must not discover the A->C->B detour"
    )
    assert remaining == pytest.approx(50.0)
    assert all(float(v[2]) == 1.0 for v in policy.flows.values())


def test_lossy_static_paths_carry_what_fits_without_equalizing():
    """Pinned routes under EQUAL_BALANCED_LOSSY: each LSP is offered its share
    and delivers what fits. A -> B direct cap 10 (edge 0); A -> C -> B cap 100
    (edges 1, 2). Demand 50 over both routes: 25 offered each, 10 + 25 = 35.
    """
    g = _graph([0, 0, 2], [1, 2, 1], [10.0, 100.0, 100.0], [1, 1, 1], 3)
    algs = ngc.Algorithms(ngc.Backend.cpu())
    handle = algs.build_graph(g)
    bundles = [ngc.PredDAG.from_edges(g, [0]), ngc.PredDAG.from_edges(g, [1, 2])]

    def policy(placement):
        cfg = ngc.FlowPolicyConfig()
        cfg.path_alg = ngc.PathAlg.SPF
        cfg.flow_placement = placement
        cfg.selection = SEL
        cfg.require_capacity = False
        cfg.min_flow_count = 1
        cfg.max_flow_count = 2
        p = ngc.FlowPolicy(algs, handle, cfg)
        p.set_static_paths(0, 1, bundles)
        return p

    fg = ngc.FlowGraph(g)
    placed, remaining = policy(ngc.FlowPlacement.EQUAL_BALANCED_LOSSY).place_demand(
        fg, 0, 1, 0, 50.0
    )
    assert placed == pytest.approx(35.0)
    assert remaining == pytest.approx(15.0)
    assert fg.edge_flow_view().tolist() == pytest.approx([10.0, 25.0, 25.0])

    fg = ngc.FlowGraph(g)
    placed, _ = policy(ngc.FlowPlacement.EQUAL_BALANCED).place_demand(fg, 0, 1, 0, 50.0)
    assert placed == pytest.approx(20.0, abs=1e-3), (
        "lossless: equal carried share, bottleneck 10"
    )
