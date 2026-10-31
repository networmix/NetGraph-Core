"""Terminal zero-cost ties retain all ECMP destinations, with exact ledgers."""

import numpy as np
import pytest

import netgraph_core as c

EB = c.FlowPlacement.EQUAL_BALANCED
LOSSY = c.FlowPlacement.EQUAL_BALANCED_LOSSY
SEL = c.EdgeSelection(multi_edge=True, require_capacity=False)


def _graph(edges, n):
    src, dst, caps, costs = zip(*edges, strict=True)
    g = c.StrictMultiDiGraph.from_arrays(
        n,
        np.array(src, dtype=np.int32),
        np.array(dst, dtype=np.int32),
        np.array(caps, dtype=np.float64),
        np.array(costs, dtype=np.int64),
        np.arange(len(edges), dtype=np.int64),
    )
    a = c.Algorithms(c.Backend.cpu())
    return g, a, a.build_graph(g)


def _input(g, values):
    return np.asarray(values)[np.argsort(g.ext_edge_ids_view())]


def _state(g, fg, flows, drops, expected_drops, source, target, delivered):
    assert _input(g, fg.edge_flow_view()) == pytest.approx(flows, rel=0, abs=1e-9)
    assert _input(g, fg.residual_view()) == pytest.approx(
        _input(g, g.capacity_view()) - flows, rel=0, abs=1e-9
    )
    dv = np.zeros(g.num_edges())
    balance = np.zeros(g.num_nodes())
    for e, value in drops:
        dv[e] += value
        balance[g.edge_src_view()[e]] -= value
    assert _input(g, dv) == pytest.approx(expected_drops, rel=0, abs=1e-9)
    for u, v, value in zip(
        g.edge_src_view(), g.edge_dst_view(), fg.edge_flow_view(), strict=True
    ):
        balance[u] -= value
        balance[v] += value
    balance[source] += delivered + sum(expected_drops)
    balance[target] -= delivered
    assert balance == pytest.approx(np.zeros(g.num_nodes()), rel=0, abs=1e-9)


@pytest.mark.parametrize("caps", [(10, 5), (10, 8), (5, 10), (10, 0)])
@pytest.mark.parametrize("placement", [EB, LOSSY])
def test_anycast_terminal_tie_placement_and_policy(caps, placement):
    c1, c2 = caps
    g, a, h = _graph(
        [(0, 1, c1, 1), (0, 2, c2, 1), (1, 3, 1000, 0), (2, 3, 1000, 0)], 4
    )
    d, dag = a.spf(h, 0, 3, selection=SEL)
    assert d.tolist() == [0, 1, 1, 1]
    assert set(map(int, g.ext_edge_ids_view()[dag.via_edges])) == {0, 1, 2, 3}
    paths = dag.resolve_to_paths(0, 3)
    assert {tuple(int(u) for u, _ in p) for p in paths} == {(0, 1, 3), (0, 2, 3)}
    sent = [c1, c2] if placement == LOSSY else [min(caps)] * 2
    expected = sum(sent)
    expected_drops = [20 - c1, 20 - c2, 0, 0] if placement == LOSSY else [0] * 4
    fg = c.FlowGraph(g)
    placed, drops = fg.place_with_drops(
        c.FlowIndex(0, 3, 0, 0), 0, 3, dag, 40, placement
    )
    assert placed == expected
    _state(g, fg, sent * 2, drops, expected_drops, 0, 3, expected)

    # The dynamic policy must consume the same complete terminal route.
    policy = c.FlowPolicy(
        a,
        h,
        c.FlowPolicyConfig(
            selection=SEL,
            flow_placement=placement,
            shortest_path=True,
            require_capacity=False,
            min_flow_count=1,
            max_flow_count=1,
            reoptimize_flows_on_each_placement=False,
        ),
    )
    pf = c.FlowGraph(g)
    assert policy.place_demand(pf, 0, 3, flowClass=0, volume=40) == (
        expected,
        40 - expected,
    )
    assert _input(g, pf.edge_flow_view()) == pytest.approx(sent * 2, rel=0, abs=1e-9)
    assert _input(g, pf.residual_view()) == pytest.approx(
        [c1 - sent[0], c2 - sent[1], 1000 - sent[0], 1000 - sent[1]],
        rel=0,
        abs=1e-9,
    )


@pytest.mark.parametrize("placement", [EB, LOSSY])
@pytest.mark.parametrize("shape", ["shared-prefix", "unequal-hops", "transit"])
def test_terminal_tie_after_multihop_or_transit(shape, placement):
    lossy = placement == LOSSY
    if shape == "shared-prefix":
        edges = [
            (0, 1, 100, 1),
            (1, 2, 10, 1),
            (1, 3, 8, 1),
            (2, 4, 1000, 0),
            (3, 4, 1000, 0),
        ]
        flows = [40, 10, 8, 10, 8] if lossy else [16, 8, 8, 8, 8]
        drops = [0, 10, 12, 0, 0] if lossy else [0] * 5
        target, cost = 4, 2
    elif shape == "unequal-hops":
        edges = [
            (0, 1, 10, 2),
            (0, 2, 100, 1),
            (2, 3, 8, 1),
            (1, 4, 1000, 0),
            (3, 4, 1000, 0),
        ]
        flows = [10, 20, 8, 10, 8] if lossy else [8] * 5
        drops = [10, 0, 12, 0, 0] if lossy else [0] * 5
        target, cost = 4, 2
    else:
        edges = [
            (0, 1, 10, 1),
            (0, 2, 8, 1),
            (1, 3, 1000, 0),
            (2, 3, 1000, 0),
            (1, 2, 100, 1),
        ]
        flows = [10, 8, 10, 8, 0] if lossy else [8, 8, 8, 8, 0]
        drops = [10, 12, 0, 0, 0] if lossy else [0] * 5
        target, cost = 3, 1
    g, a, h = _graph(edges, target + 1)
    for dst in (None, target):
        d, dag = a.spf(h, 0, dst, selection=SEL)
        assert d[target] == cost
        assert len(dag.resolve_to_paths(0, target)) == 2
        fg = c.FlowGraph(g)
        placed, actual_drops = fg.place_with_drops(
            c.FlowIndex(0, target, 0, 0), 0, target, dag, 40, placement
        )
        assert placed == (18 if lossy else 16)
        _state(g, fg, flows, actual_drops, drops, 0, target, placed)


@pytest.mark.parametrize(
    "caps,expected", [((10, 5), 10), ((10, 8), 16), ((5, 10), 10), ((10, 0), 0)]
)
def test_terminal_ecmp_max_flow_and_batch(caps, expected, monkeypatch):
    g, a, h = _graph(
        [(0, 1, caps[0], 1), (0, 2, caps[1], 1), (1, 3, 1000, 0), (2, 3, 1000, 0)], 4
    )
    kw = dict(
        flow_placement=EB,
        shortest_path=True,
        require_capacity=False,
        with_edge_flows=True,
        with_residuals=True,
    )
    total, s = a.max_flow(h, 0, 3, **kw)
    assert total == s.total_flow == expected
    assert _input(g, s.edge_flows).tolist() == [min(caps)] * 4
    assert _input(g, s.residual_capacity).tolist() == [
        caps[0] - min(caps),
        caps[1] - min(caps),
        1000 - min(caps),
        1000 - min(caps),
    ]
    assert list(s.costs) == ([1] if expected else [])
    assert list(s.flows) == ([expected] if expected else [])
    monkeypatch.setenv("NGRAPH_CORE_BATCH_THREADS", "2")
    results = a.batch_max_flow(h, np.array([[0, 3]] * 4, dtype=np.int32), **kw)
    assert len(results) == 4
    for r in results:
        assert r.total_flow == expected
        assert np.array_equal(r.edge_flows, s.edge_flows)
        assert np.array_equal(r.residual_capacity, s.residual_capacity)
        assert list(r.costs) == list(s.costs)
        assert list(r.flows) == list(s.flows)


def test_terminal_sensitivity_and_true_max_flow_certificate():
    g, a, h = _graph(
        [(0, 1, 10, 2), (0, 2, 100, 1), (2, 3, 8, 1), (1, 4, 1000, 0), (3, 4, 1000, 0)],
        5,
    )
    delta = a.sensitivity_analysis(
        h, 0, 4, flow_placement=EB, shortest_path=True, require_capacity=False
    )
    # Only x->d2 saturates: deleting it leaves the 10-unit direct route.
    assert [(int(g.ext_edge_ids_view()[e]), v) for e, v in delta] == [(2, 6)]
    for masked, expected in [(False, 18), (True, 10)]:
        em = np.array([not masked or e != 2 for e in g.ext_edge_ids_view()], dtype=bool)
        total, s = a.max_flow(
            h,
            0,
            4,
            flow_placement=c.FlowPlacement.PROPORTIONAL,
            require_capacity=True,
            with_edge_flows=True,
            with_residuals=True,
            edge_mask=em,
        )
        assert total == s.total_flow == expected
        expected_flows = [10, 0, 0, 10, 0] if masked else [10, 8, 8, 10, 8]
        assert _input(g, s.edge_flows).tolist() == expected_flows
        assert np.all(s.edge_flows >= 0)
        assert np.all(s.edge_flows <= g.capacity_view())
        balance = np.zeros(5)
        for u, v, f in zip(
            g.edge_src_view(), g.edge_dst_view(), s.edge_flows, strict=True
        ):
            balance[u] -= f
            balance[v] += f
        assert balance.tolist() == [-expected, 0, 0, 0, expected]
        cut = list(s.min_cut.edges)
        assert sum(g.capacity_view()[cut]) == expected
        assert _input(g, s.residual_capacity).tolist() == [
            e[2] - f
            for e, f in zip(
                [(0, 1, 10), (0, 2, 100), (2, 3, 8), (1, 4, 1000), (3, 4, 1000)],
                expected_flows,
                strict=True,
            )
        ]


@pytest.mark.parametrize("masked,offer", [(False, 0), (True, 40)])
def test_terminal_empty_placement(masked, offer):
    g, a, h = _graph([(0, 1, 10, 1), (0, 2, 8, 1), (1, 3, 1000, 0), (2, 3, 1000, 0)], 4)
    _, dag = a.spf(h, 0, 3, selection=SEL, edge_mask=np.full(4, not masked, dtype=bool))
    fg = c.FlowGraph(g)
    placed, drops = fg.place_with_drops(
        c.FlowIndex(0, 3, 0, 0), 0, 3, dag, offer, LOSSY
    )
    assert placed == 0
    _state(g, fg, [0] * 4, drops, [0] * 4, 0, 3, 0)


def test_reverse_terminal_fanout_does_not_duplicate_ties():
    g, a, h = _graph([(0, 1, 1000, 0), (0, 2, 1000, 0), (1, 3, 10, 1), (2, 3, 8, 1)], 4)
    fan = [int(e) for e, ext in enumerate(g.ext_edge_ids_view()) if ext in (0, 1)]
    for forced in ([], fan):
        d, dag = a.spf_to(h, 3, selection=SEL, fanout_edges=forced)
        assert d.tolist() == [1, 1, 1, 0]
        assert list(g.ext_edge_ids_view()[dag.via_edges]) == [0, 1, 2, 3]
        assert dag.parent_offsets.tolist() == [0, 0, 1, 2, 4]
        assert dag.parents.tolist() == [0, 0, 1, 2]
        fg = c.FlowGraph(g)
        placed, drops = fg.place_with_drops(
            c.FlowIndex(0, 3, 0, 0), 0, 3, dag, 40, LOSSY
        )
        assert placed == 18
        _state(g, fg, [20, 20, 10, 8], drops, [0, 0, 10, 12], 0, 3, 18)


def test_reverse_zero_cycle_still_has_only_one_simple_path():
    g, a, h = _graph([(0, 1, 10, 1), (1, 2, 10, 0), (2, 1, 10, 0), (2, 3, 10, 1)], 4)
    d, dag = a.spf_to(h, 3, selection=SEL)
    assert d.tolist() == [2, 1, 1, 0]
    assert set(map(int, g.ext_edge_ids_view()[dag.via_edges])) == {0, 1, 3}
    assert [[int(u) for u, _ in p] for p in dag.resolve_to_paths(0, 3)] == [
        [0, 1, 2, 3]
    ]


def test_terminal_ksp_spurs_keep_both_ties_with_cost_factor():
    g, a, h = _graph([(0, 1, 10, 1), (0, 2, 8, 1), (1, 3, 1000, 0), (2, 3, 1000, 0)], 4)
    for factor in (None, 1.0):
        items = a.ksp(h, 0, 3, k=3, unique=True, max_cost_factor=factor)
        assert len(items) == 2
        routes = set()
        for d, dag in items:
            assert d[3] == 1
            paths = dag.resolve_to_paths(0, 3)
            assert len(paths) == 1
            routes.add(tuple(int(u) for u, _ in paths[0]))
        assert routes == {(0, 1, 3), (0, 2, 3)}
