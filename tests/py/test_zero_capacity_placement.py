"""Fixed ECMP membership and fresh lossy offers, with hand-derived vectors."""

from __future__ import annotations

import numpy as np
import pytest

import netgraph_core as c

EB = c.FlowPlacement.EQUAL_BALANCED
LOSSY = c.FlowPlacement.EQUAL_BALANCED_LOSSY


def _graph(src, dst, caps, costs, n):
    g = c.StrictMultiDiGraph.from_arrays(
        n,
        np.array(src, dtype=np.int32),
        np.array(dst, dtype=np.int32),
        np.array(caps, dtype=np.float64),
        np.array(costs, dtype=np.int64),
        np.arange(len(src), dtype=np.int64),
    )
    a = c.Algorithms(c.Backend.cpu())
    return g, a, a.build_graph(g)


def _selection(capacity=False):
    return c.EdgeSelection(multi_edge=True, require_capacity=capacity)


def _input_order(g, values):
    return np.asarray(values)[np.argsort(g.ext_edge_ids_view())]


def _drop_vector(g, drops):
    values = np.zeros(g.num_edges())
    for e, amount in drops:
        values[e] += amount
    return _input_order(g, values)


def _state(g, fg, flows, residual):
    # Absolute tolerance also protects the large pseudo-edge residuals: a relative
    # tolerance would hide a substantial ledger error on their million-unit caps.
    assert _input_order(g, fg.edge_flow_view()) == pytest.approx(flows, rel=0, abs=1e-9)
    assert _input_order(g, fg.residual_view()) == pytest.approx(
        residual, rel=0, abs=1e-9
    )
    assert fg.edge_flow_view() + fg.residual_view() == pytest.approx(
        g.capacity_view(), rel=0, abs=1e-9
    )


def _conservation(g, fg, drops, src, dst, delivered):
    # Count a dropped share at the edge's source, where forwarding was attempted.
    balance = np.zeros(g.num_nodes())
    for u, v, amount in zip(
        g.edge_src_view(), g.edge_dst_view(), fg.edge_flow_view(), strict=True
    ):
        balance[u] -= amount
        balance[v] += amount
    lost = 0.0
    for e, amount in drops:
        balance[g.edge_src_view()[e]] -= amount
        lost += amount
    balance[src] += delivered + lost
    balance[dst] -= delivered
    assert balance == pytest.approx(np.zeros(g.num_nodes()), abs=1e-10)


# Input edge order, capacities, costs, offer, lossy delivery, F, D, R.
CASES = [
    (
        [0, 0, 1, 2],
        [1, 2, 3, 3],
        [10, 10, 10, 0],
        [1, 2, 2, 1],
        20,
        10,
        [10, 10, 10, 0],
        [0, 0, 0, 10],
        [0, 0, 0, 0],
    ),
    (
        [0, 0, 1, 2],
        [1, 2, 3, 3],
        [0, 10, 10, 10],
        [1, 2, 2, 1],
        10,
        5,
        [0, 5, 0, 5],
        [5, 0, 0, 0],
        [0, 5, 10, 5],
    ),
    # x receives nothing but must release m, which still has to forward to t.
    (
        [0, 0, 1, 2, 3],
        [1, 2, 3, 3, 4],
        [0, 10, 10, 10, 20],
        [1, 1, 1, 1, 1],
        20,
        10,
        [0, 10, 0, 10, 10],
        [10, 0, 0, 0, 0],
        [0, 0, 10, 0, 10],
    ),
    ([0, 0], [1, 1], [0, 10], [3, 3], 20, 10, [0, 10], [10, 0], [0, 0]),
    ([0, 1], [1, 2], [10, 0], [2, 5], 10, 0, [10, 0], [0, 10], [0, 0]),
    (
        [0, 0, 1, 2],
        [1, 2, 3, 3],
        [0, 0, 10, 10],
        [1, 2, 2, 1],
        20,
        0,
        [0, 0, 0, 0],
        [10, 10, 0, 0],
        [0, 0, 10, 10],
    ),
    (
        [0, 0],
        [1, 1],
        [2**-14, 10],
        [3, 3],
        20,
        10 + 2**-14,
        [2**-14, 10],
        [10 - 2**-14, 0],
        [0, 0],
    ),
]


@pytest.mark.parametrize("mode", [EB, LOSSY])
@pytest.mark.parametrize(
    "src,dst,caps,costs,offer,delivered,flows,drops,residual",
    CASES,
    ids=[
        "dead_end",
        "zero_member",
        "zero_fed_merge",
        "parallel",
        "zero_into_dst",
        "zero_source",
        "below_resolution",
    ],
)
def test_zero_members(
    src, dst, caps, costs, offer, delivered, flows, drops, residual, mode
):
    n = max(src + dst) + 1
    g, a, h = _graph(src, dst, caps, costs, n)
    _, dag = a.spf(h, 0, n - 1, selection=_selection())
    fg = c.FlowGraph(g)
    placed, actual_drops = fg.place_with_drops(
        c.FlowIndex(0, n - 1, 0, 0), 0, n - 1, dag, offer, mode
    )
    if mode == EB:
        delivered, flows, drops, residual = 0, [0] * len(caps), [0] * len(caps), caps
    assert placed == pytest.approx(delivered)
    assert _drop_vector(g, actual_drops) == pytest.approx(drops)
    _state(g, fg, flows, residual)
    _conservation(g, fg, actual_drops, 0, n - 1, placed)


@pytest.mark.parametrize("mode", [EB, LOSSY])
@pytest.mark.parametrize("last_capacity", [0, 10])
def test_source_is_also_transit(mode, last_capacity):
    # Cost-compensated root gives equal origination at a and b; a also transits b.
    g, a, h = _graph(
        [0, 0, 1, 2], [1, 2, 2, 3], [100, 100, 10, last_capacity], [1, 2, 1, 1], 4
    )
    _, dag = a.spf(h, 0, 3, selection=_selection())
    fg = c.FlowGraph(g)
    placed, drops = fg.place_with_drops(c.FlowIndex(0, 3, 0, 0), 0, 3, dag, 20, mode)
    expected = 0 if last_capacity == 0 else 10
    assert placed == expected
    if mode == LOSSY:
        _state(g, fg, [10, 10, 10, last_capacity], [90, 90, 0, 0])
        assert _drop_vector(g, drops) == pytest.approx([0, 0, 0, 20 - last_capacity])
    elif last_capacity:
        _state(g, fg, [5, 5, 5, 10], [95, 95, 5, 0])
        assert drops == []
    else:
        _state(g, fg, [0, 0, 0, 0], [100, 100, 10, 0])
        assert drops == []
    _conservation(g, fg, drops, 0, 3, placed)


@pytest.mark.parametrize("cleanup", ["remove", "reset"])
def test_zero_delivery_still_owns_upstream_flow(cleanup):
    g, a, h = _graph([0, 1], [1, 2], [10, 0], [2, 5], 3)
    _, dag = a.spf(h, 0, 2, selection=_selection())
    fg = c.FlowGraph(g)
    idx = c.FlowIndex(0, 2, 0, 0)
    placed, drops = fg.place_with_drops(idx, 0, 2, dag, 10, LOSSY)
    assert placed == 0
    assert _drop_vector(g, drops).tolist() == [0, 10]
    assert [(int(g.ext_edge_ids_view()[e]), v) for e, v in fg.get_flow_edges(idx)] == [
        (0, 10)
    ]
    _state(g, fg, [10, 0], [0, 0])
    fg.remove(idx) if cleanup == "remove" else fg.reset()
    assert fg.get_flow_edges(idx) == []
    _state(g, fg, [0, 0], [10, 0])


@pytest.mark.parametrize("mode", [EB, LOSSY, c.FlowPlacement.PROPORTIONAL])
@pytest.mark.parametrize("route", ["masked", "capacity", "zero_offer", "fully_masked"])
def test_excluded_members_and_no_offer(mode, route):
    g, a, h = _graph([0, 0], [1, 1], [0, 10], [3, 3], 2)
    mask = np.array([False, route != "fully_masked"], dtype=bool)
    _, dag = a.spf(
        h,
        0,
        1,
        selection=_selection(route == "capacity"),
        edge_mask=None if route in ("capacity", "zero_offer") else mask,
    )
    offer = 0 if route == "zero_offer" else 10
    fg = c.FlowGraph(g)
    placed, drops = fg.place_with_drops(c.FlowIndex(0, 1, 0, 0), 0, 1, dag, offer, mode)
    expected = 0 if route in ("zero_offer", "fully_masked") else 10
    assert placed == expected
    assert drops == []
    _state(g, fg, [0, expected], [0, 10 - expected])


@pytest.mark.parametrize(
    "mode,delivered,flows,drops",
    [(EB, 0, [10, 0], [0, 0]), (LOSSY, 5, [10, 5], [5, 0])],
)
def test_residual_zero_is_a_fixed_member(mode, delivered, flows, drops):
    g, a, h = _graph([0, 0], [1, 1], [10, 10], [3, 3], 2)
    fg = c.FlowGraph(g)
    idx = c.FlowIndex(0, 1, 0, 0)
    first = c.PredDAG.from_edges(g, [0])
    assert fg.place(idx, 0, 1, first, 10, EB) == 10
    _, dag = a.spf(h, 0, 1, selection=_selection())
    placed, actual = fg.place_with_drops(c.FlowIndex(0, 1, 0, 1), 0, 1, dag, 10, mode)
    assert placed == delivered
    assert _drop_vector(g, actual) == pytest.approx(drops)
    _state(g, fg, flows, [0, 10 - flows[1]])


@pytest.mark.parametrize("shortest", [True, False])
def test_cost_only_analyses_and_masked_control(shortest):
    g, a, h = _graph([0, 0, 1, 2], [1, 2, 3, 3], [10, 10, 10, 0], [1, 2, 2, 1], 4)
    options = dict(flow_placement=EB, shortest_path=shortest, require_capacity=False)
    total, summary = a.max_flow(
        h, 0, 3, with_edge_flows=True, with_residuals=True, **options
    )
    assert total == summary.total_flow == 0
    assert _input_order(g, summary.edge_flows).tolist() == [0, 0, 0, 0]
    assert _input_order(g, summary.residual_capacity).tolist() == [10, 10, 10, 0]
    assert list(summary.costs) == list(summary.flows) == []
    assert a.sensitivity_analysis(h, 0, 3, **options) == []
    assert c.FlowState(g).place_max_flow(0, 3, **options) == 0
    batch = a.batch_max_flow(h, np.array([[0, 3]], dtype=np.int32), **options)
    assert [s.total_flow for s in batch] == [0]
    # Removing b makes the surviving cost-3 route carry 10; capacity-aware
    # discovery and proportional placement also have a physical cut of 10.
    for controls in [
        dict(node_mask=np.array([True, True, False, True])),
        dict(require_capacity=True),
        dict(flow_placement=c.FlowPlacement.PROPORTIONAL),
    ]:
        opts = options | controls
        value, result = a.max_flow(
            h, 0, 3, with_edge_flows=True, with_residuals=True, **opts
        )
        assert value == result.total_flow == 10
        assert _input_order(g, result.edge_flows).tolist() == [10, 0, 10, 0]
        assert _input_order(g, result.residual_capacity).tolist() == [0, 10, 0, 0]
        assert list(result.costs) == [3]
        assert list(result.flows) == [10]
        assert sum(g.capacity_view()[result.min_cut.edges]) == 10


@pytest.mark.parametrize(
    "case,offer,total,flows,drops",
    [
        ("healthy", 49, 48, [49 / 6] * 6 + [24, 24], [0] * 6 + [0.5, 0.5]),
        ("deadwan", 48, 24, [8] * 6 + [24, 0], [0] * 7 + [24]),
        ("deadpod", 60, 40, [10, 10, 10, 10, 0, 0, 20, 20], [0, 0, 0, 0, 10, 10, 0, 0]),
        ("deadwan1", 20, 10, [10, 10, 0, 0, 0, 0, 10, 0], [0] * 7 + [10]),
        ("partial", 30, 25, [5, 0, 5, 5, 5, 5, 15, 10], [0, 5, 0, 0, 0, 0, 0, 0]),
        (
            "asymmetric",
            60,
            49,
            [10, 5, 10, 10, 6, 10, 24, 25],
            [0, 5, 0, 0, 4, 0, 2, 0],
        ),
    ],
)
@pytest.mark.parametrize("mode", [EB, LOSSY])
def test_walkthrough_vectors(case, offer, total, flows, drops, mode):
    caps = [10] * 6 + [24, 24]
    if case in ("deadwan", "deadwan1"):
        caps[7] = 0
    elif case == "deadpod":
        caps[4:6] = [0, 0]
    elif case == "partial":
        caps[1] = 0
    elif case == "asymmetric":
        caps = [10, 5, 10, 10, 6, 12, 24, 30]
    src, dst = [0, 0, 1, 1, 2, 2, 3, 4], [3, 4, 3, 4, 3, 4, 5, 5]
    costs = [1] * 8
    root = 0 if case == "deadwan1" else 6
    if root == 6:
        src += [6, 6, 6]
        dst += [0, 1, 2]
        caps += [1e6] * 3
        costs += [0] * 3
    g, a, h = _graph(src, dst, caps, costs, 7)
    if root == 6:
        fanout = [i for i, ext in enumerate(g.ext_edge_ids_view()) if ext >= 8]
        _, dag = a.spf_to(h, 5, selection=_selection(), fanout_edges=fanout)
    else:
        _, dag = a.spf(h, root, 5, selection=_selection())
    if mode == EB:
        total = 48 if case == "healthy" else (30 if case == "asymmetric" else 0)
        flows = (
            [8] * 6 + [24, 24]
            if case == "healthy"
            else [5] * 6 + [15, 15]
            if case == "asymmetric"
            else [0] * 8
        )
        drops = [0] * 8
    if root == 6:
        flows = flows + [total / 3 if mode == EB else offer / 3] * 3
        drops = drops + [0] * 3
    residual = [capacity - flow for capacity, flow in zip(caps, flows, strict=True)]
    fg = c.FlowGraph(g)
    placed, actual = fg.place_with_drops(
        c.FlowIndex(root, 5, 0, 0), root, 5, dag, 90 if mode == EB else offer, mode
    )
    assert placed == pytest.approx(total)
    _state(g, fg, flows, residual)
    assert _drop_vector(g, actual) == pytest.approx(drops)
    _conservation(g, fg, actual, root, 5, placed)


def _policy(a, h, *, mode=LOSSY, minimum=1, maximum=1, shortest=True):
    cfg = c.FlowPolicyConfig(
        flow_placement=mode,
        selection=_selection(),
        require_capacity=False,
        shortest_path=shortest,
        min_flow_count=minimum,
        max_flow_count=maximum,
    )
    return c.FlowPolicy(a, h, cfg)


@pytest.mark.parametrize(
    "caps,expected", [([0, 10], 10), ([1, 10], 11), ([100, 100], 20)]
)
@pytest.mark.parametrize("minimum,maximum", [(1, 1), (2, 2), (1, 3), (2, 4), (1, None)])
def test_policy_full_and_repeated_offers(caps, expected, minimum, maximum):
    g, a, h = _graph([0, 0], [1, 1], caps, [3, 3], 2)
    fg = c.FlowGraph(g)
    p = _policy(a, h, minimum=minimum, maximum=maximum)
    for call in (1, 2):
        total = expected if call == 1 or caps == [100, 100] else 0
        assert p.place_demand(fg, 0, 1, 0, 20) == pytest.approx((total, 20 - total))
        flows = [10 * call] * 2 if caps == [100, 100] else caps
        _state(
            g,
            fg,
            flows,
            [capacity - flow for capacity, flow in zip(caps, flows, strict=True)],
        )
        assert p.flow_count() == minimum
    p.remove_demand(fg)
    _state(g, fg, [0, 0], caps)


def test_policy_explicit_target_remains_cumulative():
    g, a, h = _graph([0, 0], [1, 1], [100, 100], [3, 3], 2)
    fg = c.FlowGraph(g)
    p = _policy(a, h)
    assert p.place_demand(fg, 0, 1, 0, 20, target_per_flow=10) == (10, 10)
    assert p.place_demand(fg, 0, 1, 0, 20, target_per_flow=10) == (0, 20)
    _state(g, fg, [5, 5], [95, 95])


@pytest.mark.parametrize("offer", [2**-12, 2**-13])
def test_policy_offers_small_shares_to_every_seeded_record(offer):
    # A share below the lossless admission threshold is still a lossy offer.
    g, a, h = _graph([0, 0], [1, 1], [1, 1], [3, 3], 2)
    p = _policy(a, h, minimum=2, maximum=2)
    fg = c.FlowGraph(g)
    assert p.place_demand(fg, 0, 1, 0, offer) == (offer, 0)
    _state(g, fg, [offer / 2, offer / 2], [1 - offer / 2, 1 - offer / 2])


@pytest.mark.parametrize("caps", [[1, 10], [100, 100]])
def test_direct_repeated_call_drop_vectors(caps):
    g, a, h = _graph([0, 0], [1, 1], caps, [3, 3], 2)
    _, dag = a.spf(h, 0, 1, selection=_selection())
    fg = c.FlowGraph(g)
    for call in (1, 2):
        placed, drops = fg.place_with_drops(
            c.FlowIndex(0, 1, 0, call), 0, 1, dag, 20, LOSSY
        )
        if caps == [100, 100]:
            assert placed == 20
            assert drops == []
            _state(g, fg, [10 * call, 10 * call], [100 - 10 * call, 100 - 10 * call])
        else:
            assert placed == (11 if call == 1 else 0)
            assert _drop_vector(g, drops).tolist() == (
                [9, 0] if call == 1 else [10, 10]
            )
            _state(g, fg, [1, 10], [0, 0])


def test_pinned_repeats_and_explicit_targets():
    # Two independent one-edge bundles, sharing each new offer at the head.
    g, a, h = _graph([0, 0], [1, 1], [100, 100], [3, 3], 2)
    for explicit in (None, 10):
        p = _policy(a, h, maximum=2, shortest=False)
        p.set_static_paths(
            0, 1, [c.PredDAG.from_edges(g, [0]), c.PredDAG.from_edges(g, [1])]
        )
        fg = c.FlowGraph(g)
        for call in (1, 2):
            placed = 0 if explicit is not None and call == 2 else 20
            assert p.place_demand(fg, 0, 1, 0, 20, target_per_flow=explicit) == (
                placed,
                20 - placed,
            )
            carried = 10 if explicit is not None else 10 * call
            _state(g, fg, [carried, carried], [100 - carried, 100 - carried])


@pytest.mark.parametrize("caps,total", [([0, 10], 10), ([1, 10], 11)])
@pytest.mark.parametrize("pinned", [False, True])
def test_policy_unlimited_offer(caps, total, pinned):
    g, a, h = _graph([0, 0], [1, 1], caps, [3, 3], 2)
    p = _policy(a, h, shortest=not pinned)
    if pinned:
        _, dag = a.spf(h, 0, 1, selection=_selection())
        p.set_static_paths(0, 1, [dag])
    fg = c.FlowGraph(g)
    placed, remaining = p.place_demand(fg, 0, 1, 0, float("inf"))
    assert placed == total
    assert np.isposinf(remaining)
    _state(g, fg, caps, [0, 0])
    direct = c.FlowGraph(g)
    _, dag = a.spf(h, 0, 1, selection=_selection())
    assert direct.place_with_drops(
        c.FlowIndex(0, 1, 0, 0), 0, 1, dag, float("inf"), LOSSY
    ) == (total, [])
    _state(g, direct, caps, [0, 0])


@pytest.mark.parametrize(
    "low,offer,total", [(1, 20, 11), (10, 50, 35), (1, float("inf"), 11)]
)
@pytest.mark.parametrize("masked", [False, True])
def test_policy_pinned_shares(low, offer, total, masked):
    high = 10 if low == 1 else 100
    g, a, h = _graph(
        [0, 1, 0, 2], [1, 3, 2, 3], [low, low, high, high], [1, 1, 2, 2], 4
    )
    ids = {int(ext): i for i, ext in enumerate(g.ext_edge_ids_view())}
    bundles = [
        c.PredDAG.from_edges(g, [ids[0], ids[1]]),
        c.PredDAG.from_edges(g, [ids[2], ids[3]]),
    ]
    cfg = c.FlowPolicyConfig(
        flow_placement=LOSSY,
        selection=_selection(),
        require_capacity=False,
        shortest_path=False,
        min_flow_count=1,
        max_flow_count=2,
    )
    mask = np.ones(g.num_edges(), dtype=bool)
    if masked:
        mask[ids[0]] = False
    p = c.FlowPolicy(a, h, cfg, edge_mask=mask)
    p.set_static_paths(0, 3, bundles)
    fg = c.FlowGraph(g)
    expected = min(offer, high) if masked else total
    assert p.place_demand(fg, 0, 3, 0, offer) == pytest.approx(
        (expected, offer - expected)
    )
    carried_high = min(offer if masked else offer / 2, high)
    flows = (
        [0, 0, carried_high, carried_high]
        if masked
        else [low, low, carried_high, carried_high]
    )
    _state(
        g,
        fg,
        flows,
        [
            capacity - flow
            for capacity, flow in zip([low, low, high, high], flows, strict=True)
        ],
    )
    assert p.flow_count() == (1 if masked else 2)


@pytest.mark.parametrize("pinned", [False, True])
@pytest.mark.parametrize("offer,masked", [(0, False), (20, True)])
def test_policy_zero_offer_and_everything_masked(pinned, offer, masked):
    g, a, h = _graph([0, 0], [1, 1], [0, 10], [3, 3], 2)
    cfg = c.FlowPolicyConfig(
        flow_placement=LOSSY,
        selection=_selection(),
        require_capacity=False,
        shortest_path=not pinned,
        min_flow_count=1,
        max_flow_count=2,
    )
    p = c.FlowPolicy(
        a, h, cfg, edge_mask=np.full(g.num_edges(), not masked, dtype=bool)
    )
    if pinned:
        p.set_static_paths(
            0, 1, [c.PredDAG.from_edges(g, [0]), c.PredDAG.from_edges(g, [1])]
        )
    fg = c.FlowGraph(g)
    assert p.place_demand(fg, 0, 1, 0, offer) == (0, offer)
    _state(g, fg, [0, 0], [0, 10])


def test_policy_lossless_and_cost_factor_controls():
    g, a, h = _graph([0, 0, 0, 2], [1, 1, 2, 1], [1, 10, 100, 100], [3, 3, 4, 4], 3)
    for mode, delivered, flows in [(EB, 2, [1, 1, 0, 0]), (LOSSY, 11, [1, 10, 0, 0])]:
        cfg = c.FlowPolicyConfig(
            flow_placement=mode,
            selection=_selection(),
            require_capacity=False,
            shortest_path=True,
            min_flow_count=1,
            max_flow_count=1,
            max_path_cost_factor=1.5,
        )
        p = c.FlowPolicy(a, h, cfg)
        fg = c.FlowGraph(g)
        assert p.place_demand(fg, 0, 1, 0, 20) == pytest.approx(
            (delivered, 20 - delivered)
        )
        _state(
            g,
            fg,
            flows,
            [
                capacity - flow
                for capacity, flow in zip([1, 10, 100, 100], flows, strict=True)
            ],
        )
        assert {v[2] for v in p.flows.values()} == {3}
