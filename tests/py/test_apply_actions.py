"""apply_actions: several policies of one ledger in one operation (docs/semantics.md, "Several
demands at once"). The action values and the call's conversions (the types first, in argument
order, then the list's checks in list order; a refusal changes nothing), ORDERED as the same
scalar calls on a private working state, SIMULTANEOUS on one view with the losses shared across
the demands, per-entry outputs, the all-or-nothing sweeps of a list's Python allocations and the
forwarding size limits counted over the whole list. Hand-derived values on small graphs (EdgeIds
equal input positions unless mapped); two seeded generators compare an ORDERED list with the same
scalar calls, and a SIMULTANEOUS participant with its scalar call, on an identical world."""

from __future__ import annotations

import gc
import random
import weakref
from fractions import Fraction

import _netgraph_core as ext
import numpy as np
import pytest
from helpers import (
    INF,
    L62,
    assert_read_only,
    build,
    centralized,
    config,
    dag_of,
    exact,
    full_ledger_state,
    ledger_state,
    plain,
    shortest,
    sweep_python_allocation,
    t1_config,
)
from private_hooks import HOOKS, needs_hooks

import netgraph_core as ngc

LY = ngc.Admission.LOSSY
LL = ngc.Admission.LOSSLESS
ORDERED = ngc.BatchMode.ORDERED
SIMULTANEOUS = ngc.BatchMode.SIMULTANEOUS
W62 = L62 - 1
OFFER = "offer must be an integer in [0, 2^62) or permitted +inf"
PERIOD_TWO = "lossy forwarding iteration entered a period-2 cycle"
ROUND_LIMIT = (
    "lossy forwarding iteration did not reach a fixed point within the round limit"
)
ALL = ngc.OutputOptions(
    edge_flows=True,
    residual=True,
    saturated_boundary_edges=True,
    edge_reservations=True,
)
UP = ngc.TunnelStatus.UP
BLACKHOLED = ngc.TunnelStatus.BLACKHOLED


def graph(n, rows):
    """A Graph from (u, v, capacity, cost) rows written in EdgeId order (checked)."""
    g = build(n, rows)
    assert [int(i) for i in g.input_positions()] == list(range(len(rows)))
    return g


def mapped(n, rows):
    """A Graph from rows in any order and the EdgeId of each input row."""
    g = build(n, rows)
    ids = [0] * len(rows)
    for e, p in enumerate(g.input_positions()):
        ids[int(p)] = e
    return g, ids


def ecmp(admission=LY):
    return config(ngc.ShortestPaths(), admission=admission)


def t1(admission=LY):
    return t1_config(1, admission)


def t2(admission=LY):
    return config(centralized(tunnel_count=ngc.TunnelCount(1)), admission=admission)


def controller(admission=LL):
    return config(centralized(), admission=admission)


def tunnels(g, sources, destinations, specs):
    """A supplied program of tunnels: specs are (head, [EdgeIds in path order], weight)."""
    routes = [
        ngc.RouteDefinition(h, ngc.Path(g, [[e] for e in es]), weight=w)
        for h, es, w in specs
    ]
    return ngc.ForwardingProgram.from_routes(g, sources, destinations, routes)


def tunnel(g, s, d, edges, weight=1):
    return tunnels(g, [s], [d], [(s, edges, weight)])


def no_edges(g, edges):
    return g.eligibility(excluded_edges=edges)


def no_nodes(g, nodes):
    return g.eligibility(excluded_nodes=nodes)


def policy_state(p):
    """A policy's state without object identities: its delivery, its records' local ids and
    flows, its program and (private hooks) its epoch."""
    out = [
        p.delivered(),
        [(r.flow_id[1], r.delivered, r.edge_flows, r.edge_drops) for r in p.records()],
        plain(p.program),
    ]
    if HOOKS:
        epoch = p._state()
        for flows in epoch["cohorts"].values():
            for flow in flows:
                flow["flow_id"] = flow["flow_id"][1]
        out.append(epoch)
    return out


class World:
    """A ledger with its policies (created in order), each policy's sources and mask."""

    def __init__(self, g, service=None, reservation=None):
        self.g = g
        self.ledger = ngc.FlowLedger(g, service, reservation_baseline=reservation)
        self.p = []
        self.sources = []
        self.masks = []

    def add(self, sources, destinations, cfg, mask=None):
        self.p.append(ngc.FlowPolicy(self.ledger, sources, destinations, cfg, mask))
        self.sources.append(sources)
        self.masks.append(mask)
        return self.p[-1]

    def pin(self, program, admission=LY, black_holes=None, mask=None):
        cfg = ngc.FlowPolicyConfig(
            routing=ngc.Pinned(program, black_holes=black_holes), admission=admission
        )
        return self.add(list(program.sources), list(program.destinations), cfg, mask)

    def state(self):
        """The ledger (private hooks: with its next FlowId and owners) and every policy."""
        ledger = full_ledger_state(self.ledger) if HOOKS else ledger_state(self.ledger)
        return [ledger, [policy_state(p) for p in self.p]]


def run(w, actions, mode=ORDERED, outputs=None):
    out = ngc.apply_actions(w.ledger, actions, mode=mode, outputs=outputs)
    assert type(out) is list and len(out) == len(actions)
    assert all(type(r) is ngc.PlacementResult for r in out)
    return out


def refused(w, actions, mode=ORDERED, outputs=None, error=ValueError):
    """The message of the error the list raises (exactly that class); nothing changed."""
    before = w.state()
    with pytest.raises(error) as info:
        ngc.apply_actions(w.ledger, actions, mode=mode, outputs=outputs)
    assert type(info.value) is error
    assert w.state() == before
    return str(info.value)


def delivered(results):
    return [r.call_delivered for r in results]


def cyc12():
    """The cycle Pc's graph: CYC(1, 2)."""
    return graph(
        6,
        [
            (0, 1, 1, 1),
            (1, 2, 1000, 1),
            (1, 5, 1000, 1),
            (2, 3, 2, 1),
            (3, 0, 1000, 1),
            (3, 5, 1000, 1),
            (4, 0, 1000, 1),
            (4, 2, 1000, 1),
        ],
    )


def pc_program(g):
    """Two tunnels crossing X (e0) and Y (e3) in opposite orders, weights 1 and 2: replace(3)
    enters a period-2 cycle."""
    return tunnels(g, [4], [5], [(4, [6, 0, 1, 3, 5], 1), (4, [7, 3, 4, 0, 2], 2)])


def g1():
    return graph(2, [(0, 1, 10, 1)])


# -- the action values -------------------------------------------------------------------------


def test_action_fields():
    g = g1()
    p = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 1, ecmp())
    m = no_edges(g, [0])
    place, rep, fwd = (
        ngc.Place(p, 3),
        ngc.Replace(p, 3, design_demand=5),
        ngc.Forward(p, 3, black_holes=m),
    )
    assert place.policy is p and rep.policy is p and fwd.policy is p
    assert [place.offer, rep.offer, fwd.offer, rep.design_demand] == [3, 3, 3, 5]
    assert all(type(x) is int for x in (place.offer, rep.offer, rep.design_demand))
    assert fwd.black_holes is m
    assert ngc.Place(p, INF).offer == INF
    assert ngc.Replace(p, 3).design_demand is None
    assert ngc.Forward(p, 3).black_holes is None
    # Kept as given, never converted to units.
    assert ngc.Place(p, 2**70).offer == 2**70
    assert ngc.Replace(p, 1.5, design_demand=2.5).design_demand == 2.5
    assert_read_only(place, "policy", "offer")
    assert_read_only(rep, "policy", "offer", "design_demand")
    assert_read_only(fwd, "policy", "offer", "black_holes")


def test_signatures():
    g = g1()
    led = ngc.FlowLedger(g)
    p = ngc.FlowPolicy(led, 0, 1, ecmp())
    m = no_edges(g, [0])
    for build_one in (
        lambda: ngc.Replace(p, 1, 5),
        lambda: ngc.Forward(p, 1, m),
        lambda: ngc.apply_actions(led, [], ORDERED),
        lambda: ngc.apply_actions(led, [], ORDERED, None),
    ):
        with pytest.raises(TypeError, match="incompatible"):
            build_one()
    out = ngc.apply_actions(
        ledger=led, actions=[ngc.Replace(p, 1)], mode=ORDERED, outputs=None
    )
    assert type(out) is list and len(out) == 1
    assert ngc.apply_actions.__doc__.startswith(
        "apply_actions(ledger: object, actions: object, *, mode: object = BatchMode.ORDERED,"
        " outputs: object = None) -> list[_netgraph_core.PlacementResult]"
    )


def test_action_construction_errors():
    g = g1()
    p = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 1, ecmp())
    cases = [
        (lambda: ngc.Place(None, 1), "policy must be a FlowPolicy (got NoneType)"),
        (lambda: ngc.Place(3, 1), "policy must be a FlowPolicy (got int)"),
        (lambda: ngc.Place(p, True), "offer must be an int or a float (got bool)"),
        (lambda: ngc.Place(p, "1"), "offer must be an int or a float (got str)"),
        (lambda: ngc.Place(p, None), "offer must be an int or a float (got NoneType)"),
        (
            lambda: ngc.Place(p, Fraction(3, 2)),
            "offer must be an int or a float (got Fraction)",
        ),
        (
            lambda: ngc.Place(p, np.bool_(True)),
            "offer must be an int or a float (got bool)",
        ),
        (
            lambda: ngc.Replace(p, 1, design_demand=True),
            "design_demand must be an int or a float (got bool)",
        ),
        (
            lambda: ngc.Replace(p, 1, design_demand="x"),
            "design_demand must be an int or a float (got str)",
        ),
        (
            lambda: ngc.Forward(p, 1, black_holes=3),
            "black_holes must be an Eligibility or None (got int)",
        ),
        (
            lambda: ngc.Forward(p, True, black_holes=3),
            "offer must be an int or a float (got bool)",
        ),
        (
            lambda: ngc.Replace(None, True),
            "policy must be a FlowPolicy (got NoneType)",
        ),
    ]
    for make, text in cases:
        with pytest.raises(TypeError, match=exact(text)) as info:
            make()
        assert type(info.value) is TypeError
    with pytest.raises(TypeError, match="incompatible"):
        ngc.Place(p, 1, design_demand=1)
    with pytest.raises(TypeError, match="incompatible"):
        ngc.Place(p)
    # NumPy integers are ints (kept as Python ints), NumPy floats are floats.
    three = ngc.Place(p, np.int64(3)).offer
    assert three == 3 and type(three) is int
    assert ngc.Place(p, np.uint8(3)).offer == 3
    assert ngc.Place(p, np.float32(1.5)).offer == 1.5
    assert ngc.Place(p, np.float64(INF)).offer == INF

    class Index:
        def __index__(self):
            return 4

    assert type(ngc.Place(p, Index()).offer) is int


def test_offers_are_checked_by_the_call():
    g = g1()
    w = World(g)
    p = w.add(0, 1, ecmp())
    w.add(0, 1, ecmp()).place(2)
    for offer in (-1, 2**62, 2**63, 2**70, 1.5, np.float32(1.5), float("nan"), -INF):
        action = ngc.Place(p, offer)  # constructs
        assert refused(w, [action]) == OFFER
        assert refused(w, [ngc.Replace(p, offer)], SIMULTANEOUS) == OFFER
    for design in (-1, 2**63, 1.5, float("nan")):
        t = w.add(0, 1, t1())
        text = "design_demand must be an integer in [0, 2^62) or permitted +inf"
        assert refused(w, [ngc.Replace(t, 1, design_demand=design)]) == text
    # +inf (a NumPy infinity too) is +inf, never a bad finite value.
    f = w.add(0, 1, ecmp(LL))
    assert run(w, [ngc.Replace(f, np.float64(INF))])[0].requested == INF


def test_list_conversion_errors():
    g = g1()
    w = World(g)
    e = w.add(0, 1, ecmp())
    other = World(g1())
    x = other.add(0, 1, ecmp())
    e.place(2)
    cases = [
        ((None, []), {}, "ledger must be a FlowLedger"),
        ((5, []), {}, "ledger must be a FlowLedger"),
        ((None, 5), {}, "ledger must be a FlowLedger"),
        (
            (w.ledger, 5),
            {},
            "actions must be a sequence of Place, Replace or Forward (got int)",
        ),
        (
            (w.ledger, "ab"),
            {},
            "actions must be a sequence of Place, Replace or Forward (got str)",
        ),
        (
            (w.ledger, [ngc.Place(e, 1), 1]),
            {},
            "actions[1] must be a Place, Replace or Forward (got int)",
        ),
        # Before any native check: the foreign first action, the bad mode.
        (
            (w.ledger, [ngc.Place(x, 1), None]),
            {"mode": ngc.BatchMode(7)},
            "actions[1] must be a Place, Replace or Forward (got NoneType)",
        ),
        ((w.ledger, []), {"mode": 0}, "mode must be a BatchMode (got int)"),
        ((w.ledger, []), {"mode": None}, "mode must be a BatchMode (got NoneType)"),
        (
            (w.ledger, [ngc.Place(x, 1)]),
            {"outputs": 1},
            "outputs must be an OutputOptions (got int)",
        ),
    ]
    for args, kw, text in cases:
        before = w.state()
        with pytest.raises(TypeError, match=exact(text)) as info:
            ngc.apply_actions(*args, **kw)
        assert type(info.value) is TypeError
        assert w.state() == before
    # A tuple, a NumPy object array and a generator-free sequence are sequences.
    assert delivered(run(w, (ngc.Replace(e, 3),))) == [3]
    assert delivered(run(w, np.array([ngc.Replace(e, 4)], dtype=object))) == [4]


def test_actions_keep_their_objects_alive():
    g = graph(2, [(0, 1, 10, 1), (0, 1, 10, 1)])
    led = ngc.FlowLedger(g)
    action = ngc.Replace(ngc.FlowPolicy(led, 0, 1, ecmp()), 4)
    gc.collect()
    assert delivered(ngc.apply_actions(led, [action])) == [4]
    assert led.carried().tolist() == [2, 2]
    held = ngc.FlowPolicy(led, 0, 1, ecmp())
    held.replace(2)
    forward = ngc.Forward(held, 6, black_holes=g.eligibility(excluded_edges=[1]))
    alive = weakref.ref(forward.black_holes)
    del held
    gc.collect()
    assert alive() is not None
    r = ngc.apply_actions(led, [forward])[0]
    assert r.call_delivered == 3 and r.call_edge_drops == [(1, 3)]
    assert len(led.records()) == 2


# -- validation --------------------------------------------------------------------------------


def test_empty_list():
    w = World(g1())
    e = w.add(0, 1, ecmp())
    e.place(3)
    before = w.state()
    assert run(w, []) == [] and run(w, [], SIMULTANEOUS) == []
    assert w.state() == before
    text = "mode: 7 is not a value of BatchMode"
    assert refused(w, [], ngc.BatchMode(7)) == text
    assert refused(w, [ngc.Place(e, 1)], ngc.BatchMode(7)) == text


def test_validation_order():
    w = World(g1())
    e = w.add(0, 1, ecmp())
    e2 = w.add(0, 1, ecmp())
    f = w.add(0, 1, ecmp(LL))
    t = w.add(0, 1, t1())
    k = w.add(0, 1, controller())
    opt = w.add(0, 1, config(ngc.ShortestPaths(), ngc.RouteAllocation.OPTIMIZED))
    tt = w.add(0, 1, t2())
    p = w.pin(tunnel(w.g, 0, 1, [0]))
    pll = w.pin(tunnel(w.g, 0, 1, [0]), LL)
    inf = w.add(0, 1, ecmp())
    other = World(g1())
    x = other.add(0, 1, ecmp())
    inf.place(INF)
    k.place(1)
    e.place(2)
    P, R, F = ngc.Place, ngc.Replace, ngc.Forward
    model = "SIMULTANEOUS Replace requires ShortestPaths EQUAL or WEIGHTED, or Pinned"
    lossy = "SIMULTANEOUS requires admission=LOSSY"
    cases = [
        # Ownership, then the duplicate rule (at the second occurrence).
        ([P(e, 1), P(e, 2)], ORDERED, "duplicate policy in list"),
        ([R(e, 1), F(e, 1)], ORDERED, "duplicate policy in list"),
        ([R(e, 1), F(e, 1)], SIMULTANEOUS, "duplicate policy in list"),
        ([P(x, 1)], ORDERED, "list policy belongs to another ledger"),
        ([P(x, 1), P(x, 1)], ORDERED, "list policy belongs to another ledger"),
        ([R(x, 1)], SIMULTANEOUS, "list policy belongs to another ledger"),
        ([P(e, 1), P(e, -1)], ORDERED, "duplicate policy in list"),
        # Every check of an action before the next action's.
        ([P(e, -1), P(e2, 1)], ORDERED, OFFER),
        ([P(e, -1), P(x, 1)], ORDERED, OFFER),
        ([P(e, -1), P(e, 1)], ORDERED, OFFER),
        # SIMULTANEOUS's rules in their order, then the scalar checks.
        ([P(e, 1)], SIMULTANEOUS, "SIMULTANEOUS does not accept Place"),
        ([R(t, 1)], SIMULTANEOUS, model),
        ([R(k, 1)], SIMULTANEOUS, model),
        ([R(opt, 1)], SIMULTANEOUS, model),
        ([R(tt, 1)], SIMULTANEOUS, model),
        ([R(f, 1)], SIMULTANEOUS, lossy),
        ([R(pll, 1)], SIMULTANEOUS, lossy),
        ([F(k, 1)], SIMULTANEOUS, lossy),
        ([R(e, INF)], SIMULTANEOUS, "SIMULTANEOUS requires finite offers"),
        ([F(p, 1)], SIMULTANEOUS, "Pinned forwards with replace"),
        ([F(e2, 1)], SIMULTANEOUS, "forward requires an installed program"),
        ([R(e, -1)], SIMULTANEOUS, OFFER),
        (
            [R(e, 1, design_demand=1)],
            SIMULTANEOUS,
            "ShortestPaths does not accept design_demand",
        ),
        # Several invalid parts in one action.
        ([F(e2, INF)], SIMULTANEOUS, "SIMULTANEOUS requires finite offers"),
        ([F(k, INF)], SIMULTANEOUS, lossy),
        ([R(t, -1)], SIMULTANEOUS, model),
        (
            [R(e2, 1), P(e, 1), R(t, 1)],
            SIMULTANEOUS,
            "SIMULTANEOUS does not accept Place",
        ),
        # ORDERED runs the scalar checks in the call's order.
        ([P(t, 1)], ORDERED, "DistributedTE does not support incremental placement"),
        (
            [P(tt, 0)],
            ORDERED,
            "CentralizedTE with tunnel_count does not support incremental placement",
        ),
        (
            [F(e, 1, black_holes=no_edges(other.g, [0]))],
            ORDERED,
            "black_holes belongs to another graph",
        ),
        (
            [P(inf, 1)],
            ORDERED,
            "finite placement after infinite placement requires a new epoch",
        ),
        ([P(t, -1)], ORDERED, "DistributedTE does not support incremental placement"),
        ([F(p, 1)], ORDERED, "Pinned forwards with replace"),
        ([F(e2, 1)], ORDERED, "forward requires an installed program"),
        (
            [R(f, 1, design_demand=2)],
            ORDERED,
            "ShortestPaths does not accept design_demand",
        ),
        (
            [R(t, 1, design_demand=INF)],
            ORDERED,
            "DistributedTE requires a finite design demand",
        ),
    ]
    for actions, mode, text in cases:
        assert refused(w, actions, mode) == text, (actions, mode)


def test_invalid_last_action_changes_nothing():
    w = World(cyc12())
    pc = w.pin(pc_program(w.g))
    t = w.add(4, 5, t1())
    fresh = [w.add(4, 5, ecmp()) for _ in range(48)]
    e = w.add(4, 5, ecmp())
    ordered = [ngc.Replace(pc, 3)] + [ngc.Place(f, 1) for f in fresh]
    simultaneous = [ngc.Replace(pc, 3)] + [ngc.Replace(f, 1) for f in fresh]
    assert (
        refused(w, [*ordered, ngc.Place(t, 1)])
        == "DistributedTE does not support incremental placement"
    )
    assert (
        refused(w, [*simultaneous, ngc.Place(e, 1)], SIMULTANEOUS)
        == "SIMULTANEOUS does not accept Place"
    )
    # Controls: without the last action, the computation's error.
    assert refused(w, ordered) == PERIOD_TWO
    assert refused(w, simultaneous, SIMULTANEOUS) == PERIOD_TWO


# -- ORDERED -----------------------------------------------------------------------------------


def test_earlier_actions_are_served_first():
    # A reservation protects nothing: on one link of 10 the first action takes the link.
    for tunnel_first in (False, True):
        w = World(g1())
        e = w.add(0, 1, ecmp())
        t = w.add(0, 1, t1())
        re, rt = ngc.Replace(e, 10), ngc.Replace(t, 10, design_demand=10)
        r = run(w, [rt, re] if tunnel_first else [re, rt])
        er, tr = (r[1], r[0]) if tunnel_first else (r[0], r[1])
        first, second = (tr, er) if tunnel_first else (er, tr)
        assert first.call_delivered == 10 and first.call_edge_drops == []
        assert second.call_delivered == 0 and second.call_edge_drops == [(0, 10)]
        assert tr.by_tunnel[0].status == UP
        assert w.ledger.carried().tolist() == [10]
        assert w.ledger.reserved().tolist() == [10]
        assert w.ledger.residual().tolist() == [0]
        assert w.ledger.unreserved().tolist() == [0]


def test_flow_ids_in_list_order():
    w = World(g1())
    a = w.add(0, 1, ecmp())
    b = w.add(0, 1, ecmp())
    run(w, [ngc.Replace(b, 1), ngc.Replace(a, 1)])
    assert b.records()[0].flow_id == (w.ledger.uid(), 0)
    assert a.records()[0].flow_id == (w.ledger.uid(), 1)


def test_ordered_outputs_after_each_action():
    w = World(g1())
    a = w.add(0, 1, ecmp())
    b = w.add(0, 1, ecmp())
    o = ngc.OutputOptions(edge_flows=True, residual=True, saturated_boundary_edges=True)
    ra, rb = run(w, [ngc.Replace(a, 7), ngc.Replace(b, 7)], outputs=o)
    assert ra.call_delivered == 7 and ra.residual.tolist() == [3]
    assert ra.saturated_boundary_edges == [] and ra.edge_flows.tolist() == [7]
    assert rb.call_delivered == 3 and rb.call_edge_drops == [(0, 4)]
    assert rb.residual.tolist() == [0] and rb.saturated_boundary_edges == [0]
    assert rb.edge_flows.tolist() == [3]


def test_ordered_holds_after_each_action():
    for reversed_order in (False, True):
        w = World(g1())
        t = w.add(0, 1, t1())
        u = w.add(0, 1, t1())
        o = ngc.OutputOptions(edge_reservations=True)
        rt = ngc.Replace(t, 10, design_demand=10)
        ru = ngc.Replace(u, 10, design_demand=10)
        first, second = run(w, [ru, rt] if reversed_order else [rt, ru], outputs=o)
        assert first.by_tunnel[0].status == UP and first.call_delivered == 10
        assert first.edge_reservations.tolist() == [10]
        assert second.by_tunnel[0].status == ngc.TunnelStatus.INSUFFICIENT_BANDWIDTH
        assert second.call_source_drops == 10 and second.call_delivered == 0
        assert second.edge_reservations.tolist() == [0]
        assert w.ledger.reserved().tolist() == [10]


def test_normal_outcomes_commit():
    w = World(g1())
    e = w.add(0, 1, ecmp())
    f = w.add(0, 1, ecmp(LL))
    e2 = w.add(0, 1, ecmp())
    re, rf = run(w, [ngc.Replace(e, 10), ngc.Replace(f, 5)])
    assert re.call_delivered == 10
    assert rf.call_delivered == 0 and rf.call_undelivered == 5
    assert rf.termination == ngc.Termination.EXHAUSTED
    assert [r.delivered for r in f.records()] == [0]
    # A zero place on a new policy: zeros, no program, no record, flags false; nothing published.
    before = w.state()
    (z,) = run(w, [ngc.Place(e2, 0)], outputs=ngc.OutputOptions(residual=True))
    assert z.call_delivered == 0 and z.program is None
    assert not (z.routes_changed or z.weights_changed or z.reservations_changed)
    assert z.residual.tolist() == [0]
    assert w.state() == before


def test_error_leaves_everything_unchanged():
    rows = [
        (0, 1, 1, 1),
        (1, 2, 1000, 1),
        (1, 5, 1000, 1),
        (2, 3, 2, 1),
        (3, 0, 1000, 1),
        (3, 5, 1000, 1),
        (4, 0, 1000, 1),
        (4, 2, 1000, 1),
        (6, 7, 10, 1),
    ]
    w = World(graph(8, rows))
    a = w.add(6, 7, ecmp())
    pc = w.pin(pc_program(w.g))
    a.place(2)
    old = a.records()[0].flow_id
    assert refused(w, [ngc.Replace(a, 5), ngc.Replace(pc, 3)]) == PERIOD_TWO
    assert a.delivered() == 2 and a.records()[0].flow_id == old


def test_infinite_ordered_actions():
    w = World(g1())
    f = w.add(0, 1, ecmp(LL))
    e = w.add(0, 1, ecmp())
    rf, re = run(w, [ngc.Replace(f, INF), ngc.Replace(e, 1)])
    assert rf.call_delivered == 10 and rf.termination == ngc.Termination.EXHAUSTED
    assert rf.requested == INF and rf.offered is None
    assert rf.unadmitted is None and rf.call_undelivered == INF
    assert re.call_delivered == 0 and re.call_edge_drops == [(0, 1)]
    assert re.termination == ngc.Termination.COMPLETE

    def par():
        return graph(2, [(0, 1, 4, 1), (0, 1, 10, 1)])

    w = World(par())
    e = w.add(0, 1, ecmp())
    (r,) = run(w, [ngc.Place(e, INF)])
    assert r.requested == INF and r.offered == 14 and r.call_delivered == 14
    assert r.call_edge_drops == [] and r.termination == ngc.Termination.EXHAUSTED
    assert (
        refused(w, [ngc.Place(e, 0)])
        == "finite placement after infinite placement requires a new epoch"
    )
    (r,) = run(w, [ngc.Replace(e, 12)])
    assert r.call_delivered == 10 and r.call_edge_drops == [(0, 2)]
    before = w.state()
    run(w, [ngc.Place(e, 0)])
    assert w.state() == before


def test_counters_cross_two_to_the_64():
    w = World(g1())
    e = w.add(0, 1, ecmp())
    for _ in range(4):
        e.place(W62)
    old = e.records()[0].flow_id
    (r,) = run(w, [ngc.Place(e, W62)])
    assert r.requested == r.offered == r.admitted == 5 * W62 == 2**64 + 2**62 - 5
    assert r.delivered == 10 and r.edge_drops == [(0, 5 * W62 - 10)]
    assert r.call_delivered == 0 and r.call_undelivered == W62
    assert r.call_edge_drops == [(0, W62)]
    assert e.records()[0].flow_id == old


# -- SIMULTANEOUS ------------------------------------------------------------------------------


def test_ordered_against_simultaneous():
    # docs/semantics.md's call: two offers of 10 on a link of 10.
    for mode, want, drops in (
        (SIMULTANEOUS, [5, 5], [[(0, 5)], [(0, 5)]]),
        (ORDERED, [10, 0], [[], [(0, 10)]]),
    ):
        w = World(g1())
        a = w.add(0, 1, ecmp())
        b = w.add(0, 1, ecmp())
        r = run(w, [ngc.Replace(a, 10), ngc.Replace(b, 10)], mode)
        assert delivered(r) == want
        assert [x.call_edge_drops for x in r] == drops


def released_world():
    """An outside owner O (4), then A (1) and B (5) on one link of 10."""
    w = World(g1())
    o, a, b = (w.add(0, 1, ecmp()) for _ in range(3))
    o.place(4)
    a.place(1)
    b.place(5)
    return w, o, a, b


def test_simultaneous_view_releases_together():
    w, o, a, b = released_world()
    outside = policy_state(o)
    r = run(w, [ngc.Replace(a, 10), ngc.Replace(b, 10)], SIMULTANEOUS)
    assert delivered(r) == [3, 3] and r[0].call_edge_drops == [(0, 7)]
    assert [x.delivery_delta for x in r] == [2, -2]
    assert w.ledger.carried().tolist() == [10]
    assert policy_state(o) == outside
    # ORDERED (the contrast): A on 0 + 1, B on 0 + 5.
    w3, _, a3, b3 = released_world()
    assert delivered(run(w3, [ngc.Replace(a3, 10), ngc.Replace(b3, 10)])) == [1, 5]


def test_permutation_changes_nothing():
    w, _, a, b = released_world()
    r = run(w, [ngc.Replace(a, 10), ngc.Replace(b, 10)], SIMULTANEOUS)
    w2, _, a2, b2 = released_world()
    r2 = run(w2, [ngc.Replace(b2, 10), ngc.Replace(a2, 10)], SIMULTANEOUS)
    assert plain(r2) == plain([r[1], r[0]])
    assert w2.ledger.carried().tolist() == w.ledger.carried().tolist()
    assert (a2.delivered(), b2.delivered()) == (a.delivered(), b.delivered())
    # The new records' FlowIds follow list order.
    assert b2.records()[0].flow_id[1] + 1 == a2.records()[0].flow_id[1]


def test_repeat_gives_the_same_state():
    w, _, a, b = released_world()
    run(w, [ngc.Replace(a, 10), ngc.Replace(b, 10)], SIMULTANEOUS)
    old = a.records()[0].flow_id
    again = run(w, [ngc.Replace(a, 10), ngc.Replace(b, 10)], SIMULTANEOUS)
    assert delivered(again) == [3, 3] and again[0].delivery_delta == 0
    assert a.records()[0].flow_id != old
    # Removing A later frees only its traffic; B does not rebalance.
    a.remove()
    assert w.ledger.carried().tolist() == [7] and b.delivered() == 3


def test_holds_and_orphans_stay_outside():
    w = World(g1())
    t = w.add(0, 1, t1())
    e = w.add(0, 1, ecmp())
    t.replace(0, design_demand=10)
    assert run(w, [ngc.Replace(e, 10)], SIMULTANEOUS)[0].call_delivered == 10
    w = World(g1())
    w.add(0, 1, ecmp()).place(4)
    w.p.pop()  # dropped without remove(): its record is an orphan
    gc.collect()
    a = w.add(0, 1, ecmp())
    b = w.add(0, 1, ecmp())
    r = run(w, [ngc.Replace(a, 10), ngc.Replace(b, 10)], SIMULTANEOUS)
    assert delivered(r) == [3, 3] and w.ledger.carried().tolist() == [10]


def test_simultaneous_preparation_reads_the_view():
    for mode in (SIMULTANEOUS, ORDERED):
        w = World(graph(2, [(0, 1, 10, 1), (0, 1, 10, 1)]))
        snap = shortest(
            allocation=ngc.RouteAllocation.WEIGHTED,
            weight_basis=ngc.WeightBasis.RESIDUAL_SNAPSHOT,
        )
        a = w.add(0, 1, config(snap, admission=LY))
        a.replace(2)
        w.pin(tunnel(w.g, 0, 1, [0])).place(8)
        b = w.add(0, 1, ecmp())
        ra, rb = run(w, [ngc.Replace(a, 12), ngc.Replace(b, 4)], mode)
        if mode == SIMULTANEOUS:
            assert ra.weights_changed and not ra.routes_changed
            assert ra.call_delivered == 10 and ra.call_edge_drops == [(0, 1), (1, 1)]
            assert rb.call_delivered == 2 and rb.call_edge_drops == [(0, 1), (1, 1)]
            assert w.ledger.carried().tolist() == [10, 10]
        else:
            assert ra.call_delivered == 12
            assert rb.call_delivered == 0 and rb.call_edge_drops == [(0, 2), (1, 2)]


def test_simultaneous_preparation_reads_the_list_view():
    # The other participant's old load is released before the preparation: B (a tunnel on edge 0)
    # carries 8, the list replaces A (12) and B (8). On the list's view A's snapshot is (10, 10): A
    # splits (6, 6) and edge 0 shares 10 over A 6 and B 8, (4, 6) (A created first); on A's own
    # view it would be (2, 10) with no drop. A single path by HIGHER_AVAILABLE_CAPACITY ties on the
    # list's view and takes edge 0 (edge 1 on A's own view): edge 0 shares 10 over 12 and 8, (6, 4).
    snap = shortest(
        allocation=ngc.RouteAllocation.WEIGHTED,
        weight_basis=ngc.WeightBasis.RESIDUAL_SNAPSHOT,
    )
    single = shortest(
        tie_break=ngc.RouteTieBreak.HIGHER_AVAILABLE_CAPACITY,
        multipath=False,
        parallel_edges=False,
    )
    for model, a_flows, want, carried in (
        (snap, [(0, 4), (1, 6)], [10, 6], [10, 6]),
        (single, [(0, 6)], [6, 4], [10, 0]),
    ):
        w = World(graph(2, [(0, 1, 10, 1), (0, 1, 10, 1)]))
        a = w.add(0, 1, config(model, admission=LY))
        b = w.pin(tunnel(w.g, 0, 1, [0]))
        b.replace(8)
        r = run(w, [ngc.Replace(a, 12), ngc.Replace(b, 8)], SIMULTANEOUS)
        assert [rec.edge_flows for rec in a.records()] == [a_flows]
        assert delivered(r) == want and w.ledger.carried().tolist() == carried


def test_ties_follow_creation_order():
    for first_created_listed_second in (True, False):
        w = World(graph(2, [(0, 1, 1, 1)]))
        prog = tunnel(w.g, 0, 1, [0])
        p0, p1 = w.pin(prog), w.pin(prog)
        if first_created_listed_second:
            r = run(w, [ngc.Replace(p1, 1), ngc.Replace(p0, 1)], SIMULTANEOUS)
            assert delivered(r) == [0, 1] and r[0].call_edge_drops == [(0, 1)]
        else:
            r = run(w, [ngc.Replace(p0, 1), ngc.Replace(p1, 1)], SIMULTANEOUS)
            assert delivered(r) == [1, 0]


def test_pooled_demand_shares_a_link():
    # Every field of both entries.
    w = World(
        graph(
            5, [(0, 2, 10, 1), (1, 2, 10, 1), (2, 3, 1, 1), (2, 3, 1, 1), (4, 2, 10, 1)]
        )
    )
    a = w.add([0, 1], 3, ecmp())
    b = w.pin(tunnel(w.g, 4, 3, [4, 2]))
    ra, rb = run(w, [ngc.Replace(a, 2), ngc.Replace(b, 2)], SIMULTANEOUS, ALL)
    assert (ra.call_delivered, ra.call_undelivered, ra.delivery_delta) == (1, 1, 1)
    assert ra.call_edge_drops == [(2, 1)] and ra.call_source_drops == 0
    assert ra.termination == ngc.Termination.COMPLETE
    assert ra.routes_changed and ra.weights_changed and not ra.reservations_changed
    assert ra.program is not None and len(ra.program.tunnels) == 0 and ra.plan is None
    assert (ra.requested, ra.offered, ra.admitted, ra.delivered) == (2, 2, 2, 1)
    assert ra.unadmitted == 0 and ra.source_drops == 0 and ra.edge_drops == [(2, 1)]
    assert [
        (s.source, s.offered, s.delivered, s.total_edge_drops) for s in ra.by_source
    ] == [
        (0, 1, 0, 1),
        (1, 1, 1, 0),
    ]
    assert len(ra.by_tunnel) == 0 and ra.total_cost is None
    assert ra.marginal_cost_curve is None and ra.min_cut is None
    assert ra.edge_flows.tolist() == [1, 1, 0, 1, 0]
    assert ra.residual.tolist() == [9, 9, 0, 0, 8]
    assert ra.saturated_boundary_edges == [2, 3]
    assert ra.edge_reservations.tolist() == [0] * 5
    recs = a.records()
    assert [r.edge_flows for r in recs] == [[(0, 1)], [(1, 1), (3, 1)]]
    assert [len(r.edge_drops) for r in recs] == [1, 0]
    assert (rb.call_delivered, rb.call_undelivered, rb.delivery_delta) == (1, 1, 1)
    assert rb.call_edge_drops == [(2, 1)] and rb.termination == ngc.Termination.COMPLETE
    assert rb.routes_changed and rb.weights_changed
    assert [
        (x.head, x.status, x.assigned, x.delivered, x.total_edge_drops)
        for x in rb.by_tunnel
    ] == [(4, UP, 2, 1, 1)]
    assert rb.edge_flows.tolist() == [0, 0, 1, 0, 2]
    assert rb.residual.tolist() == [9, 9, 0, 0, 8]
    assert w.ledger.carried().tolist() == [1, 1, 1, 1, 2]


def test_results_in_input_order():
    w = World(graph(4, [(0, 1, 10, 1), (2, 3, 10, 1)]))
    a = w.pin(tunnel(w.g, 0, 1, [0]))
    b = w.pin(tunnel(w.g, 2, 3, [1]))
    for mode in (ORDERED, SIMULTANEOUS):
        r = run(w, [ngc.Replace(b, 7), ngc.Replace(a, 3)], mode)
        assert delivered(r) == [7, 3] and r[0].program.tunnels[0].head == 2


def test_demand_level_above_tunnels():
    w = World(graph(2, [(0, 1, 1, 1)]))
    a = w.pin(tunnels(w.g, [0], [1], [(0, [0], 1), (0, [0], 1)]))
    b = w.pin(tunnel(w.g, 0, 1, [0]))
    ra, rb = run(w, [ngc.Replace(a, 2), ngc.Replace(b, 2)], SIMULTANEOUS)
    assert ra.call_delivered == 1 and [x.delivered for x in ra.by_tunnel] == [1, 0]
    assert rb.call_delivered == 0
    # Created the other way round, the tie goes to the other demand.
    w2 = World(graph(2, [(0, 1, 1, 1)]))
    b2 = w2.pin(tunnel(w2.g, 0, 1, [0]))
    a2 = w2.pin(tunnels(w2.g, [0], [1], [(0, [0], 1), (0, [0], 1)]))
    r2 = run(w2, [ngc.Replace(a2, 2), ngc.Replace(b2, 2)], SIMULTANEOUS)
    assert delivered(r2) == [0, 1]


def test_demand_level_above_sources():
    w = World(graph(6, [(1, 3, 100, 1), (2, 3, 100, 1), (3, 4, 1, 1), (5, 3, 100, 1)]))
    a = w.pin(tunnels(w.g, [1, 2], [4], [(1, [0, 2], 1), (2, [1, 2], 1)]))
    b = w.pin(tunnel(w.g, 5, 4, [3, 2]))
    ra, rb = run(w, [ngc.Replace(a, 2), ngc.Replace(b, 2)], SIMULTANEOUS)
    assert ra.call_delivered == 1 and ra.by_tunnel[0].delivered == 1
    assert rb.call_delivered == 0


def test_black_holes_of_one_participant():
    w = World(g1())
    a = w.pin(tunnel(w.g, 0, 1, [0]), black_holes=no_edges(w.g, [0]))
    b = w.add(0, 1, ecmp())
    o = ngc.OutputOptions(residual=True, saturated_boundary_edges=True)
    ra, rb = run(w, [ngc.Replace(a, 10), ngc.Replace(b, 10)], SIMULTANEOUS, o)
    assert ra.call_delivered == 0 and ra.call_edge_drops == [(0, 10)]
    assert ra.by_tunnel[0].status == BLACKHOLED
    assert rb.call_delivered == 10 and rb.call_edge_drops == []
    assert ra.residual.tolist() == [0] == rb.residual.tolist()
    assert ra.saturated_boundary_edges == [0]
    # A Forward's own black holes.
    w = World(graph(2, [(0, 1, 10, 1), (0, 1, 10, 1)]))
    a = w.add(0, 1, ecmp())
    b = w.add(0, 1, ecmp())
    a.replace(2)
    ra, rb = run(
        w,
        [ngc.Forward(a, 10, black_holes=no_edges(w.g, [1])), ngc.Replace(b, 10)],
        SIMULTANEOUS,
    )
    assert ra.call_delivered == 5 and ra.call_edge_drops == [(1, 5)]
    assert rb.call_delivered == 10 and w.ledger.carried().tolist() == [10, 5]


def test_black_hole_downstream_of_a_shared_link():
    w = World(graph(3, [(0, 1, 10, 1), (1, 2, 10, 1)]))
    a = w.pin(tunnel(w.g, 0, 2, [0, 1]), black_holes=no_edges(w.g, [1]))
    b = w.pin(tunnel(w.g, 0, 2, [0, 1]))
    o = ngc.OutputOptions(residual=True, saturated_boundary_edges=True)
    ra, rb = run(w, [ngc.Replace(a, 10), ngc.Replace(b, 10)], SIMULTANEOUS, o)
    assert ra.call_edge_drops == [(0, 5), (1, 5)]
    assert ra.by_tunnel[0].status == BLACKHOLED
    assert rb.call_delivered == 5 and rb.call_edge_drops == [(0, 5)]
    assert ra.residual.tolist() == [0, 5]
    assert ra.saturated_boundary_edges == [0] == rb.saturated_boundary_edges


def test_forward_keeps_holds():
    w = World(g1())
    t = w.add(0, 1, t1())
    e = w.add(0, 1, ecmp())
    t.replace(5, design_demand=5)
    rt, re = run(w, [ngc.Forward(t, 10), ngc.Replace(e, 10)], SIMULTANEOUS)
    assert rt.call_delivered == 5 and rt.call_edge_drops == [(0, 5)]
    assert not rt.routes_changed and not rt.reservations_changed
    assert rt.plan.total_reserved_bandwidth == 5
    assert [r.edge_reservations for r in t.records()] == [[(0, 5)]]
    assert re.call_delivered == 5
    assert w.ledger.carried().tolist() == [10] and w.ledger.reserved().tolist() == [5]


def test_empty_roster_participants():
    w = World(g1())
    mask = no_edges(w.g, [0])
    m = w.add(0, 1, ecmp(), mask)
    e = w.add(0, 1, ecmp())
    p = w.pin(tunnel(w.g, 0, 1, [0]), mask=mask)
    rm, re, rp = run(
        w, [ngc.Replace(m, 5), ngc.Replace(e, 10), ngc.Replace(p, 5)], SIMULTANEOUS
    )
    assert rm.termination == ngc.Termination.NO_ROUTE and rm.call_source_drops == 5
    assert rm.offered == 5 and rm.source_drops == 5 and m.records() == []
    assert re.call_delivered == 10
    assert rp.termination == ngc.Termination.NO_ROUTE and rp.call_source_drops == 5
    assert [r.delivered for r in p.records()] == [0]


def test_everything_masked():
    for mode in (ORDERED, SIMULTANEOUS):
        w = World(g1())
        mask = no_edges(w.g, [0])
        m1 = w.add(0, 1, ecmp(), mask)
        m2 = w.add(0, 1, ecmp(), mask)
        p = w.pin(tunnel(w.g, 0, 1, [0]), mask=mask)
        r = run(w, [ngc.Replace(x, 5) for x in (m1, m2, p)], mode)
        assert all(x.termination == ngc.Termination.NO_ROUTE for x in r)
        assert all(x.call_source_drops == 5 for x in r)
        assert m1.records() == [] and m2.records() == [] and len(p.records()) == 1
        assert w.ledger.carried().tolist() == [0]


def test_zero_offers():
    w = World(g1())
    a = w.add(0, 1, ecmp())
    b = w.add(0, 1, ecmp())
    for x in run(w, [ngc.Replace(a, 0), ngc.Replace(b, 0)], SIMULTANEOUS):
        assert x.termination == ngc.Termination.COMPLETE and x.call_delivered == 0
        assert x.program is not None and x.routes_changed and x.weights_changed
    assert len(a.records()) == len(b.records()) == 1
    assert w.ledger.carried().tolist() == [0]


def test_iteration_errors_change_nothing():
    # D-C3 as three Pinned tunnels, each link on its own nodes with wide connectors.
    g, i = mapped(
        8,
        [
            (0, 1, 1, 1),
            (2, 3, 1, 1),
            (4, 5, 2, 1),
            (1, 2, 1000, 1),
            (3, 4, 1000, 1),
            (5, 0, 1000, 1),
            (6, 0, 1000, 1),
            (6, 2, 1000, 1),
            (6, 4, 1000, 1),
            (3, 7, 1000, 1),
            (5, 7, 1000, 1),
            (1, 7, 1000, 1),
        ],
    )
    w = World(g)
    a = w.pin(tunnel(g, 6, 7, [i[k] for k in (6, 0, 3, 1, 9)]))
    b = w.pin(tunnel(g, 6, 7, [i[k] for k in (7, 1, 4, 2, 10)]))
    c = w.pin(tunnel(g, 6, 7, [i[k] for k in (8, 2, 5, 0, 11)]))
    for x in (a, b, c):
        x.replace(1)
    actions = [ngc.Replace(a, 1), ngc.Replace(b, 1), ngc.Replace(c, 2)]
    assert refused(w, actions, SIMULTANEOUS) == PERIOD_TWO


def test_dark_arrivals_add_no_barrier():
    # The round-limit example's graph plus 16 -> 6; its three demands (black holes {5}, used by
    # head 12 only) and a demand D from 16 whose black holes {16, 11} drop its offer on its first
    # link: the round-limit error, nothing published (a dark copy of link 11 as a barrier would
    # raise the limit).
    g, i = mapped(
        17,
        [
            (0, 1, 12, 1),
            (1, 6, 1000, 1),
            (2, 3, 9, 1),
            (3, 0, 1000, 1),
            (3, 10, 1000, 1),
            (4, 5, 7, 1),
            (5, 6, 1000, 1),
            (6, 7, 13, 1),
            (7, 2, 1000, 1),
            (7, 10, 1000, 1),
            (8, 9, 1, 1),
            (10, 11, 10, 1),
            (11, 2, 1000, 1),
            (12, 4, 1000, 1),
            (13, 6, 1000, 1),
            (14, 0, 1000, 1),
            (16, 6, 1000, 1),
        ],
    )
    w = World(g)

    def ids(rows):
        return [i[k] for k in rows]

    a = w.pin(
        tunnel(g, 12, 11, ids([13, 5, 6, 7, 8, 2, 4, 11])),
        black_holes=no_edges(g, ids([5])),
    )
    b = w.pin(tunnel(g, 13, 1, ids([14, 7, 9, 11, 12, 2, 3, 0])))
    c = w.pin(tunnel(g, 14, 7, ids([15, 0, 1, 7])))
    three = [ngc.Replace(a, 13), ngc.Replace(b, 11), ngc.Replace(c, 10)]
    assert refused(w, three, SIMULTANEOUS) == ROUND_LIMIT
    d = w.pin(
        tunnel(g, 16, 1, ids([16, 7, 9, 11, 12, 2, 3, 0])),
        black_holes=no_edges(g, ids([16, 11])),
    )
    assert refused(w, [*three, ngc.Replace(d, 1)], SIMULTANEOUS) == ROUND_LIMIT
    (r,) = run(w, [ngc.Replace(d, 1)], SIMULTANEOUS)
    assert r.call_edge_drops == [(i[16], 1)] and r.by_tunnel[0].status == BLACKHOLED


def test_zero_shares_add_no_arc():
    # An entry whose share is 0 can carry nothing: it makes no link wait and joins no group. Two
    # cycles bridged by two tunnels of offer 0 fail with the period-2 error with or without
    # the bridges; so does a link X before a period-2 pair with or without Z (offer 0) over L0 then
    # X, as a Pinned tunnel or a ShortestPaths demand.
    cycle = [
        (0, 1, 1, 1),
        (1, 2, 1000, 1),
        (1, 5, 1000, 1),
        (2, 3, 2, 1),
        (3, 0, 1000, 1),
        (3, 5, 1000, 1),
        (4, 0, 1000, 1),
        (4, 2, 1000, 1),
    ]
    other = [
        (0, 1, 1, 1),
        (1, 4, 1000, 1),
        (2, 3, 2, 1),
        (3, 0, 1000, 1),
        (3, 7, 1000, 1),
        (4, 5, 1, 1),
        (5, 2, 1000, 1),
        (5, 7, 1000, 1),
        (6, 0, 1000, 1),
        (6, 2, 1000, 1),
    ]
    rows = cycle + [(a + 10, b + 10, c, d) for a, b, c, d in other]
    rows += [(1, 10, 1000, 1), (11, 0, 1000, 1)]
    bridged = [
        (4, 5, [6, 0, 1, 3, 5], 1),
        (4, 5, [7, 3, 4, 0, 2], 2),
        (16, 17, [16, 8, 9, 13, 14, 10, 12], 1),
        (16, 17, [17, 10, 11, 8, 9, 13, 15], 2),
        (4, 17, [6, 0, 18, 8, 9, 13, 15], 0),
        (16, 5, [16, 8, 19, 0, 2], 0),
    ]
    for count in (4, 6):
        g, i = mapped(18, rows)
        w = World(g)
        actions = []
        for s, d, path, offer in bridged[:count]:
            p = w.pin(tunnel(g, s, d, [i[k] for k in path]))
            p.replace(1)
            actions.append(ngc.Replace(p, offer))
        assert refused(w, actions, SIMULTANEOUS) == PERIOD_TWO, count
    # sA 0, sB 1, sZ 2, n1..n6 3..8, dA 9, dB 10, dZ 11; rows X, L0, L1 first.
    rows = [
        (3, 4, 1, 1),
        (5, 6, 1, 1),
        (7, 8, 2, 1),
        (0, 3, 100, 1),
        (4, 5, 100, 1),
        (6, 7, 100, 1),
        (8, 9, 100, 1),
        (1, 7, 100, 1),
        (8, 5, 100, 1),
        (6, 10, 100, 1),
        (2, 5, 100, 1),
        (6, 3, 100, 1),
        (4, 11, 100, 1),
    ]
    for z in ("none", "tunnel", "shortest paths"):
        g, i = mapped(12, rows)
        w = World(g)
        a = w.pin(tunnel(g, 0, 9, [i[k] for k in (3, 0, 4, 1, 5, 2, 6)]))
        b = w.pin(tunnel(g, 1, 10, [i[k] for k in (7, 2, 8, 1, 9)]))
        actions = [ngc.Replace(a, 2), ngc.Replace(b, 2)]
        if z == "tunnel":
            zt = w.pin(tunnel(g, 2, 11, [i[k] for k in (10, 1, 11, 0, 12)]))
            actions.append(ngc.Replace(zt, 0))
        if z == "shortest paths":
            actions.append(ngc.Replace(w.add(2, 11, ecmp()), 0))
        assert refused(w, actions, SIMULTANEOUS) == PERIOD_TWO, z
    # A member that can carry nothing in a program that is not pooled: Z's program from sZ has a
    # supplied DAG route (sZ's exit to dZ0 below its L0 branch, L0 then X) and a tunnel over that
    # exit; offer 2 gives each entry 1 and the DAG route's unit takes the exit, so the L0 branch never
    # carries anything. sA 0, sB 1, sZ 2, dZ0 3, n1..n6 4..9, dA 10, dB 11, dZ 12; rows X, L0, L1 first.
    rows = [
        (4, 5, 1, 1),
        (6, 7, 1, 1),
        (8, 9, 2, 1),
        (0, 4, 100, 1),
        (5, 6, 100, 1),
        (7, 8, 100, 1),
        (9, 10, 100, 1),
        (1, 8, 100, 1),
        (9, 6, 100, 1),
        (7, 11, 100, 1),
        (2, 6, 100, 1),
        (7, 4, 100, 1),
        (5, 12, 100, 1),
        (2, 3, 100, 1),
    ]
    g, i = mapped(13, rows)
    assert i[13] < i[10]
    w = World(g)
    a = w.pin(tunnel(g, 0, 10, [i[k] for k in (3, 0, 4, 1, 5, 2, 6)]))
    b = w.pin(tunnel(g, 1, 11, [i[k] for k in (7, 2, 8, 1, 9)]))
    routes = [
        ngc.RouteDefinition(2, dag_of(g, [i[k] for k in (10, 1, 11, 0, 12, 13)])),
        ngc.RouteDefinition(2, ngc.Path(g, [[i[13]]])),
    ]
    z = w.pin(ngc.ForwardingProgram.from_routes(g, [2], [3, 12], routes))
    actions = [ngc.Replace(a, 2), ngc.Replace(b, 2), ngc.Replace(z, 2)]
    assert refused(w, actions, SIMULTANEOUS) == PERIOD_TWO


def test_members_that_carry_nothing_change_nothing():
    # D-cycle-subdivision-after (A over links 1, 0, 2, 3 with offer 2; B over 2, 3, 0 with offer 1;
    # C over 2, 3, 0, 1 with offer 4; links of 2, 4, 2, 2; links 2 and 3 joined directly at node 10)
    # delivers 0, 0, 1. B becomes a DAG program with a member that carries nothing at node 10: a
    # second pooled source (node 12, its own link into node 10) with share 0, or whose link is B's
    # own black hole; a branch from node 10 to a destination of its own (node 13) that round robin
    # leaves empty; a merge into node 10 from B's source over node 12 that round robin leaves empty
    # (pooled, and as a DAG route beside a tunnel of weight 0). Links 2 and 3 stay a pair: the same
    # deliveries and link flows as without the member. sA 0, sB 1, sC 2, d 3; link e from 4 + 2e.
    rows = [
        (4, 5, 2, 1),
        (6, 7, 4, 1),
        (8, 10, 2, 1),
        (10, 11, 2, 1),
        (0, 6, 100, 1),
        (7, 4, 100, 1),
        (5, 8, 100, 1),
        (11, 3, 100, 1),
        (1, 8, 100, 1),
        (11, 4, 100, 1),
        (5, 3, 100, 1),
        (2, 8, 100, 1),
        (5, 6, 100, 1),
        (7, 3, 100, 1),
        (12, 10, 9, 1),
        (10, 13, 100, 1),
        (1, 12, 100, 1),
    ]

    def call(case, pooled=True):
        g, i = mapped(14, rows)
        assert i[3] < i[15] and i[8] < i[16]
        w = World(g)
        a = w.pin(tunnel(g, 0, 3, [i[k] for k in (4, 1, 5, 0, 6, 2, 3, 7)]))
        path = [i[k] for k in (8, 2, 3, 9, 0, 10)]
        members, sources, dests = list(path), [1], [3]
        if case in ("zero share", "black hole"):
            members, sources = [*members, i[14]], [1, 12]
        if case == "branch":
            members, dests = [*members, i[15]], [3, 13]
        if case == "merge":
            members = [*members, i[16], i[14]]
        if pooled:
            program = ngc.ForwardingProgram.from_dag(
                g, sources, dests, dag_of(g, members)
            )
        else:
            routes = [
                ngc.RouteDefinition(1, dag_of(g, members)),
                ngc.RouteDefinition(1, ngc.Path(g, [[e] for e in path]), weight=0),
            ]
            program = ngc.ForwardingProgram.from_routes(g, sources, dests, routes)
        black = no_edges(g, [i[14]]) if case == "black hole" else None
        b = w.pin(program, LY, black)
        c = w.pin(tunnel(g, 2, 3, [i[k] for k in (11, 2, 3, 9, 0, 12, 1, 13)]))
        offer = 2 if case == "black hole" else 1
        actions = [ngc.Replace(a, 2), ngc.Replace(b, offer), ngc.Replace(c, 4)]
        try:
            results = run(w, actions, SIMULTANEOUS)
        except ValueError as e:
            return str(e)
        carried = w.ledger.carried().tolist()
        return [r.call_delivered for r in results], [carried[i[e]] for e in range(4)]

    for case, pooled in (
        ("zero share", True),
        ("black hole", True),
        ("branch", True),
        ("branch", False),
        ("merge", True),
        ("merge", False),
    ):
        without = call("none", pooled)
        assert without[0] == [0, 0, 1], (case, pooled)
        assert call(case, pooled) == without, (case, pooled)


def test_dead_tails_take_no_part():
    # Traffic behind a demand's own black hole can carry nothing: D-dead-tail (A over links 2, 0
    # with offer 1; B over 1, 0, 2 with offer 2; links of 2, 2, 1, 9) gives the period-2 error
    # without Z, and with Z (its own black hole on link 3) over links 3, 0, 1 with offer 1 as a
    # tunnel or as a supplied DAG route (a pooled context), or with offer 2 as that DAG route
    # beside a tunnel on its own exit (a split route of a program that is not pooled). sA 0, sB 1,
    # sZ 2, d 3; link e from 4 + 2e to 5 + 2e (rows 0 to 3), joined by links of 100; Z's exit is
    # row 15.
    rows = [
        (4, 5, 2, 1),
        (6, 7, 2, 1),
        (8, 9, 1, 1),
        (10, 11, 9, 1),
        (0, 8, 100, 1),
        (9, 4, 100, 1),
        (5, 3, 100, 1),
        (1, 6, 100, 1),
        (7, 4, 100, 1),
        (5, 8, 100, 1),
        (9, 3, 100, 1),
        (2, 10, 100, 1),
        (11, 4, 100, 1),
        (5, 6, 100, 1),
        (7, 3, 100, 1),
        (2, 3, 100, 1),
    ]
    for z in ("none", "tunnel", "DAG route", "DAG route and tunnel"):
        g, i = mapped(12, rows)
        w = World(g)
        a = w.pin(tunnel(g, 0, 3, [i[k] for k in (4, 2, 5, 0, 6)]))
        b = w.pin(tunnel(g, 1, 3, [i[k] for k in (7, 1, 8, 0, 9, 2, 10)]))
        actions = [ngc.Replace(a, 1), ngc.Replace(b, 2)]
        zr = [i[k] for k in (11, 3, 12, 0, 13, 1, 14)]
        black = no_edges(g, [i[3]])
        if z == "tunnel":
            actions.append(ngc.Replace(w.pin(tunnel(g, 2, 3, zr), LY, black), 1))
        if z.startswith("DAG route"):
            routes = [ngc.RouteDefinition(2, dag_of(g, zr))]
            if z == "DAG route and tunnel":
                routes.append(ngc.RouteDefinition(2, ngc.Path(g, [[i[15]]])))
            program = ngc.ForwardingProgram.from_routes(g, [2], [3], routes)
            offer = 2 if z == "DAG route and tunnel" else 1
            actions.append(ngc.Replace(w.pin(program, LY, black), offer))
        assert refused(w, actions, SIMULTANEOUS) == PERIOD_TWO, z


def test_maximal_offers_share_a_link():
    # Link sums beyond 2^63 and 2^64.
    for n, cap, want in (
        (2, 10, [5, 5]),
        (3, 10, [4, 3, 3]),
        (5, 10, [2, 2, 2, 2, 2]),
        (5, 11, [3, 2, 2, 2, 2]),
    ):
        w = World(graph(2, [(0, 1, cap, 1)]))
        prog = tunnel(w.g, 0, 1, [0])
        r = run(w, [ngc.Replace(w.pin(prog), W62) for _ in range(n)], SIMULTANEOUS)
        assert delivered(r) == want
        assert [x.call_edge_drops for x in r] == [[(0, W62 - k)] for k in want]


def test_disjoint_maximal_participants():
    each = W62 // 5  # the graph's capacity sum stays below 2^62
    w = World(graph(2, [(0, 1, each, 1)] * 5))
    actions = [ngc.Replace(w.pin(tunnel(w.g, 0, 1, [e])), each) for e in range(5)]
    assert delivered(run(w, actions, SIMULTANEOUS)) == [each] * 5
    assert w.ledger.carried().tolist() == [each] * 5


def test_imported_tunnels_book_no_holds():
    w = World(g1(), [8], [6])
    t = w.add(0, 1, t1())
    t.replace(0, design_demand=6)
    p = w.pin(t.program)
    o = ngc.OutputOptions(edge_reservations=True)
    rt, rp = run(w, [ngc.Forward(t, 8), ngc.Replace(p, 8)], SIMULTANEOUS, o)
    assert delivered([rt, rp]) == [4, 4]
    assert w.ledger.carried().tolist() == [8] and w.ledger.reserved().tolist() == [6]
    assert rp.edge_reservations.tolist() == [0] and rp.plan is None
    assert rt.plan is not None
    # The program of an orphaned T1 policy.
    w = World(g1(), [8], [6])
    w.add(0, 1, t1()).replace(0, design_demand=6)
    prog = w.p.pop().program
    gc.collect()
    p = w.pin(prog)
    assert run(w, [ngc.Replace(p, 8)], SIMULTANEOUS)[0].call_delivered == 8
    assert w.ledger.reserved().tolist() == [6]


def test_zero_forward():
    w = World(g1())
    t = w.add(0, 1, t1())
    e = w.add(0, 1, ecmp())
    t.replace(6, design_demand=6)
    old = t.records()[0].flow_id
    prog = t.program
    rt, re = run(w, [ngc.Forward(t, 0), ngc.Replace(e, 10)], SIMULTANEOUS)
    assert rt.termination == ngc.Termination.COMPLETE and rt.call_delivered == 0
    assert rt.delivery_delta == -6 and rt.by_tunnel[0].assigned == 0
    assert plain(t.program) == plain(prog)
    assert len(t.records()) == 1 and t.records()[0].flow_id != old
    assert t.records()[0].edge_reservations == [(0, 6)]
    assert re.call_delivered == 10
    f = t.forward(8)
    assert f.call_delivered == 0 and f.call_edge_drops == [(0, 8)]
    assert w.ledger.reserved().tolist() == [6]


def test_black_holed_sources_and_nodes():
    w = World(graph(5, [(0, 2, 10, 1), (1, 2, 10, 1), (3, 4, 10, 1)]))
    a = w.add([0, 1], 2, ecmp())
    a.replace(2)
    z = w.pin(tunnel(w.g, 3, 4, [2]))
    ra, rz = run(
        w,
        [ngc.Forward(a, 4, black_holes=no_nodes(w.g, [0])), ngc.Replace(z, 3)],
        SIMULTANEOUS,
    )
    assert ra.call_source_drops == 2 and ra.call_edge_drops == []
    assert ra.call_delivered == 2 and rz.call_delivered == 3
    # A repair mask first: source 1 takes the whole offer.
    q = w.pin(a.program, mask=no_nodes(w.g, [0]))
    rq, _ = run(w, [ngc.Replace(q, 4), ngc.Replace(z, 3)], SIMULTANEOUS)
    assert rq.call_delivered == 4 and rq.call_source_drops == 0
    # An internal black-holed node.
    w2 = World(graph(3, [(0, 1, 10, 1), (1, 2, 10, 1)]))
    n = w2.pin(tunnel(w2.g, 0, 2, [0, 1]), black_holes=no_nodes(w2.g, [1]))
    (rn,) = run(w2, [ngc.Replace(n, 6)], SIMULTANEOUS)
    assert rn.by_tunnel[0].status == BLACKHOLED and rn.call_edge_drops == [(0, 6)]


def test_zero_service_and_zero_weight():
    w = World(graph(4, [(0, 1, 0, 1), (2, 3, 10, 1)]))
    z = w.add(0, 1, ecmp())
    pos = w.pin(tunnel(w.g, 2, 3, [1]))
    rz, rp = run(w, [ngc.Replace(z, 2), ngc.Replace(pos, 4)], SIMULTANEOUS)
    assert rz.termination == ngc.Termination.COMPLETE and rz.call_edge_drops == [(0, 2)]
    assert rp.call_delivered == 4
    w2 = World(graph(4, [(0, 1, 10, 1), (2, 3, 10, 1)]))
    zw = w2.pin(tunnel(w2.g, 0, 1, [0], weight=0))
    pos2 = w2.pin(tunnel(w2.g, 2, 3, [1]))
    rz, rp = run(w2, [ngc.Replace(zw, 2), ngc.Replace(pos2, 4)], SIMULTANEOUS)
    assert rz.termination == ngc.Termination.NO_ACTIVE_WEIGHT
    assert rz.call_source_drops == 2 and len(zw.records()) == 1
    assert rp.call_delivered == 4


def test_pooled_roots_block_serial_pairs():
    def make():
        return graph(6, [(0, 1, 2, 1), (1, 2, 2, 1), (2, 3, 2, 1), (4, 5, 10, 1)])

    scalar = World(make())
    want = plain(scalar.add([0, 1, 2], 3, ecmp()).replace(6))
    for beside in (False, True):
        w = World(make())
        a = w.add([0, 1, 2], 3, ecmp())
        z = w.pin(tunnel(w.g, 4, 5, [3]))
        actions = [ngc.Replace(a, 6)] + ([ngc.Replace(z, 1)] if beside else [])
        r = run(w, actions, SIMULTANEOUS)[0]
        assert plain(r) == want
        assert r.call_delivered == 2
        assert [s.delivered for s in r.by_source] == [1, 0, 1]


def test_pooled_cycle_keeps_source_state():
    g, i = mapped(
        9,
        [
            (0, 1, 1, 1),
            (2, 3, 1, 1),
            (4, 0, 100, 1),
            (1, 2, 100, 1),
            (5, 2, 100, 1),
            (3, 6, 100, 1),
            (7, 2, 100, 1),
            (3, 0, 100, 1),
            (1, 8, 100, 1),
        ],
    )
    w = World(g)
    a = w.add([4, 5], 6, ecmp())
    b = w.pin(tunnel(g, 7, 8, [i[6], i[1], i[7], i[0], i[8]]))
    ra, rb = run(w, [ngc.Replace(a, 2), ngc.Replace(b, 1)], SIMULTANEOUS)
    assert ra.call_delivered == 1 and [s.delivered for s in ra.by_source] == [1, 0]
    assert rb.call_delivered == 0 and rb.call_edge_drops == [(i[1], 1)]


def test_outputs_per_entry():
    w = World(graph(4, [(0, 2, 10, 1), (1, 2, 1, 1), (2, 3, 10, 1)]))
    a = w.pin(tunnel(w.g, 0, 3, [0, 2]))
    b = w.pin(tunnel(w.g, 1, 3, [1, 2]))
    o = ngc.OutputOptions(residual=True, saturated_boundary_edges=True)
    ra, rb = run(w, [ngc.Replace(a, 5), ngc.Replace(b, 1)], SIMULTANEOUS, o)
    assert ra.residual.tolist() == [5, 0, 4] == rb.residual.tolist()
    assert ra.saturated_boundary_edges == [] and rb.saturated_boundary_edges == [1]


def test_captured_weights_in_a_list():
    g = graph(3, [(0, 2, 3, 1), (1, 2, 1, 1)])
    prog = ngc.FlowPolicy(ngc.FlowLedger(g), [0, 1], 2, controller()).replace(4).program
    w = World(g)
    a = w.pin(prog)
    b = w.pin(tunnel(g, 0, 2, [0]))
    r = run(w, [ngc.Replace(a, 4), ngc.Replace(b, 1)], SIMULTANEOUS)
    assert delivered(r) == [4, 0] and r[1].call_edge_drops == [(0, 1)]


# -- all or nothing ----------------------------------------------------------------------------


def list_world(mode):
    """Fresh state and a list of three actions that `mode` accepts (two pooled demands and a
    Pinned one with black holes; each policy already holds records)."""
    g = graph(5, [(0, 1, 5, 1), (0, 2, 3, 2), (1, 3, 5, 1), (2, 3, 3, 2), (3, 4, 6, 1)])
    w = World(g)
    a = w.add([0, 1], 4, ecmp())
    b = w.pin(tunnel(g, 0, 4, [1, 3, 4]), black_holes=no_edges(g, [3]))
    c = w.add(0, 4, ecmp(LL) if mode == ORDERED else ecmp())
    for p in (a, b, c):
        p.replace(1)
    return w, [ngc.Replace(a, 6), ngc.Replace(b, 4), ngc.Replace(c, 3)]


@needs_hooks
def test_lists_are_all_or_nothing():
    # Every Python allocation of the call (the list and its result objects, made before anything
    # is published) fails in turn: MemoryError, nothing changed, then the call returns the
    # reference. The call is positional (ORDERED, no outputs): pybind11's parsing of a call that
    # passes keywords does not check its allocations (limit 1 of the module docstring), and mode
    # and outputs are keyword-only. The binding's Python allocations do not depend on the mode;
    # a SIMULTANEOUS list's native allocations are swept in test_native_allocation.py.
    def make():
        w, actions = list_world(ORDERED)
        return ngc.apply_actions, (w.ledger, actions), w.state

    assert sweep_python_allocation(make) >= 4
    # The same lists with every output and in both modes return the same as their reference.
    for mode in (ORDERED, SIMULTANEOUS):
        w, actions = list_world(mode)
        twin, twin_actions = list_world(mode)
        got = plain(ngc.apply_actions(w.ledger, actions, mode=mode, outputs=ALL))
        assert got == plain(
            ngc.apply_actions(twin.ledger, twin_actions, mode=mode, outputs=ALL)
        )


@needs_hooks
def test_forwarding_bounds():
    # A list's joint network counts every participant: a Pinned program of two one-edge routes
    # on PAR(6, 4) (2 entries, 2 occurrences, 4 routers; the 4-unit link drops: 1 event, 3
    # dependencies), a Pinned program of two 2-hop tunnels on a diamond (2 entries, 4 occurrences)
    # and a pooled ShortestPaths demand of sources 6 and 7 (2 entries, 4 occurrences, 4 routers, 3
    # dependencies): 6 entries, 10 occurrences, 8 routers, 9 events, 6 dependencies.
    g = graph(
        10,
        [
            (0, 1, 6, 1),
            (0, 1, 4, 1),
            (2, 3, 100, 1),
            (2, 4, 100, 1),
            (3, 5, 100, 1),
            (4, 5, 100, 1),
            (6, 8, 100, 1),
            (7, 8, 100, 1),
            (8, 9, 100, 1),
        ],
    )
    w = World(g)
    routes = [ngc.RouteDefinition(0, dag_of(g, [e]), weight=1) for e in (0, 1)]
    a = w.pin(ngc.ForwardingProgram.from_routes(g, [0], [1], routes))
    b = w.pin(tunnels(g, [2], [5], [(2, [2, 4], 1), (2, [3, 5], 1)]))
    c = w.add([6, 7], 9, ecmp())
    actions = [ngc.Replace(a, 12), ngc.Replace(b, 2), ngc.Replace(c, 2)]

    def call(*listed):
        return ngc.apply_actions(w.ledger, list(listed), mode=SIMULTANEOUS)

    counts = [
        ("entry", 6),
        ("occurrence", 10),
        ("router", 8),
        ("event", 9),
        ("dependency", 6),
    ]
    for k, (name, n) in enumerate(counts):
        limits = [0] * 5
        limits[k] = n - 1
        before = w.state()
        value, error, structures = ext._with_forwarding_limits(
            call, tuple(actions), limits
        )
        assert value is None and type(error) is ValueError
        assert str(error) == f"forwarding {name} limit exceeded: {n} > {n - 1}"
        assert w.state() == before
        if name in ("event", "dependency"):
            assert structures == 0, name
        for one in actions:  # each participant alone stays within
            assert ext._with_forwarding_limits(call, (one,), limits)[1] is None, name
        limits[k] = n
        value, error, structures = ext._with_forwarding_limits(
            call, tuple(actions), limits
        )
        assert error is None and len(value) == 3 and structures == 1, name
    _, error, _ = ext._with_forwarding_limits(call, tuple(actions), [1] * 5)
    assert str(error) == "forwarding entry limit exceeded: 6 > 1"
    # The scalar call counts its own network: a's alone has 2 entries.
    _, error, _ = ext._with_forwarding_limits(a.replace, (12,), [1, 0, 0, 0, 0])
    assert str(error) == "forwarding entry limit exceeded: 2 > 1"


# -- the generators ----------------------------------------------------------------------------


class Generated:
    """A world built from a seed (built twice, identically): a random DAG on nodes 0..5 (sources 0
    and 1, destination 5; nonuniform costs, capacities 0 to 8) plus a separate link 6 -> 7, and one
    policy of each kind: ShortestPaths EQUAL LY (two sources), EQUAL LL, WEIGHTED LY with capacity,
    explicit and residual-snapshot weights, EQUAL LY under a mask; OPTIMIZED LL; CentralizedTE LL;
    T1 LY (two tunnels per head); T2 LY; Pinned LY over tunnels, over split routes, over from_dag
    under a mask, over a controller program with black holes; Pinned LY on the separate link. Then
    scalar calls, an orphan and an outside owner. None when a construction fails (both copies fail
    alike). g: the first copy's graph, shared by the second (records compare their routes)."""

    LOSSY_REPLACE = (
        0,
        2,
        3,
        4,
        5,
        10,
        11,
        12,
        13,
    )  # what a SIMULTANEOUS Replace accepts
    FORWARDABLE = (0, 2, 3, 4, 5, 8, 9)  # LOSSY computing policies
    SEPARATE = 14

    @classmethod
    def make(cls, seed, g=None):
        rng = random.Random(seed)
        rows = [
            (u, v, rng.randrange(9), 1 + rng.randrange(3))
            for u in range(6)
            for v in range(u + 1, 6)
            if rng.randrange(3) != 0
        ]
        rows.append((6, 7, 10, 1))
        try:
            return cls(rng, g or build(8, rows))
        except ValueError:
            return None

    def __init__(self, rng, g):
        self.w = w = World(g)
        n_edges = g.num_edges()
        sep = n_edges - 1  # 6 -> 7: the largest tail
        tails = g.edge_sources().tolist()
        heads_of = g.edge_destinations().tolist()

        def random_edge():
            return rng.randrange(n_edges - 1)

        def path_from(h):
            path, u = [], h
            while u != 5:
                outs = [e for e in range(sep) if tails[e] == u]
                if not outs:
                    return []
                e = rng.choice(outs)
                path.append(e)
                u = heads_of[e]
            return path

        S, D = [0, 1], [5]
        weighted = ngc.RouteAllocation.WEIGHTED
        w.add(S, D, ecmp())
        w.add([0], D, ecmp(LL))
        w.add(S, D, config(shortest(allocation=weighted), admission=LY))
        explicit = shortest(
            allocation=weighted,
            weight_basis=ngc.WeightBasis.EXPLICIT,
            edge_weights=[rng.randrange(4) for _ in range(n_edges)],
        )
        w.add(S, D, config(explicit, admission=LY))
        snapshot = shortest(
            allocation=weighted, weight_basis=ngc.WeightBasis.RESIDUAL_SNAPSHOT
        )
        w.add([0], D, config(snapshot, admission=LY))
        w.add(S, D, ecmp(), no_edges(g, [random_edge()]))
        w.add(S, D, config(ngc.ShortestPaths(), ngc.RouteAllocation.OPTIMIZED))
        w.add(S, D, controller())
        w.add(S, D, t1_config(2, LY))
        w.add(S, D, t2())
        specs = []
        for h in S:
            for _ in range(2):
                path = path_from(h)
                if path:
                    specs.append((h, path, 1 + rng.randrange(3)))
        if not specs:
            raise ValueError("no tunnel")
        heads = sorted({h for h, _, _ in specs})
        w.pin(tunnels(g, heads, D, specs))
        everything = dag_of(g, range(n_edges))
        split = [
            ngc.RouteDefinition(h, everything, weight=1 + rng.randrange(2))
            for h in heads
        ]
        w.pin(ngc.ForwardingProgram.from_routes(g, heads, D, split))
        w.pin(
            ngc.ForwardingProgram.from_dag(g, S, D, everything),
            mask=no_edges(g, [random_edge()]),
        )
        prog = ngc.FlowPolicy(ngc.FlowLedger(g), S, D, controller()).replace(6).program
        w.pin(prog, black_holes=no_edges(g, [random_edge()]))
        w.pin(tunnel(g, 6, 7, [sep]))
        # State before the list: scalar calls, an outside owner, an orphan.
        w.p[0].place(rng.randrange(6))
        w.p[2].replace(rng.randrange(6))
        w.p[6].place(rng.randrange(6))
        w.p[8].replace(2, design_demand=2)
        w.p[9].replace(2)
        w.p[10].replace(rng.randrange(6))
        w.add(S, D, ecmp()).place(1 + rng.randrange(3))
        w.add([0], D, ecmp()).place(1 + rng.randrange(3))
        w.p.pop()
        w.sources.pop()
        w.masks.pop()
        gc.collect()


def ordered_action(rng, w, i):
    """One random valid ORDERED action on policy i of a generated world."""
    p = w.p[i]
    routing, lossy = p.config.routing, p.config.admission == LY
    tunnel_model = isinstance(routing, ngc.DistributedTE) or (
        isinstance(routing, ngc.CentralizedTE) and routing.tunnel_count is not None
    )
    v = rng.randrange(13)
    kind = rng.randrange(3)
    if kind == 2 and not isinstance(routing, ngc.Pinned) and p.program is not None:
        bh = None
        if rng.randrange(3) == 0:
            bh = no_edges(w.g, [rng.randrange(w.g.num_edges() - 1)])
        return ngc.Forward(
            p, INF if not lossy and rng.randrange(4) == 0 else v, black_holes=bh
        )
    if kind == 0 and not tunnel_model:
        return ngc.Place(p, v)
    if isinstance(routing, ngc.DistributedTE):
        return ngc.Replace(p, v, design_demand=v + rng.randrange(3))
    return ngc.Replace(p, INF if not lossy and rng.randrange(4) == 0 else v)


def scalar(action, policy, outputs):
    """The action as a scalar call of `policy` (its own policy, or its twin's)."""
    if isinstance(action, ngc.Place):
        return policy.place(action.offer, outputs)
    if isinstance(action, ngc.Replace):
        return policy.replace(
            action.offer, design_demand=action.design_demand, outputs=outputs
        )
    return policy.forward(action.offer, black_holes=action.black_holes, outputs=outputs)


def index_of(w, policy):
    return next(k for k, p in enumerate(w.p) if p is policy)


def test_ordered_equals_scalar_calls():
    checked = forwards = holds = 0
    for seed in range(1, 201):
        a = Generated.make(seed)
        if a is None:
            continue
        b = Generated.make(seed, a.w.g)
        rng = random.Random(seed * 7919)
        order = list(range(len(a.w.p)))
        rng.shuffle(order)
        actions = [ordered_action(rng, a.w, i) for i in order[: 2 + rng.randrange(6)]]
        try:
            got, error = (
                plain(ngc.apply_actions(a.w.ledger, actions, outputs=ALL)),
                None,
            )
        except ValueError as e:
            got, error = None, str(e)
        want, scalar_error = [], None
        for x in actions:
            try:
                want.append(plain(scalar(x, b.w.p[index_of(a.w, x.policy)], ALL)))
            except ValueError as e:
                scalar_error = str(e)
                break
        assert error == scalar_error, seed
        if error is not None:
            continue
        assert got == want, seed
        assert a.w.state() == b.w.state(), seed
        forwards += sum(isinstance(x, ngc.Forward) for x in actions)
        holds += any(a.w.ledger.reserved().tolist())
        checked += 1
    assert checked > 150 and forwards > 80 and holds > 30, (checked, forwards, holds)


def test_one_participant_equals_the_scalar_call():
    checked = pooled = black = 0
    for seed in range(1, 301):
        a = Generated.make(seed)
        if a is None:
            continue
        b = Generated.make(seed, a.w.g)
        rng = random.Random(seed * 104729)
        forward = rng.randrange(2) == 0
        pool = Generated.FORWARDABLE if forward else Generated.LOSSY_REPLACE
        i = rng.choice(pool)
        v = rng.randrange(13)
        p = a.w.p[i]
        if forward:
            if p.program is None:
                p.replace(1)
                b.w.p[i].replace(1)
            bh = None
            if rng.randrange(2) == 0:
                bh = no_edges(a.w.g, [rng.randrange(a.w.g.num_edges() - 1)])
                black += 1
            actions = [ngc.Forward(p, v, black_holes=bh)]
        else:
            actions = [ngc.Replace(p, v)]
            pooled += i <= 5
        pair = rng.randrange(2) == 0
        if pair:
            separate = a.w.p[Generated.SEPARATE]
            actions.append(ngc.Replace(separate, rng.randrange(13)))
        got = plain(
            ngc.apply_actions(a.w.ledger, actions, mode=SIMULTANEOUS, outputs=ALL)
        )
        want = [plain(scalar(x, b.w.p[index_of(a.w, x.policy)], ALL)) for x in actions]
        # The scalar calls' outputs read the ledger after each; the list's after both: with the
        # pair, each entry's residual is the one after the list (the other one is disjoint).
        if pair:
            want[0]["residual"] = want[1]["residual"]
            want[0]["saturated_boundary_edges"] = b.w.ledger.saturated_boundary_edges(
                b.w.sources[i], b.w.masks[i]
            )
        assert got == want, (seed, i)
        assert a.w.state() == b.w.state(), (seed, i)
        checked += 1
    assert checked > 200 and pooled > 60 and black > 40, (checked, pooled, black)
