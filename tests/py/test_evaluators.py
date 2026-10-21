"""The evaluators' offer (docs/semantics.md, "10. Evaluators"): a finite offer for every model under
each admission it supports, in all three evaluators; EvaluationOptions keeps its offer and design
as given and the evaluators check them after the configuration, in the documented order; the
worked EVAL examples run on the library. PAR(a, b)
is two parallel links 0 -> 1 of capacities a and b, cost 1."""

from __future__ import annotations

import math

import numpy as np
import pytest
from helpers import build, exact, plain

import netgraph_core as ngc

INF = math.inf
A = ngc.RouteAllocation
S = ngc.SourceAllocation
LL = ngc.Admission.LOSSLESS
LY = ngc.Admission.LOSSY
T = ngc.Termination
L62 = 1 << 62


@pytest.mark.parametrize("workers", [1, 2, 4])
@pytest.mark.parametrize("bits", range(32))
def test_optimized_deletion_shortcuts_keep_complete_results(workers, bits):
    fixtures = [
        (
            5,
            [(0, 1, 3, 1), (0, 2, 3, 3), (1, 3, 3, 1), (2, 3, 3, 1), (3, 4, 3, 1)],
            [0],
            [4],
            3,
            9,
            [(0, 0), (2, 0), (4, 3)],
        ),
        (
            3,
            [(0, 1, 3, 1), (1, 2, 3, 2), (1, 2, 1, 5)],
            [0],
            [2],
            3,
            9,
            [(0, 3), (1, 2)],
        ),
        (3, [(0, 1, 2, 1), (1, 2, 2, 3)], [0], [2], 2, 8, [(0, 2), (1, 2)]),
        (
            4,
            [(0, 2, 1, 1), (1, 2, 1, 2), (2, 3, 1, 1)],
            [0, 1],
            [3],
            1,
            2,
            [(0, 0), (2, 1)],
        ),
        (
            3,
            [(0, 1, 3, 1), (1, 0, 1, 7), (1, 1, 1, 9), (1, 2, 3, 2), (1, 2, 1, 5)],
            [0],
            [2],
            3,
            9,
            [(0, 3), (3, 2)],
        ),
    ]
    outputs = ngc.OutputOptions(
        edge_flows=bool(bits & 1),
        residual=bool(bits & 2),
        saturated_boundary_edges=bool(bits & 4),
        min_cut=bool(bits & 8),
        edge_reservations=bool(bits & 16),
    )
    model = ngc.CentralizedTE(source_allocation=S.OPTIMIZED)
    for n, edges, sources, destinations, value, cost, reductions in fixtures:
        g = build(n, edges)
        options = opts(model, outputs=outputs)
        got = ngc.evaluate_edge_failures(g, sources, destinations, options, workers)
        baseline = ngc.evaluate_demand(g, sources, destinations, options)
        assert plain(got.baseline) == plain(baseline)
        assert (got.baseline.call_delivered, got.baseline.total_cost) == (value, cost)
        by_input = dict(reductions)
        expected = [
            (e, by_input[int(position)])
            for e, position in enumerate(g.input_positions())
            if int(position) in by_input
        ]
        assert [
            (entry.edge, entry.delivery_reduction) for entry in got.entries
        ] == expected
        for entry in got.entries:
            without = opts(
                model,
                eligibility=g.eligibility(excluded_edges=[entry.edge]),
                outputs=outputs,
            )
            fresh = ngc.evaluate_demand(g, sources, destinations, without)
            assert entry.delivery_reduction == value - fresh.call_delivered
            if fresh.edge_flows is not None:
                caps = list(g.edge_capacity_units())
                flows = list(fresh.edge_flows)
                assert all(
                    0 <= flow <= cap for flow, cap in zip(flows, caps, strict=True)
                )


def test_optimized_deletion_shortcuts_masks_finite_zero_and_wide_cost():
    model = ngc.CentralizedTE(source_allocation=S.OPTIMIZED)
    outputs = ngc.OutputOptions(
        edge_flows=True,
        residual=True,
        min_cut=True,
        saturated_boundary_edges=True,
        edge_reservations=True,
    )
    g = build(3, [(0, 1, 3, 1), (1, 2, 3, 2), (1, 2, 1, 5)])
    for eligibility in [
        None,
        g.eligibility(excluded_edges=[2]),
        g.eligibility(excluded_nodes=[1]),
    ]:
        for offer in [None, 0, 1, 3, 4]:
            options = opts(model, eligibility=eligibility, outputs=outputs, offer=offer)
            for workers in [1, 2]:
                got = ngc.evaluate_edge_failures(g, [0], [2], options, workers)
                assert plain(got.baseline) == plain(
                    ngc.evaluate_demand(g, [0], [2], options)
                )
                for entry in got.entries:
                    deleted = g.eligibility(excluded_edges=[entry.edge])
                    if eligibility is not None:
                        deleted = ngc.Eligibility.intersect(eligibility, deleted)
                    fresh = ngc.evaluate_demand(
                        g, [0], [2], opts(model, eligibility=deleted, offer=offer)
                    )
                    assert (
                        entry.delivery_reduction
                        == got.baseline.call_delivered - fresh.call_delivered
                    )
    amount = (1 << 32) + 1
    wide = build(2, [(0, 1, amount, amount)])
    got = ngc.evaluate_edge_failures(wide, 0, 1, opts(model, outputs=outputs), 2)
    assert got.baseline.total_cost == (1 << 64) + (1 << 33) + 1
    assert [(e.edge, e.delivery_reduction) for e in got.entries] == [(0, amount)]


DOMAIN = "offer must be an integer in [0, 2^62)"


def par(a, b):
    return build(2, [(0, 1, a, 1), (0, 1, b, 1)])


def opts(routing, admission=LL, **kw):
    return ngc.EvaluationOptions(
        ngc.FlowPolicyConfig(routing=routing, admission=admission), **kw
    )


def demand(g, routing, admission=LL, **kw):
    return ngc.evaluate_demand(g, 0, 1, opts(routing, admission, **kw))


def failures(g, routing, admission=LL, workers=1, **kw):
    s = ngc.evaluate_edge_failures(g, 0, 1, opts(routing, admission, **kw), workers)
    return s.baseline.call_delivered, [x.delivery_reduction for x in s.entries]


def dte(count=2):
    return ngc.DistributedTE(ngc.TunnelCount(count))


def t2(count=2):
    return ngc.CentralizedTE(
        source_allocation=S.EQUAL, tunnel_count=ngc.TunnelCount(count)
    )


def two_paths(g):
    # The program of two tunnels, one on each link.
    routes = [ngc.RouteDefinition(0, ngc.Path(g, [[e]])) for e in (0, 1)]
    return ngc.Pinned(ngc.ForwardingProgram.from_routes(g, 0, 1, routes))


# -- construction -----------------------------------------------------------------------------


def test_options_construction():
    # The offer and the design are Units kept as given: a negative, an infinite or a float value
    # constructs and reads back unchanged (the evaluator checks it); another type is a TypeError.
    cfg = ngc.FlowPolicyConfig(routing=ngc.ShortestPaths())
    for value in (-1, INF, 4.0, 2**70):
        assert ngc.EvaluationOptions(cfg, offer=value).offer == value
    assert ngc.EvaluationOptions(cfg, design_demand=-1).design_demand == -1
    assert ngc.EvaluationOptions(cfg, offer=np.int64(3)).offer == 3
    for bad in ("4", True, np.bool_(True), [4]):
        with pytest.raises(TypeError, match=r"^offer must be an int or a float"):
            ngc.EvaluationOptions(cfg, offer=bad)
        with pytest.raises(
            TypeError, match=r"^design_demand must be an int or a float"
        ):
            ngc.EvaluationOptions(cfg, design_demand=bad)
    with pytest.raises(TypeError):
        ngc.EvaluationOptions(cfg, offered=4)


# -- one model at a time (E1 to E7) ----------------------------------------------------------


def test_offer_shortest_paths_lossless():
    # ECMP on PAR(4, 10): the uncapped bound 8 (each link 4); an offer below it is delivered in
    # full, above it up to 8; a zero offer delivers nothing.
    g = par(4, 10)
    r = demand(g, ngc.ShortestPaths())
    assert (r.call_delivered, r.requested, r.termination) == (8, INF, T.EXHAUSTED)
    r = demand(g, ngc.ShortestPaths(), offer=5)
    assert (r.call_delivered, r.requested, r.call_undelivered, r.termination) == (
        5,
        5,
        0,
        T.COMPLETE,
    )
    r = demand(g, ngc.ShortestPaths(), offer=12)
    assert (r.call_delivered, r.call_undelivered, r.termination) == (8, 4, T.EXHAUSTED)
    r = demand(g, ngc.ShortestPaths(), offer=0)
    assert (r.call_delivered, r.termination) == (0, T.COMPLETE)


def test_offer_shortest_paths_lossy():
    # ECMP Lossy on PAR(4, 10): the offer 12 splits 6 and 6, link 0 drops 2. Without an offer,
    # the first-hop experiment: each first-hop link injects its own free capacity (4 and 10),
    # offered 14 against requested +inf, nothing dropped.
    g = par(4, 10)
    r = demand(g, ngc.ShortestPaths(), LY, offer=12)
    assert (r.call_delivered, r.call_edge_drops, r.by_source[0].offered) == (
        10,
        [(0, 2)],
        12,
    )
    assert r.termination == T.COMPLETE
    r = demand(g, ngc.ShortestPaths(), LY)
    assert (r.call_delivered, r.requested, r.by_source[0].offered) == (14, INF, 14)
    assert r.call_edge_drops == []
    r = demand(g, ngc.ShortestPaths(), LY, offer=0)
    assert (r.call_delivered, r.by_source[0].offered) == (0, 0)


def test_offer_weighted():
    # WCMP by capacity on PAR(4, 10): 14 uncapped; the offer 7 splits 2 and 5.
    g = par(4, 10)
    wcmp = ngc.ShortestPaths(allocation=A.WEIGHTED)
    assert demand(g, wcmp).call_delivered == 14
    r = demand(g, wcmp, offer=7, outputs=ngc.OutputOptions(edge_flows=True))
    assert r.call_delivered == 7 and r.edge_flows.tolist() == [2, 5]


def test_offer_controllers():
    # The controllers deliver the optimum up to the offer at cost 1 per unit; CentralizedTE's
    # min_cut exists only when the solve reaches the maximum; LOSSY is refused by the matrix.
    g = par(4, 10)
    cut = ngc.OutputOptions(min_cut=True)
    for model in (
        ngc.ShortestPaths(allocation=A.OPTIMIZED),
        ngc.CentralizedTE(source_allocation=S.OPTIMIZED),
    ):
        r = demand(g, model, outputs=cut)
        assert (r.call_delivered, r.total_cost) == (14, 14)
        assert r.min_cut == ([0, 1] if isinstance(model, ngc.CentralizedTE) else None)
        r = demand(g, model, offer=5, outputs=cut)
        assert (r.call_delivered, r.total_cost, r.min_cut) == (5, 5, None)
        with pytest.raises(ValueError, match="requires admission=LOSSLESS"):
            demand(g, model, LY, offer=5)


def test_offer_distributed():
    # T1, two tunnels sized for the design 8 (requests 4 and 4: one on each link): Lossless
    # bound 8, an offer 5 delivered in full; Lossy 12 splits 6 and 6, link 0 drops 2.
    g = par(4, 10)
    assert demand(g, dte(), design_demand=8).call_delivered == 8
    assert demand(g, dte(), design_demand=8, offer=5).call_delivered == 5
    r = demand(g, dte(), LY, design_demand=8, offer=12)
    assert [list(t.edges) for t in r.program.tunnels] == [[0], [1]]
    assert (r.call_delivered, r.call_edge_drops) == (10, [(0, 2)])
    with pytest.raises(
        ValueError,
        match=exact("DistributedTE LOSSY evaluation requires a finite offer"),
    ):
        demand(g, dte(), LY, design_demand=8)
    with pytest.raises(
        ValueError, match=exact("DistributedTE evaluation requires design_demand")
    ):
        demand(g, dte())


def test_offer_centralized_tunnels():
    # T2 on PAR(4, 10), N = 2: without a design U = 5, both tunnels on link 1, the bound 10; an
    # offer 5 is also the design; Lossy 12 (the design 12, u 5) drops 2 on link 1; Lossy with
    # the design 4 and no offer forwards 4 (u 2, both tunnels on link 0).
    g = par(4, 10)
    r = demand(g, t2())
    assert (r.call_delivered, r.plan.max_bandwidth_per_tunnel) == (10, 5)
    assert [list(t.edges) for t in r.program.tunnels] == [[1], [1]]
    r = demand(g, t2(), offer=5)
    assert (r.call_delivered, r.plan.design_demand) == (5, 5)
    r = demand(g, t2(), LY, offer=12)
    assert (r.call_delivered, r.call_edge_drops) == (10, [(1, 2)])
    assert (r.plan.design_demand, r.plan.planned_bandwidth_per_tunnel) == (12, 5)
    r = demand(g, t2(), LY, design_demand=4)
    assert (r.call_delivered, r.plan.planned_bandwidth_per_tunnel) == (4, 2)
    assert [list(t.edges) for t in r.program.tunnels] == [[0], [0]]
    with pytest.raises(
        ValueError,
        match=exact("CentralizedTE LOSSY evaluation requires a finite offer"),
    ):
        demand(g, t2(), LY)


def test_offer_pinned():
    # Pinned over one tunnel per link: the bound 8 (equal split, link 0 holds 4), an offer 5 in
    # full, Lossy 12 delivers 10.
    g = par(4, 10)
    pin = two_paths(g)
    assert demand(g, pin).call_delivered == 8
    assert demand(g, pin, offer=5).call_delivered == 5
    assert demand(g, pin, LY, offer=12).call_delivered == 10
    with pytest.raises(
        ValueError, match=exact("Pinned LOSSY evaluation requires a finite offer")
    ):
        demand(g, pin, LY)


# -- the three evaluators (E8 to E12) -------------------------------------------------------


def test_failures_fixed_offer():
    # EVAL-fixed-offer and EVAL-first-hop on PAR(10, 10); then PAR(4, 10) ECMP Lossy at the
    # offers 4 and 12 (every deletion answers the same offer) and without one (the first-hop
    # experiment injects each remaining link's capacity, a different question).
    for workers in (1, 2):
        assert failures(par(10, 10), ngc.ShortestPaths(), LY, workers, offer=4) == (
            4,
            [0, 0],
        )
        assert failures(par(10, 10), ngc.ShortestPaths(), LY, workers) == (20, [10, 10])
        g = par(4, 10)
        assert failures(g, ngc.ShortestPaths(), LY, workers, offer=4) == (4, [0, 0])
        assert failures(g, ngc.ShortestPaths(), LY, workers, offer=12) == (10, [0, 6])
        assert failures(g, ngc.ShortestPaths(), LY, workers) == (14, [4, 10])


def test_failures_distributed_offer():
    # T1, design 8, Lossy 12: the baseline 10; deleting link 0 sets both tunnels up on link 1
    # (10, reduction 0); deleting link 1 leaves one tunnel on link 0 carrying 4 (reduction 6).
    g = par(4, 10)
    assert failures(g, dte(), LY, design_demand=8, offer=12) == (10, [0, 6])


def test_pairs_offer():
    # The offer is the offer of every pair; nothing is shared between pairs.
    g = par(4, 10)
    rs = ngc.evaluate_pairs(g, [(0, 1), (0, 1)], opts(ngc.ShortestPaths(), offer=5))
    assert [r.call_delivered for r in rs] == [5, 5]


def test_offer_zero_and_masked():
    # Both links masked: Lossy's offer 4 is the demand's source drops (its row zero), Lossless
    # leaves it undelivered; NO_ROUTE both. A zero offer is COMPLETE for every model, T1's plan
    # still set up (both tunnels up, 8 reserved).
    g = par(4, 10)
    mask = g.eligibility(excluded_edges=[0, 1])
    r = demand(g, ngc.ShortestPaths(), LY, eligibility=mask, offer=4)
    assert (r.termination, r.call_delivered, r.call_source_drops) == (T.NO_ROUTE, 0, 4)
    assert [(s.offered, s.source_drops) for s in r.by_source] == [(0, 0)]
    r = demand(g, ngc.ShortestPaths(), eligibility=mask, offer=4)
    assert (r.termination, r.call_delivered, r.call_undelivered) == (T.NO_ROUTE, 0, 4)
    pin = two_paths(g)
    for model, admission, kw in (
        (ngc.ShortestPaths(), LL, {}),
        (ngc.ShortestPaths(), LY, {}),
        (ngc.ShortestPaths(allocation=A.OPTIMIZED), LL, {}),
        (ngc.CentralizedTE(source_allocation=S.OPTIMIZED), LL, {}),
        (dte(), LL, {"design_demand": 8}),
        (dte(), LY, {"design_demand": 8}),
        (t2(), LL, {}),
        (t2(), LY, {}),
        (pin, LL, {}),
        (pin, LY, {}),
    ):
        r = demand(g, model, admission, offer=0, **kw)
        assert (r.termination, r.call_delivered) == (T.COMPLETE, 0), model
    r = demand(g, dte(), LY, design_demand=8, offer=0)
    assert r.plan.total_reserved_bandwidth == 8
    assert [t.status for t in r.program.tunnels] == [ngc.TunnelStatus.UP] * 2


def test_offer_limits():
    # One link of 2^62 - 1: every model delivers the largest offer (T1 with that design); 2^62 is
    # refused; deleting the link takes it all.
    big = L62 - 1
    g = build(2, [(0, 1, big, 1)])
    route = [ngc.RouteDefinition(0, ngc.Path(g, [[0]]))]
    pin = ngc.Pinned(ngc.ForwardingProgram.from_routes(g, 0, 1, route))
    for model, admission, kw in (
        (ngc.ShortestPaths(), LL, {}),
        (ngc.ShortestPaths(), LY, {}),
        (ngc.ShortestPaths(allocation=A.OPTIMIZED), LL, {}),
        (ngc.CentralizedTE(source_allocation=S.OPTIMIZED), LL, {}),
        (dte(1), LL, {"design_demand": big}),
        (dte(1), LY, {"design_demand": big}),
        (t2(1), LL, {}),
        (t2(1), LY, {}),
        (pin, LL, {}),
        (pin, LY, {}),
    ):
        assert demand(g, model, admission, offer=big, **kw).call_delivered == big, model
    with pytest.raises(ValueError, match=exact(DOMAIN)):
        demand(g, ngc.ShortestPaths(), offer=L62)
    assert failures(g, ngc.ShortestPaths(), offer=big) == (big, [big])


# -- the order of the checks (V9 to V11) -----------------------------------------------------


def outcome(model, kwargs, admission="LOSSLESS", offer=None, design_demand=None):
    """One EVAL case on the library: the first error, or the offer and design it evaluates."""
    g = build(3, [(0, 1, 5, 1), (0, 2, 5, 1), (1, 2, 5, 1)])
    if model == "Pinned":
        route = [ngc.RouteDefinition(0, ngc.Path(g, [[1]]))]
        routing = ngc.Pinned(ngc.ForwardingProgram.from_routes(g, 0, 2, route))
    else:
        fields = dict(kwargs)
        if "source_allocation" in fields:
            fields["source_allocation"] = getattr(S, fields["source_allocation"])
        if "tunnel_count" in fields:
            fields["tunnel_count"] = ngc.TunnelCount(int(fields["tunnel_count"]))
        if model == "DistributedTE":
            routing = ngc.DistributedTE(fields.pop("tunnel_count"), **fields)
        else:
            routing = getattr(ngc, model)(**fields)
    try:
        options = ngc.EvaluationOptions(
            ngc.FlowPolicyConfig(
                routing=routing, admission=getattr(ngc.Admission, admission)
            ),
            offer=offer,
            design_demand=design_demand,
        )
    except TypeError:
        return "TypeError"
    try:
        r = ngc.evaluate_demand(g, 0, 2, options)
    except ValueError as e:
        return str(e)
    if r.requested == INF and r.by_source[0].offered is not None:
        evaluated = "first-hop experiment"
    else:
        evaluated = r.requested
    design = None if r.plan is None else r.plan.design_demand
    return f"offer={evaluated} design_demand={design}"


def test_first_error():
    # EVAL-first-error: the configuration, an explicit offer, the design, T1's missing design,
    # the missing Lossy offer.
    cases = [
        dict(
            model="CentralizedTE",
            kwargs={"source_allocation": "EQUAL", "prefer_headroom": True},
            offer=-1,
        ),
        dict(
            model="DistributedTE",
            kwargs={"tunnel_count": "4"},
            admission="LOSSY",
            offer=-1,
            design_demand=INF,
        ),
        dict(
            model="DistributedTE",
            kwargs={"tunnel_count": "4"},
            admission="LOSSY",
            offer=5,
            design_demand=INF,
        ),
        dict(model="DistributedTE", kwargs={"tunnel_count": "4"}, admission="LOSSY"),
        dict(
            model="DistributedTE",
            kwargs={"tunnel_count": "4"},
            admission="LOSSY",
            design_demand=10,
        ),
        dict(model="ShortestPaths", kwargs={}, offer=INF, design_demand=5),
        dict(model="ShortestPaths", kwargs={}, design_demand=5),
        dict(
            model="CentralizedTE",
            kwargs={"source_allocation": "EQUAL", "tunnel_count": "16"},
            admission="LOSSY",
        ),
        dict(
            model="CentralizedTE",
            kwargs={"source_allocation": "EQUAL", "tunnel_count": "16"},
            admission="LOSSY",
            design_demand=10,
        ),
        dict(model="Pinned", kwargs={}, admission="LOSSY"),
        dict(model="ShortestPaths", kwargs={}, admission="LOSSY"),
        dict(model="ShortestPaths", kwargs={}, admission="LOSSY", offer=4),
    ]
    assert [outcome(**c) for c in cases] == [
        "prefer_headroom requires tunnel_count",
        DOMAIN,
        "DistributedTE requires a finite design demand",
        "DistributedTE evaluation requires design_demand",
        "DistributedTE LOSSY evaluation requires a finite offer",
        DOMAIN,
        "ShortestPaths does not accept design_demand",
        "CentralizedTE LOSSY evaluation requires a finite offer",
        "offer=10 design_demand=10",
        "Pinned LOSSY evaluation requires a finite offer",
        "offer=first-hop experiment design_demand=None",
        "offer=4 design_demand=None",
    ]


def test_offer_infinite():
    # EVAL-offer-infinite: +inf constructs and is refused by the evaluator after the
    # configuration; a string is a construction TypeError.
    headroom = {"source_allocation": "EQUAL", "prefer_headroom": True}
    assert [
        outcome("CentralizedTE", headroom, offer=INF),
        outcome("CentralizedTE", headroom, offer="4"),
        outcome("ShortestPaths", {}, offer=INF),
    ] == ["prefer_headroom requires tunnel_count", "TypeError", DOMAIN]


def test_design_model_before_range():
    # A model without a design refuses one before its range is read.
    assert (
        outcome("ShortestPaths", {}, design_demand=-1)
        == "ShortestPaths does not accept design_demand"
    )
    assert (
        outcome("DistributedTE", {"tunnel_count": "4"}, design_demand=-1)
        == "design_demand must be an integer in [0, 2^62) or permitted +inf"
    )


def test_pairs_order():
    # Pinned first (its offer is never read); the offer before per_pair_eligibility's length;
    # then each pair in order.
    g = par(4, 10)
    with pytest.raises(
        ValueError, match=exact("evaluate_pairs does not accept Pinned")
    ):
        ngc.evaluate_pairs(g, [(0, 1)], opts(two_paths(g), LY, offer=-1))
    with pytest.raises(ValueError, match=exact(DOMAIN)):
        ngc.evaluate_pairs(
            g, [(0, 1)], opts(ngc.ShortestPaths(), offer=-1), per_pair_eligibility=[]
        )
    with pytest.raises(
        ValueError,
        match=exact(
            "node 1 is in both endpoint sets: sources and destinations must not overlap"
        ),
    ):
        ngc.evaluate_pairs(g, [(0, 1), (1, 1)], opts(ngc.ShortestPaths(), offer=4))


# -- every model at finite offers: the failure study against its definition -----------------


def _models(g):
    # Every model the evaluators take, with its admission and design (T1 always names one).
    route = ngc.Path(g, [[0], [2]])
    pin = ngc.Pinned(
        ngc.ForwardingProgram.from_routes(
            g,
            0,
            3,
            [
                ngc.RouteDefinition(0, route),
                ngc.RouteDefinition(0, ngc.Path(g, [[1], [3]])),
            ],
        )
    )
    return [
        ("ECMP", ngc.ShortestPaths(), (LL, LY), {}),
        ("WCMP", ngc.ShortestPaths(allocation=A.WEIGHTED), (LL, LY), {}),
        ("SP OPTIMIZED", ngc.ShortestPaths(allocation=A.OPTIMIZED), (LL,), {}),
        ("CTE OPTIMIZED", ngc.CentralizedTE(source_allocation=S.OPTIMIZED), (LL,), {}),
        ("CTE EQUAL", ngc.CentralizedTE(source_allocation=S.EQUAL), (LL,), {}),
        (
            "CTE support",
            ngc.CentralizedTE(source_allocation=S.OPTIMIZED, support=route.to_dag()),
            (LL,),
            {},
        ),
        ("T1", dte(2), (LL, LY), {"design_demand": 8}),
        ("T2", t2(2), (LL, LY), {}),
        ("Pinned", pin, (LL, LY), {}),
    ]


def test_failures_on_random_graphs_at_an_offer():
    # 300 random graphs (4 to 7 nodes, 6 to 15 edges, capacities 0..8, costs 1..3, one or two
    # sources, a random offer 0..14), eight configurations each: every deletion's reduction is the
    # baseline minus evaluate_demand with that edge excluded too (an edge left out of the
    # candidates: 0), on 1 worker for even trials and 2 for odd ones. Optimized CentralizedTE
    # repairs its baseline when the offer exceeds the maximum and solves afresh when the offer
    # binds below it: both happen here.
    rng = np.random.default_rng(7)
    repaired = fresh = rows = 0
    for trial in range(300):
        n = int(rng.integers(4, 8))
        m = int(rng.integers(6, 16))
        es, ed = [], []
        while len(es) < m:
            u, v = int(rng.integers(0, n)), int(rng.integers(0, n))
            if u != v:
                es.append(u)
                ed.append(v)
        cap = [int(x) for x in rng.integers(0, 9, m)]
        cost = [int(x) for x in rng.integers(1, 4, m)]
        g = ngc.Graph(n, es, ed, cap, cost, 1.0)
        nodes = [int(x) for x in rng.permutation(n)]
        src = sorted(nodes[: int(rng.integers(1, 3))])
        dst = [nodes[-1]]
        offer = int(rng.integers(0, 15))
        workers = 1 + trial % 2
        configs = [
            (ngc.CentralizedTE(source_allocation=S.OPTIMIZED), LL, {"offer": offer}),
            (ngc.CentralizedTE(source_allocation=S.EQUAL), LL, {"offer": offer}),
            (ngc.ShortestPaths(allocation=A.OPTIMIZED), LL, {"offer": offer}),
            (ngc.ShortestPaths(), LY, {"offer": offer}),
            (ngc.ShortestPaths(), LL, {"offer": offer}),
            (dte(2), LY, {"offer": offer, "design_demand": offer}),
            (t2(2), LY, {"offer": offer}),
            (ngc.CentralizedTE(source_allocation=S.OPTIMIZED), LL, {}),
        ]
        for k, (model, admission, kw) in enumerate(configs):
            o = opts(model, admission, **kw)
            s = ngc.evaluate_edge_failures(g, src, dst, o, workers)
            base = ngc.evaluate_demand(g, src, dst, o).call_delivered
            assert s.baseline.call_delivered == base, (trial, k)
            if k == 0 and any(x.delivery_reduction for x in s.entries):
                top = ngc.evaluate_demand(g, src, dst, opts(model, LL)).call_delivered
                repaired += base < offer
                fresh += base == offer < top
            got = {x.edge: x.delivery_reduction for x in s.entries}
            for e in range(g.num_edges()):
                deleted = opts(
                    model,
                    admission,
                    eligibility=g.eligibility(excluded_edges=[e]),
                    **kw,
                )
                want = base - ngc.evaluate_demand(g, src, dst, deleted).call_delivered
                assert got.get(e, 0) == want, (trial, k, e, workers)
                rows += 1
    assert repaired >= 10 and fresh >= 10, (repaired, fresh)
    assert rows > 20000


def test_failures_at_every_offer():
    # A diamond 0 -> {1, 2} -> 3, capacities 4, 6, 3, 9 and costs 1, 1, 2, 1: the maximum is 9.
    # For every model, admission and offer (below, at and above its maximum, and none where a
    # model takes none), evaluate_edge_failures gives the baseline of evaluate_demand and, per
    # allowed edge, the baseline minus evaluate_demand with that edge excluded too (an edge left
    # out of the candidates must have reduction 0); evaluate_pairs agrees with evaluate_demand.
    g = build(4, [(0, 1, 4, 1), (0, 2, 6, 1), (1, 3, 3, 2), (2, 3, 9, 1)])
    checked = 0
    for name, model, admissions, extra in _models(g):
        for admission in admissions:
            top = ngc.evaluate_demand(
                g,
                0,
                3,
                opts(model, LL, **extra)
                if LL in admissions
                else opts(model, admission, offer=L62 - 1, **extra),
            ).call_delivered
            offers = [1, top, top + 5] + ([None] if admission == LL else [])
            for offer in offers:
                o = opts(model, admission, offer=offer, **extra)
                base = ngc.evaluate_demand(g, 0, 3, o).call_delivered
                for workers in (1, 2):
                    s = ngc.evaluate_edge_failures(g, 0, 3, o, workers)
                    assert s.baseline.call_delivered == base, (name, admission, offer)
                    got = {x.edge: x.delivery_reduction for x in s.entries}
                    for e in range(g.num_edges()):
                        deleted = ngc.evaluate_demand(
                            g,
                            0,
                            3,
                            opts(
                                model,
                                admission,
                                offer=offer,
                                eligibility=g.eligibility(excluded_edges=[e]),
                                **extra,
                            ),
                        ).call_delivered
                        assert got.get(e, 0) == base - deleted, (
                            name,
                            admission,
                            offer,
                            workers,
                            e,
                        )
                    checked += 1
                if not isinstance(model, ngc.Pinned):
                    (pair,) = ngc.evaluate_pairs(g, [(0, 3)], o, workers=2)
                    assert pair.call_delivered == base, (name, admission, offer)
    assert checked > 100


def test_invalid_offer_precedes_frozen_cache_miss():
    g = par(5, 7)
    cache = ngc.ShortestPathCache(g)
    cache.freeze()
    invalid = opts(ngc.ShortestPaths(), cache=cache, offer=INF)
    valid = opts(ngc.ShortestPaths(), cache=cache, offer=1)
    calls = [
        lambda o: ngc.evaluate_demand(g, 0, 1, o),
        lambda o: ngc.evaluate_pairs(g, [(0, 1)], o, workers=2),
        lambda o: ngc.evaluate_edge_failures(g, 0, 1, o, workers=2),
    ]
    for call in calls:
        before = cache.stats()
        with pytest.raises(ValueError) as error:
            call(invalid)
        assert str(error.value) == DOMAIN
        after = cache.stats()
        assert (after.hits, after.misses, after.entries) == (
            before.hits,
            before.misses,
            before.entries,
        )
        with pytest.raises(ngc.ShortestPathCacheMiss):
            call(valid)
        after = cache.stats()
        assert (
            after.misses == before.misses + 1
            and after.hits == before.hits
            and after.entries == 0
        )
    cache.thaw()
    result = ngc.evaluate_demand(g, 0, 1, valid)
    assert result.call_delivered == 1


@pytest.mark.parametrize("model", ["distributed", "centralized", "support"])
def test_pair_rosters_use_the_composed_mask(model):
    """Every reused roster has the same endpoints and effective mask as its evaluation."""
    import random

    from helpers import ALL_OUT, dag_of, plain

    def routing(g):
        if model == "distributed":
            return dte(2)
        return ngc.CentralizedTE(
            source_allocation=S.EQUAL,
            tunnel_count=ngc.TunnelCount(2),
            support=dag_of(g, range(g.num_edges())) if model == "support" else None,
        )

    # A removed source (or its only edge) leaves one head: the full offer 4 goes
    # to that head. Using the unmasked two-head roster would lose its share.
    g = build(4, [(0, 3, 4, 1), (1, 3, 8, 1)])
    masks = [
        g.eligibility(excluded_nodes=[0]),
        g.eligibility(excluded_edges=[0]),
        g.eligibility(excluded_nodes=[0, 1]),
    ]
    common = g.eligibility(excluded_nodes=[2])
    o = opts(routing(g), eligibility=common, design_demand=6, offer=4, outputs=ALL_OUT)
    result = ngc.evaluate_pairs(g, [([0, 1], 3)] * 3, o, masks, 2)
    assert [r.delivered for r in result] == [4, 4, 0]
    assert [r.edge_flows.tolist() for r in result] == [[0, 4], [0, 4], [0, 0]]
    zero = opts(
        routing(g), eligibility=common, design_demand=6, offer=0, outputs=ALL_OUT
    )
    assert all(
        r.delivered == 0
        for r in ngc.evaluate_pairs(g, [([0, 1], 3)] * 3, zero, masks, 2)
    )

    rng = random.Random(7)
    compared = 0
    for _ in range(150):
        n = rng.randint(5, 9)
        edges = [
            (u, v, rng.randint(1, 9), rng.randint(1, 4))
            for u in range(n)
            for v in range(u + 1, n)
            if rng.random() < 0.45
        ]
        if not edges:
            continue
        g = build(n, edges)
        common = (
            g.eligibility(excluded_edges=[rng.randrange(g.num_edges())])
            if rng.random() < 0.5
            else None
        )
        pairs, per = [], []
        for _ in range(3):
            nodes = list(range(n - 2))
            rng.shuffle(nodes)
            heads = sorted(nodes[: rng.randint(2, 3)])
            pairs.append((heads, [rng.choice([n - 2, n - 1])]))
            per.append(
                g.eligibility(
                    excluded_nodes=[h for h in heads if rng.random() < 0.4],
                    excluded_edges=[
                        e for e in range(g.num_edges()) if rng.random() < 0.25
                    ],
                )
            )
        config = ngc.FlowPolicyConfig(routing=routing(g))
        options = ngc.EvaluationOptions(
            config, eligibility=common, design_demand=6, offer=4, outputs=ALL_OUT
        )
        batch = ngc.evaluate_pairs(g, pairs, options, per, 2)
        assert len(batch) == 3
        for i, (heads, destinations) in enumerate(pairs):
            mask = (
                per[i] if common is None else ngc.Eligibility.intersect(common, per[i])
            )
            single = ngc.EvaluationOptions(
                config, eligibility=mask, design_demand=6, offer=4, outputs=ALL_OUT
            )
            assert plain(batch[i]) == plain(
                ngc.evaluate_demand(g, heads, destinations, single)
            )
            compared += 1
    assert compared >= 400


def test_tunnel_pair_search_error_precedes_later_conversion_error():
    # floor((2^30-1)/10000)+1 reachable heads exceed the tunnel limit.
    heads = ((1 << 30) - 1) // 10000 + 1
    g = build(heads + 1, [(h, heads, 1, 1) for h in range(heads)])
    expected = f"tunnel limit exceeded: {heads} heads x 10000 tunnels > {(1 << 30) - 1}"
    for model in (dte(10000), t2(10000)):
        options = opts(model, design_demand=1, offer=1)
        for invalid in (None, (0,), ([], heads), (0, 0)):
            with pytest.raises(ValueError, match=exact(expected)):
                ngc.evaluate_pairs(
                    g, [(list(range(heads)), heads), invalid], options, workers=2
                )
        # If the earlier pair fails its linear checks, no later search takes precedence.
        with pytest.raises(ValueError, match=exact("sources: the node set is empty")):
            ngc.evaluate_pairs(
                g, [([], heads), (list(range(heads)), heads)], options, workers=2
            )
