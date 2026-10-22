"""The policy engine, the ShortestPathCache and the evaluators through the bindings (docs/semantics.md):
configuration objects and their structural validation, PlacementResult fields and outputs,
lifetimes, the evaluators with workers, and the cache statistics. The fixtures of the reference
(``tests/reference/policy/``) replay in ``test_flow_policy_fixtures.py``."""

from __future__ import annotations

import gc
import weakref

import _netgraph_core as ext
import numpy as np
import pytest
from helpers import (
    ALL_OUT,
    GLOBAL,
    HIGHER_AVAIL,
    INF,
    L62,
    SHORTEST_PATHS,
    assert_read_only,
    build,
    config,
    dag_of,
    exact,
    full_ledger_state,
    ids,
    plain,
    sweep_python_allocation,
    t1_config,
    to_ids,
)
from private_hooks import HOOKS, needs_hooks

import netgraph_core as ngc

T = ngc.Termination


def diamond(a=5, b=3, c=1):
    # input edges: 0: 0->1, 1: 1->3 (cap a, cost 1); 2: 0->2, 3: 2->3 (cap b, cost c)
    return build(4, [(0, 1, a, 1), (1, 3, a, 1), (0, 2, b, c), (2, 3, b, c)])


def test_configurations_are_structural_and_explicit():
    with pytest.raises(TypeError):
        ngc.FlowPolicyConfig(ngc.ShortestPaths(), ngc.Admission.LOSSLESS)
    # Keyword-only; the routing is required, admission LOSSLESS by default. The
    # allocations, the weights and the route choice are the model's: no keyword of the config.
    # A wrong class is TypeError.
    plain = ngc.FlowPolicyConfig(routing=ngc.ShortestPaths())
    assert plain.admission == ngc.Admission.LOSSLESS
    with pytest.raises(TypeError):
        ngc.FlowPolicyConfig(admission=ngc.Admission.LOSSLESS)
    for kw in (
        "selection",
        "placement",
        "origin",
        "weights",
        "edge_weights",
        "bundles",
    ):
        with pytest.raises(TypeError):
            ngc.FlowPolicyConfig(routing=ngc.ShortestPaths(), **{kw: None})
    with pytest.raises(TypeError):
        ngc.ShortestPaths(allocation=1)
    with pytest.raises(TypeError):
        ngc.ShortestPaths(source_allocation="EQUAL")
    with pytest.raises(TypeError):
        ngc.ShortestPaths(weight_basis="EXPLICIT")
    with pytest.raises(TypeError, match="incompatible"):
        ngc.ShortestPaths(
            tunnels=ngc.TunnelCount(2)
        )  # a field outside its routing variant: TypeError
    with pytest.raises(TypeError, match="incompatible"):
        ngc.CentralizedTE(
            source_allocation=ngc.SourceAllocation.OPTIMIZED,
            cost_ceiling=ngc.AbsoluteCost(4),
        )
    with pytest.raises(TypeError, match="routing"):
        config("ShortestPaths")
    # T1: the tunnel count (a TunnelCount), the ceiling and the search minimum, every field as built.
    ca = ngc.DistributedTE(
        ngc.TunnelCount(2),
        cost_ceiling=ngc.CostStretch(1.5),
        min_unreserved_bandwidth=3,
    )
    assert (ca.tunnel_count, ca.cost_ceiling, ca.min_unreserved_bandwidth) == (
        ngc.TunnelCount(2),
        ngc.CostStretch(1.5),
        3,
    )
    assert (ca.tunnel_count.per_head, ca.tunnel_count.up_to) == (2, False)
    assert ngc.TunnelCount(3, up_to=True).up_to is True
    d = ngc.DistributedTE(ngc.TunnelCount(1))
    assert (d.cost_ceiling, d.min_unreserved_bandwidth) == (None, 0)
    assert ngc.DistributedTE(
        ngc.TunnelCount(1), cost_ceiling=ngc.AbsoluteCost(7)
    ).cost_ceiling == ngc.AbsoluteCost(7)
    with pytest.raises(TypeError):
        ngc.DistributedTE()  # the tunnel count is required
    with pytest.raises(TypeError):
        ngc.DistributedTE(2)  # a TunnelCount, not an int
    with pytest.raises(TypeError):
        ngc.DistributedTE(
            ngc.TunnelCount(1), ngc.AbsoluteCost(1)
        )  # keyword-only after tunnels
    with pytest.raises(TypeError):
        ngc.TunnelCount(2, True)  # up_to is keyword-only
    for bad, error, message in (
        (
            2**70,
            ValueError,
            "TunnelCount per_head = 1180591620717411303424 is out of range",
        ),
        (True, TypeError, "TunnelCount per_head must be an integer (got bool)"),
        (1.5, TypeError, "TunnelCount per_head must be an integer (got float)"),
    ):
        with pytest.raises(
            error, match=message.replace("(", r"\(").replace(")", r"\)")
        ):
            ngc.TunnelCount(bad)
    with pytest.raises(
        ValueError,
        match=exact(
            "min_unreserved_bandwidth = 1208925819614629174706176 is out of range "
            "[-9223372036854775808, 9223372036854775807]"
        ),
    ):
        ngc.DistributedTE(ngc.TunnelCount(1), min_unreserved_bandwidth=2**80)
    with pytest.raises(TypeError):
        ngc.Pinned([[0, 1]])
    cfg = ngc.FlowPolicyConfig(routing=ca)
    assert isinstance(cfg.routing, ngc.DistributedTE)
    assert cfg.routing.tunnel_count.per_head == 2
    cfg = ngc.FlowPolicyConfig(
        routing=ngc.ShortestPaths(
            allocation=ngc.RouteAllocation.EQUAL,
            source_allocation=ngc.SourceAllocation.OPTIMIZED,
        ),
        admission=ngc.Admission.LOSSY,
    )  # a plain value: the combination is checked where it is used
    assert cfg.routing.allocation == ngc.RouteAllocation.EQUAL
    assert cfg.admission == ngc.Admission.LOSSY
    assert cfg.routing.source_allocation == ngc.SourceAllocation.OPTIMIZED
    # The ineffective tie rule is checked where the model is used, not at construction.
    tie = ngc.ShortestPaths(tie_break=HIGHER_AVAIL)
    assert (tie.tie_break, tie.multipath, tie.parallel_edges) == (
        HIGHER_AVAIL,
        True,
        True,
    )
    o = ngc.OutputOptions()
    names = (
        "edge_flows",
        "residual",
        "saturated_boundary_edges",
        "edge_reservations",
        "min_cut",
    )
    assert tuple(getattr(o, n) for n in names) == (False,) * 5
    for k in range(len(names)):
        flags = tuple(i == k for i in range(len(names)))
        o = ngc.OutputOptions(**dict(zip(names, flags, strict=True)))
        assert tuple(getattr(o, n) for n in names) == flags
    # Pinned keeps its program and black holes by reference (one shared, immutable value each).
    g = diamond()
    a, b = dag_of(g, ids(g, [0, 1])), dag_of(g, ids(g, [2, 3]))
    program = ngc.ForwardingProgram.from_routes(
        g, 0, 3, [ngc.RouteDefinition(0, b), ngc.RouteDefinition(0, a)]
    )
    dark = g.eligibility(excluded_edges=ids(g, [0]))
    pinned = ngc.Pinned(program, black_holes=dark)
    assert pinned.program is program and pinned.black_holes is dark
    assert ngc.Pinned(program).black_holes is None
    # A policy's config is the configuration it was built with, ShortestPaths' absent source
    # allocation and weight basis resolved to the allocation's defaults (Pinned has neither).
    pol = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 3, ngc.FlowPolicyConfig(routing=pinned))
    assert pol.config.routing.program is program
    for allocation, source_allocation, weight_basis in (
        (ngc.RouteAllocation.EQUAL, ngc.SourceAllocation.EQUAL, None),
        (
            ngc.RouteAllocation.WEIGHTED,
            ngc.SourceAllocation.EQUAL,
            ngc.WeightBasis.CAPACITY,
        ),
        (ngc.RouteAllocation.OPTIMIZED, ngc.SourceAllocation.OPTIMIZED, None),
    ):
        got = ngc.FlowPolicy(
            ngc.FlowLedger(g),
            0,
            3,
            ngc.FlowPolicyConfig(routing=ngc.ShortestPaths(allocation=allocation)),
        ).config
        assert (got.routing.source_allocation, got.routing.weight_basis) == (
            source_allocation,
            weight_basis,
        )
    got = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 3, GLOBAL).config
    assert got.routing.source_allocation == ngc.SourceAllocation.OPTIMIZED
    got = ngc.FlowPolicy(
        ngc.FlowLedger(g),
        0,
        3,
        ngc.FlowPolicyConfig(routing=ngc.DistributedTE(ngc.TunnelCount(1))),
    ).config
    assert (got.routing.tie_break, got.routing.parallel_edges) == (
        ngc.RouteTieBreak.FEWER_HOPS,
        True,
    )
    # The constructors' own checks: an integer field outside its C++ type raises ValueError; a
    # field of the wrong class and an unknown keyword raise TypeError.
    assert ngc.AbsoluteCost(2**63 - 1).cost == 2**63 - 1
    assert ngc.AbsoluteCost(-(2**63)).cost == -(2**63)
    for bad in (2**63, -(2**63) - 1):
        with pytest.raises(
            ValueError,
            match=exact(
                f"AbsoluteCost cost = {bad} is out of range [-9223372036854775808, 9223372036854775807]"
            ),
        ):
            ngc.AbsoluteCost(bad)
    with pytest.raises(
        TypeError, match=exact("cost_ceiling must be AbsoluteCost, CostStretch or None")
    ):
        ngc.DistributedTE(ngc.TunnelCount(1), cost_ceiling=3)
    with pytest.raises(TypeError):
        ngc.DistributedTE(
            ngc.TunnelCount(1), min_flow_bandwidth=2
        )  # T0's field is gone


class _ConversionFailed(Exception):
    pass


class _FailingFloat(np.float32):
    # A NumPy floating scalar whose own conversion to a float fails.
    def __float__(self):
        raise _ConversionFailed("own float")


class _FailingIndex:
    # An object whose own conversion to an index fails.
    def __index__(self):
        raise _ConversionFailed("own index")


class _OverflowingFloat(np.float32):
    # A NumPy floating scalar whose own __float__ raises OverflowError.
    def __float__(self):
        raise OverflowError("own conversion")


class _OverflowingIndex:
    # An integer-like object whose own __index__ raises OverflowError.
    def __index__(self):
        raise OverflowError("own conversion")


def test_validation_messages_before_binding():
    g = diamond()
    fg = ngc.FlowLedger(g)

    def rejects(cfg, msg, src=0):
        with pytest.raises(ValueError, match=exact(msg)):
            ngc.FlowPolicy(fg, src, 3, cfg)

    stretch = "cost_ceiling CostStretch factor must be >= 1 (got {})"
    rejects(
        config(ngc.ShortestPaths(), ngc.RouteAllocation.OPTIMIZED, ngc.Admission.LOSSY),
        LOSSY_KERNEL,
    )
    rejects(config(GLOBAL.routing, admission=ngc.Admission.LOSSY), CTE_LOSSY)
    rejects(t1_config(0), "tunnel count must be in [1, 10000]")
    rejects(t1_config(10001), "tunnel count must be in [1, 10000]")
    rejects(
        ngc.FlowPolicyConfig(routing=ngc.DistributedTE(ngc.TunnelCount(2, up_to=True))),
        "DistributedTE requires an exact tunnel count",
    )
    rejects(t1_config(1, cost_ceiling=ngc.CostStretch(0.5)), stretch.format("1/2"))
    for bad in (-1, 2**62):
        rejects(
            t1_config(1, min_unreserved_bandwidth=bad),
            "min_unreserved_bandwidth must be an integer in [0, 2^62)",
        )
    b = dag_of(g, ids(g, [0, 1]))
    program = ngc.ForwardingProgram.from_routes(g, 0, 3, [ngc.RouteDefinition(0, b)])
    rejects(
        ngc.FlowPolicyConfig(routing=ngc.Pinned(program)),
        "program endpoints do not match policy endpoints",
        src=[0, 1],
    )
    rejects(
        SHORTEST_PATHS,
        "node 3 is in both endpoint sets: sources and destinations must not overlap",
        src=[0, 3],
    )
    rejects(SHORTEST_PATHS, "sources: the node set is empty", src=[])
    rejects(SHORTEST_PATHS, f"sources: node {L62} out of node range [0, 4)", src=[L62])
    with pytest.raises(
        ValueError, match=exact("foreign Eligibility: built for another graph")
    ):
        ngc.FlowPolicy(
            fg, 0, 3, SHORTEST_PATHS, diamond().eligibility(excluded_edges=[0])
        )
    with pytest.raises(
        ValueError, match=exact("foreign ShortestPathCache: bound to another graph")
    ):
        ngc.FlowPolicy(fg, 0, 3, SHORTEST_PATHS, None, ngc.ShortestPathCache(diamond()))
    # Arguments of the wrong kind are TypeError: a ledger of None, non-integer node ids, a ledger,
    # config, eligibility or cache of the wrong class and an unknown keyword.
    with pytest.raises(TypeError, match=exact("ledger must be a FlowLedger")):
        ngc.FlowPolicy(None, 0, 3, SHORTEST_PATHS)
    with pytest.raises(
        TypeError, match=exact("sources[0] must be an integer (got float)")
    ):
        ngc.FlowPolicy(fg, [0.0], 3, SHORTEST_PATHS)
    with pytest.raises(
        TypeError, match=exact("destinations must be a sequence of integers")
    ):
        ngc.FlowPolicy(fg, 0, 3.0, SHORTEST_PATHS)
    for wrong in (
        lambda: ngc.FlowPolicy(g, 0, 3, SHORTEST_PATHS),
        lambda: ngc.FlowPolicy(fg, 0, 3, ngc.ShortestPaths()),
        lambda: ngc.FlowPolicy(fg, 0, 3, SHORTEST_PATHS, ngc.ShortestPathCache(g)),
        lambda: ngc.FlowPolicy(fg, 0, 3, SHORTEST_PATHS, None, g.eligibility()),
        lambda: ngc.FlowPolicy(fg, 0, 3, SHORTEST_PATHS, cash=None),
    ):
        with pytest.raises(TypeError):
            wrong()
    fg.reset()  # nothing bound
    pol = ngc.FlowPolicy(fg, 0, 3, SHORTEST_PATHS)
    with pytest.raises(
        TypeError, match=exact("outputs must be an OutputOptions (got ShortestPaths)")
    ):
        pol.place(1, ngc.ShortestPaths())
    assert pol.records() == []


LOSSY_KERNEL = "ShortestPaths OPTIMIZED requires admission=LOSSLESS"
CTE_LOSSY = "CentralizedTE without tunnel_count requires admission=LOSSLESS"
P = ngc.RouteAllocation
A = ngc.Admission
ORIGIN = ngc.SourceAllocation
HC = ngc.RouteTieBreak.HIGHER_CAPACITY


def T1():
    return ngc.DistributedTE(ngc.TunnelCount(1))


INEFFECTIVE = "ShortestPaths tie_break requires multipath=False or parallel_edges=False"


@pytest.mark.parametrize(
    ("routing", "admission", "msg"),
    [
        # The rows of the matrix, their checks in order: the source allocation, admission, the
        # ineffective tie rule (a HIGHER_* tie_break over every route and parallel edge),
        # prefer_headroom, the weights, the tunnel count.
        (
            lambda: ngc.CentralizedTE(
                source_allocation=ORIGIN.EQUAL, prefer_headroom=True
            ),
            A.LOSSY,
            CTE_LOSSY,
        ),
        (
            lambda: ngc.ShortestPaths(
                allocation=P.EQUAL, source_allocation=ORIGIN.OPTIMIZED, tie_break=HC
            ),
            A.LOSSY,
            "ShortestPaths with EQUAL or WEIGHTED requires source_allocation=EQUAL",
        ),
        (
            lambda: ngc.ShortestPaths(allocation=P.OPTIMIZED, tie_break=HC),
            A.LOSSY,
            LOSSY_KERNEL,
        ),
        (
            lambda: ngc.ShortestPaths(
                allocation=P.OPTIMIZED, source_allocation=ORIGIN.EQUAL, tie_break=HC
            ),
            A.LOSSLESS,
            INEFFECTIVE,
        ),
        (lambda: ngc.ShortestPaths(tie_break=HC), A.LOSSY, INEFFECTIVE),
        (
            lambda: ngc.ShortestPaths(
                tie_break=HC, weight_basis=ngc.WeightBasis.CAPACITY
            ),
            A.LOSSLESS,
            INEFFECTIVE,
        ),
        (
            lambda: ngc.ShortestPaths(edge_weights=[1, 1, 1, 1, 1]),
            A.LOSSY,
            "edge_weights require weight_basis=EXPLICIT",
        ),
        (
            lambda: ngc.DistributedTE(
                ngc.TunnelCount(1, up_to=True), cost_ceiling=ngc.CostStretch(0.5)
            ),
            A.LOSSY,
            "DistributedTE requires an exact tunnel count",
        ),
    ],
)
def test_several_defects_report_the_first_in_the_contract_order(
    routing, admission, msg
):
    # Several defects at once: one exact ValueError text, the first in the matrix's order, at
    # construction and at the three evaluators (evaluate_pairs's precheck included); nothing is
    # bound.
    g = diamond()
    fg = ngc.FlowLedger(g)
    cfg = ngc.FlowPolicyConfig(routing=routing(), admission=admission)
    with pytest.raises(ValueError) as e:
        ngc.FlowPolicy(fg, 0, 3, cfg)
    assert str(e.value) == msg
    with pytest.raises(ValueError) as e:
        ngc.evaluate_demand(g, 0, 3, ngc.EvaluationOptions(cfg))
    assert str(e.value) == msg
    with pytest.raises(ValueError) as e:
        ngc.evaluate_pairs(g, [(0, 3)], ngc.EvaluationOptions(cfg))
    assert str(e.value) == msg
    with pytest.raises(ValueError) as e:
        ngc.evaluate_edge_failures(g, 0, 3, ngc.EvaluationOptions(cfg), 1)
    assert str(e.value) == msg
    fg.reset()  # nothing bound


def test_enum_values_outside_their_enumeration_are_rejected():
    # An enum built from an integer outside the enumeration is a domain error named by its field:
    # FlowPolicyConfig and the models are plain values, checked wherever they are used, first
    # among their checks, so it never reaches a policy or record; the searches check their own.
    bad_tie = ngc.RouteTieBreak(99)
    for options, search in (
        (ngc.SpfOptions(ngc.Direction.FORWARD, tie_break=bad_tie), ngc.spf),
        (ngc.KspOptions(1, tie_break=bad_tie), lambda g, s, o: ngc.ksp(g, s, 3, o)),
    ):
        with pytest.raises(
            ValueError, match=exact("tie_break: 99 is not a value of RouteTieBreak")
        ):
            search(diamond(), 0, options)
    g = diamond()
    bad_source = "source_allocation: 99 is not a value of SourceAllocation"
    cases = [
        (
            ngc.ShortestPaths(allocation=ngc.RouteAllocation(99)),
            "allocation: 99 is not a value of RouteAllocation",
        ),
        (ngc.ShortestPaths(source_allocation=ORIGIN(99)), bad_source),
        (
            ngc.ShortestPaths(allocation=P.WEIGHTED, weight_basis=ngc.WeightBasis(99)),
            "weight_basis: 99 is not a value of WeightBasis",
        ),
        (ngc.CentralizedTE(source_allocation=ORIGIN(99)), bad_source),
    ]
    for routing in (ngc.ShortestPaths(), T1(), GLOBAL.routing):
        cases.append(
            (
                ngc.FlowPolicyConfig(routing=routing, admission=ngc.Admission(99)),
                "admission: 99 is not a value of Admission",
            )
        )
    for cfg, message in cases:
        if not isinstance(cfg, ngc.FlowPolicyConfig):
            cfg = ngc.FlowPolicyConfig(routing=cfg)
        with pytest.raises(ValueError, match=exact(message)):
            ngc.FlowPolicy(ngc.FlowLedger(g), 0, 3, cfg)
        with pytest.raises(ValueError, match=exact(message)):
            ngc.evaluate_demand(g, 0, 3, ngc.EvaluationOptions(cfg))
    for routing in (
        ngc.ShortestPaths(tie_break=bad_tie),
        ngc.DistributedTE(ngc.TunnelCount(1), tie_break=bad_tie),
    ):
        cfg = ngc.FlowPolicyConfig(routing=routing, admission=ngc.Admission(98))
        message = (
            "admission: 98 is not a value of Admission"  # the config's own enums first
        )
        with pytest.raises(ValueError, match=exact(message)):
            ngc.FlowPolicy(ngc.FlowLedger(g), 0, 3, cfg)
        cfg = ngc.FlowPolicyConfig(routing=routing)
        with pytest.raises(
            ValueError, match=exact("tie_break: 99 is not a value of RouteTieBreak")
        ):
            ngc.FlowPolicy(ngc.FlowLedger(g), 0, 3, cfg)
    assert len(ngc.RouteAllocation.__members__) == 3  # the loop below is not empty
    for valid in ngc.RouteAllocation.__members__.values():
        assert config(ngc.ShortestPaths(), valid).routing.allocation == valid


def test_placement_result_fields_and_outputs():
    g = diamond()
    fg = ngc.FlowLedger(g)
    pol = ngc.FlowPolicy(fg, [0], [3], SHORTEST_PATHS)
    r = pol.place(4)
    assert (r.call_delivered, r.call_undelivered, r.termination) == (4, 0, T.COMPLETE)
    assert not hasattr(r, "local_outcomes")
    # The epoch's first call installs the program: its routes and weights.
    assert (r.routes_changed, r.weights_changed, r.call_edge_drops) == (True, True, [])
    assert (r.delivery_delta, r.call_source_drops, r.edge_drops, r.source_drops) == (
        4,
        0,
        [],
        0,
    )
    assert (r.requested, r.offered, r.admitted, r.delivered, r.unadmitted) == (
        4,
        None,
        4,
        4,
        0,
    )
    row = r.by_source[0]
    assert isinstance(r.by_source, tuple) and len(r.by_source) == 1
    assert (row.source, row.offered, row.admitted, row.delivered) == (0, None, 4, 4)
    assert (row.source_drops, row.total_edge_drops) == (0, 0)
    with pytest.raises(TypeError):
        ngc.SourceFlow()
    assert_read_only(row, "delivered")
    assert r.marginal_cost_curve is None and r.total_cost is None
    assert (
        r.edge_flows is None
        and r.residual is None
        and r.saturated_boundary_edges is None
        and r.min_cut is None
    )
    r = pol.place(INF, ALL_OUT)
    assert (r.call_delivered, r.call_undelivered, r.termination) == (
        3,
        INF,
        T.EXHAUSTED,
    )
    assert r.edge_flows.dtype == np.int64 and r.residual.dtype == np.int64
    assert list(r.residual) == list(fg.residual())
    own = [0] * 4
    for rec in pol.records():
        for e, x in rec.edge_flows:
            own[e] += x
    assert list(r.edge_flows) == own
    assert r.saturated_boundary_edges == fg.saturated_boundary_edges(0)
    # By hand: the flow's 7 units split equally at node 0, the leftover unit to the lower
    # EdgeId: 4 on 0 -> 1 -> 3 and 3 on 0 -> 2 -> 3, the most that keeps 0 -> 2 within its 3.
    # Only 0 -> 2 is saturated where it leaves what the source still reaches.
    assert r.edge_flows.tolist() == to_ids(g, [4, 4, 3, 3])
    assert r.residual.tolist() == to_ids(g, [1, 1, 0, 0])
    assert r.saturated_boundary_edges == ids(g, [2])
    assert r.min_cut is None  # CentralizedTE only
    assert pol.delivered() == 7
    with pytest.raises(
        ValueError,
        match=exact("offer must be an integer in [0, 2^62) or permitted +inf"),
    ):
        pol.place(1.5)
    with pytest.raises(_ConversionFailed):  # a float request whose own conversion fails
        pol.place(_FailingFloat(1.0))
    for bad in (_OverflowingFloat(1.0), _OverflowingIndex()):  # its OverflowError too
        with pytest.raises(OverflowError, match=exact("own conversion")):
            pol.place(bad)
    with pytest.raises(
        ValueError,
        match=exact("finite placement after infinite placement requires a new epoch"),
    ):
        pol.place(np.int64(1))
    pol.remove()
    assert pol.records() == [] and list(fg.carried()) == [0, 0, 0, 0]
    assert pol.place(np.int64(2)).call_delivered == 2
    # CentralizedTE: the marginal curve; min_cut with bottleneck.
    fg2 = ngc.FlowLedger(diamond(5, 3, 2))
    gl = ngc.FlowPolicy(fg2, 0, 3, GLOBAL)
    r = gl.place(INF, ngc.OutputOptions(saturated_boundary_edges=True))
    assert r.marginal_cost_curve == [(2, 5), (4, 3)] and r.min_cut is None
    gl.remove()
    r = gl.place(INF, ngc.OutputOptions(min_cut=True))
    assert r.min_cut == sorted(ids(diamond(5, 3, 2), [0, 2]))
    assert [rec.delivered for rec in gl.records()] == [8]


def test_lifetimes_keep_alive_and_orphans():
    g = diamond()
    fg = ngc.FlowLedger(g)
    el = g.eligibility(excluded_edges=ids(g, [2]))
    cache = ngc.ShortestPathCache(g)
    pol = ngc.FlowPolicy(fg, 0, 3, SHORTEST_PATHS, el, cache)
    held = [weakref.ref(x) for x in (fg, el, cache)]
    del fg, el, cache
    gc.collect()
    assert all(r() is not None for r in held)
    assert (
        pol.place(10).call_delivered == 5
    )  # the ledger, the Eligibility and the cache are alive
    del pol
    gc.collect()
    assert all(r() is None for r in held)  # the policy alone held them
    fg2 = ngc.FlowLedger(g)
    p = ngc.FlowPolicy(fg2, 0, 3, SHORTEST_PATHS)
    p.place(2)
    del p
    gc.collect()
    if HOOKS:
        assert fg2._state()["orphans"] == [0]
    assert fg2.records()[0].delivered == 2
    with pytest.raises(ValueError, match="bound"):
        q = ngc.FlowPolicy(fg2, 0, 3, SHORTEST_PATHS)
        fg2.reset()
    del q
    gc.collect()
    fg2.reset()
    assert fg2.records() == []


@needs_hooks
def test_owner_hook_is_the_owner_of_the_records():
    # The private hook FlowPolicy._owner reads the ledger-local owner id: one per policy of a
    # ledger, known before the policy has a record, and carried by its records as FlowRecord.policy_id.
    g = diamond()
    fg = ngc.FlowLedger(g)
    policies = [ngc.FlowPolicy(fg, 0, 3, SHORTEST_PATHS) for _ in range(3)]
    owners = [p._owner() for p in policies]
    assert len(set(owners)) == 3
    for p in reversed(policies):
        assert p.place(1).call_delivered == 1
    for p, owner in zip(policies, owners, strict=True):
        assert [r.policy_id for r in p.records()] == [owner]


def pinned_over(g, *dags):
    """Pinned over a program of the split routes `dags` of head 0 to node 3."""
    program = ngc.ForwardingProgram.from_routes(
        g, 0, 3, [ngc.RouteDefinition(0, d) for d in dags]
    )
    return ngc.FlowPolicyConfig(routing=ngc.Pinned(program))


def test_pinned_records():
    g = diamond()
    fg = ngc.FlowLedger(g)
    r1, r2 = dag_of(g, ids(g, [0, 1])), dag_of(g, ids(g, [2, 3]))
    el = g.eligibility(excluded_edges=ids(g, [3]))
    # The replay fails the second route (its last edge excluded): the head's share goes to the
    # first; one record per program entry, the failed one with an empty route.
    pol = ngc.FlowPolicy(fg, 0, 3, pinned_over(g, r1, r2), el)
    r = pol.place(9)
    assert (r.call_delivered, r.termination) == (5, T.EXHAUSTED)
    rec, failed = pol.records()
    assert rec.route == r1 and pol.record(rec.flow_id).delivered == 5
    assert (rec.entry_id, failed.entry_id, failed.min_route_cost) == (0, 1, None)
    # record() reads the policy's own records only, not another policy's on the same ledger.
    other = ngc.FlowPolicy(fg, 0, 3, pinned_over(g, r1))
    other.place(1)
    theirs = other.records()[0].flow_id
    with pytest.raises(
        ValueError, match=exact(f"FlowId {theirs} is not a record of this policy")
    ):
        pol.record(theirs)
    assert fg.record(theirs).delivered == 0
    elsewhere = ngc.FlowLedger(g).uid()
    with pytest.raises(
        ValueError,
        match=exact(
            f"FlowId ({elsewhere}, 0) belongs to another ledger (uid {fg.uid()})"
        ),
    ):
        pol.record((elsewhere, 0))
    pol.remove()
    with pytest.raises(
        ValueError, match=exact(f"FlowId {rec.flow_id} is retired or unknown")
    ):
        pol.record(rec.flow_id)
    all_down = ngc.FlowPolicy(
        fg, 0, 3, pinned_over(g, r2), g.eligibility(excluded_edges=ids(g, [3]))
    )
    r = all_down.place(1)
    assert (r.call_delivered, r.termination, r.by_source[0].admitted) == (
        0,
        T.NO_ROUTE,
        0,
    )


def test_dag_cache_statistics_and_projection():
    g = diamond()
    cache = ngc.ShortestPathCache(g)
    fg = ngc.FlowLedger(g)
    a = ngc.FlowPolicy(fg, 0, 3, SHORTEST_PATHS, None, cache)
    b = ngc.FlowPolicy(
        fg,
        0,
        3,
        config(
            ngc.ShortestPaths(
                tie_break=HIGHER_AVAIL, multipath=False, parallel_edges=False
            )
        ),
        None,
        cache,
    )
    s = cache.stats()
    assert (s.hits, s.misses, s.entries, s.live_values) == (1, 1, 1, 1)
    assert s.bytes > 0
    assert a.place(2).call_delivered == 2 and b.place(2).call_delivered == 2
    cache.clear()
    s = cache.stats()
    assert (s.hits, s.misses, s.entries, s.bytes) == (1, 1, 0, 0)
    assert s.live_values == 1  # held by the policies
    assert cache.graph_uid() == g.uid()


def test_capacity_batch_and_workers(monkeypatch):
    g = diamond(5, 3)
    opts = ngc.EvaluationOptions(GLOBAL, outputs=ngc.OutputOptions(edge_flows=True))
    r = ngc.evaluate_demand(g, 0, 3, opts)
    assert (r.call_delivered, r.termination, r.marginal_cost_curve) == (
        8,
        T.EXHAUSTED,
        [(2, 8)],
    )
    assert int(np.sum(r.edge_flows)) == 16
    common = g.eligibility(excluded_edges=ids(g, [3]))
    pairs = [(0, 3), ([0], [3]), (1, 3), (0, [1, 2])]
    per = [None, g.eligibility(excluded_edges=ids(g, [0])), None, None]
    want = [5, 0, 5, 8]
    for workers in (1, 3, None):
        got = ngc.evaluate_pairs(
            g, pairs, ngc.EvaluationOptions(GLOBAL, eligibility=common), per, workers
        )
        assert [x.call_delivered for x in got] == want
    # workers=None reads NGRAPH_CORE_BATCH_THREADS. The count changes only how many threads run,
    # never a result, so this checks only that the default path with the variable set gives the
    # same results (its parsing: the C++ test Evaluators.DefaultWorkersParseTheEnvironment).
    monkeypatch.setenv("NGRAPH_CORE_BATCH_THREADS", "4")
    got = ngc.evaluate_pairs(
        g, pairs, ngc.EvaluationOptions(GLOBAL, eligibility=common), per
    )
    assert [x.call_delivered for x in got] == want
    assert ngc.evaluate_pairs(g, [], ngc.EvaluationOptions(GLOBAL)) == []
    with pytest.raises(ValueError, match="node 7"):
        ngc.evaluate_pairs(
            g, [(0, 3), (0, 7), (0, 8)], ngc.EvaluationOptions(GLOBAL), workers=4
        )
    with pytest.raises(
        ValueError, match=exact("per_pair_eligibility has 1 entries for 4 pairs")
    ):
        ngc.evaluate_pairs(g, pairs, ngc.EvaluationOptions(GLOBAL), [None])
    with pytest.raises(
        ValueError, match=exact("evaluate_pairs does not accept Pinned")
    ):
        ngc.evaluate_pairs(
            g,
            pairs,
            ngc.EvaluationOptions(pinned_over(g, dag_of(g, ids(g, [0, 1])))),
        )
    cache = ngc.ShortestPathCache(g)
    eb = ngc.EvaluationOptions(SHORTEST_PATHS, cache=cache)
    with pytest.raises(ValueError, match="workers == 1"):
        ngc.evaluate_pairs(g, pairs, eb, workers=2)
    got = ngc.evaluate_pairs(
        g, pairs, eb
    )  # None with a cache: 1 worker, even with the env set
    assert [x.call_delivered for x in got] == [
        7,
        7,
        5,
        7,
    ]  # EB: [4, 3] over links of 5 and 3
    with pytest.raises(ValueError, match="workers"):
        ngc.evaluate_pairs(g, pairs, ngc.EvaluationOptions(GLOBAL), workers=0)


def test_edge_failures_candidates_and_signed_reductions(monkeypatch):
    # CentralizedTE: an unsaturated critical edge appears.
    g = build(3, [(0, 1, 9, 1), (1, 2, 4, 1), (0, 2, 2, 5)])
    s = ngc.evaluate_edge_failures(g, 0, 2, ngc.EvaluationOptions(GLOBAL))
    assert s.baseline.call_delivered == 6
    pos = g.input_positions()
    assert {int(pos[e.edge]): e.delivery_reduction for e in s.entries} == {
        0: 4,
        1: 4,
        2: 2,
    }
    # workers=None reads NGRAPH_CORE_SENSITIVITY_THREADS: the same results on 3 workers (only that
    # the default path with the variable set works; its parsing is tested in C++).
    monkeypatch.setenv("NGRAPH_CORE_SENSITIVITY_THREADS", "3")
    assert plain(
        ngc.evaluate_edge_failures(g, 0, 2, ngc.EvaluationOptions(GLOBAL))
    ) == plain(s)
    # A policy model: every allowed edge is a candidate, a zero-capacity one too. Equal splitting
    # over parallel links of 1, 10 and 0 fits 2 units (a third would cross the zero link);
    # without the thin link 1 fits, without the wide link 1, without the zero link 2.
    p = build(2, [(0, 1, 1, 1), (0, 1, 10, 1), (0, 1, 0, 1)])
    opts = ngc.EvaluationOptions(SHORTEST_PATHS)
    s = ngc.evaluate_edge_failures(p, 0, 1, opts, workers=2)
    assert s.baseline.call_delivered == 2
    pos = p.input_positions()
    assert {int(pos[e.edge]): e.delivery_reduction for e in s.entries} == {
        0: 1,
        1: 1,
        2: 0,
    }
    # A negative loss on a zero-capacity link (the search graph S): the only shortest route 0 -> 1 has capacity
    # 0, so equal splitting carries 0; without it the route 0 -> 2 -> 1 carries 5.
    g8 = build(3, [(0, 1, 0, 1), (0, 2, 5, 1), (2, 1, 5, 1)])
    s8 = ngc.evaluate_edge_failures(g8, 0, 1, opts)
    assert s8.baseline.call_delivered == 0
    pos8 = g8.input_positions()
    assert {int(pos8[e.edge]): e.delivery_reduction for e in s8.entries} == {
        0: -5,
        1: 0,
        2: 0,
    }
    with pytest.raises(ValueError, match="both endpoint sets"):
        ngc.evaluate_edge_failures(p, 0, 0, opts)
    # The baseline is evaluate_demand() with the requested outputs.
    for cfg in (GLOBAL, SHORTEST_PATHS):
        s = ngc.evaluate_edge_failures(
            g, 0, 2, ngc.EvaluationOptions(cfg, outputs=ALL_OUT)
        )
        want = ngc.evaluate_demand(g, 0, 2, ngc.EvaluationOptions(cfg, outputs=ALL_OUT))
        assert want.edge_flows is not None and want.residual is not None
        assert plain(s.baseline) == plain(want)


def test_evaluator_options_and_results():
    # The options hold what they were given (None: everything allowed, no supplied cache, no
    # output); the results have no public constructor.
    g = diamond()
    el = g.eligibility(excluded_edges=ids(g, [0]))
    cache = ngc.ShortestPathCache(g)
    for cls in (ngc.EvaluationOptions, ngc.EvaluationOptions):
        o = cls(
            GLOBAL,
            eligibility=el,
            cache=cache,
            outputs=ngc.OutputOptions(residual=True),
        )
        assert (
            o.config.routing.source_allocation == GLOBAL.routing.source_allocation
            and o.eligibility is el
            and o.cache is cache
        )
        out = o.outputs
        assert (out.edge_flows, out.residual, out.saturated_boundary_edges) == (
            False,
            True,
            False,
        )
        d = cls(SHORTEST_PATHS)
        assert d.eligibility is None and d.cache is None
        assert (
            d.outputs.edge_flows,
            d.outputs.residual,
            d.outputs.saturated_boundary_edges,
        ) == (
            False,
            False,
            False,
        )
    for cls in (ngc.SensitivityEntry, ngc.SensitivityResult):
        with pytest.raises(TypeError, match="No constructor defined"):
            cls()


def test_explicit_workers_win_over_the_environment(monkeypatch):
    # The environment supplies only the default of workers=None, which is 1 whenever a cache is
    # supplied; an explicit count is always honoured (a building cache accepts only 1).
    monkeypatch.setenv("NGRAPH_CORE_BATCH_THREADS", "4")
    monkeypatch.setenv("NGRAPH_CORE_SENSITIVITY_THREADS", "4")
    g = diamond(5, 3)
    building = ngc.ShortestPathCache(g)
    bopts = ngc.EvaluationOptions(SHORTEST_PATHS, cache=building)
    sopts = ngc.EvaluationOptions(SHORTEST_PATHS, cache=building)
    want = [
        r.call_delivered
        for r in ngc.evaluate_pairs(g, [(0, 3), (1, 3)], bopts, workers=1)
    ]
    assert want == [7, 5]
    assert [
        r.call_delivered for r in ngc.evaluate_pairs(g, [(0, 3), (1, 3)], bopts)
    ] == want
    base = plain(ngc.evaluate_edge_failures(g, 0, 3, sopts, workers=1))
    assert plain(ngc.evaluate_edge_failures(g, 0, 3, sopts)) == base
    for call in (
        lambda: ngc.evaluate_pairs(g, [(0, 3)], bopts, workers=2),
        lambda: ngc.evaluate_edge_failures(g, 0, 3, sopts, workers=2),
    ):
        with pytest.raises(ValueError, match="unless workers == 1"):
            call()
    with pytest.raises(ValueError, match=exact("workers must be >= 1 (got 0)")):
        ngc.evaluate_edge_failures(g, 0, 3, ngc.EvaluationOptions(GLOBAL), workers=0)


def test_python_keeps_what_a_policy_and_an_option_need_alive():
    # Once the caller drops them, the evaluators' options still apply their Eligibility and cache,
    # the ledger still reads its graph, and a policy still places on its ledger. Graphs of the
    # same shape are built in between, so the memory of a released value would be reused.
    g = diamond()
    want = [
        g.edge_capacity_units().tolist(),
        g.edge_costs().tolist(),
        g.edge_destinations().tolist(),
    ]
    ledger = ngc.FlowLedger(g)
    policy = ngc.FlowPolicy(ledger, 0, 3, SHORTEST_PATHS)
    opts = ngc.EvaluationOptions(
        SHORTEST_PATHS,
        eligibility=g.eligibility(excluded_edges=ids(g, [2])),
        cache=ngc.ShortestPathCache(g),
    )
    gc.collect()
    reused = [diamond(7, 9, 2) for _ in range(64)]
    assert (
        ngc.evaluate_demand(g, 0, 3, opts).call_delivered == 5
    )  # without the edge 0 -> 2
    assert opts.cache.stats().misses == 1  # the call's lookup
    del g, opts
    gc.collect()
    reused += [diamond(7, 9, 2) for _ in range(64)]
    kept = ledger.graph()
    assert [
        kept.edge_capacity_units().tolist(),
        kept.edge_costs().tolist(),
        kept.edge_destinations().tolist(),
    ] == want
    del kept, ledger
    gc.collect()
    reused += [diamond(7, 9, 2) for _ in range(64)]
    assert policy.place(10).call_delivered == 7  # Equal over routes of 5 and 3: [4, 3]
    assert len(reused) == 192


def test_exception_classes_of_conversions():
    # An integer beyond the float range is an accepted type with a value outside the domain, a
    # ValueError, both where a float is expected (the resolution) and as a request (integer
    # units). A request's own conversion failure reaches the caller unchanged.
    for call in (
        lambda: ngc.Graph(2, [0], [1], [5], [1], 10**400),
        lambda: ngc.to_units(1.0, 10**400),
        lambda: ngc.as_float(1, 10**400),
    ):
        with pytest.raises(
            ValueError,
            match=exact("resolution must be a power of two in [2^-32, 2^32]"),
        ) as error:
            call()
        assert error.type is ValueError
    g = build(2, [(0, 1, 5, 1)])
    p = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 1, SHORTEST_PATHS)
    for big in (10**400, -(10**400)):
        with pytest.raises(
            ValueError,
            match=exact("offer must be an integer in [0, 2^62) or permitted +inf"),
        ):
            p.place(big)
    # An index's own failure (a float's: test_placement_result_fields_and_outputs).
    with pytest.raises(_ConversionFailed):
        p.place(_FailingIndex())
    assert p.delivered() == 0


@needs_hooks
@pytest.mark.parametrize("routing", ["shortest_paths", "global", "weighted"])
@pytest.mark.parametrize("traced", [False, True])
def test_python_publication_is_all_or_nothing(routing, traced):
    # Every Python object of a placement exists before anything is published, so the k-th
    # Python allocation failing raises with the ledger and the policy unchanged (an earlier
    # call's flow in the same epoch included), never crashes, and a retry places as a call
    # without failure does.
    cfg = GLOBAL if routing == "global" else SHORTEST_PATHS
    if routing == "weighted":
        cfg = config(
            ngc.CentralizedTE(
                source_allocation=ngc.SourceAllocation.WEIGHTED,
                source_weights=[1001, 0],
            )
        )
    g = build(2, [(0, 1, 10, 1)])

    for entry in (
        ("_place_traced",) if traced else ("place", "replace", "forward", "evaluate")
    ):
        current_probe = []

        def make(entry=entry, current_probe=current_probe):
            ledger = ngc.FlowLedger(g)
            policy = ngc.FlowPolicy(ledger, 0, 1, cfg)
            assert policy.place(1).call_delivered == 1
            probe = ext._ConcurrencyProbe()
            current_probe[:] = [probe]

            def snapshot():
                return [
                    full_ledger_state(ledger),
                    repr(policy._state()),
                    probe.counts()["native_begin"],
                ]

            return getattr(policy, entry), (3,), snapshot

        def allocation_call(fn, args, k, current_probe=current_probe):
            probe = current_probe[0]
            result = probe.run(ext._fail_python_allocation, fn, args, k)
            if not result[0]:
                assert probe.counts()["native_begin"] == 1
            return result

        # Each scalar entry owns its result allocation path. Not vacuous: the result
        # object at least; traced, also the pair and trace. Retry must return every field.
        assert sweep_python_allocation(make, allocation_call=allocation_call) >= (
            3 if traced else 1
        ), entry


@needs_hooks
def test_evaluator_results_survive_python_allocation_failures():
    # Fixed evaluator result objects exist before native entry, even without a supplied cache.
    # A failed Python allocation raises before native work; the completed call matches a clean call.
    g = diamond(5, 3)
    calls = [
        (ngc.evaluate_demand, (g, 0, 3, ngc.EvaluationOptions(GLOBAL))),
        (ngc.evaluate_pairs, (g, [(0, 3), (1, 3)], ngc.EvaluationOptions(GLOBAL))),
        # Node sets given as lists: the call makes no temporary list, so the result list itself
        # is a new allocation the sweep reaches.
        (
            ngc.evaluate_pairs,
            (g, [([0], [3]), ([1], [3])], ngc.EvaluationOptions(GLOBAL)),
        ),
        (ngc.evaluate_edge_failures, (g, 0, 3, ngc.EvaluationOptions(GLOBAL))),
    ]
    for fn, args in calls:
        current_probe = []

        def make(fn=fn, args=args, current_probe=current_probe):
            probe = ext._ConcurrencyProbe()
            current_probe[:] = [probe]
            return fn, args, lambda: probe.counts()["native_begin"]

        def allocation_call(fn, args, k, current_probe=current_probe):
            probe = current_probe[0]
            result = probe.run(ext._fail_python_allocation, fn, args, k)
            if not result[0]:
                assert probe.counts()["native_begin"] == 1
            return result

        failures = sweep_python_allocation(make, allocation_call=allocation_call)
        # Not vacuous: the result objects allocate (their fields convert when read).
        assert failures >= 1, fn.__name__


def test_batch_reports_the_lowest_failing_pair():
    # The arguments first, then each pair completely in pair order, across the
    # Python/C++ boundary.
    g = build(2, [(0, 1, 3, 7)])
    opts = ngc.EvaluationOptions(GLOBAL)
    overlap = (
        "node 0 is in both endpoint sets: sources and destinations must not overlap"
    )
    out_of_range = "sources: node 2 out of node range [0, 2)"
    with pytest.raises(ValueError, match=exact(overlap)):
        ngc.evaluate_pairs(g, [(0, 0), (2, 1)], opts, workers=2)
    with pytest.raises(ValueError, match=exact(out_of_range)):
        ngc.evaluate_pairs(g, [(0, 1), (2, 1), (0, 0)], opts)
    other = build(2, [(0, 1, 3, 7)]).eligibility()
    with pytest.raises(
        ValueError,
        match=exact(
            "foreign Eligibility: per_pair_eligibility[0] built for another graph"
        ),
    ):
        ngc.evaluate_pairs(
            g, [(0, 1), (2, 1)], opts, per_pair_eligibility=[other, None]
        )
    with pytest.raises(
        TypeError,
        match=exact("per_pair_eligibility entries must be Eligibility or None"),
    ):
        ngc.evaluate_pairs(g, [(0, 1), (2, 1)], opts, per_pair_eligibility=[3, None])
    with pytest.raises(ValueError, match=exact(out_of_range)):
        ngc.evaluate_pairs(g, [(0, 1), (2, 1)], opts, per_pair_eligibility=[None, 3])
    # An argument error comes before every pair error.
    with pytest.raises(ValueError, match=exact("workers must be >= 1 (got 0)")):
        ngc.evaluate_pairs(g, [(0, 0)], opts, workers=0)


def test_batch_checks_an_explicit_empty_eligibility_list():
    # A given per_pair_eligibility has one entry per pair, also when it is empty.
    g = build(2, [(0, 1, 3, 7)])
    opts = ngc.EvaluationOptions(GLOBAL)
    with pytest.raises(
        ValueError, match=exact("per_pair_eligibility has 0 entries for 1 pairs")
    ):
        ngc.evaluate_pairs(g, [(0, 1)], opts, per_pair_eligibility=[])
    assert ngc.evaluate_pairs(g, [], opts, per_pair_eligibility=[]) == []
    assert [
        r.call_delivered
        for r in ngc.evaluate_pairs(g, [(0, 1)], opts, per_pair_eligibility=[None])
    ] == [3]


def test_global_zero_call_reports_the_cumulative_curve():
    # CentralizedTE's curve is cumulative (its amounts sum to delivered()); a zero call
    # reports it (None before the epoch's preparation and after remove()), changes nothing and
    # claims no cut, even when asked for one.
    g = build(2, [(0, 1, 3, 7)])
    ledger = ngc.FlowLedger(g)
    p = ngc.FlowPolicy(ledger, 0, 1, GLOBAL)
    assert p.place(0).marginal_cost_curve is None
    assert p.place(2).marginal_cost_curve == [(7, 2)]
    z = p.place(0, ALL_OUT)
    assert (z.call_delivered, z.marginal_cost_curve, z.min_cut) == (0, [(7, 2)], None)
    assert p.place(0, ngc.OutputOptions(min_cut=True)).min_cut is None  # a zero call
    assert p.delivered() == 2 and [int(x) for x in ledger.carried()] == [2]
    p.remove()
    assert p.place(0).marginal_cost_curve is None
