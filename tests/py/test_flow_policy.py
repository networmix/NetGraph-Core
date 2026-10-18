"""The policy engine, the DagCache and the evaluators through the bindings (docs/semantics.md):
configuration objects and their structural validation, PlacementResult fields and outputs,
lifetimes, the evaluators with workers, and the cache statistics. The fixtures of the reference
(``tests/reference/policy/``) replay in ``test_flow_policy_fixtures.py``."""

from __future__ import annotations

import gc
import math
import sys
import weakref

import numpy as np
import pytest
from helpers import (
    ALL_OUT,
    COSTONLY,
    GLOBAL,
    INF,
    L62,
    RES,
    build,
    config,
    dag_of,
    exact,
    full_ledger_state,
    ids,
    plain,
    sel,
    sweep_python_allocation,
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
        ngc.FlowPolicyConfig(
            ngc.CostOnly(), ngc.FlowPlacement.MAX_FLOW, ngc.Admission.LOSSLESS, sel()
        )
    with pytest.raises(TypeError):
        ngc.FlowPolicyConfig(
            routing=ngc.CostOnly(), placement=ngc.FlowPlacement.MAX_FLOW
        )
    with pytest.raises(TypeError, match="incompatible"):
        ngc.CostOnly(flow_count=2)  # a field outside its routing variant: TypeError
    with pytest.raises(TypeError, match="incompatible"):
        ngc.Global(max_flow_count=4)
    with pytest.raises(TypeError, match="routing"):
        config("CostOnly")
    ca = ngc.CapacityAware(
        flow_count=2,
        max_flow_count=4,
        max_path_cost=ngc.Stretch(1.5),
        min_flow_bandwidth=3,
        reoptimize=True,
        diminishing=ngc.DiminishingReturns(2, 0.25),
    )
    assert (ca.flow_count, ca.max_flow_count, ca.min_flow_bandwidth, ca.reoptimize) == (
        2,
        4,
        3,
        True,
    )
    assert ca.max_path_cost == ngc.Stretch(1.5)
    assert (ca.diminishing.window, ca.diminishing.fraction) == (2, 0.25)
    assert ngc.CapacityAware().max_path_cost is None
    assert ngc.CapacityAware(
        max_path_cost=ngc.Absolute(7)
    ).max_path_cost == ngc.Absolute(7)
    with pytest.raises(TypeError):
        ngc.CapacityAware(2)  # keyword-only
    for bad, message in (
        (2**40, "flow_count = 1099511627776 is out of range [-2147483648, 2147483647]"),
        (True, "flow_count must be an integer (got bool)"),
        (1.5, "flow_count must be an integer (got float)"),
    ):
        with pytest.raises(ValueError, match=exact(message)):
            ngc.CapacityAware(flow_count=bad)
    with pytest.raises(
        ValueError,
        match=exact(
            "min_flow_bandwidth (units) = 1208925819614629174706176 is out of range "
            "[-9223372036854775808, 9223372036854775807]"
        ),
    ):
        ngc.CapacityAware(min_flow_bandwidth=2**80)
    with pytest.raises(TypeError, match="PredDAG"):
        ngc.Pinned([[0, 1]])
    cfg = config(ca)
    assert isinstance(cfg.routing, ngc.CapacityAware)
    assert cfg.routing.flow_count == 2
    # Every field reads back as built, with the documented defaults.
    d = ngc.CapacityAware()
    assert (
        d.flow_count,
        d.max_flow_count,
        d.max_path_cost,
        d.min_flow_bandwidth,
        d.reoptimize,
        d.diminishing,
    ) == (1, None, None, 0, False, None)
    s = sel(ngc.EdgeTieBreak.PREFER_HIGHER_PRIORITY, RES, False, True)
    cfg = config(ngc.Global(), ngc.FlowPlacement.MAX_FLOW, ngc.Admission.LOSSY, s)
    assert cfg.placement == ngc.FlowPlacement.MAX_FLOW
    assert cfg.admission == ngc.Admission.LOSSY
    assert cfg.selection == s
    assert (ngc.OutputOptions().edge_flows, ngc.OutputOptions().residual) == (
        False,
        False,
    )
    assert ngc.OutputOptions().bottleneck is False
    for flags in ((True, False, False), (False, True, False), (False, False, True)):
        o = ngc.OutputOptions(
            edge_flows=flags[0], residual=flags[1], bottleneck=flags[2]
        )
        assert (o.edge_flows, o.residual, o.bottleneck) == flags
    # Pinned keeps the bundles in supply order, copied from the given list: each read is a new
    # list of the same PredDAG objects, and changing a list read changes nothing.
    g = diamond()
    a, b = dag_of(g, ids(g, [0, 1])), dag_of(g, ids(g, [2, 3]))
    given = [b, a]
    pinned = ngc.Pinned(given)
    given.append(a)
    first = pinned.bundles
    assert first == [b, a] and first[0] is b and first[1] is a
    first.append(b)
    second = pinned.bundles
    assert second == [b, a] and second is not first
    # A policy's config is the configuration it was built with.
    pol = ngc.FlowPolicy(
        ngc.FlowLedger(g), 0, 3, config(pinned, ngc.FlowPlacement.MAX_FLOW)
    )
    assert pol.config.placement == ngc.FlowPlacement.MAX_FLOW
    assert pol.config.routing.bundles == [b, a]
    # The constructors' own checks: an integer field outside its C++ type raises ValueError; a
    # field of the wrong class, a reoptimize that is not a truth value and an unknown keyword
    # raise TypeError. reoptimize takes a bool, a number or None by its truth value.
    with pytest.raises(
        ValueError,
        match=exact(
            "max_flow_count = 2147483648 is out of range [-2147483648, 2147483647]"
        ),
    ):
        ngc.CapacityAware(max_flow_count=2**31)
    with pytest.raises(
        ValueError,
        match=exact(
            "diminishing window = 2147483648 is out of range [-2147483648, 2147483647]"
        ),
    ):
        ngc.DiminishingReturns(2**31, 0.5)
    assert ngc.Absolute(2**63 - 1).cost == 2**63 - 1
    assert ngc.Absolute(-(2**63)).cost == -(2**63)
    for bad in (2**63, -(2**63) - 1):
        with pytest.raises(
            ValueError,
            match=exact(
                f"Absolute cost = {bad} is out of range [-9223372036854775808, 9223372036854775807]"
            ),
        ):
            ngc.Absolute(bad)
    with pytest.raises(
        TypeError, match=exact("max_path_cost must be Absolute, Stretch or None")
    ):
        ngc.CapacityAware(max_path_cost=3)
    with pytest.raises(
        TypeError, match=exact("diminishing must be DiminishingReturns or None")
    ):
        ngc.CapacityAware(diminishing=0.5)
    with pytest.raises(TypeError):
        ngc.CapacityAware(reoptimize="yes")
    assert [ngc.CapacityAware(reoptimize=x).reoptimize for x in (1, 0, 2.5, None)] == [
        True,
        False,
        True,
        False,
    ]
    with pytest.raises(TypeError):
        ngc.CapacityAware(flow_counts=2)


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


class _BigInt(int):
    pass


class _OwnFloat(int):
    # An int subclass whose own __float__ gives 2.5, whatever its integer value.
    def __float__(self):
        return 2.5


class _OverflowingFloat(np.float32):
    # A NumPy floating scalar whose own __float__ raises OverflowError.
    def __float__(self):
        raise OverflowError("own conversion")


class _OverflowingIntFloat(int):
    # An int subclass whose own __float__ raises OverflowError.
    def __float__(self):
        raise OverflowError("own conversion")


class _OverflowingIndex:
    # An integer-like object whose own __index__ raises OverflowError.
    def __index__(self):
        raise OverflowError("own conversion")


def test_float_fields_hold_only_finite_floats():
    # Stretch's factor and DiminishingReturns' fraction are C++ doubles. An integer with no finite
    # double (at either sign, as a Python int or through __index__) is outside the field's domain:
    # ValueError naming the field. The last integer that rounds to the largest finite double is
    # kept as that double; the next one rounds past it.
    top = sys.float_info.max
    last = 2**1024 - 2**970 - 1
    assert ngc.Stretch(last).factor == top
    assert ngc.DiminishingReturns(1, -last).fraction == -top

    class Index:
        def __index__(self):
            return 10**400

    too_large = "must be a number in the float range (got an integer beyond the largest finite float)"
    for bad in (last + 1, -(last + 1), 10**400, Index()):
        with pytest.raises(ValueError, match=exact("Stretch factor " + too_large)):
            ngc.Stretch(bad)
        with pytest.raises(
            ValueError, match=exact("diminishing fraction " + too_large)
        ):
            ngc.DiminishingReturns(1, bad)
    # So is an int subclass that keeps int's own conversion; one with its own __float__ is
    # converted by that first, as float() does.
    with pytest.raises(ValueError, match=exact("Stretch factor " + too_large)):
        ngc.Stretch(_BigInt(10**400))
    assert ngc.Stretch(_OwnFloat(10**400)).factor == 2.5
    # Any other failure of the object's own conversion reaches the caller unchanged, an
    # OverflowError raised by its own __float__ or __index__ too: no integer was returned.
    with pytest.raises(_ConversionFailed):
        ngc.Stretch(_FailingFloat(2.0))
    with pytest.raises(_ConversionFailed):
        ngc.DiminishingReturns(1, _FailingFloat(0.5))
    for bad in (_OverflowingFloat(2.0), _OverflowingIndex(), _OverflowingIntFloat(2)):
        with pytest.raises(OverflowError, match=exact("own conversion")):
            ngc.Stretch(bad)
        with pytest.raises(OverflowError, match=exact("own conversion")):
            ngc.DiminishingReturns(1, bad)


def test_validation_messages_before_binding():
    g = diamond()
    fg = ngc.FlowLedger(g)

    def rejects(cfg, msg, src=0):
        with pytest.raises(ValueError, match=exact(msg)):
            ngc.FlowPolicy(fg, src, 3, cfg)

    stretch = "max_path_cost Stretch factor must be finite and in [1, inf) (got {})"
    rejects(
        config(ngc.CostOnly(), ngc.FlowPlacement.MAX_FLOW, ngc.Admission.LOSSY),
        LOSSY_KERNEL,
    )
    rejects(config(ngc.CapacityAware(), admission=ngc.Admission.LOSSY), LOSSY_ROUTING)
    rejects(config(ngc.Global()), "Global routing requires placement = MaxFlow")
    rejects(
        config(ngc.CapacityAware(flow_count=0)),
        "CapacityAware flow_count must be >= 1 (got 0)",
    )
    rejects(
        config(ngc.CapacityAware(flow_count=2, max_flow_count=1)),
        "CapacityAware max_flow_count (1) must be >= flow_count (2)",
    )
    rejects(
        config(ngc.CapacityAware(max_path_cost=ngc.Stretch(0.5))),
        stretch.format("0.500000"),
    )
    rejects(
        config(ngc.CapacityAware(max_path_cost=ngc.Stretch(math.nan))),
        stretch.format("nan"),
    )
    rejects(
        config(ngc.CapacityAware(diminishing=ngc.DiminishingReturns(1, 1.0))),
        "diminishing fraction must be in (0, 1) (got 1.000000)",
    )
    rejects(config(ngc.Pinned([])), "Pinned: an empty bundle list is rejected")
    b = dag_of(g, ids(g, [0, 1]))
    rejects(
        config(ngc.Pinned([b, b])),
        "Pinned: duplicate bundle 1 (equal to bundle 0 before pruning)",
    )
    rejects(
        config(ngc.Pinned([dag_of(g, ids(g, [0, 3]))])),
        "Pinned: malformed bundle 0, no src -> dst walk (src 0, dst 3)",
    )
    rejects(
        config(ngc.Pinned([b])),
        "Pinned routing requires single endpoints (bundles bind one source and one "
        "destination)",
        src=[0, 1],
    )
    rejects(
        COSTONLY,
        "node 3 is in both endpoint sets: sources and destinations must not overlap",
        src=[0, 3],
    )
    rejects(COSTONLY, "sources: the node set is empty", src=[])
    rejects(COSTONLY, f"sources: node {L62} out of node range [0, 4)", src=[L62])
    with pytest.raises(
        ValueError, match=exact("foreign Eligibility: built for another graph")
    ):
        ngc.FlowPolicy(fg, 0, 3, COSTONLY, diamond().eligibility(excluded_edges=[0]))
    with pytest.raises(
        ValueError, match=exact("foreign DagCache: bound to another graph")
    ):
        ngc.FlowPolicy(fg, 0, 3, COSTONLY, None, ngc.DagCache(diamond()))
    # Arguments of the wrong kind: a ledger of None and non-integer node ids are ValueError; a
    # ledger, config, elig or cache of the wrong class and an unknown keyword are TypeError.
    with pytest.raises(ValueError, match=exact("ledger must be a FlowLedger")):
        ngc.FlowPolicy(None, 0, 3, COSTONLY)
    with pytest.raises(
        ValueError, match=exact("sources[0] must be an integer (got float)")
    ):
        ngc.FlowPolicy(fg, [0.0], 3, COSTONLY)
    with pytest.raises(
        ValueError, match=exact("destinations must be a sequence of integers")
    ):
        ngc.FlowPolicy(fg, 0, 3.0, COSTONLY)
    for wrong in (
        lambda: ngc.FlowPolicy(g, 0, 3, COSTONLY),
        lambda: ngc.FlowPolicy(fg, 0, 3, ngc.CostOnly()),
        lambda: ngc.FlowPolicy(fg, 0, 3, COSTONLY, ngc.DagCache(g)),
        lambda: ngc.FlowPolicy(fg, 0, 3, COSTONLY, None, g.eligibility()),
        lambda: ngc.FlowPolicy(fg, 0, 3, COSTONLY, cash=None),
    ):
        with pytest.raises(TypeError):
            wrong()
    fg.reset()  # nothing bound
    pol = ngc.FlowPolicy(fg, 0, 3, COSTONLY)
    with pytest.raises(
        TypeError, match=exact("outputs must be an OutputOptions (got CostOnly)")
    ):
        pol.place(1, ngc.CostOnly())
    assert pol.records() == []


LOSSY_KERNEL = (
    "Lossy admission is rejected with the MaxFlow kernel: a controller never drops "
    "(EqualBalanced or Weighted only)"
)
HEADROOM_KERNEL = (
    "the Weighted kernel is rejected with priority Headroom (weights are units)"
)
LOSSY_ROUTING = (
    "Lossy admission is rejected with Global and CapacityAware routing (hop-by-hop "
    "forwarding without admission control neither grows nor reroutes)"
)
P = ngc.FlowPlacement
A = ngc.Admission


def zero_flow_count():
    return ngc.CapacityAware(flow_count=0)


@pytest.mark.parametrize(
    ("routing", "placement", "admission", "priority", "msg"),
    [
        (ngc.Global, P.MAX_FLOW, A.LOSSY, "CAPACITY", LOSSY_KERNEL),
        (zero_flow_count, P.MAX_FLOW, A.LOSSY, "CAPACITY", LOSSY_KERNEL),
        (ngc.Global, P.WEIGHTED, A.LOSSLESS, "HEADROOM", HEADROOM_KERNEL),
        (zero_flow_count, P.WEIGHTED, A.LOSSY, "HEADROOM", HEADROOM_KERNEL),
        (zero_flow_count, P.EQUAL_BALANCED, A.LOSSY, "CAPACITY", LOSSY_ROUTING),
        (ngc.Global, P.WEIGHTED, A.LOSSY, "CAPACITY", LOSSY_ROUTING),
    ],
)
def test_several_defects_report_the_first_in_the_contract_order(
    routing, placement, admission, priority, msg
):
    # Several defects at once: the kernel combinations (Lossy with MaxFlow,
    # Weighted with Headroom) before the routing rules (Lossy with Global or CapacityAware,
    # Global without MaxFlow, CapacityAware's counts); one exact ValueError text at
    # construction and at the three evaluators (batch_capacity's precheck included); nothing is
    # bound.
    g = diamond()
    fg = ngc.FlowLedger(g)
    cfg = config(
        routing(), placement, admission, sel(prio=getattr(ngc.EdgePriority, priority))
    )
    with pytest.raises(ValueError) as e:
        ngc.FlowPolicy(fg, 0, 3, cfg)
    assert str(e.value) == msg
    with pytest.raises(ValueError) as e:
        ngc.capacity(g, 0, 3, ngc.CapacityOptions(cfg))
    assert str(e.value) == msg
    with pytest.raises(ValueError) as e:
        ngc.batch_capacity(g, [(0, 3)], ngc.CapacityOptions(cfg))
    assert str(e.value) == msg
    with pytest.raises(ValueError) as e:
        ngc.sensitivity(g, 0, 3, ngc.SensitivityOptions(cfg), 1)
    assert str(e.value) == msg
    fg.reset()  # nothing bound


def test_enum_values_outside_their_enumeration_are_rejected():
    # An enum built from an integer outside the enumeration is a domain error named
    # by its field, at every constructor that takes it, so it never reaches a policy or record.
    with pytest.raises(
        ValueError, match="tie_break: 99 is not a value of EdgeTieBreak"
    ):
        ngc.PathSelection(ngc.EdgeTieBreak(99), ngc.EdgePriority.CAPACITY)
    with pytest.raises(ValueError, match="priority: 99 is not a value of EdgePriority"):
        ngc.PathSelection(ngc.EdgeTieBreak.DETERMINISTIC, ngc.EdgePriority(99))
    for routing in (ngc.CostOnly(), ngc.CapacityAware(), ngc.Global()):
        with pytest.raises(
            ValueError, match="placement: 99 is not a value of FlowPlacement"
        ):
            config(routing, ngc.FlowPlacement(99))
        with pytest.raises(
            ValueError, match="admission: 99 is not a value of Admission"
        ):
            config(routing, ngc.FlowPlacement.MAX_FLOW, ngc.Admission(99))
    assert len(ngc.FlowPlacement.__members__) == 3  # the loop below is not empty
    for valid in ngc.FlowPlacement.__members__.values():
        assert config(ngc.CostOnly(), valid).placement == valid


def test_placement_result_fields_and_outputs():
    g = diamond()
    fg = ngc.FlowLedger(g)
    pol = ngc.FlowPolicy(fg, [0], [3], COSTONLY)
    r = pol.place(4)
    assert (r.placed, r.remaining, r.termination) == (4, 0, T.COMPLETE)
    assert r.local_outcomes == [(0, T.COMPLETE)]
    assert (r.routes_changed, r.weights_changed, r.drops) == (True, False, [])
    assert r.cost_distribution is None
    assert (
        r.edge_flows is None
        and r.residual is None
        and r.bottleneck is None
        and r.min_cut is None
    )
    r = pol.place(INF, ALL_OUT)
    assert (r.placed, r.remaining, r.termination) == (3, INF, T.EXHAUSTED)
    assert r.edge_flows.dtype == np.int64 and r.residual.dtype == np.int64
    assert list(r.residual) == list(fg.residual())
    own = [0] * 4
    for rec in pol.records():
        for e, x in rec.edges:
            own[e] += x
    assert list(r.edge_flows) == own
    assert r.bottleneck == fg.bottleneck(0)
    # By hand: the flow's 7 units split equally at node 0, the leftover unit to the lower
    # EdgeId: 4 on 0 -> 1 -> 3 and 3 on 0 -> 2 -> 3, the most that keeps 0 -> 2 within its 3.
    # Only 0 -> 2 is saturated where it leaves what the source still reaches.
    assert r.edge_flows.tolist() == to_ids(g, [4, 4, 3, 3])
    assert r.residual.tolist() == to_ids(g, [1, 1, 0, 0])
    assert r.bottleneck == ids(g, [2])
    assert r.min_cut is None  # Global only
    assert pol.placed_demand() == 7
    with pytest.raises(
        ValueError,
        match=exact(
            "requests are integer units or +inf (a physical request goes through "
            "to_units first; got 1.5)"
        ),
    ):
        pol.place(1.5)
    with pytest.raises(_ConversionFailed):  # a float request whose own conversion fails
        pol.place(_FailingFloat(1.0))
    for bad in (_OverflowingFloat(1.0), _OverflowingIndex()):  # its OverflowError too
        with pytest.raises(OverflowError, match=exact("own conversion")):
            pol.place(bad)
    with pytest.raises(ValueError, match=r"after \+inf"):
        pol.place(np.int64(1))
    pol.remove()
    assert pol.records() == [] and list(fg.allocated()) == [0, 0, 0, 0]
    assert pol.place(np.int64(2)).placed == 2
    # Global: the marginal curve; min_cut with bottleneck.
    fg2 = ngc.FlowLedger(diamond(5, 3, 2))
    gl = ngc.FlowPolicy(fg2, 0, 3, GLOBAL)
    r = gl.place(INF, ngc.OutputOptions(bottleneck=True))
    assert r.cost_distribution == [(2, 5), (4, 3)]
    assert sorted(r.min_cut) == ids(diamond(5, 3, 2), [0, 2])
    assert [rec.delivered for rec in gl.records()] == [8]


def three_routes(a, b, c):
    # Three disjoint routes 0 -> x -> 4 (input order): A through 1 of cost 1 + 1 = 2, B through 2
    # of cost 2 + 1 = 3, C through 3 of cost 2 + 2 = 4; capacities a, b, c on both edges.
    return build(
        5,
        [
            (0, 1, a, 1),
            (1, 4, a, 1),
            (0, 2, b, 2),
            (2, 4, b, 1),
            (0, 3, c, 2),
            (3, 4, c, 2),
        ],
    )


def routes_of(g, policy):
    """Per record: (route cost, delivered, the bundle's input edges, units per input edge)."""
    pos = g.input_index()
    return [
        (
            r.route_cost,
            r.delivered,
            sorted(int(pos[e]) for e in r.bundle.via_edges()),
            sorted((int(pos[e]), x) for e, x in r.edges),
        )
        for r in policy.records()
    ]


def gated(kernel, ceiling=None, floor=0):
    return config(
        ngc.CapacityAware(max_path_cost=ceiling, min_flow_bandwidth=floor), kernel
    )


# Each route of three_routes(2, 2, 2) carrying its 2 units, as routes_of reports it.
A2 = (2, 2, [0, 1], [(0, 2), (1, 2)])
B2 = (3, 2, [2, 3], [(2, 2), (3, 2)])
C2 = (4, 2, [4, 5], [(4, 2), (5, 2)])
KERNELS = [ngc.FlowPlacement.MAX_FLOW, ngc.FlowPlacement.EQUAL_BALANCED]


@pytest.mark.parametrize("kernel", KERNELS)
@pytest.mark.parametrize(
    "ceiling, placed, routes",
    [
        (ngc.Absolute(3), 4, [A2, B2]),
        (ngc.Stretch(1.5), 4, [A2, B2]),
        (ngc.Absolute(4), 6, [A2, B2, C2]),
        (ngc.Stretch(2.0), 6, [A2, B2, C2]),
        (None, 6, [A2, B2, C2]),
        (ngc.Absolute(1), 0, []),
    ],
)
def test_path_cost_ceiling_excludes_a_costlier_route(kernel, ceiling, placed, routes):
    # CapacityAwareGates.PathCostCeilingExcludesACostlierRoute through the binding, derived by
    # hand on three_routes(2, 2, 2): A (cost 2) is created and filled to 2,
    # growth adds B (cost 3; EqualBalanced: cohort [2, 2]), and C (cost 4) only when the
    # inclusive ceiling admits it: Absolute(3) and Stretch(1.5) (1.5 * 2 = 3) exclude it
    # (placed 4, exhausted), Absolute(4) and Stretch(2) admit it (placed 6). Absolute(1) admits
    # no route: placed 0, no_route.
    g = three_routes(2, 2, 2)
    p = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 4, gated(kernel, ceiling))
    r = p.place(10)
    term = T.NO_ROUTE if placed == 0 else T.EXHAUSTED
    assert (r.placed, r.remaining, r.termination) == (placed, 10 - placed, term)
    assert r.local_outcomes == [(0, term)]
    assert (r.routes_changed, r.weights_changed) == (bool(routes), False)
    assert routes_of(g, p) == routes


@pytest.mark.parametrize("kernel", KERNELS)
@pytest.mark.parametrize("factor, placed, routes", [(1.5, 2, [B2]), (2.0, 4, [B2, C2])])
def test_stretch_baseline_is_the_static_shortest_cost(kernel, factor, placed, routes):
    # A third party holds both units of A: the static shortest cost stays 2 (the Stretch
    # baseline ignores the residual), so Stretch(1.5) admits B (3) and not C (4): placed 2. A
    # baseline on the residual (3) would admit C under 4.5. Stretch(2) admits C: placed 4.
    g = three_routes(2, 2, 2)
    ledger = ngc.FlowLedger(g)
    third = ngc.FlowPolicy(
        ledger,
        0,
        4,
        config(ngc.Pinned([dag_of(g, ids(g, [0, 1]))]), ngc.FlowPlacement.MAX_FLOW),
    )
    assert third.place(2).placed == 2
    p = ngc.FlowPolicy(ledger, 0, 4, gated(kernel, ngc.Stretch(factor)))
    r = p.place(10)
    assert (r.placed, r.termination) == (placed, T.EXHAUSTED)
    assert routes_of(g, p) == routes


@pytest.mark.parametrize(
    "kernel, floor, placed, routes",
    [
        (ngc.FlowPlacement.MAX_FLOW, 3, 4, [(2, 4, [0, 1], [(0, 4), (1, 4)])]),
        (ngc.FlowPlacement.MAX_FLOW, 2, 6, [(2, 4, [0, 1], [(0, 4), (1, 4)]), B2]),
        (ngc.FlowPlacement.MAX_FLOW, 5, 0, []),
        (ngc.FlowPlacement.EQUAL_BALANCED, 3, 4, [(2, 4, [0, 1], [(0, 4), (1, 4)])]),
        (
            ngc.FlowPlacement.EQUAL_BALANCED,
            2,
            5,
            [(2, 3, [0, 1], [(0, 3), (1, 3)]), B2],
        ),
        (ngc.FlowPlacement.EQUAL_BALANCED, 5, 0, []),
    ],
)
def test_min_flow_bandwidth_stops_growth_below_the_floor(kernel, floor, placed, routes):
    # CapacityAwareGates.MinFlowBandwidthStopsGrowthBelowTheFloor through the binding, derived
    # by hand (every search gated at residual >= max(1, floor)) on three_routes(4, 2, 0): A of
    # 4 units at cost 2, B of 2 at cost 3. Floor 3: A is filled to 4, B (residual 2) is below
    # the floor, growth stops: placed 4. Floor 2: B is added: MaxFlow 4 + 2; EqualBalanced the
    # cohort [3, 2] (round robin, the odd unit to the first flow; [3, 3] exceeds B), then A's
    # residual 1 is below the floor: placed 5. Floor 5: no route at all, no_route.
    g = three_routes(4, 2, 0)
    p = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 4, gated(kernel, floor=floor))
    r = p.place(10)
    assert (r.placed, r.termination) == (
        placed,
        T.NO_ROUTE if placed == 0 else T.EXHAUSTED,
    )
    assert routes_of(g, p) == routes


def test_lifetimes_keep_alive_and_orphans():
    g = diamond()
    fg = ngc.FlowLedger(g)
    el = g.eligibility(excluded_edges=ids(g, [2]))
    cache = ngc.DagCache(g)
    pol = ngc.FlowPolicy(fg, 0, 3, COSTONLY, el, cache)
    held = [weakref.ref(x) for x in (fg, el, cache)]
    del fg, el, cache
    gc.collect()
    assert all(r() is not None for r in held)
    assert (
        pol.place(10).placed == 5
    )  # the ledger, the Eligibility and the cache are alive
    del pol
    gc.collect()
    assert all(r() is None for r in held)  # the policy alone held them
    fg2 = ngc.FlowLedger(g)
    p = ngc.FlowPolicy(fg2, 0, 3, COSTONLY)
    p.place(2)
    del p
    gc.collect()
    if HOOKS:
        assert fg2._state()["orphans"] == [0]
    assert fg2.records()[0].delivered == 2
    with pytest.raises(ValueError, match="bound"):
        q = ngc.FlowPolicy(fg2, 0, 3, COSTONLY)
        fg2.reset()
    del q
    gc.collect()
    fg2.reset()
    assert fg2.records() == []


@needs_hooks
def test_owner_hook_is_the_owner_of_the_records():
    # The private hook FlowPolicy._owner reads the ledger-local owner id: one per policy of a
    # ledger, known before the policy has a record, and carried by its records as FlowRecord.owner.
    g = diamond()
    fg = ngc.FlowLedger(g)
    policies = [ngc.FlowPolicy(fg, 0, 3, COSTONLY) for _ in range(3)]
    owners = [p._owner() for p in policies]
    assert len(set(owners)) == 3
    for p in reversed(policies):
        assert p.place(1).placed == 1
    for p, owner in zip(policies, owners, strict=True):
        assert [r.owner for r in p.records()] == [owner]


def test_pinned_usable_bundles_and_records():
    g = diamond()
    fg = ngc.FlowLedger(g)
    r1, r2 = dag_of(g, ids(g, [0, 1])), dag_of(g, ids(g, [2, 3]))
    el = g.eligibility(excluded_edges=ids(g, [3]))
    pol = ngc.FlowPolicy(fg, 0, 3, config(ngc.Pinned([r1, r2])), el)
    assert pol.usable_bundles() == 1
    r = pol.place(9)
    assert (r.placed, r.termination) == (5, T.EXHAUSTED)
    rec = pol.records()[0]
    assert rec.bundle == r1 and pol.record(rec.flow_id).delivered == 5
    # record() reads the policy's own records only, not another policy's on the same ledger.
    other = ngc.FlowPolicy(fg, 0, 3, config(ngc.Pinned([r1])))
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
        fg, 0, 3, config(ngc.Pinned([r2])), g.eligibility(excluded_edges=ids(g, [3]))
    )
    assert all_down.usable_bundles() == 0
    r = all_down.place(1)
    assert (r.placed, r.termination, r.local_outcomes) == (0, T.NO_ROUTE, [])


def test_dag_cache_statistics_and_projection():
    g = diamond()
    cache = ngc.DagCache(g)
    fg = ngc.FlowLedger(g)
    a = ngc.FlowPolicy(fg, 0, 3, COSTONLY, None, cache)
    b = ngc.FlowPolicy(
        fg,
        0,
        3,
        config(
            ngc.CostOnly(),
            selection=sel(
                ngc.EdgeTieBreak.PREFER_HIGHER_PRIORITY,
                RES,
                False,
                False,
            ),
        ),
        None,
        cache,
    )
    s = cache.stats()
    assert (s.hits, s.misses, s.entries, s.live_values) == (1, 1, 1, 1)
    assert s.bytes > 0
    assert a.place(2).placed == 2 and b.place(2).placed == 2
    cache.clear()
    s = cache.stats()
    assert (s.hits, s.misses, s.entries, s.bytes) == (1, 1, 0, 0)
    assert s.live_values == 1  # held by the policies
    assert cache.graph_uid() == g.uid()


def test_capacity_batch_and_workers(monkeypatch):
    g = diamond(5, 3)
    opts = ngc.CapacityOptions(GLOBAL, outputs=ngc.OutputOptions(edge_flows=True))
    r = ngc.capacity(g, 0, 3, opts)
    assert (r.placed, r.termination, r.cost_distribution) == (8, T.EXHAUSTED, [(2, 8)])
    assert int(np.sum(r.edge_flows)) == 16
    common = g.eligibility(excluded_edges=ids(g, [3]))
    pairs = [(0, 3), ([0], [3]), (1, 3), (0, [1, 2])]
    per = [None, g.eligibility(excluded_edges=ids(g, [0])), None, None]
    want = [5, 0, 5, 8]
    for workers in (1, 3, None):
        got = ngc.batch_capacity(
            g, pairs, ngc.CapacityOptions(GLOBAL, elig=common), per, workers
        )
        assert [x.placed for x in got] == want
    # workers=None reads NGRAPH_CORE_BATCH_THREADS. The count changes only how many threads run,
    # never a result, so this checks only that the default path with the variable set gives the
    # same results (its parsing: the C++ test Evaluators.DefaultWorkersParseTheEnvironment).
    monkeypatch.setenv("NGRAPH_CORE_BATCH_THREADS", "4")
    got = ngc.batch_capacity(g, pairs, ngc.CapacityOptions(GLOBAL, elig=common), per)
    assert [x.placed for x in got] == want
    assert ngc.batch_capacity(g, [], ngc.CapacityOptions(GLOBAL)) == []
    with pytest.raises(ValueError, match="node 7"):
        ngc.batch_capacity(
            g, [(0, 3), (0, 7), (0, 8)], ngc.CapacityOptions(GLOBAL), workers=4
        )
    with pytest.raises(
        ValueError, match=exact("per_pair_elig has 1 entries for 4 pairs")
    ):
        ngc.batch_capacity(g, pairs, ngc.CapacityOptions(GLOBAL), [None])
    with pytest.raises(ValueError, match="Pinned"):
        ngc.batch_capacity(
            g,
            pairs,
            ngc.CapacityOptions(config(ngc.Pinned([dag_of(g, ids(g, [0, 1]))]))),
        )
    cache = ngc.DagCache(g)
    eb = ngc.CapacityOptions(COSTONLY, cache=cache)
    with pytest.raises(ValueError, match="workers == 1"):
        ngc.batch_capacity(g, pairs, eb, workers=2)
    got = ngc.batch_capacity(
        g, pairs, eb
    )  # None with a cache: 1 worker, even with the env set
    assert [x.placed for x in got] == [7, 7, 5, 7]  # EB: [4, 3] over links of 5 and 3
    with pytest.raises(ValueError, match="workers"):
        ngc.batch_capacity(g, pairs, ngc.CapacityOptions(GLOBAL), workers=0)


def test_sensitivity_candidates_and_signed_losses(monkeypatch):
    # Global: an unsaturated critical edge appears.
    g = build(3, [(0, 1, 9, 1), (1, 2, 4, 1), (0, 2, 2, 5)])
    s = ngc.sensitivity(g, 0, 2, ngc.SensitivityOptions(GLOBAL))
    assert s.baseline.placed == 6
    pos = g.input_index()
    assert {int(pos[e.edge]): e.loss for e in s.entries} == {0: 4, 1: 4, 2: 2}
    # workers=None reads NGRAPH_CORE_SENSITIVITY_THREADS: the same results on 3 workers (only that
    # the default path with the variable set works; its parsing is tested in C++).
    monkeypatch.setenv("NGRAPH_CORE_SENSITIVITY_THREADS", "3")
    assert plain(ngc.sensitivity(g, 0, 2, ngc.SensitivityOptions(GLOBAL))) == plain(s)
    # A policy model: every allowed edge is a candidate, a zero-capacity one too. Equal splitting
    # over parallel links of 1, 10 and 0 fits 2 units (a third would cross the zero link);
    # without the thin link 1 fits, without the wide link 1, without the zero link 2.
    p = build(2, [(0, 1, 1, 1), (0, 1, 10, 1), (0, 1, 0, 1)])
    opts = ngc.SensitivityOptions(COSTONLY)
    s = ngc.sensitivity(p, 0, 1, opts, workers=2)
    assert s.baseline.placed == 2
    pos = p.input_index()
    assert {int(pos[e.edge]): e.loss for e in s.entries} == {0: 1, 1: 1, 2: 0}
    # A negative loss on a zero-capacity link (G8): the only shortest route 0 -> 1 has capacity
    # 0, so equal splitting carries 0; without it the route 0 -> 2 -> 1 carries 5.
    g8 = build(3, [(0, 1, 0, 1), (0, 2, 5, 1), (2, 1, 5, 1)])
    s8 = ngc.sensitivity(g8, 0, 1, opts)
    assert s8.baseline.placed == 0
    pos8 = g8.input_index()
    assert {int(pos8[e.edge]): e.loss for e in s8.entries} == {0: -5, 1: 0, 2: 0}
    with pytest.raises(ValueError, match="both endpoint sets"):
        ngc.sensitivity(p, 0, 0, opts)
    # The baseline is capacity() with the requested outputs.
    for cfg in (GLOBAL, COSTONLY):
        s = ngc.sensitivity(g, 0, 2, ngc.SensitivityOptions(cfg, outputs=ALL_OUT))
        want = ngc.capacity(g, 0, 2, ngc.CapacityOptions(cfg, outputs=ALL_OUT))
        assert want.edge_flows is not None and want.residual is not None
        assert plain(s.baseline) == plain(want)


def test_evaluator_options_and_results():
    # The options hold what they were given (None: everything allowed, no supplied cache, no
    # output); the results have no public constructor.
    g = diamond()
    el = g.eligibility(excluded_edges=ids(g, [0]))
    cache = ngc.DagCache(g)
    for cls in (ngc.CapacityOptions, ngc.SensitivityOptions):
        o = cls(GLOBAL, elig=el, cache=cache, outputs=ngc.OutputOptions(residual=True))
        assert (
            o.config.placement == GLOBAL.placement and o.elig is el and o.cache is cache
        )
        out = o.outputs
        assert (out.edge_flows, out.residual, out.bottleneck) == (False, True, False)
        d = cls(COSTONLY)
        assert d.elig is None and d.cache is None
        assert (d.outputs.edge_flows, d.outputs.residual, d.outputs.bottleneck) == (
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
    building = ngc.DagCache(g)
    bopts = ngc.CapacityOptions(COSTONLY, cache=building)
    sopts = ngc.SensitivityOptions(COSTONLY, cache=building)
    want = [r.placed for r in ngc.batch_capacity(g, [(0, 3), (1, 3)], bopts, workers=1)]
    assert want == [7, 5]
    assert [r.placed for r in ngc.batch_capacity(g, [(0, 3), (1, 3)], bopts)] == want
    base = plain(ngc.sensitivity(g, 0, 3, sopts, workers=1))
    assert plain(ngc.sensitivity(g, 0, 3, sopts)) == base
    for call in (
        lambda: ngc.batch_capacity(g, [(0, 3)], bopts, workers=2),
        lambda: ngc.sensitivity(g, 0, 3, sopts, workers=2),
    ):
        with pytest.raises(ValueError, match="unless workers == 1"):
            call()
    with pytest.raises(ValueError, match=exact("workers must be >= 1 (got 0)")):
        ngc.sensitivity(g, 0, 3, ngc.SensitivityOptions(GLOBAL), workers=0)


def test_python_keeps_what_a_policy_and_an_option_need_alive():
    # Once the caller drops them, the evaluators' options still apply their Eligibility and cache,
    # the ledger still reads its graph, and a policy still places on its ledger. Graphs of the
    # same shape are built in between, so the memory of a released value would be reused.
    g = diamond()
    want = [g.capacity_units().tolist(), g.cost().tolist(), g.edge_dst().tolist()]
    ledger = ngc.FlowLedger(g)
    policy = ngc.FlowPolicy(ledger, 0, 3, COSTONLY)
    opts = ngc.CapacityOptions(
        COSTONLY, elig=g.eligibility(excluded_edges=ids(g, [2])), cache=ngc.DagCache(g)
    )
    gc.collect()
    reused = [diamond(7, 9, 2) for _ in range(64)]
    assert ngc.capacity(g, 0, 3, opts).placed == 5  # without the edge 0 -> 2
    assert opts.cache.stats().misses == 1  # the call's lookup
    del g, opts
    gc.collect()
    reused += [diamond(7, 9, 2) for _ in range(64)]
    kept = ledger.graph()
    assert [
        kept.capacity_units().tolist(),
        kept.cost().tolist(),
        kept.edge_dst().tolist(),
    ] == want
    del kept, ledger
    gc.collect()
    reused += [diamond(7, 9, 2) for _ in range(64)]
    assert policy.place(10).placed == 7  # EqualBalanced over routes of 5 and 3: [4, 3]
    assert len(reused) == 192


def test_exception_classes_of_conversions():
    # An integer beyond the float range where a float is expected: the resolution (a plain float
    # parameter) cannot be converted, a TypeError; a request (integer units) is outside its
    # domain, a ValueError. A request's own conversion failure reaches the caller unchanged.
    for call in (
        lambda: ngc.Graph(2, [0], [1], [5], [1], 10**400),
        lambda: ngc.to_units(1.0, 10**400),
        lambda: ngc.as_float(1, 10**400),
    ):
        with pytest.raises(TypeError, match="incompatible"):
            call()
    g = build(2, [(0, 1, 5, 1)])
    p = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 1, COSTONLY)
    for big in (10**400, -(10**400)):
        with pytest.raises(
            ValueError,
            match=exact(f"request (units) = {big} is out of range [0, {L62 - 1}]"),
        ):
            p.place(big)
    # An index's own failure (a float's: test_placement_result_fields_and_outputs).
    with pytest.raises(_ConversionFailed):
        p.place(_FailingIndex())
    assert p.placed_demand() == 0


@needs_hooks
@pytest.mark.parametrize("routing", ["costonly", "global"])
@pytest.mark.parametrize("traced", [False, True])
def test_python_publication_is_all_or_nothing(routing, traced):
    # Every Python object of a placement exists before anything is published, so the k-th
    # Python allocation failing raises with the ledger and the policy unchanged (an earlier
    # call's flow in the same epoch included), never crashes, and a retry places as a call
    # without failure does.
    cfg = GLOBAL if routing == "global" else COSTONLY
    g = build(2, [(0, 1, 10, 1)])

    def make():
        ledger = ngc.FlowLedger(g)
        policy = ngc.FlowPolicy(ledger, 0, 1, cfg)
        assert policy.place(1).placed == 1

        def snapshot():
            return [full_ledger_state(ledger), repr(policy._state())]

        return (policy._place_traced if traced else policy.place), (3,), snapshot

    assert sweep_python_allocation(make) >= 3  # not vacuous


@needs_hooks
def test_evaluator_results_survive_python_allocation_failures():
    # The evaluators publish nothing; their result objects are made through the checked path
    # too, so a failed Python allocation raises instead of crashing, and the call that
    # completes returns what a clean call returns.
    g = diamond(5, 3)
    calls = [
        (ngc.capacity, (g, 0, 3, ngc.CapacityOptions(GLOBAL))),
        (ngc.batch_capacity, (g, [(0, 3), (1, 3)], ngc.CapacityOptions(GLOBAL))),
        (ngc.sensitivity, (g, 0, 3, ngc.SensitivityOptions(GLOBAL))),
    ]
    for fn, args in calls:
        # An evaluator keeps no state between calls: nothing to snapshot.
        failures = sweep_python_allocation(lambda fn=fn, args=args: (fn, args, list))
        assert failures >= 3, fn.__name__  # not vacuous: every result object allocates


def test_batch_reports_the_lowest_failing_pair():
    # The arguments first, then each pair completely in pair order, across the
    # Python/C++ boundary.
    g = build(2, [(0, 1, 3, 7)])
    opts = ngc.CapacityOptions(GLOBAL)
    overlap = (
        "node 0 is in both endpoint sets: sources and destinations must not overlap"
    )
    out_of_range = "sources: node 2 out of node range [0, 2)"
    with pytest.raises(ValueError, match=exact(overlap)):
        ngc.batch_capacity(g, [(0, 0), (2, 1)], opts, workers=2)
    with pytest.raises(ValueError, match=exact(out_of_range)):
        ngc.batch_capacity(g, [(0, 1), (2, 1), (0, 0)], opts)
    other = build(2, [(0, 1, 3, 7)]).eligibility()
    with pytest.raises(
        ValueError,
        match=exact("foreign Eligibility: per_pair_elig[0] built for another graph"),
    ):
        ngc.batch_capacity(g, [(0, 1), (2, 1)], opts, per_pair_elig=[other, None])
    with pytest.raises(
        TypeError, match=exact("per_pair_elig entries must be Eligibility or None")
    ):
        ngc.batch_capacity(g, [(0, 1), (2, 1)], opts, per_pair_elig=[3, None])
    with pytest.raises(ValueError, match=exact(out_of_range)):
        ngc.batch_capacity(g, [(0, 1), (2, 1)], opts, per_pair_elig=[None, 3])
    # An argument error comes before every pair error.
    with pytest.raises(ValueError, match=exact("workers must be >= 1 (got 0)")):
        ngc.batch_capacity(g, [(0, 0)], opts, workers=0)


def test_batch_checks_an_explicit_empty_eligibility_list():
    # A given per_pair_elig has one entry per pair, also when it is empty.
    g = build(2, [(0, 1, 3, 7)])
    opts = ngc.CapacityOptions(GLOBAL)
    with pytest.raises(
        ValueError, match=exact("per_pair_elig has 0 entries for 1 pairs")
    ):
        ngc.batch_capacity(g, [(0, 1)], opts, per_pair_elig=[])
    assert ngc.batch_capacity(g, [], opts, per_pair_elig=[]) == []
    assert [
        r.placed for r in ngc.batch_capacity(g, [(0, 1)], opts, per_pair_elig=[None])
    ] == [3]


def test_global_zero_call_reports_the_cumulative_curve():
    # Global's curve is cumulative (its amounts sum to placed_demand()); a zero call
    # reports it ([] before the first placement and after remove()), changes nothing and claims
    # no cut.
    g = build(2, [(0, 1, 3, 7)])
    ledger = ngc.FlowLedger(g)
    p = ngc.FlowPolicy(ledger, 0, 1, GLOBAL)
    assert p.place(0).cost_distribution == []
    assert p.place(2).cost_distribution == [(7, 2)]
    z = p.place(0, ALL_OUT)
    assert (z.placed, z.cost_distribution, z.min_cut) == (0, [(7, 2)], None)
    assert p.placed_demand() == 2 and [int(x) for x in ledger.allocated()] == [2]
    p.remove()
    assert p.place(0).cost_distribution == []
