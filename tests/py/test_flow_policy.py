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
    CAP,
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
        ngc.FlowPolicyConfig(
            ngc.CostOnly(), ngc.FlowPlacement.MAX_FLOW, ngc.Admission.LOSSLESS, sel()
        )
    # Keyword-only; every field but the routing has a default: placement, selection and origin
    # None, admission LOSSLESS (contract 3.1). A wrong class is TypeError.
    plain = ngc.FlowPolicyConfig(routing=ngc.CostOnly())
    assert (plain.placement, plain.admission, plain.selection, plain.origin) == (
        None,
        ngc.Admission.LOSSLESS,
        None,
        None,
    )
    with pytest.raises(TypeError):
        ngc.FlowPolicyConfig(placement=ngc.FlowPlacement.MAX_FLOW)
    with pytest.raises(TypeError):
        ngc.FlowPolicyConfig(routing=ngc.CostOnly(), origin="EQUAL")
    with pytest.raises(TypeError):
        ngc.FlowPolicyConfig(routing=ngc.CostOnly(), placement=1)
    with pytest.raises(TypeError):
        ngc.FlowPolicyConfig(routing=ngc.Global(), bundles=None)  # an unknown keyword
    with pytest.raises(TypeError):
        ngc.FlowPolicyConfig(routing=ngc.CostOnly(), weights="EXPLICIT")
    with pytest.raises(TypeError, match="incompatible"):
        ngc.CostOnly(
            tunnels=ngc.Tunnels(2)
        )  # a field outside its routing variant: TypeError
    with pytest.raises(TypeError, match="incompatible"):
        ngc.Global(max_path_cost=ngc.Absolute(4))
    with pytest.raises(TypeError, match="routing"):
        config("CostOnly")
    # T1: the tunnel count (a Tunnels), the ceiling and the search minimum, every field as built.
    ca = ngc.CapacityAware(
        ngc.Tunnels(2), max_path_cost=ngc.Stretch(1.5), min_unreserved_bandwidth=3
    )
    assert (ca.tunnels, ca.max_path_cost, ca.min_unreserved_bandwidth) == (
        ngc.Tunnels(2),
        ngc.Stretch(1.5),
        3,
    )
    assert (ca.tunnels.count, ca.tunnels.up_to) == (2, False)
    assert ngc.Tunnels(3, up_to=True).up_to is True
    d = ngc.CapacityAware(ngc.Tunnels(1))
    assert (d.max_path_cost, d.min_unreserved_bandwidth) == (None, 0)
    assert ngc.CapacityAware(
        ngc.Tunnels(1), max_path_cost=ngc.Absolute(7)
    ).max_path_cost == ngc.Absolute(7)
    with pytest.raises(TypeError):
        ngc.CapacityAware()  # the tunnel count is required
    with pytest.raises(TypeError):
        ngc.CapacityAware(2)  # a Tunnels, not an int
    with pytest.raises(TypeError):
        ngc.CapacityAware(ngc.Tunnels(1), ngc.Absolute(1))  # keyword-only after tunnels
    with pytest.raises(TypeError):
        ngc.Tunnels(2, True)  # up_to is keyword-only
    for bad, message in (
        (2**70, "tunnel count = 1180591620717411303424 is out of range"),
        (True, "tunnel count must be an integer (got bool)"),
        (1.5, "tunnel count must be an integer (got float)"),
    ):
        with pytest.raises(
            ValueError, match=message.replace("(", r"\(").replace(")", r"\)")
        ):
            ngc.Tunnels(bad)
    with pytest.raises(
        ValueError,
        match=exact(
            "min_unreserved_bandwidth = 1208925819614629174706176 is out of range "
            "[-9223372036854775808, 9223372036854775807]"
        ),
    ):
        ngc.CapacityAware(ngc.Tunnels(1), min_unreserved_bandwidth=2**80)
    with pytest.raises(TypeError):
        ngc.Pinned([[0, 1]])
    cfg = ngc.FlowPolicyConfig(routing=ca)
    assert isinstance(cfg.routing, ngc.CapacityAware)
    assert cfg.routing.tunnels.count == 2
    s = sel(ngc.EdgeTieBreak.PREFER_HIGHER_PRIORITY, RES, False, True)
    cfg = ngc.FlowPolicyConfig(
        routing=ngc.Global(),
        placement=ngc.FlowPlacement.MAX_FLOW,
        admission=ngc.Admission.LOSSY,
        selection=s,
        origin=ngc.Origin.EQUAL,
    )  # a plain value: the combination is checked where it is used
    assert cfg.placement == ngc.FlowPlacement.MAX_FLOW
    assert cfg.admission == ngc.Admission.LOSSY
    assert cfg.selection == s
    assert cfg.origin == ngc.Origin.EQUAL
    assert ngc.PathSelection() == sel(ngc.EdgeTieBreak.DETERMINISTIC, CAP, None, True)
    assert ngc.PathSelection().multipath is None
    assert ngc.PathSelection() != sel(ngc.EdgeTieBreak.DETERMINISTIC, CAP, True, True)
    o = ngc.OutputOptions()
    names = (
        "edge_flows",
        "residual",
        "bottleneck",
        "reservation_edge_flows",
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
        g, 0, 3, [ngc.Route(0, b), ngc.Route(0, a)]
    )
    dark = g.eligibility(excluded_edges=ids(g, [0]))
    pinned = ngc.Pinned(program, black_holes=dark)
    assert pinned.program is program and pinned.black_holes is dark
    assert ngc.Pinned(program).black_holes is None
    # A policy's config is the configuration it was built with, its absent origin and selection
    # resolved to the model's defaults (Pinned takes none).
    pol = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 3, ngc.FlowPolicyConfig(routing=pinned))
    assert pol.config.placement is None and pol.config.selection is None
    assert pol.config.routing.program is program
    assert pol.config.origin is None
    for routing, placement, origin in (
        (ngc.CostOnly(), ngc.FlowPlacement.EQUAL_BALANCED, ngc.Origin.EQUAL),
        (ngc.CostOnly(), ngc.FlowPlacement.WEIGHTED, ngc.Origin.EQUAL),
        (ngc.CostOnly(), ngc.FlowPlacement.MAX_FLOW, ngc.Origin.OPTIMIZED),
        (ngc.Global(), None, ngc.Origin.OPTIMIZED),
    ):
        got = ngc.FlowPolicy(
            ngc.FlowLedger(g),
            0,
            3,
            ngc.FlowPolicyConfig(routing=routing, placement=placement),
        ).config
        assert got.origin == origin
        assert got.selection == (None if placement is None else ngc.PathSelection())
    got = ngc.FlowPolicy(
        ngc.FlowLedger(g),
        0,
        3,
        ngc.FlowPolicyConfig(routing=ngc.CapacityAware(ngc.Tunnels(1))),
    ).config
    assert (got.origin, got.selection, got.placement) == (
        ngc.Origin.EQUAL,
        ngc.PathSelection(),
        None,
    )
    # The constructors' own checks: an integer field outside its C++ type raises ValueError; a
    # field of the wrong class and an unknown keyword raise TypeError.
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
        ngc.CapacityAware(ngc.Tunnels(1), max_path_cost=3)
    with pytest.raises(TypeError):
        ngc.CapacityAware(ngc.Tunnels(1), min_flow_bandwidth=2)  # T0's field is gone


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
    # Stretch's factor is a C++ double. An integer with no finite
    # double (at either sign, as a Python int or through __index__) is outside the field's domain:
    # ValueError naming the field. The last integer that rounds to the largest finite double is
    # kept as that double; the next one rounds past it.
    top = sys.float_info.max
    last = 2**1024 - 2**970 - 1
    assert ngc.Stretch(last).factor == top
    assert ngc.Stretch(-last).factor == -top

    class Index:
        def __index__(self):
            return 10**400

    too_large = "must be a number in the float range (got an integer beyond the largest finite float)"
    for bad in (last + 1, -(last + 1), 10**400, Index()):
        with pytest.raises(ValueError, match=exact("Stretch factor " + too_large)):
            ngc.Stretch(bad)
    # So is an int subclass that keeps int's own conversion; one with its own __float__ is
    # converted by that first, as float() does.
    with pytest.raises(ValueError, match=exact("Stretch factor " + too_large)):
        ngc.Stretch(_BigInt(10**400))
    assert ngc.Stretch(_OwnFloat(10**400)).factor == 2.5
    # Any other failure of the object's own conversion reaches the caller unchanged, an
    # OverflowError raised by its own __float__ or __index__ too: no integer was returned.
    with pytest.raises(_ConversionFailed):
        ngc.Stretch(_FailingFloat(2.0))
    for bad in (_OverflowingFloat(2.0), _OverflowingIndex(), _OverflowingIntFloat(2)):
        with pytest.raises(OverflowError, match=exact("own conversion")):
            ngc.Stretch(bad)


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
    rejects(
        ngc.FlowPolicyConfig(
            routing=ngc.Global(), placement=ngc.FlowPlacement.MAX_FLOW
        ),
        "Global does not accept placement",
    )
    rejects(ngc.FlowPolicyConfig(routing=ngc.CostOnly()), "CostOnly requires placement")
    rejects(t1_config(0), "tunnel count must be in [1, 10000]")
    rejects(t1_config(10001), "tunnel count must be in [1, 10000]")
    rejects(
        ngc.FlowPolicyConfig(routing=ngc.CapacityAware(ngc.Tunnels(2, up_to=True))),
        "T1 requires an exact tunnel count",
    )
    rejects(t1_config(1, max_path_cost=ngc.Stretch(0.5)), stretch.format("0.500000"))
    rejects(t1_config(1, max_path_cost=ngc.Stretch(math.nan)), stretch.format("nan"))
    for bad in (-1, 2**62):
        rejects(
            t1_config(1, min_unreserved_bandwidth=bad),
            "min_unreserved_bandwidth must be an integer in [0, 2^62)",
        )
    b = dag_of(g, ids(g, [0, 1]))
    program = ngc.ForwardingProgram.from_routes(g, 0, 3, [ngc.Route(0, b)])
    rejects(
        ngc.FlowPolicyConfig(
            routing=ngc.Pinned(program), placement=ngc.FlowPlacement.MAX_FLOW
        ),
        "Pinned does not accept placement",
    )
    rejects(
        ngc.FlowPolicyConfig(routing=ngc.Pinned(program)),
        "program endpoints do not match policy endpoints",
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


LOSSY_KERNEL = "MAX_FLOW requires admission=LOSSLESS"
P = ngc.FlowPlacement
A = ngc.Admission
ORIGIN = ngc.Origin


def T1():
    return ngc.CapacityAware(ngc.Tunnels(1))


@pytest.mark.parametrize(
    ("routing", "placement", "admission", "selection", "origin", "msg"),
    [
        # The rows of the matrix, their columns left to right: placement, origin, admission,
        # selection (contract 3.2, 14.3, 14.4).
        # Global rejects any selection before its row.
        (
            ngc.Global,
            P.MAX_FLOW,
            A.LOSSY,
            True,
            ORIGIN.EQUAL,
            "Global does not accept selection",
        ),
        (
            ngc.Global,
            P.MAX_FLOW,
            A.LOSSY,
            False,
            ORIGIN.EQUAL,
            "Global does not accept placement",
        ),
        (ngc.Global, None, A.LOSSY, False, ORIGIN.EQUAL, LOSSY_KERNEL),
        (
            ngc.Global,
            None,
            A.LOSSLESS,
            True,
            ORIGIN.EQUAL,
            "Global does not accept selection",
        ),
        (
            ngc.CostOnly,
            None,
            A.LOSSY,
            True,
            ORIGIN.OPTIMIZED,
            "CostOnly requires placement",
        ),
        (
            ngc.CostOnly,
            P.EQUAL_BALANCED,
            A.LOSSY,
            True,
            ORIGIN.OPTIMIZED,
            "CostOnly with EQUAL_BALANCED or WEIGHTED requires origin=EQUAL",
        ),
        (ngc.CostOnly, P.MAX_FLOW, A.LOSSY, True, ORIGIN.EQUAL, LOSSY_KERNEL),
        (
            T1,
            P.MAX_FLOW,
            A.LOSSY,
            True,
            ORIGIN.OPTIMIZED,
            "T1 does not accept placement",
        ),
        (T1, None, A.LOSSY, True, ORIGIN.OPTIMIZED, "T1 requires origin=EQUAL"),
    ],
)
def test_several_defects_report_the_first_in_the_contract_order(
    routing, placement, admission, selection, origin, msg
):
    # Several defects at once: one exact ValueError text, the first in the matrix's order, at
    # construction and at the three evaluators (batch_capacity's precheck included); nothing is
    # bound.
    g = diamond()
    fg = ngc.FlowLedger(g)
    cfg = ngc.FlowPolicyConfig(
        routing=routing(),
        placement=placement,
        admission=admission,
        selection=sel() if selection else None,
        origin=origin,
    )
    with pytest.raises(ValueError) as e:
        ngc.FlowPolicy(fg, 0, 3, cfg)
    assert str(e.value) == msg
    with pytest.raises(ValueError) as e:
        ngc.capacity(g, 0, 3, ngc.QueryOptions(cfg))
    assert str(e.value) == msg
    with pytest.raises(ValueError) as e:
        ngc.batch_capacity(g, [(0, 3)], ngc.QueryOptions(cfg))
    assert str(e.value) == msg
    with pytest.raises(ValueError) as e:
        ngc.sensitivity(g, 0, 3, ngc.QueryOptions(cfg), 1)
    assert str(e.value) == msg
    fg.reset()  # nothing bound


def test_enum_values_outside_their_enumeration_are_rejected():
    # An enum built from an integer outside the enumeration is a domain error named by its field:
    # at PathSelection's constructor, and for FlowPolicyConfig (a plain value) wherever it is used,
    # first among its checks, so it never reaches a policy or record.
    with pytest.raises(
        ValueError, match="tie_break: 99 is not a value of EdgeTieBreak"
    ):
        ngc.PathSelection(ngc.EdgeTieBreak(99), ngc.EdgePriority.CAPACITY)
    with pytest.raises(ValueError, match="priority: 99 is not a value of EdgePriority"):
        ngc.PathSelection(ngc.EdgeTieBreak.DETERMINISTIC, ngc.EdgePriority(99))
    g = diamond()
    for routing in (ngc.CostOnly(), T1(), ngc.Global()):
        for cfg, message in (
            (
                ngc.FlowPolicyConfig(routing=routing, placement=ngc.FlowPlacement(99)),
                "placement: 99 is not a value of FlowPlacement",
            ),
            (
                ngc.FlowPolicyConfig(routing=routing, admission=ngc.Admission(99)),
                "admission: 99 is not a value of Admission",
            ),
            (
                ngc.FlowPolicyConfig(routing=routing, origin=ngc.Origin(99)),
                "origin: 99 is not a value of Origin",
            ),
        ):
            with pytest.raises(ValueError, match=exact(message)):
                ngc.FlowPolicy(ngc.FlowLedger(g), 0, 3, cfg)
            with pytest.raises(ValueError, match=exact(message)):
                ngc.capacity(g, 0, 3, ngc.QueryOptions(cfg))
    assert len(ngc.FlowPlacement.__members__) == 3  # the loop below is not empty
    for valid in ngc.FlowPlacement.__members__.values():
        assert config(ngc.CostOnly(), valid).placement == valid


def test_placement_result_fields_and_outputs():
    g = diamond()
    fg = ngc.FlowLedger(g)
    pol = ngc.FlowPolicy(fg, [0], [3], COSTONLY)
    r = pol.place(4)
    assert (r.placed, r.remaining, r.termination) == (4, 0, T.COMPLETE)
    assert not hasattr(r, "local_outcomes")
    # The epoch's first call installs the program: its routes and weights.
    assert (r.routes_changed, r.weights_changed, r.call_drops) == (True, True, [])
    assert (r.delivery_delta, r.call_source_loss, r.drops, r.source_loss) == (
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
    assert (row.source_loss, row.dropped) == (0, 0)
    with pytest.raises(TypeError):
        ngc.SourceFlow()
    with pytest.raises(AttributeError):
        row.delivered = 1
    assert r.cost_distribution is None and r.total_cost is None
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
    assert pol.records() == [] and list(fg.allocated()) == [0, 0, 0, 0]
    assert pol.place(np.int64(2)).placed == 2
    # Global: the marginal curve; min_cut with bottleneck.
    fg2 = ngc.FlowLedger(diamond(5, 3, 2))
    gl = ngc.FlowPolicy(fg2, 0, 3, GLOBAL)
    r = gl.place(INF, ngc.OutputOptions(bottleneck=True))
    assert r.cost_distribution == [(2, 5), (4, 3)] and r.min_cut is None
    gl.remove()
    r = gl.place(INF, ngc.OutputOptions(min_cut=True))
    assert r.min_cut == sorted(ids(diamond(5, 3, 2), [0, 2]))
    assert [rec.delivered for rec in gl.records()] == [8]


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


def pinned_over(g, *dags):
    """Pinned over a program of the split routes `dags` of head 0 to node 3."""
    program = ngc.ForwardingProgram.from_routes(
        g, 0, 3, [ngc.Route(0, d) for d in dags]
    )
    return ngc.FlowPolicyConfig(routing=ngc.Pinned(program))


def test_pinned_records():
    g = diamond()
    fg = ngc.FlowLedger(g)
    r1, r2 = dag_of(g, ids(g, [0, 1])), dag_of(g, ids(g, [2, 3]))
    el = g.eligibility(excluded_edges=ids(g, [3]))
    # The replay fails the second route (its last edge excluded): the head's share goes to the
    # first; one record per program entry, the failed one with an empty bundle.
    pol = ngc.FlowPolicy(fg, 0, 3, pinned_over(g, r1, r2), el)
    r = pol.place(9)
    assert (r.placed, r.termination) == (5, T.EXHAUSTED)
    rec, failed = pol.records()
    assert rec.bundle == r1 and pol.record(rec.flow_id).delivered == 5
    assert (rec.entry, failed.entry, failed.route_cost) == (0, 1, None)
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
    assert (r.placed, r.termination, r.by_source[0].admitted) == (0, T.NO_ROUTE, 0)


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
    opts = ngc.QueryOptions(GLOBAL, outputs=ngc.OutputOptions(edge_flows=True))
    r = ngc.capacity(g, 0, 3, opts)
    assert (r.placed, r.termination, r.cost_distribution) == (8, T.EXHAUSTED, [(2, 8)])
    assert int(np.sum(r.edge_flows)) == 16
    common = g.eligibility(excluded_edges=ids(g, [3]))
    pairs = [(0, 3), ([0], [3]), (1, 3), (0, [1, 2])]
    per = [None, g.eligibility(excluded_edges=ids(g, [0])), None, None]
    want = [5, 0, 5, 8]
    for workers in (1, 3, None):
        got = ngc.batch_capacity(
            g, pairs, ngc.QueryOptions(GLOBAL, elig=common), per, workers
        )
        assert [x.placed for x in got] == want
    # workers=None reads NGRAPH_CORE_BATCH_THREADS. The count changes only how many threads run,
    # never a result, so this checks only that the default path with the variable set gives the
    # same results (its parsing: the C++ test Evaluators.DefaultWorkersParseTheEnvironment).
    monkeypatch.setenv("NGRAPH_CORE_BATCH_THREADS", "4")
    got = ngc.batch_capacity(g, pairs, ngc.QueryOptions(GLOBAL, elig=common), per)
    assert [x.placed for x in got] == want
    assert ngc.batch_capacity(g, [], ngc.QueryOptions(GLOBAL)) == []
    with pytest.raises(ValueError, match="node 7"):
        ngc.batch_capacity(
            g, [(0, 3), (0, 7), (0, 8)], ngc.QueryOptions(GLOBAL), workers=4
        )
    with pytest.raises(
        ValueError, match=exact("per_pair_elig has 1 entries for 4 pairs")
    ):
        ngc.batch_capacity(g, pairs, ngc.QueryOptions(GLOBAL), [None])
    with pytest.raises(
        ValueError, match=exact("batch_capacity does not accept Pinned")
    ):
        ngc.batch_capacity(
            g,
            pairs,
            ngc.QueryOptions(pinned_over(g, dag_of(g, ids(g, [0, 1])))),
        )
    cache = ngc.DagCache(g)
    eb = ngc.QueryOptions(COSTONLY, cache=cache)
    with pytest.raises(ValueError, match="workers == 1"):
        ngc.batch_capacity(g, pairs, eb, workers=2)
    got = ngc.batch_capacity(
        g, pairs, eb
    )  # None with a cache: 1 worker, even with the env set
    assert [x.placed for x in got] == [7, 7, 5, 7]  # EB: [4, 3] over links of 5 and 3
    with pytest.raises(ValueError, match="workers"):
        ngc.batch_capacity(g, pairs, ngc.QueryOptions(GLOBAL), workers=0)


def test_sensitivity_candidates_and_signed_losses(monkeypatch):
    # Global: an unsaturated critical edge appears.
    g = build(3, [(0, 1, 9, 1), (1, 2, 4, 1), (0, 2, 2, 5)])
    s = ngc.sensitivity(g, 0, 2, ngc.QueryOptions(GLOBAL))
    assert s.baseline.placed == 6
    pos = g.input_index()
    assert {int(pos[e.edge]): e.loss for e in s.entries} == {0: 4, 1: 4, 2: 2}
    # workers=None reads NGRAPH_CORE_SENSITIVITY_THREADS: the same results on 3 workers (only that
    # the default path with the variable set works; its parsing is tested in C++).
    monkeypatch.setenv("NGRAPH_CORE_SENSITIVITY_THREADS", "3")
    assert plain(ngc.sensitivity(g, 0, 2, ngc.QueryOptions(GLOBAL))) == plain(s)
    # A policy model: every allowed edge is a candidate, a zero-capacity one too. Equal splitting
    # over parallel links of 1, 10 and 0 fits 2 units (a third would cross the zero link);
    # without the thin link 1 fits, without the wide link 1, without the zero link 2.
    p = build(2, [(0, 1, 1, 1), (0, 1, 10, 1), (0, 1, 0, 1)])
    opts = ngc.QueryOptions(COSTONLY)
    s = ngc.sensitivity(p, 0, 1, opts, workers=2)
    assert s.baseline.placed == 2
    pos = p.input_index()
    assert {int(pos[e.edge]): e.loss for e in s.entries} == {0: 1, 1: 1, 2: 0}
    # A negative loss on a zero-capacity link (the search graph S): the only shortest route 0 -> 1 has capacity
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
        s = ngc.sensitivity(g, 0, 2, ngc.QueryOptions(cfg, outputs=ALL_OUT))
        want = ngc.capacity(g, 0, 2, ngc.QueryOptions(cfg, outputs=ALL_OUT))
        assert want.edge_flows is not None and want.residual is not None
        assert plain(s.baseline) == plain(want)


def test_evaluator_options_and_results():
    # The options hold what they were given (None: everything allowed, no supplied cache, no
    # output); the results have no public constructor.
    g = diamond()
    el = g.eligibility(excluded_edges=ids(g, [0]))
    cache = ngc.DagCache(g)
    for cls in (ngc.QueryOptions, ngc.QueryOptions):
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
    bopts = ngc.QueryOptions(COSTONLY, cache=building)
    sopts = ngc.QueryOptions(COSTONLY, cache=building)
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
        ngc.sensitivity(g, 0, 3, ngc.QueryOptions(GLOBAL), workers=0)


def test_python_keeps_what_a_policy_and_an_option_need_alive():
    # Once the caller drops them, the evaluators' options still apply their Eligibility and cache,
    # the ledger still reads its graph, and a policy still places on its ledger. Graphs of the
    # same shape are built in between, so the memory of a released value would be reused.
    g = diamond()
    want = [g.capacity_units().tolist(), g.cost().tolist(), g.edge_dst().tolist()]
    ledger = ngc.FlowLedger(g)
    policy = ngc.FlowPolicy(ledger, 0, 3, COSTONLY)
    opts = ngc.QueryOptions(
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
            match=exact("offer must be an integer in [0, 2^62) or permitted +inf"),
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
        (ngc.capacity, (g, 0, 3, ngc.QueryOptions(GLOBAL))),
        (ngc.batch_capacity, (g, [(0, 3), (1, 3)], ngc.QueryOptions(GLOBAL))),
        (ngc.sensitivity, (g, 0, 3, ngc.QueryOptions(GLOBAL))),
    ]
    for fn, args in calls:
        # An evaluator keeps no state between calls: nothing to snapshot.
        failures = sweep_python_allocation(lambda fn=fn, args=args: (fn, args, list))
        assert failures >= 3, fn.__name__  # not vacuous: every result object allocates


def test_batch_reports_the_lowest_failing_pair():
    # The arguments first, then each pair completely in pair order, across the
    # Python/C++ boundary.
    g = build(2, [(0, 1, 3, 7)])
    opts = ngc.QueryOptions(GLOBAL)
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
    opts = ngc.QueryOptions(GLOBAL)
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
    # reports it (None before the epoch's preparation and after remove()), changes nothing and
    # claims no cut, even when asked for one.
    g = build(2, [(0, 1, 3, 7)])
    ledger = ngc.FlowLedger(g)
    p = ngc.FlowPolicy(ledger, 0, 1, GLOBAL)
    assert p.place(0).cost_distribution is None
    assert p.place(2).cost_distribution == [(7, 2)]
    z = p.place(0, ALL_OUT)
    assert (z.placed, z.cost_distribution, z.min_cut) == (0, [(7, 2)], None)
    assert p.place(0, ngc.OutputOptions(min_cut=True)).min_cut is None  # a zero call
    assert p.placed_demand() == 2 and [int(x) for x in ledger.allocated()] == [2]
    p.remove()
    assert p.place(0).cost_distribution is None
