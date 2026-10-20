"""The integer ledger FlowLedger through the bindings (docs/semantics.md#6-the-ledger):
baselines on unnarrowed values, snapshots and records, FlowIds, bottleneck, and the lifetime
table driven by manual placements on one fixed route (a manual flow: CentralizedTE with the route as its
support, or Pinned over a program of the one split route), whose PlacementResult the tests read."""

from __future__ import annotations

import copy
import gc
import itertools

import numpy as np
import pytest
from helpers import (
    COSTONLY,
    GLOBAL,
    INF,
    L62,
    build,
    dag_of,
    exact,
    ids,
    ledger_state,
)
from private_hooks import HOOKS, needs_hooks

import netgraph_core as ngc

T = ngc.Termination


def pinned(fg, src, dst, edges, placement, admission=None):
    """A manual flow on the edges' DAG: OPTIMIZED is CentralizedTE with that DAG as its support (a
    controller restricted to it); EQUAL and WEIGHTED are Pinned over a program of the one
    split route (WEIGHTED's weights: the ledger's residual now, frozen, as the old flow read them
    at its creation)."""
    g = fg.graph()
    admission = admission or ngc.Admission.LOSSLESS
    if placement == ngc.RouteAllocation.OPTIMIZED:
        cfg = ngc.FlowPolicyConfig(
            routing=ngc.CentralizedTE(
                source_allocation=ngc.SourceAllocation.OPTIMIZED,
                support=dag_of(g, edges),
            ),
            admission=admission,
        )
        return ngc.FlowPolicy(fg, src, dst, cfg)
    snapshot = placement == ngc.RouteAllocation.WEIGHTED
    program = ngc.ForwardingProgram.from_routes(
        g,
        src,
        dst,
        [ngc.RouteDefinition(src, dag_of(g, edges))],
        allocation=placement,
        weight_basis=ngc.WeightBasis.RESIDUAL_SNAPSHOT if snapshot else None,
        residual=fg.residual() if snapshot else None,
    )
    cfg = ngc.FlowPolicyConfig(routing=ngc.Pinned(program), admission=admission)
    return ngc.FlowPolicy(fg, src, dst, cfg)


def flow_id(policy):
    recs = policy.records()
    return recs[0].flow_id if recs else None


EB = ngc.RouteAllocation.EQUAL
MF = ngc.RouteAllocation.OPTIMIZED
W = ngc.RouteAllocation.WEIGHTED
LOSSY = ngc.Admission.LOSSY


def test_baseline_validation_order_and_snapshots():
    g = build(2, [(0, 1, 5, 1), (1, 0, 3, 1)])
    fg = ngc.FlowLedger(g)
    assert list(fg.residual()) == [5, 3]
    assert fg.residual().dtype == np.int64
    assert list(fg.carried()) == [0, 0]
    assert list(fg.edge_capacity_units()) == [5, 3]
    assert fg.records() == []
    assert fg.graph() is g
    for bad, msg in (
        (
            [6, 0],
            r"service_baseline\[0\] = 6 is outside \[0, edge_capacity_units\[0\] = 5\]",
        ),
        ([1], "service_baseline has length 1, expected 2"),
        ([-1, 0], r"service_baseline\[0\] = -1"),
        ([L62 * 4, 0], r"service_baseline\[0\] = 18446744073709551616"),
        ([2.5, 0], "integer"),
        ([True, 0], "integer"),
        (np.zeros((1, 2), dtype=np.int64), "1-D"),
    ):
        # An entry or a rank of a wrong type is a TypeError, a length or a value a ValueError.
        with pytest.raises((TypeError, ValueError), match=msg) as error:
            ngc.FlowLedger(g, bad)
        assert error.type is (TypeError if msg in ("integer", "1-D") else ValueError)
    fg2 = ngc.FlowLedger(graph=g, service_baseline=np.array([1, 2], dtype=np.uint64))
    assert list(fg2.residual()) == [1, 2]
    assert fg2.uid() != fg.uid()
    # residual() and carried() are copies: writing one changes neither the ledger nor the next
    # snapshot. edge_capacity_units() is a read-only view of the graph's own capacities.
    for snapshot in (fg.residual, fg.carried):
        a = snapshot()
        assert a.dtype == np.int64
        before = a.tolist()
        a[0] = 7
        assert snapshot().tolist() == before
    assert list(fg.residual()) == [5, 3] and list(fg.carried()) == [0, 0]
    view = fg.edge_capacity_units()
    assert view.dtype == np.int64 and not view.flags.writeable
    assert np.shares_memory(view, g.edge_capacity_units())


def test_the_ledger_has_no_public_placement_or_removal():
    # Only policies place or remove flows: the public members are the readers, reset and copy.
    public = {n for n in dir(ngc.FlowLedger) if not n.startswith("_")}
    assert public == {
        "carried",
        "copy",
        "edge_capacity_units",
        "graph",
        "record",
        "records",
        "reserved",
        "reset",
        "residual",
        "saturated_boundary_edges",
        "uid",
        "unreserved",
    }


def test_lifetime_bound_rows_and_epochs():
    g = build(2, [(0, 1, 5, 1)])
    fg = ngc.FlowLedger(g)
    p = pinned(fg, 0, 1, [0], EB)
    r = p.place(0)
    assert (r.call_delivered, r.call_undelivered, r.termination, r.delivery_delta) == (
        0,
        0,
        T.COMPLETE,
        0,
    )
    assert (r.routes_changed, r.weights_changed, r.call_edge_drops, fg.records()) == (
        False,
        False,
        [],
        [],
    )
    if HOOKS:
        assert fg._state()["next_flow_id"] == 0
    with pytest.raises(ValueError, match="bound"):
        fg.reset()
    with pytest.raises(ValueError, match="bound"):
        fg.copy()
    with pytest.raises(ValueError, match="bound"):
        copy.copy(fg)
    with pytest.raises(ValueError, match="bound"):
        copy.deepcopy(fg)
    r = p.place(2)
    assert (r.call_delivered, r.termination, r.routes_changed, flow_id(p)[1]) == (
        2,
        T.COMPLETE,
        True,
        0,
    )
    r = p.place(INF)
    assert (r.call_delivered, r.termination, r.call_undelivered, r.requested) == (
        3,
        T.EXHAUSTED,
        INF,
        INF,
    )
    assert (r.admitted, r.delivered, r.unadmitted) == (5, 5, None)
    assert (r.routes_changed, r.weights_changed, r.marginal_cost_curve) == (
        False,
        False,
        None,
    )
    for finite in (1, 0):
        with pytest.raises(
            ValueError,
            match=exact(
                "finite placement after infinite placement requires a new epoch"
            ),
        ):
            p.place(finite)
    old = flow_id(p)
    p.remove()
    assert list(fg.carried()) == [0]
    with pytest.raises(ValueError, match="retired or unknown"):
        fg.record(old)
    r = p.place(1)
    assert (r.call_delivered, r.routes_changed, flow_id(p)[1]) == (1, True, 1)
    p.remove()
    p.remove()
    assert p.records() == []
    if HOOKS:
        assert (p._state()["requested"], p._state()["infinite"]) == (0, False)
        assert fg._state() == {"next_flow_id": 2, "bound": [0], "orphans": []}


def test_orphans_from_destruction_reset_and_copy():
    g = build(2, [(0, 1, 5, 1)])
    fg = ngc.FlowLedger(g)
    q = pinned(fg, 0, 1, [0], EB)
    q.place(2)
    del q
    gc.collect()  # destroyed without remove(): an orphan
    if HOOKS:
        assert fg._state()["orphans"] == [0]
    assert fg.records()[0].delivered == 2
    p = pinned(fg, 0, 1, [0], EB)
    assert p.place(5).call_delivered == 3
    with pytest.raises(ValueError, match="bound"):
        fg.copy()
    del p
    gc.collect()  # the second orphan; the ledger has no bound policy left
    c = fg.copy()
    assert c.uid() != fg.uid()
    assert [r.flow_id for r in c.records()] == [(c.uid(), 0), (c.uid(), 1)]
    assert list(c.carried()) == [5]
    d = copy.deepcopy(fg)
    assert list(d.carried()) == [5]
    fg.reset([4])
    assert (fg.records(), list(fg.residual())) == ([], [4])
    if HOOKS:
        assert fg._state() == {"next_flow_id": 2, "bound": [], "orphans": []}
    assert len(c.records()) == 2  # the copy is independent
    fg.reset()  # no baseline: back to the capacities
    assert list(fg.residual()) == [5]


def test_policy_keeps_the_ledger_alive():
    g = build(2, [(0, 1, 5, 1)])
    fg = ngc.FlowLedger(g)
    p = pinned(fg, 0, 1, [0], MF)
    del fg
    gc.collect()
    assert p.place(2).call_delivered == 2
    assert p.delivered() == 2


def test_records_fields_and_counters_beyond_64_bits():
    g = build(4, [(0, 1, 10, 1), (0, 2, 10, 1), (1, 3, 1, 1), (2, 3, 10, 1)])
    fg = ngc.FlowLedger(g)
    p = pinned(fg, 0, 3, range(4), EB, LOSSY)
    r = p.place(10)
    assert (r.call_delivered, r.call_undelivered, r.termination, r.edge_drops) == (
        6,
        4,
        T.COMPLETE,
        [(2, 4)],
    )
    rec = fg.records()[0]
    assert rec.flow_id == (fg.uid(), 0)
    assert (rec.policy_id, rec.sources, rec.destinations) == (0, [0], [3])
    # The route is the ledger's own PredecessorDAG, shared by every snapshot (not a copy).
    assert rec.route is fg.records()[0].route is fg.record(rec.flow_id).route
    assert rec.admission == LOSSY
    assert list(rec.route.via_edges()) == [0, 1, 2, 3]
    assert (rec.min_route_cost, rec.offered, rec.delivered) == (2, 10, 6)
    assert rec.edge_flows == [(0, 5), (1, 5), (2, 1), (3, 5)]
    assert rec.edge_drops == [(2, 4)]
    assert fg.record(rec.flow_id).delivered == 6
    # Records are snapshots: a later call changes the ledger's record, not one read before it.
    held = fg.record(rec.flow_id)
    p.place(1)
    assert (rec.offered, held.offered, fg.record(rec.flow_id).offered) == (10, 10, 11)
    # Lossless records: no offered.
    g2 = build(2, [(0, 1, 1, 1)])
    fg2 = ngc.FlowLedger(g2)
    pinned(fg2, 0, 1, [0], MF).place(1)
    assert fg2.records()[0].offered is None
    # An empty route: min_route_cost None. CentralizedTE routing keeps a record per source; one unit from
    # sources 0 and 1 takes the cheaper edge 0 -> 2, so source 1 has no share.
    g3 = build(3, [(0, 2, 1, 1), (1, 2, 1, 5)])
    fg3 = ngc.FlowLedger(g3)
    assert ngc.FlowPolicy(fg3, [0, 1], 2, GLOBAL).place(1).call_delivered == 1
    rec_used, rec_idle = fg3.records()
    assert (rec_used.sources, rec_used.min_route_cost, rec_used.delivered) == (
        [0],
        1,
        1,
    )
    assert (rec_idle.sources, rec_idle.min_route_cost, rec_idle.delivered) == (
        [1],
        None,
        0,
    )
    assert rec_idle.edge_flows == []
    assert list(rec_idle.route.via_edges()) == []
    # Every unit dropped on a saturated link, nine times: the counters pass 2^64.
    sat = pinned(fg2, 0, 1, [0], EB, LOSSY)
    for _ in range(9):
        sat.place(L62 - 1)
    lossy = fg2.records()[1]
    assert lossy.offered == 9 * (L62 - 1)
    assert lossy.edge_drops == [(0, 9 * (L62 - 1))]
    if HOOKS:
        assert sat._state()["requested"] == 9 * (L62 - 1)
    assert lossy.offered > 1 << 64


def test_record_endpoint_sets_are_ascending():
    # A destination set given out of order and with a repeat: every record lists it ascending,
    # once; a record per source under CentralizedTE, each with its one source.
    g = build(4, [(0, 2, 5, 1), (0, 3, 5, 1), (1, 2, 5, 1), (1, 3, 5, 1)])
    fg = ngc.FlowLedger(g)
    ngc.FlowPolicy(fg, 0, [3, 2, 3], COSTONLY).place(2)
    ngc.FlowPolicy(fg, [1, 0, 1], [3, 2], GLOBAL).place(4)
    records = fg.records()
    assert [r.destinations for r in records] == [[2, 3]] * len(records)
    assert sorted(r.sources for r in records[1:]) == [[0], [1]]
    assert records[0].sources == [0] and len(records) == 3


def test_empty_split_graph_keeps_the_route_cost():
    # Weighted with a residual snapshot of 0 on every route edge under another owner's flow:
    # the record exists at zero delivery with an empty split graph (no positive weight:
    # NO_ACTIVE_WEIGHT); min_route_cost is the route's minimum walk cost, not a split-graph value.
    # The first preparation installs the program: routes and weights changed.
    g = build(2, [(0, 1, 2, 3)])
    fg = ngc.FlowLedger(g)
    pinned(fg, 0, 1, [0], MF).place(2)
    w = pinned(fg, 0, 1, [0], W)
    r = w.place(1)
    assert (r.call_delivered, r.termination, r.routes_changed, r.weights_changed) == (
        0,
        T.NO_ACTIVE_WEIGHT,
        True,
        True,
    )
    rec = fg.record(flow_id(w))
    assert (rec.min_route_cost, rec.delivered, rec.edge_flows) == (3, 0, [])
    if HOOKS:
        # The flow's route is the split route before weight pruning (its replay basis).
        flows = next(iter(w._state()["cohorts"].values()))
        assert flows[0]["route"] == [0]


def test_flow_ids_of_other_ledgers_and_bad_shapes():
    g = build(2, [(0, 1, 5, 1)])
    a, b = ngc.FlowLedger(g), ngc.FlowLedger(g)
    p = pinned(a, 0, 1, [0], MF)
    p.place(1)
    with pytest.raises(
        ValueError,
        match=exact(f"FlowId {flow_id(p)} belongs to another ledger (uid {b.uid()})"),
    ):
        b.record(flow_id(p))
    # Not a tuple is a wrong type, a tuple of another length a wrong value: one text.
    for bad, cls in (
        (0, TypeError),
        ([a.uid(), 0], TypeError),
        ((a.uid(),), ValueError),
    ):
        with pytest.raises(
            (TypeError, ValueError), match=exact("a FlowId is a tuple (ledger uid, id)")
        ) as error:
            a.record(bad)
        assert error.type is cls
    for bad, what in ((True, "ledger uid"), (1.0, "id")):
        fid = (bad, 0) if what == "ledger uid" else (a.uid(), bad)
        with pytest.raises(
            TypeError,
            match=exact(f"FlowId {what} must be an integer (got {type(bad).__name__})"),
        ):
            a.record(fid)
    with pytest.raises(
        ValueError, match=exact(f"FlowId ({a.uid()}, {L62 * 8}) is retired or unknown")
    ):
        a.record((a.uid(), L62 * 8))


@needs_hooks
def test_a_flow_id_read_passes_an_allocation_failure_through():
    # A ledger uid outside [0, 2^64) names no ledger: the read raises "retired or unknown". When a
    # Python allocation fails while that range error is raised, the call raises MemoryError, never
    # the "retired or unknown" answer (which would hide the failure). Both readers of a FlowId.
    import _netgraph_core as ext

    g = build(2, [(0, 1, 5, 1)])
    fg = ngc.FlowLedger(g)
    policy = pinned(fg, 0, 1, [0], MF)
    for uid, read in itertools.product(
        (-1, 2**64, -(2**70)), (fg.record, policy.record)
    ):
        failures = 0
        for k in range(100):
            hit, _, value, error = ext._fail_python_allocation(read, ((uid, 0),), k)
            assert value is None, (k, value)
            if not hit:
                assert isinstance(error, ValueError), (k, repr(error))
                assert str(error) == f"FlowId ({uid}, 0) is retired or unknown"
                break
            failures += 1
            assert type(error) is MemoryError, (k, repr(error))
        else:
            pytest.fail("the read never completed")
        assert failures > 0


@pytest.mark.parametrize("routing", ["Pinned", "ShortestPaths", "CentralizedTE"])
def test_rejected_requests_leave_the_policy_and_the_ledger(routing):
    """Every rejected request (a value that is not an integer unit or +inf, outside [0, 2^62),
    or finite after +inf), on a fresh policy, after a finite call and after +inf, leaves the
    ledger and the policy as they were: the records, allocations and delivered amount always, and
    with the private hooks the epoch state itself (the requested total, the +inf flag and the
    flows). A NumPy integer is an integer unit."""
    g = build(4, [(0, 1, 5, 1), (1, 3, 5, 1), (0, 2, 3, 2), (2, 3, 3, 2)])
    fg = ngc.FlowLedger(g)
    if routing == "Pinned":
        p = pinned(fg, 0, 3, [0, 2], EB)  # EdgeIds of 0 -> 1 -> 3
    else:
        p = ngc.FlowPolicy(fg, 0, 3, GLOBAL if routing == "CentralizedTE" else COSTONLY)

    def state():
        own = p._state() if HOOKS else None
        return ledger_state(fg), p.delivered(), own

    def rejected(requests):
        before = state()
        for value, cls, msg in requests:
            with pytest.raises((TypeError, ValueError), match=exact(msg)) as error:
                p.place(value)
            assert error.type is cls, value
        assert state() == before
        return before

    # One text for every offer outside the domain: a type an offer never takes (a bool, a
    # string) TypeError, a value outside it ValueError.
    offer = "offer must be an integer in [0, 2^62) or permitted +inf"
    bad = [
        (value, ValueError, offer)
        for value in (1.5, -1, L62, -INF, float("nan"), np.float64(2.0))
    ] + [(value, TypeError, offer) for value in (True, "1")]
    rejected(bad)
    assert fg.records() == []
    assert p.place(np.int64(4)).call_delivered == 4
    rejected(bad)
    p.place(INF)
    after_inf = "finite placement after infinite placement requires a new epoch"
    before = rejected([(1, ValueError, after_inf), (0, ValueError, after_inf), *bad])
    if HOOKS:
        assert (before[2]["requested"], before[2]["infinite"]) == (4, True)


def test_ledger_bottleneck_with_eligibility_and_errors():
    g = build(4, [(0, 1, 1, 1), (0, 2, 2, 3), (1, 2, 1, 1), (2, 3, 2, 1)])
    fg = ngc.FlowLedger(g)
    p = pinned(fg, 0, 3, range(4), MF)
    r = p.place(INF)
    assert (r.call_delivered, r.termination) == (2, T.EXHAUSTED)
    assert fg.saturated_boundary_edges(0) == ids(g, [0, 3])
    assert fg.saturated_boundary_edges([0, 0]) == ids(g, [0, 3])
    assert fg.saturated_boundary_edges([0], g.eligibility(excluded_nodes=[0])) == []
    # An excluded edge is neither crossed nor listed: without 0 -> 2 (the one edge with residual)
    # the sources reach only {0}; of its leaving edges 0 -> 1 is listed and 0 -> 2 is not.
    assert fg.saturated_boundary_edges(
        [0], g.eligibility(excluded_edges=ids(g, [1]))
    ) == ids(g, [0])
    with pytest.raises(ValueError, match="empty"):
        fg.saturated_boundary_edges([])
    with pytest.raises(ValueError, match="out of node range"):
        fg.saturated_boundary_edges([L62 * 8])
    with pytest.raises(ValueError, match="foreign"):
        fg.saturated_boundary_edges(
            [0], build(4, [(0, 1, 1, 1)]).eligibility(excluded_edges=[0])
        )


def test_manual_flow_construction_errors():
    g = build(3, [(0, 1, 5, 1), (2, 1, 5, 1)])
    fg = ngc.FlowLedger(g)
    with pytest.raises(
        ValueError, match=exact("routes[0]: no walk from head 0 to a destination")
    ):
        pinned(fg, 0, 2, [0, 1], EB)
    # As a support, a route without a walk leaves the roster empty.
    assert pinned(fg, 0, 2, [0, 1], MF).place(1).termination == T.NO_ROUTE
    with pytest.raises(
        ValueError,
        match=exact("CentralizedTE without tunnel_count requires admission=LOSSLESS"),
    ):
        pinned(fg, 0, 1, [0], MF, LOSSY)
    with pytest.raises(ValueError, match="both endpoint sets"):
        pinned(fg, 1, 1, [0], MF)
    other = build(3, [(0, 1, 5, 1), (2, 1, 5, 1)])
    foreign = ngc.FlowPolicyConfig(
        routing=ngc.CentralizedTE(
            source_allocation=ngc.SourceAllocation.OPTIMIZED,
            support=dag_of(other, [0]),
        )
    )
    with pytest.raises(ValueError, match=exact("support belongs to another graph")):
        ngc.FlowPolicy(fg, 0, 1, foreign)
    foreign_route = ngc.ForwardingProgram.from_routes(
        other, 0, 1, [ngc.RouteDefinition(0, dag_of(other, [0]))]
    )
    with pytest.raises(ValueError, match=exact("program belongs to another graph")):
        ngc.FlowPolicy(
            fg, 0, 1, ngc.FlowPolicyConfig(routing=ngc.Pinned(foreign_route))
        )
    fg.reset()  # no policy is bound
