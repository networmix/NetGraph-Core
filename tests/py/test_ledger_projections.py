"""The ledger's two projections through the binding (contract 10.1, 10.3): carried traffic
(allocated, residual) and the strict holds of T1 tunnels (reserved, unreserved), never added;
their baselines, texts, removal, orphans, reset, copy, the policy's reservation output and
records, imported programs and the bottleneck report."""

from __future__ import annotations

import numpy as np
import pytest
from helpers import COSTONLY, build, exact, ledger_state, t1_config

import netgraph_core as ngc

LOSSY = ngc.Admission.LOSSY
HOLDS = ngc.OutputOptions(reservation_edge_flows=True)


def one10():
    return build(2, [(0, 1, 10, 1)])


def twoh(a, b):
    return build(3, [(0, 2, a, 1), (1, 2, b, 1)])


def t1(ledger, sources, dest, n, admission=None):
    return ngc.FlowPolicy(ledger, sources, dest, t1_config(n, admission))


def test_two_projections_of_one_link():
    # Carried 4, held 6 on a link of 10: residual 6, unreserved 4; then replace(2,
    # 6): placed 2, delivery_delta -2.
    fg = ngc.FlowLedger(one10())
    p = t1(fg, 0, 1, 1)
    p.replace(4, design_demand=6)
    assert [
        x.tolist()
        for x in (fg.allocated(), fg.residual(), fg.reserved(), fg.unreserved())
    ] == [
        [4],
        [6],
        [6],
        [4],
    ]
    r = p.replace(2, design_demand=6)
    assert (r.placed, r.delivery_delta) == (2, -2)


def test_holds_are_never_added_to_carried_load():
    # After the T1 call (carried 4, held 6) a CostOnly policy places the residual 6 (not 10 - 4 - 6).
    fg = ngc.FlowLedger(one10())
    t1(fg, 0, 1, 1).replace(4, design_demand=6)
    assert ngc.FlowPolicy(fg, 0, 1, COSTONLY).place(10).placed == 6


def test_unreserved_texts():
    # The reservation baseline is keyword-only, checked after the
    # residual's, with the residual baseline's messages under its own name.
    g = one10()
    fg = ngc.FlowLedger(g, unreserved=[3])
    p = t1(fg, 0, 1, 1)
    r = p.replace(4, design_demand=4)
    assert r.program.tunnels[0].status == ngc.TunnelStatus.NO_BANDWIDTH
    assert (fg.residual().tolist(), fg.unreserved().tolist()) == ([10], [3])
    for bad, text in (
        ([3, 1], "unreserved has length 2, expected 1 (one entry per edge)"),
        ([11], "unreserved[0] = 11 is outside [0, capacity_units[0] = 10]"),
        ([-1], "unreserved[0] = -1 is outside [0, capacity_units[0] = 10]"),
    ):
        with pytest.raises(ValueError, match=exact(text)):
            ngc.FlowLedger(g, unreserved=bad)
    with pytest.raises(
        ValueError,
        match=exact("residual[0] = 11 is outside [0, capacity_units[0] = 10]"),
    ):
        ngc.FlowLedger(g, [11], unreserved=[11])
    with pytest.raises(TypeError):
        ngc.FlowLedger(g, [10], [3])  # unreserved is keyword-only
    # The contract's two-edge cases.
    two = build(2, [(0, 1, 10, 1), (0, 1, 10, 1)])
    with pytest.raises(
        ValueError,
        match=exact("unreserved has length 3, expected 2 (one entry per edge)"),
    ):
        ngc.FlowLedger(two, unreserved=[1, 2, 3])
    with pytest.raises(
        ValueError,
        match=exact("unreserved[1] = 11 is outside [0, capacity_units[1] = 10]"),
    ):
        ngc.FlowLedger(two, unreserved=[5, 11])
    assert ngc.FlowLedger(two, unreserved=[10, 0]).unreserved().tolist() == [10, 0]


def test_snapshots_are_int64_copies():
    fg = ngc.FlowLedger(one10(), unreserved=np.array([7], dtype=np.int64))
    for snap in (fg.reserved, fg.unreserved):
        a = snap()
        assert a.dtype == np.int64
        before = a.tolist()
        a[0] = 99
        assert snap().tolist() == before


def test_remove_releases_both():
    fg = ngc.FlowLedger(one10())
    p = t1(fg, 0, 1, 1)
    p.replace(4, design_demand=6)
    p.remove()
    assert (fg.reserved().tolist(), fg.allocated().tolist()) == ([0], [0])
    assert p.program is None and p.records() == []
    p.remove()
    assert fg.reserved().tolist() == [0]


def test_orphans_reset_and_copy():
    # Holds of a policy dropped without remove() stay as orphans, readable; reset() clears
    # them; reset(unreserved=...) sets the reservation baseline; a copy keeps holds and records.
    g = one10()
    fg = ngc.FlowLedger(g)
    p = t1(fg, 0, 1, 1)
    p.replace(4, design_demand=6)
    with pytest.raises(ValueError, match="cannot be copied"):
        fg.copy()
    with pytest.raises(ValueError, match="rejected while a policy is bound"):
        fg.reset()
    del p
    assert fg.reserved().tolist() == [6]
    assert [r.reserved_edges for r in fg.records()] == [[(0, 6)]]
    copy = fg.copy()
    assert copy.uid() != fg.uid()
    assert ledger_state(copy) == ledger_state(fg)
    fg.reset()
    assert (fg.reserved().tolist(), fg.records()) == ([0], [])
    assert copy.reserved().tolist() == [6]
    fg.reset(unreserved=[3])
    assert (fg.unreserved().tolist(), fg.residual().tolist()) == ([3], [10])
    fg.reset([7])
    assert (fg.unreserved().tolist(), fg.residual().tolist()) == ([10], [7])
    with pytest.raises(
        ValueError,
        match=exact("unreserved[0] = 11 is outside [0, capacity_units[0] = 10]"),
    ):
        fg.reset(unreserved=[11])
    assert fg.residual().tolist() == [7]  # a failed reset changes nothing


def test_reservation_output_and_record_holds():
    # The 6.6 policy's holds per edge [0, 6] and per record (tunnels 1 and 3 hold 3 on e1);
    # a CostOnly policy's reservation output is zeros.
    fg = ngc.FlowLedger(twoh(2, 10))
    p = t1(fg, [0, 1], 2, 2, LOSSY)
    r = p.replace(12, design_demand=12, outputs=HOLDS)
    assert r.reservation_edge_flows.tolist() == [0, 6]
    assert r.reservation_edge_flows.dtype == np.int64
    recs = p.records()
    assert [rec.reserved_edges for rec in recs] == [[], [(1, 3)], [], [(1, 3)]]
    held = [0, 0]
    for rec in recs:
        for e, x in rec.reserved_edges:
            held[e] += x
    assert held == fg.reserved().tolist()
    c = ngc.FlowPolicy(ngc.FlowLedger(twoh(2, 10)), [0, 1], 2, COSTONLY)
    assert c.place(4, HOLDS).reservation_edge_flows.tolist() == [0, 0]
    assert r.residual is None  # only what is requested
    assert p.replace(12, design_demand=12).reservation_edge_flows is None


def test_imported_programs_book_nothing():
    # Pinned over T1's program forwards it without a hold.
    g = twoh(2, 10)
    prog = (
        t1(ngc.FlowLedger(g), [0, 1], 2, 2, LOSSY).replace(12, design_demand=12).program
    )
    fg = ngc.FlowLedger(g)
    pin = ngc.FlowPolicy(
        fg, [0, 1], 2, ngc.FlowPolicyConfig(routing=ngc.Pinned(prog), admission=LOSSY)
    )
    r = pin.replace(12, outputs=HOLDS)
    assert (r.reservation_edge_flows.tolist(), fg.reserved().tolist()) == (
        [0, 0],
        [0, 0],
    )
    assert all(rec.reserved_edges == [] for rec in pin.records())


def test_bottleneck_reads_the_service_view_only():
    # One link held whole (unreserved 0) carries nothing: no saturation is reported.
    fg = ngc.FlowLedger(one10())
    p = t1(fg, 0, 1, 1)
    r = p.replace(0, design_demand=10, outputs=ngc.OutputOptions(bottleneck=True))
    assert fg.unreserved().tolist() == [0]
    assert fg.bottleneck(0) == [] and r.bottleneck == []
