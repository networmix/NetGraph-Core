"""The per-edge amount arrays at the binding: a ledger's baseline (at
construction and reset), a search's residual (the public search and its K shortest paths) and
the flow core's view (the private _place_global_units) share one validation of the unnarrowed Python integers: the length first,
then every entry in EdgeId order, the first one outside [0, capacity_units[e]] reported with its
exact value. Each user keeps its own names: ``residual`` / ``residual (units)`` and ``view`` /
``view (capacity units)``.
"""

from __future__ import annotations

import _netgraph_core as ext
import numpy as np
import pytest
from helpers import GLOBAL, build, exact, ledger_state, sel

import netgraph_core as ngc


def graph():
    # EdgeIds 0: 0->1 (capacity 5), 1: 1->2 (capacity 3)
    return build(3, [(0, 1, 5, 1), (1, 2, 3, 1)])


# Each user: (name, label, call(graph, amounts, later_defect)). With later_defect the call also
# carries an argument checked after the amounts, which must not be the error reported.
def ledger(g, a, later):
    return ngc.FlowLedger(g, a)


def reset(g, a, later):
    fg = ngc.FlowLedger(g, [4, 2])
    fg.reset(a)
    return fg


def search(g, a, later):
    opts = ngc.SpfOptions(
        ngc.Direction.FORWARD, sel(), residual=a, min_residual=0 if later else 1
    )
    return ngc.spf(g, 0, opts)


def paths(g, a, later):
    elig = graph().eligibility() if later else None  # a foreign Eligibility
    return ngc.ksp(g, 0, 2, ngc.KspOptions(1, sel(), residual=a, elig=elig))


def view(g, a, later):
    return ext._place_global_units(g, a, [9] if later else [0], [2], None, None)


USERS = {
    "ledger": ("residual", "residual (units)", ledger),
    "reset": ("residual", "residual (units)", reset),
    "search": ("residual", "residual (units)", search),
    "paths": ("residual", "residual (units)", paths),
    "view": ("view", "view (capacity units)", view),
}

# (amounts, later defect, the message with {name} and {label}).
CASES = [
    # the length (the shape step) before any entry is read
    ([1.5], False, "{name} has length 1, expected 2 (one entry per edge)"),
    ([1, 1, 1], False, "{name} has length 3, expected 2 (one entry per edge)"),
    (
        np.zeros((1, 2), dtype=np.int64),
        False,
        "{name} must be a 1-D array (got 2 dimensions)",
    ),
    (
        [5, [1]],
        False,
        "{name} must be a 1-D sequence of integers (entry 1 is a sequence)",
    ),
    ("ab", False, "{name} must be a 1-D sequence of integers (got str)"),
    # the entries' types, all of them before any value
    ([1, 2.5], False, "{label}[1] must be an integer (got float)"),
    ([True, 0], False, "{label}[0] must be an integer (got bool)"),
    ([6, 2.5], False, "{label}[1] must be an integer (got float)"),
    (
        np.array([1.0, 2.0]),
        False,
        "{label} entries must be integers (got dtype float64)",
    ),
    # the values: the first entry outside, exact also beyond int64
    ([-1, 99], False, "{name}[0] = -1 is outside [0, capacity_units[0] = 5]"),
    ([5, 4], False, "{name}[1] = 4 is outside [0, capacity_units[1] = 3]"),
    (
        [5, 2**70],
        False,
        "{name}[1] = 1180591620717411303424 is outside [0, capacity_units[1] = 3]",
    ),
    (
        [-(2**70), 9],
        False,
        "{name}[0] = -1180591620717411303424 is outside [0, capacity_units[0] = 5]",
    ),
    (
        np.array([5, 2**64 - 1], dtype=np.uint64),
        False,
        "{name}[1] = 18446744073709551615 is outside [0, capacity_units[1] = 3]",
    ),
    # an amount above its capacity (one above; beyond 64 bits) wins over an argument checked
    # after the amounts
    ([6, 0], True, "{name}[0] = 6 is outside [0, capacity_units[0] = 5]"),
    (
        [5, 2**70],
        True,
        "{name}[1] = 1180591620717411303424 is outside [0, capacity_units[1] = 3]",
    ),
]


# No argument of a ledger is checked after its baseline: its users take no later defect.
PAIRS = [
    (user, case)
    for user in sorted(USERS)
    for case in range(len(CASES))
    if not (CASES[case][1] and user in ("ledger", "reset"))
]


@pytest.mark.parametrize(("user", "case"), PAIRS)
def test_every_user_reports_the_same_first_defect(user, case):
    name, label, call = USERS[user]
    amounts, later, text = CASES[case]
    want = text.format(name=name, label=label)
    with pytest.raises(ValueError, match=exact(want)):
        call(graph(), amounts, later)


@pytest.mark.parametrize(
    ("user", "text"),
    [
        ("paths", "foreign Eligibility: built for another graph"),
        ("search", "min_residual must be at least 1 unit (got 0)"),
        ("view", "sources: node 9 out of node range [0, 3)"),
    ],
)
def test_the_later_defect_alone_is_rejected(user, text):
    # The later defect of the cases above is a defect of its own: with valid amounts the call
    # still fails, so those cases test an order.
    with pytest.raises(ValueError, match=exact(text)):
        USERS[user][2](graph(), [5, 3], True)


# What each user's result shows of the amounts it took, for full amounts (the capacities) and
# for zero ones: the ledger's residual (its baseline), the search's costs (an edge with 0 units
# left is not used: its head stays unreached), the number of paths and the flow placed.
BOUND_OUTCOMES = {
    "ledger": (lambda r: r.residual().tolist(), [5, 3], [0, 0]),
    "reset": (lambda r: r.residual().tolist(), [5, 3], [0, 0]),
    "search": (
        lambda r: r.costs.tolist(),
        [0, 1, 2],
        [0, ngc.UNREACHABLE, ngc.UNREACHABLE],
    ),
    "paths": (len, 1, 0),
    "view": (lambda r: (r["placed"], r["edge_flows"]), (3, [3, 3]), (0, [0, 0])),
}


@pytest.mark.parametrize("user", sorted(USERS))
def test_every_user_accepts_the_bounds(user):
    _name, _label, call = USERS[user]
    shown, full, zero = BOUND_OUTCOMES[user]
    for amounts, want in (
        ([5, 3], full),
        ([0, 0], zero),
        (np.array([5, 3], dtype=np.uint64), full),
    ):
        assert shown(call(graph(), amounts, False)) == want, amounts


def test_a_failed_reset_leaves_the_ledger_unchanged():
    g = graph()
    fg = ngc.FlowLedger(g, [4, 2])
    policy = ngc.FlowPolicy(fg, 0, 2, GLOBAL)
    assert policy.place(2).placed == 2
    del policy  # its records stay as orphans: reset() is allowed and would clear them

    uid, before = fg.uid(), ledger_state(fg)
    assert before[1] == [2, 2] and len(before[4]) == 1
    for bad in ([1], [1, 2.5], [6, 0], [4, 2**70]):
        with pytest.raises(ValueError):
            fg.reset(bad)
        assert (fg.uid(), ledger_state(fg)) == (uid, before)
    fg.reset([1, 1])
    assert (fg.residual().tolist(), fg.allocated().tolist(), fg.records()) == (
        [1, 1],
        [0, 0],
        [],
    )
