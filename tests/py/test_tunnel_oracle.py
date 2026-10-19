"""The tunnel oracle (tunnel_oracle.py) against hand derivations, two per rule it restates: one
test per rule, each with two cases derived by hand. These tests import nothing of the library; a
copy of the oracle with one misreading of a rule must fail that rule's test."""

from __future__ import annotations

from fractions import Fraction

import tunnel_oracle as orc

UP, NB, CC, FAILED = orc.UP, orc.NB, orc.CC, orc.FAILED

# The graphs of the hand derivations: (u, v, capacity, cost) in EdgeId order.
ONE10 = [(0, 1, 10, 1)]
PARX = [(0, 1, 4, 1), (0, 1, 10, 1)]
PVIEW = [(0, 1, 10, 1), (0, 1, 8, 1)]
FILT = [(0, 1, 2, 1), (0, 1, 10, 1)]
PARC = [(0, 1, 4, 1), (0, 1, 10, 2)]
PAR612 = [(0, 1, 6, 1), (0, 1, 6, 2)]
ZL = [(0, 1, 0, 1), (0, 1, 10, 2)]
X8 = [(0, 2, 0, 1), (1, 2, 8, 3)]
LASTDOWN = [(0, 1, 2, 1), (0, 1, 10, 2)]  # the last tunnel without bandwidth
SHARED = [(0, 2, 6, 1), (1, 2, 6, 1), (2, 3, 6, 1)]  # two heads sharing a link
DIAM5 = [(0, 1, 6, 1), (0, 2, 6, 3), (1, 2, 4, 1), (1, 3, 6, 3), (2, 3, 6, 1)]
REPLAY = [(0, 2, 10, 1), (0, 3, 10, 5), (1, 2, 10, 1), (3, 2, 10, 5)]
SX = [(0, 1, 0, 5), (0, 1, 1, 6)]  # Stretch(1.2) just below the next path


def twoh(a, b):
    return [(0, 2, a, 1), (1, 2, b, 1)]


def caps(edges):
    return [c for _, _, c, _ in edges]


def statuses(s):
    return [t.status for t in s.tunnels]


def routes(s):
    return [t.path or [] for t in s.tunnels]


def requests(s):
    return [t.requested for t in s.tunnels]


# The requests RR(RR(D, H, h), N, j) in round order (6.2)
def test_round_order_requests():
    # RR(10, 2) = (5, 5), RR(5, 2) = (3, 2): round order 3, 3, 2, 2.
    assert requests(orc.setup(twoh(10, 10), [0, 1], [2], 10, 2)) == [3, 3, 2, 2]
    # D = 5 over two heads, N = 2: shares (3, 2), requests 2, 1, 1, 1.
    assert requests(orc.setup(twoh(10, 10), [0, 1], [2], 5, 2)) == [2, 1, 1, 1]


# Every attempt on U minus this setup's earlier up tunnels (6.4)
def test_attempts_see_earlier_tunnels():
    # Two heads sharing a link: the first tunnel of each head fills node 2's exit (6 of 3 + 3).
    assert statuses(orc.setup(SHARED, [0, 1], [3], 12, 2)) == [UP, UP, NB, NB]
    # The last tunnel without bandwidth: three requests of 3 fit the 10-unit link, the fourth
    # does not.
    assert statuses(orc.setup(LASTDOWN, [0], [1], 12, 4)) == [UP, UP, UP, NB]


# The filter U >= max(b, minimum) is inclusive (6.4)
def test_inclusive_filter():
    # PAR612, D = 12, N = 4: two requests of 3 fill a link of 6.
    assert routes(orc.setup(PAR612, [0], [1], 12, 4)) == [[0], [0], [1], [1]]
    # Minimum 5 filters the 4-unit link: both on e1, holds [0, 2]; minimum 0: both on e0.
    s = orc.setup(PARC, [0], [1], 2, 2, minimum=5)
    assert (routes(s), s.holds) == ([[1], [1]], [0, 2])
    s = orc.setup(PARC, [0], [1], 2, 2)
    assert (routes(s), s.holds) == ([[0], [0]], [2, 0])


# A request of 0 has no filter unless the minimum is positive (6.3)
def test_zero_requests():
    # An idle tunnel on a zero-capacity link: request 1 on e1 (e0 has 0), request 0 over the
    # zero link e0 (cost 1).
    s = orc.setup(ZL, [0], [1], 1, 2)
    assert (requests(s), routes(s), statuses(s)) == ([1, 0], [[1], [0]], [UP, UP])
    # Minimum 1 filters the zero link for both.
    s = orc.setup(ZL, [0], [1], 1, 2, minimum=1)
    assert (routes(s), s.holds) == ([[1], [1]], [0, 1])


# The filter first, then the collapse of parallel members (6.4)
def test_filter_before_collapse():
    sel = orc.Selection(multi_edge=False)
    # e0 (2 < 6) is filtered before the collapse would keep it (lowest EdgeId).
    assert routes(orc.setup(FILT, [0], [1], 6, 1, sel=sel)) == [[1]]
    # The last tunnel without bandwidth with multi_edge False: as with True.
    assert statuses(orc.setup(LASTDOWN, [0], [1], 12, 4, sel=sel)) == [UP, UP, UP, NB]


# The single-path order: cost, PHP bottleneck, hops, EdgeId sequence (3.3)
def test_single_path_order():
    # DETERMINISTIC: equal costs and hops: the smaller EdgeId, then e0 is filtered.
    assert routes(orc.setup(PARX, [0], [1], 8, 2)) == [[0], [1]]
    # The diamond, design 8: the cheapest path 0-1-2-3 (cost 3, three hops) before the two-hop
    # paths of cost 4.
    assert routes(orc.setup(DIAM5, [0], [3], 8, 2)) == [[0, 2, 4], []]


# The priority view: RESIDUAL reads U, HEADROOM U / C, CAPACITY C (6.4)
def test_priority_view():
    # Unreserved [5, 8], service residual the capacities [10, 8]: RESIDUAL reads U: e1.
    sel = orc.Selection(php=True, priority="RESIDUAL")
    assert routes(orc.setup(PVIEW, [0], [1], 4, 1, sel=sel, unreserved=[5, 8])) == [[1]]
    # PHP HEADROOM: 4/4 = 10/10, a tie: e0, then e0 is filtered: e1.
    sel = orc.Selection(php=True, priority="HEADROOM")
    assert routes(orc.setup(PARX, [0], [1], 8, 2, sel=sel)) == [[0], [1]]


# Absolute(c) is inclusive (6.4)
def test_absolute_inclusive():
    # X8 with Absolute(3): head 1's path costs 3: up.
    assert statuses(orc.setup(X8, [0, 1], [2], 4, 1, absolute=3)) == [NB, UP]
    # Absolute(2): COST_CEILING.
    assert statuses(orc.setup(X8, [0, 1], [2], 4, 1, absolute=2)) == [NB, CC]


# Stretch(f): the exact f times the pool minimum c0 (6.4)
def test_stretch_pool_baseline():
    # X8 with Stretch(2): c0 = 1 (head 0's zero-capacity link), ceiling 2 < 3.
    s = orc.setup(X8, [0, 1], [2], 4, 1, stretch=2.0)
    assert (s.c0, statuses(s)) == (1, [NB, CC])
    # Stretch(1.2) on SX: c0 = 5, the stored 1.2 times 5 is below 6.
    s = orc.setup(SX, [0], [1], 1, 1, stretch=1.2)
    assert Fraction(1.2) * 5 < 6
    assert (s.c0, statuses(s)) == (5, [CC])


# An up tunnel holds b on every link, a down one nothing (6.4)
def test_holds():
    # The last tunnel without bandwidth: three up tunnels of 3 on e1, the fourth down.
    s = orc.setup(LASTDOWN, [0], [1], 12, 4)
    assert (s.holds, sum(t.reserved for t in s.tunnels)) == ([0, 9], 9)
    # ONE10, N = 2, D = 12: 6 up, the second 6 down (U = 4).
    s = orc.setup(ONE10, [0], [1], 12, 2)
    assert (statuses(s), s.holds, sum(t.reserved for t in s.tunnels)) == (
        [UP, NB],
        [6],
        6,
    )


# Forwarding: RR over the heads, RR over a head's up tunnels, no renormalization (6.5)
def test_forwarding():
    # TWOH(2, 10), N = 2: head 0 has no up tunnel and loses its 6; head 1 delivers 6.
    g = twoh(2, 10)
    s = orc.setup(g, [0, 1], [2], 12, 2)
    delivered, lost, _, _ = orc.lossy(g, s.roster, s.tunnels, 12, caps(g))
    assert (delivered, lost) == (6, 6)
    # Shares RR(13, 4) = (4, 3, 3, 3); e0 carries H(6, [4, 3]) = (4, 2).
    s = orc.setup(PAR612, [0], [1], 12, 4)
    delivered, lost, per, drops = orc.lossy(
        PAR612, s.roster, s.tunnels, 13, caps(PAR612)
    )
    assert [a for a, _, _ in per] == [4, 3, 3, 3]
    assert ([d for _, d, _ in per], delivered, lost, drops) == (
        [4, 2, 3, 3],
        12,
        0,
        [1, 0],
    )


# Lossless admits only totals that give a head without an up tunnel nothing (3.4, 8.3)
def test_lossless_dead_head():
    # TWOH(2, 10): head 0 is first and down: 0.
    g = twoh(2, 10)
    s = orc.setup(g, [0, 1], [2], 12, 2)
    assert orc.lossless_bound(g, s.roster, s.tunnels, caps(g)) == 0
    # X8: the same with head 0 behind a zero-capacity link.
    s = orc.setup(X8, [0, 1], [2], 4, 1)
    assert orc.lossless_bound(X8, s.roster, s.tunnels, caps(X8)) == 0


# The replay roster from the program's basis, not from surviving entries (8.4)
def test_replay_roster():
    # REPLAY: head 0's tunnel e0 fails, head 0 still reaches 2 over e1, e3: it keeps 5 and
    # loses it; head 1 delivers its 5.
    s = orc.setup(REPLAY, [0, 1], [2], 10, 1)
    assert routes(s) == [[0], [2]]
    roster, st, alive = orc.replay(REPLAY, s.roster, s.tunnels, [2], ([], [0]))
    assert (roster, st) == ([0, 1], [FAILED, UP])
    delivered, lost, _, _ = orc.lossy(
        REPLAY, roster, s.tunnels, 10, caps(REPLAY), alive
    )
    assert (delivered, lost) == (5, 5)
    # The 6.6 program with e1 masked: head 1 leaves, head 0 keeps the whole 12 and loses it;
    # the program has up tunnels none can use: NO_ACTIVE_WEIGHT.
    g = twoh(2, 10)
    s = orc.setup(g, [0, 1], [2], 12, 2)
    roster, st, alive = orc.replay(g, s.roster, s.tunnels, [2], ([], [1]))
    assert (roster, st) == ([0], [NB, FAILED, NB, FAILED])
    delivered, lost, _, _ = orc.lossy(g, roster, s.tunnels, 12, caps(g), alive)
    assert (delivered, lost) == (0, 12)
    assert orc.termination(12, roster, s.tunnels, alive, True, 0) == "NO_ACTIVE_WEIGHT"


# One record per entry, down and failed ones included (10.4)
def test_records_per_entry():
    # DIAM5 at D = 14: both down, two zero records.
    s = orc.setup(DIAM5, [0], [3], 14, 2)
    assert statuses(s) == [NB, NB]
    assert orc.records(s.tunnels) == [(0, 0, []), (1, 0, [])]
    # The 6.6 program: four records, the up tunnels holding 3 on e1.
    s = orc.setup(twoh(2, 10), [0, 1], [2], 12, 2)
    assert orc.records(s.tunnels) == [
        (0, 0, []),
        (1, 0, [(1, 3)]),
        (2, 0, []),
        (3, 0, [(1, 3)]),
    ]
