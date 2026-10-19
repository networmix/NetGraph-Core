"""The tunnel oracle (tunnel_oracle.py) against hand derivations, two per rule it restates: one
test per rule, each with two cases derived by hand (a T2 test adds a third where the first two
cannot see a misreading of its rule, and says so). These tests import nothing of the library; a
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
SX = [(0, 1, 0, 5), (0, 1, 1, 6)]  # CostStretch(6/5) admits the next path at equality


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


# The requests RR(RR(D, H, h), N, j) in round order
def test_round_order_requests():
    # RR(10, 2) = (5, 5), RR(5, 2) = (3, 2): round order 3, 3, 2, 2.
    assert requests(orc.setup(twoh(10, 10), [0, 1], [2], 10, 2)) == [3, 3, 2, 2]
    # D = 5 over two heads, N = 2: shares (3, 2), requests 2, 1, 1, 1.
    assert requests(orc.setup(twoh(10, 10), [0, 1], [2], 5, 2)) == [2, 1, 1, 1]


# Every attempt on U minus this setup's earlier up tunnels
def test_attempts_see_earlier_tunnels():
    # Two heads sharing a link: the first tunnel of each head fills node 2's exit (6 of 3 + 3).
    assert statuses(orc.setup(SHARED, [0, 1], [3], 12, 2)) == [UP, UP, NB, NB]
    # The last tunnel without bandwidth: three requests of 3 fit the 10-unit link, the fourth
    # does not.
    assert statuses(orc.setup(LASTDOWN, [0], [1], 12, 4)) == [UP, UP, UP, NB]


# The filter U >= max(b, minimum) is inclusive
def test_inclusive_filter():
    # PAR612, D = 12, N = 4: two requests of 3 fill a link of 6.
    assert routes(orc.setup(PAR612, [0], [1], 12, 4)) == [[0], [0], [1], [1]]
    # Minimum 5 filters the 4-unit link: both on e1, holds [0, 2]; minimum 0: both on e0.
    s = orc.setup(PARC, [0], [1], 2, 2, minimum=5)
    assert (routes(s), s.holds) == ([[1], [1]], [0, 2])
    s = orc.setup(PARC, [0], [1], 2, 2)
    assert (routes(s), s.holds) == ([[0], [0]], [2, 0])


# A request of 0 has no filter unless the minimum is positive
def test_zero_requests():
    # An idle tunnel on a zero-capacity link: request 1 on e1 (e0 has 0), request 0 over the
    # zero link e0 (cost 1).
    s = orc.setup(ZL, [0], [1], 1, 2)
    assert (requests(s), routes(s), statuses(s)) == ([1, 0], [[1], [0]], [UP, UP])
    # Minimum 1 filters the zero link for both.
    s = orc.setup(ZL, [0], [1], 1, 2, minimum=1)
    assert (routes(s), s.holds) == ([[1], [1]], [0, 1])


# The filter first, then the collapse of parallel members
def test_filter_before_collapse():
    sel = orc.Selection(parallel_edges=False)
    # e0 (2 < 6) is filtered before the collapse would keep it (lowest EdgeId).
    assert routes(orc.setup(FILT, [0], [1], 6, 1, sel=sel)) == [[1]]
    # The last tunnel without bandwidth with parallel_edges False: as with True.
    assert statuses(orc.setup(LASTDOWN, [0], [1], 12, 4, sel=sel)) == [UP, UP, UP, NB]


# The single-path order: cost, HIGHER_* bottleneck, hops, EdgeId sequence
def test_single_path_order():
    # FEWER_HOPS: equal costs and hops: the smaller EdgeId, then e0 is filtered.
    assert routes(orc.setup(PARX, [0], [1], 8, 2)) == [[0], [1]]
    # The diamond, design 8: the cheapest path 0-1-2-3 (cost 3, three hops) before the two-hop
    # paths of cost 4.
    assert routes(orc.setup(DIAM5, [0], [3], 8, 2)) == [[0, 2, 4], []]


# The tie scores: HIGHER_AVAILABLE_CAPACITY reads U, HIGHER_HEADROOM U / C, HIGHER_CAPACITY C
def test_tie_score_view():
    # Unreserved [5, 8], service residual the capacities [10, 8]: HIGHER_AVAILABLE_CAPACITY reads
    # U: e1.
    sel = orc.Selection("HIGHER_AVAILABLE_CAPACITY")
    assert routes(orc.setup(PVIEW, [0], [1], 4, 1, sel=sel, unreserved=[5, 8])) == [[1]]
    # HIGHER_HEADROOM: 4/4 = 10/10, a tie: e0, then e0 is filtered: e1.
    sel = orc.Selection("HIGHER_HEADROOM")
    assert routes(orc.setup(PARX, [0], [1], 8, 2, sel=sel)) == [[0], [1]]


# AbsoluteCost(c) is inclusive
def test_absolute_inclusive():
    # X8 with AbsoluteCost(3): head 1's path costs 3: up.
    assert statuses(orc.setup(X8, [0, 1], [2], 4, 1, absolute=3)) == [NB, UP]
    # AbsoluteCost(2): COST_CEILING_EXCEEDED.
    assert statuses(orc.setup(X8, [0, 1], [2], 4, 1, absolute=2)) == [NB, CC]


# CostStretch(f): the exact f times the pool minimum c0
def test_stretch_pool_baseline():
    # X8 with CostStretch(2): c0 = 1 (head 0's zero-capacity link), ceiling 2 < 3.
    s = orc.setup(X8, [0, 1], [2], 4, 1, stretch=2.0)
    assert (s.c0, statuses(s)) == (1, [NB, CC])
    # CostStretch(1.2) on SX: c0 = 5 and 1.2 is exactly 6/5 (the decimal it prints, not its binary
    # value, which times 5 is below 6), so the cost-6 path is admitted at equality; 11/10 is not.
    s = orc.setup(SX, [0], [1], 1, 1, stretch=1.2)
    assert Fraction(1.2) * 5 < 6 == Fraction("1.2") * 5
    assert (s.c0, statuses(s)) == (5, [UP])
    s = orc.setup(SX, [0], [1], 1, 1, stretch=Fraction(11, 10))
    assert (s.c0, statuses(s)) == (5, [CC])


# An up tunnel holds b on every link, a down one nothing
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


# Forwarding: RR over the heads, RR over a head's up tunnels, no renormalization
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


# Lossless admits only totals that give a head without an up tunnel nothing
def test_lossless_dead_head():
    # TWOH(2, 10): head 0 is first and down: 0.
    g = twoh(2, 10)
    s = orc.setup(g, [0, 1], [2], 12, 2)
    assert orc.lossless_bound(g, s.roster, s.tunnels, caps(g)) == 0
    # X8: the same with head 0 behind a zero-capacity link.
    s = orc.setup(X8, [0, 1], [2], 4, 1)
    assert orc.lossless_bound(X8, s.roster, s.tunnels, caps(X8)) == 0


# The replay roster from the program's basis, not from surviving entries
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
    # The TWOH call's program with e1 masked: head 1 leaves, head 0 keeps the whole 12 and loses it;
    # the program has up tunnels none can use: NO_ACTIVE_WEIGHT.
    g = twoh(2, 10)
    s = orc.setup(g, [0, 1], [2], 12, 2)
    roster, st, alive = orc.replay(g, s.roster, s.tunnels, [2], ([], [1]))
    assert (roster, st) == ([0], [NB, FAILED, NB, FAILED])
    delivered, lost, _, _ = orc.lossy(g, roster, s.tunnels, 12, caps(g), alive)
    assert (delivered, lost) == (0, 12)
    assert orc.termination(12, roster, s.tunnels, alive, True, 0) == "NO_ACTIVE_WEIGHT"


# One record per entry, down and failed ones included
def test_records_per_entry():
    # DIAM5 at D = 14: both down, two zero records.
    s = orc.setup(DIAM5, [0], [3], 14, 2)
    assert statuses(s) == [NB, NB]
    assert orc.records(s.tunnels) == [(0, 0, []), (1, 0, [])]
    # The TWOH call's program: four records, the up tunnels holding 3 on e1.
    s = orc.setup(twoh(2, 10), [0, 1], [2], 12, 2)
    assert orc.records(s.tunnels) == [
        (0, 0, []),
        (1, 0, [(1, 3)]),
        (2, 0, []),
        (3, 0, [(1, 3)]),
    ]


# ---------------------------------------------------------------- T2

G1 = [(0, 1, 8192, 1), (0, 2, 40960, 1), (1, 3, 8192, 1), (2, 3, 40960, 1)]
PAR412 = [(0, 1, 4, 1), (0, 1, 12, 1)]
# A reverse arc: the first unit's cheapest walk 0-1-2-3 (cost 3) blocks every second walk unless the
# second unit cancels its middle link: 0-2, 2->1 reversed, 1-3 (5 - 1 + 5).
REVERSE = [(0, 1, 1, 1), (0, 2, 1, 5), (1, 2, 1, 1), (1, 3, 1, 5), (2, 3, 1, 1)]
DEAD100 = [(h, 101, 0, 1) for h in range(100)] + [(100, 101, 2**61, 1)]
ROOT_TIE_A = [
    (0, 2, 1, 1),
    (0, 4, 1, 1),
    (1, 2, 1, 1),
    (1, 5, 1, 1),
    (2, 3, 1, 1),
    (4, 3, 1, 1),
    (5, 3, 1, 1),
]
ROOT_TIE_B = [
    (0, 2, 1, 1),
    (0, 5, 1, 1),
    (1, 2, 1, 1),
    (1, 4, 1, 1),
    (2, 3, 1, 1),
    (4, 3, 1, 1),
    (5, 3, 1, 1),
]
FEWEST_A = [(0, 1, 1, 1), (0, 2, 1, 2), (1, 2, 1, 1)]
FEWEST_B = [(0, 1, 1, 1), (0, 3, 1, 3), (1, 2, 1, 1), (2, 3, 1, 1)]
CUT_A = [(0, 1, 2, 2), (1, 2, 2, 2), (1, 2, 3, 2)]
CUT_B = [
    (0, 2, 1, 1),
    (1, 2, 1, 1),
    (2, 3, 1, 1),
    (2, 4, 1, 1),
    (3, 5, 1, 1),
    (4, 5, 1, 1),
]


def x_of(plan, edges):
    out = [0] * len(edges)
    for t in plan.tunnels:
        for e in t.path:
            out[e] += 1
    return out


# Live heads reach a destination over links with V >= 1
def test_heads_with_capacity():
    # A dead head: head 0's only link is 0: dead.
    assert orc.t2(twoh(0, 10), [0, 1], [2], 2, 10).live == [1]
    # A tiny head: links 1 and 100: both live (U = 0 does not make a head dead).
    assert orc.t2(twoh(1, 100), [0, 1], [2], 16).live == [0, 1]


# U, the largest integer size passing the quota flow of the live heads
def test_max_bandwidth_per_tunnel():
    assert orc.t2(G1, [0], [3], 16, 49152, headroom=True).U == 2925  # G1, design 49152
    assert orc.t2(twoh(3, 8), [0, 1], [2], 2, 6).U == 1  # Heads behind 3 and 8
    # A dead head (U = 5 from head 1 alone): the two above have no dead head, which the
    # misreading "dead heads given quotas" needs.
    assert orc.t2(twoh(0, 10), [0, 1], [2], 2, 10).U == 5


# The fluid share s*, the exact breakpoint V[e] / j
def test_max_fractional_bandwidth_per_tunnel():
    # G1, design 49152
    assert orc.t2(G1, [0], [3], 16, 49152).fluid == Fraction(20480, 7)
    p = orc.t2(twoh(1, 100), [0, 1], [2], 16)  # A tiny head
    assert (p.fluid, p.U, p.insufficient, p.tunnels) == (Fraction(1, 16), 0, True, [])


# The working size u from the design shares of the whole roster, dead heads included
def test_working_size():
    assert orc.t2(twoh(0, 10), [0, 1], [2], 2, 10).u == 3  # A dead head: ceil(5 / 2)
    assert orc.t2(G1, [0], [3], 16, 8192, headroom=True).u == 512  # G1 design 8192


# The minimum route-cost sum J first
def test_route_cost():
    # Both tunnels on a cheap link of 2 (u = 1), J = 2.
    p = orc.t2([(0, 1, 2, 1), (0, 1, 10, 2)], [0], [1], 2, 2, headroom=True)
    assert (x_of(p, [0, 1]), p.J) == ([2, 0], 2)
    # G4 design 61440: 15 branches of 4096, 16 tunnels at 2048: J = 32, branch 1 (the smallest keys) gets 2.
    g4 = [(0, i, 4096, 1) for i in range(1, 16)] + [
        (i, 16, 4096, 1) for i in range(1, 16)
    ]
    p = orc.t2(g4, [0], [16], 16, 61440, headroom=True)
    assert (p.J, x_of(p, g4)[:3], x_of(p, g4)[15:18]) == (32, [2, 1, 1], [2, 1, 1])


# The slot prices floor(2^32 (B - V + j u) / B)
def test_headroom_prices():
    # PAR412, N = 2, D = 2 (u = 1): both on the 12-link, 2^32 / 12 + 2^33 / 12 rounded down.
    p = orc.t2(PAR412, [0], [1], 2, 2, headroom=True)
    assert (x_of(p, PAR412), p.phi) == ([0, 2], 1073741823)
    # G1, design 49152
    assert orc.t2(G1, [0], [3], 16, 49152, headroom=True).phi == 73610035200
    # A preload: ONE10 with 6 units free, N = 2, D = 4: U = 3, u = 2, the slots pay
    # floor(2^32 (10 - 6 + 2) / 10) + floor(2^32 (10 - 6 + 4) / 10) = 2576980377 + 3435973836. The
    # two above have V = B, where the misreading "V in the denominator" gives the same
    # prices (so has G10 head 0, whose preload equals the static capacities).
    p = orc.t2(ONE10, [0], [1], 2, 4, V=[6], headroom=True)
    assert (p.U, p.u, x_of(p, ONE10), p.phi) == (3, 2, [2], 6012954213)


# The canonical augmentation: fewest arcs, then the smallest keys, reverse arcs allowed
def test_canonical_augmentation():
    g4 = [(0, i, 4096, 1) for i in range(1, 16)] + [
        (i, 16, 4096, 1) for i in range(1, 16)
    ]
    p = orc.t2(g4, [0], [16], 16, 61440, headroom=True)  # G4 design 61440
    assert [t.path for t in p.tunnels][:3] == [[0, 15], [0, 15], [1, 16]]
    # G1 without prefer_headroom at D = 8192: the cost tie goes to the smallest keys.
    assert x_of(orc.t2(G1, [0], [3], 16, 8192), G1) == [16, 0, 16, 0]
    # A reverse arc: the second unit cancels the first unit's middle link (the two above never need
    # a reverse arc).
    p = orc.t2(REVERSE, [0], [3], 2)
    assert (p.U, x_of(p, REVERSE), [t.path for t in p.tunnels], p.J) == (
        1,
        [1, 1, 0, 1, 1],
        [[0, 3], [1, 4]],
        12,
    )


# The root keys (1, head NodeId, direction): on a tie at the root the smaller head goes first
def test_root_key_tie():
    # Capacities and costs 1, N = 1: every first walk costs 2 with four arcs and every slot pays
    # 2^32. Head 0 goes first and takes e0 e4 (its smaller link key); head 1 then takes e3 e6 (four
    # arcs), not e2, e0 reversed, e1, e5 (six). Heads in descending order: e2 e4 and e1 e5.
    for headroom in (False, True):
        p = orc.t2(ROOT_TIE_A, [0, 1], [3], 1, headroom=headroom)
        assert [(t.head, t.path) for t in p.tunnels] == [(0, [0, 4]), (1, [3, 6])]
    # The same with node 4 and 5 swapped behind the heads: head 0 e0 e4, head 1 e3 e5
    # (descending: e2 e4 and e1 e6).
    p = orc.t2(ROOT_TIE_B, [0, 1], [3], 1)
    assert [(t.head, t.path) for t in p.tunnels] == [(0, [0, 4]), (1, [3, 5])]


# Fewest arcs before the keys
def test_fewest_arcs_before_keys():
    # e1 (cost 2, one link) ties e0 e2 (1 + 1); three arcs beat four, the keys alone take e0 e2.
    assert [t.path for t in orc.t2(FEWEST_A, [0], [2], 1).tunnels] == [[1]]
    # e1 (cost 3) ties e0 e2 e3 (1 + 1 + 1): three arcs beat five.
    assert [t.path for t in orc.t2(FEWEST_B, [0], [3], 1).tunnels] == [[1]]


# The cut into unit walks, sorted, indices j H + i with holes
def test_program():
    p = orc.t2(twoh(0, 10), [0, 1], [2], 2, 10)  # A dead head
    assert [(t.index, t.head, t.path) for t in p.tunnels] == [(1, 1, [1]), (3, 1, [1])]
    p = orc.t2(twoh(3, 8), [0, 1], [2], 2, 6)  # Heads behind 3 and 8
    assert [(t.index, t.head, t.round, t.path) for t in p.tunnels] == [
        (0, 0, 0, [0]),
        (1, 1, 0, [1]),
        (2, 0, 1, [0]),
        (3, 1, 1, [1]),
    ]
    # The walks follow the smallest positive EdgeId, live heads ascending. Head 1's units took e1
    # (cost 2 against head 0's 4), head 0's e0 e2: x = (2, 2, 2); head 0 cuts first, e0 e1 twice,
    # leaving e2 to head 1 (the largest EdgeId: head 0 e0 e2, head 1 e1).
    p = orc.t2(CUT_A, [0, 1], [2], 2)
    assert [(t.head, t.path) for t in p.tunnels] == [
        (0, [0, 1]),
        (1, [2]),
        (0, [0, 1]),
        (1, [2]),
    ]
    # x all ones through node 2: head 0 e0 e2 e4, head 1 e1 e3 e5 (the largest EdgeId swaps the
    # suffixes).
    p = orc.t2(CUT_B, [0, 1], [5], 1)
    assert [(t.head, t.path) for t in p.tunnels] == [(0, [0, 2, 4]), (1, [1, 3, 5])]


# F, the program's lossless bound, a dead head's share blocking
def test_max_lossless_offer():
    assert orc.t2(twoh(0, 10), [0, 1], [2], 2, 10).F == 0  # A dead head
    assert (
        orc.t2([(0, 1, 3, 1), (0, 1, 2, 1)], [0], [1], 2, 5).F == 5
    )  # The remainder on the first link
    assert (
        orc.t2([(0, 1, 2, 1), (0, 1, 3, 1)], [0], [1], 2, 5).F == 4
    )  # The remainder on the second link


# Up to L: the smallest k reaching the live part of D, else the best score, smaller k
def test_up_to():
    p = orc.t2(G1, [0], [3], 6, 49152, up_to=True, headroom=True)  # G1 up to 6
    assert (p.scores, p.count) == ([40960] * 5 + [49152], 6)
    p = orc.t2(twoh(0, 10), [0, 1], [2], 2, 10, up_to=True)  # Up to, a dead head
    assert (p.scores, p.count, [t.path for t in p.tunnels]) == ([5, 5], 1, [[1]])
    # DEAD100, unlimited: equal scores 2^61, no target: the smaller count (the two above reach
    # their target, which the misreading "the largest k on a tie" never meets).
    p = orc.t2(DEAD100, list(range(101)), [101], 2, None, up_to=True)
    assert (p.scores, p.count, p.F) == ([2**61, 2**61], 1, 0)


# No tunnels is an installed program keeping the roster
def test_no_tunnels_keep_roster():
    # Forwarding without tunnels: a tiny head (U = 0) forwarding 10: 5 lost at each head, COMPLETE.
    p = orc.t2(twoh(1, 100), [0, 1], [2], 16, 200)
    roster = orc.program_roster(p)
    amounts, lost = orc.entry_shares(10, roster, p.tunnels)
    assert (roster, lost, [orc.rr(10, len(roster), i) for i in range(len(roster))]) == (
        [0, 1],
        10,
        [5, 5],
    )
    assert orc.termination(10, roster, p.tunnels, [], True, 0) == "COMPLETE"
    # Every head dead: TWOH(0, 0): no live head, no tunnel; Lossless EXHAUSTED.
    p = orc.t2(twoh(0, 0), [0, 1], [2], 2, 10)
    roster = orc.program_roster(p)
    assert (p.live, p.tunnels, roster) == ([], [], [0, 1])
    assert orc.termination(10, roster, p.tunnels, [], False, 0) == "EXHAUSTED"


# A dead head's T2 program replayed under masks
def test_t2_replay_roster():
    g = twoh(0, 10)
    p = orc.t2(g, [0, 1], [2], 2, 10)
    roster, st, alive = orc.replay(g, orc.program_roster(p), p.tunnels, [2], ([], []))
    delivered, lost, _, _ = orc.lossy(g, roster, p.tunnels, 10, caps(g), alive)
    assert (roster, st, delivered, lost) == ([0, 1], [UP, UP], 5, 5)
    roster, st, alive = orc.replay(g, orc.program_roster(p), p.tunnels, [2], ([], [0]))
    delivered, lost, _, _ = orc.lossy(g, roster, p.tunnels, 10, caps(g), alive)
    assert (roster, st, delivered, lost) == ([1], [UP, UP], 10, 0)
