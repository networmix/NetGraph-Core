"""An independent reference for T1 (CapacityAware(Tunnels(N)), contract 6.2 to 6.5, 8.3,
8.4), the property tests of test_tunnels_t1.py. It restates the contract in plain Python and calls
nothing of the library: the requests in round order, every attempt a brute force over the simple
paths (a path stops at the first destination) with the bandwidth filter, the collapse of parallel
members and the single-path order, the cost ceilings with the pool baseline c0, the forwarding
(round robin over the heads, then over a head's up tunnels), the lossless bound, the replay under
a mask and the records per entry.

The Lossy loss rule is program_oracle.lossy (the synchronous rounds from zero, the nested
split at each link): the tunnels' paths here are drawn over small graphs, and the property tests
use only programs whose links the rounds settle on.

A graph is a list of edges (u, v, capacity, cost) in EdgeId order. A mask is (excluded nodes,
excluded edges); an edge is allowed when neither it nor an endpoint is excluded."""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

import program_oracle

UP, NB, CC, FAILED, BLACKHOLED = (
    "UP",
    "NO_BANDWIDTH",
    "COST_CEILING",
    "FAILED",
    "BLACKHOLED",
)


def rr(x, k, i):
    """RR(x, k, i) = x // k + [i < x % k] (2.3)."""
    return x // k + (1 if i < x % k else 0)


def allowed(edges, mask):
    xn, xe = (set(mask[0]), set(mask[1])) if mask else (set(), set())
    return [
        e not in xe and u not in xn and v not in xn
        for e, (u, v, _, _) in enumerate(edges)
    ]


def simple_paths(edges, ok, src, dests):
    """Every simple path (EdgeIds in walk order) from src that ends at the first destination it
    meets, over the edges with ok[e]."""
    dests = set(dests)
    out = []

    def walk(u, seen, path):
        for e, (a, b, _, _) in enumerate(edges):
            if a != u or not ok[e] or b in seen:
                continue
            if b in dests:
                out.append([*path, e])
            else:
                walk(b, seen | {b}, [*path, e])

    walk(src, {src}, [])
    return out


def cost_of(edges, path):
    return sum(edges[e][3] for e in path)


def priority(edges, sel_priority, U, e, service):
    """The priority view of 6.4: RESIDUAL reads U, HEADROOM U / C (0 when C = 0), CAPACITY C;
    `service` (the service residual) is not read."""
    cap = edges[e][2]
    if sel_priority == "RESIDUAL":
        return Fraction(U[e])
    if sel_priority == "HEADROOM":
        return Fraction(U[e], cap) if cap > 0 else Fraction(0)
    return Fraction(cap)


@dataclass
class Selection:
    php: bool = False
    priority: str = "CAPACITY"
    multi_edge: bool = True


@dataclass
class Tunnel:
    index: int
    head: int
    slot: int
    round: int
    requested: int
    status: str
    path: list | None = None
    cost: int | None = None

    @property
    def up(self):
        return self.status == UP

    @property
    def reserved(self):
        return self.requested if self.up else 0


@dataclass
class Setup:
    roster: list
    tunnels: list
    holds: list
    c0: int | None


def roster_of(edges, sources, dests, mask=None):
    """The sources whose node is allowed that statically reach an allowed destination (2.2)."""
    ok = allowed(edges, mask)
    xn = set(mask[0]) if mask else set()
    live = [d for d in dests if d not in xn]
    return [
        s
        for s in sorted(set(sources))
        if s not in xn and simple_paths(edges, ok, s, live)
    ]


def static_cost(edges, ok, src, dests):
    paths = simple_paths(edges, ok, src, dests)
    return min((cost_of(edges, p) for p in paths), default=None)


def choose(edges, U, candidate_ok, sel, src, dests, threshold, service):
    """One attempt's path: the filter, then the collapse, then the single-path order."""
    # The filter U >= max(b, minimum); a threshold of 0 keeps every allowed edge.
    filtered = [ok and U[e] >= threshold for e, ok in enumerate(candidate_ok)]
    kept = list(filtered)
    if not sel.multi_edge:
        best = {}
        for e, (u, v, _, c) in enumerate(edges):
            if not filtered[e]:
                continue
            key = (
                c,
                -priority(edges, sel.priority, U, e, service) if sel.php else 0,
                e,
            )
            if (u, v) not in best or key < best[(u, v)][0]:
                best[(u, v)] = (key, e)
        chosen = {e for _, e in best.values()}
        kept = [filtered[e] and e in chosen for e in range(len(edges))]
    paths = simple_paths(edges, kept, src, dests)
    if not paths:
        return None

    def order(p):
        bottleneck = min(priority(edges, sel.priority, U, e, service) for e in p)
        cost, hops = cost_of(edges, p), len(p)
        return (cost, -bottleneck if sel.php else 0, hops, tuple(p))

    return min(paths, key=order)


def exact_fraction(f):
    return Fraction(f) if not isinstance(f, Fraction) else f


def setup(
    edges,
    sources,
    dests,
    D,
    N,
    *,
    sel=None,
    minimum=0,
    absolute=None,
    stretch=None,
    mask=None,
    unreserved=None,
    held=None,
    service=None,
):
    """T1's setup (6.2 to 6.4). unreserved: the reservation baseline R0 (default: the
    capacities); held: other owners' holds per edge (default 0); U = R0 - held. service: the
    service residual (default: the capacities), which setup never reads (rule 7's misreading
    does)."""
    sel = sel or Selection()
    ok = allowed(edges, mask)
    xn = set(mask[0]) if mask else set()
    live = [d for d in dests if d not in xn]
    roster = roster_of(edges, sources, dests, mask)
    H = len(roster)
    base = list(unreserved) if unreserved is not None else [c for _, _, c, _ in edges]
    U = [r - (held[e] if held else 0) for e, r in enumerate(base)]
    service = list(service) if service is not None else [c for _, _, c, _ in edges]
    # Stretch: c0 is the smallest static shortest cost over the roster (before any filter).
    c0 = min((static_cost(edges, ok, h, live) for h in roster), default=None)
    tunnels = []
    for j in range(N):
        for i, head in enumerate(roster):
            b = rr(rr(D, H, i), N, j)
            threshold = max(b, minimum)
            path = choose(edges, U, ok, sel, head, live, threshold, service)
            t = Tunnel(j * H + i, head, i, j, b, NB)
            if path is not None:
                c = cost_of(edges, path)
                limit_ok = True
                if absolute is not None:
                    limit_ok = c <= absolute
                elif stretch is not None:
                    limit_ok = c <= exact_fraction(stretch) * c0
                if limit_ok:
                    t.status, t.path, t.cost = UP, path, c
                    for e in path:
                        U[e] -= b
                else:
                    t.status = CC
            tunnels.append(t)
    holds = [0] * len(edges)
    for t in tunnels:
        for e in t.path or []:
            holds[e] += t.reserved
    return Setup(roster, tunnels, holds, c0)


# ---------------------------------------------------------------- forwarding (6.5, 8.3)


def entry_shares(x, roster, tunnels, alive=None):
    """The shares of a total x: RR over the roster, then RR over each head's up tunnels in local
    order. Returns (per tunnel amount, lost at the heads)."""
    alive = alive if alive is not None else [t.up for t in tunnels]
    H = len(roster)
    if H == 0:
        return [0] * len(tunnels), x
    amounts = [0] * len(tunnels)
    lost = 0
    for i, h in enumerate(roster):
        share = rr(x, H, i)
        ids = [
            k
            for k, (t, a) in enumerate(zip(tunnels, alive, strict=True))
            if a and t.head == h
        ]
        if not ids:
            lost += share
            continue
        for m, k in enumerate(ids):
            amounts[k] = rr(share, len(ids), m)
    return amounts, lost


def loads(edges, tunnels, amounts):
    out = [0] * len(edges)
    for t, a in zip(tunnels, amounts, strict=True):
        for e in t.path or []:
            out[e] += a
    return out


def lossless_bound(edges, roster, tunnels, V, alive=None, top=None):
    """8.3: the largest total whose no-drop loads fit V with no share lost at a head."""
    top = top if top is not None else sum(V) + 1

    def fits(x):
        amounts, lost = entry_shares(x, roster, tunnels, alive)
        return lost == 0 and all(
            u <= c for u, c in zip(loads(edges, tunnels, amounts), V, strict=True)
        )

    lo, hi = 0, top
    if fits(hi):
        return hi
    while hi - lo > 1:
        mid = (lo + hi) // 2
        lo, hi = (mid, hi) if fits(mid) else (lo, mid)
    return lo


def lossy(edges, roster, tunnels, x, V, alive=None):
    """The Lossy forwarding of x on the free capacity V: (delivered, lost at the heads, per
    tunnel (assigned, delivered, dropped), per edge drops). The loss rule at each link is
    program_oracle.nested_split (highest averages over the heads' arrivals, then within a head over
    its tunnels'), evaluated by synchronous rounds from zero until they settle."""
    amounts, lost = entry_shares(x, roster, tunnels, alive)
    alive = alive if alive is not None else [t.up for t in tunnels]
    live = [k for k, a in enumerate(alive) if a and amounts[k] > 0]
    occ = {}
    for k in live:
        for j, e in enumerate(tunnels[k].path):
            occ.setdefault(e, []).append((k, j))
    carried = {k: [0] * len(tunnels[k].path) for k in live}
    for _ in range(4 + max((len(tunnels[k].path) for k in live), default=0)):
        new = {k: [0] * len(tunnels[k].path) for k in live}
        for e, members in occ.items():
            arrival = [amounts[k] if j == 0 else carried[k][j - 1] for k, j in members]
            heads = [tunnels[k].head for k, _ in members]
            split = program_oracle.nested_split(min(V[e], sum(arrival)), arrival, heads)
            for (k, j), c in zip(members, split, strict=True):
                new[k][j] = c
        if new == carried:
            break
        carried = new
    else:
        raise AssertionError("the rounds did not settle")
    per = [(amounts[k], 0, 0) for k in range(len(tunnels))]
    drops = [0] * len(edges)
    for k in live:
        path = tunnels[k].path
        for j, e in enumerate(path):
            arrival = amounts[k] if j == 0 else carried[k][j - 1]
            drops[e] += arrival - carried[k][j]
        per[k] = (amounts[k], carried[k][-1], amounts[k] - carried[k][-1])
    delivered = sum(d for _, d, _ in per)
    return delivered, lost, per, drops


# ---------------------------------------------------------------- replay (8.4) and records (10.4)


def replay(edges, program_roster, tunnels, dests, mask, producer_mask=None):
    """A Pinned replay of a T1 program under `mask`: the roster rebuilt from the topology basis
    (the program's exclusions and the mask together), a tunnel FAILED when the mask touches it or
    its head left, setup-down tunnels keeping their status. Returns (replay roster, statuses,
    alive)."""
    xn = set(mask[0]) | (set(producer_mask[0]) if producer_mask else set())
    xe = set(mask[1]) | (set(producer_mask[1]) if producer_mask else set())
    both = (sorted(xn), sorted(xe))
    ok = allowed(edges, both)
    live_dests = [d for d in dests if d not in xn]
    roster = [
        h
        for h in program_roster
        if h not in xn and simple_paths(edges, ok, h, live_dests)
    ]
    statuses, alive = [], []
    for t in tunnels:
        if not t.up:
            statuses.append(t.status)
            alive.append(False)
        elif t.head not in roster or not all(ok[e] for e in t.path):
            statuses.append(FAILED)
            alive.append(False)
        else:
            statuses.append(UP)
            alive.append(True)
    return roster, statuses, alive


def records(tunnels, amounts=None):
    """One record per program entry, in index order: (entry, delivered, holds) (10.4)."""
    amounts = amounts or [0] * len(tunnels)
    out = []
    for t, a in zip(tunnels, amounts, strict=True):
        out.append(
            (
                t.index,
                a,
                [(e, t.reserved) for e in sorted(t.path or [])] if t.up else [],
            )
        )
    return out


def termination(offer, roster, tunnels, alive, lossy_admission, admitted):
    """14.2 as amended: zero offer, empty roster, rule (3) only for a program with an up tunnel."""
    if offer == 0:
        return "COMPLETE"
    if not roster:
        return "NO_ROUTE"
    has_path = any(t.up for t in tunnels)
    usable = any(a and t.head in roster for t, a in zip(tunnels, alive, strict=True))
    if has_path and not usable:
        return "NO_ACTIVE_WEIGHT"
    if lossy_admission:
        return "COMPLETE"
    return "COMPLETE" if admitted == offer else "EXHAUSTED"
