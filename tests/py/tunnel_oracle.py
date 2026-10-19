"""An independent reference for T1 (DistributedTE(TunnelCount(N))) and T2
(CentralizedTE(tunnel_count=TunnelCount(N))) as docs/semantics.md states them ("T1: tunnels set up
by their heads", "T2: the controller's equal tunnels", "Forwarding programs"), for the property
tests of test_distributed_te.py and test_centralized_te_tunnels.py. It restates those rules in plain
Python and calls nothing of the library: for
T1 the requests in round order, every attempt a brute force over the simple paths (a path stops at
the first destination) with the bandwidth filter, the collapse of parallel members and the
single-path order, the cost ceilings with the pool baseline c0; for T2 the live heads by simple
paths, the quota flow by augmenting paths (Edmonds-Karp), U by bisection over the integer sizes,
s* by trying every candidate fraction V[e] / j, the objectives by enumerating the count vectors of
path multisets, the canonical augmentation with Bellman-Ford, the cut into unit walks, and "up to"
over every count without pruning; for both the forwarding (round robin over the heads, then over a
head's up tunnels), the lossless bound, the replay under a mask and the records per entry.

The Lossy loss rule is program_oracle.lossy (the synchronous rounds from zero, the nested
split at each link): the tunnels' paths here are drawn over small graphs, and the property tests
use only programs whose links the rounds settle on.

A graph is a list of edges (u, v, capacity, cost) in EdgeId order. A mask is (excluded nodes,
excluded edges); an edge is allowed when neither it nor an endpoint is excluded."""

from __future__ import annotations

import itertools
from dataclasses import dataclass, field
from fractions import Fraction

import program_oracle

UP, NB, CC, FAILED, BLACKHOLED = (
    "UP",
    "INSUFFICIENT_BANDWIDTH",
    "COST_CEILING_EXCEEDED",
    "FAILED",
    "BLACKHOLED",
)


def rr(x, k, i):
    """RR(x, k, i) = x // k + [i < x % k]."""
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


def priority(edges, tie_break, U, e, service):
    """The tie rule's edge score on T1's view: HIGHER_AVAILABLE_CAPACITY reads U,
    HIGHER_HEADROOM U / C (0 when C = 0), HIGHER_CAPACITY C; `service` (the service residual)
    is not read."""
    cap = edges[e][2]
    if tie_break == "HIGHER_AVAILABLE_CAPACITY":
        return Fraction(U[e])
    if tie_break == "HIGHER_HEADROOM":
        return Fraction(U[e], cap) if cap > 0 else Fraction(0)
    return Fraction(cap)


@dataclass
class Selection:
    tie_break: str = "FEWER_HOPS"
    parallel_edges: bool = True

    @property
    def higher(self):
        return self.tie_break != "FEWER_HOPS"


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
    planned: int = 0

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
    """The sources whose node is allowed that statically reach an allowed destination."""
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
    if not sel.parallel_edges:
        best = {}
        for e, (u, v, _, c) in enumerate(edges):
            if not filtered[e]:
                continue
            key = (
                c,
                -priority(edges, sel.tie_break, U, e, service) if sel.higher else 0,
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
        bottleneck = min(priority(edges, sel.tie_break, U, e, service) for e in p)
        cost, hops = cost_of(edges, p), len(p)
        return (cost, -bottleneck if sel.higher else 0, hops, tuple(p))

    return min(paths, key=order)


def exact_fraction(f):
    """A CostStretch factor as the library reads it: a Fraction as it is, a float as the decimal it prints
    (Fraction(str(f)): 1.2 is 6/5), an integer n/1."""
    if isinstance(f, Fraction):
        return f
    return Fraction(str(f)) if isinstance(f, float) else Fraction(f)


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
    """T1's setup. unreserved: the reservation baseline R0 (default: the
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
    # CostStretch: c0 is the smallest static shortest cost over the roster (before any filter).
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


# ---------------------------------------------------------------- forwarding


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
    """The lossless bound: the largest total whose no-drop loads fit V with no share lost at a
    head."""
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


# ---------------------------------------------------------------- replay and records


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
    """One record per program entry, in index order: (entry, delivered, holds)."""
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
    """The termination rule: zero offer, empty roster, rule (3) only for a program with an up
    tunnel."""
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


# ---------------------------------------------------------------- T2

SCALE = 2**32


def t2_live(edges, roster, dests, V, mask=None):
    """The live heads: those with a simple path whose every link has V >= 1."""
    ok = allowed(edges, mask)
    xn = set(mask[0]) if mask else set()
    live_dests = [d for d in dests if d not in xn]
    return [
        h
        for h in roster
        if any(
            all(V[e] >= 1 for e in p) for p in simple_paths(edges, ok, h, live_dests)
        )
    ]


def _max_flow(n, arcs, s, t):
    """Edmonds-Karp on (u, v, capacity) arcs; exact integers."""
    cap = {}
    adj = [set() for _ in range(n)]
    for u, v, c in arcs:
        cap[(u, v)] = cap.get((u, v), 0) + c
        cap.setdefault((v, u), 0)
        adj[u].add(v)
        adj[v].add(u)
    flow = 0
    while True:
        prev = {s: None}
        queue = [s]
        for u in queue:
            for v in sorted(adj[u]):
                if v not in prev and cap[(u, v)] > 0:
                    prev[v] = u
                    queue.append(v)
        if t not in prev:
            return flow
        step, v = None, t
        while prev[v] is not None:
            step = cap[(prev[v], v)] if step is None else min(step, cap[(prev[v], v)])
            v = prev[v]
        v = t
        while prev[v] is not None:
            cap[(prev[v], v)] -= step
            cap[(v, prev[v])] += step
            v = prev[v]
        flow += step


def t2_fits(edges, ok, roster, live, dests, k, size, V):
    """Do k tunnels of the common size `size` fit at every live head? The quota flow: link
    capacity min(K, floor(V / size)), an arc of k from a root to each live head, K from each
    destination to a sink; links leaving a destination are not in it."""
    heads = live
    K = len(heads) * k
    n = 1 + max(max(u, v) for u, v, _, _ in edges)
    root, sink = n, n + 1
    arcs = [
        (u, v, min(K, int(Fraction(V[e]) / size)))
        for e, (u, v, _, _) in enumerate(edges)
        if ok[e] and u not in dests
    ]
    arcs += [(root, h, k) for h in heads]
    arcs += [(d, sink, K) for d in dests]
    return _max_flow(n + 2, arcs, root, sink) == K


def t2_common(edges, ok, roster, live, dests, k, V):
    """U: the largest integer size that fits (feasibility decreases in the size: bisection over
    [1, max V]); 0 when none."""
    lo, hi = 0, max(V, default=0) + 1  # lo fits (0 by convention), hi does not
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if t2_fits(edges, ok, roster, live, dests, k, Fraction(mid), V):
            lo = mid
        else:
            hi = mid
    return lo


def t2_fluid(edges, ok, roster, live, dests, k, V):
    """s*: the largest fraction V[e] / j (j <= K) that fits."""
    K = len(live) * k
    candidates = sorted(
        {Fraction(v, j) for v in V if v for j in range(1, K + 1)}, reverse=True
    )
    for c in candidates:
        if t2_fits(edges, ok, roster, live, dests, k, c, V):
            return c
    return Fraction(0)


def t2_working(U, roster, live, k, D):
    """u: for a finite design, min(U, max(1, max over live heads of ceil(RR(D, H, i) / k))),
    the shares over the whole roster; unlimited: U."""
    if D is None:
        return U
    H = len(roster)
    need = max(-(-rr(D, H, roster.index(h)) // k) for h in live)
    return min(U, max(1, need))


def slot_price(B, V, j, u):
    """The headroom price of slot j of an edge."""
    return SCALE * (B - V + j * u) // B


def t2_objectives(edges, x, V, u, headroom):
    J = sum(edges[e][3] * x[e] for e in range(len(edges)))
    phi = (
        sum(
            slot_price(edges[e][2], V[e], j, u)
            for e in range(len(edges))
            if edges[e][2] > 0
            for j in range(1, x[e] + 1)
        )
        if headroom
        else 0
    )
    return J, phi


def t2_best_objectives(edges, ok, live, dests, k, u, V, headroom):
    """The minimum (J, Phi) over every count vector of k simple paths per live head with x[e] <=
    min(K, floor(V[e] / u)), by enumerating the path multisets (small graphs only)."""
    K = len(live) * k
    bound = [min(K, v // u) for v in V]
    choices = []
    for h in live:
        paths = simple_paths(edges, ok, h, dests)
        choices.append(list(itertools.combinations_with_replacement(paths, k)))
    best = None
    for pick in itertools.product(*choices):
        x = [0] * len(edges)
        for multiset in pick:
            for p in multiset:
                for e in p:
                    x[e] += 1
        if any(a > b for a, b in zip(x, bound, strict=True)):
            continue
        value = t2_objectives(edges, x, V, u, headroom)
        best = value if best is None or value < best else best
    return best


def t2_canonical(edges, ok, live, dests, k, u, V, headroom):
    """The canonical unit augmentation: K single units from zero, each on the residual path
    from the root to the sink of minimum (J, Phi) price, then the fewest arcs, then the smallest
    sequence of arc keys: (0, EdgeId, direction) for links, (1, head, direction) for root arcs,
    (2, destination, direction) for sink arcs (forward 0, reverse 1). Bellman-Ford distances."""
    n = 1 + max(max(a, b) for a, b, _, _ in edges)
    root, sink = n, n + 1
    K = len(live) * k
    bound = [min(K, v // u) for v in V]
    x = [0] * len(edges)
    rootf = {h: 0 for h in live}
    sinkf = {d: 0 for d in dests}

    def price(e, j):
        B = edges[e][2]
        return slot_price(B, V[e], j, u) if headroom and B > 0 else 0

    for _ in range(K):
        arcs = []  # (tail, head, price, key, action)
        for e, (a, b, _, c) in enumerate(edges):
            if not ok[e] or a in dests:
                continue
            cj = c
            if x[e] < bound[e]:
                arcs.append((a, b, (cj, price(e, x[e] + 1)), (0, e, 0), ("x", e, 1)))
            if x[e] > 0:
                arcs.append((b, a, (-cj, -price(e, x[e])), (0, e, 1), ("x", e, -1)))
        for h in live:
            if rootf[h] < k:
                arcs.append((root, h, (0, 0), (1, h, 0), ("r", h, 1)))
            if rootf[h] > 0:
                arcs.append((h, root, (0, 0), (1, h, 1), ("r", h, -1)))
        for d in dests:
            if sinkf[d] < K:
                arcs.append((d, sink, (0, 0), (2, d, 0), ("s", d, 1)))
            if sinkf[d] > 0:
                arcs.append((sink, d, (0, 0), (2, d, 1), ("s", d, -1)))
        dist = {root: (0, 0)}
        for _round in range(n + 3):
            changed = False
            for a, b, p, _, _ in arcs:
                if a in dist:
                    cand = (dist[a][0] + p[0], dist[a][1] + p[1])
                    if b not in dist or cand < dist[b]:
                        dist[b], changed = cand, True
            if not changed:
                break
        assert not changed and sink in dist, "a negative cycle or no path"

        on = [
            (a, b, p, key, act)
            for a, b, p, key, act in arcs
            if a in dist
            and b in dist
            and (dist[a][0] + p[0], dist[a][1] + p[1]) == dist[b]
        ]
        hop = {sink: 0}
        frontier = [sink]
        while frontier:
            nxt = []
            for v in frontier:
                for arc in on:
                    if arc[1] == v and arc[0] not in hop:
                        hop[arc[0]] = hop[v] + 1
                        nxt.append(arc[0])
            frontier = nxt
        v = root
        while v != sink:
            arc = min(
                (a for a in on if a[0] == v and hop.get(a[1]) == hop[v] - 1),
                key=lambda a: a[3],
            )
            kind, ident, step = arc[4]
            if kind == "x":
                x[ident] += step
            elif kind == "r":
                rootf[ident] += step
            else:
                sinkf[ident] += step
            v = arc[1]
    return x


def t2_program(edges, roster, live, dests, k, x, u):
    """The tunnels: x cut, live head by live head in ascending NodeId, into unit walks along
    the smallest positive EdgeId to a destination; each head's walks sorted by EdgeId sequence;
    in round order (local index j, then roster position i), index j * H + i."""
    rem = list(x)
    H = len(roster)
    per_head = {}
    for h in live:
        walks = []
        for _ in range(k):
            cur, walk = h, []
            while cur not in dests:
                e = next(
                    e for e, (a, _, _, _) in enumerate(edges) if a == cur and rem[e] > 0
                )
                rem[e] -= 1
                walk.append(e)
                cur = edges[e][1]
            walks.append(walk)
        per_head[h] = sorted(walks)
    assert not any(rem), "the walks do not cover x"
    tunnels = []
    for j in range(k):
        for i, h in enumerate(roster):
            if h not in per_head:
                continue
            path = per_head[h][j]
            index = j * H + i
            t = Tunnel(index, h, i, j, 0, UP, path, cost_of(edges, path), u)
            tunnels.append(t)
    return tunnels


@dataclass
class T2Plan:
    roster: list
    live: list
    count: int
    U: int = 0
    fluid: Fraction = Fraction(0)
    u: int = 0
    x: list = field(default_factory=list)
    tunnels: list = field(default_factory=list)
    J: int = 0
    phi: int = 0
    F: int = 0
    scores: list = field(default_factory=list)

    @property
    def insufficient(self):
        return bool(self.live) and self.U == 0


def t2_bound(edges, plan, V):
    """F: the program's lossless maximum: every roster head's share counts, a dead head's
    included (it blocks every total that gives it a unit)."""
    if not plan.roster:
        return 0
    return lossless_bound(edges, plan.roster, plan.tunnels, V)


def t2_at(edges, ok, roster, live, dests, k, V, D, headroom, U=None):
    """The T2 program at count k for the design D (None: unlimited)."""
    plan = T2Plan(roster, live, k)
    if not live:
        return plan
    plan.U = U if U is not None else t2_common(edges, ok, roster, live, dests, k, V)
    plan.fluid = t2_fluid(edges, ok, roster, live, dests, k, V)
    if plan.U == 0:
        return plan  # no tunnels; program_roster keeps the roster
    plan.u = t2_working(plan.U, roster, live, k, D)
    plan.x = t2_canonical(edges, ok, live, dests, k, plan.u, V, headroom)
    plan.tunnels = t2_program(edges, roster, live, dests, k, plan.x, plan.u)
    plan.J, plan.phi = t2_objectives(edges, plan.x, V, plan.u, headroom)
    return plan


def live_share(roster, live, z):
    """G(z): the live heads' part of the round-robin total z."""
    H = len(roster)
    return sum(rr(z, H, i) for i, h in enumerate(roster) if h in live)


def live_budget(edges, plan, V):
    """F'(k): the largest total whose live shares fit on the program (dead shares
    ignored)."""
    roster = plan.roster

    def fits(z):
        amounts, _ = entry_shares(z, roster, plan.tunnels)
        return all(
            a <= c for a, c in zip(loads(edges, plan.tunnels, amounts), V, strict=True)
        )

    M = len(roster) * plan.count
    lo, hi = 0, M * (plan.U + 1)
    while hi - lo > 1:
        mid = (lo + hi) // 2
        lo, hi = (mid, hi) if fits(mid) else (lo, mid)
    return lo


def t2(
    edges, sources, dests, N, D=None, *, up_to=False, V=None, mask=None, headroom=False
):
    """T2 for TunnelCount(N) or TunnelCount(N, up_to=True) on the free capacity V (default: the
    capacities). Returns the T2Plan with its forwarding bound, and for "up to" the per-count
    scores."""
    V = list(V) if V is not None else [c for _, _, c, _ in edges]
    ok = allowed(edges, mask)
    xn = set(mask[0]) if mask else set()
    live_dests = [d for d in dests if d not in xn]
    roster = roster_of(edges, sources, dests, mask)
    live = t2_live(edges, roster, dests, V, mask)
    scores = []
    if up_to:
        target = live_share(roster, live, D) if D is not None else None
        for k in range(1, N + 1):
            at = t2_at(edges, ok, roster, live, live_dests, k, V, None, headroom)
            budget = live_budget(edges, at, V) if at.tunnels else 0
            z = budget if D is None else min(D, budget)
            scores.append(live_share(roster, live, z) if at.U else 0)
        hit = [
            k for k in range(1, N + 1) if target is not None and scores[k - 1] == target
        ]
        chosen = (
            hit[0] if hit else min(range(1, N + 1), key=lambda k: (-scores[k - 1], k))
        )
    else:
        chosen = N
    plan = t2_at(edges, ok, roster, live, live_dests, chosen, V, D, headroom)
    if not live:
        plan.count = 1 if up_to else N
    plan.F = t2_bound(edges, plan, V)
    plan.scores = scores
    return plan


def program_roster(plan):
    """The roster a T2 program installs: the policy's whole roster, also without tunnels."""
    return plan.roster
