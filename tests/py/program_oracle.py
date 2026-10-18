"""An independent reference for forwarding programs of tunnels (the property tests of
test_programs.py). It restates the contract's rules in plain Python and calls nothing of the
library: round robin and highest averages (2.3), the injection over the roster and the split of
each head's share over its entries (8.1), the lossless bound of a fixed program (8.3), the loss
rule inside one program with the nested split at each link, sources first and then entries
(9.1, 9.2), and a maximum flow with its residual-graph certificates (4.1).

The loss rule is evaluated by synchronous rounds from zero, not by the library's ordering of
events: for tunnels that visit links in one topological order (the generators draw paths over a
DAG) the rounds settle on the unique solution of the link equations within the longest path's
number of hops."""

from __future__ import annotations

from fractions import Fraction


def rr(n, k):
    """Round robin: n units over k members; the first n % k members get one more."""
    return [n // k + (i < n % k) for i in range(k)] if k else []


def ha(n, weights):
    """Highest averages: each unit to the largest weight / (assigned + 1), ties to the lower
    position; all weights 0 assign nothing."""
    total = sum(weights)
    if total == 0:
        return [0] * len(weights)
    got = [n * w // total for w in weights]
    for _ in range(n - sum(got)):
        best = max(
            range(len(weights)), key=lambda i: (Fraction(weights[i], got[i] + 1), -i)
        )
        got[best] += 1
    return got


def replay_roster(specs, rows, excluded_edges, destinations):
    """The replay roster (8.4): the heads (ascending) that still reach a destination over their
    original routes' edges (the union of their tunnels, the saved basis) without the excluded
    edges, whether or not one of their tunnels survives whole. rows: (u, v, ...) per EdgeId."""
    heads = sorted({h for h, _, _ in specs})
    out = []
    for h in heads:
        basis = {e for head, es, _ in specs if head == h for e in es} - set(
            excluded_edges
        )
        seen, stack = {h}, [h]
        while stack:
            u = stack.pop()
            for e in basis:
                if rows[e][0] == u and rows[e][1] not in seen:
                    seen.add(rows[e][1])
                    stack.append(rows[e][1])
        if seen & set(destinations):
            out.append(h)
    return out


def entry_amounts(x, specs, alive=None, source_weights=None, roster=None):
    """The amount of each tunnel of a program for a total offer x.

    specs: (head, edges, weight) in index order; alive: per tunnel, False when it is FAILED;
    roster: the replay roster (default: every head with a tunnel, the program's own). Equal
    sources split x by round robin over the roster, captured source weights by highest averages;
    a head then splits its share by highest averages over its surviving tunnels' weights. What no
    tunnel carries is lost at the source: an empty roster's whole offer, the offer captured
    weights of 0 do not inject, a head whose surviving tunnels have weight 0 or that has none.
    Returns (amounts, source_loss)."""
    alive = alive if alive is not None else [True] * len(specs)
    roster = roster if roster is not None else sorted({h for h, _, _ in specs})
    if not roster:
        return [0] * len(specs), x
    if source_weights is None:
        shares = rr(x, len(roster))
    else:
        weights = dict(source_weights)
        shares = ha(x, [weights[h] for h in roster])
    amounts = [0] * len(specs)
    lost = x - sum(shares)
    for h, share in zip(roster, shares, strict=True):
        ids = [i for i, (head, _, _) in enumerate(specs) if head == h and alive[i]]
        split = ha(share, [specs[i][2] for i in ids])
        lost += share - sum(split)
        for i, a in zip(ids, split, strict=True):
            amounts[i] = a
    return amounts, lost


def loads(specs, amounts, num_edges):
    """Per edge, the sum of the amounts of the tunnels that use it (once per use)."""
    out = [0] * num_edges
    for (_, edges, _), a in zip(specs, amounts, strict=True):
        for e in edges:
            out[e] += a
    return out


def lossless_bound(
    specs, caps, alive=None, source_weights=None, roster=None, top=1 << 20
):
    """The largest x whose amounts fit every capacity and lose nothing at the sources (the
    split rules are house-monotone, so the fitting offers are a prefix: bisection)."""

    def fits(x):
        amounts, lost = entry_amounts(x, specs, alive, source_weights, roster)
        load = loads(specs, amounts, len(caps))
        return lost == 0 and all(u <= c for u, c in zip(load, caps, strict=True))

    lo, hi = 0, top
    if fits(hi):
        return hi
    while hi - lo > 1:
        mid = (lo + hi) // 2
        lo, hi = (mid, hi) if fits(mid) else (lo, mid)
    return lo


def nested_split(service, arrivals, heads):
    """One link's service split over its occurrences: highest averages over the heads' arrival
    totals (ascending head), then within each head over its occurrences' arrivals (in order)."""
    out = [0] * len(arrivals)
    order = sorted(set(heads))
    totals = [
        sum(a for a, h in zip(arrivals, heads, strict=True) if h == k) for k in order
    ]
    for k, quota in zip(order, ha(service, totals), strict=True):
        ids = [i for i, h in enumerate(heads) if h == k]
        for i, c in zip(ids, ha(quota, [arrivals[i] for i in ids]), strict=True):
            out[i] = c
    return out


def lossy(specs, amounts, caps):
    """The loss rule by synchronous rounds from zero. Returns per tunnel (delivered, dropped)
    and the carried amount of each edge; raises when the rounds do not settle."""
    occ = {}
    for i, (_, edges, _) in enumerate(specs):
        for j, e in enumerate(edges):
            occ.setdefault(e, []).append((i, j))
    carried = [[0] * len(edges) for _, edges, _ in specs]
    limit = 2 + max((len(edges) for _, edges, _ in specs), default=0)
    for _ in range(limit + 1):
        arrival = [
            [a if j == 0 else carried[i][j - 1] for j in range(len(edges))]
            for i, ((_, edges, _), a) in enumerate(zip(specs, amounts, strict=True))
        ]
        new = [[0] * len(edges) for _, edges, _ in specs]
        for e, members in occ.items():
            arr = [arrival[i][j] for i, j in members]
            heads = [specs[i][0] for i, _ in members]
            for (i, j), c in zip(
                members, nested_split(min(caps[e], sum(arr)), arr, heads), strict=True
            ):
                new[i][j] = c
        if new == carried:
            break
        carried = new
    else:
        raise AssertionError("the rounds did not settle")
    rows = []
    for i, a in enumerate(amounts):
        delivered = carried[i][-1] if carried[i] else a
        rows.append((delivered, a - delivered))
    edge = [0] * len(caps)
    for e, members in occ.items():
        edge[e] = sum(carried[i][j] for i, j in members)
    return rows, edge


def max_flow(num_nodes, arcs, sources, sinks):
    """Edmonds-Karp over arcs (u, v, cap) from a super source over `sources` to a super sink
    over `sinks`. Returns (value, flow per arc)."""
    s, t = num_nodes, num_nodes + 1
    big = sum(c for _, _, c in arcs) + 1
    adj = [[] for _ in range(num_nodes + 2)]
    head, cap = [], []

    def add(u, v, c):
        adj[u].append(len(head))
        head.append(v)
        cap.append(c)
        adj[v].append(len(head))
        head.append(u)
        cap.append(0)

    for u, v, c in arcs:
        add(u, v, c)
    for u in sources:
        add(s, u, big)
    for v in sinks:
        add(v, t, big)
    value = 0
    while True:
        parent = [-1] * (num_nodes + 2)
        parent[s] = -2
        queue = [s]
        for u in queue:
            for a in adj[u]:
                if cap[a] > 0 and parent[head[a]] == -1:
                    parent[head[a]] = a
                    queue.append(head[a])
        if parent[t] == -1:
            break
        push, v = big, t
        while v != s:
            a = parent[v]
            push = min(push, cap[a])
            v = head[a ^ 1]
        v = t
        while v != s:
            a = parent[v]
            cap[a] -= push
            cap[a ^ 1] += push
            v = head[a ^ 1]
        value += push
    return value, [cap[2 * k + 1] for k in range(len(arcs))]


def residual_certificates(num_nodes, arcs, flow, sources, sinks):
    """For a flow on arcs (u, v, cap, cost): (no augmenting path from the sources to the sinks,
    no negative-cost cycle in the residual graph). Both true: a minimum-cost maximum flow."""
    res = []
    for (u, v, c, w), f in zip(arcs, flow, strict=True):
        if f < c:
            res.append((u, v, w))
        if f > 0:
            res.append((v, u, -w))
    reach, stack = set(sources), list(sources)
    while stack:
        u = stack.pop()
        for a, b, _ in res:
            if a == u and b not in reach:
                reach.add(b)
                stack.append(b)
    maximal = not (reach & set(sinks))
    dist = [
        0
    ] * num_nodes  # Bellman-Ford from a virtual root joined to every node at cost 0
    for _ in range(num_nodes):
        changed = False
        for a, b, w in res:
            if dist[a] + w < dist[b]:
                dist[b] = dist[a] + w
                changed = True
        if not changed:
            return maximal, True
    return maximal, False
