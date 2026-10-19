"""Independent references for pooled sources and the source allocation (docs/semantics.md, section 7).

Nothing here calls the extension. ``rr`` and ``ha`` are the round robin and the highest averages
(the latter one unit at a time, with exact cross-products). ``PooledReference`` replays a ShortestPaths
policy with ``EQUAL`` or ``WEIGHTED`` routing: each roster source's static shortest-path
DAG (or canonical path) toward the destinations, their union, one split per router, attribution
to sources member by member (highest averages over what each source still holds there) and, on an
edge that drops, highest averages over what each source brought to it; Lossless admits the
largest total whose round robin forwards with no drop, Lossy forwards hop by hop on the free
capacity, and an infinite Lossy request is the first-hop experiment. ``equal_source_allocation`` is the
equal-source-allocation controller: the largest total ``T`` whose round-robin supplies a simple
augmenting-path maximum flow can carry, and the minimum cost at those supplies by successive
shortest paths (Bellman-Ford). ``optimized_source_allocation`` is the same without supply caps.
"""

from __future__ import annotations

INF = float("inf")


def rr(x, k):
    """Round robin of x units over k members: x // k each, one more for the first x % k."""
    return [x // k + (1 if i < x % k else 0) for i in range(k)]


def ha(x, w):
    """Highest averages: x units one at a time to the member with the largest
    w[i] / (given[i] + 1), ties to the earlier member; a zero weight gets nothing and
    all-zero weights give nothing."""
    out = [0] * len(w)
    for _ in range(x):
        best = None
        for i, wi in enumerate(w):
            if wi <= 0:
                continue
            if best is None or wi * (out[best] + 1) > w[best] * (out[i] + 1):
                best = i
        if best is None:
            break
        out[best] += 1
    return out


def _graph_lists(g):
    src = [int(x) for x in g.edge_sources()]
    dst = [int(x) for x in g.edge_destinations()]
    cost = [int(x) for x in g.edge_costs()]
    cap = [int(x) for x in g.edge_capacity_units()]
    return src, dst, cost, cap


def distances_to(n, src, dst, cost, dests, node_allowed, edge_allowed):
    """Static distance of every node to the nearest allowed destination over allowed edges
    (Bellman-Ford; capacity-blind), None where none is reachable."""
    dist = [None] * n
    for d in dests:
        if node_allowed[d]:
            dist[d] = 0
    for _ in range(n):
        changed = False
        for e in range(len(src)):
            u, v = src[e], dst[e]
            if (
                not edge_allowed[e]
                or not node_allowed[u]
                or not node_allowed[v]
                or dist[v] is None
            ):
                continue
            if u in dests:
                continue  # a path stops at the first destination it meets
            if dist[u] is None or dist[v] + cost[e] < dist[u]:
                dist[u] = dist[v] + cost[e]
                changed = True
        if not changed:
            break
    return dist


def tight_edges(src, dst, cost, dests, node_allowed, edge_allowed, dist):
    """Per node, its tight out-edges in ascending EdgeId (none out of a destination)."""
    out = {}
    for e in range(len(src)):
        u, v = src[e], dst[e]
        if (
            not edge_allowed[e]
            or not node_allowed[u]
            or not node_allowed[v]
            or u in dests
        ):
            continue
        if dist[u] is not None and dist[v] is not None and dist[u] == cost[e] + dist[v]:
            out.setdefault(u, []).append(e)
    return out


def canonical_next(dst, dests, tight, dist):
    """Per node, the next edge of its canonical path under FEWER_HOPS: the lowest cost (tight
    edges), then the fewest hops, then the lexicographically smallest EdgeId sequence."""
    hops = {d: 0 for d in dests}
    order = sorted((u for u in tight), key=lambda u: dist[u])
    for u in order:  # every tight edge leads to a strictly smaller distance
        best = None
        for e in tight[u]:
            h = hops.get(dst[e])
            if h is not None and (best is None or h + 1 < best):
                best = h + 1
        if best is not None:
            hops[u] = best
    nxt = {}
    for u in order:
        if u not in hops:
            continue
        for e in tight[u]:  # ascending EdgeId: the first edge of the optimal hop count
            if hops.get(dst[e]) == hops[u] - 1:
                nxt[u] = e
                break
    return nxt


class PooledReference:
    """A ShortestPaths EQUAL or WEIGHTED policy, replayed call by call."""

    def __init__(
        self,
        g,
        sources,
        dests,
        kernel,  # "EB" or "W"
        lossy,
        residual,  # the ledger's residual before the first call (per EdgeId)
        multipath=True,
        weight_basis="CAPACITY",  # W's weights: "CAPACITY" or "RESIDUAL_SNAPSHOT"
        excluded_nodes=(),
        excluded_edges=(),
    ):
        self.n = g.num_nodes()
        self.src, self.dst, self.cost, self.cap = _graph_lists(g)
        self.m = len(self.src)
        self.node_allowed = [v not in set(excluded_nodes) for v in range(self.n)]
        self.edge_allowed = [e not in set(excluded_edges) for e in range(self.m)]
        self.dests = set(dests)
        self.inputs = sorted(set(sources))
        self.kernel = kernel
        self.lossy = lossy
        self.multipath = multipath
        self.weight_basis = weight_basis
        self.free = list(residual)  # the free capacity (the ledger's residual)
        dist = distances_to(
            self.n,
            self.src,
            self.dst,
            self.cost,
            self.dests,
            self.node_allowed,
            self.edge_allowed,
        )
        self.roster = [
            s for s in self.inputs if self.node_allowed[s] and dist[s] is not None
        ]
        self.dist = dist
        self.prepared = False
        self.requested = 0
        self.infinite = False
        self.offered = 0  # Lossy: the round-robin position
        self.records = {}  # source -> dict(edges, drops, offered, delivered)
        self.source_drops = {s: 0 for s in self.inputs}

    # -- preparation ----------------------------------------------------------------------
    def _prepare(self):
        tight = tight_edges(
            self.src,
            self.dst,
            self.cost,
            self.dests,
            self.node_allowed,
            self.edge_allowed,
            self.dist,
        )
        if not self.multipath:
            nxt = canonical_next(self.dst, self.dests, tight, self.dist)
            tight = {u: [e] for u, e in nxt.items()}
        # Each source's route: the tight edges it reaches.
        self.route = {}
        union_nodes = set()
        for s in self.roster:
            seen, stack, edges = {s}, [s], set()
            while stack:
                u = stack.pop()
                for e in tight.get(u, []):
                    edges.add(e)
                    if self.dst[e] not in seen:
                        seen.add(self.dst[e])
                        stack.append(self.dst[e])
            self.route[s] = edges
            union_nodes |= seen
        if self.kernel == "W":
            weight = [
                (self.cap[e] if self.weight_basis == "CAPACITY" else self.free[e])
                for e in range(self.m)
            ]
            # Usable nodes: a destination, or a node with a member of positive weight to a
            # usable node (processed by increasing distance).
            usable = set(self.dests)
            members = {}
            for u in sorted(union_nodes, key=lambda x: self.dist[x]):
                if u in self.dests:
                    continue
                ms = [
                    e
                    for e in tight.get(u, [])
                    if weight[e] > 0 and self.dst[e] in usable
                ]
                if ms:
                    members[u] = ms
                    usable.add(u)
            self.weight = weight
        else:
            members = {u: list(tight[u]) for u in union_nodes if u in tight}
            self.weight = None
        self.members = members
        self.order = sorted(
            (u for u in union_nodes if u not in self.dests),
            key=lambda x: -self.dist[x],
        )
        self.records = {
            s: {"edges": {}, "drops": {}, "offered": 0, "delivered": 0}
            for s in self.roster
        }
        self.prepared = True

    def _split(self, u, x):
        ms = self.members[u]
        if self.kernel == "EB":
            return rr(x, len(ms))
        return ha(x, [self.weight[e] for e in ms])

    # -- forwarding -----------------------------------------------------------------------
    def _forward(self, inject, view, first_hop=False):
        """One forwarding from the roots' injections. Returns per source: carried per edge,
        dropped per edge, delivered, source loss, offered; with view None, nothing drops and
        the result says whether every share fitted (fits) and no share was memberless."""
        k = len(self.roster)
        index = {s: i for i, s in enumerate(self.roster)}
        held = {u: [0] * k for u in self.order}  # what each source holds at a node
        delivered = [0] * k
        loss = [0] * k
        offered = [0] * k
        carried = [dict() for _ in range(k)]
        dropped = [dict() for _ in range(k)]
        fits = True
        for s in self.roster:
            i = index[s]
            if not first_hop:
                held[s][i] += inject[i]
                offered[i] += inject[i]
        for u in self.order:
            have = held[u]
            x = sum(have)
            own_hop = (
                first_hop and u in index and u in self.members
            )  # the first-hop experiment's own injection, per member
            if u not in self.members:
                if x > 0:  # only a root can hold traffic without a member
                    for i in range(k):
                        loss[i] += have[i]
                    fits = False
                continue
            ms = self.members[u]
            shares = self._split(u, x)
            remaining = list(have)
            for e, share in zip(ms, shares, strict=True):
                per = ha(share, remaining)
                remaining = [a - b for a, b in zip(remaining, per, strict=True)]
                if own_hop:
                    own = view[e] if view is not None else 0
                    per[index[u]] += own
                    offered[index[u]] += own
                # Lossless (no view) carries every share; its caller checks the edge bounds.
                if view is None or sum(per) <= view[e]:
                    take = per
                else:
                    take = ha(view[e], per)
                for i in range(k):
                    if take[i]:
                        carried[i][e] = carried[i].get(e, 0) + take[i]
                    if per[i] - take[i]:
                        dropped[i][e] = dropped[i].get(e, 0) + per[i] - take[i]
                head = self.dst[e]
                for i in range(k):
                    if head in self.dests:
                        delivered[i] += take[i]
                    else:
                        held[head][i] += take[i]
        return carried, dropped, delivered, loss, offered, fits

    def _fits(self, total, view):
        inject = rr(total, len(self.roster))
        carried, _, _, _, _, ok = self._forward(inject, None)
        if not ok:
            return False
        for e in range(self.m):
            if sum(c.get(e, 0) for c in carried) > view[e]:
                return False
        return True

    # -- the call -------------------------------------------------------------------------
    def place(self, offer):
        """One place(offer) call; returns the result fields the tests compare and applies it
        to the free capacity and the records."""
        if offer == 0:
            return self._result(0, {}, 0, offer)
        if not self.roster:
            if offer != INF:
                self.requested += offer
            else:
                self.infinite = True
            # Nothing can enter: a finite Lossy offer is lost before it enters (offered =
            # requested, all of it source loss); every row stays zero.
            lost = offer if self.lossy and offer != INF else 0
            self.offered += lost
            return self._result(0, {}, lost, offer, empty=True)
        if not self.prepared:
            self._prepare()
        k = len(self.roster)
        if offer == INF:
            self.infinite = True
        else:
            self.requested += offer
        if not self.lossy:
            own = [0] * self.m
            admitted = 0
            for r in self.records.values():
                admitted += r["delivered"]
                for e, x in r["edges"].items():
                    own[e] += x
            view = [self.free[e] + own[e] for e in range(self.m)]
            lo = admitted
            hi = admitted + offer if offer != INF else admitted + sum(view) + 1
            while lo < hi:  # the largest total that fits (a prefix)
                mid = (lo + hi + 1) // 2
                if self._fits(mid, view):
                    lo = mid
                else:
                    hi = mid - 1
            carried, _, delivered, _, _, _ = self._forward(rr(lo, k), None)
            for e in range(self.m):
                self.free[e] = view[e] - sum(c.get(e, 0) for c in carried)
            for i, s in enumerate(self.roster):
                self.records[s]["edges"] = dict(carried[i])
                self.records[s]["delivered"] = delivered[i]
            return self._result(lo - admitted, {}, 0, offer)
        if offer == INF:
            carried, dropped, delivered, loss, offered, _ = self._forward(
                [0] * k, list(self.free), first_hop=True
            )
        else:
            before = rr(self.offered, k)
            self.offered += offer
            after = rr(self.offered, k)
            inject = [a - b for a, b in zip(after, before, strict=True)]
            carried, dropped, delivered, loss, offered, _ = self._forward(
                inject, list(self.free)
            )
        call_edge_drops = {}
        for i, s in enumerate(self.roster):
            r = self.records[s]
            for e, x in carried[i].items():
                r["edges"][e] = r["edges"].get(e, 0) + x
                self.free[e] -= x
            for e, x in dropped[i].items():
                r["drops"][e] = r["drops"].get(e, 0) + x
                call_edge_drops[e] = call_edge_drops.get(e, 0) + x
            r["offered"] += offered[i]
            r["delivered"] += delivered[i]
            self.source_drops[s] += loss[i]
        return self._result(sum(delivered), call_edge_drops, sum(loss), offer)

    def _result(self, call_delivered, call_edge_drops, call_loss, offer, empty=False):
        rows = []
        for s in self.inputs:
            r = self.records.get(s)
            if r is None:
                rows.append((s, 0 if self.lossy else None, 0, 0, 0, 0))
                continue
            dropped = sum(r["drops"].values())
            if self.lossy:
                rows.append(
                    (
                        s,
                        r["offered"],
                        r["offered"],
                        r["delivered"],
                        self.source_drops[s],
                        dropped,
                    )
                )
            else:
                rows.append((s, None, r["delivered"], r["delivered"], 0, 0))
        drops = {}
        for r in self.records.values():
            for e, x in r["drops"].items():
                drops[e] = drops.get(e, 0) + x
        delivered = sum(r["delivered"] for r in self.records.values())
        requested = INF if self.infinite else self.requested
        if self.lossy:
            offered = sum(r["offered"] for r in self.records.values())
            source_drops = sum(self.source_drops.values())
            if empty or not self.roster:
                offered, source_drops = self.offered, self.offered
            admitted, unadmitted = offered, 0
        else:
            offered, source_drops, admitted = None, 0, delivered
            unadmitted = None if self.infinite else self.requested - delivered
        if offer == 0:
            term = "COMPLETE"
        elif not self.roster:
            term = "NO_ROUTE"
        elif not self.members:
            term = "NO_ACTIVE_WEIGHT"
        elif self.lossy:
            term = "COMPLETE" if offer != INF else "EXHAUSTED"
        else:
            term = (
                "COMPLETE" if offer != INF and call_delivered == offer else "EXHAUSTED"
            )
        return {
            "call_delivered": call_delivered,
            "termination": term,
            "call_edge_drops": sorted((e, x) for e, x in call_edge_drops.items() if x),
            "call_source_drops": call_loss,
            "requested": requested,
            "offered": offered,
            "admitted": admitted,
            "delivered": delivered,
            "unadmitted": unadmitted,
            "source_drops": source_drops,
            "drops": sorted((e, x) for e, x in drops.items() if x),
            "by_source": rows,
            "records": {
                s: (sorted((e, x) for e, x in r["edges"].items() if x), r["delivered"])
                for s, r in self.records.items()
            },
        }


# -- the controllers ------------------------------------------------------------------------


def _max_flow(n, src, dst, view, allowed, sources, caps, dests):
    """Maximum flow (augmenting paths, BFS) from a virtual root with an arc of capacity caps[i]
    to sources[i] (None: unbounded) to a virtual sink fed by every destination."""
    root, sink = n, n + 1
    arcs = []  # [u, v, cap, flow]
    for e in range(len(src)):
        if allowed[e]:
            arcs.append([src[e], dst[e], view[e], 0])
    big = sum(view) + 1
    for s, c in zip(sources, caps, strict=True):
        arcs.append([root, s, big if c is None else c, 0])
    for d in dests:
        arcs.append([d, sink, big, 0])
    total = 0
    while True:
        parent = {root: None}
        queue = [root]
        for u in queue:
            for i, (a, b, c, f) in enumerate(arcs):
                if a == u and f < c and b not in parent:
                    parent[b] = (i, 1)
                    queue.append(b)
                if b == u and f > 0 and a not in parent:
                    parent[a] = (i, -1)
                    queue.append(a)
        if sink not in parent:
            return total
        path, v = [], sink
        while v != root:
            i, sign = parent[v]
            path.append((i, sign))
            v = arcs[i][0] if sign == 1 else arcs[i][1]
        amount = min(arcs[i][2] - arcs[i][3] if s == 1 else arcs[i][3] for i, s in path)
        for i, s in path:
            arcs[i][3] += s * amount
        total += amount


def _min_cost(n, src, dst, cost, view, allowed, sources, caps, dests, value):
    """The minimum cost of a flow of `value` with root arcs capped at caps (None: unbounded),
    by successive shortest paths with Bellman-Ford; None when value is not reachable."""
    root, sink = n, n + 1
    big = sum(view) + 1
    arcs = []
    for e in range(len(src)):
        if allowed[e]:
            arcs.append([src[e], dst[e], view[e], cost[e], 0])
    for s, c in zip(sources, caps, strict=True):
        arcs.append([root, s, big if c is None else c, 0, 0])
    for d in dests:
        arcs.append([d, sink, big, 0, 0])
    sent, total = 0, 0
    while sent < value:
        dist = [None] * (n + 2)
        parent = [None] * (n + 2)
        dist[root] = 0
        for _ in range(n + 2):
            changed = False
            for i, (u, v, c, w, f) in enumerate(arcs):
                if (
                    f < c
                    and dist[u] is not None
                    and (dist[v] is None or dist[u] + w < dist[v])
                ):
                    dist[v], parent[v], changed = dist[u] + w, (i, 1), True
                if (
                    f > 0
                    and dist[v] is not None
                    and (dist[u] is None or dist[v] - w < dist[u])
                ):
                    dist[u], parent[u], changed = dist[v] - w, (i, -1), True
            if not changed:
                break
        if dist[sink] is None:
            return None
        path, v = [], sink
        while v != root:
            i, s = parent[v]
            path.append((i, s))
            v = arcs[i][0] if s == 1 else arcs[i][1]
        amount = min(arcs[i][2] - arcs[i][4] if s == 1 else arcs[i][4] for i, s in path)
        amount = min(amount, value - sent)
        for i, s in path:
            arcs[i][4] += s * amount
        sent += amount
        total += dist[sink] * amount
    return total


def equal_source_allocation(g, roster, dests, view, allowed, limit=None):
    """(T, cost): the largest total up to `limit` whose round-robin supplies over the roster
    fit (the solve reaches T), and the minimum cost at those supplies."""
    n = g.num_nodes()
    src, dst, cost, _ = _graph_lists(g)
    k = len(roster)
    if k == 0:
        return 0, 0
    top = sum(view) if limit is None else limit
    best = 0
    for t in range(1, top + 1):  # feasible totals form a prefix
        if _max_flow(n, src, dst, view, allowed, roster, rr(t, k), dests) != t:
            break
        best = t
    return best, _min_cost(
        n, src, dst, cost, view, allowed, roster, rr(best, k), dests, best
    )


def optimized_source_allocation(g, roster, dests, view, allowed, limit=None):
    """(T, cost): the maximum flow up to `limit` from the roster and its minimum cost."""
    n = g.num_nodes()
    src, dst, cost, _ = _graph_lists(g)
    if not roster:
        return 0, 0
    top = _max_flow(n, src, dst, view, allowed, roster, [None] * len(roster), dests)
    if limit is not None:
        top = min(top, limit)
    return top, _min_cost(
        n, src, dst, cost, view, allowed, roster, [None] * len(roster), dests, top
    )
