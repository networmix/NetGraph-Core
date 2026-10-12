"""Independent checks of the flow-core payload (``_place_global_units``) for the Python tests.

Nothing here calls the extension's solver: an exhaustive minimum cost of every integer
value on small graphs, an independent successive-shortest-path solver (Bellman-Ford) for the
others, and the certificate in its two ordered branches (the payload shape, feasibility,
conservation, attained value, the full residual closure of the real sources and its cut,
potentials with the NEVER_REACHED sentinel tested first (never below it, never on a closure
vertex, and exactly on the vertices that no allowed source reaches over eligible edges of
positive view at zero flow), the virtual-terminal convention,
pi(t*) equal to the last curve cost, the curve shape, total and moment, and the complete
curve: against the exhaustive optimum of every value when the graph has at most 7 edges and
every view is at most 3, else against the independent solver at every breakpoint of either
curve; the C++ certificate checks the same clauses and makes the same two comparisons, its
branch chosen by each test).
"""

from __future__ import annotations

from collections import deque
from fractions import Fraction
from itertools import product

NEVER_REACHED = -1


def exhaustive(num_nodes, src, dst, cost, view, sources, dests, allowed=None):
    """Minimum cost of every integer flow value over the allowed edges (id coordinates)."""
    edges = [e for e in range(len(src)) if allowed is None or allowed[e]]
    s_set, d_set = set(sources), set(dests)
    best: dict[int, int] = {}
    for flow in product(*(range(view[e] + 1) for e in edges)):
        bal = [0] * num_nodes
        for e, x in zip(edges, flow, strict=True):
            bal[src[e]] -= x
            bal[dst[e]] += x
        ok = all(
            (bal[v] <= 0)
            if v in s_set
            else (bal[v] >= 0)
            if v in d_set
            else bal[v] == 0
            for v in range(num_nodes)
        )
        if not ok:
            continue
        value = -sum(bal[s] for s in s_set)
        c = sum(cost[e] * x for e, x in zip(edges, flow, strict=True))
        if value not in best or c < best[value]:
            best[value] = c
    return best


def independent_curve(num_nodes, src, dst, cost, view, sources, dests, allowed=None):
    """The minimum-cost curve as breakpoints [(value, cost)], from (0, 0): one shortest
    augmenting path at a time, each found by Bellman-Ford over the residual arcs (a virtual
    source feeds every source, every destination drains into a virtual sink)."""
    big = sum(view) + 1  # more than any flow, so the virtual arcs never bind
    arcs = [
        [src[e], dst[e], view[e], cost[e], 0]
        for e in range(len(src))
        if allowed is None or allowed[e]
    ]
    s_v, t_v = num_nodes, num_nodes + 1
    arcs += [[s_v, s, big, 0, 0] for s in sources]
    arcs += [[d, t_v, big, 0, 0] for d in dests]
    points = [(0, 0)]
    while True:
        dist: list[int | None] = [None] * (num_nodes + 2)
        parent: list[tuple[int, int] | None] = [None] * (num_nodes + 2)
        dist[s_v] = 0
        for _ in range(num_nodes + 2):
            changed = False
            for i, (u, v, cap, w, flow) in enumerate(arcs):
                du, dv = dist[u], dist[v]
                if flow < cap and du is not None and (dv is None or du + w < dv):
                    dist[v], parent[v], changed = du + w, (i, 1), True
                du, dv = dist[u], dist[v]
                if flow > 0 and dv is not None and (du is None or dv - w < du):
                    dist[u], parent[u], changed = dv - w, (i, -1), True
            if not changed:
                break
        if dist[t_v] is None:
            return points
        path, v = [], t_v
        while v != s_v:
            i, sign = parent[v]
            path.append((i, sign))
            v = arcs[i][0] if sign == 1 else arcs[i][1]
        amount = min(
            arcs[i][2] - arcs[i][4] if sign == 1 else arcs[i][4] for i, sign in path
        )
        for i, sign in path:
            arcs[i][4] += sign * amount
        value, total = points[-1]
        points.append((value + amount, total + dist[t_v] * amount))


def points_value(points, x):
    """Cost of the piecewise-linear curve through the breakpoints at a value in its domain
    (an empty or backward segment, which the curve-shape clause rejects, is skipped)."""
    for (a, ca), (b, cb) in zip(points, points[1:], strict=False):
        if a <= x <= b and a < b:
            return ca + Fraction(cb - ca, b - a) * (x - a)
    return Fraction(points[-1][1])


def curve_value(curve, x):
    """Cost of the convex function given by its marginal curve at the integer value x."""
    total = Fraction(0)
    left = x
    for c, u in curve:
        take = min(left, u)
        total += c * take
        left -= take
        if left == 0:
            break
    return total


def certify(
    g, view, sources, dests, cap, payload, node_allowed=None, edge_allowed=None
):
    """Failed clauses of the certificate (empty when the payload is certified)."""
    n, m = g.num_nodes(), g.num_edges()
    src = [int(x) for x in g.edge_src()]
    dst = [int(x) for x in g.edge_dst()]
    cost = [int(x) for x in g.cost()]
    node_ok = [True] * n if node_allowed is None else [bool(x) for x in node_allowed]
    edge_ok = [True] * m if edge_allowed is None else [bool(x) for x in edge_allowed]
    s_all, d_all = sorted(set(sources)), sorted(set(dests))
    f = payload["edge_flows"]
    pi = payload["potentials"]
    placed = payload["placed"]
    curve = [tuple(t) for t in payload["cost_distribution"]]
    fails: list[str] = []
    if len(f) != m or len(pi) != n:
        fails.append("shape")
        return fails
    if cap == 0:
        zero = {
            "placed": 0,
            "cost_distribution": [],
            "edge_flows": [0] * m,
            "min_cut": None,
            "potentials": [NEVER_REACHED] * n,
        }
        if any(payload[k] != v for k, v in zero.items()) or any(
            payload["counters"].values()
        ):
            fails.append("zero cap: the prescribed zero payload")
        return fails
    bal = [0] * n
    for e in range(m):
        if not 0 <= f[e] <= view[e]:
            fails.append(f"feasibility at {e}")
        if not edge_ok[e] and f[e]:
            fails.append(f"flow on ineligible edge {e}")
        bal[src[e]] -= f[e]
        bal[dst[e]] += f[e]
    supply = {s: -bal[s] for s in s_all}
    inflow = {d: bal[d] for d in d_all}
    for v in range(n):
        if v in supply:
            if supply[v] < 0 or (not node_ok[v] and supply[v]):
                fails.append(f"source {v}")
        elif v in inflow:
            if inflow[v] < 0 or (not node_ok[v] and inflow[v]):
                fails.append(f"destination {v}")
        elif bal[v]:
            fails.append(f"conservation at {v}")
    if sum(supply.values()) != placed or sum(inflow.values()) != placed:
        fails.append("value differs from supply or inflow")
    if placed < 0 or (cap is not None and placed > cap):
        fails.append("placed outside [0, cap]")
    s_star, t_star = n, n + 1
    out: list[list[tuple[int, int]]] = [[] for _ in range(n + 2)]
    for e in range(m):
        if not edge_ok[e]:
            continue
        if f[e] < view[e]:
            out[src[e]].append((dst[e], cost[e]))
        if f[e] > 0:
            out[dst[e]].append((src[e], -cost[e]))
    roots = [s for s in s_all if node_ok[s]]
    for s in roots:
        out[s_star].append((s, 0))
        if supply[s] > 0:
            out[s].append((s_star, 0))
    for d in d_all:
        if node_ok[d]:
            out[d].append((t_star, 0))
            if inflow[d] > 0:
                out[t_star].append((d, 0))
    seen = set(roots)
    queue = deque(roots)
    while queue:
        u = queue.popleft()
        for v, _w in out[u]:
            if v not in seen:
                seen.add(v)
                queue.append(v)
    connected = t_star in seen
    if cap is None and connected:
        fails.append("uncapped solve left the destinations connected")
    if cap is not None and placed != cap and connected:
        fails.append("attained value: early stop")
    if connected:
        if payload["min_cut"] is not None:
            fails.append("min_cut on a connected solve")
    else:
        expect = [
            e
            for e in range(m)
            if edge_ok[e] and src[e] in seen and dst[e] not in seen and view[e] > 0
        ]
        if payload["min_cut"] != expect:
            fails.append(f"min_cut {payload['min_cut']} != closure cut {expect}")
        if sum(view[e] for e in expect) != placed:
            fails.append("cut value differs from placed")
    # The vertices some search settled: the first search runs on the zero flow and settles
    # every vertex an allowed source reaches over eligible edges of positive view; a later
    # search adds only the reverses of carrying edges, whose ends that one already settled.
    settled = set(roots)
    queue = deque(roots)
    while queue:
        u = queue.popleft()
        for e in range(m):
            if src[e] == u and edge_ok[e] and view[e] > 0 and dst[e] not in settled:
                settled.add(dst[e])
                queue.append(dst[e])
    for v in range(n):
        if pi[v] < NEVER_REACHED:
            fails.append("potential below the sentinel")
        if v in seen and pi[v] == NEVER_REACHED:
            fails.append(f"closure vertex {v} holds NEVER_REACHED")
        if (pi[v] == NEVER_REACHED) == (v in settled):
            fails.append(
                f"NEVER_REACHED differs from reachability at zero flow at node {v}"
            )
    for s in roots:
        if pi[s] != 0:
            fails.append(f"root {s} potential {pi[s]}")
    for u in range(n):
        if u not in seen or pi[u] == NEVER_REACHED:
            continue
        for v, w in out[u]:
            if v < n and v in seen and pi[v] != NEVER_REACHED and w + pi[u] - pi[v] < 0:
                fails.append(f"negative reduced cost {u}->{v}")
    reached = [d for d in d_all if pi[d] != NEVER_REACHED]
    if reached:
        pt = min(pi[d] for d in reached)
        for d in reached:
            if inflow[d] > 0 and pi[d] != pt:
                fails.append(f"inflow destination {d} potential differs from pi(t*)")
        if curve and curve[-1][0] != pt:
            fails.append("pi(t*) differs from the last curve cost")
    elif placed:
        fails.append("flow delivered but no destination reached")
    if any(u <= 0 for _c, u in curve):
        fails.append("curve volume not positive")
    if any(b[0] <= a[0] for a, b in zip(curve, curve[1:], strict=False)):
        fails.append("curve costs not strictly increasing")
    if sum(u for _c, u in curve) != placed:
        fails.append("curve total differs from placed")
    if sum(c * u for c, u in curve) != sum(cost[e] * f[e] for e in range(m)):
        fails.append("curve moment differs from the flow cost")
    allowed_dests = [d for d in d_all if node_ok[d]]
    if m <= 7 and max(view, default=0) <= 3:
        best = exhaustive(n, src, dst, cost, view, roots, allowed_dests, edge_ok)
        fmax = max(best)
        top = fmax if cap is None else min(cap, fmax)
        if placed != top:
            fails.append(f"placed {placed} != min(cap, F_max) {top} (exhaustive)")
        for x in range(min(top, placed) + 1):
            if curve_value(curve, x) != best[x]:
                fails.append(f"curve differs from the exhaustive optimum at {x}")
                break
    else:
        ref = independent_curve(n, src, dst, cost, view, roots, allowed_dests, edge_ok)
        fmax = ref[-1][0]
        top = fmax if cap is None else min(cap, fmax)
        if placed != top:
            fails.append(f"placed {placed} != min(cap, F_max) {top} (independent)")
        mine = [(0, 0)]
        for c, u in curve:
            mine.append((mine[-1][0] + u, mine[-1][1] + c * u))
        for x in sorted({x for x, _c in mine + ref}):
            if x <= min(top, placed) and points_value(mine, x) != points_value(ref, x):
                fails.append(f"curve differs from the independent solver at {x}")
                break
    return fails
