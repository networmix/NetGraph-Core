"""Independent checks of the flow core (APPROACH revision 24, section 2.2).

Nothing here imports ``flowcore``: the exhaustive minimum-cost function (tier a),
an independent complete solver (tier b: successive shortest paths by Bellman-Ford on
the residual, no potentials, no Dinic, its own arc order), and ``check_solve``, the
certificate: feasibility, conservation with supply and demand, the attained value,
the full residual closure, ``min_cut`` presence and value, dual feasibility inside
the closure with the ``NEVER_REACHED`` sentinel tested first, the virtual-terminal
convention, the curve shape and moment, and the complete curve on the independent
domain ``0..min(cap, F_max)`` at the union of both breakpoint sets.
"""

from __future__ import annotations

from collections import deque
from fractions import Fraction
from itertools import product

NEVER_REACHED = -1
VIRTUAL = (1 << 63) - 2


def exhaustive(g, view, sources, dests, edges=None, want_flows=False):
    """Minimum cost of every integer value; optionally every optimal flow."""
    edges = list(range(g.num_edges)) if edges is None else sorted(edges)
    S, D = set(sources), set(dests)
    best: dict[int, int] = {}
    argmin: dict[int, list[tuple[int, ...]]] = {}
    for flow in product(*(range(view[e] + 1) for e in edges)):
        bal = [0] * g.num_nodes
        for e, x in zip(edges, flow, strict=True):
            bal[g.src[e]] -= x
            bal[g.dst[e]] += x
        ok = True
        for v in range(g.num_nodes):
            if v in S:
                ok = bal[v] <= 0
            elif v in D:
                ok = bal[v] >= 0
            else:
                ok = bal[v] == 0
            if not ok:
                break
        if not ok:
            continue
        value = -sum(bal[s] for s in S)
        cost = sum(g.cost[e] * x for e, x in zip(edges, flow, strict=True))
        full = [0] * g.num_edges
        for e, x in zip(edges, flow, strict=True):
            full[e] = x
        if value not in best or cost < best[value]:
            best[value] = cost
            argmin[value] = [tuple(full)]
        elif cost == best[value] and want_flows:
            argmin[value].append(tuple(full))
    assert set(best) == set(range(max(best) + 1))
    return (best, argmin) if want_flows else best


def independent_curve(g, view, sources, dests, edges=None):
    """Convex breakpoints ``[(value, cost)]`` by Bellman-Ford augmentation."""
    edges = list(range(g.num_edges)) if edges is None else sorted(edges)
    n = g.num_nodes
    s_star, t_star = n, n + 1
    arcs = []  # [u, v, cap, cost]
    for e in edges:
        arcs.append([g.src[e], g.dst[e], view[e], g.cost[e]])
    for s in sources:
        arcs.append([s_star, s, VIRTUAL, 0])
    for d in dests:
        arcs.append([d, t_star, VIRTUAL, 0])
    flow = [0] * len(arcs)
    points = [(0, 0)]
    while True:
        dist: list[int | None] = [None] * (n + 2)
        dist[s_star] = 0
        parent: dict[int, tuple[int, int]] = {}
        for _ in range(n + 2):
            changed = False
            for i, (u, v, c, w) in enumerate(arcs):
                if flow[i] < c and dist[u] is not None:
                    nd = dist[u] + w
                    if dist[v] is None or nd < dist[v]:
                        dist[v] = nd
                        parent[v] = (i, 1)
                        changed = True
                if flow[i] > 0 and dist[v] is not None:
                    nd = dist[v] - w
                    if dist[u] is None or nd < dist[u]:
                        dist[u] = nd
                        parent[u] = (i, -1)
                        changed = True
            if not changed:
                break
        if dist[t_star] is None:
            return points
        path = []
        v = t_star
        while v != s_star:
            i, sign = parent[v]
            path.append((i, sign))
            v = arcs[i][0] if sign == 1 else arcs[i][1]
        amount = min(arcs[i][2] - flow[i] if sign == 1 else flow[i] for i, sign in path)
        for i, sign in path:
            flow[i] += sign * amount
        x, c = points[-1]
        points.append((x + amount, c + amount * dist[t_star]))


def value_at(points, x):
    for (a, c), (b, d) in zip(points, points[1:], strict=False):
        if a <= x <= b:
            return Fraction(c) + Fraction(d - c, b - a) * (x - a)
    if points and x == points[-1][0]:
        return Fraction(points[-1][1])
    raise ValueError(f"value {x} outside the curve domain")


def curve_points(curve):
    pts = [(0, 0)]
    for c, u in curve:
        pts.append((pts[-1][0] + u, pts[-1][1] + c * u))
    return pts


def check_solve(g, view, sources, dests, cap, sol, *, edges=None, tier="b"):
    """Return the list of failed certificate clauses (empty when certified)."""
    fails: list[str] = []
    active = set(range(g.num_edges)) if edges is None else set(edges)
    S, D = sorted(set(sources)), sorted(set(dests))
    f = sol["edge_flows"]
    placed = sol["placed"]
    pi = sol["potentials"]
    curve = [tuple(t) for t in sol["cost_distribution"]]
    for e in range(g.num_edges):
        if e not in active and f[e]:
            fails.append(f"flow on edge {e} outside the view")
        if not 0 <= f[e] <= view[e]:
            fails.append(f"edge {e}: flow {f[e]} outside [0, {view[e]}]")
    bal = [0] * g.num_nodes
    for e in range(g.num_edges):
        bal[g.src[e]] -= f[e]
        bal[g.dst[e]] += f[e]
    supply = {s: -bal[s] for s in S}
    inflow = {d: bal[d] for d in D}
    for v in range(g.num_nodes):
        if v in supply:
            if supply[v] < 0:
                fails.append(f"source {v} absorbs flow")
        elif v in inflow:
            if inflow[v] < 0:
                fails.append(f"destination {v} emits flow")
        elif bal[v] != 0:
            fails.append(f"conservation at node {v}")
    if sum(supply.values()) != placed or sum(inflow.values()) != placed:
        fails.append("value differs from supply or inflow")
    if cap is not None and not 0 <= placed <= cap:
        fails.append("placed outside [0, cap]")
    if cap == 0:
        if pi != [NEVER_REACHED] * g.num_nodes:
            fails.append("zero cap: potentials must all be NEVER_REACHED")
        if curve or sol["min_cut"] is not None or any(sol["counters"].values()):
            fails.append("zero cap: empty curve, min_cut None, counters 0")
        return fails
    # full residual closure of the super-source
    n = g.num_nodes
    s_star, t_star = n, n + 1
    out: list[list[tuple[int, int, int]]] = [[] for _ in range(n + 2)]
    for e in active:
        u, v = g.src[e], g.dst[e]
        if f[e] < view[e]:
            out[u].append((v, g.cost[e], e))
        if f[e] > 0:
            out[v].append((u, -g.cost[e], e))
    for s in S:
        out[s_star].append((s, 0, -1))
        if supply[s] > 0:
            out[s].append((s_star, 0, -1))
    for d in D:
        out[d].append((t_star, 0, -1))
        if inflow[d] > 0:
            out[t_star].append((d, 0, -1))
    seen = set(S)  # seeded at the real sources (revision 20)
    q = deque(S)
    while q:
        u = q.popleft()
        for v, _w, _e in out[u]:
            if v not in seen:
                seen.add(v)
                q.append(v)
    connected = t_star in seen
    if cap is None and connected:
        fails.append("uncapped solve left the destinations connected")
    if cap is not None and placed != cap and connected:
        fails.append("attained value: placed < cap with a residual path left")
    mc = sol["min_cut"]
    if connected:
        if mc is not None:
            fails.append("min_cut present on a connected (non-maximal) solve")
    else:
        expect = sorted(
            e
            for e in active
            if g.src[e] in seen and g.dst[e] not in seen and view[e] > 0
        )
        if mc is None:
            fails.append("min_cut missing on a maximal solve")
        elif sorted(mc) != expect or list(mc) != sorted(mc):
            fails.append(f"min_cut {mc} differs from the closure cut {expect}")
        if sum(view[e] for e in expect) != placed:
            fails.append("cut value differs from placed")
    # dual feasibility inside the closure (sentinel tested before arithmetic)
    for v in range(n):
        if v in seen and pi[v] == NEVER_REACHED:
            fails.append(f"closure vertex {v} holds NEVER_REACHED")
    for s in S:
        if s in seen and pi[s] != NEVER_REACHED and pi[s] != 0:
            fails.append(f"reached root {s} has potential {pi[s]}")
    for u in range(n):
        if u not in seen or pi[u] == NEVER_REACHED:
            continue
        for v, w, _e in out[u]:
            if v >= n or v not in seen or pi[v] == NEVER_REACHED:
                continue
            if w + pi[u] - pi[v] < 0:
                fails.append(f"negative reduced cost on arc {u}->{v}")
    reached = [d for d in D if pi[d] != NEVER_REACHED]
    if reached:
        pt = min(pi[d] for d in reached)
        for d in reached:
            if inflow[d] > 0 and pi[d] != pt:
                fails.append(f"destination {d} with inflow: potential {pi[d]} != {pt}")
            if pi[d] < pt:
                fails.append("destination potential below pi(t*)")
    elif placed:
        fails.append("flow delivered but no destination reached")
    # curve shape and moment
    if any(u <= 0 for _c, u in curve):
        fails.append("curve volume not positive")
    if any(b[0] <= a[0] for a, b in zip(curve, curve[1:], strict=False)):
        fails.append("curve costs not strictly increasing")
    if sum(u for _c, u in curve) != placed:
        fails.append("curve total differs from placed")
    if sum(c * u for c, u in curve) != sum(
        g.cost[e] * f[e] for e in range(g.num_edges)
    ):
        fails.append("curve moment differs from the flow cost")
    # complete curve on the independent domain
    pts = curve_points(curve)
    if tier == "a":
        best = exhaustive(g, view, S, D, edges=active)
        fmax = max(best)
        top = fmax if cap is None else min(cap, fmax)
        if placed != top:
            fails.append(f"placed {placed} != min(cap, F_max) = {top}")
        for x in range(min(top, placed) + 1):
            if value_at(pts, x) != best[x]:
                fails.append(f"curve differs from the exhaustive optimum at {x}")
                break
    else:
        ref = independent_curve(g, view, S, D, edges=active)
        fmax = ref[-1][0]
        top = fmax if cap is None else min(cap, fmax)
        if placed != top:
            fails.append(f"placed {placed} != min(cap, F_max) = {top}")
        xs = {x for x, _ in pts} | {x for x, _ in ref}
        for x in sorted(xs):
            if x > min(top, placed):
                continue
            if value_at(pts, x) != value_at(ref, x):
                fails.append(f"curve differs from the independent solver at {x}")
                break
    return fails
