"""The one integer minimum-cost flow core (APPROACH revision 24, section 2.2).

One objective on every view: minimum cost for the maximum value attainable under the
value cap. Successive shortest paths with potentials and blocking flows: each tier
runs Dijkstra on reduced costs to completion over the reachable residual (from the
super-source, through the super-sink and onward), adds the distances to the
potentials of the vertices it settles, then drains the zero-reduced-cost admissible
subgraph by Dinic phases (BFS levels, DFS with current-arc pointers) until no
admissible augmenting path remains or the cap is reached; then the next tier.

Arc identity: every edge is its own arc. At a node the forward arcs of its
out-edges and the reverse arcs of its in-edges are visited in one sequence ordered
by ascending EdgeId (a self-loop's forward arc before its reverse arc); at a real
endpoint node its virtual arc comes first; the super-source visits its arcs by
ascending source NodeId and the super-sink by ascending destination NodeId.
Virtual arcs cost 0 and carry ``2^63 - 2`` units (never saturated). Every solve starts from zero flow
on the view, so potentials 0 are valid.

``arc_order`` and ``arc_rank`` are mutant switches for the adjacency-order controls
only: ``"interleaved"`` is the contract and ``"forward_first"`` visits every forward
arc before any reverse arc; ``arc_rank[e]`` replaces the EdgeId as the sort key
(for example the old global ``(cost, src, dst)`` rank of layout B).
"""

from __future__ import annotations

import heapq
from collections import deque
from dataclasses import dataclass, field

from graph import Eligibility, Graph, check_disjoint, check_node_set
from units import (
    INTERNAL_CAP_MAX,
    NEVER_REACHED,
    VIRTUAL_CAP,
    ContractError,
)

COUNTER_NAMES = (
    "dijkstra_runs",
    "dijkstra_pops",
    "settled",
    "arc_scans",
    "dinic_phases",
    "augmentations",
)


@dataclass
class FlowSolve:
    """The private payload of ``_place_global_units`` plus diagnostics."""

    placed: int
    cost_distribution: list[tuple[int, int]]
    edge_flows: list[int]
    min_cut: list[int] | None
    potentials: list[int]
    counters: dict[str, int]
    # Diagnostics (not payload fields): per-source supply, per-destination inflow,
    # the internal super-sink potential, the full residual closure, tier traces.
    supply: dict[int, int] = field(default_factory=dict)
    inflow: dict[int, int] = field(default_factory=dict)
    pi_sink: int | None = None
    closure: list[int] = field(default_factory=list)
    sink_in_closure: bool = False
    tiers: list[dict] = field(default_factory=list)

    def payload(self) -> dict:
        return {
            "placed": self.placed,
            "cost_distribution": [list(t) for t in self.cost_distribution],
            "edge_flows": list(self.edge_flows),
            "min_cut": None if self.min_cut is None else list(self.min_cut),
            "potentials": list(self.potentials),
            "counters": dict(self.counters),
        }


def place_global_units(
    g: Graph,
    view,
    sources,
    destinations,
    elig: Eligibility | None,
    cap: int | None,
) -> FlowSolve:
    """The recorder's private entry point: validated, whole eligible graph."""
    return solve(g, view, sources, destinations, cap, elig=elig)


def _validate(g, view, sources, destinations, cap, elig):
    if len(view) != g.num_edges:
        raise ContractError("view length differs from the edge count")
    for e, x in enumerate(view):
        if not isinstance(x, int) or isinstance(x, bool):
            raise ContractError("view entries must be integer units")
        if not 0 <= x <= g.capacity[e]:
            raise ContractError(f"view[{e}] = {x} outside [0, capacity]")
    s = check_node_set(g, sources, "sources")
    d = check_node_set(g, destinations, "destinations")
    check_disjoint(s, d)
    if cap is not None:
        if not isinstance(cap, int) or isinstance(cap, bool):
            raise ContractError("cap must be an int or None")
        if not 0 <= cap <= INTERNAL_CAP_MAX:
            raise ContractError("cap outside [0, 2^63 - 2]")
    if elig is not None:
        elig.check(g)
    return s, d


def solve(
    g: Graph,
    view,
    sources,
    destinations,
    cap: int | None,
    *,
    edges=None,
    elig: Eligibility | None = None,
    arc_order: str = "interleaved",
    arc_rank=None,
    trace: bool = False,
) -> FlowSolve:
    """Minimum-cost flow of value ``min(cap, max)`` on ``view`` restricted to ``edges``."""
    srcs, dsts = _validate(g, view, sources, destinations, cap, elig)
    n = g.num_nodes
    S, T = n, n + 1
    allowed = None if elig is None else elig.edge_allowed
    node_ok = None if elig is None else elig.node_allowed
    if edges is None:
        edges = range(g.num_edges)
    active = sorted(
        {e for e in edges if allowed is None or allowed[e]},
    )
    srcs = [s for s in srcs if node_ok is None or node_ok[s]]
    dsts = [d for d in dsts if node_ok is None or node_ok[d]]

    counters = dict.fromkeys(COUNTER_NAMES, 0)
    edge_flows = [0] * g.num_edges
    if cap == 0:
        # Zero cap: no Dijkstra, every potential NEVER_REACHED, an empty curve,
        # min_cut None (the cap was attained), every counter 0.
        return FlowSolve(
            placed=0,
            cost_distribution=[],
            edge_flows=edge_flows,
            min_cut=None,
            potentials=[NEVER_REACHED] * n,
            counters=counters,
            supply=dict.fromkeys(srcs, 0),
            inflow=dict.fromkeys(dsts, 0),
        )

    # ---- arcs: index j and j ^ 1 are a forward/reverse pair
    head: list[int] = []
    res: list[int] = []
    acost: list[int] = []
    aedge: list[int] = []  # edge id, or -1 for virtual arcs

    def add(u: int, v: int, c: int, w: int, e: int) -> int:
        j = len(head)
        head.extend((v, u))
        res.extend((c, 0))
        acost.extend((w, -w))
        aedge.extend((e, e))
        return j

    tail: list[int] = []
    real_arcs: list[list[tuple]] = [[] for _ in range(n)]
    edge_arc: dict[int, int] = {}
    if arc_order not in ("interleaved", "forward_first"):
        raise ValueError(f"unknown arc order {arc_order!r}")
    for e in active:
        u, v = g.src[e], g.dst[e]
        j = add(u, v, view[e], g.cost[e], e)
        tail.extend((u, v))
        edge_arc[e] = j
        key = e if arc_rank is None else arc_rank[e]
        if arc_order == "forward_first":
            real_arcs[u].append((0, key, 0, j))
            real_arcs[v].append((1, key, 1, j + 1))
        else:
            real_arcs[u].append((key, 0, j))
            real_arcs[v].append((key, 1, j + 1))
    src_arc: dict[int, int] = {}
    for s in srcs:
        j = add(S, s, VIRTUAL_CAP, 0, -1)
        tail.extend((S, s))
        src_arc[s] = j
    dst_arc: dict[int, int] = {}
    for d in dsts:
        j = add(d, T, VIRTUAL_CAP, 0, -1)
        tail.extend((d, T))
        dst_arc[d] = j
    adj: list[list[int]] = [[] for _ in range(n + 2)]
    adj[S] = [src_arc[s] for s in srcs]
    adj[T] = [dst_arc[d] + 1 for d in dsts]
    for u in range(n):
        first = []
        if u in src_arc:
            first.append(src_arc[u] + 1)
        if u in dst_arc:
            first.append(dst_arc[u])
        adj[u] = first + [item[-1] for item in sorted(real_arcs[u])]

    # A positive cap always runs the solver; the clamp of the internal domain at the
    # sources' view out-capacity only bounds the arithmetic (Python ints are exact),
    # so the loop runs on the input cap (revision 20).
    limit = cap

    pi = [0] * (n + 2)
    ever = [False] * (n + 2)
    total = 0
    curve: list[tuple[int, int]] = []
    tiers: list[dict] = []

    def reduced(j: int) -> int:
        return acost[j] + pi[tail[j]] - pi[head[j]]

    while limit is None or total < limit:
        # ---- Dijkstra on reduced costs, to completion over the reachable residual
        counters["dijkstra_runs"] += 1
        dist: list[int | None] = [None] * (n + 2)
        dist[S] = 0
        heap = [(0, S)]
        settled_now: list[int] = []
        while heap:
            du, u = heapq.heappop(heap)
            counters["dijkstra_pops"] += 1
            if du != dist[u]:
                continue
            settled_now.append(u)
            counters["settled"] += 1
            for j in adj[u]:
                counters["arc_scans"] += 1
                if res[j] <= 0:
                    continue
                v = head[j]
                rc = reduced(j)
                assert rc >= 0, ("negative reduced cost", j, rc)
                nd = du + rc
                dv = dist[v]
                if dv is None or nd < dv:
                    dist[v] = nd
                    heapq.heappush(heap, (nd, v))
        for u in settled_now:
            pi[u] += dist[u]
            ever[u] = True
        if dist[T] is None:
            break
        tier_cost = pi[T] - pi[S]
        before = total
        # ---- drain the admissible subgraph by Dinic phases
        while limit is None or total < limit:
            level = [-1] * (n + 2)
            level[S] = 0
            q = deque([S])
            while q:
                u = q.popleft()
                for j in adj[u]:
                    counters["arc_scans"] += 1
                    v = head[j]
                    if res[j] > 0 and level[v] < 0 and reduced(j) == 0:
                        level[v] = level[u] + 1
                        q.append(v)
            # every BFS level build counts, the last unsuccessful one of a tier included
            # (revision 23)
            counters["dinic_phases"] += 1
            if level[T] < 0:
                break
            cur = [0] * (n + 2)
            while limit is None or total < limit:
                path = _find_path(S, T, adj, cur, res, head, level, reduced, counters)
                if path is None:
                    break
                amount = min(res[j] for j in path)
                if limit is not None:
                    amount = min(amount, limit - total)
                for j in path:
                    res[j] -= amount
                    res[j ^ 1] += amount
                total += amount
                counters["augmentations"] += 1
        drained = total - before
        assert drained > 0, "a successful Dijkstra drains at least one unit"
        curve.append((tier_cost, drained))
        if trace:
            tiers.append(
                {
                    "cost": tier_cost,
                    "units": drained,
                    "edge_flows": {
                        e: res[edge_arc[e] + 1] for e in active if res[edge_arc[e] + 1]
                    },
                    "counters": dict(counters),
                }
            )

    for e in active:
        edge_flows[e] = res[edge_arc[e] + 1]
    # ---- full residual closure seeded at the real sources (reverse and virtual
    # arcs included; revision 20)
    seen = [False] * (n + 2)
    for s in srcs:
        seen[s] = True
    q = deque(srcs)
    while q:
        u = q.popleft()
        for j in adj[u]:
            v = head[j]
            if res[j] > 0 and not seen[v]:
                seen[v] = True
                q.append(v)
    if seen[T]:
        min_cut = None
    else:
        min_cut = sorted(
            e for e in active if seen[g.src[e]] and not seen[g.dst[e]] and view[e] > 0
        )
    return FlowSolve(
        placed=total,
        cost_distribution=curve,
        edge_flows=edge_flows,
        min_cut=min_cut,
        potentials=[pi[v] if ever[v] else NEVER_REACHED for v in range(n)],
        counters=counters,
        supply={s: res[src_arc[s] + 1] for s in srcs},
        inflow={d: res[dst_arc[d] + 1] for d in dsts},
        pi_sink=pi[T] if ever[T] else None,
        closure=[v for v in range(n) if seen[v]],
        sink_in_closure=seen[T],
        tiers=tiers,
    )


def _find_path(S, T, adj, cur, res, head, level, reduced, counters):
    """One augmenting path in the level graph (current-arc pointers kept)."""
    stack = [S]
    path: list[int] = []
    while stack:
        u = stack[-1]
        if u == T:
            return path
        advanced = False
        arcs = adj[u]
        while cur[u] < len(arcs):
            j = arcs[cur[u]]
            counters["arc_scans"] += 1
            v = head[j]
            if res[j] > 0 and level[v] == level[u] + 1 and reduced(j) == 0:
                stack.append(v)
                path.append(j)
                advanced = True
                break
            cur[u] += 1
        if not advanced:
            stack.pop()
            if path:
                j = path.pop()
                # The arc into the dead end fails: advance the parent's pointer.
                cur[stack[-1]] += 1
    return None


def flow_cost(g: Graph, flows) -> int:
    return sum(g.cost[e] * x for e, x in enumerate(flows))
