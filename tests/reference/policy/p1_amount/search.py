"""Search semantics of APPROACH revision 24, section 6, as internal routines.

* member graph: static eligibility (capacity is never read: an edge of capacity 0 is a
  member like any other), optional residual gating at ``min_residual``; ``multi_edge = False`` keeps one member per ``(u, v)`` by
  ``select_adjacency``: minimum cost, then (``PreferHigherPriority``) the highest
  priority value, then the lowest EdgeId;
* the unique shortest-path DAG (positive costs, 2.4): tight edges; with targets the
  walks from the roots to the nearest targets (anycast), one edge set for both
  directions;
* the staged canonical single path: (1) the tight DAG of the pair, (2) under
  ``PreferHigherPriority`` the best attainable bottleneck ``B*`` (capped by a Yen
  prefix's bottleneck) and the DAG restricted to priority ``>= B*``, (3) one reverse
  hop count ``h`` and a greedy walk taking the smallest tight edge whose head has
  ``h(v) = h(u) - 1``;

Priority values: ``Residual`` = residual units, ``Capacity`` = static capacity
units, ``Headroom`` = the exact rational ``residual / capacity`` (0 on an edge of
capacity 0). A search without a residual reads the static capacity as the residual:
``Residual`` and ``Capacity`` the capacity, ``Headroom`` 1 (0 on an edge of capacity 0)
(contract r8 3.5, step 1 of the next Core change).
"""

from __future__ import annotations

import heapq
import math
from collections import deque
from dataclasses import dataclass
from fractions import Fraction

from graph import Eligibility, Graph, check_disjoint, check_node_set
from units import ContractError

DETERMINISTIC = "Deterministic"
PREFER = "PreferHigherPriority"
PRIORITIES = ("Residual", "Capacity", "Headroom")


@dataclass(frozen=True)
class Selection:
    """PathSelection {multipath, multi_edge, tie_break, priority}."""

    multipath: bool = True
    multi_edge: bool = True
    tie_break: str = DETERMINISTIC
    priority: str = "Residual"

    def __post_init__(self) -> None:
        if self.tie_break not in (DETERMINISTIC, PREFER):
            raise ContractError(f"unknown tie_break {self.tie_break!r}")
        if self.priority not in PRIORITIES:
            raise ContractError(f"unknown priority {self.priority!r}")


def priority_fn(g: Graph, priority: str, residual):
    """Per-edge priority value (exact: ints or Fractions)."""
    free = (lambda e: g.capacity[e]) if residual is None else (lambda e: residual[e])
    if priority == "Capacity":
        return lambda e: g.capacity[e]
    if priority == "Residual":
        return free
    if priority == "Headroom":
        return lambda e: Fraction(free(e), g.capacity[e]) if g.capacity[e] else 0
    raise ContractError(f"unknown priority {priority!r}")


def check_residual(g: Graph, residual) -> None:
    if residual is None:
        return
    if len(residual) != g.num_edges:
        raise ContractError("residual length differs from the edge count")
    for e, r in enumerate(residual):
        if not isinstance(r, int) or isinstance(r, bool):
            raise ContractError("residual entries must be integer units")
        if not 0 <= r <= g.capacity[e]:
            raise ContractError(f"residual[{e}] = {r} outside [0, capacity]")


def member_graph(
    g: Graph,
    elig: Eligibility | None,
    selection: Selection,
    residual=None,
    min_residual: int = 1,
) -> list[int]:
    """Eligible routable edges, collapsed per ``(u, v)`` when ``multi_edge`` is off."""
    if elig is not None:
        elig.check(g)
    if min_residual < 1:
        raise ContractError("min_residual is at least 1 unit")
    ok = []
    for e in range(g.num_edges):
        if elig is not None and not elig.edge_allowed[e]:
            continue
        if residual is not None and residual[e] < min_residual:
            continue
        ok.append(e)
    if selection.multi_edge:
        return ok
    return collapse(g, ok, selection, residual)


def collapse(g: Graph, edges, selection: Selection, residual) -> list[int]:
    prio = priority_fn(g, selection.priority, residual)
    best: dict[tuple[int, int], int] = {}
    for e in edges:
        key = (g.src[e], g.dst[e])
        cur = best.get(key)
        if cur is None or _select_key(g, e, selection, prio) < _select_key(
            g, cur, selection, prio
        ):
            best[key] = e
    return sorted(best.values())


def _select_key(g, e, selection, prio):
    if selection.tie_break == PREFER:
        return (g.cost[e], -prio(e), e)
    return (g.cost[e], 0, e)


def distances(
    g: Graph, edges, roots, reverse: bool = False, excluded_nodes=frozenset()
) -> list[int | None]:
    """Dijkstra over ``edges`` from (or, reverse, toward) the root set."""
    adj: list[list[tuple[int, int]]] = [[] for _ in range(g.num_nodes)]
    for e in edges:
        u, v = g.src[e], g.dst[e]
        if u in excluded_nodes or v in excluded_nodes:
            continue
        if reverse:
            adj[v].append((u, g.cost[e]))
        else:
            adj[u].append((v, g.cost[e]))
    dist: list[int | None] = [None] * g.num_nodes
    heap = []
    for r in roots:
        if r in excluded_nodes:
            continue
        dist[r] = 0
        heap.append((0, r))
    heapq.heapify(heap)
    while heap:
        d, u = heapq.heappop(heap)
        if d != dist[u]:
            continue
        for v, c in adj[u]:
            nd = d + c
            dv = dist[v]
            if dv is None or nd < dv:
                dist[v] = nd
                heapq.heappush(heap, (nd, v))
    return dist


def full_dag(g: Graph, edges, roots, reverse: bool = False) -> list[int]:
    """All tight edges over the nodes reached from (or toward) the roots."""
    d = distances(g, edges, roots, reverse)
    out = []
    for e in edges:
        u, v = g.src[e], g.dst[e]
        if reverse:
            if d[v] is not None and d[u] is not None and d[v] + g.cost[e] == d[u]:
                out.append(e)
        elif d[u] is not None and d[v] is not None and d[u] + g.cost[e] == d[v]:
            out.append(e)
    return sorted(out)


def tight_dag(
    g: Graph, edges, sources, destinations, excluded_nodes=frozenset()
) -> tuple[int | None, list[int]]:
    """Edges on the minimum-cost walks from the source set to the destination set.

    The walks of cost ``d* = min dist(S, D)``; destinations and sources are terminal
    (positive costs). Identical for a forward search with targets and a reverse
    search with the roles swapped (2.4).
    """
    ds = distances(g, edges, sources, False, excluded_nodes)
    dd = distances(g, edges, destinations, True, excluded_nodes)
    reach = [ds[t] for t in destinations if ds[t] is not None]
    if not reach:
        return None, []
    dstar = min(reach)
    out = []
    for e in edges:
        u, v = g.src[e], g.dst[e]
        if u in excluded_nodes or v in excluded_nodes:
            continue
        du, dv = ds[u], dd[v]
        if du is not None and dv is not None and du + g.cost[e] + dv == dstar:
            out.append(e)
    return dstar, sorted(out)


@dataclass
class SinglePath:
    edges: list[int]
    cost: int
    bottleneck: object  # min priority value along the walk (None: Deterministic)
    b_star: object  # the stage-2 threshold (None: Deterministic)


def canonical_path(
    g: Graph,
    edges,
    sources,
    destinations,
    prio=None,
    prefix_cap=None,
    excluded_nodes=frozenset(),
) -> SinglePath | None:
    """The staged canonical single path (section 6) on the member edges ``edges``.

    ``prio`` is the priority function under ``PreferHigherPriority`` (``None`` for
    ``Deterministic``); ``prefix_cap`` caps ``B*`` by a Yen root prefix's bottleneck.
    """
    dstar, dag = tight_dag(g, edges, sources, destinations, excluded_nodes)
    if dstar is None:
        return None
    dests = set(destinations)
    b_star = None
    if prio is not None:
        ds = distances(g, dag, sources, False, excluded_nodes)
        best: dict[int, object] = {}
        for s in sources:
            if s not in excluded_nodes and ds[s] == 0:
                best[s] = math.inf
        order = sorted(
            (v for v in range(g.num_nodes) if ds[v] is not None), key=ds.__getitem__
        )
        out_adj: dict[int, list[int]] = {}
        for e in dag:
            out_adj.setdefault(g.src[e], []).append(e)
        for u in order:
            if u not in best:
                continue
            for e in out_adj.get(u, ()):
                cand = min(best[u], prio(e))
                v = g.dst[e]
                if v not in best or cand > best[v]:
                    best[v] = cand
        b_star = max(best[t] for t in destinations if t in best and ds[t] == dstar)
        if prefix_cap is not None and prefix_cap < b_star:
            b_star = prefix_cap
        dag = [e for e in dag if prio(e) >= b_star]
    # hop counts toward the destinations inside the (restricted) DAG
    radj: dict[int, list[int]] = {}
    for e in dag:
        radj.setdefault(g.dst[e], []).append(e)
    h: dict[int, int] = {t: 0 for t in destinations if t not in excluded_nodes}
    q = deque(h)
    while q:
        v = q.popleft()
        for e in radj.get(v, ()):
            u = g.src[e]
            if u not in h:
                h[u] = h[v] + 1
                q.append(u)
    starts = [s for s in sources if s in h]
    if not starts:
        return None
    hmin = min(h[s] for s in starts)
    out_adj2: dict[int, list[int]] = {}
    for e in dag:
        out_adj2.setdefault(g.src[e], []).append(e)
    frontier = [s for s in starts if h[s] == hmin]
    walk: list[int] = []
    u = None
    while True:
        cands = []
        for x in [u] if u is not None else frontier:
            for e in out_adj2.get(x, ()):
                v = g.dst[e]
                if v in h and h[v] == h[x] - 1:
                    cands.append(e)
        if not cands:
            break
        e = min(cands)
        walk.append(e)
        u = g.dst[e]
        if u in dests:
            break
    if not walk:
        return None
    cost = sum(g.cost[e] for e in walk)
    assert cost == dstar
    bottleneck = None if prio is None else min(prio(e) for e in walk)
    return SinglePath(walk, cost, bottleneck, b_star)


# ------------------------------------------------------------------ public-shaped spf


@dataclass
class SpfResult:
    costs: list[int | None] | None
    dag: list[int]
    target_cost: int | None
    path: SinglePath | None = None


def spf(
    g: Graph,
    roots,
    *,
    direction: str = "Forward",
    targets=None,
    selection: Selection | None = None,
    residual=None,
    min_residual: int = 1,
    elig: Eligibility | None = None,
) -> SpfResult:
    selection = Selection() if selection is None else selection
    roots = check_node_set(g, roots, "roots")
    if targets is not None:
        targets = check_node_set(g, targets, "targets")
        check_disjoint(roots, targets)
    elif not selection.multipath:
        raise ContractError("multipath = false requires targets (rule 13)")
    if direction not in ("Forward", "Reverse"):
        raise ContractError(f"unknown direction {direction!r}")
    check_residual(g, residual)
    members = member_graph(g, elig, selection, residual, min_residual)
    node_ok = (lambda v: True) if elig is None else (lambda v: elig.node_allowed[v])
    live_roots = [r for r in roots if node_ok(r)]
    if targets is None:
        d = distances(g, members, live_roots, direction == "Reverse")
        return SpfResult(
            d, full_dag(g, members, live_roots, direction == "Reverse"), None
        )
    live_targets = [t for t in targets if node_ok(t)]
    src, dst = (
        (live_roots, live_targets)
        if direction == "Forward"
        else (live_targets, live_roots)
    )
    if not src or not dst:
        return SpfResult(None, [], None)
    if selection.multipath:
        dstar, dag = tight_dag(g, members, src, dst)
        return SpfResult(None, dag, dstar)
    prio = (
        priority_fn(g, selection.priority, residual)
        if selection.tie_break == PREFER
        else None
    )
    path = canonical_path(g, members, src, dst, prio)
    if path is None:
        return SpfResult(None, [], None)
    return SpfResult(None, sorted(path.edges), path.cost, path)
