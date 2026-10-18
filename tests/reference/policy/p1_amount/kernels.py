"""Placement kernels in units (APPROACH revision 24, sections 2.1, 7 and 16).

A plan is frozen at creation: the flow's route (a DAG of edge ids), its split graph
(members pruned to the walks from the source to the destination set) and, for
Weighted, the weights. Membership: MaxFlow and EqualBalanced take the route (an edge of
capacity 0 is a member that carries nothing; EqualBalanced still gives it its share);
Weighted the edges of positive weight by its ``priority`` axis (``Residual``: residual
>= 1 unit on the creation view; ``Capacity``: static capacity >= 1 unit); Lossy keeps
its kernel's membership.

* ``alloc(plan, T)``: node by node in topological order, the source's inflow being
  ``T``; round robin (EqualBalanced) or highest averages over the frozen weights
  (Weighted) at every node; member-monotone in ``T`` (house monotonicity).
* ``t_max(plan, view)``: the largest ``T`` whose allocation fits every member,
  by binary search over units (exact because the allocation is monotone).
* ``maxflow(...)``: the flow core restricted to the route on the view (one
  objective: minimum cost for the maximum value under the cap).
* ``lossy(plan, offer, W)``: the offer split at every node by the kernel's rule,
  each edge carrying ``min(share, residual)``, the rest a transit drop on that edge;
  ``+inf`` injects every source member to its residual.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import flowcore
import search
from apportion import split
from graph import Graph
from units import ContractError

MAX_FLOW = "MaxFlow"
WEIGHTED = "Weighted"
EQUAL_BALANCED = "EqualBalanced"
KERNELS = (MAX_FLOW, WEIGHTED, EQUAL_BALANCED)
LOSSLESS = "Lossless"
LOSSY = "Lossy"


def check_kernel(kernel: str, admission: str, priority: str) -> None:
    if kernel not in KERNELS:
        raise ContractError(f"unknown placement {kernel!r}")
    if admission not in (LOSSLESS, LOSSY):
        raise ContractError(f"unknown admission {admission!r}")
    if admission == LOSSY and kernel == MAX_FLOW:
        raise ContractError("Lossy is rejected with MaxFlow (rule 13)")
    if kernel == WEIGHTED and priority == "Headroom":
        raise ContractError("Weighted with priority = Headroom is rejected")


@dataclass
class Plan:
    kernel: str
    admission: str
    src: int
    dsts: tuple[int, ...]
    route: tuple[int, ...]
    members: tuple[int, ...]
    weights: dict[int, int] | None
    out: dict[int, list[int]] = field(default_factory=dict)
    topo: list[int] = field(default_factory=list)

    @property
    def rule(self) -> str:
        return "ha" if self.kernel == WEIGHTED else "rr"

    def to_json(self) -> dict:
        return {
            "kernel": self.kernel,
            "admission": self.admission,
            "route": list(self.route),
            "members": list(self.members),
            "weights": None
            if self.weights is None
            else [[e, w] for e, w in sorted(self.weights.items())],
        }


def prune(g: Graph, edges, src: int, dsts) -> tuple[int, ...]:
    """Edges on some walk from ``src`` to the destination set (destinations terminal)."""
    dset = set(dsts)
    edges = [e for e in edges if g.src[e] not in dset]
    fwd = {src}
    changed = True
    while changed:
        changed = False
        for e in edges:
            if g.src[e] in fwd and g.dst[e] not in fwd:
                fwd.add(g.dst[e])
                changed = True
    back = set(dset)
    changed = True
    while changed:
        changed = False
        for e in edges:
            if g.dst[e] in back and g.src[e] not in back:
                back.add(g.src[e])
                changed = True
    return tuple(sorted(e for e in edges if g.src[e] in fwd and g.dst[e] in back))


def make_plan(
    g: Graph,
    kernel: str,
    admission: str,
    src: int,
    dsts,
    route,
    priority: str,
    view,
) -> Plan:
    """Capture membership and weights on the creation view (frozen until reroute)."""
    check_kernel(kernel, admission, priority)
    route = tuple(sorted(route))
    dsts = tuple(sorted(dsts))
    weights = None
    if kernel != WEIGHTED:
        cand = route
    elif priority == "Residual":
        cand = tuple(e for e in route if view[e] >= 1)
    else:
        cand = tuple(e for e in route if g.capacity[e] >= 1)
    members = prune(g, cand, src, dsts)
    if kernel == WEIGHTED:
        weights = {
            e: (view[e] if priority == "Residual" else g.capacity[e]) for e in members
        }
    plan = Plan(kernel, admission, src, dsts, route, members, weights)
    out: dict[int, list[int]] = {}
    indeg: dict[int, int] = {src: 0}
    for e in members:
        out.setdefault(g.src[e], []).append(e)
        indeg[g.dst[e]] = indeg.get(g.dst[e], 0) + 1
        indeg.setdefault(g.src[e], 0)
    for lst in out.values():
        lst.sort()
    topo = []
    ready = sorted(v for v, d in indeg.items() if d == 0)
    while ready:
        u = ready.pop(0)
        topo.append(u)
        for e in out.get(u, ()):
            v = g.dst[e]
            indeg[v] -= 1
            if indeg[v] == 0:
                ready.append(v)
                ready.sort()
    assert len(topo) == len(indeg), "split graph is cyclic"
    plan.out, plan.topo = out, topo
    return plan


def alloc(g: Graph, plan: Plan, total: int) -> dict[int, int]:
    """Cumulative allocation ``alloc_i(T)`` of a Weighted or EqualBalanced plan."""
    if plan.kernel == MAX_FLOW:
        raise ValueError("MaxFlow allocations come from the flow core")
    if total == 0:
        return {}
    if not plan.members:
        raise ValueError("an empty split graph carries nothing")
    inflow = {plan.src: total}
    out: dict[int, int] = {}
    dset = set(plan.dsts)
    for u in plan.topo:
        amount = inflow.get(u, 0)
        if u in dset or amount == 0:
            continue
        outs = plan.out.get(u, [])
        assert outs, "a non-destination node of the split graph has members"
        weights = None if plan.weights is None else [plan.weights[e] for e in outs]
        for e, x in zip(
            outs, split(plan.rule, amount, weights, len(outs)), strict=True
        ):
            if x:
                out[e] = x
                inflow[g.dst[e]] = inflow.get(g.dst[e], 0) + x
    assert sum(inflow.get(d, 0) for d in dset) == total
    return out


def fits(allocation: dict[int, int], view) -> bool:
    return all(x <= view[e] for e, x in allocation.items())


@dataclass
class TMax:
    value: int
    evaluations: int
    upper: int


def t_max(g: Graph, plan: Plan, view) -> TMax:
    """Largest ``T`` with ``alloc(T)`` within ``view`` on every member (exact)."""
    if not plan.members:
        return TMax(0, 0, 0)
    hi = sum(view[e] for e in plan.out.get(plan.src, []))
    upper = hi
    lo = 0
    evals = 0
    while lo < hi:
        mid = (lo + hi + 1) // 2
        evals += 1
        if fits(alloc(g, plan, mid), view):
            lo = mid
        else:
            hi = mid - 1
    return TMax(lo, evals, upper)


def maxflow(
    g: Graph,
    plan: Plan,
    view,
    cap: int | None,
    *,
    arc_order="interleaved",
    arc_rank=None,
    trace=False,
):
    """The MaxFlow kernel: a fresh flow-core solve on the view over the route."""
    return flowcore.solve(
        g,
        view,
        [plan.src],
        plan.dsts,
        cap,
        edges=plan.members,
        arc_order=arc_order,
        arc_rank=arc_rank,
        trace=trace,
    )


@dataclass
class LossyCall:
    offered: int
    injected: int
    delivered: int
    carried: dict[int, int]
    drops: dict[int, int]
    shares: dict[int, int]

    @property
    def transit_drops(self) -> int:
        return sum(self.drops.values())

    @property
    def source_drops(self) -> int:
        return self.offered - self.delivered - self.transit_drops


def lossy(g: Graph, plan: Plan, offer: int | None, W) -> LossyCall:
    """One Lossy offer on the working residual ``W`` (``None`` = ``+inf``)."""
    carried: dict[int, int] = {}
    drops: dict[int, int] = {}
    shares: dict[int, int] = {}
    src_outs = plan.out.get(plan.src, [])
    if not src_outs:
        offered = 0 if offer is None else offer
        return LossyCall(offered, 0, 0, carried, drops, shares)
    dset = set(plan.dsts)
    inflow: dict[int, int] = {}
    injected = 0
    for u in plan.topo:
        if u in dset:
            continue
        outs = plan.out.get(u, [])
        if u == plan.src:
            if offer is None:
                part = [W[e] for e in outs]
            else:
                weights = (
                    None if plan.weights is None else [plan.weights[e] for e in outs]
                )
                part = split(plan.rule, offer, weights, len(outs))
            injected = sum(part)
        else:
            amount = inflow.get(u, 0)
            if amount == 0:
                continue
            weights = None if plan.weights is None else [plan.weights[e] for e in outs]
            part = split(plan.rule, amount, weights, len(outs))
        for e, x in zip(outs, part, strict=True):
            if not x:
                continue
            c = min(x, W[e])
            shares[e] = x
            if c:
                carried[e] = c
                inflow[g.dst[e]] = inflow.get(g.dst[e], 0) + c
            if x - c:
                drops[e] = x - c
    delivered = sum(inflow.get(d, 0) for d in dset)
    offered = injected if offer is None else offer
    return LossyCall(offered, injected, delivered, carried, drops, shares)


def route_cost(g: Graph, edges, src: int, dsts) -> int | None:
    """Minimum walk cost of a route from ``src`` to the destination set."""
    d = search.distances(g, edges, [src])
    reach = [d[t] for t in dsts if d[t] is not None]
    return min(reach) if reach else None
