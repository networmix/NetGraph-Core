"""Placement kernels in units (APPROACH revision 24, sections 2.1, 7 and 16).

The kernel names and the admission rule, a flow's plan as a record (its route, its members and,
for Weighted, its weights), the split graph of a route (``prune``: the edges on the walks from a
source to the destination set) and a route's cost. The per-flow kernels of the per-source
transaction (allocation, the T_max search, the per-flow MaxFlow solve, the Lossy split) left with
``CapacityAware`` when T0 was deleted (README "T1 replaces CapacityAware"): the pooled
CostOnly program and the controllers of ``policy.py`` forward and solve on their own.
"""

from __future__ import annotations

from dataclasses import dataclass

import search
from graph import Graph
from units import ContractError

MAX_FLOW = "MaxFlow"
WEIGHTED = "Weighted"
EQUAL_BALANCED = "EqualBalanced"
KERNELS = (MAX_FLOW, WEIGHTED, EQUAL_BALANCED)
LOSSLESS = "Lossless"
LOSSY = "Lossy"


def check_kernel(kernel: str, admission: str) -> None:
    if kernel not in KERNELS:
        raise ContractError(f"unknown placement {kernel!r}")
    if admission not in (LOSSLESS, LOSSY):
        raise ContractError(f"unknown admission {admission!r}")
    if admission == LOSSY and kernel == MAX_FLOW:
        raise ContractError("Lossy is rejected with MaxFlow (rule 13)")


@dataclass
class Plan:
    kernel: str
    admission: str
    src: int
    dsts: tuple[int, ...]
    route: tuple[int, ...]
    members: tuple[int, ...]
    weights: dict[int, int] | None

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


def route_cost(g: Graph, edges, src: int, dsts) -> int | None:
    """Minimum walk cost of a route from ``src`` to the destination set."""
    d = search.distances(g, edges, [src])
    reach = [d[t] for t in dsts if d[t] is not None]
    return min(reach) if reach else None
