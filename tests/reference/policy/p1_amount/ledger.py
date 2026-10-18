"""The integer ledger ``FlowGraph`` (APPROACH revision 24, section 8).

The ledger is mutated only through policies: ``_bind`` registers a policy,
``_commit`` exchanges a policy's records atomically (feasibility checked before any
mutation, new FlowIds assigned in creation order), ``_detach`` reverts a policy's
allocations and retires its ids (``remove()``), ``_orphan`` marks a destroyed
policy's records as orphans (readable, never changed again). ``reset`` and
``copy`` are rejected while a policy is bound. Ids are never reused.
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass, field, replace

from graph import Eligibility, Graph, check_node_set
from units import COUNTER_LIMIT, UNREACHABLE, ContractError

_UID = itertools.count(1)


@dataclass
class FlowRecord:
    flow_id: tuple[int, int]
    owner: int
    sources: tuple[int, ...]
    destinations: tuple[int, ...]
    admission: str
    bundle: tuple[int, ...]
    route_cost: int
    offered: int | None
    delivered: int
    edges: dict[int, int] = field(default_factory=dict)
    drops: dict[int, int] = field(default_factory=dict)

    def snapshot(self) -> FlowRecord:
        return replace(self, edges=dict(self.edges), drops=dict(self.drops))

    def to_json(self, g: Graph | None = None) -> dict:
        pos = None if g is None else g.input_index()

        def pairs(d):
            return [
                [e, x] if pos is None else [e, pos[e], x] for e, x in sorted(d.items())
            ]

        return {
            "flow_id": list(self.flow_id),
            "owner": self.owner,
            "sources": list(self.sources),
            "destinations": list(self.destinations),
            "admission": self.admission,
            "bundle": list(self.bundle),
            "bundle_input": None
            if pos is None
            else sorted(pos[e] for e in self.bundle),
            "route_cost": None if self.route_cost == UNREACHABLE else self.route_cost,
            "offered": self.offered,
            "delivered": self.delivered,
            "edges": pairs(self.edges),
            "drops": pairs(self.drops),
        }


@dataclass
class Draft:
    """A record proposed by a policy for the commit (``flow_id`` None: new)."""

    flow_id: tuple[int, int] | None
    sources: tuple[int, ...]
    destinations: tuple[int, ...]
    admission: str
    bundle: tuple[int, ...]
    route_cost: int
    offered: int | None
    delivered: int
    edges: dict[int, int]
    drops: dict[int, int]


class FlowGraph:
    def __init__(self, g: Graph, residual=None) -> None:
        self.graph = g
        self.uid = next(_UID)
        self._baseline = self._check_baseline(residual)
        self._records: dict[int, FlowRecord] = {}
        self._next_id = 0
        self._next_owner = 0
        self._bound: set[int] = set()
        self._orphans: set[int] = set()
        self._retired: set[int] = set()

    def _check_baseline(self, residual) -> list[int]:
        g = self.graph
        if residual is None:
            return list(g.capacity)
        if len(residual) != g.num_edges:
            raise ContractError("baseline length differs from the edge count")
        out = []
        for e, r in enumerate(residual):
            if not isinstance(r, int) or isinstance(r, bool):
                raise ContractError("baselines are integer units (to_units first)")
            if not 0 <= r <= g.capacity[e]:
                raise ContractError(f"baseline[{e}] = {r} outside [0, capacity]")
            out.append(r)
        return out

    # ------------------------------------------------------------ observers

    def capacity_units(self) -> list[int]:
        return list(self.graph.capacity)

    def baseline(self) -> list[int]:
        return list(self._baseline)

    def allocated(self) -> list[int]:
        out = [0] * self.graph.num_edges
        for r in self._records.values():
            for e, x in r.edges.items():
                out[e] += x
        return out

    def residual(self) -> list[int]:
        return [b - a for b, a in zip(self._baseline, self.allocated(), strict=True)]

    def records(self) -> list[FlowRecord]:
        return [self._records[i].snapshot() for i in sorted(self._records)]

    def record(self, flow_id) -> FlowRecord:
        uid, i = flow_id
        if uid != self.uid:
            raise ContractError("FlowId of another ledger")
        if i not in self._records:
            raise ContractError(f"retired or unknown FlowId {flow_id}")
        return self._records[i].snapshot()

    def next_flow_id(self) -> int:
        return self._next_id

    def owner_records(self, owner: int) -> list[FlowRecord]:
        return [r for i, r in sorted(self._records.items()) if r.owner == owner]

    def bottleneck(self, sources, elig: Eligibility | None = None) -> list[int]:
        """Saturation diagnostic: allowed edges of capacity >= 1 leaving the set
        reachable from the sources over forward arcs with residual >= 1."""
        g = self.graph
        srcs = check_node_set(g, sources, "sources")
        if elig is not None:
            elig.check(g)
        ok = (lambda e: True) if elig is None else (lambda e: elig.edge_allowed[e])
        node_ok = (lambda v: True) if elig is None else (lambda v: elig.node_allowed[v])
        res = self.residual()
        seen = {s for s in srcs if node_ok(s)}
        stack = list(seen)
        while stack:
            u = stack.pop()
            for e in g.out_ids(u):
                v = g.dst[e]
                if ok(e) and res[e] >= 1 and v not in seen:
                    seen.add(v)
                    stack.append(v)
        return [
            e
            for e in range(g.num_edges)
            if g.src[e] in seen
            and g.dst[e] not in seen
            and ok(e)
            and g.capacity[e] >= 1
        ]

    # ------------------------------------------------------------ lifetime

    def bound_policies(self) -> int:
        return len(self._bound)

    def reset(self, residual=None) -> None:
        if self._bound:
            raise ContractError("reset is rejected while a policy is bound")
        baseline = self._check_baseline(residual)
        self._retired.update(self._records)
        self._records.clear()
        self._orphans.clear()
        self._baseline = baseline

    def copy(self) -> FlowGraph:
        if self._bound:
            raise ContractError("a ledger with bound policies cannot be copied")
        other = FlowGraph.__new__(FlowGraph)
        other.graph = self.graph
        other.uid = next(_UID)
        other._baseline = list(self._baseline)
        other._records = {
            i: replace(r.snapshot(), flow_id=(other.uid, i))
            for i, r in self._records.items()
        }
        other._next_id = self._next_id
        other._next_owner = self._next_owner
        other._bound = set()
        other._orphans = set(self._orphans)
        other._retired = set(self._retired)
        return other

    # ------------------------------------------------------------ policy-private

    def _bind(self) -> int:
        owner = self._next_owner
        self._next_owner += 1
        self._bound.add(owner)
        return owner

    def _orphan(self, owner: int) -> None:
        """The policy was destroyed without ``remove()``: its records stay."""
        self._bound.discard(owner)
        self._orphans.add(owner)

    def _detach(self, owner: int) -> list[int]:
        """``remove()``: revert the owner's allocations and retire its ids."""
        ids = [i for i, r in self._records.items() if r.owner == owner]
        for i in ids:
            del self._records[i]
            self._retired.add(i)
        return ids

    def _commit(self, owner: int, drafts: list[Draft]) -> list[tuple[int, int]]:
        """One atomic exchange of the owner's records; returns the new FlowIds."""
        if owner not in self._bound:
            raise ContractError("commit by an unbound policy")
        g = self.graph
        current = {i for i, r in self._records.items() if r.owner == owner}
        kept = {d.flow_id[1] for d in drafts if d.flow_id is not None}
        if kept != current:
            raise ContractError("a commit must carry every live record of the policy")
        alloc = [0] * g.num_edges
        for r in self._records.values():
            if r.owner != owner:
                for e, x in r.edges.items():
                    alloc[e] += x
        for d in drafts:
            for e, x in d.edges.items():
                if x < 0:
                    raise AssertionError("negative allocation in a commit")
                alloc[e] += x
            for x in d.drops.values():
                if not 0 <= x < COUNTER_LIMIT:
                    raise ContractError("drop counter outside the 128-bit domain")
            if d.offered is not None and not 0 <= d.offered < COUNTER_LIMIT:
                raise ContractError("offered counter outside the 128-bit domain")
        for e in range(g.num_edges):
            if alloc[e] > self._baseline[e]:
                raise AssertionError(f"infeasible commit on edge {e}")
        new_ids: list[tuple[int, int]] = []
        for d in drafts:
            if d.flow_id is None:
                i = self._next_id
                self._next_id += 1
                fid = (self.uid, i)
                new_ids.append(fid)
            else:
                fid = d.flow_id
                i = fid[1]
            self._records[i] = FlowRecord(
                flow_id=fid,
                owner=owner,
                sources=d.sources,
                destinations=d.destinations,
                admission=d.admission,
                bundle=d.bundle,
                route_cost=d.route_cost,
                offered=d.offered,
                delivered=d.delivered,
                edges={e: x for e, x in d.edges.items() if x},
                drops={e: x for e, x in d.drops.items() if x},
            )
        return new_ids

    # ------------------------------------------------------------ state digest

    def state(self) -> dict:
        g = self.graph
        pos = g.input_index()
        alloc = self.allocated()
        res = self.residual()
        by_input_alloc = [0] * g.num_edges
        by_input_res = [0] * g.num_edges
        for e in range(g.num_edges):
            by_input_alloc[pos[e]] = alloc[e]
            by_input_res[pos[e]] = res[e]
        return {
            "ledger_uid": self.uid,
            "allocated_input": by_input_alloc,
            "residual_input": by_input_res,
            "records": [r.to_json(g) for r in self.records()],
            "bound": sorted(self._bound),
            "orphan_owners": sorted(self._orphans),
            "next_flow_id": self._next_id,
        }
