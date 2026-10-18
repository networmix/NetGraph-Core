"""The policy engine ``FlowPolicy`` in integer units (APPROACH revision 24, 8 and 9).

One placement engine for every routing mode: ``CostOnly``, ``CapacityAware``,
``Global`` and ``Pinned``, the kernels ``MaxFlow``, ``Weighted`` and
``EqualBalanced`` and the admissions ``Lossless`` and ``Lossy`` (prototype 1's
kernels). The transaction of 9.2 runs per source (9.3) for ``CostOnly`` and
``CapacityAware``: the roster is fixed at construction (eligible sources that
statically reach a destination, ascending NodeId), each source is offered the
increment of its round-robin share of the cumulative request, and the stages run on
the shared working residual ``W``; ``created`` is one counter per call; there is one
commit. ``Pinned`` and ``Global`` have one "source", the policy. ``Global`` is a
replacement kernel over the whole eligible graph with the canonical per-source
decomposition. Local outcomes, the termination reduction and the flag table follow
9.4 (revision 20: flags accumulate committed events). Epochs follow section 8.

Readings where the contract is silent are listed in ``NOTES.md`` (P-numbers).
"""

from __future__ import annotations

import copy
from dataclasses import dataclass, field
from fractions import Fraction

import _p1  # noqa: F401
import flowcore
import kernels as K
import search
from apportion import round_robin
from graph import Eligibility, check_dag_edges, check_disjoint, check_node_set
from ledger import Draft, FlowGraph
from units import (
    INF,
    K_MAX_FLOWS_PER_CALL,
    UNREACHABLE,
    ContractError,
    check_request,
)

MF, W_, EB = K.MAX_FLOW, K.WEIGHTED, K.EQUAL_BALANCED
LOSSLESS, LOSSY = K.LOSSLESS, K.LOSSY
OUTCOME_RANK = {"limit_flows": 4, "limit_diminishing": 3, "exhausted": 2, "no_route": 1}

# Test hook only: the engine constant kMaxFlowsPerCall is 10,000 (9.2); fixtures
# that exercise the call-wide cap override this module attribute (NOTES P9).
MAX_FLOWS_PER_CALL = K_MAX_FLOWS_PER_CALL


# ====================================================================== configuration


@dataclass(frozen=True)
class CostOnly:
    pass


@dataclass(frozen=True)
class Diminishing:
    window: int
    fraction: float


@dataclass(frozen=True)
class CapacityAware:
    flow_count: int = 1
    max_flow_count: int | None = None
    max_path_cost: tuple | None = None
    min_flow_bandwidth: int = 0
    reoptimize: bool = False
    diminishing: Diminishing | None = None


@dataclass(frozen=True)
class Global:
    pass


@dataclass(frozen=True)
class Pinned:
    bundles: tuple[tuple[int, ...], ...]


@dataclass(frozen=True)
class FlowPolicyConfig:
    routing: object
    placement: str
    admission: str
    selection: search.Selection


@dataclass
class Outputs:
    edge_flows: bool = False
    residual: bool = False
    bottleneck: bool = False


@dataclass
class PlacementResult:
    placed: int
    remaining: int | None
    termination: str
    routes_changed: bool
    weights_changed: bool
    drops: list[tuple[int, int]]
    cost_distribution: list[tuple[int, int]] | None = None
    edge_flows: list[int] | None = None
    residual: list[int] | None = None
    bottleneck: list[int] | None = None
    min_cut: list[int] | None = None
    local_outcomes: list[tuple[int, str]] = field(default_factory=list)
    trace: dict = field(default_factory=dict)

    def to_json(self) -> dict:
        out = {
            "placed": self.placed,
            "remaining": "inf" if self.remaining is None else self.remaining,
            "termination": self.termination,
            "routes_changed": self.routes_changed,
            "weights_changed": self.weights_changed,
            "drops": [list(x) for x in self.drops],
            "local_outcomes": [list(x) for x in self.local_outcomes],
        }
        if self.cost_distribution is not None:
            out["cost_distribution"] = [list(x) for x in self.cost_distribution]
        for name in ("edge_flows", "residual", "bottleneck", "min_cut"):
            value = getattr(self, name)
            if value is not None:
                out[name] = list(value)
        return out


# ====================================================================== working flows


@dataclass
class Flow:
    src: int
    route: tuple[int, ...]
    plan: K.Plan
    cost: int
    flow_id: tuple[int, int] | None = None
    T: int = 0
    alloc: dict[int, int] = field(default_factory=dict)
    offered: int | None = None
    delivered: int = 0
    drops: dict[int, int] = field(default_factory=dict)
    new: bool = False
    seq: int = 0
    bundle_index: int | None = None


def add_view(W: list[int], alloc: dict[int, int]) -> list[int]:
    out = list(W)
    for e, x in alloc.items():
        out[e] += x
    return out


def sub_view(V: list[int], alloc: dict[int, int]) -> list[int]:
    out = list(V)
    for e, x in alloc.items():
        out[e] -= x
        assert out[e] >= 0, ("overfilled view", e)
    return out


def merge(allocs) -> dict[int, int]:
    out: dict[int, int] = {}
    for a in allocs:
        for e, x in a.items():
            out[e] = out.get(e, 0) + x
    return out


@dataclass
class SourceState:
    """Per-source scope of one call (9.2 binding table)."""

    src: int
    offer: int | None
    P: int
    B: int | None
    T_entry: int
    flows: list[Flow]
    created_here: int = 0
    stage1_route: bool = True
    capped: bool = False
    limit: str | None = None
    gains: list[int] = field(default_factory=list)
    trace: dict = field(default_factory=dict)


class FlowPolicy:
    def __init__(
        self,
        fg: FlowGraph,
        sources,
        destinations,
        cfg: FlowPolicyConfig,
        elig: Eligibility | None = None,
        cache=None,
    ) -> None:
        g = fg.graph
        self.fg, self.g, self.cfg = fg, g, cfg
        if cache is not None:
            raise ContractError("DagCache is not modelled by prototype 2")
        if elig is not None:
            elig.check(g)
        self.elig = elig
        self.sources = check_node_set(g, sources, "sources")
        self.destinations = check_node_set(g, destinations, "destinations")
        check_disjoint(self.sources, self.destinations)
        self._validate_cfg()
        self.mode = type(cfg.routing).__name__
        self.kernel, self.admission = cfg.placement, cfg.admission
        self.selection = cfg.selection
        static = search.member_graph(g, elig, search.Selection())
        self._static = static
        node_ok = (lambda v: True) if elig is None else (lambda v: elig.node_allowed[v])
        self._dests_ok = tuple(d for d in self.destinations if node_ok(d))
        self.baseline: dict[int, int | None] = {}
        if self.mode == "Pinned":
            self.bundles = self._validate_bundles()
            self.usable = [(i, route) for i, route in enumerate(self._pruned) if route]
            self.roster = [self.sources[0]] if self.usable else []
        else:
            roster = []
            for s in self.sources:
                if not node_ok(s):
                    continue
                d = search.distances(g, static, [s])
                if any(d[t] is not None for t in self._dests_ok):
                    roster.append(s)
                    self.baseline[s] = min(
                        d[t] for t in self._dests_ok if d[t] is not None
                    )
            self.roster = roster
        self.owner = fg._bind()
        self._reset_epoch()

    # ------------------------------------------------------------ validation

    def _validate_cfg(self) -> None:
        cfg = self.cfg
        r = cfg.routing
        if not isinstance(r, (CostOnly, CapacityAware, Global, Pinned)):
            raise ContractError("unknown routing variant")
        if not isinstance(cfg.selection, search.Selection):
            raise ContractError("selection is explicit (no default)")
        K.check_kernel(cfg.placement, cfg.admission, cfg.selection.priority)
        if cfg.admission == LOSSY and isinstance(r, (Global, CapacityAware)):
            raise ContractError("Lossy is rejected with Global and CapacityAware")
        if isinstance(r, Global) and cfg.placement != MF:
            raise ContractError("Global requires placement = MaxFlow")
        if isinstance(r, CapacityAware):
            if not isinstance(r.flow_count, int) or r.flow_count < 1:
                raise ContractError("flow_count >= 1")
            if r.max_flow_count is not None and r.max_flow_count < r.flow_count:
                raise ContractError("max_flow_count >= flow_count")
            if not isinstance(r.min_flow_bandwidth, int) or r.min_flow_bandwidth < 0:
                raise ContractError("min_flow_bandwidth is a non-negative integer")
            search.check_ceiling(r.max_path_cost)
            if r.diminishing is not None:
                d = r.diminishing
                if not isinstance(d.window, int) or d.window < 1:
                    raise ContractError("diminishing.window >= 1")
                if not 0 < d.fraction < 1:
                    raise ContractError("diminishing.fraction in (0, 1)")
        if isinstance(r, Pinned):
            if len(self.sources) != 1 or len(self.destinations) != 1:
                raise ContractError("Pinned requires single endpoints")

    def _validate_bundles(self):
        g = self.g
        bundles = self.cfg.routing.bundles
        if not bundles:
            raise ContractError("an empty bundle list is rejected")
        src, dst = self.sources[0], self.destinations[0]
        seen = set()
        out = []
        self._pruned = []
        for b in bundles:
            edges = check_dag_edges(g, b)
            if edges in seen:
                raise ContractError("duplicate bundle (pre-prune)")
            seen.add(edges)
            if not K.prune(g, edges, src, (dst,)):
                raise ContractError("malformed bundle: no src -> dst walk")
            out.append(edges)
            ok = [e for e in edges if self.elig is None or self.elig.edge_allowed[e]]
            self._pruned.append(K.prune(g, ok, src, (dst,)))
        return out

    # ------------------------------------------------------------ state

    def _reset_epoch(self) -> None:
        self.cohorts: dict[int, list[Flow]] = {}
        self.requested = 0
        self.infinite = False
        self.epoch_started = False

    def placed_demand(self) -> int:
        return sum(f.delivered for c in self.cohorts.values() for f in c)

    def records(self):
        return [r for r in self.fg.records() if r.owner == self.owner]

    def record(self, flow_id):
        return self.fg.record(flow_id)

    def usable_bundles(self) -> int:
        return len(self.usable) if self.mode == "Pinned" else 0

    def flow_ids(self) -> list[tuple[int, int]]:
        return sorted(f.flow_id for c in self.cohorts.values() for f in c)

    def state(self) -> dict:
        pos = self.g.input_index()
        return {
            "mode": self.mode,
            "roster": list(self.roster),
            "requested": self.requested,
            "infinite": self.infinite,
            "placed_demand": self.placed_demand(),
            "cohorts": {
                s: [
                    {
                        "flow_id": None if f.flow_id is None else list(f.flow_id),
                        "route": list(f.route),
                        "route_input": sorted(pos[e] for e in f.route),
                        "T": f.T,
                        "delivered": f.delivered,
                        "alloc": sorted(f.alloc.items()),
                        "plan": f.plan.to_json(),
                        "bundle_index": f.bundle_index,
                    }
                    for f in c
                ]
                for s, c in sorted(self.cohorts.items())
            },
        }

    def remove(self) -> None:
        """Detach every flow (allocations reverted, ids retired), new epoch."""
        self.fg._detach(self.owner)
        self._reset_epoch()

    def destroy(self) -> None:
        self.fg._orphan(self.owner)

    # ------------------------------------------------------------ routes

    def _route_cost(self, route, src) -> int:
        c = K.route_cost(self.g, route, src, self._dests_ok)
        return UNREACHABLE if c is None else c

    def _new_flow(self, src, route, view, ctx) -> Flow:
        plan = K.make_plan(
            self.g,
            self.kernel,
            self.admission,
            src,
            self._dests_ok,
            route,
            self.selection.priority,
            view,
        )
        ctx["seq"] += 1
        return Flow(
            src=src,
            route=plan.route,
            plan=plan,
            cost=self._route_cost(plan.route, src),
            offered=0 if self.admission == LOSSY else None,
            new=True,
            seq=ctx["seq"],
        )

    def _cost_only_route(self, src, W) -> tuple[int, ...]:
        """Static support restricted to the source, selection projected on W."""
        g, sel = self.g, self.selection
        _, support = search.tight_dag(g, self._static, [src], self._dests_ok)
        members = support
        if not sel.multi_edge:
            members = search.collapse(g, support, sel, W)
        if sel.multipath:
            _, dag = search.tight_dag(g, members, [src], self._dests_ok)
            return tuple(dag)
        prio = (
            search.priority_fn(g, sel.priority, W)
            if sel.tie_break == search.PREFER
            else None
        )
        p = search.canonical_path(g, members, [src], self._dests_ok, prio)
        return tuple() if p is None else tuple(sorted(p.edges))

    def _search(self, src, view):
        """CapacityAware route search on a view (residual gating, ceiling)."""
        g, sel, r = self.g, self.selection, self.cfg.routing
        threshold = max(1, r.min_flow_bandwidth)
        members = search.member_graph(g, self.elig, sel, view, threshold)
        if sel.multipath:
            dstar, dag = search.tight_dag(g, members, [src], self._dests_ok)
            if dstar is None:
                return None
            route, cost = tuple(dag), dstar
        else:
            prio = (
                search.priority_fn(g, sel.priority, view)
                if sel.tie_break == search.PREFER
                else None
            )
            p = search.canonical_path(g, members, [src], self._dests_ok, prio)
            if p is None:
                return None
            route, cost = tuple(sorted(p.edges)), p.cost
        if not search.within_ceiling(cost, r.max_path_cost, self.baseline.get(src)):
            return None
        return route

    # ------------------------------------------------------------ amounts

    def _t_max(self, flow: Flow, view, evals) -> int:
        if self.kernel == MF:
            sol = K.maxflow(self.g, flow.plan, view, None)
            return sol.placed
        tm = K.t_max(self.g, flow.plan, view)
        evals.append(tm.evaluations)
        return tm.value

    def _alloc(self, flow: Flow, view, T) -> dict[int, int]:
        if self.kernel == MF:
            sol = K.maxflow(self.g, flow.plan, view, T)
            assert sol.placed == T, "MaxFlow fill below its target"
            return {e: x for e, x in enumerate(sol.edge_flows) if x}
        return K.alloc(self.g, flow.plan, T)

    def _cohort_alloc(self, flows: list[Flow], T: int) -> list[tuple[int, dict]]:
        parts = round_robin(T, len(flows))
        return [
            (t, K.alloc(self.g, f.plan, t)) for f, t in zip(flows, parts, strict=True)
        ]

    def _cohort_t_max(self, flows: list[Flow], view, evals) -> int:
        src_members = set()
        for f in flows:
            src_members.update(f.plan.out.get(f.src, []))
        hi = sum(view[e] for e in src_members)
        lo, n = 0, 0
        while lo < hi:
            mid = (lo + hi + 1) // 2
            n += 1
            total = merge(a for _, a in self._cohort_alloc(flows, mid))
            if K.fits(total, view):
                lo = mid
            else:
                hi = mid - 1
        evals.append(n)
        return lo

    def _set_cohort(self, flows: list[Flow], T: int) -> None:
        for f, (t, a) in zip(flows, self._cohort_alloc(flows, T), strict=True):
            if t != f.T or f.new:
                f.alloc = a
            f.T = t
            f.delivered = t

    # ------------------------------------------------------------ place

    def place(self, volume, outputs: Outputs | None = None) -> PlacementResult:
        req = check_request(volume)
        if self.infinite and req is not None:
            # validation precedes the zero shortcut (revision 20, section 8)
            raise ContractError("a finite request after +inf is rejected (rule 13)")
        if req is not None and self.requested + req >= 1 << 127:
            raise ContractError("cumulative request counter overflow")
        outputs = outputs or Outputs()
        if req == 0:
            served = [(src, "complete") for src in self.roster]
            return self._finish(
                PlacementResult(
                    0, 0, "complete", False, False, [], local_outcomes=served
                ),
                outputs,
                None,
            )
        if not self.roster:
            # the epoch's finite-or-infinite state and `requested` record every positive
            # request, also one that finds no roster (section 8)
            if req is None:
                self.infinite = True
            else:
                self.requested += req
            return self._finish(
                PlacementResult(
                    0, req, "no_route", False, False, [], trace={"roster": []}
                ),
                outputs,
                None,
            )
        if self.mode == "Global":
            return self._place_global(req, outputs)
        return self._place_flows(req, outputs)

    def _finish(
        self, result: PlacementResult, outputs: Outputs, sol
    ) -> PlacementResult:
        if outputs.edge_flows:
            alloc = [0] * self.g.num_edges
            for r in self.records():
                for e, x in r.edges.items():
                    alloc[e] += x
            result.edge_flows = alloc
        if outputs.residual:
            result.residual = self.fg.residual()
        if outputs.bottleneck:
            result.bottleneck = self.fg.bottleneck(self.sources, self.elig)
            if self.mode == "Global":
                result.min_cut = None if sol is None else sol.min_cut
        return result

    # ------------------------------------------------------------ Global

    def _place_global(self, req, outputs) -> PlacementResult:
        g, fg = self.g, self.fg
        W = fg.residual()
        old = self.cohorts
        created = not old
        own = merge(f.alloc for c in old.values() for f in c)
        V = add_view(W, own)
        P = self.placed_demand()
        cap = None if req is None else P + req
        sol = flowcore.solve(g, V, self.roster, self._dests_ok, cap, elig=self.elig)
        assert sol.placed >= P, "a Global call never reduces placed_demand()"
        shares = self._decompose(sol)
        changed = created
        new_cohorts: dict[int, list[Flow]] = {}
        for s in self.roster:
            alloc = shares[s]
            route = tuple(sorted(alloc))
            prev = old.get(s, [None])[0]
            if prev is not None and prev.route != route:
                changed = True
            plan = K.Plan(MF, LOSSLESS, s, self._dests_ok, route, route, None)
            f = Flow(
                src=s,
                route=route,
                plan=plan,
                cost=self._route_cost(route, s),
                flow_id=None if prev is None else prev.flow_id,
                T=sol.supply.get(s, 0),
                alloc=alloc,
                delivered=sol.supply.get(s, 0),
                new=prev is None,
                seq=self.roster.index(s),
            )
            new_cohorts[s] = [f]
        placed = sol.placed - P
        drafts, order = self._drafts(new_cohorts)
        new_ids = fg._commit(self.owner, drafts)
        self._assign_ids(order, new_ids)
        self.cohorts = new_cohorts
        if req is None:
            self.infinite = True
        else:
            self.requested += req
        if req is not None and placed == req:
            term = "complete"
        elif sol.placed == 0:
            term = "no_route"
        else:
            term = "exhausted"
        trace = {
            "V": V,
            "cap": cap,
            "flow_core": {**sol.payload(), "supply": sol.supply},
            "shares": {s: sorted(a.items()) for s, a in shares.items()},
        }
        result = PlacementResult(
            placed,
            None if req is None else req - placed,
            term,
            changed,
            False,
            [],
            cost_distribution=list(sol.cost_distribution),
            local_outcomes=[(src, term) for src in self.roster],
            trace=trace,
        )
        return self._finish(result, outputs, sol)

    def _decompose(self, sol) -> dict[int, dict[int, int]]:
        """Canonical per-source decomposition (9.1): sources ascending; the smallest
        EdgeId walk with positive remaining flow to the first destination reached,
        saturated by min(walk bottleneck, remaining virtual supply)."""
        g = self.g
        remaining = dict(enumerate(sol.edge_flows))
        dests = set(self._dests_ok)
        out: dict[int, dict[int, int]] = {}
        for s in self.roster:
            mine: dict[int, int] = {}
            left = sol.supply.get(s, 0)
            while left > 0:
                u, walk = s, []
                while not (walk and u in dests):
                    cands = [e for e in g.out_ids(u) if remaining.get(e, 0) > 0]
                    assert cands, "flow conservation guarantees a continuation"
                    e = min(cands)
                    walk.append(e)
                    u = g.dst[e]
                    assert len(walk) <= g.num_nodes, "acyclic support"
                amount = min(min(remaining[e] for e in walk), left)
                for e in walk:
                    remaining[e] -= amount
                    mine[e] = mine.get(e, 0) + amount
                left -= amount
            out[s] = mine
        assert all(x == 0 for x in remaining.values()), "decomposition covers the flow"
        return out

    # ------------------------------------------------------------ route-search modes

    def _offers(self, req) -> list[int | None]:
        k = len(self.roster)
        if self.mode == "Pinned":
            return [req]
        if req is None:
            return [None] * k
        old = round_robin(self.requested, k)
        new = round_robin(self.requested + req, k)
        return [b - a for a, b in zip(old, new, strict=True)]

    def _place_flows(self, req, outputs) -> PlacementResult:
        fg = self.fg
        ctx = {
            "W": fg.residual(),
            "created": 0,
            "routes_changed": False,
            "weights_changed": False,
            "seq": 0,
            "drops": {},
        }
        cohorts = copy.deepcopy(self.cohorts)
        entry = {s: sum(f.delivered for f in cohorts.get(s, [])) for s in self.roster}
        offers = self._offers(req)
        states = []
        for s, offer in zip(self.roster, offers, strict=True):
            flows = cohorts.setdefault(s, [])
            st = SourceState(
                src=s,
                offer=offer,
                P=entry[s],
                B=offer,
                T_entry=sum(f.T for f in flows),
                flows=flows,
            )
            self._run_source(st, ctx)
            states.append(st)
        drafts, order = self._drafts(cohorts)
        new_ids = fg._commit(self.owner, drafts)
        self._assign_ids(order, new_ids)
        self.cohorts = {s: c for s, c in cohorts.items() if c}
        if req is None:
            self.infinite = True
        else:
            self.requested += req
        placed = sum(
            sum(f.delivered for f in cohorts.get(s, [])) - entry[s] for s in self.roster
        )
        assert placed >= 0
        if req is not None:
            assert placed <= req
        outcomes = {st.src: self._local_outcome(st) for st in states}
        unserved = [o for o in outcomes.values() if o != "complete"]
        term = max(unserved, key=OUTCOME_RANK.__getitem__) if unserved else "complete"
        trace = {
            "offers": {
                st.src: "inf" if st.offer is None else st.offer for st in states
            },
            "local_outcomes": outcomes,
            "created": ctx["created"],
            "sources": {st.src: st.trace for st in states},
            "W_final": ctx["W"],
        }
        result = PlacementResult(
            placed,
            None if req is None else req - placed,
            term,
            ctx["routes_changed"],
            ctx["weights_changed"],
            sorted(ctx["drops"].items()),
            local_outcomes=[(st.src, outcomes[st.src]) for st in states],
            trace=trace,
        )
        return self._finish(result, outputs, None)

    def _local_outcome(self, st: SourceState) -> str:
        delivered = sum(f.delivered for f in st.flows) - st.P
        if st.offer == 0:
            return "complete"
        if self.admission == LOSSY:
            if st.offer is not None and st.trace.get("injected") == st.offer:
                return "complete"
        elif st.offer is not None and delivered == st.offer:
            return "complete"
        if st.limit == "limit_flows" or st.capped:
            return "limit_flows"
        if st.limit == "limit_diminishing":
            return "limit_diminishing"
        if not st.flows:
            return "no_route"
        return "exhausted"

    def _run_source(self, st: SourceState, ctx) -> None:
        if not st.flows:
            self._stage1(st, ctx)
        if not st.flows:
            return
        if self.mode == "CapacityAware" and self.cfg.routing.reoptimize:
            self._stage2(st, ctx)
        if self.admission == LOSSY:
            self._stage3_lossy(st, ctx)
        elif self.kernel == EB:
            self._stage3_eb(st, ctx)
        else:
            self._stage3_seq(st, ctx)
        if self.mode == "CapacityAware":
            self._stage4(st, ctx)

    # ---- stage 1

    def _stage1(self, st: SourceState, ctx) -> None:
        W = ctx["W"]
        tr: dict = {}
        if self.mode == "CostOnly":
            route = self._cost_only_route(st.src, W)
            tr["route"] = list(route)
            st.flows.append(self._new_flow(st.src, route, W, ctx))
        elif self.mode == "Pinned":
            for i, route in self.usable:
                f = self._new_flow(st.src, route, W, ctx)
                f.bundle_index = i
                st.flows.append(f)
            tr["usable"] = [i for i, _ in self.usable]
        else:
            r = self.cfg.routing
            S = list(W)
            rem = st.offer
            target = None
            if self.kernel == EB and st.offer is not None:
                target = (st.P + st.offer) // r.flow_count
            tr.update(eb_target=target, reservations=[])
            budget = min(r.flow_count, MAX_FLOWS_PER_CALL - ctx["created"])
            stopped = False
            for _ in range(budget):
                route = self._search(st.src, S)
                if route is None:
                    stopped = True
                    break
                f = self._new_flow(st.src, route, S, ctx)
                evals: list[int] = []
                bound = self._t_max(f, S, evals)
                if self.kernel == EB:
                    t = bound if target is None else min(target, bound)
                else:
                    t = bound if rem is None else min(rem, bound)
                    if rem is not None:
                        rem -= t
                a = self._alloc(f, S, t) if t else {}
                tr["reservations"].append(
                    {
                        "route": list(route),
                        "S": list(S),
                        "T_bound": bound,
                        "reserved": t,
                    }
                )
                S = sub_view(S, a)
                st.flows.append(f)
                ctx["created"] += 1
                st.created_here += 1
            # the call-wide cap blocked (part of) the creation, not a missing route
            st.capped = budget < r.flow_count and not stopped
            tr["created_cap_blocked"] = st.capped
        if st.flows:
            ctx["routes_changed"] = True
        st.trace["stage1"] = tr

    # ---- stage 2

    def _stage2(self, st: SourceState, ctx) -> None:
        rows = []
        for f in st.flows:
            if f.new:
                continue
            W_i = add_view(ctx["W"], f.alloc)
            cand = self._search(st.src, W_i)
            row = {
                "flow": list(f.flow_id),
                "T": f.T,
                "candidate": None if cand is None else list(cand),
            }
            if cand is None or cand == f.route:
                row["decision"] = "kept"
                rows.append(row)
                continue
            trial = self._new_flow(st.src, cand, W_i, {"seq": 0})
            evals: list[int] = []
            bound = self._t_max(trial, W_i, evals)
            if bound < f.T:
                row.update(decision="rejected: delivers less than T_i", bound=bound)
                rows.append(row)
                continue
            a = self._alloc(trial, W_i, f.T) if f.T else {}
            ctx["W"] = sub_view(W_i, a)
            f.route, f.plan, f.cost, f.alloc = trial.route, trial.plan, trial.cost, a
            ctx["routes_changed"] = True
            if self.kernel == W_:
                ctx["weights_changed"] = True
            row.update(decision="adopted", alloc=sorted(a.items()))
            rows.append(row)
        st.trace["stage2"] = rows

    # ---- stage 3

    def _stage3_seq(self, st: SourceState, ctx) -> None:
        rows = []
        for f in st.flows:
            V = add_view(ctx["W"], f.alloc)
            evals: list[int] = []
            if self.kernel == MF:
                cap = None if st.B is None else f.T + st.B
                sol = K.maxflow(self.g, f.plan, V, cap)
                t_new = sol.placed
                assert t_new >= f.T
                f.alloc = {e: x for e, x in enumerate(sol.edge_flows) if x}
                curve = sol.cost_distribution
            else:
                bound = self._t_max(f, V, evals)
                target = bound if st.B is None else min(f.T + st.B, bound)
                t_new = max(f.T, target)
                if t_new != f.T or f.new:
                    f.alloc = K.alloc(self.g, f.plan, t_new)
                curve = None
            ctx["W"] = sub_view(V, f.alloc)
            gain = t_new - f.T
            rows.append(
                {
                    "T": [f.T, t_new],
                    "evals": evals,
                    "curve": curve,
                    "alloc": sorted(f.alloc.items()),
                }
            )
            f.T = t_new
            f.delivered = t_new
            if st.B is not None:
                st.B -= gain
        st.trace["stage3"] = rows

    def _stage3_eb(self, st: SourceState, ctx) -> None:
        flows = st.flows
        current = merge(f.alloc for f in flows)
        V = add_view(ctx["W"], current)
        T = sum(f.T for f in flows)
        evals: list[int] = []
        bound = self._cohort_t_max(flows, V, evals)
        target = bound if st.B is None else min(T + st.B, bound)
        t_new = max(T, target)
        self._set_cohort(flows, t_new)
        ctx["W"] = sub_view(V, merge(f.alloc for f in flows))
        if st.B is not None:
            st.B -= t_new - T
        st.trace["stage3"] = {
            "T": [T, t_new],
            "bound": bound,
            "evals": evals,
            "shares": [f.T for f in flows],
        }

    def _stage3_lossy(self, st: SourceState, ctx) -> None:
        flows = st.flows
        offers = (
            [None] * len(flows)
            if st.offer is None
            else round_robin(st.offer, len(flows))
        )
        injected = 0
        rows = []
        for f, offer in zip(flows, offers, strict=True):
            call = K.lossy(self.g, f.plan, offer, ctx["W"])
            ctx["W"] = sub_view(ctx["W"], call.carried)
            for e, x in call.carried.items():
                f.alloc[e] = f.alloc.get(e, 0) + x
            for e, x in call.drops.items():
                f.drops[e] = f.drops.get(e, 0) + x
                ctx["drops"][e] = ctx["drops"].get(e, 0) + x
            f.offered = (f.offered or 0) + call.offered
            f.delivered += call.delivered
            injected += call.injected
            rows.append(
                {
                    "offer": "inf" if offer is None else offer,
                    "offered": call.offered,
                    "delivered": call.delivered,
                    "carried": sorted(call.carried.items()),
                    "drops": sorted(call.drops.items()),
                }
            )
        st.trace["stage3"] = rows
        st.trace["injected"] = injected

    # ---- stage 4

    def _stage4(self, st: SourceState, ctx) -> None:
        r = self.cfg.routing
        steps = []
        v = st.offer
        threshold = None
        if r.diminishing is not None and v is not None:
            threshold = Fraction(r.diminishing.fraction) * v
        accepted = 0
        while st.B is None or st.B > 0:
            if ctx["created"] >= MAX_FLOWS_PER_CALL:
                st.limit = "limit_flows"
                steps.append({"stop": "kMaxFlowsPerCall"})
                break
            if r.max_flow_count is not None and len(st.flows) >= r.max_flow_count:
                if r.max_flow_count > r.flow_count:
                    st.limit = "limit_flows"
                steps.append({"stop": "max_flow_count", "k": len(st.flows)})
                break
            cand = self._search(st.src, ctx["W"])
            if cand is None:
                steps.append({"stop": "no candidate"})
                break
            if self.selection.multipath and any(f.route == cand for f in st.flows):
                # revision 24 (owner review): with multipath routes a DAG already spans
                # every equal-cost path, so a second flow on a structurally equal DAG
                # gains only round-robin remainders: growth stops. A repeated single path
                # is a legitimate tunnel and stays a growth candidate.
                steps.append(
                    {"stop": "candidate equals a cohort route", "route": list(cand)}
                )
                break
            f = self._new_flow(st.src, cand, ctx["W"], ctx)
            step: dict = {"route": list(cand), "W": list(ctx["W"])}
            evals: list[int] = []
            if self.kernel == EB:
                current = merge(x.alloc for x in st.flows)
                V = add_view(ctx["W"], current)
                T = sum(x.T for x in st.flows)
                trial = st.flows + [f]
                bound = self._cohort_t_max(trial, V, evals)
                t_cand = bound if v is None else min(st.T_entry + v, bound)
                step.update(T=T, T_candidate=t_cand, bound=bound)
                if t_cand <= T:
                    step["decision"] = "rejected: cohort total not larger"
                    ctx["seq"] -= 1
                    steps.append(step)
                    break
                st.flows.append(f)
                self._set_cohort(st.flows, t_cand)
                ctx["W"] = sub_view(V, merge(x.alloc for x in st.flows))
                gain = t_cand - T
                step["shares"] = [x.T for x in st.flows]
            else:
                bound = self._t_max(f, ctx["W"], evals)
                fill = bound if st.B is None else min(st.B, bound)
                step.update(bound=bound, fill=fill)
                if fill == 0:
                    step["decision"] = "rejected: zero fill"
                    ctx["seq"] -= 1
                    steps.append(step)
                    break
                f.alloc = self._alloc(f, ctx["W"], fill)
                f.T = f.delivered = fill
                ctx["W"] = sub_view(ctx["W"], f.alloc)
                st.flows.append(f)
                gain = fill
            ctx["created"] += 1
            ctx["routes_changed"] = True
            accepted += 1
            st.gains.append(gain)
            if st.B is not None:
                st.B -= gain
            step.update(decision="accepted", gain=gain, evals=evals)
            steps.append(step)
            if threshold is not None and accepted >= r.diminishing.window:
                recent = sum(st.gains[-r.diminishing.window :])
                if recent < threshold and st.B > 0:
                    st.limit = "limit_diminishing"
                    steps.append(
                        {
                            "stop": "diminishing",
                            "recent": recent,
                            "threshold": str(threshold),
                        }
                    )
                    break
        st.trace["stage4"] = steps

    # ------------------------------------------------------------ commit

    def _drafts(self, cohorts) -> tuple[list[Draft], list[Flow]]:
        drafts: list[Draft] = []
        order: list[Flow] = []
        flows = [f for s in sorted(cohorts) for f in cohorts[s]]
        old = [f for f in flows if not f.new]
        new = sorted(
            (f for f in flows if f.new),
            key=lambda f: (f.seq if self.mode != "Global" else f.src),
        )
        for f in old + new:
            drafts.append(
                Draft(
                    flow_id=None if f.new else f.flow_id,
                    sources=(f.src,),
                    destinations=self._dests_ok,
                    placement=self.kernel,
                    admission=self.admission,
                    bundle=f.route,
                    route_cost=f.cost,
                    offered=f.offered if self.admission == LOSSY else None,
                    delivered=f.delivered,
                    edges=dict(f.alloc),
                    drops=dict(f.drops),
                )
            )
            order.append(f)
        return drafts, [f for f in order if f.new]

    def _assign_ids(self, new_flows: list[Flow], new_ids) -> None:
        for f, fid in zip(new_flows, new_ids, strict=True):
            f.flow_id = fid
            f.new = False


def capacity(g, sources, destinations, cfg, elig=None, outputs=None) -> PlacementResult:
    """``capacity`` = ``place(+inf)`` on a fresh ledger (section 11)."""
    fg = FlowGraph(g)
    return FlowPolicy(fg, sources, destinations, cfg, elig).place(INF, outputs)
