"""The policy engine ``FlowPolicy`` in integer units (APPROACH revision 24, 8 and 9).

One placement engine for the routing modes ``CostOnly`` and ``Global`` (with or without
a support), the kernels ``MaxFlow``, ``Weighted`` and ``EqualBalanced`` and the admissions
``Lossless`` and ``Lossy``. The roster is fixed at construction (eligible sources that
statically reach a destination, ascending NodeId); there is one commit per call. ``Global``
is a replacement kernel over the whole eligible graph (or exactly over its support) with the
canonical per-source decomposition. Epochs follow section 8.

Step 2 of the next Core change (contract r8, sections 2.2, 4.1, 5.1, 9.1, 10.2, 14.1,
14.2; README "Pooled sources and the origin"): ``CostOnly`` pools its sources. Its
routes are prepared once per epoch on the residual at the first positive call's entry
and joined into one program (``make_program``); EqualBalanced and Weighted forward the
equal origin over it (``forward``: one split per router, attribution by highest
averages); MaxFlow and Global are controllers (``_place_controller``) with the optimized
origin. The result has the 14.1 shape (``PlacementResult``) and the 14.2 termination.

Step 3 (contract r8, sections 4.1, 5.2, 8, 10.4; README "Programs and Global with a
support"): ``Pinned`` is gone from the reference (its programs are checked by the library's
own tests); ``Global(support)`` is the controller restricted to a supplied DAG, trimmed to the
walks from the sources to the allowed destinations; CostOnly Weighted takes its weights from a
weight basis (``weights``: capacity, or the residual snapshot at preparation) instead of the
selection's priority; the change flags of CostOnly and the controllers compare programs (a
CostOnly program is created at the epoch's first positive call; a controller's program is its
solve cut into paths, compared path by path, with the optimized origin's supplies).

When T0 was deleted (the contract's section 6; README "T1 replaces CapacityAware"), the per-source
``CapacityAware`` transaction is gone from the reference with the library's (its tunnel
form, T1, is checked by the library's own tests).

Readings where the contract is silent are listed in ``NOTES.md`` (P-numbers).
"""

from __future__ import annotations

import copy
from dataclasses import dataclass, field

import _p1  # noqa: F401
import flowcore
import kernels as K
import search
from apportion import highest_averages, round_robin, split
from graph import Eligibility, check_dag_edges, check_disjoint, check_node_set
from ledger import Draft, FlowGraph
from units import (
    INF,
    UNREACHABLE,
    ContractError,
    check_request,
)

MF, W_, EB = K.MAX_FLOW, K.WEIGHTED, K.EQUAL_BALANCED
LOSSLESS, LOSSY = K.LOSSLESS, K.LOSSY


# ====================================================================== configuration


@dataclass(frozen=True)
class CostOnly:
    pass


@dataclass(frozen=True)
class Global:
    support: tuple[int, ...] | None = None


@dataclass(frozen=True)
class FlowPolicyConfig:
    """Global takes no placement and no selection (both None); the origin is the model's
    default (step 2: EQUAL for CostOnly EqualBalanced and Weighted, OPTIMIZED for CostOnly
    MaxFlow and Global; the equal origin of the controllers is not modelled here)."""

    routing: object
    placement: str | None
    admission: str
    selection: search.Selection | None
    weights: str | None = (
        None  # CostOnly Weighted: "Capacity" (None) or "ResidualSnapshot"
    )


@dataclass
class Outputs:
    edge_flows: bool = False
    residual: bool = False
    bottleneck: bool = False
    min_cut: bool = False


@dataclass
class PlacementResult:
    """The result of one call (step 2 shape): this call's values, then the demand's state
    after it (epoch totals), then the optional outputs."""

    placed: int
    remaining: int | None
    termination: str
    routes_changed: bool
    weights_changed: bool
    call_drops: list[tuple[int, int]]
    call_source_loss: int = 0
    requested: int | None = 0  # None: +inf
    offered: int | None = None  # Lossy only
    admitted: int = 0
    delivered: int = 0
    unadmitted: int | None = 0
    source_loss: int = 0
    drops: list[tuple[int, int]] = field(default_factory=list)
    by_source: list[tuple] = field(default_factory=list)
    total_cost: int | None = None
    cost_distribution: list[tuple[int, int]] | None = None
    edge_flows: list[int] | None = None
    residual: list[int] | None = None
    bottleneck: list[int] | None = None
    min_cut: list[int] | None = None
    trace: dict = field(default_factory=dict)

    def to_json(self) -> dict:
        out = {
            "placed": self.placed,
            "remaining": "inf" if self.remaining is None else self.remaining,
            "delivery_delta": self.placed,
            "call_drops": [list(x) for x in self.call_drops],
            "call_source_loss": self.call_source_loss,
            "termination": self.termination,
            "routes_changed": self.routes_changed,
            "weights_changed": self.weights_changed,
            "requested": "inf" if self.requested is None else self.requested,
            "admitted": self.admitted,
            "delivered": self.delivered,
            "source_loss": self.source_loss,
            "drops": [list(x) for x in self.drops],
            "by_source": [list(x) for x in self.by_source],
        }
        for name in ("offered", "unadmitted", "total_cost"):
            value = getattr(self, name)
            if value is not None:
                out[name] = value
        if self.cost_distribution is not None:
            out["cost_distribution"] = [list(x) for x in self.cost_distribution]
        for name in ("edge_flows", "residual", "bottleneck", "min_cut"):
            value = getattr(self, name)
            if value is not None:
                out[name] = list(value)
        return out


def ha(total: int, weights: list[int]) -> list[int]:
    """Highest averages with zero weights allowed: a zero weight gets nothing, all-zero
    weights give nothing (step 2: the attribution of pooled traffic to its sources)."""
    out = [0] * len(weights)
    pos = [i for i, w in enumerate(weights) if w > 0]
    if pos and total:
        for i, x in zip(
            pos, highest_averages(total, [weights[i] for i in pos]), strict=True
        ):
            out[i] = x
    return out


@dataclass
class Program:
    """A CostOnly demand's pooled program (step 2, contract r8 5.1): the union of the roster's
    routes, its members (EqualBalanced: every route edge; Weighted: the edges of positive weight
    on the walks from a root to a destination), the frozen weights, per node its members in
    ascending EdgeId and a topological order of every node."""

    kernel: str
    roots: list[int]
    dsts: tuple[int, ...]
    members: tuple[int, ...]
    weights: dict[int, int] | None
    out: dict[int, list[int]]
    topo: list[int]


def make_program(g, kernel, roots, dsts, union, basis, view) -> Program:
    """Weighted weights by the weight basis (contract r8 5.2): the static capacity, or the free
    capacity `view` at preparation (ResidualSnapshot)."""
    union = tuple(sorted(union))
    snapshot = basis == "ResidualSnapshot"
    if kernel != W_:
        cand = union
    elif snapshot:
        cand = tuple(e for e in union if view[e] >= 1)
    else:
        cand = tuple(e for e in union if g.capacity[e] >= 1)
    members = tuple(sorted({e for r in roots for e in K.prune(g, cand, r, dsts)}))
    weights = None
    if kernel == W_:
        weights = {e: (view[e] if snapshot else g.capacity[e]) for e in members}
    out: dict[int, list[int]] = {}
    nodes = set(roots) | set(dsts)
    for e in members:
        out.setdefault(g.src[e], []).append(e)
        nodes.update((g.src[e], g.dst[e]))
    indeg = dict.fromkeys(nodes, 0)
    for e in members:
        indeg[g.dst[e]] += 1
    topo, ready = [], sorted(v for v in nodes if indeg[v] == 0)
    while ready:
        u = ready.pop(0)
        topo.append(u)
        for e in out.get(u, ()):
            indeg[g.dst[e]] -= 1
            if indeg[g.dst[e]] == 0:
                ready.append(g.dst[e])
                ready.sort()
    assert len(topo) == len(nodes), "the pooled program is cyclic"
    return Program(kernel, list(roots), tuple(dsts), members, weights, out, topo)


def forward(g, prog: Program, inject, W):
    """One forwarding of a pooled program: inject[i] at roots[i] (None: the first-hop
    experiment, every root's out-members each injected with exactly their W), each node's
    inflow (own injection plus carried arrivals) split once by round robin or highest
    averages; a member's amount attributed to the roots by highest averages over their
    remaining inflows (member by member), then (W given) the member carries min(arrivals, W[e])
    attributed by highest averages over the roots' arrivals (contract 5.1, 9.1). Returns per
    root: offered, source loss, delivered and {edge: carried}, {edge: dropped}; and the totals."""
    k = len(prog.roots)
    index = {r: i for i, r in enumerate(prog.roots)}
    offered, lost, delivered = [0] * k, [0] * k, [0] * k
    carried = [dict() for _ in range(k)]
    dropped = [dict() for _ in range(k)]
    inflow: dict[int, dict[int, int]] = {}
    own_at: dict[int, int] = {}
    for i, r in enumerate(prog.roots):
        if inject is None:
            if prog.out.get(r):
                own_at[r] = i
                offered[i] = sum(W[e] for e in prog.out[r])
            continue
        offered[i] = inject[i]
        if not inject[i]:
            continue
        if not prog.out.get(r):
            lost[i] = inject[i]
        else:
            inflow.setdefault(r, {})[i] = inject[i]
    dset = set(prog.dsts)
    for u in prog.topo:
        rem = dict(inflow.get(u, {}))
        if u in dset:
            for i, x in rem.items():
                delivered[i] += x
            continue
        own = own_at.get(u)
        amount = sum(rem.values())
        if amount == 0 and own is None:
            continue
        outs = prog.out.get(u, [])
        weights = None if prog.weights is None else [prog.weights[e] for e in outs]
        part = (
            split("ha" if prog.kernel == W_ else "rr", amount, weights, len(outs))
            if amount
            else [0] * len(outs)
        )
        order = sorted(set(rem) | ({own} if own is not None else set()))
        for e, x in zip(outs, part, strict=True):
            row = ha(x, [rem.get(i, 0) for i in order])
            for i, y in zip(order, row, strict=True):
                rem[i] = rem.get(i, 0) - y
            arrive = [
                y + (W[e] if i == own else 0) for i, y in zip(order, row, strict=True)
            ]
            total = sum(arrive)
            if total == 0:
                continue
            c = total if W is None else min(total, W[e])
            keep = arrive if c == total else ha(c, arrive)
            for i, a, kk in zip(order, arrive, keep, strict=True):
                if kk:
                    carried[i][e] = carried[i].get(e, 0) + kk
                    inflow.setdefault(g.dst[e], {})
                    inflow[g.dst[e]][i] = inflow[g.dst[e]].get(i, 0) + kk
                if a - kk:
                    dropped[i][e] = dropped[i].get(e, 0) + a - kk
        assert all(v == 0 for v in rem.values()), (
            "a router's inflow is attributed in full"
        )
    return offered, lost, delivered, carried, dropped, index


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


def add_view(W: list[int], alloc: dict[int, int]) -> list[int]:
    out = list(W)
    for e, x in alloc.items():
        out[e] += x
    return out


def merge(allocs) -> dict[int, int]:
    out: dict[int, int] = {}
    for a in allocs:
        for e, x in a.items():
            out[e] = out.get(e, 0) + x
    return out


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
        self.kernel = MF if self.mode == "Global" else cfg.placement
        self.support = cfg.routing.support if self.mode == "Global" else None
        self.admission = cfg.admission
        self.selection = cfg.selection or search.Selection()
        static = search.member_graph(g, elig, search.Selection())
        self._static = static
        node_ok = (lambda v: True) if elig is None else (lambda v: elig.node_allowed[v])
        self._dests_ok = tuple(d for d in self.destinations if node_ok(d))
        if self.support is not None:
            # The support trimmed to the walks from each allowed source to an allowed destination
            # over its allowed edges; the roster is the sources that reach one inside it.
            ok = [
                e
                for e in check_dag_edges(g, self.support)
                if elig is None or elig.edge_allowed[e]
            ]
            trimmed: set[int] = set()
            roster = []
            for s in self.sources:
                if not node_ok(s):
                    continue
                part = K.prune(g, ok, s, self._dests_ok)
                if part:
                    roster.append(s)
                    trimmed.update(part)
            self.roster = roster
            self.trimmed = tuple(sorted(trimmed))
        else:
            roster = []
            for s in self.sources:
                if not node_ok(s):
                    continue
                d = search.distances(g, static, [s])
                if any(d[t] is not None for t in self._dests_ok):
                    roster.append(s)
            self.roster = roster
        self.owner = fg._bind()
        self._reset_epoch()

    # ------------------------------------------------------------ validation

    def _validate_cfg(self) -> None:
        cfg = self.cfg
        r = cfg.routing
        if not isinstance(r, (CostOnly, Global)):
            raise ContractError("unknown routing variant")
        if isinstance(r, Global):
            # Global takes no placement and no selection; a controller is lossless.
            if cfg.placement is not None or cfg.selection is not None:
                raise ContractError("Global takes no placement and no selection")
            if cfg.admission != LOSSLESS:
                raise ContractError("MAX_FLOW requires admission=LOSSLESS")
        else:
            if cfg.placement is None:
                raise ContractError("placement is required")
            sel = cfg.selection or search.Selection()
            K.check_kernel(cfg.placement, cfg.admission)
            if cfg.weights is not None and cfg.placement != W_:
                raise ContractError("weights require CostOnly placement=WEIGHTED")
            if (
                len(self.sources) > 1
                and not sel.multipath
                and sel.tie_break == search.PREFER
            ):
                raise ContractError("the pooled single-path rule (Q12)")

    # ------------------------------------------------------------ state

    def _reset_epoch(self) -> None:
        self.cohorts: dict[int, list[Flow]] = {}
        self.requested = 0
        self.offered = 0
        self.infinite = False
        self.epoch_started = False
        self.prepared = False
        self.program: Program | None = None
        self.union: tuple[int, ...] = ()
        self.curve: list[tuple[int, int]] = []
        # A controller's program (10.4): its paths (source, walk) and its weights (the paths'
        # amounts, the optimized origin's supplies); None before the first solve.
        self.paths: tuple | None = None

    def placed_demand(self) -> int:
        return sum(f.delivered for c in self.cohorts.values() for f in c)

    def records(self):
        return [r for r in self.fg.records() if r.owner == self.owner]

    def record(self, flow_id):
        return self.fg.record(flow_id)

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

    # ------------------------------------------------------------ place

    def place(self, volume, outputs: Outputs | None = None) -> PlacementResult:
        req = check_request(volume)
        if self.infinite and req is not None:
            # validation precedes the zero shortcut (revision 20, section 8)
            raise ContractError("a finite request after +inf is rejected (rule 13)")
        if req is not None and self.requested + req >= 1 << 127:
            raise ContractError("cumulative request counter overflow")
        outputs = outputs or Outputs()
        controller = self.mode == "Global" or (
            self.mode == "CostOnly" and self.kernel == MF
        )
        if req == 0:
            # The current state, every flag false; before the epoch's first positive call it is
            # all zeros (contract r8 10.3). Global reports its curve once prepared.
            r = PlacementResult(0, 0, "complete", False, False, [])
            if self._dense_global() and self.prepared:
                r.cost_distribution = list(self.curve)
            return self._finish(r, outputs, None, controller)
        if not self.roster:
            # the epoch's finite-or-infinite state and `requested` record every positive
            # request, also one that finds no roster (section 8); a Lossy offer is lost before
            # it enters (the demand's source loss; every source's row stays zero)
            if req is None:
                self.infinite = True
            else:
                self.requested += req
                if self.admission == LOSSY:
                    self.offered += req
            self.prepared = True
            r = PlacementResult(
                0,
                req,
                "no_route",
                False,
                False,
                [],
                call_source_loss=req
                if self.admission == LOSSY and req is not None
                else 0,
                trace={"roster": []},
            )
            if self._dense_global():
                r.cost_distribution = []
                if outputs.min_cut:
                    r.min_cut = []
            return self._finish(r, outputs, None, controller)
        if controller:
            return self._place_controller(req, outputs)
        return self._place_pooled(req, outputs)

    def _finish(
        self, result: PlacementResult, outputs: Outputs, sol, controller=False
    ) -> PlacementResult:
        """The demand's state after the call (contract r8 14.1) from the committed records,
        then the optional outputs."""
        lossy = self.admission == LOSSY
        recs = self.records()
        result.requested = None if self.infinite else self.requested
        result.delivered = sum(r.delivered for r in recs)
        rows = []
        for s in self.sources:
            mine = [r for r in recs if r.sources[0] == s]
            delivered = sum(r.delivered for r in mine)
            if lossy:
                off = sum(r.offered or 0 for r in mine)
                dropped = sum(x for r in mine for x in r.drops.values())
                rows.append(
                    (s, off, off, delivered, off - delivered - dropped, dropped)
                )
            else:
                rows.append((s, None, delivered, delivered, 0, 0))
        result.by_source = rows
        if lossy:
            drops: dict[int, int] = {}
            for r in recs:
                for e, x in r.drops.items():
                    drops[e] = drops.get(e, 0) + x
            result.drops = sorted((e, x) for e, x in drops.items() if x)
            result.offered = self.offered
            result.admitted = self.offered
            result.unadmitted = 0
            result.source_loss = self.offered - result.delivered - sum(drops.values())
        else:
            result.admitted = result.delivered
            result.unadmitted = (
                None if self.infinite else self.requested - result.delivered
            )
        if controller:
            result.total_cost = sum(
                self.g.cost[e] * x for r in recs for e, x in r.edges.items()
            )
        if outputs.edge_flows:
            alloc = [0] * self.g.num_edges
            for r in recs:
                for e, x in r.edges.items():
                    alloc[e] += x
            result.edge_flows = alloc
        if outputs.residual:
            result.residual = self.fg.residual()
        if outputs.bottleneck:
            result.bottleneck = self.fg.bottleneck(self.sources, self.elig)
        if outputs.min_cut and self._dense_global() and sol is not None:
            result.min_cut = sol.min_cut
        return result

    def _dense_global(self) -> bool:
        """Global over every eligible edge (no support): the curve and the cut."""
        return self.mode == "Global" and self.support is None

    # ------------------------------------------------------------ CostOnly (step 2)

    def _prepare(self, W):
        """The epoch's routes of a CostOnly demand on the residual at the first positive call's
        entry: per roster source its route (record bundle, route cost), their union, and one
        continuation per router (contract r8 5.1)."""
        routes = {s: self._cost_only_route(s, W) for s in self.roster}
        union = tuple(sorted({e for r in routes.values() for e in r}))
        for r in routes.values():
            for e in r:
                u = self.g.src[e]
                assert sorted(x for x in r if self.g.src[x] == u) == sorted(
                    x for x in union if self.g.src[x] == u
                ), "two roster sources continue differently at a router"
        return routes, union

    def _source_flow(self, s, route, prev) -> Flow:
        plan = K.Plan(
            self.kernel, self.admission, s, self._dests_ok, route, route, None
        )
        return Flow(
            src=s,
            route=route,
            plan=plan,
            cost=self._route_cost(route, s),
            flow_id=None if prev is None else prev.flow_id,
            offered=0 if self.admission == LOSSY else None,
            new=prev is None,
            seq=self.roster.index(s),
        )

    def _place_pooled(self, req, outputs) -> PlacementResult:
        g, fg = self.g, self.fg
        W = fg.residual()
        created = not self.cohorts
        trace: dict = {"prepare": None}
        if created:
            routes, union = self._prepare(W)
            prog = make_program(
                g,
                self.kernel,
                self.roster,
                self._dests_ok,
                union,
                self.cfg.weights or "Capacity",
                W,
            )
            flows = {s: self._source_flow(s, routes[s], None) for s in self.roster}
            trace["prepare"] = {s: list(routes[s]) for s in self.roster}
        else:
            prog = self.program
            flows = {s: c[0] for s, c in copy.deepcopy(self.cohorts).items()}
        k = len(self.roster)
        if self.admission == LOSSLESS:
            admitted = sum(f.delivered for f in flows.values())
            own = merge(f.alloc for f in flows.values())
            V = add_view(W, own)
            # the largest total in [admitted, admitted + req] whose no-drop forwarding of the
            # round-robin shares fits V and gives no root without a member a share
            bound = 0
            for i, r in enumerate(self.roster):
                if not prog.out.get(r):
                    bound = min(bound, i)
                    break
                bound += sum(V[e] for e in prog.out[r])
            hi = bound if req is None else min(admitted + req, bound)

            def fits(t):
                inj = round_robin(t, k)
                if any(
                    x and not prog.out.get(r)
                    for r, x in zip(self.roster, inj, strict=True)
                ):
                    return False
                _, _, _, carried, _, _ = forward(g, prog, inj, None)
                return all(x <= V[e] for e, x in merge(carried).items())

            lo = admitted
            if hi > lo and fits(hi):
                lo = hi
            else:
                while hi - lo > 1:
                    mid = (lo + hi) // 2
                    if fits(mid):
                        lo = mid
                    else:
                        hi = mid
            _, _, delivered, carried, _, _ = forward(g, prog, round_robin(lo, k), None)
            for i, s in enumerate(self.roster):
                flows[s].alloc = dict(carried[i])
                flows[s].T = flows[s].delivered = delivered[i]
            placed = lo - admitted
            call_drops, source_loss = [], 0
            trace["T"] = [admitted, lo]
        else:
            if req is None:
                inj = None
            else:
                old = round_robin(self.offered, k)
                new = round_robin(self.offered + req, k)
                inj = [b - a for a, b in zip(old, new, strict=True)]
            offered, lost, delivered, carried, dropped, _ = forward(g, prog, inj, W)
            for i, s in enumerate(self.roster):
                f = flows[s]
                for e, x in carried[i].items():
                    f.alloc[e] = f.alloc.get(e, 0) + x
                for e, x in dropped[i].items():
                    f.drops[e] = f.drops.get(e, 0) + x
                f.offered = (f.offered or 0) + offered[i]
                f.delivered += delivered[i]
            placed = sum(delivered)
            call_drops = sorted(merge(dropped).items())
            source_loss = sum(lost)
            self.offered += sum(offered)
            trace["offers"] = {
                s: "inf" if inj is None else inj[i] for i, s in enumerate(self.roster)
            }
        cohorts = {s: [f] for s, f in flows.items()}
        drafts, order = self._drafts(cohorts)
        new_ids = fg._commit(self.owner, drafts)
        self._assign_ids(order, new_ids)
        self.cohorts = cohorts
        self.program = prog
        self.prepared = True
        if req is None:
            self.infinite = True
        else:
            self.requested += req
        if not prog.members:
            term = "no_active_weight"
        elif self.admission == LOSSY:
            term = "complete" if req is not None else "exhausted"
        else:
            term = "complete" if req is not None and placed == req else "exhausted"
        # The flags compare programs (10.4): the epoch's first call creates one (its routes and
        # its weights); an increment changes neither.
        result = PlacementResult(
            placed,
            None if req is None else req - placed,
            term,
            created,
            created,
            call_drops,
            call_source_loss=source_loss,
            trace=trace,
        )
        return self._finish(result, outputs, None)

    # ------------------------------------------------------------ the controllers

    def _place_controller(self, req, outputs) -> PlacementResult:
        """Global over every eligible edge, CostOnly MaxFlow over the union of the roster's
        routes: a minimum-cost flow of the cumulative value from one virtual source (the
        optimized origin)."""
        g, fg = self.g, self.fg
        W = fg.residual()
        old = self.cohorts
        created = not old
        glob = self.mode == "Global"
        dense = self._dense_global()
        trace: dict = {}
        if glob:
            routes, union = None, self.trimmed if not dense else None
        elif created:
            routes, union = self._prepare(W)
            trace["prepare"] = {s: list(routes[s]) for s in self.roster}
        else:
            routes, union = None, self.union
            trace["prepare"] = None
        own = merge(f.alloc for c in old.values() for f in c)
        V = add_view(W, own)
        P = self.placed_demand()
        cap = None if req is None else P + req
        if dense:
            sol = flowcore.solve(g, V, self.roster, self._dests_ok, cap, elig=self.elig)
        else:
            sol = flowcore.solve(g, V, self.roster, self._dests_ok, cap, edges=union)
        assert sol.placed >= P, "a controller call never reduces placed_demand()"
        shares, paths = self._decompose(sol)
        new_cohorts: dict[int, list[Flow]] = {}
        for s in self.roster:
            alloc = shares[s]
            prev = old.get(s, [None])[0]
            if glob:
                route = tuple(sorted(alloc))
                plan = K.Plan(MF, LOSSLESS, s, self._dests_ok, route, route, None)
                f = Flow(
                    src=s,
                    route=route,
                    plan=plan,
                    cost=self._route_cost(route, s),
                    flow_id=None if prev is None else prev.flow_id,
                    new=prev is None,
                    seq=self.roster.index(s),
                )
            else:
                f = (
                    copy.deepcopy(prev)
                    if prev is not None
                    else self._source_flow(s, routes[s], None)
                )
            f.T = f.delivered = sol.supply.get(s, 0)
            f.alloc = alloc
            new_cohorts[s] = [f]
        placed = sol.placed - P
        drafts, order = self._drafts(new_cohorts)
        new_ids = fg._commit(self.owner, drafts)
        self._assign_ids(order, new_ids)
        self.cohorts = new_cohorts
        # The program (10.4): the paths in decomposition order with their amounts, and the
        # optimized origin's supplies (every roster source, zero included). No program and a
        # program without a path are equal.
        program = (
            tuple((s, walk) for s, walk, _ in paths),
            (
                tuple(x for _, _, x in paths),
                tuple(sol.supply.get(s, 0) for s in self.roster),
            ),
        )
        before = self.paths or ((), None)
        if not before[0] and not program[0]:
            routes_changed = weights_changed = False
        else:
            routes_changed = before[0] != program[0]
            weights_changed = before[1] != program[1]
        self.paths = program
        if not glob and created:
            self.union = union
        self.prepared = True
        if req is None:
            self.infinite = True
        else:
            self.requested += req
        term = "complete" if req is not None and placed == req else "exhausted"
        if dense:
            trace.update(
                {
                    "V": V,
                    "cap": cap,
                    "flow_core": {**sol.payload(), "supply": sol.supply},
                    "shares": {s: sorted(a.items()) for s, a in shares.items()},
                }
            )
            self.curve = list(sol.cost_distribution)
        else:
            trace.update(
                {
                    "cap": cap,
                    "flow_core": {
                        "placed": sol.placed,
                        "edge_flows": list(sol.edge_flows),
                        "supply": sol.supply,
                    },
                    "shares": {s: sorted(a.items()) for s, a in shares.items()},
                }
            )
        result = PlacementResult(
            placed,
            None if req is None else req - placed,
            term,
            routes_changed,
            weights_changed,
            [],
            cost_distribution=list(sol.cost_distribution) if dense else None,
            trace=trace,
        )
        return self._finish(result, outputs, sol, True)

    def _decompose(self, sol):
        """Canonical per-source decomposition (9.1, contract r8 4.1): sources ascending; the
        smallest EdgeId walk with positive remaining flow to the first destination reached,
        saturated by min(walk bottleneck, remaining virtual supply). Returns the per-source shares
        and the paths (source, walk, amount) in that order."""
        g = self.g
        remaining = dict(enumerate(sol.edge_flows))
        dests = set(self._dests_ok)
        out: dict[int, dict[int, int]] = {}
        paths: list[tuple[int, tuple[int, ...], int]] = []
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
                paths.append((s, tuple(walk), amount))
            out[s] = mine
        assert all(x == 0 for x in remaining.values()), "decomposition covers the flow"
        return out, paths

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
