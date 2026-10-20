"""The worked examples, run on the library: each example's own input and expected
values (copied as written), one adapter per kind of example that builds the input with the library
and returns what the example names, compared field by field (an expected mapping may name fewer
fields than the adapter returns; an int is not a float or a bool; a Fraction compares by its exact
string). The examples of several demands at once run as one apply_actions list (each route or head
its own demand); the examples of concurrent publication are not here.

Examples given as routes of link indices (no nodes) are laid out as one head, a feeder link to each
route's first link, a connector between consecutive links and an exit link to one destination;
these added links are wider than anything offered, so they never bind and the example's links
decide every split. An example's per-route offers are the routes' weights, so the head's split by
highest averages hands each route exactly its offer."""

from __future__ import annotations

import copy
import itertools
import math
from decimal import Decimal
from fractions import Fraction

import _netgraph_core as ext
import numpy as np
import pytest

import netgraph_core as ngc

INF = math.inf
LL = ngc.Admission.LOSSLESS
LY = ngc.Admission.LOSSY
S = ngc.SourceAllocation
WIDE = 1 << 50  # wider than any amount an example offers


class NotExposed(Exception):
    """The example names a value the library does not report (an internal quantity)."""


# -- comparison -------------------------------------------------------------------------------


def mismatch(got, want, where="value"):
    """None if `got` matches `want`, else the first difference. Typed and recursive: an int must be
    an int (not a float or a bool), a bool a bool, None None; an expected mapping names the keys it
    checks; a sequence must have the expected length; a Fraction compares by its exact string."""
    if isinstance(got, Fraction):
        got = f"{got.numerator}/{got.denominator}" if got.denominator != 1 else str(got)
    if isinstance(want, dict):
        if not isinstance(got, dict):
            return f"{where}: got {got!r}, expected a mapping"
        for key, value in want.items():
            if key not in got:
                return f"{where}.{key}: missing"
            sub = mismatch(got[key], value, f"{where}.{key}")
            if sub:
                return sub
        return None
    if isinstance(want, (list, tuple)):
        if not isinstance(got, (list, tuple)) or len(got) != len(want):
            return f"{where}: got {got!r}, expected {want!r}"
        for n, (g, w) in enumerate(zip(got, want, strict=True)):
            sub = mismatch(g, w, f"{where}[{n}]")
            if sub:
                return sub
        return None
    if type(got) is not type(want) or got != want:
        return f"{where}: got {got!r}, expected {want!r}"
    return None


# -- building blocks --------------------------------------------------------------------------


def graph_of(edges):
    """The example's graph: [source, destination, capacity, cost] rows in EdgeId order."""
    n = 1 + max(max(a, b) for a, b, _, _ in edges)
    g = ngc.Graph(
        n,
        [e[0] for e in edges],
        [e[1] for e in edges],
        [e[2] for e in edges],
        [e[3] for e in edges],
        1.0,
    )
    assert g.input_positions().tolist() == list(range(len(edges)))
    return g


def dag_of(g, edges):
    """The PredecessorDAG of the given links (each link's tail is its predecessor)."""
    n = g.num_nodes()
    tails = g.edge_sources().tolist()
    heads = g.edge_destinations().tolist()
    preds = [[] for _ in range(n)]
    for e in sorted(edges):
        preds[heads[e]].append((tails[e], e))
    offsets, parents, via = [0], [], []
    for v in range(n):
        for u, e in preds[v]:
            parents.append(u)
            via.append(e)
        offsets.append(len(parents))
    return ngc.PredecessorDAG.from_csr(g, offsets, parents, via)


def path_of(g, edges):
    return ngc.Path(g, [[e] for e in edges])


def per_edge(pairs, m):
    out = [0] * m
    for e, x in pairs:
        out[e] += x
    return out


def config(routing, admission=LL):
    return ngc.FlowPolicyConfig(routing=routing, admission=admission)


def policy(g, sources, destinations, cfg, mask=None, ledger=None):
    return ngc.FlowPolicy(
        ledger if ledger is not None else ngc.FlowLedger(g),
        sources,
        destinations,
        cfg,
        mask,
    )


def pinned(g, program, admission=LL, mask=None, black_holes=None, ledger=None):
    cfg = config(ngc.Pinned(program, black_holes=black_holes), admission)
    return policy(
        g, list(program.sources), list(program.destinations), cfg, mask, ledger
    )


def edge_mask(g, links=(), nodes=()):
    """The Eligibility excluding the given links and nodes (None: nothing excluded)."""
    if not links and not nodes:
        return None
    return g.eligibility(excluded_nodes=list(nodes), excluded_edges=list(links))


def laid_out(cap, routes):
    """The layout of routes given as link indices: one head (node 0), one destination (node 1),
    per link its own two nodes, a feeder from the head to each first link, a connector between
    consecutive links and an exit from each last link. Returns the graph and per example link its
    EdgeId, and per route its path of EdgeIds."""
    rows = []
    node = itertools.count(2)
    ends = {e: (next(node), next(node)) for e in range(len(cap))}
    for e, c in enumerate(cap):
        rows.append((ends[e][0], ends[e][1], c, 1, ("link", e)))
    added = {}

    def link(a, b, key):
        if key not in added:
            added[key] = len(rows)
            rows.append((a, b, WIDE, 1, key))
        return key

    walks = []
    for route in routes:
        walk = [link(0, ends[route[0]][0], ("feed", route[0]))]
        for e, f in itertools.pairwise(route):
            walk += [("link", e), link(ends[e][1], ends[f][0], ("join", e, f))]
        walk += [("link", route[-1]), link(ends[route[-1]][1], 1, ("exit", route[-1]))]
        walks.append(walk)
    n = 2 + 2 * len(cap)
    g = ngc.Graph(
        n,
        [r[0] for r in rows],
        [r[1] for r in rows],
        [r[2] for r in rows],
        [r[3] for r in rows],
        1.0,
    )
    position = {r[4]: i for i, r in enumerate(rows)}
    eid = {pos: e for e, pos in enumerate(g.input_positions().tolist())}
    edge = {key: eid[pos] for key, pos in position.items()}
    links = [edge[("link", e)] for e in range(len(cap))]
    paths = [[edge[k] for k in walk] for walk in walks]
    return g, links, paths


def laid_out_tunnels(cap, routes, weights, groups=None):
    """The program of tunnels along laid-out routes: one head per group (default: one head for
    every route), each route a tunnel of its weight. Returns the graph, the example links'
    EdgeIds and the program."""
    if groups is None:
        g, links, paths = laid_out(cap, routes)
        tunnels = [
            ngc.RouteDefinition(0, path_of(g, p), weight=w)
            for p, w in zip(paths, weights, strict=True)
        ]
        return g, links, ngc.ForwardingProgram.from_routes(g, [0], [1], tunnels)
    # One head per group: a node per group feeding the shared head-side layout.
    g, links, paths = laid_out_groups(cap, routes, groups)
    tunnels = [
        ngc.RouteDefinition(head, path_of(g, p), weight=w)
        for (head, p), w in zip(paths, weights, strict=True)
    ]
    heads = sorted({h for h, _ in paths})
    return g, links, ngc.ForwardingProgram.from_routes(g, heads, [1], tunnels)


def laid_out_groups(cap, routes, groups):
    """As laid_out, with route i fed from its group's own head (heads are nodes 2 + 2L ...)."""
    rows = []
    ends = {e: (2 + 2 * e, 3 + 2 * e) for e in range(len(cap))}
    for e, c in enumerate(cap):
        rows.append((ends[e][0], ends[e][1], c, 1, ("link", e)))
    first_head = 2 + 2 * len(cap)
    head_of = {grp: first_head + i for i, grp in enumerate(sorted(set(groups)))}
    added = set()

    def link(a, b, key):
        if key not in added:
            added.add(key)
            rows.append((a, b, WIDE, 1, key))
        return key

    walks = []
    for route, grp in zip(routes, groups, strict=True):
        h = head_of[grp]
        walk = [link(h, ends[route[0]][0], ("feed", grp, route[0]))]
        for e, f in itertools.pairwise(route):
            walk += [("link", e), link(ends[e][1], ends[f][0], ("join", e, f))]
        walk += [("link", route[-1]), link(ends[route[-1]][1], 1, ("exit", route[-1]))]
        walks.append((h, walk))
    n = first_head + len(head_of)
    g = ngc.Graph(
        n,
        [r[0] for r in rows],
        [r[1] for r in rows],
        [r[2] for r in rows],
        [r[3] for r in rows],
        1.0,
    )
    position = {r[4]: i for i, r in enumerate(rows)}
    eid = {pos: e for e, pos in enumerate(g.input_positions().tolist())}
    edge = {key: eid[pos] for key, pos in position.items()}
    links = [edge[("link", e)] for e in range(len(cap))]
    return g, links, [(h, [edge[k] for k in walk]) for h, walk in walks]


def factor_of(spec):
    """An example's CostStretch factor: itself, or ["float", x], ["fraction", n, d], ["int", n],
    ["bool", b], ["index", n] (an object known only by __index__), ["numpy", type, value],
    ["str", s] or ["other", x] (a value of another type)."""
    if not isinstance(spec, list):
        return spec
    kind = spec[0]
    if kind == "float":
        return float(spec[1])
    if kind == "fraction":
        return Fraction(spec[1], spec[2])
    if kind == "index":
        return _Index(spec[1])
    if kind == "numpy":
        return getattr(np, spec[1])(spec[2])
    if kind in ("str", "other"):
        return spec[1]
    return bool(spec[1]) if kind == "bool" else int(spec[1])


class _Index:
    def __init__(self, n):
        self.n = n

    def __index__(self):
        return self.n


def ceiling_of(spec):
    """A ceiling argument: a CostStretch of the factor, or the raw value of another type."""
    if isinstance(spec, list) and spec[0] == "other":
        return spec[1]
    return ngc.CostStretch(factor_of(spec))


def error_text(error):
    return f"TypeError: {error}" if isinstance(error, TypeError) else str(error)


# -- splitting rules ----------------------------------------------------------------------------


def split_on_parallel_links(n, weights):
    """H(n, w) as the library splits it: one head, a tunnel of weight w_i on its own link."""
    if n >= 1 << 62:
        raise NotExposed("an amount at or above 2^62 is outside the library's domain")
    m = len(weights)
    g = ngc.Graph(2, [0] * m, [1] * m, [max(n, 1)] * m, [1] * m, 1.0)
    tunnels = [
        ngc.RouteDefinition(0, path_of(g, [e]), weight=w) for e, w in enumerate(weights)
    ]
    program = ngc.ForwardingProgram.from_routes(g, [0], [1], tunnels)
    return [t.assigned for t in pinned(g, program, LY).place(n).by_tunnel]


def run_ha(n, weights):
    return split_on_parallel_links(n, weights)


def split_report(n, weights, assigned):
    total = sum(weights)
    exact = [Fraction(n * w, total) for w in weights]
    return dict(
        assigned=assigned,
        max_excess=max(Fraction(a) - x for a, x in zip(assigned, exact, strict=True)),
        max_deficit=max(x - Fraction(a) for a, x in zip(assigned, exact, strict=True)),
        below_k=all(
            abs(a - x) < len(weights) for a, x in zip(assigned, exact, strict=True)
        ),
        below_one=all(abs(a - x) < 1 for a, x in zip(assigned, exact, strict=True)),
    )


def run_split_error(n, weights):
    return split_report(n, weights, split_on_parallel_links(n, weights))


def run_bound_sweep(max_total, max_members, max_weight):
    # One program per weight vector; every total is an evaluation of it (nothing published).
    cases = violations = equal_violations = 0
    for k in range(1, max_members + 1):
        g = ngc.Graph(2, [0] * k, [1] * k, [max(max_total, 1)] * k, [1] * k, 1.0)
        for weights in itertools.product(range(max_weight + 1), repeat=k):
            if not any(weights):
                continue
            tunnels = [
                ngc.RouteDefinition(0, path_of(g, [e]), weight=w)
                for e, w in enumerate(weights)
            ]
            p = pinned(g, ngc.ForwardingProgram.from_routes(g, [0], [1], tunnels), LY)
            for n in range(max_total + 1):
                cases += 1
                assigned = [t.assigned for t in p.evaluate(n).by_tunnel]
                r = split_report(n, list(weights), assigned)
                violations += not r["below_k"]
                if len(set(weights)) == 1:
                    equal_violations += not r["below_one"]
    return dict(cases=cases, violations=violations, equal_violations=equal_violations)


def run_member_order(amount, order):
    # The members of one split, taken in prepared (EdgeId) order: the member order[j] is link j.
    k = len(order)
    g = ngc.Graph(2, [0] * k, [1] * k, [max(amount, 1)] * k, [1] * k, 1.0)
    program = ngc.ForwardingProgram.from_dag(g, [0], [1], dag_of(g, range(k)))
    flows = (
        pinned(g, program, LY)
        .place(amount, ngc.OutputOptions(edge_flows=True))
        .edge_flows.tolist()
    )
    result = [0] * k
    for j, member in enumerate(order):
        result[member] = flows[j]
    return result


def run_pooled(sources, members):
    # Sources 0..s-1 each feed router s by a wide link; the router reaches s + 1 over `members`
    # wide parallel links. Pooled: from_dag; independent: one split route per source.
    s = len(sources)
    if len(set(sources)) != 1:
        raise NotExposed("unequal source amounts: a program's sources are equal")
    rows = [(i, s, WIDE) for i in range(s)] + [(s, s + 1, WIDE)] * members
    g = ngc.Graph(
        s + 2,
        [r[0] for r in rows],
        [r[1] for r in rows],
        [r[2] for r in rows],
        [1] * len(rows),
        1.0,
    )
    member_links = list(range(s, s + members))
    whole = dag_of(g, range(len(rows)))
    pooled_policy = pinned(
        g, ngc.ForwardingProgram.from_dag(g, list(range(s)), [s + 1], whole), LY
    )
    flows = pooled_policy.place(
        sum(sources), ngc.OutputOptions(edge_flows=True)
    ).edge_flows.tolist()
    records = sorted(pooled_policy.records(), key=lambda r: r.sources[0])
    by_column = [
        [dict(rec.edge_flows).get(e, 0) for rec in records] for e in member_links
    ]
    separate = ngc.ForwardingProgram.from_routes(
        g, list(range(s)), [s + 1], [ngc.RouteDefinition(i, whole) for i in range(s)]
    )
    independent = (
        pinned(g, separate, LY)
        .place(sum(sources), ngc.OutputOptions(edge_flows=True))
        .edge_flows.tolist()
    )
    return dict(
        columns=[flows[e] for e in member_links],
        source_by_column=by_column,
        independent_columns=[independent[e] for e in member_links],
    )


# -- the ledger and the policy's operations ------------------------------------------------------


def run_ledger(service, reservation, carried, strict, new_delivery, old_delivery):
    # One link: a ShortestPaths policy carries `old_delivery` (then is replaced by `new_delivery`)
    # and a DistributedTE policy holds `strict` for its design.
    g = ngc.Graph(2, [0], [1], [max(service, reservation)], [1], 1.0)
    ledger = ngc.FlowLedger(g, [service], reservation_baseline=[reservation])
    carrier = policy(g, 0, 1, config(ngc.ShortestPaths(), LY), ledger=ledger)
    assert carrier.place(old_delivery).call_delivered == carried
    holder = policy(
        g, 0, 1, config(ngc.DistributedTE(ngc.TunnelCount(1)), LY), ledger=ledger
    )
    holder.replace(0, design_demand=strict)
    out = dict(
        carried=int(ledger.carried()[0]),
        residual=int(ledger.residual()[0]),
        reserved=int(ledger.reserved()[0]),
        unreserved=int(ledger.unreserved()[0]),
    )
    out["delivery_delta"] = carrier.replace(new_delivery).delivery_delta
    return out


def run_refusal(capacity, outside, calls):
    # One link; another owner holds outside[c] at the entry of call c; a Lossless policy places
    # each call; its delivery after each call. The rejected rule (re-offering refused units) is
    # what the library must not do: only the correct sequence is reported.
    g = ngc.Graph(2, [0], [1], [capacity], [1], 1.0)
    ledger = ngc.FlowLedger(g)
    other = policy(g, 0, 1, config(ngc.ShortestPaths()), ledger=ledger)
    mine = policy(g, 0, 1, config(ngc.ShortestPaths()), ledger=ledger)
    out = []
    for load, offer in zip(outside, calls, strict=True):
        other.replace(load)
        out.append(mine.place(offer).delivered)
    return dict(correct=out)


def private_links(caps, costs=None):
    """Sources 0..n-1, each with its own link to node n."""
    n = len(caps)
    return ngc.Graph(
        n + 1, list(range(n)), [n] * n, list(caps), list(costs or [1] * n), 1.0
    )


def run_increments(caps, calls, admission, outside=None):
    n = len(caps)
    g = private_links(caps)
    ledger = ngc.FlowLedger(g)
    adm = getattr(ngc.Admission, admission)
    others = [
        policy(g, i, n, config(ngc.ShortestPaths()), ledger=ledger) for i in range(n)
    ]
    mine = policy(g, list(range(n)), n, config(ngc.ShortestPaths(), adm), ledger=ledger)
    delivered, r = [], None
    for c, v in enumerate(calls):
        if outside is not None:
            for i, load in enumerate(outside[c]):
                others[i].replace(load)
        r = mine.place(v)
        delivered.append(r.call_delivered)
    out = dict(
        call_delivered=delivered,
        requested=r.requested,
        delivered=r.delivered,
        unadmitted=r.unadmitted,
        offered=r.offered,
    )
    if admission == "LOSSLESS":
        out["admitted"] = r.admitted
    else:
        out["edge_drops"] = per_edge(r.edge_drops, n)
    return out


def run_epoch(cap, calls, full=False, admission="LOSSY"):
    # ShortestPaths ECMP over parallel links.
    k = len(cap)
    g = ngc.Graph(2, [0] * k, [1] * k, list(cap), [1] * k, 1.0)
    ledger = ngc.FlowLedger(g)
    p = policy(
        g,
        0,
        1,
        config(ngc.ShortestPaths(), getattr(ngc.Admission, admission)),
        ledger=ledger,
    )
    per_call, r = [], None
    for v in [sum(calls)] if full else calls:
        r = p.place(v)
        per_call.append(r.call_delivered)
    return dict(
        requested=sum(calls),
        delivered=r.delivered,
        edge_drops=sum(x for _, x in r.edge_drops),
        edge_flows=ledger.carried().tolist(),
        residual=ledger.residual().tolist(),
        per_call_delivery=per_call,
    )


def run_equal_increment(capacities, costs, calls):
    # The equal-source controller: CentralizedTE EQUAL, sources with private links.
    k = len(capacities)
    g = private_links(capacities, costs)
    ledger = ngc.FlowLedger(g)
    p = policy(
        g,
        list(range(k)),
        k,
        config(ngc.CentralizedTE(source_allocation=S.EQUAL)),
        ledger=ledger,
    )
    deliveries, r = [], None
    for v in calls:
        r = p.place(v)
        deliveries.append(r.call_delivered)
    return dict(
        call_deliveries=deliveries,
        allocation=ledger.carried().tolist(),
        delivered=r.delivered,
        cost=r.total_cost,
    )


def run_zero_state(service, outside):
    # place(0) on a ShortestPaths policy without a prepared program, new and after remove(); another
    # owner's tunnels hold `outside` (weights = its loads) on parallel links of capacities `service`.
    k = len(service)
    g = ngc.Graph(2, [0] * k, [1] * k, list(service), [1] * k, 1.0)
    ledger = ngc.FlowLedger(g)
    held = [
        ngc.RouteDefinition(0, path_of(g, [e]), weight=w) for e, w in enumerate(outside)
    ]
    pinned(
        g, ngc.ForwardingProgram.from_routes(g, [0], [1], held), ledger=ledger
    ).place(sum(outside))
    outputs = ngc.OutputOptions(
        edge_flows=True, residual=True, edge_reservations=True, min_cut=True
    )

    def state(r):
        return dict(
            call_delivered=r.call_delivered,
            call_undelivered=r.call_undelivered,
            delivery_delta=r.delivery_delta,
            call_edge_drops=r.call_edge_drops,
            call_source_drops=r.call_source_drops,
            termination=r.termination.name,
            routes_changed=r.routes_changed,
            weights_changed=r.weights_changed,
            reservations_changed=r.reservations_changed,
            program=r.program,
            edge_flows=r.edge_flows.tolist(),
            edge_reservations=r.edge_reservations.tolist(),
            residual=r.residual.tolist(),
            min_cut=r.min_cut,
            marginal_cost_curve=r.marginal_cost_curve,
        )

    fresh = state(
        policy(g, 0, 1, config(ngc.ShortestPaths()), ledger=ledger).place(0, outputs)
    )
    used = policy(g, 0, 1, config(ngc.ShortestPaths()), ledger=ledger)
    used.place(1)
    used.remove()
    again = state(used.place(0, outputs))
    assert again == fresh, (fresh, again)
    return fresh


def run_change_flags(path_cap, calls, recomputed):
    # One link; a recomputed program (CentralizedTE) or a fixed one (Pinned over one tunnel).
    g = ngc.Graph(2, [0], [1], [path_cap], [1], 1.0)
    if recomputed:
        p = policy(g, 0, 1, config(ngc.CentralizedTE(source_allocation=S.OPTIMIZED)))
    else:
        tunnel = [ngc.RouteDefinition(0, path_of(g, [0]))]
        p = pinned(g, ngc.ForwardingProgram.from_routes(g, [0], [1], tunnel))
    flags = []
    for i, (op, v) in enumerate(calls):
        r = getattr(p, op)(v)
        if i:
            flags.append(r)
    return dict(
        routes_changed=[r.routes_changed for r in flags],
        weights_changed=[r.weights_changed for r in flags],
        reservations_changed=[r.reservations_changed for r in flags],
    )


# -- programs: tunnels, split routes, the LOSSY rule ------------------------------------------------


def run_lossy(offers, routes, cap, labels=None):
    """Routes of link indices carrying their offers. Each demand of the labels is one source of
    the program (a demand's routes are its tunnels, weights = offers), so a link splits by highest
    averages over sources, then tunnels: the labels' demand-then-tunnel order. Several demands must
    offer equal totals (a program's sources are equal)."""
    labels = labels or [[i, 0, 0] for i in range(len(offers))]
    demands = [lab[0] for lab in labels]
    if len(set(demands)) == len(demands):
        groups = (
            None  # one tunnel per demand: one head, split by H over the routes' offers
        )
        total = sum(offers)
    else:
        totals = {}
        for d, o in zip(demands, offers, strict=True):
            totals[d] = totals.get(d, 0) + o
        if len(set(totals.values())) != 1:
            raise NotExposed("demands with unequal totals")
        groups = demands
        total = sum(offers)
    g, links, program = laid_out_tunnels(cap, routes, offers, groups)
    p = pinned(g, program, LY)
    try:
        r = p.place(total, ngc.OutputOptions(edge_flows=True))
    except ValueError as e:
        return dict(error=str(e))
    flows = r.edge_flows.tolist()
    drops = per_edge(r.call_edge_drops, g.num_edges())
    carried, edge_drops = [], []
    records = {rec.entry_id: rec for rec in p.records()}
    for i, route in enumerate(routes):
        rec_flows = dict(records[i].edge_flows) if i in records else {}
        rec_drops = dict(records[i].edge_drops) if i in records else {}
        carried.append([rec_flows.get(links[e], 0) for e in route])
        edge_drops.append([rec_drops.get(links[e], 0) for e in route])
    return dict(
        carried=carried,
        edge_drops=edge_drops,
        delivered=[t.delivered for t in r.by_tunnel],
        edge_flows=[flows[links[e]] for e in range(len(cap))],
        link_drops=[drops[links[e]] for e in range(len(cap))],
    )


def run_bound(cap, routes, weights):
    # One head, a tunnel per route of the given weight: the Lossless bound of the program.
    g, _, program = laid_out_tunnels(cap, routes, weights)
    return pinned(g, program).place(INF).call_delivered


def run_split_contexts(cap, routes, weights, pool):
    # One head, split routes whose members are parallel links to the destination. Separate: one
    # split route per entry (from_routes); pooled: the one shared context of from_dag.
    k = len(cap)
    g = ngc.Graph(2, [0] * k, [1] * k, list(cap), [1] * k, 1.0)
    if pool:
        if len({tuple(r) for r in routes}) != 1:
            raise NotExposed("pooled routes with different members")
        program = ngc.ForwardingProgram.from_dag(g, [0], [1], dag_of(g, routes[0]))
    else:
        program = ngc.ForwardingProgram.from_routes(
            g,
            [0],
            [1],
            [
                ngc.RouteDefinition(0, dag_of(g, r), weight=w)
                for r, w in zip(routes, weights, strict=True)
            ],
        )
    bound = pinned(g, program).place(INF).call_delivered
    r = pinned(g, program, LY).place(2, ngc.OutputOptions(edge_flows=True))
    drops = per_edge(r.call_edge_drops, k)
    return dict(
        bound=bound,
        loads_at_2=[f + d for f, d in zip(r.edge_flows.tolist(), drops, strict=True)],
    )


def run_dag_contexts(edges, sources, destination, offer, pooled):
    g = graph_of(edges)
    whole = dag_of(g, range(len(edges)))
    if pooled:
        program = ngc.ForwardingProgram.from_dag(g, sources, [destination], whole)
    else:
        program = ngc.ForwardingProgram.from_routes(
            g, sources, [destination], [ngc.RouteDefinition(s, whole) for s in sources]
        )
    return dict(
        lossless_bound=pinned(g, program).place(INF).call_delivered,
        lossy_delivered=pinned(g, program, LY).place(offer).call_delivered,
    )


def run_stale_dag(
    edges, sources, destinations, offer, weights,
    black_links=(), black_nodes=(), mask_links=(), mask_nodes=(),
):  # fmt: skip
    # One split route (from_dag of the member links, WEIGHTED by the given frozen weights) forwarded
    # by Pinned under the repair mask, with the black holes laid over.
    g = graph_of(edges)
    members = [e for e, w in enumerate(weights) if w is not None]
    program = ngc.ForwardingProgram.from_dag(
        g,
        sources,
        [destinations],
        dag_of(g, members),
        allocation=ngc.RouteAllocation.WEIGHTED,
        weight_basis=ngc.WeightBasis.EXPLICIT,
        edge_weights=[w or 0 for w in weights],
    )
    p = pinned(
        g,
        program,
        LY,
        mask=edge_mask(g, mask_links, mask_nodes),
        black_holes=edge_mask(g, black_links, black_nodes),
    )
    r = p.place(offer, ngc.OutputOptions(edge_flows=True))
    return dict(
        delivered=r.call_delivered,
        edge_flows=r.edge_flows.tolist(),
        edge_drops=per_edge(r.call_edge_drops, len(edges)),
        source_drops=r.call_source_drops,
    )


def run_stale_paths(
    edges, tunnels, offers, black_links=(), black_nodes=(), mask_links=(), mask_nodes=()
):
    # Tunnels from their heads (one program; each head its own source, offers equal), forwarded by
    # Pinned under the repair mask with the black holes laid over.
    g = graph_of(edges)
    heads = sorted({h for h, _ in tunnels})
    if len({offers[h] for h in heads}) != 1:
        raise NotExposed("heads with unequal offers")
    program = ngc.ForwardingProgram.from_routes(
        g,
        heads,
        sorted({edges[p[-1]][1] for _, p in tunnels}),
        [ngc.RouteDefinition(h, path_of(g, p)) for h, p in tunnels],
    )
    mask = edge_mask(g, mask_links, mask_nodes)
    holes = edge_mask(g, black_links, black_nodes)
    total = sum(offers[h] for h in heads)
    try:
        r = pinned(g, program, LY, mask=mask, black_holes=holes).place(total)
    except ValueError as e:
        return dict(error=str(e))
    out = dict(
        status=[t.status.name for t in r.by_tunnel],
        assigned=[t.assigned for t in r.by_tunnel],
        delivered=[t.delivered for t in r.by_tunnel],
        edge_drops=per_edge(r.call_edge_drops, len(edges)),
        source_drops=r.call_source_drops,
    )
    if len(heads) == 1 and not r.call_source_drops:
        out["lossless_admitted"] = (
            pinned(g, program, mask=mask, black_holes=holes).place(total).call_delivered
        )
    return out


def run_termination(offer, model=None, args=None, program=None):
    # The installed program decides: a T1 or T2 setup, or a program of tunnels (head, path,
    # weight) on one link per head.
    if model is not None:
        g = graph_of(args["edges"])
        count = ngc.TunnelCount(args["per_head"])
        if model == "t1":
            routing = ngc.DistributedTE(count)
        else:
            routing = ngc.CentralizedTE(source_allocation=S.EQUAL, tunnel_count=count)
        p = policy(g, args["heads"], args["destinations"], config(routing, LY))
        r = p.replace(offer, design_demand=args.get("design_demand"))
        installed = p.program
        forward_error = None
        try:
            p.forward(offer)
        except ValueError as e:
            forward_error = str(e)
        lossless = pinned(g, installed).place(offer).termination.name
    else:
        roster, entries = program
        n = 1 + max(roster)
        g = ngc.Graph(n + 1, list(range(n)), [n] * n, [WIDE] * n, [1] * n, 1.0)
        routes = []
        for head, has_path, weight in entries:
            if not has_path:
                raise NotExposed("an entry without a path")
            routes.append(ngc.RouteDefinition(head, path_of(g, [head]), weight=weight))
        installed = ngc.ForwardingProgram.from_routes(g, roster, [n], routes)
        p = pinned(g, installed, LY)
        r = p.place(offer)
        forward_error = None
        lossless = pinned(g, installed).place(offer).termination.name
    return dict(
        tunnels=len(installed.tunnels),
        lossy=r.termination.name,
        lossless=lossless,
        source_drops=[row.source_drops for row in r.by_source if row.source_drops],
        forward_error=forward_error,
    )


# -- the tunnel models ------------------------------------------------------------------------------


def t1_routing(per_head, absolute=None, stretch=None, minimum=0):
    ceiling = None
    if absolute is not None:
        ceiling = ngc.AbsoluteCost(absolute)
    elif stretch is not None:
        ceiling = ceiling_of(stretch)
    return ngc.DistributedTE(
        ngc.TunnelCount(per_head),
        cost_ceiling=ceiling,
        min_unreserved_bandwidth=minimum,
    )


def run_t1(
    edges, heads, destinations, design_demand, per_head, offer=None,
    excluded=(), absolute=None, stretch=None, minimum=0,
):  # fmt: skip
    g = graph_of(edges)
    mask = edge_mask(g, excluded)
    try:
        routing = t1_routing(per_head, absolute, stretch, minimum)
        ledger = ngc.FlowLedger(g)
        setup = policy(g, heads, destinations, config(routing), mask, ledger)
    except (TypeError, ValueError) as e:
        return dict(error=error_text(e))
    r = setup.replace(INF, design_demand=design_demand)
    tunnels = r.program.tunnels
    out = dict(
        requested_bandwidth=[t.requested_bandwidth for t in tunnels],
        routes=[list(t.edges) if t.edges is not None else [] for t in tunnels],
        status=[t.status.name for t in tunnels],
        up=sum(t.status == ngc.TunnelStatus.UP for t in tunnels),
        total_reserved_bandwidth=r.plan.total_reserved_bandwidth,
        edge_reservations=ledger.reserved().tolist(),
        bound=r.call_delivered,
    )
    if offer is not None:
        lossy = policy(g, heads, destinations, config(routing, LY), mask).replace(
            offer, design_demand=design_demand
        )
        assigned = [t.assigned for t in lossy.by_tunnel]
        load = [0] * len(edges)
        for t, a in zip(lossy.program.tunnels, assigned, strict=True):
            for e in t.edges or ():
                load[e] += a
        lossless = policy(g, heads, destinations, config(routing), mask).replace(
            offer, design_demand=design_demand
        )
        out.update(
            offered_by_entry=assigned,
            source_drops=lossy.call_source_drops,
            no_drop_load=load,
            lossy_delivery=lossy.call_delivered,
            lossy_edge_drops=per_edge(lossy.call_edge_drops, len(edges)),
            lossless_admitted=lossless.call_delivered,
            lossless_unadmitted=lossless.call_undelivered,
        )
    return out


def run_head_requests(heads, per_head, design_demand):
    # T1's requests in program order, heads with one wide link each.
    g = ngc.Graph(
        heads + 1, list(range(heads)), [heads] * heads, [WIDE] * heads, [1] * heads, 1.0
    )
    p = policy(
        g,
        list(range(heads)),
        heads,
        config(ngc.DistributedTE(ngc.TunnelCount(per_head)), LY),
    )
    r = p.replace(0, design_demand=design_demand)
    return [t.requested_bandwidth for t in r.program.tunnels]


def t2_policy(
    edges,
    heads,
    destinations,
    count,
    prefer_headroom,
    admission,
    cap=None,
    support=None,
    excluded=(),
):
    g = graph_of(edges)
    routing = ngc.CentralizedTE(
        source_allocation=S.EQUAL,
        support=None if support is None else dag_of(g, support),
        tunnel_count=count,
        prefer_headroom=prefer_headroom,
    )
    ledger = ngc.FlowLedger(g, cap)
    return policy(
        g,
        heads,
        destinations,
        config(routing, admission),
        edge_mask(g, excluded),
        ledger,
    )


def run_t2(
    edges, heads, destinations, per_head, design_demand=None, cap=None,
    prefer_headroom=False, lossy_offer=None, support=None, excluded=(),
):  # fmt: skip
    count = ngc.TunnelCount(per_head)
    p = t2_policy(
        edges, heads, destinations, count, prefer_headroom, LL, cap, support, excluded
    )
    r = p.replace(INF, design_demand=design_demand)
    plan = r.plan
    tunnels = r.program.tunnels
    counts = [0] * len(edges)
    for t in tunnels:
        for e in t.edges:
            counts[e] += 1
    out = dict(
        heads_with_capacity=list(plan.heads_with_capacity),
        max_fractional_bandwidth_per_tunnel=plan.max_fractional_bandwidth_per_tunnel,
        max_bandwidth_per_tunnel=plan.max_bandwidth_per_tunnel,
        planned_bandwidth_per_tunnel=plan.planned_bandwidth_per_tunnel,
        total_planned_bandwidth=plan.total_planned_bandwidth,
        total_tunnel_cost=plan.total_tunnel_cost,
        headroom_penalty=plan.headroom_penalty,
        max_lossless_offer=plan.max_lossless_offer,
        insufficient_resolution=plan.insufficient_resolution,
        design_by_head=[list(x) for x in plan.design_by_head],
        tunnels_per_head=plan.tunnels_per_head,
        counts=counts,
        paths=[list(t.edges) for t in tunnels],
        entry_ids=[t.entry_id for t in tunnels],
    )
    if lossy_offer is not None:
        q = t2_policy(
            edges,
            heads,
            destinations,
            count,
            prefer_headroom,
            LY,
            cap,
            support,
            excluded,
        )
        lr = q.replace(lossy_offer, design_demand=design_demand)
        out.update(
            lossy_delivery=lr.call_delivered, lossy_source_drops=lr.call_source_drops
        )
    return out


def run_upto(
    edges, heads, destinations, limit, design_demand=None, prefer_headroom=False
):
    count = ngc.TunnelCount(limit, up_to=True)
    p = t2_policy(edges, heads, destinations, count, prefer_headroom, LL)
    offer = INF if design_demand is None else design_demand
    r = p.replace(offer, design_demand=design_demand)
    return dict(
        tunnels_per_head=r.plan.tunnels_per_head,
        delivered=r.call_delivered,
        paths=[list(t.edges) for t in r.program.tunnels],
    )


# -- the controllers ------------------------------------------------------------------------------


def run_flow(edges, heads, destinations, limit, source_allocation="EQUAL"):
    g = graph_of(edges)
    routing = ngc.CentralizedTE(source_allocation=getattr(S, source_allocation))
    r = policy(g, heads, destinations, config(routing)).place(
        limit, ngc.OutputOptions(edge_flows=True)
    )
    return dict(
        delivered=r.call_delivered, cost=r.total_cost, edge_flows=r.edge_flows.tolist()
    )


def run_min_cut(edges, heads, destinations):
    g = graph_of(edges)
    routing = ngc.CentralizedTE(source_allocation=S.OPTIMIZED)
    r = policy(g, heads, destinations, config(routing)).place(
        INF, ngc.OutputOptions(min_cut=True)
    )
    return dict(value=r.call_delivered, cut=r.min_cut)


def run_residual_witness(edges, flow, source, sink, value, potentials):
    # The library's minimum cost at the value (CentralizedTE up to `value`), against the given
    # flow's cost; the cycle costs are the given flow's residual structure (a property of the
    # example's flow, read here, not of a library result).
    g = graph_of(edges)
    r = policy(
        g, source, sink, config(ngc.CentralizedTE(source_allocation=S.OPTIMIZED))
    ).place(value)
    assert r.call_delivered == value
    cost = sum(f * e[3] for f, e in zip(flow, edges, strict=True))
    residual, positive = [], []
    for (a, b, cap, c), f in zip(edges, flow, strict=True):
        if f < cap:
            residual.append((a, b, c))
        if f:
            residual.append((b, a, -c))
            positive.append((a, b, c))
    return dict(
        value=value,
        cost=cost,
        minimum_cost_at_value=r.total_cost,
        is_minimum_cost=cost == r.total_cost,
        positive_flow_cycle_costs=cycle_costs(positive),
        residual_cycle_costs=cycle_costs(residual),
    )


def cycle_costs(arcs):
    """The costs of the simple cycles of a small arc list (u, v, cost), each cycle once."""
    nodes = {x for a, b, _ in arcs for x in (a, b)}
    found = set()

    def walk(start, u, seen, ids):
        for i, (a, b, _) in enumerate(arcs):
            if a != u:
                continue
            if b == start:
                cycle = ids + [i]
                found.add(min(tuple(cycle[k:] + cycle[:k]) for k in range(len(cycle))))
            elif b not in seen:
                walk(start, b, seen | {b}, ids + [i])

    for start in sorted(nodes):
        walk(start, start, {start}, [])
    return sorted(sum(arcs[i][2] for i in cycle) for cycle in found)


def run_evaluate_offer(capacities, offer=None):
    # ECMP LOSSY from node 0 over parallel links to node 1 (equal costs: EdgeId = input position):
    # the baseline, each deletion, and the baseline's accounting (requested, offered, the drop on
    # each link).
    k = len(capacities)
    g = ngc.Graph(2, [0] * k, [1] * k, list(capacities), [1] * k, 1.0)
    opts = ngc.EvaluationOptions(config(ngc.ShortestPaths(), LY), offer=offer)
    s = ngc.evaluate_edge_failures(g, 0, 1, opts, 1)
    return dict(
        baseline=s.baseline.call_delivered,
        delivery_reduction=per_edge(
            [(x.edge, x.delivery_reduction) for x in s.entries], k
        ),
        requested=s.baseline.requested,
        offered=s.baseline.offered,
        edge_drops=per_edge(s.baseline.edge_drops, k),
    )


# -- searches and the route choice ------------------------------------------------------------------


def run_search(
    edges, source, destinations, residual=None, min_residual=1, k=None, stretch=None
):
    g = graph_of(edges)
    try:
        if k is not None:
            options = ngc.KspOptions(
                k,
                cost_ceiling=None if stretch is None else ceiling_of(stretch),
                residual=residual,
            )
            found = ngc.ksp(g, source, destinations, options)
            return dict(paths=[[cost, [h[0] for h in p.hops()]] for cost, p in found])
        options = ngc.SpfOptions(
            ngc.Direction.FORWARD,
            targets=[destinations],
            residual=residual,
            min_residual=min_residual,
        )
        r = ngc.spf(g, source, options)
    except (TypeError, ValueError) as e:
        return dict(error=error_text(e))
    if r.target_cost is None or r.target_cost == ngc.UNREACHABLE:
        return dict(distance=None, dag=[])
    return dict(distance=r.target_cost, dag=sorted(r.dag.via_edges().tolist()))


def run_tie_choice(edges, source, destination, available=None):
    # The single path spf chooses under each RouteTieBreak (availability: the residual).
    g = graph_of(edges)
    out = {}
    for tie in ngc.RouteTieBreak.__members__.values():
        options = ngc.SpfOptions(
            ngc.Direction.FORWARD,
            targets=[destination],
            tie_break=tie,
            multipath=False,
            residual=available,
        )
        out[tie.name] = sorted(ngc.spf(g, source, options).dag.via_edges().tolist())
    return out


def run_stretch_values(factors, equal_pairs=()):
    values = []
    for spec in factors:
        try:
            f = ngc.CostStretch(factor_of(spec)).factor
            values.append(f"{f.numerator}/{f.denominator}")
        except TypeError:
            values.append("TypeError")
        except ValueError as e:
            values.append(str(e))
    equal = [
        ngc.CostStretch(factor_of(a)) == ngc.CostStretch(factor_of(b))
        for a, b in equal_pairs
    ]
    return dict(values=values, equal=equal)


# -- configuration, evaluation and texts ------------------------------------------------------------

# 0 -> 1 (e0), 0 -> 2 (e1), 1 -> 2 (e2): one source is [0], two are [0, 1]; the destination is 2.
CONFIG_EDGES = [[0, 1, 5, 1], [0, 2, 5, 1], [1, 2, 5, 1]]
ENUMS = {
    "allocation": ngc.RouteAllocation,
    "source_allocation": ngc.SourceAllocation,
    "weight_basis": ngc.WeightBasis,
    "tie_break": ngc.RouteTieBreak,
}
MODELS = {
    "ShortestPaths": ngc.ShortestPaths,
    "DistributedTE": lambda tunnel_count, **kw: ngc.DistributedTE(tunnel_count, **kw),
    "CentralizedTE": ngc.CentralizedTE,
    "Pinned": lambda program, **kw: ngc.Pinned(program, **kw),
    "FlowPolicyConfig": ngc.FlowPolicyConfig,
}


def config_program(g, spec):
    """A program of the configuration graph: "p" from 0 to 2, or {"sources", "destinations"}."""
    if spec == "p":
        spec = {"sources": [0], "destinations": [2]}
    (d,) = spec["destinations"]
    route = path_of(g, [0 if d == 1 else 1])
    return ngc.ForwardingProgram.from_routes(
        g, spec["sources"], [d], [ngc.RouteDefinition(0, route)]
    )


def config_argument(g, key, value):
    if key in ENUMS and isinstance(value, str):
        return getattr(ENUMS[key], value)
    if key == "tunnel_count":
        count, *up_to = value.split()
        return ngc.TunnelCount(int(count), up_to=up_to == ["up_to"])
    if key == "cost_ceiling":
        return ceiling_of(value)
    if key == "support":
        return dag_of(g, [1])
    if key == "program":
        return config_program(g, value)
    if key == "routing":
        return ngc.ShortestPaths()
    return value


def resolved(routing):
    def name(v):
        return "None" if v is None else v.name

    if isinstance(routing, ngc.ShortestPaths):
        return (
            f"allocation={name(routing.allocation)}"
            f" source_allocation={name(routing.source_allocation)}"
            f" weight_basis={name(routing.weight_basis)}"
            f" tie_break={name(routing.tie_break)}"
            f" multipath={routing.multipath} parallel_edges={routing.parallel_edges}"
        )
    if isinstance(routing, ngc.CentralizedTE):
        t = routing.tunnel_count
        count = "None" if t is None else f"{t.per_head}" + (" up_to" if t.up_to else "")
        return (
            f"source_allocation={name(routing.source_allocation)}"
            f" support={'None' if routing.support is None else 'dag'}"
            f" tunnel_count={count} prefer_headroom={routing.prefer_headroom}"
        )
    return "valid"


def config_case(
    model, kwargs, sources=1, admission="LOSSLESS", cache=False, endpoints=None
):
    g = graph_of(CONFIG_EDGES)
    try:
        value = MODELS[model](
            **{k: config_argument(g, k, v) for k, v in kwargs.items()}
        )
    except TypeError as e:
        text = str(e)
        return (
            f"TypeError: {text}"
            if text == "cost_ceiling must be AbsoluteCost, CostStretch or None"
            else "TypeError"
        )
    if model == "FlowPolicyConfig":
        return "constructed"
    cfg = ngc.FlowPolicyConfig(
        routing=value, admission=getattr(ngc.Admission, admission)
    )
    src, dst = ([0], [2]) if sources == 1 else ([0, 1], [2])
    if endpoints is not None:
        src, dst = endpoints
    tunnels = isinstance(value, ngc.DistributedTE) or (
        isinstance(value, ngc.CentralizedTE) and value.tunnel_count is not None
    )
    try:
        p = ngc.FlowPolicy(
            ngc.FlowLedger(g),
            src,
            dst,
            cfg,
            None,
            ngc.ShortestPathCache(g) if cache else None,
        )
    except ValueError as e:
        return str(e)
    (p.replace if tunnels else p.place)(0)
    return resolved(p.config.routing)


def run_model_config(cases):
    return [config_case(**c) for c in cases]


def evaluation_case(
    model, kwargs, admission="LOSSLESS", offer=None, design_demand=None, sources=1
):
    g = graph_of(CONFIG_EDGES)
    fields = {k: config_argument(g, k, v) for k, v in kwargs.items()}
    if model == "Pinned":
        routing = ngc.Pinned(fields["program"])
    elif model == "DistributedTE":
        routing = ngc.DistributedTE(fields.pop("tunnel_count"), **fields)
    else:
        routing = MODELS[model](**fields)
    try:
        options = ngc.EvaluationOptions(
            config(routing, getattr(ngc.Admission, admission)),
            offer=offer,
            design_demand=design_demand,
        )
    except TypeError:
        return "TypeError"
    src = [0] if sources == 1 else [0, 1]
    try:
        r = ngc.evaluate_demand(g, src, 2, options)
    except ValueError as e:
        return str(e)
    if r.requested == INF and r.by_source[0].offered is not None:
        evaluated = "first-hop experiment"
    else:
        evaluated = r.requested
    design = None if r.plan is None else r.plan.design_demand
    return f"offer={evaluated} design_demand={design}"


def run_evaluation_validation(cases):
    return [evaluation_case(**c) for c in cases]


def run_ownership_texts(objects):
    # Each object of another graph where it is used.
    g = graph_of(CONFIG_EDGES)
    other = graph_of(CONFIG_EDGES)
    mine = config_program(g, "p")
    out = []
    for name in objects:
        try:
            if name == "program":
                policy(g, 0, 2, config(ngc.Pinned(config_program(other, "p"))))
            elif name == "support":
                cte = ngc.CentralizedTE(
                    source_allocation=S.OPTIMIZED, support=dag_of(other, [1])
                )
                policy(g, 0, 2, config(cte))
            elif name == "black_holes":
                policy(
                    g, 0, 2, config(ngc.Pinned(mine, black_holes=other.eligibility()))
                )
            elif name == "Eligibility":
                policy(
                    g,
                    0,
                    2,
                    config(ngc.ShortestPaths()),
                    other.eligibility(excluded_edges=[0]),
                )
            elif name == "ShortestPathCache":
                ngc.FlowPolicy(
                    ngc.FlowLedger(g),
                    0,
                    2,
                    config(ngc.ShortestPaths()),
                    None,
                    ngc.ShortestPathCache(other),
                )
            elif name == "dag":
                ngc.ForwardingProgram.from_dag(g, [0], [2], dag_of(other, [1]))
            out.append(None)
        except ValueError as e:
            out.append(str(e))
    return out


def run_input_texts(cases):
    out = []
    for c in cases:
        caps = c["capacities"]
        k = len(caps)
        g = ngc.Graph(2, [0] * k, [1] * k, list(caps), [1] * k, 1.0)
        row = []
        try:
            routing = ngc.DistributedTE(
                ngc.TunnelCount(1),
                min_unreserved_bandwidth=c["min_unreserved_bandwidth"],
            )
            policy(g, 0, 1, config(routing))
            row.append(None)
        except ValueError as e:
            row.append(str(e))
        for baseline in c["reservation_baseline"]:
            try:
                ngc.FlowLedger(g, reservation_baseline=baseline)
                row.append(None)
            except ValueError as e:
                row.append(str(e))
        out.append(row)
    return out


class IndexOnly:
    """An integer known only through __index__."""

    def __init__(self, n):
        self.n = n

    def __index__(self):
        return self.n


class FloatOnly:
    """A number known only through __float__."""

    def __init__(self, x):
        self.x = x

    def __float__(self):
        return float(self.x)


def input_value(spec):
    """An example's argument: itself, or ["bool", b], ["float", x] ("inf", "-inf" and "nan" as
    text), ["str", s], ["none"], ["fraction", n, d], ["decimal", s], ["floatonly", x] (FloatOnly),
    ["index", n] (IndexOnly), ["pow2", n] (the int 2**n), ["numpy", type name, value], ["list",
    spec, ...] or ["array", dtype, nested list]."""
    if not isinstance(spec, list):
        return spec
    kind, rest = spec[0], spec[1:]
    if kind == "list":
        return [input_value(x) for x in rest]
    if kind == "array":
        return np.array(rest[1], dtype=rest[0])
    if kind == "numpy":
        return getattr(np, rest[0])(rest[1])
    return {
        "bool": lambda: bool(rest[0]),
        "float": lambda: float(rest[0]),
        "str": lambda: rest[0],
        "none": lambda: None,
        "fraction": lambda: Fraction(rest[0], rest[1]),
        "decimal": lambda: Decimal(rest[0]),
        "floatonly": lambda: FloatOnly(rest[0]),
        "index": lambda: IndexOnly(rest[0]),
        "pow2": lambda: 2 ** rest[0],
    }[kind]()


def input_class(call, *args):
    """ "valid", or the input error of call(*args) as "<class>: <text>" (TypeError or ValueError
    exactly)."""
    try:
        call(*args)
    except (TypeError, ValueError) as e:
        assert type(e) in (TypeError, ValueError), type(e)
        return f"{type(e).__name__}: {e}"
    return "valid"


def run_input_classes(cases):
    """Each value at its position: an spf or ksp search option (read at the call), place's offer
    (on a fresh policy), a ledger's service baseline over links of capacities `caps`, or one Place
    of a LOSSY ShortestPaths policy of a fresh ledger in an ORDERED apply_actions list, with the
    case's ledger (the policy's own, or None)."""
    g = ngc.Graph(2, [0], [1], [10], [1], 1.0)
    forward = ngc.Direction.FORWARD

    def baseline(v, caps):
        k = len(caps)
        return ngc.FlowLedger(ngc.Graph(2, [0] * k, [1] * k, caps, [1] * k, 1.0), v)

    def listed(v, _):
        ledger = ngc.FlowLedger(g)
        p = policy(g, 0, 1, config(ngc.ShortestPaths(), LY), ledger=ledger)
        action = ngc.Place(p, input_value(v["offer"]))  # its construction comes first
        own = None if input_value(v["ledger"]) is None else ledger
        return ngc.apply_actions(own, [action])

    calls = {
        "spf min_residual": lambda v, _: ngc.spf(
            g, 0, ngc.SpfOptions(forward, min_residual=v)
        ),
        "ksp k": lambda v, _: ngc.ksp(g, 0, 1, ngc.KspOptions(v)),
        "place offer": lambda v, _: policy(g, 0, 1, config(ngc.ShortestPaths())).place(
            v
        ),
        "service_baseline": baseline,
        "apply_actions": listed,
    }
    out = []
    for case in cases:
        call, caps = calls[case["entry"]], case.get("caps")
        out.append([input_class(call, input_value(v), caps) for v in case["values"]])
    return out


def number_outcome(call, report=False):
    """ "valid" ("valid <result>" when reported), or the error as "<class>: <text>", pybind11's
    refusal of an argument's type as "TypeError: pybind11 conversion"."""
    try:
        got = call()
    except (TypeError, ValueError) as e:
        assert type(e) in (TypeError, ValueError), type(e)
        if type(e) is TypeError and "incompatible" in str(e):
            return "TypeError: pybind11 conversion"
        return f"{type(e).__name__}: {e}"
    if not report:
        return "valid"
    return f"valid {got.tolist() if isinstance(got, np.ndarray) else got}"


def run_number_types(cases):
    """Each value at its position: a resolution (Graph, to_units and as_float, which must agree),
    to_units' physical value at the case's resolution (reported: the units), as_float's integer
    amount, or place's offer (on a fresh policy)."""
    g = ngc.Graph(2, [0], [1], [10], [1], 1.0)

    def resolution(q):
        outcomes = {
            number_outcome(lambda: ngc.Graph(2, [0], [1], [1], [1], q)),
            number_outcome(lambda: ngc.to_units(1.0, q)),
            number_outcome(lambda: ngc.as_float(1, q)),
        }
        assert len(outcomes) == 1, outcomes
        return outcomes.pop()

    out = []
    for case in cases:
        entry = case["entry"]
        row = []
        for spec in case["values"]:
            v = input_value(spec)
            if entry == "resolution":
                row.append(resolution(v))
            elif entry == "to_units":
                q = input_value(case.get("resolution", 1.0))
                row.append(number_outcome(lambda v=v, q=q: ngc.to_units(v, q), True))
            elif entry == "as_float x":
                row.append(number_outcome(lambda v=v: ngc.as_float(v, 1.0)))
            else:
                p = policy(g, 0, 1, config(ngc.ShortestPaths()))
                row.append(number_outcome(lambda v=v, p=p: p.place(v)))
        out.append(row)
    return out


# -- several demands at once (apply_actions) ---------------------------------------------------


def several_demands_layout(cap, routes):
    """The list layout of routes given as link indices: each route is its own demand from its own
    head (node i) to one destination (node n); each example link gets its own two nodes, joined to
    the next link of a route by a link that never limits anything, except where every occurrence
    of e continues on f and every occurrence of f comes from e: there e's head is f's tail (no
    join), so the examples' serial pairs stay pairs. Returns the graph, per example link its
    EdgeId and per route its path of EdgeIds."""
    nxt = [set() for _ in cap]
    prv = [set() for _ in cap]
    for r in routes:
        for j, e in enumerate(r):
            nxt[e].add(r[j + 1] if j + 1 < len(r) else -1)
            prv[e].add(r[j - 1] if j > 0 else -1)

    def fused(e, f):
        return nxt[e] == {f} and prv[f] == {e}

    n = len(routes)
    tail = [n + 1 + e for e in range(len(cap))]
    node = n + 1 + len(cap)
    head = []
    for e in range(len(cap)):
        f = next(iter(nxt[e])) if len(nxt[e]) == 1 else -1
        if f >= 0 and fused(e, f):
            head.append(tail[f])
        else:
            head.append(node)
            node += 1
    rows = [(tail[e], head[e], c, 1) for e, c in enumerate(cap)]
    positions = []
    for i, r in enumerate(routes):
        path = [len(rows)]
        rows.append((i, tail[r[0]], WIDE, 1))
        for j, e in enumerate(r):
            path.append(e)
            if j + 1 < len(r) and not fused(e, r[j + 1]):
                path.append(len(rows))
                rows.append((head[e], tail[r[j + 1]], WIDE, 1))
        path.append(len(rows))
        rows.append((head[r[-1]], n, WIDE, 1))
        positions.append(path)
    g = ngc.Graph(node, *([row[k] for row in rows] for k in range(4)), 1.0)
    ids = [0] * len(rows)
    for e, position in enumerate(g.input_positions().tolist()):
        ids[position] = e
    return g, ids, [[ids[k] for k in path] for path in positions]


def run_list_lossy(offers, routes, cap):
    """Each route its own demand: a Pinned LOSSY policy over one tunnel, created in route order,
    all replaced in one SIMULTANEOUS list. Per demand its arrivals, carried and drops on each link
    of its route (from its records), its delivery, and per example link the total carried."""
    g, ids, paths = several_demands_layout(cap, routes)
    n = len(routes)
    ledger = ngc.FlowLedger(g)
    policies = []
    for i, path in enumerate(paths):
        program = ngc.ForwardingProgram.from_routes(
            g, [i], [n], [ngc.RouteDefinition(i, path_of(g, path))]
        )
        policies.append(pinned(g, program, LY, ledger=ledger))
    actions = [ngc.Replace(p, o) for p, o in zip(policies, offers, strict=True)]
    try:
        results = ngc.apply_actions(ledger, actions, mode=ngc.BatchMode.SIMULTANEOUS)
    except ValueError as e:
        return dict(error=str(e))
    arrivals, carried, edge_drops = [], [], []
    for p, route in zip(policies, routes, strict=True):
        flows, drops = {}, {}
        for rec in p.records():
            for e, x in rec.edge_flows:
                flows[e] = flows.get(e, 0) + x
            for e, x in rec.edge_drops:
                drops[e] = drops.get(e, 0) + x
        carried.append([flows.get(ids[e], 0) for e in route])
        edge_drops.append([drops.get(ids[e], 0) for e in route])
        arrivals.append(
            [c + d for c, d in zip(carried[-1], edge_drops[-1], strict=True)]
        )
    total = ledger.carried().tolist()
    return dict(
        arrivals=arrivals,
        carried=carried,
        edge_drops=edge_drops,
        delivered=[r.call_delivered for r in results],
        edge_flows=[total[ids[e]] for e in range(len(cap))],
    )


def run_stale_list(edges, tunnels, offers, black_links=(), black_nodes=()):
    """Each head its own demand (the examples' labels): a Pinned LOSSY policy over its tunnels with
    the example's black holes, created in head order, all replaced in one SIMULTANEOUS list. The
    statuses come from each policy's evaluation of 0 (they depend on the program and the black
    holes only)."""
    g = graph_of(edges)
    ledger = ngc.FlowLedger(g)
    holes = edge_mask(g, black_links, black_nodes)
    policies, actions = [], []
    for h in sorted(offers):
        mine = [p for th, p in tunnels if th == h]
        program = ngc.ForwardingProgram.from_routes(
            g,
            [h],
            [edges[mine[0][-1]][1]],
            [ngc.RouteDefinition(h, path_of(g, p)) for p in mine],
        )
        policies.append(pinned(g, program, LY, black_holes=holes, ledger=ledger))
        actions.append(ngc.Replace(policies[-1], offers[h]))
    status = [t.status.name for p in policies for t in p.evaluate(0).by_tunnel]
    try:
        results = ngc.apply_actions(ledger, actions, mode=ngc.BatchMode.SIMULTANEOUS)
    except ValueError as e:
        return dict(error=str(e), status=status)
    return dict(status=status, delivered=[r.call_delivered for r in results])


def run_action_list(edges, actions, mode):
    """One-tunnel LOSSY Pinned participants created in list order on a fresh ledger, replaced in
    one list of the mode: per entry its delivery, drops, tunnel status, residual and
    saturated_boundary_edges (from the entry's own source)."""
    g = graph_of(edges)
    ledger = ngc.FlowLedger(g)
    listed = []
    for a in actions:
        route = path_of(g, a["route"])
        program = ngc.ForwardingProgram.from_routes(
            g,
            [a["source"]],
            [route.destination()],
            [ngc.RouteDefinition(a["source"], route)],
        )
        holes = edge_mask(g, a.get("black_holes", ()))
        listed.append(
            ngc.Replace(
                pinned(g, program, LY, black_holes=holes, ledger=ledger), a["offer"]
            )
        )
    outputs = ngc.OutputOptions(residual=True, saturated_boundary_edges=True)
    try:
        results = ngc.apply_actions(
            ledger, listed, mode=getattr(ngc.BatchMode, mode), outputs=outputs
        )
    except ValueError as e:
        return dict(error=str(e))
    entries = [
        dict(
            delivered=r.call_delivered,
            call_edge_drops=[list(x) for x in r.call_edge_drops],
            status=r.by_tunnel[0].status.name,
            residual=r.residual.tolist(),
            saturated_boundary_edges=r.saturated_boundary_edges,
        )
        for r in results
    ]
    return dict(entries=entries)


def run_list_validation(cases):
    """Each case's list with real policies, on one ledger per case over CYC(1, 2) plus a separate
    link 6 -> 7 of 10: an action's policy (by name, created at its first use) has the case's model
    and admission on 6 -> 7, a foreign one is of another ledger; a Pinned action that `computes`
    the period-2 error is the step-3 cycle (two crossing tunnels from 4, weights 1 and 2) with that
    cycle's offer 3 in place of the case's, another Pinned one a tunnel over 6 -> 7; `program` installs the policy's program first (a replace or a
    place of 1), `after_infinite` places +inf first, `foreign_black_holes` gives black holes of
    another graph. The list's first error, or "valid"."""
    rows = [
        (0, 1, 1, 1),
        (1, 2, 1000, 1),
        (1, 5, 1000, 1),
        (2, 3, 2, 1),
        (3, 0, 1000, 1),
        (3, 5, 1000, 1),
        (4, 0, 1000, 1),
        (4, 2, 1000, 1),
        (6, 7, 10, 1),
    ]
    models = {
        "ShortestPaths": lambda: ngc.ShortestPaths(),
        "DistributedTE": lambda: ngc.DistributedTE(ngc.TunnelCount(1)),
        "CentralizedTE": lambda: ngc.CentralizedTE(source_allocation=S.OPTIMIZED),
    }
    out = []
    for case in cases:
        g = graph_of(rows)
        ledger = ngc.FlowLedger(g)
        policies, listed = {}, []
        for a in case["actions"]:
            admission = getattr(ngc.Admission, a["admission"])
            if a["policy"] not in policies:
                own = ngc.FlowLedger(g) if a.get("foreign") else ledger
                if a["model"] != "Pinned":
                    cfg = config(models[a["model"]](), admission)
                    p = policy(g, 6, 7, cfg, ledger=own)
                elif a.get("computes"):
                    cyc = [(4, [6, 0, 1, 3, 5], 1), (4, [7, 3, 4, 0, 2], 2)]
                    routes = [
                        ngc.RouteDefinition(h, path_of(g, es), weight=w)
                        for h, es, w in cyc
                    ]
                    program = ngc.ForwardingProgram.from_routes(g, [4], [5], routes)
                    p = pinned(g, program, admission, ledger=own)
                else:
                    program = ngc.ForwardingProgram.from_routes(
                        g, [6], [7], [ngc.RouteDefinition(6, path_of(g, [8]))]
                    )
                    p = pinned(g, program, admission, ledger=own)
                if a.get("program"):
                    (p.place if admission == LL else p.replace)(1)
                if a.get("after_infinite"):
                    p.place(INF)
                policies[a["policy"]] = p
            p = policies[a["policy"]]
            offer = INF if a["offer"] == "inf" else a["offer"]
            if a.get("computes"):
                offer = 3
            if a["kind"] == "Place":
                listed.append(ngc.Place(p, offer))
            elif a["kind"] == "Replace":
                listed.append(
                    ngc.Replace(p, offer, design_demand=a.get("design_demand"))
                )
            else:
                holes = None
                if a.get("foreign_black_holes"):
                    holes = graph_of(rows).eligibility(excluded_edges=[8])
                listed.append(ngc.Forward(p, offer, black_holes=holes))
        try:
            ngc.apply_actions(ledger, listed, mode=getattr(ngc.BatchMode, case["mode"]))
            out.append("valid")
        except ValueError as e:
            out.append(str(e))
    return out


def run_forwarding_size(cases):
    """The forwarding sizes' counting step on each case's sizes (per participant, groups of
    (tunnels, path length)): the entry count, then the occurrence count, each exact against
    2^31 - 1 (the networks cannot be built at these sizes; test_apply_actions.py shows the same
    refusals end to end under lowered limits). Needs the private hooks."""
    if not hasattr(ext, "_check_forwarding_count"):
        pytest.skip("needs the private hooks (NETGRAPH_CORE_PRIVATE_HOOKS)")
    out = []
    for participants in cases:
        groups = [group for participant in participants for group in participant]
        entries = sum(count for count, _ in groups)
        occurrences = sum(count * length for count, length in groups)
        try:
            ext._check_forwarding_count(0, entries)
            ext._check_forwarding_count(1, occurrences)
            out.append("within the entry and occurrence limits")
        except ValueError as e:
            out.append(str(e))
    return out


def own_highest_averages(x, w):
    """Unit-at-a-time highest averages by exact cross-multiplication, ties to the earlier label:
    the test's own, sharing nothing with the library."""
    a = [0] * len(w)
    for _ in range(x):
        best = None
        for i, wi in enumerate(w):
            if wi > 0 and (best is None or wi * (a[best] + 1) > w[best] * (a[i] + 1)):
                best = i
        if best is None:
            return [0] * len(w)
        a[best] += 1
    return a


def run_c3(offers, routes, cap):
    """The C3 instance (A: links 0, 1; B: 1, 2; C: 2, 0): every integer carried record within the
    continuation bounds (a demand's carried on a link at most its arrival there) and work
    conservation (each link carries min(capacity, arrivals)), each checked against the link
    equations with the test's own highest averages; the library's highest averages (private hooks)
    agrees with it on every equation evaluated; and the list itself raises the period-2 error."""
    assert offers == [1, 1, 2] and routes == [[0, 1], [1, 2], [2, 0]]
    records = solutions = 0
    violations = []
    for a0, a1, b1, b2, c2, c0 in itertools.product(
        range(2), range(2), range(2), range(2), range(3), range(3)
    ):
        if a1 > a0 or b2 > b1 or c0 > c2:
            continue
        # Arrivals: A 1 then a0; B 1 then b1; C 2 then c2. Links: 0 (A, C), 1 (A, B), 2 (B, C).
        links = [([1, c2], [a0, c0]), ([a0, 1], [a1, b1]), ([b1, 2], [b2, c2])]
        if any(
            sum(car) != min(v, sum(arr))
            for (arr, car), v in zip(links, cap, strict=True)
        ):
            continue
        records += 1
        bad = []
        for k, ((arr, car), v) in enumerate(zip(links, cap, strict=True)):
            t = min(v, sum(arr))
            mine = own_highest_averages(t, arr)
            if hasattr(ext, "_highest_averages") and all(arr):
                assert ext._highest_averages(t, arr) == mine, (arr, t)
            if mine != car:
                bad.append(k)
        solutions += not bad
        violations.append(bad)
    assert violations == [[2], [0], [0, 1, 2], [1], [0]], violations
    assert run_list_lossy(offers, routes, cap) == dict(
        error="lossy forwarding iteration entered a period-2 cycle"
    )
    return dict(feasible_work_conserving_records=records, local_H_solutions=solutions)


# The examples of several demands at once that share an oracle with a single program's: each
# route (or head) is its own demand, all in one SIMULTANEOUS list.
SEVERAL_DEMANDS = {
    "D-crossing": run_list_lossy,
    "D-simultaneous": run_list_lossy,
    "D-C3": run_list_lossy,
    "D-cycle-prefix-suffix": run_list_lossy,
    "D-cycle-subdivision-before": run_list_lossy,
    "D-cycle-subdivision-after": run_list_lossy,
    "STALE-cycle-period-2": run_stale_list,
    "STALE-cycle-round-limit-baseline": run_stale_list,
    "STALE-cycle-round-limit": run_stale_list,
    "STALE-black-hole-in-cycle": run_stale_list,
}


ADAPTERS = {
    "ha": run_ha,
    "split_error": run_split_error,
    "bound_sweep": run_bound_sweep,
    "member_order": run_member_order,
    "pooled": run_pooled,
    "ledger": run_ledger,
    "refusal": run_refusal,
    "increments": run_increments,
    "epoch": run_epoch,
    "equal_increment": run_equal_increment,
    "zero_state": run_zero_state,
    "change_flags": run_change_flags,
    "lossy": run_lossy,
    "bound": run_bound,
    "split_contexts": run_split_contexts,
    "dag_contexts": run_dag_contexts,
    "stale_dag": run_stale_dag,
    "stale_paths": run_stale_paths,
    "termination": run_termination,
    "t1": run_t1,
    "head_requests": run_head_requests,
    "t2": run_t2,
    "upto": run_upto,
    "flow": run_flow,
    "min_cut": run_min_cut,
    "residual_witness": run_residual_witness,
    "evaluate_offer": run_evaluate_offer,
    "search": run_search,
    "tie_choice": run_tie_choice,
    "stretch_values": run_stretch_values,
    "model_config": run_model_config,
    "evaluation_validation": run_evaluation_validation,
    "ownership_texts": run_ownership_texts,
    "input_texts": run_input_texts,
    "input_classes": run_input_classes,
    "action_list": run_action_list,
    "list_validation": run_list_validation,
    "forwarding_size": run_forwarding_size,
    "c3": run_c3,
    "number_types": run_number_types,
}


EXAMPLES = [
    {
        "id": "E-configured",
        "oracle": "ha",
        "input": {"n": 12, "weights": [1, 2]},
        "expected": [4, 8],
    },
    {
        "id": "E-snapshot",
        "oracle": "ha",
        "input": {"n": 12, "weights": [2, 10]},
        "expected": [2, 10],
    },
    {
        "id": "E-source-survivor",
        "oracle": "ha",
        "input": {"n": 16, "weights": [3, 1]},
        "expected": [12, 4],
    },
    {
        "id": "H-not-one-unit",
        "oracle": "ha",
        "input": {"n": 50, "weights": [100, 1, 1, 1]},
        "expected": [50, 0, 0, 0],
    },
    {
        "id": "R3-H-large",
        "oracle": "ha",
        "input": {
            "n": 85070591730234615865843651857942052863,
            "weights": [
                42535295865117307932921825928971026432,
                42535295865117307932921825928971026432,
            ],
        },
        "expected": [
            42535295865117307932921825928971026432,
            42535295865117307932921825928971026431,
        ],
    },
    {
        "id": "H-local-bound",
        "oracle": "split_error",
        "input": {"n": 1000000, "weights": [1000000, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]},
        "expected": {
            "assigned": [1000000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            "max_excess": "1000000/100001",
            "max_deficit": "100000/100001",
            "below_k": True,
        },
    },
    {
        "id": "RR-local-bound",
        "oracle": "split_error",
        "input": {"n": 10, "weights": [1, 1, 1]},
        "expected": {
            "assigned": [4, 3, 3],
            "max_excess": "2/3",
            "max_deficit": "1/3",
            "below_one": True,
        },
    },
    {
        "id": "H-bound-sweep",
        "oracle": "bound_sweep",
        "input": {"max_total": 24, "max_members": 4, "max_weight": 4},
        "expected": {"cases": 19400, "violations": 0, "equal_violations": 0},
    },
    {
        "id": "R3-member-order",
        "oracle": "member_order",
        "input": {"amount": 1, "order": [1, 0]},
        "expected": [0, 1],
    },
    {
        "id": "B-pooled",
        "oracle": "pooled",
        "input": {"sources": [1, 1], "members": 2},
        "expected": {
            "columns": [1, 1],
            "source_by_column": [[1, 0], [0, 1]],
            "independent_columns": [2, 0],
        },
    },
    {
        "id": "SPF-zero-edge",
        "oracle": "search",
        "input": {
            "edges": [[0, 1, 0, 1], [0, 2, 5, 1], [2, 1, 5, 1]],
            "source": 0,
            "destinations": 1,
        },
        "expected": {"distance": 1, "dag": [0]},
    },
    {
        "id": "SPF-zero-edge-residual",
        "oracle": "search",
        "input": {
            "edges": [[0, 1, 0, 1], [0, 2, 5, 1], [2, 1, 5, 1]],
            "source": 0,
            "destinations": 1,
            "residual": [0, 5, 5],
        },
        "expected": {"distance": 2, "dag": [1, 2]},
    },
    {
        "id": "KSP-zero-edge",
        "oracle": "search",
        "input": {
            "edges": [[0, 1, 0, 1], [0, 2, 5, 1], [2, 1, 5, 1]],
            "source": 0,
            "destinations": 1,
            "k": 2,
        },
        "expected": {"paths": [[1, [0]], [2, [1, 2]]]},
    },
    {
        "id": "KSP-zero-edge-residual",
        "oracle": "search",
        "input": {
            "edges": [[0, 1, 0, 1], [0, 2, 5, 1], [2, 1, 5, 1]],
            "source": 0,
            "destinations": 1,
            "k": 2,
            "residual": [0, 5, 5],
        },
        "expected": {"paths": [[2, [1, 2]]]},
    },
    {
        "id": "B1-equal",
        "oracle": "flow",
        "input": {
            "edges": [[0, 3, 2, 1], [1, 2, 10, 1], [2, 3, 10, 2]],
            "heads": [0, 1],
            "destinations": [3],
            "limit": 12,
            "source_allocation": "EQUAL",
        },
        "expected": {"cost": 8, "delivered": 4, "edge_flows": [2, 2, 2]},
    },
    {
        "id": "B1-optimized",
        "oracle": "flow",
        "input": {
            "edges": [[0, 3, 2, 1], [1, 2, 10, 1], [2, 3, 10, 2]],
            "heads": [0, 1],
            "destinations": [3],
            "limit": 12,
            "source_allocation": "OPTIMIZED",
        },
        "expected": {"cost": 32, "delivered": 12, "edge_flows": [2, 10, 10]},
    },
    {
        "id": "MINCUT-B1",
        "oracle": "min_cut",
        "input": {
            "edges": [[0, 3, 2, 1], [1, 2, 10, 1], [2, 3, 10, 2]],
            "heads": [0, 1],
            "destinations": [3],
        },
        "expected": {"value": 12, "cut": [0, 1]},
    },
    {
        "id": "MINCUT-zero-edge",
        "oracle": "min_cut",
        "input": {
            "edges": [[0, 2, 0, 1], [1, 2, 10, 1]],
            "heads": [0, 1],
            "destinations": [2],
        },
        "expected": {"value": 10, "cut": [0, 1]},
    },
    {
        "id": "B2-equal-global",
        "oracle": "flow",
        "input": {
            "edges": [[0, 1, 4, 1], [0, 2, 2, 1], [1, 2, 10, 1]],
            "heads": [0, 1],
            "destinations": [2],
            "limit": 12,
        },
        "expected": {"cost": 16, "delivered": 12, "edge_flows": [4, 2, 10]},
    },
    {
        "id": "B2-shortest",
        "oracle": "flow",
        "input": {
            "edges": [[0, 2, 2, 1], [1, 2, 10, 1]],
            "heads": [0, 1],
            "destinations": [2],
            "limit": 12,
            "source_allocation": "EQUAL",
        },
        "expected": {"cost": 4, "delivered": 4, "edge_flows": [2, 2]},
    },
    {
        "id": "B2-optimized",
        "oracle": "flow",
        "input": {
            "edges": [[0, 1, 4, 1], [0, 2, 2, 1], [1, 2, 10, 1]],
            "heads": [0, 1],
            "destinations": [2],
            "limit": 12,
            "source_allocation": "OPTIMIZED",
        },
        "expected": {"cost": 12, "delivered": 12, "edge_flows": [0, 2, 10]},
    },
    {
        "id": "B3-zero-first",
        "oracle": "flow",
        "input": {
            "edges": [[0, 2, 0, 1], [1, 2, 10, 1]],
            "heads": [0, 1],
            "destinations": [2],
            "limit": 10,
        },
        "expected": {"delivered": 0, "cost": 0, "edge_flows": [0, 0]},
    },
    {
        "id": "B3-positive-first",
        "oracle": "flow",
        "input": {
            "edges": [[0, 2, 10, 1], [1, 2, 0, 1]],
            "heads": [0, 1],
            "destinations": [2],
            "limit": 10,
        },
        "expected": {"delivered": 1, "cost": 1, "edge_flows": [1, 0]},
    },
    {
        "id": "B3-delete-zero",
        "oracle": "flow",
        "input": {
            "edges": [[0, 2, 0, 1], [1, 2, 10, 1]],
            "heads": [1],
            "destinations": [2],
            "limit": 10,
        },
        "expected": {"delivered": 10, "cost": 10, "edge_flows": [0, 10]},
    },
    {
        "id": "E-support",
        "oracle": "flow",
        "input": {
            "edges": [[0, 1, 10, 1], [0, 2, 10, 1], [1, 3, 2, 1], [2, 3, 10, 2]],
            "heads": [0],
            "destinations": [3],
            "limit": 12,
        },
        "expected": {"cost": 34, "delivered": 12, "edge_flows": [2, 10, 2, 10]},
    },
    {
        "id": "E-support-shortest",
        "oracle": "flow",
        "input": {
            "edges": [[0, 1, 10, 1], [1, 3, 2, 1]],
            "heads": [0],
            "destinations": [3],
            "limit": 12,
        },
        "expected": {"delivered": 2, "cost": 4, "edge_flows": [2, 2]},
    },
    {
        "id": "T2-OPT-global",
        "oracle": "flow",
        "input": {
            "edges": [[0, 2, 2, 1], [1, 2, 10, 1]],
            "heads": [0, 1],
            "destinations": [2],
            "limit": 12,
            "source_allocation": "OPTIMIZED",
        },
        "expected": {"delivered": 12, "cost": 12, "edge_flows": [2, 10]},
    },
    {
        "id": "R4-lossless-equal",
        "oracle": "equal_increment",
        "input": {"capacities": [2, 10], "costs": [9, 1], "calls": [3, 3]},
        "expected": {
            "call_deliveries": [3, 1],
            "allocation": [2, 2],
            "delivered": 4,
            "cost": 20,
        },
    },
    {
        "id": "R3-residual-optimal",
        "oracle": "residual_witness",
        "input": {
            "edges": [[0, 1, 1, 1], [0, 2, 1, 1], [1, 0, 1, 1]],
            "flow": [0, 1, 0],
            "source": 0,
            "sink": 2,
            "value": 1,
            "potentials": [0, 0, 1],
        },
        "expected": {
            "is_minimum_cost": True,
            "cost": 1,
            "residual_cycle_costs": [2],
            "positive_flow_cycle_costs": [],
        },
    },
    {
        "id": "R3-residual-cyclic",
        "oracle": "residual_witness",
        "input": {
            "edges": [[0, 1, 1, 1], [0, 2, 1, 1], [1, 0, 1, 1]],
            "flow": [1, 1, 1],
            "source": 0,
            "sink": 2,
            "value": 1,
            "potentials": [0, 0, 1],
        },
        "expected": {
            "is_minimum_cost": False,
            "cost": 3,
            "minimum_cost_at_value": 1,
            "residual_cycle_costs": [-2],
            "positive_flow_cycle_costs": [2],
        },
    },
    {
        "id": "A-X1",
        "oracle": "t1",
        "input": {
            "edges": [[0, 1, 6, 1], [0, 1, 6, 2]],
            "heads": [0],
            "destinations": [1],
            "design_demand": 12,
            "per_head": 4,
            "offer": 13,
        },
        "expected": {
            "lossy_delivery": 12,
            "lossy_edge_drops": [1, 0],
            "lossless_admitted": 12,
            "lossless_unadmitted": 1,
            "bound": 12,
            "no_drop_load": [7, 6],
            "offered_by_entry": [4, 3, 3, 3],
            "requested_bandwidth": [3, 3, 3, 3],
            "total_reserved_bandwidth": 12,
            "edge_reservations": [6, 6],
            "routes": [[0], [0], [1], [1]],
            "source_drops": 0,
            "up": 4,
        },
    },
    {
        "id": "A-X2",
        "oracle": "t1",
        "input": {
            "edges": [[0, 1, 3, 1], [0, 1, 2, 2]],
            "heads": [0],
            "destinations": [1],
            "design_demand": 5,
            "per_head": 2,
            "offer": 5,
        },
        "expected": {
            "bound": 5,
            "no_drop_load": [3, 2],
            "offered_by_entry": [3, 2],
            "requested_bandwidth": [3, 2],
            "total_reserved_bandwidth": 5,
            "edge_reservations": [3, 2],
            "routes": [[0], [1]],
            "source_drops": 0,
            "up": 2,
        },
    },
    {
        "id": "A-X2-small-share",
        "oracle": "t1",
        "input": {
            "edges": [[0, 1, 6, 1], [0, 1, 6, 2]],
            "heads": [0],
            "destinations": [1],
            "design_demand": 1,
            "per_head": 2,
            "offer": 1,
        },
        "expected": {
            "bound": 6,
            "no_drop_load": [1, 0],
            "offered_by_entry": [1, 0],
            "requested_bandwidth": [1, 0],
            "total_reserved_bandwidth": 1,
            "edge_reservations": [1, 0],
            "routes": [[0], [0]],
            "status": ["UP", "UP"],
            "source_drops": 0,
            "up": 2,
        },
    },
    {
        "id": "T1-over-reservation",
        "oracle": "t1",
        "input": {
            "edges": [[0, 1, 5, 1], [0, 1, 5, 1]],
            "heads": [0],
            "destinations": [1],
            "design_demand": 8,
            "per_head": 2,
            "offer": 10,
        },
        "expected": {
            "requested_bandwidth": [4, 4],
            "routes": [[0], [1]],
            "edge_reservations": [4, 4],
            "offered_by_entry": [5, 5],
            "no_drop_load": [5, 5],
            "lossy_delivery": 10,
            "lossy_edge_drops": [0, 0],
            "bound": 10,
        },
    },
    {
        "id": "T1-zero-link-drop",
        "oracle": "t1",
        "input": {
            "edges": [[0, 1, 0, 1]],
            "heads": [0],
            "destinations": [1],
            "design_demand": 0,
            "per_head": 1,
            "offer": 1,
        },
        "expected": {
            "requested_bandwidth": [0],
            "routes": [[0]],
            "status": ["UP"],
            "offered_by_entry": [1],
            "source_drops": 0,
            "lossy_delivery": 0,
            "lossy_edge_drops": [1],
            "bound": 0,
            "lossless_admitted": 0,
            "lossless_unadmitted": 1,
        },
    },
    {
        "id": "A-X3",
        "oracle": "t1",
        "input": {
            "edges": [[0, 1, 2, 1], [0, 1, 10, 2]],
            "heads": [0],
            "destinations": [1],
            "design_demand": 12,
            "per_head": 4,
            "offer": 12,
        },
        "expected": {
            "bound": 10,
            "no_drop_load": [0, 12],
            "offered_by_entry": [4, 4, 4, 0],
            "requested_bandwidth": [3, 3, 3, 3],
            "total_reserved_bandwidth": 9,
            "edge_reservations": [0, 9],
            "routes": [[1], [1], [1], []],
            "source_drops": 0,
            "up": 3,
            "status": ["UP", "UP", "UP", "INSUFFICIENT_BANDWIDTH"],
        },
    },
    {
        "id": "A-X4",
        "oracle": "t1",
        "input": {
            "edges": [[0, 2, 6, 1], [1, 2, 6, 1], [2, 3, 6, 1]],
            "heads": [0, 1],
            "destinations": [3],
            "design_demand": 12,
            "per_head": 2,
            "offer": 12,
        },
        "expected": {
            "bound": 6,
            "no_drop_load": [6, 6, 12],
            "offered_by_entry": [6, 6, 0, 0],
            "requested_bandwidth": [3, 3, 3, 3],
            "total_reserved_bandwidth": 6,
            "edge_reservations": [3, 3, 6],
            "routes": [[0, 2], [1, 2], [], []],
            "source_drops": 0,
            "up": 2,
        },
    },
    {
        "id": "A-X5",
        "oracle": "t1",
        "input": {
            "edges": [[0, 3, 20, 1], [1, 3, 20, 1], [2, 3, 20, 1]],
            "heads": [0, 1, 2],
            "destinations": [3],
            "design_demand": 2,
            "per_head": 1,
            "offer": 3,
        },
        "expected": {
            "bound": 60,
            "no_drop_load": [1, 1, 1],
            "offered_by_entry": [1, 1, 1],
            "requested_bandwidth": [1, 1, 0],
            "total_reserved_bandwidth": 2,
            "edge_reservations": [1, 1, 0],
            "routes": [[0], [1], [2]],
            "source_drops": 0,
            "up": 3,
        },
    },
    {
        "id": "T1-split",
        "oracle": "head_requests",
        "input": {"heads": 2, "per_head": 2, "design_demand": 5},
        "expected": [2, 1, 1, 1],
    },
    {
        "id": "A-X6",
        "oracle": "t1",
        "input": {
            "edges": [[0, 1, 5, 1], [0, 1, 5, 1]],
            "heads": [0],
            "destinations": [1],
            "design_demand": 8,
            "per_head": 2,
            "offer": 8,
        },
        "expected": {
            "bound": 10,
            "no_drop_load": [4, 4],
            "offered_by_entry": [4, 4],
            "requested_bandwidth": [4, 4],
            "total_reserved_bandwidth": 8,
            "edge_reservations": [4, 4],
            "routes": [[0], [1]],
            "source_drops": 0,
            "up": 2,
        },
    },
    {
        "id": "A-X8",
        "oracle": "t1",
        "input": {
            "edges": [[0, 2, 0, 1], [1, 2, 8, 3]],
            "heads": [0, 1],
            "destinations": [2],
            "design_demand": 4,
            "per_head": 1,
            "offer": 4,
        },
        "expected": {
            "bound": 0,
            "no_drop_load": [0, 2],
            "offered_by_entry": [0, 2],
            "requested_bandwidth": [2, 2],
            "total_reserved_bandwidth": 2,
            "edge_reservations": [0, 2],
            "routes": [[], [1]],
            "source_drops": 2,
            "up": 1,
        },
    },
    {
        "id": "A-X9",
        "oracle": "t1",
        "input": {
            "edges": [
                [0, 1, 8192, 1],
                [0, 2, 40960, 1],
                [1, 3, 8192, 1],
                [2, 3, 40960, 1],
            ],
            "heads": [0],
            "destinations": [3],
            "design_demand": 49152,
            "per_head": 16,
            "offer": 49152,
        },
        "expected": {
            "bound": 47262,
            "no_drop_load": [6554, 42598, 6554, 42598],
            "offered_by_entry": [
                3277,
                3277,
                3277,
                3277,
                3277,
                3277,
                3277,
                3277,
                3277,
                3277,
                3277,
                3277,
                3276,
                3276,
                3276,
                0,
            ],
            "requested_bandwidth": [
                3072,
                3072,
                3072,
                3072,
                3072,
                3072,
                3072,
                3072,
                3072,
                3072,
                3072,
                3072,
                3072,
                3072,
                3072,
                3072,
            ],
            "total_reserved_bandwidth": 46080,
            "edge_reservations": [6144, 39936, 6144, 39936],
            "routes": [
                [0, 2],
                [0, 2],
                [1, 3],
                [1, 3],
                [1, 3],
                [1, 3],
                [1, 3],
                [1, 3],
                [1, 3],
                [1, 3],
                [1, 3],
                [1, 3],
                [1, 3],
                [1, 3],
                [1, 3],
                [],
            ],
            "source_drops": 0,
            "up": 15,
        },
    },
    {
        "id": "A-X9-large-design",
        "oracle": "t1",
        "input": {
            "edges": [
                [0, 1, 8192, 1],
                [0, 2, 40960, 1],
                [1, 3, 8192, 1],
                [2, 3, 40960, 1],
            ],
            "heads": [0],
            "destinations": [3],
            "design_demand": 409600,
            "per_head": 16,
        },
        "expected": {
            "bound": 40960,
            "requested_bandwidth": [
                25600,
                25600,
                25600,
                25600,
                25600,
                25600,
                25600,
                25600,
                25600,
                25600,
                25600,
                25600,
                25600,
                25600,
                25600,
                25600,
            ],
            "total_reserved_bandwidth": 25600,
            "edge_reservations": [0, 25600, 0, 25600],
            "routes": [
                [1, 3],
                [],
                [],
                [],
                [],
                [],
                [],
                [],
                [],
                [],
                [],
                [],
                [],
                [],
                [],
                [],
            ],
            "up": 1,
        },
    },
    {
        "id": "A-X8-delete-zero-edge",
        "oracle": "t1",
        "input": {
            "edges": [[0, 2, 0, 1], [1, 2, 8, 3]],
            "heads": [0, 1],
            "destinations": [2],
            "design_demand": 4,
            "per_head": 1,
            "excluded": [0],
        },
        "expected": {
            "bound": 8,
            "requested_bandwidth": [4],
            "total_reserved_bandwidth": 4,
            "edge_reservations": [0, 4],
            "routes": [[1]],
            "up": 1,
        },
    },
    {
        "id": "A-X8-Stretch2",
        "oracle": "t1",
        "input": {
            "edges": [[0, 2, 0, 1], [1, 2, 8, 3]],
            "heads": [0, 1],
            "destinations": [2],
            "design_demand": 4,
            "per_head": 1,
            "stretch": Fraction(2, 1),
        },
        "expected": {
            "bound": 0,
            "requested_bandwidth": [2, 2],
            "total_reserved_bandwidth": 0,
            "edge_reservations": [0, 0],
            "routes": [[], []],
            "up": 0,
            "status": ["INSUFFICIENT_BANDWIDTH", "COST_CEILING_EXCEEDED"],
        },
    },
    {
        "id": "A-X8-Absolute3",
        "oracle": "t1",
        "input": {
            "edges": [[0, 2, 0, 1], [1, 2, 8, 3]],
            "heads": [0, 1],
            "destinations": [2],
            "design_demand": 4,
            "per_head": 1,
            "absolute": 3,
        },
        "expected": {
            "bound": 0,
            "requested_bandwidth": [2, 2],
            "total_reserved_bandwidth": 2,
            "edge_reservations": [0, 2],
            "routes": [[], [1]],
            "up": 1,
            "status": ["INSUFFICIENT_BANDWIDTH", "UP"],
        },
    },
    {
        "id": "A-X8-Stretch-exact",
        "oracle": "t1",
        "input": {
            "edges": [[0, 1, 0, 5], [0, 1, 1, 6]],
            "heads": [0],
            "destinations": [1],
            "design_demand": 1,
            "per_head": 1,
            "offer": 1,
            "stretch": 1.2,
        },
        "expected": {
            "status": ["UP"],
            "routes": [[1]],
            "up": 1,
            "lossy_delivery": 1,
            "source_drops": 0,
            "lossless_admitted": 1,
        },
    },
    {
        "id": "A-X8-Stretch-cost7",
        "oracle": "t1",
        "input": {
            "edges": [[0, 1, 0, 5], [0, 1, 1, 7]],
            "heads": [0],
            "destinations": [1],
            "design_demand": 1,
            "per_head": 1,
            "offer": 1,
            "stretch": Fraction(6, 5),
        },
        "expected": {
            "status": ["COST_CEILING_EXCEEDED"],
            "routes": [[]],
            "up": 0,
            "lossy_delivery": 0,
            "source_drops": 1,
            "lossless_admitted": 0,
        },
    },
    {
        "id": "A-X7-8",
        "oracle": "t1",
        "input": {
            "edges": [
                [0, 1, 6, 1],
                [0, 2, 6, 3],
                [1, 2, 4, 1],
                [1, 3, 6, 3],
                [2, 3, 6, 1],
            ],
            "heads": [0],
            "destinations": [3],
            "design_demand": 8,
            "per_head": 2,
        },
        "expected": {
            "requested_bandwidth": [4, 4],
            "routes": [[0, 2, 4], []],
            "up": 1,
            "total_reserved_bandwidth": 4,
            "edge_reservations": [4, 0, 4, 0, 4],
            "bound": 4,
        },
    },
    {
        "id": "A-X7-10",
        "oracle": "t1",
        "input": {
            "edges": [
                [0, 1, 6, 1],
                [0, 2, 6, 3],
                [1, 2, 4, 1],
                [1, 3, 6, 3],
                [2, 3, 6, 1],
            ],
            "heads": [0],
            "destinations": [3],
            "design_demand": 10,
            "per_head": 2,
        },
        "expected": {
            "requested_bandwidth": [5, 5],
            "routes": [[0, 3], [1, 4]],
            "up": 2,
            "total_reserved_bandwidth": 10,
            "edge_reservations": [5, 5, 0, 5, 5],
            "bound": 12,
        },
    },
    {
        "id": "A-X7-14",
        "oracle": "t1",
        "input": {
            "edges": [
                [0, 1, 6, 1],
                [0, 2, 6, 3],
                [1, 2, 4, 1],
                [1, 3, 6, 3],
                [2, 3, 6, 1],
            ],
            "heads": [0],
            "destinations": [3],
            "design_demand": 14,
            "per_head": 2,
        },
        "expected": {
            "requested_bandwidth": [7, 7],
            "routes": [[], []],
            "up": 0,
            "total_reserved_bandwidth": 0,
            "edge_reservations": [0, 0, 0, 0, 0],
            "bound": 0,
        },
    },
    {
        "id": "T1-two-heads",
        "oracle": "t1",
        "input": {
            "edges": [[0, 2, 2, 1], [1, 2, 10, 1]],
            "heads": [0, 1],
            "destinations": [2],
            "design_demand": 12,
            "per_head": 2,
            "offer": 12,
        },
        "expected": {
            "requested_bandwidth": [3, 3, 3, 3],
            "routes": [[], [1], [], [1]],
            "status": ["INSUFFICIENT_BANDWIDTH", "UP", "INSUFFICIENT_BANDWIDTH", "UP"],
            "up": 2,
            "total_reserved_bandwidth": 6,
            "edge_reservations": [0, 6],
            "offered_by_entry": [0, 3, 0, 3],
            "source_drops": 6,
            "no_drop_load": [0, 6],
            "bound": 0,
            "lossy_delivery": 6,
        },
    },
    {
        "id": "T1-idle",
        "oracle": "t1",
        "input": {
            "edges": [[0, 1, 10, 1]],
            "heads": [0],
            "destinations": [1],
            "design_demand": 3,
            "per_head": 4,
            "offer": 3,
        },
        "expected": {
            "requested_bandwidth": [1, 1, 1, 0],
            "routes": [[0], [0], [0], [0]],
            "up": 4,
            "total_reserved_bandwidth": 3,
            "edge_reservations": [3],
            "offered_by_entry": [1, 1, 1, 0],
            "source_drops": 0,
            "no_drop_load": [3],
            "bound": 10,
            "lossy_delivery": 3,
        },
    },
    {
        "id": "T1-idle-zero-link",
        "oracle": "t1",
        "input": {
            "edges": [[0, 1, 0, 1], [0, 1, 10, 2]],
            "heads": [0],
            "destinations": [1],
            "design_demand": 1,
            "per_head": 2,
            "offer": 2,
        },
        "expected": {
            "requested_bandwidth": [1, 0],
            "routes": [[1], [0]],
            "status": ["UP", "UP"],
            "up": 2,
            "total_reserved_bandwidth": 1,
            "edge_reservations": [0, 1],
            "offered_by_entry": [1, 1],
            "source_drops": 0,
            "no_drop_load": [1, 1],
            "bound": 1,
            "lossy_delivery": 1,
        },
    },
    {
        "id": "T1-zero-design",
        "oracle": "t1",
        "input": {
            "edges": [[0, 1, 10, 1]],
            "heads": [0],
            "destinations": [1],
            "design_demand": 0,
            "per_head": 2,
            "offer": 4,
        },
        "expected": {
            "requested_bandwidth": [0, 0],
            "routes": [[0], [0]],
            "up": 2,
            "total_reserved_bandwidth": 0,
            "edge_reservations": [0],
            "offered_by_entry": [2, 2],
            "source_drops": 0,
            "no_drop_load": [4],
            "bound": 10,
            "lossy_delivery": 4,
        },
    },
    {
        "id": "C-G1",
        "oracle": "t2",
        "input": {
            "prefer_headroom": True,
            "edges": [
                [0, 1, 8192, 1],
                [0, 2, 40960, 1],
                [1, 3, 8192, 1],
                [2, 3, 40960, 1],
            ],
            "heads": [0],
            "destinations": [3],
            "per_head": 16,
            "design_demand": 49152,
            "lossy_offer": 49152,
        },
        "expected": {
            "max_bandwidth_per_tunnel": 2925,
            "counts": [2, 14, 2, 14],
            "max_fractional_bandwidth_per_tunnel": "20480/7",
            "max_lossless_offer": 46812,
            "lossy_delivery": 47104,
            "paths": [
                [0, 2],
                [0, 2],
                [1, 3],
                [1, 3],
                [1, 3],
                [1, 3],
                [1, 3],
                [1, 3],
                [1, 3],
                [1, 3],
                [1, 3],
                [1, 3],
                [1, 3],
                [1, 3],
                [1, 3],
                [1, 3],
            ],
            "headroom_penalty": 73610035200,
            "total_planned_bandwidth": 46800,
            "total_tunnel_cost": 32,
            "planned_bandwidth_per_tunnel": 2925,
        },
    },
    {
        "id": "C-G1-finite",
        "oracle": "t2",
        "input": {
            "prefer_headroom": True,
            "edges": [
                [0, 1, 8192, 1],
                [0, 2, 40960, 1],
                [1, 3, 8192, 1],
                [2, 3, 40960, 1],
            ],
            "heads": [0],
            "destinations": [3],
            "per_head": 16,
            "design_demand": 8192,
        },
        "expected": {
            "max_bandwidth_per_tunnel": 2925,
            "counts": [2, 14, 2, 14],
            "max_fractional_bandwidth_per_tunnel": "20480/7",
            "max_lossless_offer": 46812,
            "paths": [
                [0, 2],
                [0, 2],
                [1, 3],
                [1, 3],
                [1, 3],
                [1, 3],
                [1, 3],
                [1, 3],
                [1, 3],
                [1, 3],
                [1, 3],
                [1, 3],
                [1, 3],
                [1, 3],
                [1, 3],
                [1, 3],
            ],
            "headroom_penalty": 12884901876,
            "total_planned_bandwidth": 8192,
            "total_tunnel_cost": 32,
            "planned_bandwidth_per_tunnel": 512,
        },
    },
    {
        "id": "C-G4",
        "oracle": "t2",
        "input": {
            "prefer_headroom": True,
            "edges": [
                [0, 1, 4096, 1],
                [0, 2, 4096, 1],
                [0, 3, 4096, 1],
                [0, 4, 4096, 1],
                [0, 5, 4096, 1],
                [0, 6, 4096, 1],
                [0, 7, 4096, 1],
                [0, 8, 4096, 1],
                [0, 9, 4096, 1],
                [0, 10, 4096, 1],
                [0, 11, 4096, 1],
                [0, 12, 4096, 1],
                [0, 13, 4096, 1],
                [0, 14, 4096, 1],
                [0, 15, 4096, 1],
                [1, 16, 4096, 1],
                [2, 16, 4096, 1],
                [3, 16, 4096, 1],
                [4, 16, 4096, 1],
                [5, 16, 4096, 1],
                [6, 16, 4096, 1],
                [7, 16, 4096, 1],
                [8, 16, 4096, 1],
                [9, 16, 4096, 1],
                [10, 16, 4096, 1],
                [11, 16, 4096, 1],
                [12, 16, 4096, 1],
                [13, 16, 4096, 1],
                [14, 16, 4096, 1],
                [15, 16, 4096, 1],
            ],
            "heads": [0],
            "destinations": [16],
            "per_head": 16,
            "design_demand": 61440,
            "lossy_offer": 61440,
        },
        "expected": {
            "max_bandwidth_per_tunnel": 2048,
            "planned_bandwidth_per_tunnel": 2048,
            "max_fractional_bandwidth_per_tunnel": "2048",
            "counts": [
                2,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                2,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
            ],
            "paths": [
                [0, 15],
                [0, 15],
                [1, 16],
                [2, 17],
                [3, 18],
                [4, 19],
                [5, 20],
                [6, 21],
                [7, 22],
                [8, 23],
                [9, 24],
                [10, 25],
                [11, 26],
                [12, 27],
                [13, 28],
                [14, 29],
            ],
            "max_lossless_offer": 32768,
            "lossy_delivery": 57856,
            "headroom_penalty": 73014444032,
            "total_planned_bandwidth": 32768,
            "total_tunnel_cost": 32,
        },
    },
    {
        "id": "C-equal",
        "oracle": "t2",
        "input": {
            "prefer_headroom": True,
            "edges": [[0, 2, 3, 1], [1, 2, 8, 1]],
            "heads": [0, 1],
            "destinations": [2],
            "per_head": 2,
            "design_demand": 6,
        },
        "expected": {
            "max_bandwidth_per_tunnel": 1,
            "counts": [2, 2],
            "max_fractional_bandwidth_per_tunnel": "3/2",
            "max_lossless_offer": 6,
            "paths": [[0], [1], [0], [1]],
            "headroom_penalty": 5905580031,
            "total_planned_bandwidth": 4,
            "total_tunnel_cost": 4,
            "planned_bandwidth_per_tunnel": 1,
        },
    },
    {
        "id": "C-finite-cost",
        "oracle": "t2",
        "input": {
            "prefer_headroom": True,
            "edges": [[0, 1, 2, 1], [0, 1, 10, 2]],
            "heads": [0],
            "destinations": [1],
            "per_head": 2,
            "design_demand": 2,
        },
        "expected": {
            "max_bandwidth_per_tunnel": 5,
            "counts": [2, 0],
            "max_fractional_bandwidth_per_tunnel": "5",
            "max_lossless_offer": 2,
            "paths": [[0], [0]],
            "headroom_penalty": 6442450944,
            "total_planned_bandwidth": 2,
            "total_tunnel_cost": 2,
            "planned_bandwidth_per_tunnel": 1,
        },
    },
    {
        "id": "C-G10-0",
        "oracle": "t2",
        "input": {
            "prefer_headroom": True,
            "edges": [
                [0, 1, 40960, 1],
                [0, 5, 40960, 1],
                [0, 3, 24576, 2],
                [1, 0, 40960, 1],
                [1, 2, 40960, 1],
                [1, 4, 16384, 2],
                [2, 1, 40960, 1],
                [2, 3, 40960, 1],
                [2, 5, 20480, 2],
                [3, 2, 40960, 1],
                [3, 4, 40960, 1],
                [3, 0, 24576, 2],
                [4, 3, 40960, 1],
                [4, 5, 40960, 1],
                [4, 1, 16384, 2],
                [5, 0, 40960, 1],
                [5, 4, 40960, 1],
                [5, 2, 20480, 2],
            ],
            "heads": [0],
            "destinations": [3],
            "per_head": 16,
            "design_demand": 49152,
            "cap": [
                40960,
                40960,
                24576,
                40960,
                40960,
                16384,
                40960,
                40960,
                20480,
                40960,
                40960,
                24576,
                40960,
                40960,
                16384,
                40960,
                40960,
                20480,
            ],
        },
        "expected": {
            "max_bandwidth_per_tunnel": 6144,
            "counts": [4, 4, 8, 0, 4, 0, 0, 4, 0, 0, 0, 0, 4, 0, 0, 0, 4, 0],
            "max_fractional_bandwidth_per_tunnel": "6144",
            "max_lossless_offer": 49160,
            "paths": [
                [0, 4, 7],
                [0, 4, 7],
                [0, 4, 7],
                [0, 4, 7],
                [1, 16, 12],
                [1, 16, 12],
                [1, 16, 12],
                [1, 16, 12],
                [2],
                [2],
                [2],
                [2],
                [2],
                [2],
                [2],
                [2],
            ],
            "headroom_penalty": 38654705652,
            "total_planned_bandwidth": 49152,
            "total_tunnel_cost": 40,
            "planned_bandwidth_per_tunnel": 3072,
        },
    },
    {
        "id": "C-G10-1",
        "oracle": "t2",
        "input": {
            "prefer_headroom": True,
            "edges": [
                [0, 1, 40960, 1],
                [0, 5, 40960, 1],
                [0, 3, 24576, 2],
                [1, 0, 40960, 1],
                [1, 2, 40960, 1],
                [1, 4, 16384, 2],
                [2, 1, 40960, 1],
                [2, 3, 40960, 1],
                [2, 5, 20480, 2],
                [3, 2, 40960, 1],
                [3, 4, 40960, 1],
                [3, 0, 24576, 2],
                [4, 3, 40960, 1],
                [4, 5, 40960, 1],
                [4, 1, 16384, 2],
                [5, 0, 40960, 1],
                [5, 4, 40960, 1],
                [5, 2, 20480, 2],
            ],
            "heads": [1],
            "destinations": [4],
            "per_head": 16,
            "design_demand": 49152,
            "cap": [
                28672,
                28672,
                0,
                40960,
                28672,
                16384,
                40960,
                28672,
                20480,
                40960,
                40960,
                24576,
                28672,
                40960,
                16384,
                40960,
                28672,
                20480,
            ],
        },
        "expected": {
            "max_bandwidth_per_tunnel": 4096,
            "counts": [0, 6, 0, 6, 5, 5, 0, 5, 0, 0, 5, 0, 0, 0, 0, 0, 6, 0],
            "max_fractional_bandwidth_per_tunnel": "4096",
            "max_lossless_offer": 52431,
            "paths": [
                [3, 1, 16],
                [3, 1, 16],
                [3, 1, 16],
                [3, 1, 16],
                [3, 1, 16],
                [3, 1, 16],
                [4, 7, 10],
                [4, 7, 10],
                [4, 7, 10],
                [4, 7, 10],
                [4, 7, 10],
                [5],
                [5],
                [5],
                [5],
                [5],
            ],
            "total_planned_bandwidth": 49152,
            "total_tunnel_cost": 43,
            "planned_bandwidth_per_tunnel": 3072,
        },
    },
    {
        "id": "C-G10-2",
        "oracle": "t2",
        "input": {
            "prefer_headroom": True,
            "edges": [
                [0, 1, 40960, 1],
                [0, 5, 40960, 1],
                [0, 3, 24576, 2],
                [1, 0, 40960, 1],
                [1, 2, 40960, 1],
                [1, 4, 16384, 2],
                [2, 1, 40960, 1],
                [2, 3, 40960, 1],
                [2, 5, 20480, 2],
                [3, 2, 40960, 1],
                [3, 4, 40960, 1],
                [3, 0, 24576, 2],
                [4, 3, 40960, 1],
                [4, 5, 40960, 1],
                [4, 1, 16384, 2],
                [5, 0, 40960, 1],
                [5, 4, 40960, 1],
                [5, 2, 20480, 2],
            ],
            "heads": [2],
            "destinations": [5],
            "per_head": 16,
            "design_demand": 49152,
            "cap": [
                28672,
                10240,
                0,
                22528,
                13312,
                1024,
                40960,
                13312,
                20480,
                40960,
                25600,
                24576,
                28672,
                40960,
                16384,
                40960,
                10240,
                20480,
            ],
        },
        "expected": {
            "max_bandwidth_per_tunnel": 2560,
            "counts": [0, 4, 0, 4, 0, 0, 4, 4, 8, 0, 4, 0, 0, 4, 0, 0, 0, 0],
            "max_fractional_bandwidth_per_tunnel": "2560",
            "max_lossless_offer": 40960,
            "paths": [
                [6, 3, 1],
                [6, 3, 1],
                [6, 3, 1],
                [6, 3, 1],
                [7, 10, 13],
                [7, 10, 13],
                [7, 10, 13],
                [7, 10, 13],
                [8],
                [8],
                [8],
                [8],
                [8],
                [8],
                [8],
                [8],
            ],
            "total_planned_bandwidth": 40960,
            "total_tunnel_cost": 40,
            "planned_bandwidth_per_tunnel": 2560,
        },
    },
    {
        "id": "T2-dead-head",
        "oracle": "t2",
        "input": {
            "prefer_headroom": True,
            "edges": [[0, 2, 0, 1], [1, 2, 10, 1]],
            "heads": [0, 1],
            "destinations": [2],
            "per_head": 2,
            "design_demand": 10,
            "lossy_offer": 10,
        },
        "expected": {
            "entry_ids": [1, 3],
            "design_by_head": [[0, 5], [1, 5]],
            "heads_with_capacity": [1],
            "max_fractional_bandwidth_per_tunnel": "5",
            "max_bandwidth_per_tunnel": 5,
            "planned_bandwidth_per_tunnel": 3,
            "counts": [0, 2],
            "paths": [[1], [1]],
            "total_planned_bandwidth": 6,
            "total_tunnel_cost": 2,
            "max_lossless_offer": 0,
            "budget": 21,
            "headroom_penalty": 3865470565,
            "lossy_delivery": 5,
            "lossy_source_drops": 5,
        },
    },
    {
        "id": "T2-two-heads-equal",
        "oracle": "t2",
        "input": {
            "prefer_headroom": True,
            "edges": [[0, 2, 2, 1], [1, 2, 10, 1]],
            "heads": [0, 1],
            "destinations": [2],
            "per_head": 2,
        },
        "expected": {
            "entry_ids": [0, 1, 2, 3],
            "design_by_head": [],
            "max_fractional_bandwidth_per_tunnel": "1",
            "max_bandwidth_per_tunnel": 1,
            "planned_bandwidth_per_tunnel": 1,
            "counts": [2, 2],
            "paths": [[0], [1], [0], [1]],
            "total_planned_bandwidth": 4,
            "total_tunnel_cost": 4,
            "max_lossless_offer": 4,
            "headroom_penalty": 7730941132,
        },
    },
    {
        "id": "T2-idle",
        "oracle": "t2",
        "input": {
            "prefer_headroom": True,
            "edges": [[0, 1, 10, 1]],
            "heads": [0],
            "destinations": [1],
            "per_head": 4,
            "design_demand": 3,
            "lossy_offer": 3,
        },
        "expected": {
            "max_fractional_bandwidth_per_tunnel": "5/2",
            "max_bandwidth_per_tunnel": 2,
            "planned_bandwidth_per_tunnel": 1,
            "counts": [4],
            "total_planned_bandwidth": 4,
            "max_lossless_offer": 10,
            "headroom_penalty": 4294967294,
            "lossy_delivery": 3,
        },
    },
    {
        "id": "T2-zero-design",
        "oracle": "t2",
        "input": {
            "prefer_headroom": True,
            "edges": [[0, 1, 10, 1]],
            "heads": [0],
            "destinations": [1],
            "per_head": 2,
            "design_demand": 0,
            "lossy_offer": 4,
        },
        "expected": {
            "max_fractional_bandwidth_per_tunnel": "5",
            "max_bandwidth_per_tunnel": 5,
            "planned_bandwidth_per_tunnel": 1,
            "counts": [2],
            "total_planned_bandwidth": 2,
            "max_lossless_offer": 10,
            "lossy_delivery": 4,
        },
    },
    {
        "id": "T2-remainder-first",
        "oracle": "t2",
        "input": {
            "prefer_headroom": True,
            "edges": [[0, 1, 3, 1], [0, 1, 2, 1]],
            "heads": [0],
            "destinations": [1],
            "per_head": 2,
            "design_demand": 5,
        },
        "expected": {
            "max_fractional_bandwidth_per_tunnel": "2",
            "max_bandwidth_per_tunnel": 2,
            "planned_bandwidth_per_tunnel": 2,
            "counts": [1, 1],
            "paths": [[0], [1]],
            "total_planned_bandwidth": 4,
            "max_lossless_offer": 5,
        },
    },
    {
        "id": "T2-remainder-second",
        "oracle": "t2",
        "input": {
            "prefer_headroom": True,
            "edges": [[0, 1, 2, 1], [0, 1, 3, 1]],
            "heads": [0],
            "destinations": [1],
            "per_head": 2,
            "design_demand": 5,
        },
        "expected": {
            "max_fractional_bandwidth_per_tunnel": "2",
            "max_bandwidth_per_tunnel": 2,
            "planned_bandwidth_per_tunnel": 2,
            "counts": [1, 1],
            "paths": [[0], [1]],
            "total_planned_bandwidth": 4,
            "max_lossless_offer": 4,
        },
    },
    {
        "id": "C-G1-UP_TO",
        "oracle": "upto",
        "input": {
            "prefer_headroom": True,
            "edges": [
                [0, 1, 8192, 1],
                [0, 2, 40960, 1],
                [1, 3, 8192, 1],
                [2, 3, 40960, 1],
            ],
            "heads": [0],
            "destinations": [3],
            "limit": 6,
            "design_demand": 49152,
        },
        "expected": {
            "bounds": [40960, 40960, 40960, 40960, 40960, 49152],
            "scores": [40960, 40960, 40960, 40960, 40960, 49152],
            "tunnels_per_head": 6,
            "delivered": 49152,
        },
    },
    {
        "id": "T2-UP_TO-common",
        "oracle": "upto",
        "input": {
            "prefer_headroom": True,
            "edges": [[0, 3, 6, 1], [1, 3, 3, 1], [1, 3, 3, 1]],
            "heads": [0, 1],
            "destinations": [3],
            "limit": 2,
            "design_demand": 12,
        },
        "expected": {
            "bounds": [7, 12],
            "scores": [7, 12],
            "tunnels_per_head": 2,
            "delivered": 12,
            "paths": [[0], [1], [0], [2]],
        },
    },
    {
        "id": "T2-UP_TO-dead-head",
        "oracle": "upto",
        "input": {
            "prefer_headroom": True,
            "edges": [[0, 2, 0, 1], [1, 2, 10, 1]],
            "heads": [0, 1],
            "destinations": [2],
            "limit": 2,
            "design_demand": 10,
        },
        "expected": {
            "bounds": [21, 21],
            "scores": [5, 5],
            "tunnels_per_head": 1,
            "delivered": 0,
            "paths": [[1]],
        },
    },
    {
        "id": "T2-tiny-head",
        "oracle": "t2",
        "input": {
            "prefer_headroom": True,
            "edges": [[0, 2, 1, 1], [1, 2, 100, 1]],
            "heads": [0, 1],
            "destinations": [2],
            "per_head": 16,
        },
        "expected": {
            "heads_with_capacity": [0, 1],
            "max_fractional_bandwidth_per_tunnel": "1/16",
            "max_bandwidth_per_tunnel": 0,
            "counts": [0, 0],
            "paths": [],
            "max_lossless_offer": 0,
            "insufficient_resolution": True,
        },
    },
    {
        "id": "T2-tiny-head-UP_TO",
        "oracle": "upto",
        "input": {
            "prefer_headroom": True,
            "edges": [[0, 2, 1, 1], [1, 2, 100, 1]],
            "heads": [0, 1],
            "destinations": [2],
            "limit": 3,
        },
        "expected": {
            "bounds": [2, 0, 0],
            "scores": [2, 0, 0],
            "tunnels_per_head": 1,
            "delivered": 2,
            "paths": [[0], [1]],
        },
    },
    {
        "id": "E-duplicates",
        "oracle": "bound",
        "input": {"cap": [4, 8], "routes": [[0], [1], [1]], "weights": [1, 1, 1]},
        "expected": 12,
    },
    {
        "id": "E-deduplicated",
        "oracle": "bound",
        "input": {"cap": [4, 8], "routes": [[0], [1]], "weights": [1, 1]},
        "expected": 8,
    },
    {
        "id": "E-split-contexts",
        "oracle": "split_contexts",
        "input": {
            "cap": [1, 1],
            "routes": [[0, 1], [0, 1]],
            "weights": [1, 1],
            "pool": False,
        },
        "expected": {"bound": 1, "loads_at_2": [2, 0]},
    },
    {
        "id": "E-split-contexts-pooled",
        "oracle": "split_contexts",
        "input": {
            "cap": [1, 1],
            "routes": [[0, 1], [0, 1]],
            "weights": [1, 1],
            "pool": True,
        },
        "expected": {"bound": 2, "loads_at_2": [1, 1]},
    },
    {
        "id": "E-static-path-fill",
        "oracle": "bound",
        "input": {"cap": [6, 4], "routes": [[0], [1]], "weights": [1, 1]},
        "expected": 9,
    },
    {
        "id": "E-zero-first",
        "oracle": "bound",
        "input": {"cap": [0, 4], "routes": [[0], [1]], "weights": [1, 1]},
        "expected": 0,
    },
    {
        "id": "E-zero-last",
        "oracle": "bound",
        "input": {"cap": [4, 0], "routes": [[0], [1]], "weights": [1, 1]},
        "expected": 1,
    },
    {
        "id": "E-support-pinned-bound",
        "oracle": "bound",
        "input": {
            "cap": [10, 2, 10, 10],
            "routes": [[0, 1], [2, 3]],
            "weights": [1, 1],
        },
        "expected": 4,
    },
    {
        "id": "STALE-remote-tail",
        "oracle": "stale_dag",
        "input": {
            "edges": [[0, 1, 100, 1], [0, 2, 100, 1], [1, 3, 100, 1], [2, 3, 100, 1]],
            "sources": [0],
            "destinations": 3,
            "offer": 100,
            "weights": [1, 1, 1, 1],
            "black_links": [3],
        },
        "expected": {
            "delivered": 50,
            "edge_flows": [50, 50, 50, 0],
            "edge_drops": [0, 0, 0, 50],
        },
    },
    {
        "id": "STALE-remote-tail-repair",
        "oracle": "stale_dag",
        "input": {
            "edges": [[0, 1, 100, 1], [0, 2, 100, 1], [1, 3, 100, 1], [2, 3, 100, 1]],
            "sources": [0],
            "destinations": 3,
            "offer": 100,
            "weights": [1, 1, 1, 1],
            "mask_links": [3],
        },
        "expected": {
            "delivered": 100,
            "edge_flows": [100, 0, 100, 0],
            "edge_drops": [0, 0, 0, 0],
        },
    },
    {
        "id": "STALE-wcmp-member",
        "oracle": "stale_dag",
        "input": {
            "edges": [[0, 1, 10, 1], [0, 1, 40, 1]],
            "sources": [0],
            "destinations": 1,
            "offer": 40,
            "weights": [10, 40],
            "black_links": [0],
        },
        "expected": {"delivered": 32, "edge_flows": [0, 32], "edge_drops": [8, 0]},
    },
    {
        "id": "STALE-wcmp-zero-capacity",
        "oracle": "stale_dag",
        "input": {
            "edges": [[0, 1, 0, 1], [0, 1, 40, 1]],
            "sources": [0],
            "destinations": 1,
            "offer": 40,
            "weights": [0, 40],
        },
        "expected": {"delivered": 40, "edge_flows": [0, 40], "edge_drops": [0, 0]},
    },
    {
        "id": "STALE-relief",
        "oracle": "stale_paths",
        "input": {
            "edges": [[0, 2, 10, 1], [1, 2, 10, 1], [2, 3, 10, 1]],
            "tunnels": [[0, [0, 2]], [1, [1, 2]]],
            "offers": {0: 10, 1: 10},
            "black_links": [0],
        },
        "expected": {
            "status": ["BLACKHOLED", "UP"],
            "delivered": [0, 10],
            "edge_drops": [10, 0, 0],
        },
    },
    {
        "id": "STALE-relief-original",
        "oracle": "stale_paths",
        "input": {
            "edges": [[0, 2, 10, 1], [1, 2, 10, 1], [2, 3, 10, 1]],
            "tunnels": [[0, [0, 2]], [1, [1, 2]]],
            "offers": {0: 10, 1: 10},
        },
        "expected": {
            "status": ["UP", "UP"],
            "delivered": [5, 5],
            "edge_drops": [0, 0, 10],
        },
    },
    {
        "id": "STALE-tunnel-node",
        "oracle": "stale_paths",
        "input": {
            "edges": [[0, 1, 10, 1], [0, 2, 10, 1], [1, 3, 10, 1], [2, 3, 10, 1]],
            "tunnels": [[0, [0, 2]], [0, [1, 3]]],
            "offers": {0: 10},
            "black_nodes": [2],
        },
        "expected": {
            "status": ["UP", "BLACKHOLED"],
            "assigned": [5, 5],
            "delivered": [5, 0],
            "edge_drops": [0, 5, 0, 0],
            "lossless_admitted": 1,
        },
    },
    {
        "id": "STALE-tunnel-node-repair",
        "oracle": "stale_paths",
        "input": {
            "edges": [[0, 1, 10, 1], [0, 2, 10, 1], [1, 3, 10, 1], [2, 3, 10, 1]],
            "tunnels": [[0, [0, 2]], [0, [1, 3]]],
            "offers": {0: 10},
            "mask_nodes": [2],
        },
        "expected": {
            "status": ["UP", "FAILED"],
            "assigned": [10, 0],
            "delivered": [10, 0],
            "edge_drops": [0, 0, 0, 0],
            "lossless_admitted": 10,
        },
    },
    {
        "id": "STALE-pooled-original",
        "oracle": "stale_dag",
        "input": {
            "edges": [[0, 2, 10, 1], [1, 2, 10, 1], [2, 3, 10, 1], [2, 3, 10, 1]],
            "sources": [0, 1],
            "destinations": 3,
            "offer": 2,
            "weights": [1, 1, 1, 1],
        },
        "expected": {
            "delivered": 2,
            "edge_flows": [1, 1, 1, 1],
            "edge_drops": [0, 0, 0, 0],
            "source_drops": 0,
        },
    },
    {
        "id": "STALE-pooled-remainder",
        "oracle": "stale_dag",
        "input": {
            "edges": [[0, 2, 10, 1], [1, 2, 10, 1], [2, 3, 10, 1], [2, 3, 10, 1]],
            "sources": [0, 1],
            "destinations": 3,
            "offer": 2,
            "weights": [1, 1, 1, 1],
            "black_links": [0, 3],
        },
        "expected": {
            "delivered": 1,
            "edge_flows": [0, 1, 1, 0],
            "edge_drops": [1, 0, 0, 0],
            "source_drops": 0,
        },
    },
    {
        "id": "STALE-all-zero-weight",
        "oracle": "stale_dag",
        "input": {
            "edges": [[0, 1, 0, 1]],
            "sources": [0],
            "destinations": 1,
            "offer": 2,
            "weights": [0],
        },
        "expected": {
            "delivered": 0,
            "edge_flows": [0],
            "edge_drops": [0],
            "source_drops": 2,
        },
    },
    {
        "id": "STALE-black-source",
        "oracle": "stale_dag",
        "input": {
            "edges": [[0, 1, 10, 1]],
            "sources": [0],
            "destinations": 1,
            "offer": 10,
            "weights": [1],
            "black_nodes": [0],
        },
        "expected": {
            "delivered": 0,
            "source_drops": 10,
            "edge_flows": [0],
            "edge_drops": [0],
        },
    },
    {
        "id": "STALE-black-source-tunnel",
        "oracle": "stale_paths",
        "input": {
            "edges": [[0, 1, 10, 1]],
            "tunnels": [[0, [0]]],
            "offers": {0: 10},
            "black_nodes": [0],
        },
        "expected": {
            "status": ["BLACKHOLED"],
            "assigned": [0],
            "delivered": [0],
            "source_drops": 10,
            "edge_drops": [0],
        },
    },
    {
        "id": "STALE-mask-and-black-hole",
        "oracle": "stale_paths",
        "input": {
            "edges": [
                [0, 1, 10, 1],
                [0, 2, 10, 1],
                [0, 3, 10, 1],
                [1, 4, 10, 1],
                [2, 4, 10, 1],
                [3, 4, 10, 1],
            ],
            "tunnels": [[0, [0, 3]], [0, [1, 4]], [0, [2, 5]]],
            "offers": {0: 10},
            "mask_links": [3],
            "black_links": [3, 4],
        },
        "expected": {
            "status": ["FAILED", "BLACKHOLED", "UP"],
            "assigned": [0, 5, 5],
            "delivered": [0, 0, 5],
            "edge_drops": [0, 0, 0, 0, 5, 0],
            "source_drops": 0,
        },
    },
    {
        "id": "D-merge",
        "oracle": "lossy",
        "input": {"offers": [8, 4], "routes": [[0], [0]], "cap": [9]},
        "expected": {
            "carried": [[6], [3]],
            "edge_drops": [[2], [1]],
            "delivered": [6, 3],
            "edge_flows": [9],
        },
    },
    {
        "id": "D-key-tie",
        "oracle": "lossy",
        "input": {"offers": [1, 1], "routes": [[0], [0]], "cap": [1]},
        "expected": {
            "carried": [[1], [0]],
            "edge_drops": [[0], [1]],
            "delivered": [1, 0],
            "edge_flows": [1],
        },
    },
    {
        "id": "D-rounding",
        "oracle": "lossy",
        "input": {"offers": [5, 1, 1], "routes": [[0], [0], [0]], "cap": [3]},
        "expected": {
            "carried": [[3], [0], [0]],
            "edge_drops": [[2], [1], [1]],
            "delivered": [3, 0, 0],
            "edge_flows": [3],
        },
    },
    {
        "id": "D-hierarchy",
        "oracle": "lossy",
        "input": {
            "offers": [1, 1, 2],
            "routes": [[0], [0], [0]],
            "cap": [1],
            "labels": [[0, 0, 0], [0, 0, 1], [1, 0, 0]],
        },
        "expected": {
            "carried": [[1], [0], [0]],
            "edge_drops": [[0], [1], [2]],
            "delivered": [1, 0, 0],
            "edge_flows": [1],
        },
    },
    {
        "id": "D-cycle-compressed",
        "oracle": "lossy",
        "input": {"offers": [1, 1], "routes": [[0, 1], [1, 0]], "cap": [2, 2]},
        "expected": {
            "carried": [[1, 1], [1, 1]],
            "edge_drops": [[0, 0], [0, 0]],
            "delivered": [1, 1],
            "edge_flows": [2, 2],
        },
    },
    {
        "id": "D-zero-service",
        "oracle": "lossy",
        "input": {"offers": [2], "routes": [[0]], "cap": [0]},
        "expected": {
            "carried": [[0]],
            "edge_drops": [[2]],
            "delivered": [0],
            "edge_flows": [0],
        },
    },
    {
        "id": "D-subdivision-before",
        "oracle": "lossy",
        "input": {
            "offers": [1, 1, 1],
            "routes": [[0, 1], [0, 2], [0, 3]],
            "cap": [1, 1, 1, 1],
            "labels": [[0, 0, 0], [0, 0, 1], [0, 0, 2]],
        },
        "expected": {
            "carried": [[1, 1], [0, 0], [0, 0]],
            "delivered": [1, 0, 0],
            "edge_flows": [1, 1, 0, 0],
        },
    },
    {
        "id": "D-subdivision-after",
        "oracle": "lossy",
        "input": {
            "offers": [1, 1, 1],
            "routes": [[0, 1, 2], [0, 3], [0, 4]],
            "cap": [1, 1, 1, 1, 1],
            "labels": [[0, 0, 0], [0, 0, 1], [0, 0, 2]],
        },
        "expected": {
            "carried": [[1, 1, 1], [0, 0], [0, 0]],
            "delivered": [1, 0, 0],
            "edge_flows": [1, 1, 1, 0, 0],
        },
    },
    {
        "id": "E-support-pinned-loss",
        "oracle": "lossy",
        "input": {
            "offers": [6, 6],
            "routes": [[0, 1], [2, 3]],
            "cap": [10, 2, 10, 10],
            "labels": [[0, 0, 0], [0, 0, 1]],
        },
        "expected": {
            "carried": [[6, 2], [6, 6]],
            "edge_drops": [[0, 4], [0, 0]],
            "delivered": [2, 6],
            "edge_flows": [6, 2, 6, 6],
        },
    },
    {
        "id": "D-C2",
        "oracle": "lossy",
        "input": {"offers": [1, 2], "routes": [[0, 1], [1, 0]], "cap": [1, 2]},
        "expected": {
            "error": "lossy forwarding iteration entered a period-2 cycle",
            "rounds": 3,
        },
    },
    {
        "id": "E-ledger",
        "oracle": "ledger",
        "input": {
            "service": 10,
            "reservation": 10,
            "carried": 4,
            "strict": 6,
            "new_delivery": 2,
            "old_delivery": 4,
        },
        "expected": {
            "carried": 4,
            "residual": 6,
            "reserved": 6,
            "unreserved": 4,
            "delivery_delta": -2,
        },
    },
    {
        "id": "RC02-increment",
        "oracle": "epoch",
        "input": {"cap": [1, 1], "calls": [1, 1]},
        "expected": {
            "requested": 2,
            "delivered": 1,
            "edge_drops": 1,
            "edge_flows": [1, 0],
            "residual": [0, 1],
            "per_call_delivery": [1, 0],
        },
    },
    {
        "id": "RC02-increment-lossless",
        "oracle": "epoch",
        "input": {"cap": [1, 1], "calls": [1, 1], "admission": "LOSSLESS"},
        "expected": {
            "requested": 2,
            "delivered": 2,
            "edge_drops": 0,
            "edge_flows": [1, 1],
            "residual": [0, 0],
            "per_call_delivery": [1, 1],
        },
    },
    {
        "id": "RC02-full",
        "oracle": "epoch",
        "input": {"cap": [1, 1], "calls": [1, 1], "full": True},
        "expected": {
            "requested": 2,
            "delivered": 2,
            "edge_drops": 0,
            "edge_flows": [1, 1],
            "residual": [0, 0],
            "per_call_delivery": [2],
        },
    },
    {
        "id": "INC-lossless",
        "oracle": "increments",
        "input": {"caps": [1, 3], "calls": [3, 1], "admission": "LOSSLESS"},
        "expected": {
            "call_delivered": [2, 0],
            "requested": 4,
            "admitted": 2,
            "unadmitted": 2,
            "offered": None,
            "delivered": 2,
        },
    },
    {
        "id": "INC-lossy",
        "oracle": "increments",
        "input": {"caps": [1, 3], "calls": [3, 1], "admission": "LOSSY"},
        "expected": {
            "call_delivered": [2, 1],
            "requested": 4,
            "offered": 4,
            "edge_drops": [1, 0],
            "unadmitted": 0,
            "delivered": 3,
        },
    },
    {
        "id": "INC-lossless-refusal",
        "oracle": "increments",
        "input": {
            "caps": [10],
            "calls": [5, 1],
            "admission": "LOSSLESS",
            "outside": [[8], [5]],
        },
        "expected": {
            "call_delivered": [2, 1],
            "requested": 6,
            "admitted": 3,
            "unadmitted": 3,
            "delivered": 3,
        },
    },
    {
        "id": "R4-refusal",
        "oracle": "refusal",
        "input": {"capacity": 10, "outside": [8, 5], "calls": [5, 1]},
        "expected": {"correct": [2, 3], "requeue": [2, 5]},
    },
    {
        "id": "RC02-zero-unprepared",
        "oracle": "zero_state",
        "input": {"service": [5, 8], "outside": [1, 2]},
        "expected": {
            "call_delivered": 0,
            "call_undelivered": 0,
            "delivery_delta": 0,
            "call_edge_drops": [],
            "call_source_drops": 0,
            "termination": "COMPLETE",
            "routes_changed": False,
            "weights_changed": False,
            "reservations_changed": False,
            "program": None,
            "edge_flows": [0, 0],
            "edge_reservations": [0, 0],
            "residual": [4, 6],
            "min_cut": None,
            "marginal_cost_curve": None,
        },
    },
    {
        "id": "FLAGS-controller-increment",
        "oracle": "change_flags",
        "input": {
            "path_cap": 20,
            "calls": [["place", 5], ["place", 5], ["replace", 10]],
            "recomputed": True,
        },
        "expected": {
            "routes_changed": [False, False],
            "weights_changed": [True, False],
            "reservations_changed": [False, False],
        },
    },
    {
        "id": "FLAGS-fixed-offer",
        "oracle": "change_flags",
        "input": {
            "path_cap": 20,
            "calls": [["replace", 5], ["replace", 10], ["replace", 10]],
            "recomputed": False,
        },
        "expected": {
            "routes_changed": [False, False],
            "weights_changed": [False, False],
            "reservations_changed": [False, False],
        },
    },
    {
        "id": "ERR-ownership-texts",
        "oracle": "ownership_texts",
        "input": {
            "objects": [
                "program",
                "support",
                "black_holes",
                "Eligibility",
                "ShortestPathCache",
                "dag",
            ]
        },
        "expected": [
            "program belongs to another graph",
            "support belongs to another graph",
            "black_holes belongs to another graph",
            "foreign Eligibility: built for another graph",
            "foreign ShortestPathCache: bound to another graph",
            "dag belongs to another graph",
        ],
    },
    {
        "id": "T2-prefer-headroom-off",
        "oracle": "t2",
        "input": {
            "edges": [[0, 1, 10, 1], [0, 1, 10, 1]],
            "heads": [0],
            "destinations": [1],
            "per_head": 2,
            "design_demand": 2,
            "prefer_headroom": False,
            "lossy_offer": 12,
        },
        "expected": {
            "max_bandwidth_per_tunnel": 10,
            "planned_bandwidth_per_tunnel": 1,
            "counts": [2, 0],
            "max_lossless_offer": 10,
            "lossy_delivery": 10,
        },
    },
    {
        "id": "T2-prefer-headroom-on",
        "oracle": "t2",
        "input": {
            "edges": [[0, 1, 10, 1], [0, 1, 10, 1]],
            "heads": [0],
            "destinations": [1],
            "per_head": 2,
            "design_demand": 2,
            "prefer_headroom": True,
            "lossy_offer": 12,
        },
        "expected": {
            "max_bandwidth_per_tunnel": 10,
            "planned_bandwidth_per_tunnel": 1,
            "counts": [1, 1],
            "headroom_penalty": 858993458,
            "max_lossless_offer": 20,
            "lossy_delivery": 12,
        },
    },
    {
        "id": "T2-support",
        "oracle": "t2",
        "input": {
            "edges": [[0, 1, 6, 1], [0, 1, 6, 1], [0, 1, 100, 2]],
            "heads": [0],
            "destinations": [1],
            "per_head": 2,
            "support": [0, 1],
        },
        "expected": {
            "max_bandwidth_per_tunnel": 6,
            "counts": [1, 1, 0],
            "max_lossless_offer": 12,
        },
    },
    {
        "id": "T2-support-masked",
        "oracle": "t2",
        "input": {
            "edges": [[0, 1, 6, 1], [0, 1, 6, 1], [0, 1, 100, 2]],
            "heads": [0],
            "destinations": [1],
            "per_head": 2,
            "support": [0, 1],
            "excluded": [1],
        },
        "expected": {
            "max_bandwidth_per_tunnel": 3,
            "counts": [2, 0, 0],
            "paths": [[0], [0]],
            "entry_ids": [0, 1],
            "max_lossless_offer": 6,
        },
    },
    {
        "id": "T2-support-mask-outside",
        "oracle": "t2",
        "input": {
            "edges": [[0, 1, 6, 1], [0, 1, 6, 1], [0, 1, 100, 2]],
            "heads": [0],
            "destinations": [1],
            "per_head": 2,
            "support": [0, 1],
            "excluded": [2],
        },
        "expected": {
            "max_bandwidth_per_tunnel": 6,
            "counts": [1, 1, 0],
            "paths": [[0], [1]],
            "entry_ids": [0, 1],
            "max_lossless_offer": 12,
        },
    },
    {
        "id": "CFG-defaults",
        "oracle": "model_config",
        "input": {
            "cases": [
                {"model": "ShortestPaths", "kwargs": {}},
                {"model": "ShortestPaths", "kwargs": {"allocation": "EQUAL"}},
                {"model": "ShortestPaths", "kwargs": {"allocation": "WEIGHTED"}},
                {"model": "ShortestPaths", "kwargs": {"allocation": "OPTIMIZED"}},
                {
                    "model": "CentralizedTE",
                    "kwargs": {"source_allocation": "EQUAL", "tunnel_count": "16"},
                },
            ]
        },
        "expected": [
            "allocation=EQUAL source_allocation=EQUAL weight_basis=None tie_break=FEWER_HOPS multipath=True parallel_edges=True",
            "allocation=EQUAL source_allocation=EQUAL weight_basis=None tie_break=FEWER_HOPS multipath=True parallel_edges=True",
            "allocation=WEIGHTED source_allocation=EQUAL weight_basis=CAPACITY tie_break=FEWER_HOPS multipath=True parallel_edges=True",
            "allocation=OPTIMIZED source_allocation=OPTIMIZED weight_basis=None tie_break=FEWER_HOPS multipath=True parallel_edges=True",
            "source_allocation=EQUAL support=None tunnel_count=16 prefer_headroom=False",
        ],
    },
    {
        "id": "CFG-construction",
        "oracle": "model_config",
        "input": {
            "cases": [
                {"model": "CentralizedTE", "kwargs": {}},
                {"model": "CentralizedTE", "kwargs": {"source_allocation": None}},
                {"model": "CentralizedTE", "kwargs": {"tunnel_count": "16"}},
                {
                    "model": "CentralizedTE",
                    "kwargs": {"source_allocation": "EQUAL", "selection": "default"},
                },
                {
                    "model": "DistributedTE",
                    "kwargs": {"tunnel_count": "4", "source_allocation": "EQUAL"},
                },
                {
                    "model": "DistributedTE",
                    "kwargs": {"tunnel_count": "4", "multipath": False},
                },
                {"model": "Pinned", "kwargs": {"program": "p", "allocation": "EQUAL"}},
                {
                    "model": "FlowPolicyConfig",
                    "kwargs": {"routing": "ShortestPaths()", "allocation": "EQUAL"},
                },
                {"model": "ShortestPaths", "kwargs": {"multipath": None}},
                {"model": "ShortestPaths", "kwargs": {"placement": "EQUAL"}},
            ]
        },
        "expected": [
            "TypeError",
            "TypeError",
            "TypeError",
            "TypeError",
            "TypeError",
            "TypeError",
            "TypeError",
            "TypeError",
            "TypeError",
            "TypeError",
        ],
    },
    {
        "id": "CFG-matrix",
        "oracle": "model_config",
        "input": {
            "cases": [
                {
                    "model": "ShortestPaths",
                    "kwargs": {"source_allocation": "OPTIMIZED"},
                },
                {
                    "model": "ShortestPaths",
                    "kwargs": {"allocation": "OPTIMIZED"},
                    "admission": "LOSSY",
                },
                {"model": "ShortestPaths", "kwargs": {"tie_break": "HIGHER_CAPACITY"}},
                {
                    "model": "ShortestPaths",
                    "kwargs": {"tie_break": "HIGHER_CAPACITY", "multipath": False},
                    "sources": 2,
                },
                {
                    "model": "ShortestPaths",
                    "kwargs": {"tie_break": "HIGHER_CAPACITY", "multipath": False},
                },
                {
                    "model": "ShortestPaths",
                    "kwargs": {"tie_break": "HIGHER_HEADROOM", "parallel_edges": False},
                    "sources": 2,
                },
                {"model": "ShortestPaths", "kwargs": {"weight_basis": "CAPACITY"}},
                {
                    "model": "ShortestPaths",
                    "kwargs": {"allocation": "WEIGHTED", "weight_basis": "EXPLICIT"},
                },
                {
                    "model": "CentralizedTE",
                    "kwargs": {"source_allocation": "OPTIMIZED", "tunnel_count": "2"},
                    "sources": 2,
                },
                {
                    "model": "CentralizedTE",
                    "kwargs": {"source_allocation": "OPTIMIZED", "tunnel_count": "2"},
                },
                {
                    "model": "CentralizedTE",
                    "kwargs": {"source_allocation": "EQUAL"},
                    "admission": "LOSSY",
                },
                {
                    "model": "CentralizedTE",
                    "kwargs": {"source_allocation": "EQUAL", "prefer_headroom": True},
                },
                {
                    "model": "CentralizedTE",
                    "kwargs": {
                        "source_allocation": "EQUAL",
                        "support": "dag",
                        "tunnel_count": "2",
                    },
                    "admission": "LOSSY",
                },
                {"model": "DistributedTE", "kwargs": {"tunnel_count": "4 up_to"}},
            ]
        },
        "expected": [
            "ShortestPaths with EQUAL or WEIGHTED requires source_allocation=EQUAL",
            "ShortestPaths OPTIMIZED requires admission=LOSSLESS",
            "ShortestPaths tie_break requires multipath=False or parallel_edges=False",
            "ShortestPaths with multiple input sources requires tie_break=FEWER_HOPS when multipath=False",
            "allocation=EQUAL source_allocation=EQUAL weight_basis=None tie_break=HIGHER_CAPACITY multipath=False parallel_edges=True",
            "allocation=EQUAL source_allocation=EQUAL weight_basis=None tie_break=HIGHER_HEADROOM multipath=True parallel_edges=False",
            "weight_basis requires allocation=WEIGHTED",
            "weight_basis=EXPLICIT requires edge_weights",
            "CentralizedTE with tunnel_count and several sources requires source_allocation=EQUAL",
            "source_allocation=OPTIMIZED support=None tunnel_count=2 prefer_headroom=False",
            "CentralizedTE without tunnel_count requires admission=LOSSLESS",
            "prefer_headroom requires tunnel_count",
            "source_allocation=EQUAL support=dag tunnel_count=2 prefer_headroom=False",
            "DistributedTE requires an exact tunnel count",
        ],
    },
    {
        "id": "CFG-first-error",
        "oracle": "model_config",
        "input": {
            "cases": [
                {
                    "model": "CentralizedTE",
                    "kwargs": {"source_allocation": "EQUAL", "prefer_headroom": True},
                    "admission": "LOSSY",
                },
                {
                    "model": "ShortestPaths",
                    "kwargs": {
                        "source_allocation": "OPTIMIZED",
                        "weight_basis": "CAPACITY",
                    },
                },
                {
                    "model": "ShortestPaths",
                    "kwargs": {
                        "allocation": "OPTIMIZED",
                        "tie_break": "HIGHER_CAPACITY",
                    },
                    "admission": "LOSSY",
                },
                {
                    "model": "ShortestPaths",
                    "kwargs": {
                        "tie_break": "HIGHER_CAPACITY",
                        "weight_basis": "CAPACITY",
                    },
                },
                {
                    "model": "ShortestPaths",
                    "kwargs": {
                        "tie_break": "HIGHER_CAPACITY",
                        "multipath": False,
                        "weight_basis": "CAPACITY",
                    },
                    "sources": 2,
                },
                {
                    "model": "CentralizedTE",
                    "kwargs": {
                        "source_allocation": "OPTIMIZED",
                        "tunnel_count": "0",
                        "prefer_headroom": True,
                    },
                    "sources": 2,
                },
                {
                    "model": "DistributedTE",
                    "kwargs": {
                        "tunnel_count": "0 up_to",
                        "cost_ceiling": ["fraction", 1, 2],
                    },
                },
                {
                    "model": "DistributedTE",
                    "kwargs": {
                        "tunnel_count": "4 up_to",
                        "cost_ceiling": ["fraction", 1, 2],
                    },
                },
                {
                    "model": "DistributedTE",
                    "kwargs": {
                        "tunnel_count": "4",
                        "cost_ceiling": ["fraction", 1, 2],
                        "min_unreserved_bandwidth": -1,
                    },
                },
                {
                    "model": "DistributedTE",
                    "kwargs": {"tunnel_count": "4", "min_unreserved_bandwidth": -1},
                },
            ]
        },
        "expected": [
            "CentralizedTE without tunnel_count requires admission=LOSSLESS",
            "ShortestPaths with EQUAL or WEIGHTED requires source_allocation=EQUAL",
            "ShortestPaths OPTIMIZED requires admission=LOSSLESS",
            "ShortestPaths tie_break requires multipath=False or parallel_edges=False",
            "ShortestPaths with multiple input sources requires tie_break=FEWER_HOPS when multipath=False",
            "CentralizedTE with tunnel_count and several sources requires source_allocation=EQUAL",
            "tunnel count must be in [1, 10000]",
            "DistributedTE requires an exact tunnel count",
            "cost_ceiling CostStretch factor must be >= 1 (got 1/2)",
            "min_unreserved_bandwidth must be an integer in [0, 2^62)",
        ],
    },
    {
        "id": "EVAL-first-error",
        "oracle": "evaluation_validation",
        "input": {
            "cases": [
                {
                    "model": "CentralizedTE",
                    "kwargs": {"source_allocation": "EQUAL", "prefer_headroom": True},
                    "offer": -1,
                },
                {
                    "model": "DistributedTE",
                    "kwargs": {"tunnel_count": "4"},
                    "admission": "LOSSY",
                    "offer": -1,
                    "design_demand": math.inf,
                },
                {
                    "model": "DistributedTE",
                    "kwargs": {"tunnel_count": "4"},
                    "admission": "LOSSY",
                    "offer": 5,
                    "design_demand": math.inf,
                },
                {
                    "model": "DistributedTE",
                    "kwargs": {"tunnel_count": "4"},
                    "admission": "LOSSY",
                },
                {
                    "model": "DistributedTE",
                    "kwargs": {"tunnel_count": "4"},
                    "admission": "LOSSY",
                    "design_demand": 10,
                },
                {
                    "model": "ShortestPaths",
                    "kwargs": {},
                    "offer": math.inf,
                    "design_demand": 5,
                },
                {"model": "ShortestPaths", "kwargs": {}, "design_demand": 5},
                {
                    "model": "CentralizedTE",
                    "kwargs": {"source_allocation": "EQUAL", "tunnel_count": "16"},
                    "admission": "LOSSY",
                },
                {
                    "model": "CentralizedTE",
                    "kwargs": {"source_allocation": "EQUAL", "tunnel_count": "16"},
                    "admission": "LOSSY",
                    "design_demand": 10,
                },
                {"model": "Pinned", "kwargs": {"program": "p"}, "admission": "LOSSY"},
                {"model": "ShortestPaths", "kwargs": {}, "admission": "LOSSY"},
                {
                    "model": "ShortestPaths",
                    "kwargs": {},
                    "admission": "LOSSY",
                    "offer": 4,
                },
            ]
        },
        "expected": [
            "prefer_headroom requires tunnel_count",
            "offer must be an integer in [0, 2^62)",
            "DistributedTE requires a finite design demand",
            "DistributedTE evaluation requires design_demand",
            "DistributedTE LOSSY evaluation requires a finite offer",
            "offer must be an integer in [0, 2^62)",
            "ShortestPaths does not accept design_demand",
            "CentralizedTE LOSSY evaluation requires a finite offer",
            "offer=10 design_demand=10",
            "Pinned LOSSY evaluation requires a finite offer",
            "offer=first-hop experiment design_demand=None",
            "offer=4 design_demand=None",
        ],
    },
    {
        "id": "CFG-cache-order",
        "oracle": "model_config",
        "input": {
            "cases": [
                {
                    "model": "DistributedTE",
                    "kwargs": {"tunnel_count": "0"},
                    "cache": True,
                },
                {
                    "model": "Pinned",
                    "kwargs": {"program": {"sources": [0], "destinations": [2]}},
                    "cache": True,
                    "endpoints": [[0], [1]],
                },
                {
                    "model": "CentralizedTE",
                    "kwargs": {"source_allocation": "EQUAL", "tunnel_count": "0"},
                    "cache": True,
                },
                {
                    "model": "CentralizedTE",
                    "kwargs": {"source_allocation": "EQUAL", "prefer_headroom": True},
                    "cache": True,
                },
                {
                    "model": "ShortestPaths",
                    "kwargs": {"allocation": "WEIGHTED"},
                    "cache": True,
                },
                {"model": "DistributedTE", "kwargs": {"tunnel_count": "0"}},
                {
                    "model": "Pinned",
                    "kwargs": {"program": {"sources": [0], "destinations": [2]}},
                    "endpoints": [[0], [1]],
                },
                {
                    "model": "Pinned",
                    "kwargs": {"program": {"sources": [0], "destinations": [1]}},
                    "endpoints": [[0], [1]],
                },
            ]
        },
        "expected": [
            "DistributedTE does not use ShortestPathCache",
            "Pinned does not use ShortestPathCache",
            "CentralizedTE does not use ShortestPathCache",
            "prefer_headroom requires tunnel_count",
            "allocation=WEIGHTED source_allocation=EQUAL weight_basis=CAPACITY tie_break=FEWER_HOPS multipath=True parallel_edges=True",
            "tunnel count must be in [1, 10000]",
            "program endpoints do not match policy endpoints",
            "valid",
        ],
    },
    {
        "id": "EVAL-offer-infinite",
        "oracle": "evaluation_validation",
        "input": {
            "cases": [
                {
                    "model": "CentralizedTE",
                    "kwargs": {"source_allocation": "EQUAL", "prefer_headroom": True},
                    "offer": math.inf,
                },
                {
                    "model": "CentralizedTE",
                    "kwargs": {"source_allocation": "EQUAL", "prefer_headroom": True},
                    "offer": "4",
                },
                {"model": "ShortestPaths", "kwargs": {}, "offer": math.inf},
            ]
        },
        "expected": [
            "prefer_headroom requires tunnel_count",
            "TypeError",
            "offer must be an integer in [0, 2^62)",
        ],
    },
    {
        "id": "CFG-ceiling-type",
        "oracle": "model_config",
        "input": {
            "cases": [
                {
                    "model": "DistributedTE",
                    "kwargs": {"tunnel_count": "0", "cost_ceiling": ["other", 1.5]},
                },
                {
                    "model": "DistributedTE",
                    "kwargs": {"tunnel_count": "4", "cost_ceiling": ["other", 2]},
                },
                {
                    "model": "DistributedTE",
                    "kwargs": {"tunnel_count": "4", "cost_ceiling": ["fraction", 3, 2]},
                },
            ]
        },
        "expected": [
            "TypeError: cost_ceiling must be AbsoluteCost, CostStretch or None",
            "TypeError: cost_ceiling must be AbsoluteCost, CostStretch or None",
            "valid",
        ],
    },
    {
        "id": "KSP-ceiling-type",
        "oracle": "search",
        "input": {
            "edges": [[0, 1, 10, 5]],
            "source": 0,
            "destinations": 1,
            "k": 2,
            "residual": [10, 10],
            "stretch": ["other", 1.5],
        },
        "expected": {
            "error": "TypeError: cost_ceiling must be AbsoluteCost, CostStretch or None"
        },
    },
    {
        "id": "T1-stretch-below-one",
        "oracle": "t1",
        "input": {
            "edges": [[0, 1, 10, 5]],
            "heads": [0],
            "destinations": [1],
            "design_demand": 1,
            "per_head": 1,
            "offer": 1,
            "stretch": ["fraction", 1, 2],
        },
        "expected": {"error": "cost_ceiling CostStretch factor must be >= 1 (got 1/2)"},
    },
    {
        "id": "KSP-stretch-below-one",
        "oracle": "search",
        "input": {
            "edges": [[0, 1, 10, 5]],
            "source": 0,
            "destinations": 1,
            "k": 2,
            "stretch": ["fraction", 1, 2],
        },
        "expected": {"error": "cost_ceiling CostStretch factor must be >= 1 (got 1/2)"},
    },
    {
        "id": "KSP-ceiling-before-residual",
        "oracle": "search",
        "input": {
            "edges": [[0, 1, 10, 5]],
            "source": 0,
            "destinations": 1,
            "k": 2,
            "residual": [10, 10],
            "stretch": ["fraction", 1, 2],
        },
        "expected": {"error": "cost_ceiling CostStretch factor must be >= 1 (got 1/2)"},
    },
    {
        "id": "KSP-residual-length",
        "oracle": "search",
        "input": {
            "edges": [[0, 1, 10, 5]],
            "source": 0,
            "destinations": 1,
            "k": 2,
            "residual": [10, 10],
            "stretch": ["fraction", 3, 2],
        },
        "expected": {"error": "residual has length 2, expected 1 (one entry per edge)"},
    },
    {
        "id": "TIE-parallel",
        "oracle": "tie_choice",
        "input": {
            "edges": [[0, 1, 4, 1], [0, 1, 10, 1], [0, 1, 100, 1], [0, 1, 1000, 2]],
            "source": 0,
            "destination": 1,
            "available": [4, 8, 5, 1000],
        },
        "expected": {
            "FEWER_HOPS": [0],
            "HIGHER_CAPACITY": [2],
            "HIGHER_AVAILABLE_CAPACITY": [1],
            "HIGHER_HEADROOM": [0],
        },
    },
    {
        "id": "TIE-routes",
        "oracle": "tie_choice",
        "input": {
            "edges": [[0, 1, 10, 1], [0, 2, 3, 2], [1, 2, 10, 1]],
            "source": 0,
            "destination": 2,
        },
        "expected": {
            "FEWER_HOPS": [1],
            "HIGHER_CAPACITY": [0, 2],
            "HIGHER_AVAILABLE_CAPACITY": [0, 2],
            "HIGHER_HEADROOM": [1],
        },
    },
    {
        "id": "TIE-bottleneck",
        "oracle": "tie_choice",
        "input": {
            "edges": [[0, 1, 10, 1], [0, 2, 3, 2], [1, 2, 1, 1]],
            "source": 0,
            "destination": 2,
        },
        "expected": {
            "FEWER_HOPS": [1],
            "HIGHER_CAPACITY": [1],
            "HIGHER_AVAILABLE_CAPACITY": [1],
            "HIGHER_HEADROOM": [1],
        },
    },
    {
        "id": "EVAL-fixed-offer",
        "oracle": "evaluate_offer",
        "input": {"capacities": [10, 10], "offer": 4},
        "expected": {"baseline": 4, "delivery_reduction": [0, 0]},
    },
    {
        "id": "EVAL-first-hop",
        "oracle": "evaluate_offer",
        "input": {"capacities": [10, 10]},
        "expected": {"baseline": 20, "delivery_reduction": [10, 10]},
    },
    {
        "id": "EVAL-first-hop-unequal",
        "oracle": "evaluate_offer",
        "input": {"capacities": [4, 10]},
        "expected": {
            "baseline": 14,
            "delivery_reduction": [4, 10],
            "requested": math.inf,
            "offered": 14,
            "edge_drops": [0, 0],
        },
    },
    {
        "id": "EVAL-fixed-offer-unequal",
        "oracle": "evaluate_offer",
        "input": {"capacities": [4, 10], "offer": 12},
        "expected": {
            "baseline": 10,
            "delivery_reduction": [0, 6],
            "requested": 12,
            "offered": 12,
            "edge_drops": [2, 0],
        },
    },
    {
        "id": "DAG-pooled",
        "oracle": "dag_contexts",
        "input": {
            "edges": [[0, 2, 10, 1], [1, 2, 10, 1], [2, 3, 1, 1], [2, 3, 1, 1]],
            "sources": [0, 1],
            "destination": 3,
            "offer": 2,
            "pooled": True,
        },
        "expected": {"lossless_bound": 2, "lossy_delivered": 2},
    },
    {
        "id": "DAG-separate",
        "oracle": "dag_contexts",
        "input": {
            "edges": [[0, 2, 10, 1], [1, 2, 10, 1], [2, 3, 1, 1], [2, 3, 1, 1]],
            "sources": [0, 1],
            "destination": 3,
            "offer": 2,
            "pooled": False,
        },
        "expected": {"lossless_bound": 1, "lossy_delivered": 1},
    },
    {
        "id": "STRETCH-values",
        "oracle": "stretch_values",
        "input": {
            "factors": [
                ["float", 1.2],
                ["float", 1.1],
                ["float", 1.5],
                ["float", 0.30000000000000004],
                ["fraction", 6, 5],
                ["fraction", 18446744073709551616, 4611686018427387904],
                ["int", 2],
                ["float", math.inf],
                ["float", 1e19],
                ["fraction", 1, 9223372036854775808],
                ["bool", True],
                ["index", 3],
                ["str", "1.2"],
            ],
            "equal_pairs": [
                [["float", 1.5], ["fraction", 3, 2]],
                [["float", 1.2], ["fraction", 6, 5]],
                [["float", 1.2], ["fraction", 5404319552844595, 4503599627370496]],
            ],
        },
        "expected": {
            "values": [
                "6/5",
                "11/10",
                "3/2",
                "7500000000000001/25000000000000000",
                "6/5",
                "4/1",
                "2/1",
                "CostStretch factor must be finite (got inf)",
                "CostStretch fraction numerator must be in [-2^63, 2^63)",
                "CostStretch fraction denominator must be in [1, 2^63)",
                "TypeError",
                "3/1",
                "TypeError",
            ],
            "equal": [True, True, False],
        },
    },
    {
        "id": "STRETCH-numpy",
        "oracle": "stretch_values",
        "input": {
            "factors": [
                ["numpy", "float64", 1.2],
                ["numpy", "float32", 1.2],
                ["numpy", "float16", 1.2],
                ["numpy", "float32", 0.30000001192092896],
                ["numpy", "int64", 2],
                ["numpy", "float64", math.nan],
                ["numpy", "float32", math.inf],
                ["numpy", "float64", 1e19],
                ["numpy", "bool_", True],
            ],
            "equal_pairs": [
                [["numpy", "float32", 1.2], ["fraction", 6, 5]],
                [["numpy", "float64", 1.2], ["float", 1.2]],
            ],
        },
        "expected": {
            "values": [
                "6/5",
                "6/5",
                "6/5",
                "3/10",
                "2/1",
                "CostStretch factor must be finite (got nan)",
                "CostStretch factor must be finite (got inf)",
                "CostStretch fraction numerator must be in [-2^63, 2^63)",
                "TypeError",
            ],
            "equal": [True, True],
        },
    },
    {
        "id": "KSP-stretch-numpy",
        "oracle": "search",
        "input": {
            "edges": [[0, 1, 0, 100000000], [0, 1, 10, 120000001]],
            "source": 0,
            "destinations": 1,
            "k": 2,
            "residual": [0, 10],
            "stretch": ["numpy", "float32", 1.2],
        },
        "expected": {"paths": []},
    },
    {
        "id": "KSP-stretch-numpy-equal",
        "oracle": "search",
        "input": {
            "edges": [[0, 1, 0, 100000000], [0, 1, 10, 120000000]],
            "source": 0,
            "destinations": 1,
            "k": 2,
            "residual": [0, 10],
            "stretch": ["numpy", "float32", 1.2],
        },
        "expected": {"paths": [[120000000, [1]]]},
    },
    {
        "id": "KSP-stretch-boundary",
        "oracle": "search",
        "input": {
            "edges": [[0, 1, 0, 5], [0, 1, 10, 6]],
            "source": 0,
            "destinations": 1,
            "k": 2,
            "residual": [0, 10],
            "stretch": ["float", 1.2],
        },
        "expected": {"paths": [[6, [1]]]},
    },
    {
        "id": "KSP-stretch-below",
        "oracle": "search",
        "input": {
            "edges": [[0, 1, 0, 5], [0, 1, 10, 6]],
            "source": 0,
            "destinations": 1,
            "k": 2,
            "residual": [0, 10],
            "stretch": ["fraction", 11, 10],
        },
        "expected": {"paths": []},
    },
    {
        "id": "TERM-T1-all-down",
        "oracle": "termination",
        "input": {
            "offer": 14,
            "model": "t1",
            "args": {
                "edges": [
                    [0, 1, 6, 1],
                    [0, 2, 6, 3],
                    [1, 2, 4, 1],
                    [1, 3, 6, 3],
                    [2, 3, 6, 1],
                ],
                "heads": [0],
                "destinations": [3],
                "design_demand": 14,
                "per_head": 2,
            },
        },
        "expected": {
            "tunnels": 2,
            "lossy": "COMPLETE",
            "lossless": "EXHAUSTED",
            "source_drops": [14],
            "forward_error": None,
        },
    },
    {
        "id": "TERM-T2-all-dead",
        "oracle": "termination",
        "input": {
            "offer": 10,
            "model": "t2",
            "args": {
                "edges": [[0, 2, 0, 1], [1, 2, 0, 1]],
                "heads": [0, 1],
                "destinations": [2],
                "per_head": 2,
                "design_demand": 10,
            },
        },
        "expected": {
            "tunnels": 0,
            "lossy": "COMPLETE",
            "lossless": "EXHAUSTED",
            "source_drops": [5, 5],
            "forward_error": None,
        },
    },
    {
        "id": "TERM-zero-weight-route",
        "oracle": "termination",
        "input": {"offer": 4, "program": [[0], [[0, True, 0]]]},
        "expected": {
            "tunnels": 1,
            "lossy": "NO_ACTIVE_WEIGHT",
            "lossless": "NO_ACTIVE_WEIGHT",
            "source_drops": [4],
            "forward_error": None,
        },
    },
    {
        "id": "T2-no-tunnels-forward",
        "oracle": "termination",
        "input": {
            "offer": 10,
            "model": "t2",
            "args": {
                "edges": [[0, 2, 1, 1], [1, 2, 100, 1]],
                "heads": [0, 1],
                "destinations": [2],
                "per_head": 16,
            },
        },
        "expected": {
            "tunnels": 0,
            "lossy": "COMPLETE",
            "lossless": "EXHAUSTED",
            "source_drops": [5, 5],
            "forward_error": None,
        },
    },
    {
        "id": "ERR-unreserved-texts",
        "oracle": "input_texts",
        "input": {
            "cases": [
                {
                    "min_unreserved_bandwidth": 4611686018427387904,
                    "capacities": [10, 10],
                    "reservation_baseline": [[1, 2, 3]],
                },
                {
                    "min_unreserved_bandwidth": -1,
                    "capacities": [10, 10],
                    "reservation_baseline": [[5, 11]],
                },
                {
                    "min_unreserved_bandwidth": 0,
                    "capacities": [10, 10],
                    "reservation_baseline": [[10, 0]],
                },
            ]
        },
        "expected": [
            [
                "min_unreserved_bandwidth must be an integer in [0, 2^62)",
                "reservation_baseline has length 3, expected 2 (one entry per edge)",
            ],
            [
                "min_unreserved_bandwidth must be an integer in [0, 2^62)",
                "reservation_baseline[1] = 11 is outside [0, edge_capacity_units[1] = 10]",
            ],
            [None, None],
        ],
    },
    {
        "id": "ERR-class-search",
        "oracle": "input_classes",
        "input": {
            "cases": [
                {
                    "entry": "spf min_residual",
                    "values": [
                        ["bool", 1],
                        ["float", 1.0],
                        ["str", "1"],
                        ["none"],
                        ["float", "inf"],
                        0,
                        -3,
                        1,
                    ],
                },
                {
                    "entry": "ksp k",
                    "values": [
                        ["bool", 1],
                        ["float", "inf"],
                        ["str", "2"],
                        0,
                        -5,
                        1180591620717411303424,
                        3,
                    ],
                },
            ]
        },
        "expected": [
            [
                "TypeError: min_residual must be an integer (got bool)",
                "TypeError: min_residual must be an integer (got float)",
                "TypeError: min_residual must be an integer (got str)",
                "TypeError: min_residual must be an integer (got NoneType)",
                "TypeError: min_residual must be an integer (got float)",
                "ValueError: min_residual must be at least 1 unit (got 0)",
                "ValueError: min_residual must be at least 1 unit (got -3)",
                "valid",
            ],
            [
                "TypeError: k must be an integer (got bool)",
                "TypeError: k must be an integer (got float)",
                "TypeError: k must be an integer (got str)",
                "ValueError: k must be >= 1 (got 0)",
                "ValueError: k must be >= 1 (got -5)",
                "valid",
                "valid",
            ],
        ],
    },
    {
        "id": "ERR-class-call",
        "oracle": "input_classes",
        "input": {
            "cases": [
                {
                    "entry": "place offer",
                    "values": [
                        ["bool", 1],
                        ["bool", 0],
                        ["str", "1"],
                        ["none"],
                        ["fraction", 1, 1],
                        ["float", 1.5],
                        ["float", 0.0],
                        ["float", "-inf"],
                        ["float", "nan"],
                        -1,
                        4611686018427387904,
                        9223372036854775808,
                        ["float", "inf"],
                        0,
                        4611686018427387903,
                    ],
                }
            ]
        },
        "expected": [
            [
                "TypeError: offer must be an integer in [0, 2^62) or permitted +inf",
                "TypeError: offer must be an integer in [0, 2^62) or permitted +inf",
                "TypeError: offer must be an integer in [0, 2^62) or permitted +inf",
                "TypeError: offer must be an integer in [0, 2^62) or permitted +inf",
                "TypeError: offer must be an integer in [0, 2^62) or permitted +inf",
                "ValueError: offer must be an integer in [0, 2^62) or permitted +inf",
                "ValueError: offer must be an integer in [0, 2^62) or permitted +inf",
                "ValueError: offer must be an integer in [0, 2^62) or permitted +inf",
                "ValueError: offer must be an integer in [0, 2^62) or permitted +inf",
                "ValueError: offer must be an integer in [0, 2^62) or permitted +inf",
                "ValueError: offer must be an integer in [0, 2^62) or permitted +inf",
                "ValueError: offer must be an integer in [0, 2^62) or permitted +inf",
                "valid",
                "valid",
                "valid",
            ]
        ],
    },
    {
        "id": "ERR-class-vector",
        "oracle": "input_classes",
        "input": {
            "cases": [
                {
                    "entry": "service_baseline",
                    "caps": [10, 10],
                    "values": [
                        ["list", ["bool", 1], 0],
                        ["list", ["float", 1.0], 0],
                        ["list", 0, ["str", "0"]],
                        ["str", "ab"],
                        5,
                        ["list", ["list", 0], 0],
                        ["list", 0],
                        ["list", ["bool", 1]],
                        ["list", 0, 11],
                        ["list", -1, 0],
                        ["list", 1180591620717411303424, 0],
                        ["list", 3, 10],
                    ],
                }
            ]
        },
        "expected": [
            [
                "TypeError: service_baseline (units)[0] must be an integer (got bool)",
                "TypeError: service_baseline (units)[0] must be an integer (got float)",
                "TypeError: service_baseline (units)[1] must be an integer (got str)",
                "TypeError: service_baseline must be a 1-D sequence of integers (got str)",
                "TypeError: service_baseline must be a 1-D sequence of integers (got int)",
                "TypeError: service_baseline must be a 1-D sequence of integers (entry 0 is a sequence)",
                "ValueError: service_baseline has length 1, expected 2 (one entry per edge)",
                "ValueError: service_baseline has length 1, expected 2 (one entry per edge)",
                "ValueError: service_baseline[1] = 11 is outside [0, edge_capacity_units[1] = 10]",
                "ValueError: service_baseline[0] = -1 is outside [0, edge_capacity_units[0] = 10]",
                "ValueError: service_baseline[0] = 1180591620717411303424 is outside [0, edge_capacity_units[0] = 10]",
                "valid",
            ]
        ],
    },
    {
        "id": "ERR-class-numpy",
        "oracle": "input_classes",
        "input": {
            "cases": [
                {
                    "entry": "place offer",
                    "values": [
                        ["numpy", "bool_", True],
                        ["numpy", "float64", 1.5],
                        ["numpy", "float32", "-inf"],
                        ["numpy", "float64", "inf"],
                        ["numpy", "int64", 2],
                        ["numpy", "uint64", 9223372036854775808],
                    ],
                },
                {
                    "entry": "ksp k",
                    "values": [
                        ["numpy", "int64", 1],
                        ["numpy", "float64", 1.0],
                        ["numpy", "bool_", True],
                    ],
                },
                {
                    "entry": "service_baseline",
                    "caps": [10, 10],
                    "values": [
                        ["array", "float64", [1, 1]],
                        ["array", "int64", [[1, 1]]],
                        ["array", "float64", [1]],
                        ["array", "int64", [1, 11]],
                        ["array", "uint8", [1, 10]],
                    ],
                },
                {
                    "entry": "apply_actions",
                    "values": [
                        {"ledger": "own", "offer": ["numpy", "bool_", True]},
                        {"ledger": "own", "offer": ["numpy", "float32", 1.5]},
                        {"ledger": "own", "offer": ["numpy", "int64", 5]},
                    ],
                },
            ]
        },
        "expected": [
            [
                "TypeError: offer must be an integer in [0, 2^62) or permitted +inf",
                "ValueError: offer must be an integer in [0, 2^62) or permitted +inf",
                "ValueError: offer must be an integer in [0, 2^62) or permitted +inf",
                "valid",
                "valid",
                "ValueError: offer must be an integer in [0, 2^62) or permitted +inf",
            ],
            [
                "valid",
                "TypeError: k must be an integer (got float64)",
                "TypeError: k must be an integer (got bool)",
            ],
            [
                "TypeError: service_baseline (units) entries must be integers (got dtype float64)",
                "TypeError: service_baseline must be a 1-D array (got 2 dimensions)",
                "ValueError: service_baseline has length 1, expected 2 (one entry per edge)",
                "ValueError: service_baseline[1] = 11 is outside [0, edge_capacity_units[1] = 10]",
                "valid",
            ],
            [
                "TypeError: offer must be an int or a float (got bool)",
                "ValueError: offer must be an integer in [0, 2^62) or permitted +inf",
                "valid",
            ],
        ],
    },
    {
        "id": "ERR-number-resolution",
        "oracle": "number_types",
        "input": {
            "cases": [
                {
                    "entry": "resolution",
                    "values": [
                        1,
                        ["float", 0.5],
                        ["float", 2.3283064365386963e-10],
                        4294967296,
                        ["decimal", "1"],
                        ["fraction", 1, 1],
                        ["floatonly", 1.0],
                        ["index", 1],
                        ["bool", 1],
                        ["str", "1"],
                        ["none"],
                        ["float", 0.75],
                        ["float", "nan"],
                        ["float", "inf"],
                        8589934592,
                        ["pow2", 1024],
                        0,
                        -1,
                    ],
                }
            ]
        },
        "expected": [
            [
                "valid",
                "valid",
                "valid",
                "valid",
                "TypeError: pybind11 conversion",
                "TypeError: pybind11 conversion",
                "TypeError: pybind11 conversion",
                "TypeError: pybind11 conversion",
                "TypeError: pybind11 conversion",
                "TypeError: pybind11 conversion",
                "TypeError: pybind11 conversion",
                "ValueError: resolution must be a power of two in [2^-32, 2^32]",
                "ValueError: resolution must be a power of two in [2^-32, 2^32]",
                "ValueError: resolution must be a power of two in [2^-32, 2^32]",
                "ValueError: resolution must be a power of two in [2^-32, 2^32]",
                "ValueError: resolution must be a power of two in [2^-32, 2^32]",
                "ValueError: resolution must be a power of two in [2^-32, 2^32]",
                "ValueError: resolution must be a power of two in [2^-32, 2^32]",
            ]
        ],
    },
    {
        "id": "ERR-number-list",
        "oracle": "number_types",
        "input": {
            "cases": [
                {
                    "entry": "to_units",
                    "resolution": 1.0,
                    "values": [
                        ["list", ["pow2", 64], 1.5],
                        ["list", ["decimal", "1"], 1.5],
                        ["list", ["fraction", 3, 2], 1],
                        ["list", ["floatonly", 1.0]],
                        ["list", ["index", 1], 2],
                        ["list", ["decimal", "1"], ["pow2", 64]],
                        ["list", -1.5, ["pow2", 64]],
                        ["list", ["bool", 1], 1.5],
                        ["list", ["bool", 1], ["bool", 0]],
                        ["list", ["list", 1, 2], ["list", 3, 4]],
                        ["list", ["pow2", 63], 1.5],
                        ["list", -1, ["pow2", 63]],
                        ["list", ["float", "nan"], 1],
                        ["list", 1, 2.5, 0.49],
                        ["list"],
                    ],
                },
                {
                    "entry": "to_units",
                    "resolution": 0.5,
                    "values": [["list", 1, 2.25], ["list", ["decimal", "1"]]],
                },
                {
                    "entry": "to_units",
                    "resolution": ["decimal", "1"],
                    "values": [["list", ["pow2", 64], 1.5], ["list", 1.0]],
                },
                {
                    "entry": "to_units",
                    "resolution": 0.75,
                    "values": [["list", ["decimal", "1"]], ["list", ["pow2", 64], 1.5]],
                },
            ]
        },
        "expected": [
            [
                "ValueError: to_units takes physical numbers (got dtype object)",
                "TypeError: to_units takes physical numbers (got dtype object)",
                "TypeError: to_units takes physical numbers (got dtype object)",
                "TypeError: to_units takes physical numbers (got dtype object)",
                "TypeError: to_units takes physical numbers (got dtype object)",
                "TypeError: to_units takes physical numbers (got dtype object)",
                "ValueError: to_units takes physical numbers (got dtype object)",
                "TypeError: to_units takes physical numbers (a float or an array of floats)",
                "TypeError: to_units takes physical numbers (got dtype bool)",
                "TypeError: to_units takes a float or a 1-D array of floats",
                "ValueError: to_units: the amount is at or above the bound of the units domain (2^62 units)",
                "ValueError: to_units: a negative amount has no unit value",
                "ValueError: to_units: NaN has no unit value",
                "valid [1, 3, 0]",
                "valid []",
            ],
            [
                "valid [2, 5]",
                "TypeError: to_units takes physical numbers (got dtype object)",
            ],
            ["TypeError: pybind11 conversion", "TypeError: pybind11 conversion"],
            [
                "ValueError: resolution must be a power of two in [2^-32, 2^32]",
                "ValueError: resolution must be a power of two in [2^-32, 2^32]",
            ],
        ],
    },
    {
        "id": "ERR-number-scalar",
        "oracle": "number_types",
        "input": {
            "cases": [
                {
                    "entry": "to_units",
                    "resolution": 1.0,
                    "values": [
                        ["decimal", "1"],
                        ["pow2", 64],
                        ["fraction", 1, 1],
                        ["floatonly", 1.0],
                        ["index", 1],
                        ["bool", 1],
                        ["str", "1"],
                        ["none"],
                        ["float", "-inf"],
                        4611686018427387904,
                        4611686018427387903,
                        9007199254740993,
                        2.5,
                        ["float", 0.5],
                    ],
                },
                {
                    "entry": "place offer",
                    "values": [
                        ["index", 3],
                        ["decimal", "1"],
                        ["floatonly", 1.0],
                        ["pow2", 64],
                        3,
                    ],
                },
                {
                    "entry": "as_float x",
                    "values": [
                        ["index", 3],
                        ["decimal", "1"],
                        ["floatonly", 1.0],
                        ["float", 1.5],
                        ["pow2", 64],
                    ],
                },
            ]
        },
        "expected": [
            [
                "TypeError: to_units takes physical numbers (got dtype object)",
                "ValueError: to_units takes physical numbers (got dtype object)",
                "TypeError: to_units takes physical numbers (got dtype object)",
                "TypeError: to_units takes physical numbers (got dtype object)",
                "TypeError: to_units takes physical numbers (got dtype object)",
                "TypeError: to_units takes physical numbers (a float or an array of floats)",
                "TypeError: to_units takes physical numbers (a float or an array of floats)",
                "TypeError: to_units takes physical numbers (got dtype object)",
                "ValueError: to_units: an infinite amount is not finite and has no unit value",
                "ValueError: to_units: the amount is at or above the bound of the units domain (2^62 units)",
                "ValueError: to_units: the amount is at or above the bound of the units domain (2^62 units)",
                "valid 9007199254740992",
                "valid 3",
                "valid 1",
            ],
            [
                "valid",
                "TypeError: offer must be an integer in [0, 2^62) or permitted +inf",
                "TypeError: offer must be an integer in [0, 2^62) or permitted +inf",
                "ValueError: offer must be an integer in [0, 2^62) or permitted +inf",
                "valid",
            ],
            [
                "valid",
                "TypeError: as_float: x must be an integer (got Decimal)",
                "TypeError: as_float: x must be an integer (got FloatOnly)",
                "TypeError: as_float: x must be an integer (got float)",
                "valid",
            ],
        ],
    },
    {
        "id": "ERR-number-numpy",
        "oracle": "number_types",
        "input": {
            "cases": [
                {
                    "entry": "resolution",
                    "values": [
                        ["numpy", "float32", 0.5],
                        ["numpy", "int64", 2],
                        ["numpy", "uint64", 1],
                        ["numpy", "bool_", True],
                        ["array", "float64", 1.0],
                        ["numpy", "float64", "nan"],
                        ["numpy", "float16", 0.75],
                    ],
                },
                {
                    "entry": "to_units",
                    "resolution": 1.0,
                    "values": [
                        ["array", "float32", [1.5, 2.5]],
                        ["array", "object", [1, 2]],
                        ["array", "bool", [True, False]],
                        ["list", ["numpy", "bool_", True], 1.5],
                        ["list", ["numpy", "bool_", True], ["numpy", "bool_", False]],
                        ["array", "int64", [[1, 2]]],
                        ["array", "float64", 2.5],
                        ["numpy", "float32", 2.5],
                        ["list", ["numpy", "uint64", 9223372036854775808], -1],
                    ],
                },
            ]
        },
        "expected": [
            [
                "valid",
                "valid",
                "valid",
                "TypeError: pybind11 conversion",
                "TypeError: pybind11 conversion",
                "ValueError: resolution must be a power of two in [2^-32, 2^32]",
                "ValueError: resolution must be a power of two in [2^-32, 2^32]",
            ],
            [
                "valid [2, 3]",
                "TypeError: to_units takes physical numbers (got dtype object)",
                "TypeError: to_units takes physical numbers (got dtype bool)",
                "TypeError: to_units takes physical numbers (a float or an array of floats)",
                "TypeError: to_units takes physical numbers (got dtype bool)",
                "TypeError: to_units takes a float or a 1-D array of floats",
                "valid 3",
                "valid 3",
                "ValueError: to_units: the amount is at or above the bound of the units domain (2^62 units)",
            ],
        ],
    },
    # Several demands at once (apply_actions).
    {
        "id": "D-crossing",
        "oracle": "lossy",
        "input": {"offers": [8, 4, 4], "routes": [[0, 1], [0], [1]], "cap": [9, 7]},
        "expected": {
            "arrivals": [[8, 6], [4], [4]],
            "carried": [[6, 4], [3], [3]],
            "edge_drops": [[2, 2], [1], [1]],
            "delivered": [4, 3, 3],
            "edge_flows": [9, 7],
        },
    },
    {
        "id": "D-simultaneous",
        "oracle": "lossy",
        "input": {"offers": [10, 10], "routes": [[0], [0]], "cap": [10]},
        "expected": {
            "carried": [[5], [5]],
            "edge_drops": [[5], [5]],
            "delivered": [5, 5],
            "edge_flows": [10],
        },
    },
    {
        "id": "D-C3",
        "oracle": "lossy",
        "input": {
            "offers": [1, 1, 2],
            "routes": [[0, 1], [1, 2], [2, 0]],
            "cap": [1, 1, 2],
        },
        "expected": {
            "error": "lossy forwarding iteration entered a period-2 cycle",
            "rounds": 3,
        },
    },
    {
        "id": "D-cycle-prefix-suffix",
        "oracle": "lossy",
        "input": {
            "offers": [8, 4],
            "routes": [[0, 1, 2], [2, 1, 3]],
            "cap": [2, 5, 5, 10],
        },
        "expected": {
            "arrivals": [[8, 2, 2], [4, 3, 3]],
            "carried": [[2, 2, 2], [3, 3, 3]],
            "edge_drops": [[6, 0, 0], [1, 0, 0]],
            "delivered": [2, 3],
            "edge_flows": [2, 5, 5, 3],
        },
    },
    {
        "id": "D-cycle-subdivision-before",
        "oracle": "lossy",
        "input": {
            "offers": [2, 1, 4],
            "routes": [[1, 0, 2], [2, 0], [2, 0, 1]],
            "cap": [2, 4, 2],
        },
        "expected": {
            "arrivals": [[2, 2, 1], [1, 0], [4, 2, 1]],
            "carried": [[2, 1, 0], [0, 0], [2, 1, 1]],
            "delivered": [0, 0, 1],
            "edge_flows": [2, 3, 2],
        },
    },
    {
        "id": "D-cycle-subdivision-after",
        "oracle": "lossy",
        "input": {
            "offers": [2, 1, 4],
            "routes": [[1, 0, 2, 3], [2, 3, 0], [2, 3, 0, 1]],
            "cap": [2, 4, 2, 2],
        },
        "expected": {
            "arrivals": [[2, 2, 1, 0], [1, 0, 0], [4, 2, 2, 1]],
            "carried": [[2, 1, 0, 0], [0, 0, 0], [2, 2, 1, 1]],
            "delivered": [0, 0, 1],
            "edge_flows": [2, 3, 2, 2],
        },
    },
    {
        "id": "D-no-H-fixed-point",
        "oracle": "c3",
        "input": {
            "offers": [1, 1, 2],
            "routes": [[0, 1], [1, 2], [2, 0]],
            "cap": [1, 1, 2],
        },
        "expected": {"feasible_work_conserving_records": 5, "local_H_solutions": 0},
    },
    {
        "id": "STALE-cycle-period-2",
        "oracle": "stale_paths",
        "input": {
            "edges": [
                [0, 1, 768, 1],
                [1, 2, 25600, 1],
                [2, 3, 768, 1],
                [3, 0, 25600, 1],
                [4, 2, 25600, 1],
                [5, 0, 25600, 1],
                [6, 2, 25600, 1],
                [7, 8, 25600, 1],
            ],
            "tunnels": [[4, [4, 2, 3, 0]], [5, [5, 0, 1, 2]], [6, [6, 2]], [7, [7]]],
            "offers": {4: 1024, 5: 256, 6: 256, 7: 512},
            "black_links": [7],
        },
        "expected": {
            "error": "lossy forwarding iteration entered a period-2 cycle",
            "rounds": 8,
        },
    },
    {
        "id": "STALE-cycle-round-limit-baseline",
        "oracle": "stale_paths",
        "input": {
            "edges": [
                [0, 1, 12, 1],
                [1, 6, 1000, 1],
                [2, 3, 9, 1],
                [3, 0, 1000, 1],
                [3, 10, 1000, 1],
                [4, 5, 7, 1],
                [5, 6, 1000, 1],
                [6, 7, 13, 1],
                [7, 2, 1000, 1],
                [7, 10, 1000, 1],
                [8, 9, 1, 1],
                [10, 11, 10, 1],
                [11, 2, 1000, 1],
                [12, 4, 1000, 1],
                [13, 6, 1000, 1],
                [14, 0, 1000, 1],
            ],
            "tunnels": [
                [12, [13, 5, 6, 7, 8, 2, 4, 11]],
                [13, [14, 7, 9, 11, 12, 2, 3, 0]],
                [14, [15, 0, 1, 7]],
            ],
            "offers": {12: 13, 13: 11, 14: 10},
        },
        "expected": {"delivered": [3, 4, 4]},
    },
    {
        "id": "STALE-cycle-round-limit",
        "oracle": "stale_paths",
        "input": {
            "edges": [
                [0, 1, 12, 1],
                [1, 6, 1000, 1],
                [2, 3, 9, 1],
                [3, 0, 1000, 1],
                [3, 10, 1000, 1],
                [4, 5, 7, 1],
                [5, 6, 1000, 1],
                [6, 7, 13, 1],
                [7, 2, 1000, 1],
                [7, 10, 1000, 1],
                [8, 9, 1, 1],
                [10, 11, 10, 1],
                [11, 2, 1000, 1],
                [12, 4, 1000, 1],
                [13, 6, 1000, 1],
                [14, 0, 1000, 1],
            ],
            "tunnels": [
                [12, [13, 5, 6, 7, 8, 2, 4, 11]],
                [13, [14, 7, 9, 11, 12, 2, 3, 0]],
                [14, [15, 0, 1, 7]],
            ],
            "offers": {12: 13, 13: 11, 14: 10},
            "black_links": [5],
        },
        "expected": {
            "error": "lossy forwarding iteration did not reach a fixed point within the round limit",
            "rounds": 256,
        },
    },
    {
        "id": "STALE-black-hole-in-cycle",
        "oracle": "stale_paths",
        "input": {
            "edges": [
                [0, 1, 12, 1],
                [1, 6, 1000, 1],
                [1, 16, 1000, 1],
                [2, 3, 9, 1],
                [3, 0, 1000, 1],
                [3, 10, 1000, 1],
                [4, 5, 7, 1],
                [5, 6, 1000, 1],
                [6, 7, 13, 1],
                [7, 2, 1000, 1],
                [7, 10, 1000, 1],
                [8, 9, 1, 1],
                [10, 11, 10, 1],
                [11, 2, 1000, 1],
                [12, 4, 1000, 1],
                [13, 6, 1000, 1],
                [14, 0, 1000, 1],
                [15, 0, 1000, 1],
                [16, 2, 1000, 1],
            ],
            "tunnels": [
                [12, [14, 6, 7, 8, 9, 3, 5, 12]],
                [13, [15, 8, 10, 12, 13, 3, 4, 0]],
                [14, [16, 0, 1, 8]],
                [15, [17, 0, 2, 18, 3, 5, 12]],
            ],
            "offers": {12: 13, 13: 11, 14: 10, 15: 5},
            "black_links": [2, 6, 17],
        },
        "expected": {
            "status": ["BLACKHOLED", "UP", "UP", "BLACKHOLED"],
            "error": "lossy forwarding iteration did not reach a fixed point within the round limit",
            "rounds": 256,
        },
    },
    {
        "id": "LIST-black-hole-in-cycle",
        "oracle": "action_list",
        "input": {
            "edges": [
                [0, 1, 12, 1],
                [1, 6, 1000, 1],
                [1, 16, 1000, 1],
                [2, 3, 9, 1],
                [3, 0, 1000, 1],
                [3, 10, 1000, 1],
                [4, 5, 7, 1],
                [5, 6, 1000, 1],
                [6, 7, 13, 1],
                [7, 2, 1000, 1],
                [7, 10, 1000, 1],
                [8, 9, 1, 1],
                [10, 11, 10, 1],
                [11, 2, 1000, 1],
                [12, 4, 1000, 1],
                [13, 6, 1000, 1],
                [14, 0, 1000, 1],
                [15, 0, 1000, 1],
                [16, 2, 1000, 1],
            ],
            "actions": [
                {
                    "source": 12,
                    "route": [14, 6, 7, 8, 9, 3, 5, 12],
                    "offer": 13,
                    "black_holes": [2, 6, 17],
                },
                {
                    "source": 13,
                    "route": [15, 8, 10, 12, 13, 3, 4, 0],
                    "offer": 11,
                    "black_holes": [2, 6, 17],
                },
                {
                    "source": 14,
                    "route": [16, 0, 1, 8],
                    "offer": 10,
                    "black_holes": [2, 6, 17],
                },
                {
                    "source": 15,
                    "route": [17, 0, 2, 18, 3, 5, 12],
                    "offer": 5,
                    "black_holes": [2, 6, 17],
                },
            ],
            "mode": "SIMULTANEOUS",
        },
        "expected": {
            "error": "lossy forwarding iteration did not reach a fixed point within the round limit",
            "rounds": 256,
        },
    },
    {
        "id": "LIST-black-hole-own",
        "oracle": "action_list",
        "input": {
            "edges": [[0, 1, 10, 1]],
            "actions": [
                {"source": 0, "route": [0], "offer": 10, "black_holes": [0]},
                {"source": 0, "route": [0], "offer": 10},
            ],
            "mode": "SIMULTANEOUS",
        },
        "expected": {
            "entries": [
                {
                    "delivered": 0,
                    "call_edge_drops": [[0, 10]],
                    "status": "BLACKHOLED",
                    "residual": [0],
                    "saturated_boundary_edges": [0],
                },
                {
                    "delivered": 10,
                    "call_edge_drops": [],
                    "status": "UP",
                    "residual": [0],
                    "saturated_boundary_edges": [0],
                },
            ]
        },
    },
    {
        "id": "LIST-black-hole-upstream",
        "oracle": "action_list",
        "input": {
            "edges": [[0, 1, 10, 1], [1, 2, 10, 1]],
            "actions": [
                {"source": 0, "route": [0, 1], "offer": 10, "black_holes": [1]},
                {"source": 0, "route": [0, 1], "offer": 10},
            ],
            "mode": "SIMULTANEOUS",
        },
        "expected": {
            "entries": [
                {
                    "delivered": 0,
                    "call_edge_drops": [[0, 5], [1, 5]],
                    "status": "BLACKHOLED",
                    "residual": [0, 5],
                    "saturated_boundary_edges": [0],
                },
                {
                    "delivered": 5,
                    "call_edge_drops": [[0, 5]],
                    "status": "UP",
                    "residual": [0, 5],
                    "saturated_boundary_edges": [0],
                },
            ]
        },
    },
    {
        "id": "LIST-outputs-ordered",
        "oracle": "action_list",
        "input": {
            "edges": [[0, 1, 10, 1]],
            "actions": [
                {"source": 0, "route": [0], "offer": 7},
                {"source": 0, "route": [0], "offer": 7},
            ],
            "mode": "ORDERED",
        },
        "expected": {
            "entries": [
                {
                    "delivered": 7,
                    "call_edge_drops": [],
                    "residual": [3],
                    "saturated_boundary_edges": [],
                },
                {
                    "delivered": 3,
                    "call_edge_drops": [[0, 4]],
                    "residual": [0],
                    "saturated_boundary_edges": [0],
                },
            ]
        },
    },
    {
        "id": "LIST-outputs-simultaneous",
        "oracle": "action_list",
        "input": {
            "edges": [[0, 2, 5, 1], [1, 2, 10, 1], [2, 3, 10, 1]],
            "actions": [
                {"source": 0, "route": [0, 2], "offer": 5},
                {"source": 1, "route": [1, 2], "offer": 5},
            ],
            "mode": "SIMULTANEOUS",
        },
        "expected": {
            "entries": [
                {
                    "delivered": 5,
                    "call_edge_drops": [],
                    "residual": [0, 5, 0],
                    "saturated_boundary_edges": [0],
                },
                {
                    "delivered": 5,
                    "call_edge_drops": [],
                    "residual": [0, 5, 0],
                    "saturated_boundary_edges": [2],
                },
            ]
        },
    },
    {
        "id": "LIST-first-error",
        "oracle": "list_validation",
        "input": {
            "cases": [
                {
                    "mode": "ORDERED",
                    "actions": [
                        {
                            "kind": "Replace",
                            "policy": "a",
                            "model": "ShortestPaths",
                            "admission": "LOSSY",
                            "offer": -1,
                        },
                        {
                            "kind": "Replace",
                            "policy": "x",
                            "model": "ShortestPaths",
                            "admission": "LOSSY",
                            "offer": 1,
                            "foreign": True,
                        },
                    ],
                },
                {
                    "mode": "ORDERED",
                    "actions": [
                        {
                            "kind": "Replace",
                            "policy": "a",
                            "model": "ShortestPaths",
                            "admission": "LOSSY",
                            "offer": 1,
                        },
                        {
                            "kind": "Replace",
                            "policy": "a",
                            "model": "ShortestPaths",
                            "admission": "LOSSY",
                            "offer": -1,
                        },
                    ],
                },
                {
                    "mode": "SIMULTANEOUS",
                    "actions": [
                        {
                            "kind": "Forward",
                            "policy": "a",
                            "model": "ShortestPaths",
                            "admission": "LOSSY",
                            "offer": "inf",
                        }
                    ],
                },
                {
                    "mode": "SIMULTANEOUS",
                    "actions": [
                        {
                            "kind": "Place",
                            "policy": "a",
                            "model": "ShortestPaths",
                            "admission": "LOSSLESS",
                            "offer": 1,
                        }
                    ],
                },
                {
                    "mode": "SIMULTANEOUS",
                    "actions": [
                        {
                            "kind": "Replace",
                            "policy": "t",
                            "model": "DistributedTE",
                            "admission": "LOSSLESS",
                            "offer": 1,
                        }
                    ],
                },
                {
                    "mode": "SIMULTANEOUS",
                    "actions": [
                        {
                            "kind": "Forward",
                            "policy": "c",
                            "model": "CentralizedTE",
                            "admission": "LOSSLESS",
                            "offer": "inf",
                            "program": True,
                        }
                    ],
                },
                {
                    "mode": "SIMULTANEOUS",
                    "actions": [
                        {
                            "kind": "Replace",
                            "policy": "a",
                            "model": "Pinned",
                            "admission": "LOSSY",
                            "offer": 1,
                            "computes": "lossy forwarding iteration entered a period-2 cycle",
                        },
                        {
                            "kind": "Replace",
                            "policy": "b",
                            "model": "Pinned",
                            "admission": "LOSSY",
                            "offer": 4611686018427387904,
                        },
                    ],
                },
                {
                    "mode": "SIMULTANEOUS",
                    "actions": [
                        {
                            "kind": "Replace",
                            "policy": "a",
                            "model": "Pinned",
                            "admission": "LOSSY",
                            "offer": 1,
                            "computes": "lossy forwarding iteration entered a period-2 cycle",
                        },
                        {
                            "kind": "Replace",
                            "policy": "b",
                            "model": "Pinned",
                            "admission": "LOSSY",
                            "offer": 1,
                        },
                    ],
                },
                {
                    "mode": "ORDERED",
                    "actions": [
                        {
                            "kind": "Place",
                            "policy": "t",
                            "model": "DistributedTE",
                            "admission": "LOSSLESS",
                            "offer": -1,
                        }
                    ],
                },
                {
                    "mode": "ORDERED",
                    "actions": [
                        {
                            "kind": "Forward",
                            "policy": "p",
                            "model": "Pinned",
                            "admission": "LOSSY",
                            "offer": "inf",
                        }
                    ],
                },
                {
                    "mode": "ORDERED",
                    "actions": [
                        {
                            "kind": "Place",
                            "policy": "a",
                            "model": "ShortestPaths",
                            "admission": "LOSSY",
                            "offer": 5,
                            "after_infinite": True,
                        },
                        {
                            "kind": "Replace",
                            "policy": "x",
                            "model": "ShortestPaths",
                            "admission": "LOSSY",
                            "offer": 1,
                            "foreign": True,
                        },
                    ],
                },
                {
                    "mode": "ORDERED",
                    "actions": [
                        {
                            "kind": "Replace",
                            "policy": "a",
                            "model": "ShortestPaths",
                            "admission": "LOSSY",
                            "offer": 1,
                            "design_demand": 5,
                        },
                        {
                            "kind": "Replace",
                            "policy": "x",
                            "model": "ShortestPaths",
                            "admission": "LOSSY",
                            "offer": 1,
                            "foreign": True,
                        },
                    ],
                },
                {
                    "mode": "ORDERED",
                    "actions": [
                        {
                            "kind": "Forward",
                            "policy": "a",
                            "model": "ShortestPaths",
                            "admission": "LOSSY",
                            "offer": 1,
                            "program": True,
                            "foreign_black_holes": True,
                        },
                        {
                            "kind": "Replace",
                            "policy": "x",
                            "model": "ShortestPaths",
                            "admission": "LOSSY",
                            "offer": 1,
                            "foreign": True,
                        },
                    ],
                },
                {"mode": "SIMULTANEOUS", "actions": []},
            ]
        },
        "expected": [
            "offer must be an integer in [0, 2^62) or permitted +inf",
            "duplicate policy in list",
            "SIMULTANEOUS requires finite offers",
            "SIMULTANEOUS does not accept Place",
            "SIMULTANEOUS Replace requires ShortestPaths EQUAL or WEIGHTED, or Pinned",
            "SIMULTANEOUS requires admission=LOSSY",
            "offer must be an integer in [0, 2^62) or permitted +inf",
            "lossy forwarding iteration entered a period-2 cycle",
            "DistributedTE does not support incremental placement",
            "Pinned forwards with replace",
            "finite placement after infinite placement requires a new epoch",
            "ShortestPaths does not accept design_demand",
            "black_holes belongs to another graph",
            "valid",
        ],
    },
    {
        "id": "ERR-forwarding-limit",
        "oracle": "forwarding_size",
        "input": {
            "cases": [
                [[[10000, 214749]]],
                [[[10000, 214748], [1, 3647]]],
                [[[10000, 214748], [1, 3648]]],
                [[[1, 1073741824]]],
                [[[1, 1073741824]], [[1, 1073741824]]],
                [[[2147483648, 1]]],
            ]
        },
        "expected": [
            "forwarding occurrence limit exceeded: 2147490000 > 2147483647",
            "within the entry and occurrence limits",
            "forwarding occurrence limit exceeded: 2147483648 > 2147483647",
            "within the entry and occurrence limits",
            "forwarding occurrence limit exceeded: 2147483648 > 2147483647",
            "forwarding entry limit exceeded: 2147483648 > 2147483647",
        ],
    },
    {
        "id": "ERR-class-list",
        "oracle": "input_classes",
        "input": {
            "cases": [
                {
                    "entry": "apply_actions",
                    "values": [
                        {"ledger": "own", "offer": ["bool", 1]},
                        {"ledger": "own", "offer": ["str", "1"]},
                        {"ledger": "own", "offer": ["fraction", 3, 2]},
                        {"ledger": "own", "offer": ["float", 1.5]},
                        {"ledger": "own", "offer": 9223372036854775808},
                        {"ledger": "own", "offer": -1},
                        {"ledger": ["none"], "offer": 5},
                        {"ledger": ["none"], "offer": ["bool", 1]},
                        {"ledger": ["none"], "offer": -1},
                        {"ledger": "own", "offer": ["float", "inf"]},
                        {"ledger": "own", "offer": 5},
                    ],
                }
            ]
        },
        "expected": [
            [
                "TypeError: offer must be an int or a float (got bool)",
                "TypeError: offer must be an int or a float (got str)",
                "TypeError: offer must be an int or a float (got Fraction)",
                "ValueError: offer must be an integer in [0, 2^62) or permitted +inf",
                "ValueError: offer must be an integer in [0, 2^62) or permitted +inf",
                "ValueError: offer must be an integer in [0, 2^62) or permitted +inf",
                "TypeError: ledger must be a FlowLedger",
                "TypeError: offer must be an int or a float (got bool)",
                "TypeError: ledger must be a FlowLedger",
                "valid",
                "valid",
            ]
        ],
    },
]


# Values an example names that the library does not report, by example: left out of the comparison.
NOT_REPORTED = {
    "T2-dead-head": {
        "budget": "the supply budget of T2's count choice is internal to its up-to rule",
    },
    "D-C2": {
        "rounds": "the error names the period-2 cycle, not the round it was found in"
    },
    "D-C3": {
        "rounds": "the error names the period-2 cycle, not the round it was found in"
    },
    "STALE-cycle-period-2": {
        "rounds": "the error names the period-2 cycle, not the round it was found in"
    },
    "R4-refusal": {
        "requeue": "the rejected rule (a refused unit offered again), which the library never applies",
    },
}
for _name in (
    "STALE-cycle-round-limit",
    "STALE-black-hole-in-cycle",
    "LIST-black-hole-in-cycle",
):
    NOT_REPORTED[_name] = {"rounds": "the error names the round limit, not its value"}
for _name in (
    "C-G1-UP_TO",
    "T2-UP_TO-common",
    "T2-UP_TO-dead-head",
    "T2-tiny-head-UP_TO",
):
    NOT_REPORTED[_name] = {
        "bounds": "the per-count supply budgets of the up-to choice are internal",
        "scores": "the per-count scores of the up-to choice are internal",
    }

# Inputs outside the library's domain (amounts and weights are below 2^62). Such an example is
# not run on the library; it is checked three ways: its expected value is the rule's
# (reference_ha), the library refuses its input, and the same rounding at the largest inputs of
# the domain (its analogues) gives the library the rule's value.
OUTSIDE_THE_DOMAIN = {"R3-H-large": "a split of about 2^126 units"}
TOP = (1 << 62) - 1
# Each analogue: (input, expected); the C++ test (contract_examples_tests.cpp) uses the same values.
ANALOGUES = {
    "R3-H-large": [
        # Equal powers of two summing to n + 1, as in the example: the odd unit to the first.
        ({"n": TOP, "weights": [1 << 61, 1 << 61]}, [1 << 61, (1 << 61) - 1]),
        # The largest weights, equal and unequal (products of about 2^124).
        ({"n": TOP, "weights": [TOP, TOP]}, [1 << 61, (1 << 61) - 1]),
        ({"n": TOP, "weights": [TOP, TOP - 1, 1]}, [1 << 61, (1 << 61) - 1, 0]),
        (
            {"n": TOP - 1, "weights": [TOP, 1 << 61, 3]},
            [3074457345618258601, 1537228672809129300, 1],
        ),
        (
            {"n": TOP, "weights": [3, 5, 7]},
            [922337203685477580, 1537228672809129301, 2152120141932781022],
        ),
    ]
}


def reference_ha(n, weights):
    """H(n, w) as docs/semantics.md defines it ("Terms"), in exact integers: the floor of each
    share, then each remaining unit to the largest w / (assigned + 1), ties to the earlier member."""
    total = sum(weights)
    if not total:
        return [0] * len(weights)
    ans = [n * w // total for w in weights]
    for _ in range(n - sum(ans)):
        j = max(
            range(len(weights)), key=lambda i: (Fraction(weights[i], ans[i] + 1), -i)
        )
        ans[j] += 1
    return ans


def split_offer_on_unit_links(n, weights):
    """H(n, w) as the library splits a LOSSY offer: one head, a tunnel of weight w_i on its own
    link of capacity 1 (the split is of the offer; the links only limit the delivery)."""
    m = len(weights)
    g = ngc.Graph(2, [0] * m, [1] * m, [1] * m, [1] * m, 1.0)
    tunnels = [
        ngc.RouteDefinition(0, path_of(g, [e]), weight=w) for e, w in enumerate(weights)
    ]
    program = ngc.ForwardingProgram.from_routes(g, [0], [1], tunnels)
    return [t.assigned for t in pinned(g, program, LY).place(n).by_tunnel]


def refuses_outside_the_domain(n, weights):
    # A route weight and an offer beyond the domain are each refused.
    g = ngc.Graph(
        2,
        [0] * len(weights),
        [1] * len(weights),
        [1] * len(weights),
        [1] * len(weights),
        1.0,
    )
    for w in weights:
        with pytest.raises(ValueError, match="out of range|must be an integer in"):
            ngc.RouteDefinition(0, path_of(g, [0]), weight=w)
    program = ngc.ForwardingProgram.from_routes(
        g, [0], [1], [ngc.RouteDefinition(0, path_of(g, [0]))]
    )
    with pytest.raises(ValueError, match="offer must be an integer in"):
        pinned(g, program, LY).place(n)


def check_outside_the_domain(example):
    assert mismatch(reference_ha(**example["input"]), example["expected"]) is None
    refuses_outside_the_domain(**example["input"])
    for case, want in ANALOGUES[example["id"]]:
        assert mismatch(reference_ha(**case), want) is None, case
        assert mismatch(split_offer_on_unit_links(**case), want) is None, case


def _param(example):
    return pytest.param(example, id=example["id"])


def test_every_example_once():
    ids = [e["id"] for e in EXAMPLES]
    assert len(ids) == len(set(ids)) == 183
    assert {e["oracle"] for e in EXAMPLES} <= set(ADAPTERS)


@pytest.mark.parametrize("example", [_param(e) for e in EXAMPLES])
def test_contract_example(example):
    if example["id"] in OUTSIDE_THE_DOMAIN:
        check_outside_the_domain(example)
        return
    expected = example["expected"]
    left_out = NOT_REPORTED.get(example["id"], {})
    if left_out:
        expected = {k: v for k, v in expected.items() if k not in left_out}
    adapter = SEVERAL_DEMANDS.get(example["id"], ADAPTERS[example["oracle"]])
    got = adapter(**copy.deepcopy(example["input"]))
    assert mismatch(got, expected) is None
