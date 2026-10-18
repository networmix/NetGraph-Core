"""Prototype 2 fixtures in integer units (APPROACH revision 24, sections 8, 9 and 16).

The 44 fixtures of the grid prototype regenerated in units (the old names mapped in
``FIXTURE_MAP``), plus the review witnesses of review 13 to review 17. Graphs are
stated in input coordinates; bundles and expected per-edge arrays are input
positions translated through the inverse of ``input_index()``. Expected values are
derived in the fixture text or bodies (round robin, highest averages, hand-traced
transactions), never read back from the engine. "Manual" owners are ``Pinned``
policies with one bundle.

Graph shapes (input positions):
* ``diamond(s0, t0, s1, t1)``: 0 -[0]-> 1 -[1]-> 3 and 0 -[2]-> 2 -[3]-> 3; route 0
  = {0, 1}, route 1 = {2, 3};
* ``routes([(cap, cost), ...])``: route i is 0 -[2i]-> i+1 -[2i+1]-> d;
* ``parallel(caps)``: links 0 -> 1.
"""

from __future__ import annotations

from dataclasses import dataclass

import _p1  # noqa: F401
import certify
import policy as P
import search
from ledger import FlowGraph
from recorder import PSnap
from snap import dict_to_input
from units import INF, K_MAX_FLOWS_PER_CALL

MF, W, EB = P.MF, P.W_, P.EB
LOSSLESS, LOSSY = P.LOSSLESS, P.LOSSY

# Old fixture name -> regenerated fixture (the 44 grid fixtures).
FIXTURE_MAP = {
    "admission_zero_src_dst_no_route": "admission_zero_and_no_route",
    "cancellation_dag": "cancellation_dag",
    "combine_custom_edges_explicit": "sources_augmentation_edges",
    "combine_eb_pseudo_root": "sources_eb_per_source",
    "combine_forward_vs_reverse_3_vs_5": "sources_lossy_per_source",
    "combine_inf_every_kernel": "sources_inf_every_kernel",
    "combine_limits_and_validation": "sources_limits_and_validation",
    "combine_multi_source_attribution": "sources_attribution",
    "combine_proportional_2_vs_3": "sources_maxflow_vs_global",
    "combine_source_split_additivity": "sources_split_cumulative",
    "creation_fewer_routes": "creation_fewer_routes",
    "creation_retry_after_release": "creation_retry_after_release",
    "eb_1_plus_2_vs_3": "eb_1_plus_2_vs_3",
    "eb_three_members_g_1_1": "eb_correlated_remainders",
    "eb_unequal_1_10": "eb_unequal_1_10",
    "eight_versus_ten_capacity": "eight_versus_ten_capacity",
    "eight_versus_ten_headroom": "eight_versus_ten_headroom",
    "eight_versus_ten_residual": "eight_versus_ten_residual",
    "flag_table": "flag_table",
    "four_calls_EB": "four_calls_EB",
    "four_calls_Lossy": "four_calls_Lossy",
    "four_calls_P": "four_calls_MaxFlow",
    "four_calls_W": "four_calls_W",
    "growth_256_request_1278": "growth_256_request_1278",
    "growth_request_6": "growth_request_6",
    "harmful_growth_10_1": "harmful_growth_10_1",
    "hidden_beneficial_route": "hidden_beneficial_route",
    "inf_every_kernel": "inf_every_kernel",
    "invalidation_and_detach": "lifetime_remove_orphan_reset",
    "limits_diminishing": "limits_diminishing",
    "limits_max_flows_created": "limits_kmaxflowspercall",
    "lossy_diamond_10_then_1": "lossy_diamond_10_then_1",
    "lossy_shared_ingress": "lossy_shared_ingress",
    "pinned_bundle_rules": "pinned_bundle_rules",
    "pinned_kernels_5_1": "pinned_kernels_5_1",
    "reopt_double_spend": "reopt_double_spend",
    "request_12_on_5_5": "request_12_on_5_5",
    "retention_third_party_removal": "retention_third_party_removal",
    "rounded_retention_3_plus_g": "retention_same_structure",
    "shared_capacity_view": "shared_capacity_view",
    "ten_versus_five": "ten_versus_five",
    "third_party_occupancy": "third_party_occupancy",
    "weighted_1_plus_5_vs_6": "weighted_1_plus_5_vs_6",
    "weighted_frozen_third_party": "weighted_frozen_third_party",
}


@dataclass
class Fixture:
    name: str
    family: str
    rules: tuple[str, ...]
    description: str
    run: object


FIXTURES: list[Fixture] = []


def fixture(name, family, rules, description):
    def wrap(fn):
        FIXTURES.append(Fixture(name, family, tuple(rules), description, fn))
        return fn

    return wrap


# ------------------------------------------------------------------ helpers


def diamond(s0, t0, s1, t1, c0=1, c1=1):
    return 4, [(0, 1, s0, c0), (1, 3, t0, c0), (0, 2, s1, c1), (2, 3, t1, c1)]


def routes(specs):
    d = len(specs) + 1
    edges = []
    for i, (c, w) in enumerate(specs):
        edges += [(0, i + 1, c, w), (i + 1, d, c, w)]
    return d + 1, edges


def parallel(caps, costs=None):
    costs = costs or [1] * len(caps)
    return 2, [(0, 1, c, w) for c, w in zip(caps, costs, strict=True)]


def tunnel(priority="Capacity", tie=search.PREFER):
    return search.Selection(
        multipath=False, multi_edge=False, tie_break=tie, priority=priority
    )


def ca(kernel, priority="Capacity", fc=1, mx=None, reopt=False, sel=None, **kw):
    return P.FlowPolicyConfig(
        P.CapacityAware(flow_count=fc, max_flow_count=mx, reoptimize=reopt, **kw),
        kernel,
        LOSSLESS,
        sel or tunnel(priority),
    )


def co(kernel, adm=LOSSLESS, sel=None):
    return P.FlowPolicyConfig(P.CostOnly(), kernel, adm, sel or search.Selection())


def glob():
    return P.FlowPolicyConfig(P.Global(), MF, LOSSLESS, search.Selection())


class Ctx:
    """One graph, its caller-side inverse map, and helpers in input coordinates."""

    def __init__(self, s: PSnap, label, spec, resolution=1):
        n, edges = spec
        self.s, self.n, self.edges = s, n, edges
        self.g, self.inv = s.graph(label, n, edges, resolution)

    def ids(self, positions):
        return tuple(sorted(self.inv[p] for p in positions))

    def ledger(self, residual_input=None):
        if residual_input is None:
            return FlowGraph(self.g)
        return FlowGraph(self.g, [residual_input[p] for p in self.g.input_index()])

    def pin(self, kernel, bundles, adm=LOSSLESS, sel=None):
        return P.FlowPolicyConfig(
            P.Pinned(tuple(self.ids(b) for b in bundles)),
            kernel,
            adm,
            sel or search.Selection(),
        )

    def manual(self, fg, src, dst, positions, amount, label="manual owner"):
        """A fixed-structure owner: a Pinned MaxFlow policy with one bundle."""
        pol = P.FlowPolicy(fg, src, dst, self.pin(MF, [positions]))
        r = self.s.place(label, pol, amount, fg=fg)
        return pol, r

    def alloc(self, fg):
        return fg.state()["allocated_input"]

    def own(self, fg, pol):
        out = [0] * self.g.num_edges
        for r in pol.records():
            for e, x in r.edges.items():
                out[self.g.input_index()[e]] += x
        return out

    def routes_of(self, pol, src=None):
        st = pol.state()["cohorts"]
        src = next(iter(st)) if src is None else src
        return [f["route_input"] for f in st[src]]

    def shares(self, pol, src=None):
        st = pol.state()["cohorts"]
        src = next(iter(st)) if src is None else src
        return [f["T"] for f in st[src]]


# ====================================================================== creation


@fixture(
    "creation_fewer_routes",
    "creation",
    ["8b", "9.2"],
    "One capacity-5 link, CapacityAware flow_count 2, max 2, request 20: the "
    "first reservation exhausts the link, the second search finds nothing, so one "
    "flow is created (created flows are the ones that received a route); placed 5, "
    "exhausted (growth finds no candidate), routes_changed true.",
)
def creation_fewer_routes(s: PSnap):
    c = Ctx(s, "one_link", parallel([5]))
    for kernel in (EB, MF, W):
        fg = c.ledger()
        pol = P.FlowPolicy(fg, 0, 1, ca(kernel, fc=2, mx=2))
        r = s.place(f"{kernel}: request 20", pol, 20)
        s.expect(
            f"{kernel} result",
            (r.placed, r.termination, r.routes_changed),
            (5, "exhausted", True),
        )
        s.expect(f"{kernel} cohort size", len(pol.flow_ids()), 1)


@fixture(
    "creation_retry_after_release",
    "creation",
    ["8", "9.1", "9.4"],
    "Another owner holds the only link: a CapacityAware call finds no route, ends "
    "no_route, creates nothing and consumes no id (flags false). After the owner's "
    "remove() the next call creates the cohort (creation is attempted on every call "
    "while the cohort is empty), complete, routes_changed true.",
)
def creation_retry(s: PSnap):
    c = Ctx(s, "one_link", parallel([5]))
    fg = c.ledger()
    other, _ = c.manual(fg, 0, 1, [0], 5)
    pol = P.FlowPolicy(fg, 0, 1, ca(MF))
    nxt = fg.next_flow_id()
    r = s.place("first call", pol, 1)
    s.expect(
        "no route",
        (r.termination, r.routes_changed, r.weights_changed),
        ("no_route", False, False),
    )
    s.expect("no id consumed", fg.next_flow_id(), nxt)
    other.remove()
    r = s.place("second call", pol, 1)
    s.expect(
        "created", (r.placed, r.termination, r.routes_changed), (1, "complete", True)
    )
    s.expect("created id", pol.flow_ids()[0][1], nxt)


# ====================================================================== reservations


@fixture(
    "ten_versus_five",
    "reservations",
    ["8b", "9.2", "9.4"],
    "Two disjoint capacity-5 routes, EB CapacityAware, request 10: flow_count 2 "
    "places 10 (the second route is found after the first reservation, complete); "
    "flow_count 1 with max 1 places 5 and ends exhausted (a fixed-count cohort full; "
    "revision 20, the grid prototype said limit_flows).",
)
def ten_versus_five(s: PSnap):
    c = Ctx(s, "diamond_5", diamond(5, 5, 5, 5))
    fg = c.ledger()
    r = s.place("flow_count 2", P.FlowPolicy(fg, 0, 3, ca(EB, fc=2, mx=2)), 10, fg=fg)
    s.expect("two flows", (r.placed, r.termination), (10, "complete"))
    fg = c.ledger()
    r = s.place(
        "flow_count 1, max 1", P.FlowPolicy(fg, 0, 3, ca(EB, fc=1, mx=1)), 10, fg=fg
    )
    s.expect("one flow", (r.placed, r.termination), (5, "exhausted"))


def _eight_versus_ten(s: PSnap, priority, at8, term8, routes8):
    c = Ctx(s, "diamond_5", diamond(5, 5, 5, 5))
    for volume, expect, term in ((8, at8, term8), (10, 10, "complete")):
        fg = c.ledger()
        pol = P.FlowPolicy(fg, 0, 3, ca(EB, priority=priority, fc=2, mx=2))
        r = s.place(f"{priority}: request {volume}", pol, volume, fg=fg)
        s.expect(f"{priority} {volume}", (r.placed, r.termination), (expect, term))
        s.expect(
            f"{priority} reservations at {volume}",
            [x["reserved"] for x in r.trace["sources"][0]["stage1"]["reservations"]],
            [volume // 2] * 2
            if volume == 10
            else ([4, 1] if routes8 == "shared" else [4, 4]),
        )
        if volume == 8:
            want = [[0, 1], [0, 1]] if routes8 == "shared" else [[0, 1], [2, 3]]
            s.expect(f"{priority} routes at 8", c.routes_of(pol), want)


@fixture(
    "eight_versus_ten_capacity",
    "reservations",
    ["8b", "9.2"],
    "EB reservation target (P + v) div k with Capacity priority: request 8 reserves "
    "4 on route 0, whose residual 1 still ties on static capacity, so the second "
    "search takes route 0 again (reservation 1): both flows share route 0 and place "
    "5 (exhausted: fixed count); request 10 reserves 5, route 0 is gated, the second "
    "flow takes route 1: 10. The non-monotone CapacityAware fixture (9.2).",
)
def eight_versus_ten_capacity(s: PSnap):
    _eight_versus_ten(s, "Capacity", 5, "exhausted", "shared")


@fixture(
    "eight_versus_ten_residual",
    "reservations",
    ["8b", "9.2"],
    "Same with Residual priority: the second search prefers the untouched route "
    "(bottleneck 5 against 1): 8 and 10.",
)
def eight_versus_ten_residual(s: PSnap):
    _eight_versus_ten(s, "Residual", 8, "complete", "split")


@fixture(
    "eight_versus_ten_headroom",
    "reservations",
    ["8b", "9.2"],
    "Same with Headroom priority (exact rationals 1/5 against 5/5): 8 and 10.",
)
def eight_versus_ten_headroom(s: PSnap):
    _eight_versus_ten(s, "Headroom", 8, "complete", "split")


@fixture(
    "request_12_on_5_5",
    "reservations",
    ["8b", "9.2", "9.4"],
    "Infeasible reservations: request 12, flow_count 2 = max 2, two capacity-5 "
    "routes. EB targets 6 clip to 5 and 5; MaxFlow and Weighted sequential "
    "reservations 5 and 5; placed 10, remaining 2, exhausted (fixed count).",
)
def request_12(s: PSnap):
    c = Ctx(s, "diamond_5", diamond(5, 5, 5, 5))
    for kernel in (EB, MF, W):
        fg = c.ledger()
        pol = P.FlowPolicy(fg, 0, 3, ca(kernel, fc=2, mx=2))
        r = s.place(f"{kernel}: request 12", pol, 12, fg=fg)
        s.expect(
            f"{kernel} result",
            (r.placed, r.remaining, r.termination),
            (10, 2, "exhausted"),
        )
        res = r.trace["sources"][0]["stage1"]["reservations"]
        s.expect(f"{kernel} reservations", [x["reserved"] for x in res], [5, 5])


# ====================================================================== growth


@fixture(
    "growth_request_6",
    "growth",
    ["8b", "9.2"],
    "EB CapacityAware flow_count 1, max 2, two capacity-5 routes, request 6: the "
    "first flow takes 5, growth to two flows gives the cohort 6 split [3, 3] "
    "(no over-admission); complete, routes_changed true.",
)
def growth_6(s: PSnap):
    c = Ctx(s, "diamond_5", diamond(5, 5, 5, 5))
    fg = c.ledger()
    pol = P.FlowPolicy(fg, 0, 3, ca(EB, fc=1, mx=2))
    r = s.place("request 6", pol, 6)
    s.expect(
        "result",
        (r.placed, r.termination, r.routes_changed, r.weights_changed),
        (6, "complete", True, False),
    )
    s.expect("shares", c.shares(pol), [3, 3])


@fixture(
    "growth_256_request_1278",
    "growth",
    ["8b", "9.2"],
    "256 parallel capacity-5 links, EB tunnels (Residual), flow_count 1, max 256, "
    "request 1278: 255 accepted growth steps, k = 256, the cohort total 1278 split "
    "by round robin: 254 flows at 5 and 2 at 4; complete.",
)
def growth_256(s: PSnap):
    c = Ctx(s, "parallel_256", parallel([5] * 256))
    fg = c.ledger()
    pol = P.FlowPolicy(fg, 0, 1, ca(EB, priority="Residual", fc=1, mx=256))
    r = pol.place(1278)
    steps = r.trace["sources"][0]["stage4"]
    for st in steps:
        st.pop("W", None)
        st.pop("shares", None)
    s.step("request 1278", result=r.to_json(), growth=steps, created=r.trace["created"])
    s.expect("result", (r.placed, r.termination), (1278, "complete"))
    shares = c.shares(pol)
    s.expect(
        "k and shares", (len(shares), shares.count(5), shares.count(4)), (256, 254, 2)
    )
    s.expect(
        "accepted growth steps",
        sum(1 for x in steps if x.get("decision") == "accepted"),
        255,
    )
    s.expect("created", r.trace["created"], 256)


@fixture(
    "harmful_growth_10_1",
    "growth",
    ["8b", "9.2", "9.4"],
    "Routes of capacity 10 (cost 2) and 1 (cost 4), EB, max 2. P = 10 on route 0; "
    "place(1): the candidate on route 1 gives a cohort bound of 3 (round robin [2, 1] "
    "over capacities [10, 1]) < 10 and is rejected: delta 0, exhausted, flags false, "
    "no id consumed.",
)
def harmful_growth(s: PSnap):
    c = Ctx(s, "harmful", diamond(10, 10, 1, 1, c0=1, c1=2))
    fg = c.ledger()
    pol = P.FlowPolicy(fg, 0, 3, ca(EB, fc=1, mx=2))
    s.place("request 10", pol, 10)
    nxt = fg.next_flow_id()
    r = s.place("request 1", pol, 1)
    s.expect(
        "result",
        (r.placed, r.termination, r.routes_changed, r.weights_changed),
        (0, "exhausted", False, False),
    )
    rej = [
        x
        for x in r.trace["sources"][0]["stage4"]
        if x.get("decision", "").startswith("rejected")
    ]
    s.expect(
        "rejected candidate",
        [(x["route"], x["T_candidate"]) for x in rej],
        [(list(c.ids([2, 3])), 3)],
    )
    s.expect("no id consumed", fg.next_flow_id(), nxt)
    s.expect("total", pol.placed_demand(), 10)


@fixture(
    "hidden_beneficial_route",
    "growth",
    ["8b", "9.4"],
    "Routes (capacity, cost per edge) (10, 1), (1, 2), (10, 3); P = 10 on the "
    "first; place(1) evaluates the cheapest candidate (route 1), rejects it (bound "
    "3) and stops: exhausted under the greedy rule, although routes 0 and 2 jointly "
    "would total 11 (round robin [6, 5] fits [10, 10]; documented, not searched).",
)
def hidden_beneficial(s: PSnap):
    c = Ctx(s, "three_routes", routes([(10, 1), (1, 2), (10, 3)]))
    fg = c.ledger()
    pol = P.FlowPolicy(fg, 0, 4, ca(EB, fc=1, mx=2))
    s.place("request 10", pol, 10)
    r = s.place("request 1", pol, 1)
    s.expect("result", (r.placed, r.termination), (0, "exhausted"))
    tried = [x["route"] for x in r.trace["sources"][0]["stage4"] if "route" in x]
    s.expect("candidate tried", tried, [list(c.ids([2, 3]))])
    rr = P.round_robin(11, 2)
    s.expect(
        "alternative (not searched) fits",
        (rr, rr[0] <= 10 and rr[1] <= 10),
        ([6, 5], True),
    )


# ====================================================================== four calls


def _four_calls(s: PSnap, kernel, adm, expected):
    for shape, caps in (("equal", (5, 5)), ("unequal", (1, 10))):
        c = Ctx(s, f"diamond_{shape}", diamond(caps[0], caps[0], caps[1], caps[1]))
        for mode in ("CostOnly", "Pinned"):
            fg = c.ledger()
            cfg = (
                co(kernel, adm)
                if mode == "CostOnly"
                else c.pin(kernel, [[0, 1], [2, 3]], adm)
            )
            pol = P.FlowPolicy(fg, 0, 3, cfg)
            got = [
                s.place(f"{shape} {mode} call {i}", pol, 1, fg=fg).placed
                for i in range(4)
            ]
            s.expect(f"{shape} {mode}", got, expected[shape])


@fixture(
    "four_calls_MaxFlow",
    "four_calls",
    ["8a", "16.4"],
    "MaxFlow, four calls of 1 on equal (5/5) and unequal (1/10) routes, one "
    "CostOnly flow on the diamond and two Pinned flows: 1, 1, 1, 1.",
)
def four_calls_mf(s: PSnap):
    _four_calls(s, MF, LOSSLESS, {"equal": [1] * 4, "unequal": [1] * 4})


@fixture(
    "four_calls_W",
    "four_calls",
    ["8a", "16.4"],
    "Weighted (Residual weights at creation, 1:10 on the unequal diamond; highest "
    "averages gives [0, T] up to T = 9): 1, 1, 1, 1 everywhere.",
)
def four_calls_w(s: PSnap):
    _four_calls(s, W, LOSSLESS, {"equal": [1] * 4, "unequal": [1] * 4})


@fixture(
    "four_calls_EB",
    "four_calls",
    ["8a", "16.4"],
    "EqualBalanced: equal 1, 1, 1, 1; unequal 1/10: 1, 1, 0, 0 (round robin [2, 1] "
    "overfills the 1-route at T = 3).",
)
def four_calls_eb(s: PSnap):
    _four_calls(s, EB, LOSSLESS, {"equal": [1] * 4, "unequal": [1, 1, 0, 0]})


@fixture(
    "four_calls_Lossy",
    "four_calls",
    ["19", "16.6"],
    "Lossy EB (per-call offers, round robin from the lowest EdgeId every call): "
    "equal 1, 1, 1, 1; unequal 1, 0, 0, 0 (every one-unit offer goes to the full "
    "1-route and is dropped there; the grid prototype's halves do not exist in "
    "units).",
)
def four_calls_lossy(s: PSnap):
    _four_calls(s, EB, LOSSY, {"equal": [1] * 4, "unequal": [1, 0, 0, 0]})


@fixture(
    "eb_unequal_1_10",
    "four_calls",
    ["8a", "9.4"],
    "EB on the unequal 1/10 diamond (one CostOnly flow): 1, 1, 0, 0; the third and "
    "fourth calls are exhausted with nothing placed.",
)
def eb_unequal(s: PSnap):
    c = Ctx(s, "diamond_1_10", diamond(1, 1, 10, 10))
    fg = c.ledger()
    pol = P.FlowPolicy(fg, 0, 3, co(EB))
    rs = [s.place(f"call {i}", pol, 1) for i in range(4)]
    s.expect("placed", [r.placed for r in rs], [1, 1, 0, 0])
    s.expect(
        "terminations",
        [r.termination for r in rs],
        ["complete", "complete", "exhausted", "exhausted"],
    )


@fixture(
    "cancellation_dag",
    "four_calls",
    ["20", "16.4"],
    "CostOnly MaxFlow on the cancellation DAG (all costs 1): calls 1 then 1 give 1 "
    "and 1 (edge 1->3 goes 1 -> 0), a single call of 2 gives the same allocation "
    "[1,1,0,1,1,1,1].",
)
def cancellation(s: PSnap):
    pairs = [(0, 1), (0, 2), (1, 3), (1, 4), (2, 3), (3, 5), (4, 5)]
    c = Ctx(s, "cancellation", (6, [(u, v, 1, 1) for u, v in pairs]))
    fg = c.ledger()
    pol = P.FlowPolicy(fg, 0, 5, co(MF))
    a = s.place("call 1", pol, 1)
    s.expect("after call 1", c.alloc(fg), [1, 0, 1, 0, 0, 1, 0])
    b = s.place("call 2", pol, 1)
    s.expect("two calls", (a.placed, b.placed), (1, 1))
    s.expect("after call 2", c.alloc(fg), [1, 1, 0, 1, 1, 1, 1])
    fg2 = c.ledger()
    s.place("single call 2", P.FlowPolicy(fg2, 0, 5, co(MF)), 2, fg=fg2)
    s.expect("1 + 1 equals 2", c.alloc(fg2), c.alloc(fg))


# ====================================================================== Lossy


@fixture(
    "lossy_diamond_10_then_1",
    "lossy",
    ["19", "16.6", "15"],
    "CostOnly Lossy EB on source links 10/10 over tails 1/10: request 10 delivers "
    "6 (complete, 4 dropped on the 1-tail, remaining 4); request 1 is offered to the "
    "full branch and dropped (complete, delivered 0); the record keeps offered 11, "
    "delivered 6 and 5 drops on the tail; drops are never re-offered.",
)
def lossy_diamond(s: PSnap):
    c = Ctx(s, "diamond_10_1_10_10", diamond(10, 1, 10, 10))
    fg = c.ledger()
    pol = P.FlowPolicy(fg, 0, 3, co(EB, LOSSY))
    tail = c.ids([1])[0]
    a = s.place("request 10", pol, 10)
    s.expect(
        "first",
        (a.placed, a.termination, a.remaining, a.drops),
        (6, "complete", 4, [(tail, 4)]),
    )
    b = s.place("request 1", pol, 1)
    s.expect(
        "second",
        (b.placed, b.termination, b.remaining, b.drops),
        (0, "complete", 1, [(tail, 1)]),
    )
    rec = pol.records()[0]
    s.expect("record", (rec.offered, rec.delivered, rec.drops), (11, 6, {tail: 5}))


@fixture(
    "lossy_shared_ingress",
    "lossy",
    ["19", "16.6"],
    "Two Pinned Lossy EB flows share an ingress of capacity 2 (tails 1 and 10): "
    "requests 2 then 2 deliver 2 then 0 (the second offers are dropped on the full "
    "ingress, complete); one request 4 delivers 1 (the first flow spends the ingress "
    "on its 1-tail, the second flow's offer of 2 is dropped on the ingress). Not "
    "additive (rule 19).",
)
def lossy_shared(s: PSnap):
    c = Ctx(
        s,
        "shared_ingress",
        (4, [(0, 1, 2, 1), (1, 3, 1, 1), (1, 2, 10, 1), (2, 3, 10, 1)]),
    )
    bundles = [[0, 1], [0, 2, 3]]
    fg = c.ledger()
    pol = P.FlowPolicy(fg, 0, 3, c.pin(EB, bundles, LOSSY))
    a = s.place("request 2", pol, 2)
    b = s.place("request 2 again", pol, 2)
    s.expect(
        "2 then 2",
        ((a.placed, a.termination), (b.placed, b.termination)),
        ((2, "complete"), (0, "complete")),
    )
    fg2 = c.ledger()
    pol2 = P.FlowPolicy(fg2, 0, 3, c.pin(EB, bundles, LOSSY))
    r = s.place("request 4", pol2, 4, fg=fg2)
    s.expect("single 4", (r.placed, r.termination), (1, "complete"))
    recs = pol2.records()
    s.expect(
        "records (offered, delivered, drops)",
        [(x.offered, x.delivered, dict_to_input(c.g, x.drops)) for x in recs],
        [(2, 1, [0, 1, 0, 0]), (2, 0, [2, 0, 0, 0])],
    )


# ====================================================================== +inf


@fixture(
    "inf_every_kernel",
    "infinite",
    ["15", "8a", "8b"],
    "+inf for every kernel: CostOnly MaxFlow, Weighted and EB on 5/5 place 10 "
    "(exhausted, remaining +inf); EB and MaxFlow CapacityAware grow to two flows "
    "(10); Lossy injects every source member to its residual: 11 on the 1/10 "
    "diamond; on 10/10 over 1/10 delivered 11 with 20 offered and 9 dropped.",
)
def inf_every_kernel(s: PSnap):
    c = Ctx(s, "diamond_5", diamond(5, 5, 5, 5))
    for kernel in (MF, W, EB):
        fg = c.ledger()
        r = s.place(
            f"CostOnly {kernel} +inf", P.FlowPolicy(fg, 0, 3, co(kernel)), INF, fg=fg
        )
        s.expect(
            f"CostOnly {kernel}",
            (r.placed, r.remaining, r.termination),
            (10, None, "exhausted"),
        )
    for kernel in (EB, MF):
        fg = c.ledger()
        pol = P.FlowPolicy(fg, 0, 3, ca(kernel, fc=1))
        r = s.place(f"CapacityAware {kernel} +inf", pol, INF, fg=fg)
        s.expect(
            f"CapacityAware {kernel}",
            (r.placed, r.termination, len(pol.flow_ids())),
            (10, "exhausted", 2),
        )
    c2 = Ctx(s, "diamond_1_10", diamond(1, 1, 10, 10))
    fg = c2.ledger()
    r = s.place("Lossy +inf 1/10", P.FlowPolicy(fg, 0, 3, co(EB, LOSSY)), INF, fg=fg)
    s.expect("Lossy 1/10", (r.placed, r.termination), (11, "exhausted"))
    c3 = Ctx(s, "diamond_10_1_10_10", diamond(10, 1, 10, 10))
    fg = c3.ledger()
    pol = P.FlowPolicy(fg, 0, 3, co(EB, LOSSY))
    r = s.place("Lossy +inf transit", pol, INF, fg=fg)
    rec = pol.records()[0]
    s.expect(
        "Lossy transit", (r.placed, rec.offered, sum(rec.drops.values())), (11, 20, 9)
    )


# ====================================================================== reoptimisation


@fixture(
    "retention_third_party_removal",
    "reoptimisation",
    ["8b", "9.2"],
    "Two capacity-5 routes, another owner holds 2 on route 0. EB CapacityAware "
    "count 2 (max 2), Capacity priority, reoptimisation. Request 8 reserves 3 and 4 "
    "and places 6 ([3, 3]; exhausted). The owner leaves; request 1: flow 0 keeps "
    "route 0, flow 1's candidate (route 0, view 2) cannot carry its 3 and is "
    "rejected (retention); stage 3 raises the cohort to 7 ([4, 3]); flags false.",
)
def retention(s: PSnap):
    c = Ctx(s, "diamond_5", diamond(5, 5, 5, 5))
    fg = c.ledger()
    other, _ = c.manual(fg, 0, 3, [0, 1], 2)
    pol = P.FlowPolicy(fg, 0, 3, ca(EB, fc=2, mx=2, reopt=True))
    r = s.place("request 8", pol, 8)
    s.expect("first", (r.placed, r.termination), (6, "exhausted"))
    res = r.trace["sources"][0]["stage1"]["reservations"]
    s.expect("reservations", [x["reserved"] for x in res], [3, 4])
    other.remove()
    r = s.place("request 1", pol, 1)
    rows = r.trace["sources"][0]["stage2"]
    s.expect(
        "reopt decisions",
        [x["decision"].split(":")[0] for x in rows],
        ["kept", "rejected"],
    )
    s.expect("second", (r.placed, pol.placed_demand(), c.shares(pol)), (1, 7, [4, 3]))
    s.expect("flags", (r.routes_changed, r.weights_changed), (False, False))


@fixture(
    "shared_capacity_view",
    "reoptimisation",
    ["8b", "9.2"],
    "Shared ingress e0 (10) and parallel egresses e1, e2 (10 each); a third party "
    "holds 6 on e1. Two EB tunnels (Residual, Collapse picks the egress with the "
    "larger residual) both take e2: [4, 4]. After the third party leaves, request "
    "1: flow 0 sees e0 = 2 + 4 = 6 on its view W + own, moves to e1 at T_old = 4; "
    "the cohort reaches 9 ([5, 4]); routes_changed true, weights_changed false.",
)
def shared_view(s: PSnap):
    c = Ctx(s, "shared", (3, [(0, 1, 10, 1), (1, 2, 10, 1), (1, 2, 10, 1)]))
    fg = c.ledger()
    other, _ = c.manual(fg, 1, 2, [1], 6)
    pol = P.FlowPolicy(fg, 0, 2, ca(EB, priority="Residual", fc=2, mx=2, reopt=True))
    r = s.place("request 8", pol, 8)
    s.expect("first", (r.placed, c.routes_of(pol)), (8, [[0, 2], [0, 2]]))
    other.remove()
    r = s.place("request 1", pol, 1)
    row = r.trace["sources"][0]["stage2"][0]
    s.expect(
        "flow 0 adopted e1",
        (row["decision"], row["candidate"]),
        ("adopted", list(c.ids([0, 1]))),
    )
    s.expect("second", (r.placed, pol.placed_demand(), c.shares(pol)), (1, 9, [5, 4]))
    s.expect("flags", (r.routes_changed, r.weights_changed), (True, False))


@fixture(
    "reopt_double_spend",
    "reoptimisation",
    ["8b", "9.2"],
    "MaxFlow tunnel with reoptimisation on two parallel capacity-10 links: T = 4 on "
    "e0; a third party fills e0. Request 1: the reroute to e1 happens at T_old = 4 "
    "(no budget), stage 3 adds 1: placed 1, total 5. The single-link variant (no "
    "third party) places 1 as well.",
)
def double_spend(s: PSnap):
    c = Ctx(s, "parallel_10_10", parallel([10, 10]))
    fg = c.ledger()
    pol = P.FlowPolicy(fg, 0, 1, ca(MF, priority="Residual", fc=1, mx=1, reopt=True))
    s.place("request 4", pol, 4)
    c.manual(fg, 0, 1, [0], 6, "third party fills e0")
    r = s.place("request 1", pol, 1)
    s.expect("result", (r.placed, pol.placed_demand(), c.routes_of(pol)), (1, 5, [[1]]))
    s.expect("flags", (r.routes_changed, r.weights_changed), (True, False))
    c1 = Ctx(s, "one_link", parallel([10]))
    fg1 = c1.ledger()
    p1 = P.FlowPolicy(fg1, 0, 1, ca(MF, priority="Residual", fc=1, mx=1, reopt=True))
    s.place("single link request 4", p1, 4, fg=fg1)
    r = s.place("single link request 1", p1, 1, fg=fg1)
    s.expect("single link", (r.placed, p1.placed_demand()), (1, 5))


@fixture(
    "retention_same_structure",
    "reoptimisation",
    ["8b", "9.4"],
    "Three parallel capacity-2 links, EB CapacityAware multipath with "
    "reoptimisation: request 3 gives [1, 1, 1]; third parties fill the leftovers. "
    "Request 1: the candidate DAG equals the flow's route (a structurally equal "
    "candidate is not a reroute): kept, nothing shed, placed 0, exhausted, flags "
    "false. (Replaces the grid fixture 3 + g.)",
)
def retention_same(s: PSnap):
    c = Ctx(s, "parallel_2", parallel([2, 2, 2]))
    fg = c.ledger()
    pol = P.FlowPolicy(fg, 0, 1, ca(EB, fc=1, reopt=True, sel=search.Selection()))
    s.place("request 3", pol, 3)
    s.expect("allocation", c.own(fg, pol), [1, 1, 1])
    for p in range(3):
        c.manual(fg, 0, 1, [p], 1, f"third party fills link {p}")
    r = s.place("request 1", pol, 1)
    s.expect(
        "result",
        (r.placed, r.termination, r.routes_changed, r.weights_changed),
        (0, "exhausted", False, False),
    )
    s.expect(
        "reopt", [x["decision"] for x in r.trace["sources"][0]["stage2"]], ["kept"]
    )
    s.expect("nothing shed", pol.placed_demand(), 3)


# ====================================================================== third parties


@fixture(
    "weighted_frozen_third_party",
    "third_party",
    ["8a", "16.2", "9.4"],
    "CostOnly Weighted on source links 10/10 over tails t/10: weights [10, 10] "
    "captured at creation and frozen. (a) t = 1: 2 then 1 give 2 then 0 (T = 3 "
    "needs 2 on the 1-tail). (b) t = 2: request 2 ([1, 1]); a third party takes 1 "
    "on the tail; request 1 places 0; after the third party leaves, request 1 places "
    "1. The plan never changes; weights_changed stays false.",
)
def weighted_frozen(s: PSnap):
    c = Ctx(s, "diamond_10_1", diamond(10, 1, 10, 10))
    fg = c.ledger()
    pol = P.FlowPolicy(fg, 0, 3, co(W))
    a = s.place("(a) request 2", pol, 2)
    b = s.place("(a) request 1", pol, 1)
    s.expect("(a)", (a.placed, b.placed), (2, 0))
    c2 = Ctx(s, "diamond_10_2", diamond(10, 2, 10, 10))
    fg = c2.ledger()
    pol = P.FlowPolicy(fg, 0, 3, co(W))
    s.place("(b) request 2", pol, 2, fg=fg)
    plan = pol.state()["cohorts"][0][0]["plan"]
    other, _ = c2.manual(fg, 1, 3, [1], 1, "(b) third party on the tail")
    r1 = s.place("(b) request 1 blocked", pol, 1, fg=fg)
    other.remove()
    r2 = s.place("(b) request 1 after removal", pol, 1, fg=fg)
    s.expect("(b)", (r1.placed, r2.placed), (0, 1))
    s.expect("(b) plan frozen", pol.state()["cohorts"][0][0]["plan"], plan)
    s.expect(
        "(b) weights_changed", (r1.weights_changed, r2.weights_changed), (False, False)
    )


@fixture(
    "third_party_occupancy",
    "third_party",
    ["8a", "16.2", "18"],
    "Capacity 10, baseline 8, own 4, another owner 3: the residual is 1 and the "
    "view 5; request 2 adds exactly 1 for MaxFlow, Weighted and EB (exhausted) and "
    "never touches the other owner's record.",
)
def third_party(s: PSnap):
    c = Ctx(s, "one_link", parallel([10]))
    for kernel in (MF, W, EB):
        fg = c.ledger([8])
        pol = P.FlowPolicy(fg, 0, 1, co(kernel))
        s.place(f"{kernel} own 4", pol, 4, fg=fg)
        other, _ = c.manual(fg, 0, 1, [0], 3)
        r = s.place(f"{kernel} request 2", pol, 2, fg=fg)
        s.expect(f"{kernel}", (r.placed, r.termination), (1, "exhausted"))
        s.expect(
            f"{kernel} other untouched", [x.delivered for x in other.records()], [3]
        )


# ====================================================================== cohorts in units


@fixture(
    "eb_correlated_remainders",
    "units",
    ["16.5", "9.2"],
    "Round-robin remainders correlate across flows sharing a split: three parallel "
    "1-unit links. A Pinned EB flow on the three links places 3. Two Pinned bundles "
    "that become equal after pruning (the second carries a dead branch) keep two "
    "flows on the same split: T = 2 gives each flow its unit on the lowest-id link "
    "(2 > 1), so the cohort places 1 (exact and feasible; no start rotation, "
    "decided).",
)
def eb_correlated(s: PSnap):
    c = Ctx(
        s,
        "parallel_1_dead",
        (3, [(0, 1, 1, 1), (0, 1, 1, 1), (0, 1, 1, 1), (0, 2, 1, 1)]),
    )
    fg = c.ledger()
    r = s.place("one flow", P.FlowPolicy(fg, 0, 1, c.pin(EB, [[0, 1, 2]])), 3, fg=fg)
    s.expect("one flow", (r.placed, r.termination), (3, "complete"))
    fg = c.ledger()
    pol = P.FlowPolicy(fg, 0, 1, c.pin(EB, [[0, 1, 2], [0, 1, 2, 3]]))
    r = s.place("two flows on the same split", pol, 3, fg=fg)
    s.expect(
        "cohort", (r.placed, r.termination, c.shares(pol)), (1, "exhausted", [1, 0])
    )
    s.expect("routes equal after pruning", c.routes_of(pol), [[0, 1, 2], [0, 1, 2]])


@fixture(
    "eb_1_plus_2_vs_3",
    "units",
    ["16.4", "16.5"],
    "Three ample parallel links: EB place(1); place(2) equals place(3) exactly "
    "([1, 1, 1]) for one CostOnly flow and for a Pinned cohort of three single-link "
    "flows (round robin over the cohort in FlowId order).",
)
def eb_additivity(s: PSnap):
    c = Ctx(s, "parallel_10", parallel([10, 10, 10]))
    for mode in ("CostOnly", "Pinned"):
        finals = []
        for calls in ([1, 2], [3]):
            fg = c.ledger()
            cfg = co(EB) if mode == "CostOnly" else c.pin(EB, [[0], [1], [2]])
            pol = P.FlowPolicy(fg, 0, 1, cfg)
            got = [s.place(f"{mode} {v}", pol, v, fg=fg).placed for v in calls]
            s.expect(f"{mode} {calls}", got, calls)
            finals.append(c.alloc(fg))
        s.expect(f"{mode} 1 + 2 equals 3", finals[0], finals[1])
        s.expect(f"{mode} allocation", finals[1], [1, 1, 1])


@fixture(
    "weighted_1_plus_5_vs_6",
    "units",
    ["16.4"],
    "Three equal parallel capacity-2 links, CostOnly Weighted: place(1); place(5) "
    "equals place(6) ([2, 2, 2]).",
)
def weighted_additivity(s: PSnap):
    c = Ctx(s, "parallel_2", parallel([2, 2, 2]))
    finals = []
    for calls in ([1, 5], [6]):
        fg = c.ledger()
        pol = P.FlowPolicy(fg, 0, 1, co(W))
        got = [s.place(f"{calls} {v}", pol, v, fg=fg).placed for v in calls]
        s.expect(f"{calls}", got, calls)
        finals.append(c.alloc(fg))
    s.expect("1 + 5 equals 6", finals[0], finals[1])
    s.expect("allocation", finals[1], [2, 2, 2])


# ====================================================================== admission


@fixture(
    "admission_zero_and_no_route",
    "admission",
    ["13", "9.4", "8"],
    "Zero request: complete, nothing created, flags false. Overlapping endpoint "
    "sets (src == dst) are rejected at construction (rule 13). Zero capacity is "
    "not a missing route (contract r8 3.5, step 1): CostOnly over zero-capacity "
    "links has the source in its roster and creates its zero-delivery record on "
    "the static support (exhausted: MaxFlow carries 0 there). No route: "
    "CapacityAware on a zero baseline has a static roster but no residual "
    "route (no_route, no id consumed); Pinned with every bundle excluded "
    "(usable_bundles 0, no_route).",
)
def admission(s: PSnap):
    c = Ctx(s, "diamond_5", diamond(5, 5, 5, 5))
    fg = c.ledger()
    pol = P.FlowPolicy(fg, 0, 3, ca(MF))
    r = s.place("request 0", pol, 0)
    s.expect(
        "zero",
        (
            r.termination,
            r.placed,
            r.routes_changed,
            r.weights_changed,
            fg.next_flow_id(),
        ),
        ("complete", 0, False, False, 0),
    )
    s.raises("src == dst rejected", P.FlowPolicy, fg, 0, 0, ca(MF))
    s.raises("overlapping sets rejected", P.FlowPolicy, fg, [0, 1], [1, 3], co(MF))
    c0 = Ctx(s, "diamond_0", diamond(0, 0, 0, 0))
    fg0 = c0.ledger()
    p0 = P.FlowPolicy(fg0, 0, 3, co(MF))
    r = s.place("CostOnly over zero-capacity links", p0, 1, fg=fg0)
    # Both routes of the diamond cost 2, so the static support is all four links and the
    # MaxFlow route is that support; every link has capacity 0, so the flow carries 0.
    s.expect(
        "CostOnly zero capacity",
        (
            r.termination,
            r.placed,
            p0.roster,
            [(x.delivered, x.route_cost, len(x.bundle)) for x in fg0.records()],
        ),
        ("exhausted", 0, [0], [(0, 2, 4)]),
    )
    fgz = c.ledger([0, 0, 0, 0])
    pz = P.FlowPolicy(fgz, 0, 3, ca(EB, fc=2))
    r = s.place("CapacityAware on a zero residual", pz, 1, fg=fgz)
    s.expect(
        "CapacityAware no route", (r.termination, fgz.next_flow_id()), ("no_route", 0)
    )
    fgd = c.ledger()
    elig = c.g.eligibility(excluded_edges=c.ids([0, 2]))
    pd = P.FlowPolicy(fgd, 0, 3, c.pin(MF, [[0, 1], [2, 3]]), elig)
    r = s.place("Pinned all down", pd, 1, fg=fgd)
    s.expect("Pinned no route", (r.termination, pd.usable_bundles()), ("no_route", 0))


# ====================================================================== limits


@fixture(
    "limits_kmaxflowspercall",
    "limits",
    ["8b", "9.2", "9.4"],
    "kMaxFlowsPerCall is an engine constant (10,000); the fixture lowers the "
    "module's test hook to 3. Growth over five 1-unit links stops after 3 creations "
    "(limit_flows, placed 3); a Pinned policy with 5 bundles still creates 5 flows "
    "(the cap counts CapacityAware creations); two CapacityAware sources share one "
    "created counter: source 0 creates 3, source 1's creation is blocked by the cap "
    "(local limit_flows, nothing committed for it, no id consumed).",
)
def limits_created(s: PSnap):
    s.expect("engine constant", K_MAX_FLOWS_PER_CALL, 10000)
    saved = P.MAX_FLOWS_PER_CALL
    P.MAX_FLOWS_PER_CALL = 3
    try:
        c = Ctx(s, "parallel_1x5", parallel([1] * 5))
        fg = c.ledger()
        pol = P.FlowPolicy(fg, 0, 1, ca(MF, priority="Residual"))
        r = s.place("request 5", pol, 5)
        s.expect(
            "growth capped",
            (r.placed, r.termination, len(pol.flow_ids())),
            (3, "limit_flows", 3),
        )
        fg = c.ledger()
        pin = P.FlowPolicy(fg, 0, 1, c.pin(MF, [[i] for i in range(5)]))
        r = s.place("Pinned 5 bundles", pin, 5, fg=fg)
        s.expect("Pinned not capped", (r.placed, len(pin.flow_ids())), (5, 5))
        c2 = Ctx(s, "two_sources", (3, [(0, 2, 1, 1)] * 4 + [(1, 2, 1, 1)] * 4))
        fg = c2.ledger()
        two = P.FlowPolicy(fg, [0, 1], 2, ca(MF, priority="Residual"))
        r = s.place("two sources request 8", two, 8, fg=fg)
        s.expect("offers", r.trace["offers"], {0: 4, 1: 4})
        s.expect(
            "local outcomes",
            r.trace["local_outcomes"],
            {0: "limit_flows", 1: "limit_flows"},
        )
        s.expect(
            "result",
            (r.placed, r.termination, len(two.flow_ids())),
            (3, "limit_flows", 3),
        )
    finally:
        P.MAX_FLOWS_PER_CALL = saved


@fixture(
    "limits_diminishing",
    "limits",
    ["8b", "9.2", "9.4"],
    "q = 2^-12: one 50-capacity link (204,800 units) and twenty 0.01 links (41 "
    "units each). MaxFlow tunnels (Residual), Diminishing(window 8, fraction 1e-3), "
    "request 100 (409,600 units): growth adds 41 units per step; after 8 steps the "
    "window sum 328 < 409.6 stops growth (limit_diminishing). Without the guard all "
    "20 links are used (exhausted); +inf never arms it.",
)
def limits_diminishing(s: PSnap):
    q = 2.0**-12
    from units import to_units

    small = to_units(0.01, q)
    big = to_units(50, q)
    s.expect("units", (small, big), (41, 204800))
    c = Ctx(s, "fan", parallel([big] + [small] * 20), resolution=q)
    dim = P.Diminishing(8, 1e-3)
    for label, d, volume, steps, term in (
        ("window8", dim, to_units(100, q), 8, "limit_diminishing"),
        ("no guard", None, to_units(100, q), 20, "exhausted"),
        ("inf", dim, INF, 20, "exhausted"),
    ):
        fg = c.ledger()
        pol = P.FlowPolicy(fg, 0, 1, ca(MF, priority="Residual", diminishing=d))
        r = s.place(f"{label}", pol, volume, fg=fg)
        s.expect(
            f"{label}",
            (r.placed, r.termination, len(pol.flow_ids())),
            (big + steps * small, term, steps + 1),
        )


# ====================================================================== lifetime


@fixture(
    "lifetime_remove_orphan_reset",
    "lifetime",
    ["8", "18", "13"],
    "Section 8 with policies (replaces the invalidation fixture: there is no "
    "external removal): reset and copy are rejected while a policy is bound; "
    "remove() reverts and retires the ids; a later place opens a new cohort with a "
    "new id (routes_changed true); remove() on an empty cohort resets the epoch "
    "state; a destroyed policy's records stay as orphans until reset; retired ids "
    "throw; ids are never reused.",
)
def lifetime(s: PSnap):
    c = Ctx(s, "diamond_5", diamond(5, 5, 5, 5))
    fg = c.ledger()
    pol = P.FlowPolicy(fg, 0, 3, ca(MF))
    s.place("request 3", pol, 3)
    fid = pol.flow_ids()[0]
    s.raises("reset while bound", fg.reset)
    s.raises("copy while bound", fg.copy)
    pol.remove()
    s.expect("reverted", (c.alloc(fg), pol.placed_demand()), ([0, 0, 0, 0], 0))
    s.raises("retired id", fg.record, fid)
    r = s.place("restart", pol, 2)
    s.expect(
        "restart",
        (r.routes_changed, r.weights_changed, pol.flow_ids()[0][1]),
        (True, False, 1),
    )
    pol.remove()
    pol.remove()
    s.expect(
        "empty-cohort remove resets the epoch",
        (pol.requested, pol.infinite),
        (0, False),
    )
    other = P.FlowPolicy(fg, 0, 3, co(EB))
    s.place("other 2", other, 2)
    other.destroy()
    s.expect("orphan readable", [x.delivered for x in fg.records()], [2])
    s.raises("reset while pol bound", fg.reset)
    pol.destroy()
    fg.reset()
    s.expect("reset clears orphans", (fg.records(), fg.next_flow_id()), ([], 3))


# ====================================================================== Pinned


@fixture(
    "pinned_kernels_5_1",
    "pinned",
    ["8a", "19", "16.5"],
    "Pinned routes of capacity 5 and 1, request 8: MaxFlow 6 and Weighted 6 "
    "(sequential fill, exhausted, remaining 2); EqualBalanced 3 (round robin over "
    "the cohort: [2, 1]); Lossy 5 (offers [4, 4]: the 1-route carries 1 and drops 3 "
    "on its first link; complete).",
)
def pinned_kernels(s: PSnap):
    c = Ctx(s, "diamond_5_1", diamond(5, 5, 1, 1))
    for kernel, adm, placed, term in (
        (MF, LOSSLESS, 6, "exhausted"),
        (W, LOSSLESS, 6, "exhausted"),
        (EB, LOSSLESS, 3, "exhausted"),
        (EB, LOSSY, 5, "complete"),
    ):
        fg = c.ledger()
        pol = P.FlowPolicy(fg, 0, 3, c.pin(kernel, [[0, 1], [2, 3]], adm))
        r = s.place(f"{kernel}/{adm} request 8", pol, 8, fg=fg)
        s.expect(f"{kernel}/{adm}", (r.placed, r.termination), (placed, term))
        if adm == LOSSY:
            s.expect("Lossy drops", r.drops, [(c.ids([2])[0], 3)])


@fixture(
    "pinned_bundle_rules",
    "pinned",
    ["13", "12", "9.1"],
    "Pinned rules: an empty list, duplicate bundles (pre-prune), a bundle without a "
    "src -> dst walk and a source set are rejected; bundles that become equal after "
    "pruning keep independent flows; a multi-walk DAG bundle with mixed costs "
    "reports its minimum-cost walk as route_cost (3, not 4).",
)
def pinned_rules(s: PSnap):
    c = Ctx(
        s,
        "mixed",
        (4, [(0, 1, 5, 1), (1, 3, 5, 2), (0, 2, 5, 1), (2, 3, 5, 3), (0, 3, 5, 9)]),
    )
    fg = c.ledger()
    s.raises(
        "empty list",
        P.FlowPolicy,
        fg,
        0,
        3,
        P.FlowPolicyConfig(P.Pinned(()), MF, LOSSLESS, search.Selection()),
    )
    s.raises("duplicate", P.FlowPolicy, fg, 0, 3, c.pin(MF, [[0, 1], [0, 1]]))
    s.raises("malformed", P.FlowPolicy, fg, 0, 3, c.pin(MF, [[0, 2]]))
    s.raises("source set", P.FlowPolicy, fg, [0, 1], 3, c.pin(MF, [[0, 1]]))
    elig = c.g.eligibility(excluded_edges=c.ids([4]))
    pol = P.FlowPolicy(fg, 0, 3, c.pin(MF, [[0, 1], [0, 1, 4]]), elig)
    s.expect("two usable bundles", pol.usable_bundles(), 2)
    s.place("request 1", pol, 1)
    s.expect(
        "independent flows on equal pruned routes", c.routes_of(pol), [[0, 1], [0, 1]]
    )
    fg2 = c.ledger()
    multi = P.FlowPolicy(fg2, 0, 3, c.pin(MF, [[0, 1, 2, 3]]))
    s.place("multi-walk bundle", multi, 7, fg=fg2)
    s.expect("route_cost is the min-cost walk", multi.records()[0].route_cost, 3)
    s.expect("cheapest walk fills first", c.alloc(fg2), [5, 5, 2, 2, 0])


# ====================================================================== flags


@fixture(
    "flag_table",
    "flags",
    ["14", "9.4"],
    "The authoritative flag table of revision 20 (flags accumulate committed "
    "events): creation (T, F), also at zero delivery; no-op (F, F); saturated call "
    "(F, F); adopted Weighted reroute that admits zero (T, T; review-11 C4); "
    "rejected growth only (F, F); accepted growth (T, F); EB routes [3, 2, 1] "
    "holding 3, request 3: growth onto the second route accepted (gain 2), onto the "
    "third rejected: routes_changed stays true; restart after remove() (T, F).",
)
def flag_table(s: PSnap):
    c = Ctx(s, "diamond_5", diamond(5, 5, 5, 5))
    fg = c.ledger()
    pol = P.FlowPolicy(fg, 0, 3, ca(MF, mx=1))
    rows = {}
    rows["creation"] = s.place("creation", pol, 5)
    rows["no-op"] = s.place("no-op", pol, 0)
    rows["saturated without growth candidate change"] = s.place("saturated", pol, 1)
    pol.remove()
    rows["restart"] = s.place("restart", pol, 1)
    cw = Ctx(s, "parallel_10_20", parallel([10, 20]))
    fgw = cw.ledger()
    blocking, _ = cw.manual(fgw, 0, 1, [1], 20, "blocking owner")
    wsel = tunnel("Capacity")
    wpol = P.FlowPolicy(
        fgw,
        0,
        1,
        P.FlowPolicyConfig(
            P.CapacityAware(max_flow_count=1, reoptimize=True), W, LOSSLESS, wsel
        ),
    )
    s.place("Weighted place(4)", wpol, 4, fg=fgw)
    cw.manual(fgw, 0, 1, [0], 6, "owner fills e0")
    blocking.remove()
    cw.manual(fgw, 0, 1, [1], 16, "owner takes 16 of e1")
    rows["adopted Weighted reroute at zero"] = s.place(
        "Weighted reroute", wpol, 1, fg=fgw
    )
    s.expect(
        "reroute admitted zero", rows["adopted Weighted reroute at zero"].placed, 0
    )
    ch = Ctx(s, "harmful", diamond(10, 10, 1, 1, c0=1, c1=2))
    fgh = ch.ledger()
    ph = P.FlowPolicy(fgh, 0, 3, ca(EB, fc=1, mx=2))
    s.place("P = 10", ph, 10, fg=fgh)
    rows["rejected growth only"] = s.place("rejected growth", ph, 1, fg=fgh)
    cg = Ctx(s, "diamond_5b", diamond(5, 5, 5, 5))
    fgg = cg.ledger()
    pg = P.FlowPolicy(fgg, 0, 3, ca(EB, fc=1, mx=2))
    s.place("first 5", pg, 5, fg=fgg)
    rows["accepted growth"] = s.place("accepted growth", pg, 5, fg=fgg)
    c3 = Ctx(s, "routes_3_2_1", routes([(3, 1), (2, 2), (1, 3)]))
    fg3 = c3.ledger()
    p3 = P.FlowPolicy(fg3, 0, 4, ca(EB, fc=1, mx=3))
    s.place("hold 3", p3, 3, fg=fg3)
    r = s.place("request 3", p3, 3, fg=fg3)
    decisions = [
        x.get("decision", x.get("stop")) for x in r.trace["sources"][0]["stage4"]
    ]
    s.expect(
        "growth decisions",
        [d.split(":")[0] for d in decisions],
        ["accepted", "rejected"],
    )
    s.expect(
        "[3, 2, 1] witness",
        (r.placed, r.termination, c3.shares(p3)),
        (2, "exhausted", [3, 2]),
    )
    rows["accepted then rejected growth"] = r
    want = {
        "creation": (True, False),
        "no-op": (False, False),
        "saturated without growth candidate change": (False, False),
        "restart": (True, False),
        "adopted Weighted reroute at zero": (True, True),
        "rejected growth only": (False, False),
        "accepted growth": (True, False),
        "accepted then rejected growth": (True, False),
    }
    for name, (rc, wc) in want.items():
        res = rows[name]
        s.expect(f"flags: {name}", (res.routes_changed, res.weights_changed), (rc, wc))


@fixture(
    "flag_zero_delivery_creation",
    "flags",
    ["9.4", "8"],
    "Revision-20 witness: a first place(1) on a static route whose unit another "
    "owner holds creates a zero-delivery record: routes_changed true, termination "
    "exhausted, the record stays (CostOnly and Pinned); a second zero-delivery call "
    "keeps the same ids and sets no flag (revision 21).",
)
def flag_zero_delivery(s: PSnap):
    c = Ctx(s, "one_unit", parallel([1]))
    for mode in ("CostOnly", "Pinned"):
        fg = c.ledger()
        c.manual(fg, 0, 1, [0], 1)
        cfg = co(EB) if mode == "CostOnly" else c.pin(EB, [[0]])
        pol = P.FlowPolicy(fg, 0, 1, cfg)
        r = s.place(f"{mode} place(1)", pol, 1, fg=fg)
        s.expect(
            f"{mode}",
            (
                r.placed,
                r.termination,
                r.routes_changed,
                [x.delivered for x in pol.records()],
            ),
            (0, "exhausted", True, [0]),
        )
        ids = pol.flow_ids()
        r = s.place(f"{mode} place(1) again", pol, 1, fg=fg)
        s.expect(
            f"{mode} records stay pinned (revision 21)",
            (r.placed, r.routes_changed, pol.flow_ids()),
            (0, False, ids),
        )


# ====================================================================== source sets


def two_direct(caps, costs=(1, 1)):
    """Sources 0 and 1, each with one link to the destination 2."""
    return 3, [(0, 2, caps[0], costs[0]), (1, 2, caps[1], costs[1])]


@fixture(
    "sources_lossy_per_source",
    "sources",
    ["24", "26", "19", "9.3"],
    "Lossy per source (replaces the zero-cost forward-versus-reverse witness): "
    "sources 0 and 1 reach 3 through X (0->X 4, 1->X 4, X->3 3) and source 1 also "
    "through Y (1->Y 4, Y->3 4, cost 2 each). Request 8 offers [4, 4]; source 0 "
    "places 3 and drops 1 on X->3; source 1's DAG is its own shortest walks (1->X->3 "
    "only, cost 2 < 4) and the full X->3 drops its 4. Delivered 3, complete (both "
    "offers injected).",
)
def sources_lossy(s: PSnap):
    edges = [(0, 3, 4, 1), (1, 3, 4, 1), (3, 5, 3, 1), (1, 4, 4, 2), (4, 5, 4, 2)]
    c = Ctx(s, "per_source", (6, edges))
    fg = c.ledger()
    pol = P.FlowPolicy(fg, [0, 1], 5, co(EB, LOSSY))
    r = s.place("request 8", pol, 8)
    s.expect("offers", r.trace["offers"], {0: 4, 1: 4})
    s.expect("result", (r.placed, r.termination), (3, "complete"))
    recs = {x.sources[0]: x for x in pol.records()}
    s.expect(
        "source 0",
        (recs[0].delivered, dict_to_input(c.g, recs[0].drops)),
        (3, [0, 0, 1, 0, 0]),
    )
    s.expect(
        "source 1",
        (recs[1].delivered, dict_to_input(c.g, recs[1].drops)),
        (0, [0, 0, 4, 0, 0]),
    )


@fixture(
    "sources_attribution",
    "sources",
    ["24", "26", "12", "9.3"],
    "Sources A (A->T 1, cost 1) and B (B->M 10, M->T 2, cost 1 each). Lossy "
    "request 20 (offers [10, 10]): A delivers 1 and drops 9 on A->T; B delivers 2 "
    "and drops 8 on M->T; placed 3, complete. MaxFlow, Weighted and EB per source: "
    "1 + 2 = 3, exhausted (the old pseudo-root EB gave 2). Route costs per record: "
    "1 and 2.",
)
def sources_attribution(s: PSnap):
    c = Ctx(s, "attribution", (4, [(0, 3, 1, 1), (1, 2, 10, 1), (2, 3, 2, 1)]))
    fg = c.ledger()
    pol = P.FlowPolicy(fg, [0, 1], 3, co(EB, LOSSY))
    r = s.place("Lossy 20", pol, 20)
    recs = {x.sources[0]: x for x in pol.records()}
    s.expect("Lossy", (r.placed, r.termination), (3, "complete"))
    s.expect(
        "A",
        (
            recs[0].offered,
            recs[0].delivered,
            sum(recs[0].drops.values()),
            recs[0].route_cost,
        ),
        (10, 1, 9, 1),
    )
    s.expect(
        "B",
        (
            recs[1].offered,
            recs[1].delivered,
            sum(recs[1].drops.values()),
            recs[1].route_cost,
        ),
        (10, 2, 8, 2),
    )
    for kernel in (MF, W, EB):
        fg = c.ledger()
        pol = P.FlowPolicy(fg, [0, 1], 3, co(kernel))
        r = s.place(f"{kernel} 20", pol, 20, fg=fg)
        s.expect(
            f"{kernel}",
            (r.placed, r.termination, [x.delivered for x in pol.records()]),
            (3, "exhausted", [1, 2]),
        )


@fixture(
    "sources_eb_per_source",
    "sources",
    ["26", "9.3"],
    "EqualBalanced per source (independent routers, rule 26): A (A->T 4, cost 1) "
    "and B (B->T 10, cost 2). Request 20 offers [10, 10]: A places 4, B places 10, "
    "total 14 (the old pseudo-root flow with joint equalisation capped both at 4: "
    "8). +inf also 14.",
)
def sources_eb(s: PSnap):
    c = Ctx(s, "two_direct", two_direct((4, 10), (1, 2)))
    for v, term in ((20, "exhausted"), (INF, "exhausted")):
        fg = c.ledger()
        pol = P.FlowPolicy(fg, [0, 1], 2, co(EB))
        r = s.place(f"EB {v}", pol, v, fg=fg)
        s.expect(
            f"EB {v}",
            (r.placed, r.termination, [x.delivered for x in pol.records()]),
            (14, term, [4, 10]),
        )


@fixture(
    "sources_augmentation_edges",
    "sources",
    ["21", "9.3"],
    "User augmentation links are ordinary positive-cost edges (no custom fanout): "
    "source A with a real link A->T (1, cost 2) and an augmentation path A->X->T "
    "(9 each, cost 1): the CostOnly DAG holds both cost-2 walks; MaxFlow places 10; "
    "EB splits [ceil, floor] at A over A->X (lower id) and A->T: 3.",
)
def sources_augmentation(s: PSnap):
    c = Ctx(s, "augmentation", (3, [(0, 2, 1, 2), (0, 1, 9, 1), (1, 2, 9, 1)]))
    for kernel, placed in ((MF, 10), (EB, 3)):
        fg = c.ledger()
        r = s.place(f"{kernel} +inf", P.FlowPolicy(fg, 0, 2, co(kernel)), INF, fg=fg)
        s.expect(f"{kernel}", r.placed, placed)


@fixture(
    "sources_inf_every_kernel",
    "sources",
    ["24", "15", "9.3"],
    "+inf per source on the attribution graph: MaxFlow 3, Weighted 3, EB 3 "
    "(exhausted); Lossy injects every source member to its residual (A 1, B 10): "
    "delivered 3, offered 11, drops 8, exhausted.",
)
def sources_inf(s: PSnap):
    c = Ctx(s, "attribution", (4, [(0, 3, 1, 1), (1, 2, 10, 1), (2, 3, 2, 1)]))
    for kernel in (MF, W, EB):
        fg = c.ledger()
        r = s.place(
            f"{kernel} +inf", P.FlowPolicy(fg, [0, 1], 3, co(kernel)), INF, fg=fg
        )
        s.expect(f"{kernel}", (r.placed, r.termination), (3, "exhausted"))
    fg = c.ledger()
    pol = P.FlowPolicy(fg, [0, 1], 3, co(EB, LOSSY))
    r = s.place("Lossy +inf", pol, INF, fg=fg)
    s.expect(
        "Lossy",
        (
            r.placed,
            r.termination,
            sum(x.offered for x in pol.records()),
            sum(sum(x.drops.values()) for x in pol.records()),
        ),
        (3, "exhausted", 11, 8),
    )


@fixture(
    "sources_limits_and_validation",
    "sources",
    ["13", "9.3"],
    "Endpoint-set validation: Pinned with a source set, an empty source set and "
    "overlapping sets are rejected; an excluded source is skipped from the roster "
    "(roster [1]); a source that does not reach a destination statically is not in "
    "the roster; three sources of a CostOnly policy create three records.",
)
def sources_validation(s: PSnap):
    c = Ctx(
        s,
        "three_sources",
        (5, [(0, 3, 5, 1), (1, 3, 5, 1), (2, 3, 5, 1), (4, 0, 1, 1)]),
    )
    fg = c.ledger()
    s.raises("Pinned source set", P.FlowPolicy, fg, [0, 1], 3, c.pin(MF, [[0]]))
    s.raises("empty source set", P.FlowPolicy, fg, [], 3, co(MF))
    s.raises("overlap", P.FlowPolicy, fg, [0, 3], 3, co(MF))
    elig = c.g.eligibility(excluded_nodes=[0])
    pe = P.FlowPolicy(fg, [0, 1], 3, co(MF), elig)
    s.expect("excluded source skipped", pe.roster, [1])
    fg2 = c.ledger()
    p3 = P.FlowPolicy(fg2, [0, 1, 2], 3, co(MF))
    r = s.place("three sources request 3", p3, 3, fg=fg2)
    s.expect("three records", (r.placed, len(p3.records())), (3, 3))
    c2 = Ctx(s, "unreachable_source", (3, [(0, 2, 5, 1), (2, 1, 5, 1)]))
    p = P.FlowPolicy(c2.ledger(), [0, 1], 2, co(MF))
    s.expect("roster excludes a source without a static route", p.roster, [0])


@fixture(
    "sources_maxflow_vs_global",
    "sources",
    ["24", "27", "9.3"],
    "Independent sources versus the controller: 0->2 (cost 2), 1->2 (cost 1), "
    "2->3 (cost 1) and 0->3 (cost 5), unit capacities. CostOnly MaxFlow per source: "
    "source 0's DAG is its shortest walk 0->2->3 (cost 3) and takes the shared "
    "edge; source 1 then finds 2->3 full: +inf places 1. Global's unique optimum "
    "(cost 7 against 8) routes source 0 over 0->3 and source 1 over 1->2->3: 2 (the "
    "per-source DAGs cannot cancel across sources).",
)
def sources_vs_global(s: PSnap):
    c = Ctx(
        s, "shared_edge", (4, [(0, 2, 1, 2), (1, 2, 1, 1), (2, 3, 1, 1), (0, 3, 1, 5)])
    )
    fg = c.ledger()
    r = s.place(
        "CostOnly MaxFlow +inf", P.FlowPolicy(fg, [0, 1], 3, co(MF)), INF, fg=fg
    )
    s.expect("per source", r.placed, 1)
    fg = c.ledger()
    pg = P.FlowPolicy(fg, [0, 1], 3, glob())
    r = s.place("Global +inf", pg, INF, fg=fg)
    s.expect("Global", (r.placed, c.alloc(fg)), (2, [0, 1, 1, 1]))
    s.expect(
        "Global per-source records",
        [(x.sources[0], x.delivered) for x in pg.records()],
        [(0, 1), (1, 1)],
    )


@fixture(
    "sources_split_cumulative",
    "sources",
    ["24", "16.7", "9.3"],
    "The source split is over the cumulative request (16.7): sources A, B, C with "
    "direct links of 10. place(1); place(1) offers [1, 0, 0] then [0, 1, 0] (round "
    "robin of 2 is [1, 1, 0]), identical to place(2); a per-call split would offer "
    "[1, 0, 0] twice (control). Lossy follows the same cumulative split: two calls "
    "of 1 over two sources offer [1, 1].",
)
def sources_split(s: PSnap):
    c = Ctx(s, "three_direct", (4, [(0, 3, 10, 1), (1, 3, 10, 1), (2, 3, 10, 1)]))
    fg = c.ledger()
    pol = P.FlowPolicy(fg, [0, 1, 2], 3, co(MF))
    a = s.place("call 1", pol, 1)
    b = s.place("call 2", pol, 1)
    s.expect(
        "offers",
        (a.trace["offers"], b.trace["offers"]),
        ({0: 1, 1: 0, 2: 0}, {0: 0, 1: 1, 2: 0}),
    )
    fg2 = c.ledger()
    s.place("one call of 2", P.FlowPolicy(fg2, [0, 1, 2], 3, co(MF)), 2, fg=fg2)
    s.expect("1 + 1 equals 2", c.alloc(fg), c.alloc(fg2))
    s.expect("per-call split control differs", [1, 0, 0] + [1, 0, 0] != [1, 1, 0], True)
    c2 = Ctx(s, "two_direct", two_direct((5, 5)))
    fg3 = c2.ledger()
    lp = P.FlowPolicy(fg3, [0, 1], 2, co(EB, LOSSY))
    x = s.place("Lossy 1", lp, 1, fg=fg3)
    y = s.place("Lossy 1", lp, 1, fg=fg3)
    s.expect(
        "Lossy cumulative offers",
        (x.trace["offers"], y.trace["offers"]),
        ({0: 1, 1: 0}, {0: 0, 1: 1}),
    )


@fixture(
    "sources_1_9_request_4",
    "sources",
    ["24", "9.3", "9.4"],
    "Review-13 source budget witness: sources with capacities [1, 9], request 4 "
    "offers [2, 2]; deliveries [1, 2]; unserved remainders are not redistributed: "
    "placed 3, remaining 1; local outcomes exhausted and complete; call exhausted.",
)
def sources_1_9(s: PSnap):
    c = Ctx(s, "two_direct_1_9", two_direct((1, 9)))
    for kernel in (MF, W, EB):
        fg = c.ledger()
        pol = P.FlowPolicy(fg, [0, 1], 2, co(kernel))
        r = s.place(f"{kernel} request 4", pol, 4, fg=fg)
        s.expect(
            f"{kernel}",
            (
                r.trace["offers"],
                [x.delivered for x in pol.records()],
                r.placed,
                r.remaining,
                r.trace["local_outcomes"],
                r.termination,
            ),
            ({0: 2, 1: 2}, [1, 2], 3, 1, {0: "exhausted", 1: "complete"}, "exhausted"),
        )


@fixture(
    "sources_partition_1_1_vs_2_0",
    "sources",
    ["24", "16.4", "9.2"],
    "Source-set allocation is partition-dependent (review-14 R3; 16.4 excludes "
    "competing sources): 0->2 and 1->2 (2 each) share 2->3 (2). place(2); place(2) "
    "gives deliveries [1, 1] (the second call finds 2->3 full); place(4) gives "
    "[2, 0] (source 0 is placed first in ascending NodeId). Same total, different "
    "records: independent sources, no reassignment.",
)
def sources_partition(s: PSnap):
    c = Ctx(s, "shared_2", (4, [(0, 2, 2, 1), (1, 2, 2, 1), (2, 3, 2, 1)]))
    out = []
    for calls in ([2, 2], [4]):
        fg = c.ledger()
        pol = P.FlowPolicy(fg, [0, 1], 3, co(MF))
        for v in calls:
            s.place(f"{calls} call {v}", pol, v, fg=fg)
        out.append([x.delivered for x in pol.records()])
    s.expect("partition dependence", out, [[1, 1], [2, 0]])


@fixture(
    "eb_two_flows_per_source_reservation",
    "sources",
    ["9.2", "9.3", "8b"],
    "Per-source binding of the EB reservation target (review-16 R4): two sources, "
    "EB CapacityAware flow_count 2 each, request 4: offers [2, 2], each source's "
    "target (P_s + v_s) div flow_count = 1 per flow (not the policy-wide 2).",
)
def eb_reservation_binding(s: PSnap):
    c = Ctx(
        s,
        "two_sources_two_links",
        (3, [(0, 2, 5, 1), (0, 2, 5, 1), (1, 2, 5, 1), (1, 2, 5, 1)]),
    )
    fg = c.ledger()
    pol = P.FlowPolicy(fg, [0, 1], 2, ca(EB, fc=2, mx=2, sel=tunnel("Residual")))
    r = s.place("request 4", pol, 4)
    targets = {src: tr["stage1"]["eb_target"] for src, tr in r.trace["sources"].items()}
    reserved = {
        src: [x["reserved"] for x in tr["stage1"]["reservations"]]
        for src, tr in r.trace["sources"].items()
    }
    s.expect("targets", targets, {0: 1, 1: 1})
    s.expect("reservations", reserved, {0: [1, 1], 1: [1, 1]})
    s.expect("result", (r.placed, r.termination), (4, "complete"))


# ====================================================================== termination


@fixture(
    "termination_reduction",
    "termination",
    ["9.4", "14"],
    "Local outcomes and the reduction limit_flows > limit_diminishing > exhausted > "
    "no_route over unserved sources (review-16 R2): one limit_flows and one "
    "exhausted give limit_flows; two sources without routes give no_route; one "
    "exhausted and one no_route give exhausted; a fixed-count cohort full is "
    "exhausted; a Lossy finite offer consumed with drops is complete; v = 0 is "
    "complete; an empty static roster (no walk from the source; a zero-capacity "
    "walk would count) is no_route.",
)
def termination(s: PSnap):
    # limit_flows + exhausted: source 0 grows over two 1-unit links with max 3 >
    # flow_count 1 (limit_flows needs the guard to fire: max 2 reached with
    # remaining offer); source 1 has one 1-unit link and no growth candidate.
    c = Ctx(s, "mixed", (3, [(0, 2, 1, 1), (0, 2, 1, 1), (0, 2, 1, 1), (1, 2, 1, 1)]))
    fg = c.ledger()
    pol = P.FlowPolicy(fg, [0, 1], 2, ca(MF, priority="Residual", fc=1, mx=2))
    r = s.place("request 6 (offers [3, 3])", pol, 6)
    s.expect(
        "limit_flows + exhausted",
        (r.trace["local_outcomes"], r.termination),
        ({0: "limit_flows", 1: "exhausted"}, "limit_flows"),
    )
    cz = Ctx(s, "no_routes", (3, [(0, 2, 1, 1), (1, 2, 1, 1)]))
    fgz = cz.ledger([0, 0])
    pz = P.FlowPolicy(fgz, [0, 1], 2, ca(MF))
    r = s.place("two sources without routes", pz, 2, fg=fgz)
    s.expect(
        "no_route + no_route",
        (r.trace["local_outcomes"], r.termination),
        ({0: "no_route", 1: "no_route"}, "no_route"),
    )
    fge = cz.ledger([1, 0])
    pe = P.FlowPolicy(fge, [0, 1], 2, ca(MF, mx=1))
    r = s.place("exhausted + no_route", pe, 4, fg=fge)
    s.expect(
        "exhausted + no_route",
        (r.trace["local_outcomes"], r.termination),
        ({0: "exhausted", 1: "no_route"}, "exhausted"),
    )
    cf = Ctx(s, "fixed", parallel([1]))
    fgf = cf.ledger()
    r = s.place("fixed count", P.FlowPolicy(fgf, 0, 1, ca(MF, fc=1, mx=1)), 2, fg=fgf)
    s.expect("fixed-count cohort full", (r.placed, r.termination), (1, "exhausted"))
    cl = Ctx(s, "lossy", parallel([1, 10]))
    fgl = cl.ledger()
    r = s.place(
        "Lossy finite with drops", P.FlowPolicy(fgl, 0, 1, co(EB, LOSSY)), 4, fg=fgl
    )
    s.expect(
        "Lossy complete with drops",
        (r.placed, r.termination, r.remaining),
        (3, "complete", 1),
    )
    r = s.place("v = 0", P.FlowPolicy(cl.ledger(), 0, 1, co(EB)), 0)
    s.expect("v = 0", r.termination, "complete")
    c0 = Ctx(s, "no_roster", (2, [(1, 0, 1, 1)]))
    p0 = P.FlowPolicy(c0.ledger(), 0, 1, co(EB))
    r = s.place("empty roster", p0, 1)
    s.expect("empty roster", (r.termination, r.local_outcomes), ("no_route", []))
    s.place("empty roster +inf", p0, INF)
    s.raises("finite after +inf with an empty roster", p0.place, 1)


@fixture(
    "eb_resplit_3_3_to_3_2_2",
    "growth",
    ["16.3", "9.2"],
    "Review-16 R7: CapacityAware EB, flow_count 2, max 3, three disjoint routes of "
    "capacities [3, 3, 2] and increasing cost. Request 6 reserves 3 on each of the "
    "first two ([3, 3]). Request 1: stage 3 cannot grow the pair; growth adds the "
    "third route and the cohort total 7 is resplit [3, 2, 2]: the second flow loses "
    "1 while the cohort and placed_demand() gain 1 (individual T_i and edge deltas "
    "may decrease; the cohort total never does).",
)
def eb_resplit(s: PSnap):
    c = Ctx(s, "routes_3_3_2", routes([(3, 1), (3, 2), (2, 3)]))
    fg = c.ledger()
    pol = P.FlowPolicy(fg, 0, 4, ca(EB, fc=2, mx=3))
    r = s.place("request 6", pol, 6)
    s.expect("first", (r.placed, c.shares(pol)), (6, [3, 3]))
    before = c.alloc(fg)
    r = s.place("request 1", pol, 1)
    s.expect(
        "resplit",
        (r.placed, r.termination, c.shares(pol), pol.placed_demand()),
        (1, "complete", [3, 2, 2], 7),
    )
    delta = [a - b for a, b in zip(c.alloc(fg), before, strict=True)]
    s.expect("edge deltas (signed)", delta, [0, 0, -1, -1, 2, 2])


# ====================================================================== Global


@fixture(
    "global_chain_plus_bypasses",
    "global",
    ["27", "25", "16.4"],
    "Global is a replacement kernel (review-13 D5): chain 0-1-2-3 (cost 1) with "
    "bypasses 0->2 and 1->3 (cost 3). place(1) takes the chain (cost 3); the next "
    "place(1) re-solves the cumulative value 2 on W + own: both bypasses, the middle "
    "edge released, placed 1 (total 2, not 1), cost_distribution [(3, 1), (5, 1)], "
    "routes_changed true (support changed). place(2) on a fresh ledger gives the "
    "same allocation and curve (additive on a fixed view).",
)
def global_chain(s: PSnap):
    c = Ctx(
        s,
        "chain",
        (4, [(0, 1, 1, 1), (1, 2, 1, 1), (2, 3, 1, 1), (0, 2, 1, 3), (1, 3, 1, 3)]),
    )
    fg = c.ledger()
    pol = P.FlowPolicy(fg, 0, 3, glob())
    a = s.place("place(1)", pol, 1)
    s.expect(
        "first",
        (a.placed, a.cost_distribution, c.alloc(fg), a.routes_changed),
        (1, [(3, 1)], [1, 1, 1, 0, 0], True),
    )
    b = s.place("place(1)", pol, 1)
    s.expect(
        "second",
        (b.placed, b.cost_distribution, c.alloc(fg), b.routes_changed),
        (1, [(3, 1), (5, 1)], [1, 0, 1, 1, 1], True),
    )
    fg2 = c.ledger()
    p2 = P.FlowPolicy(fg2, 0, 3, glob())
    r = s.place("place(2)", p2, 2, fg=fg2)
    s.expect(
        "additive on a fixed view",
        (r.cost_distribution, c.alloc(fg2)),
        (b.cost_distribution, c.alloc(fg)),
    )
    sol = {k: v for k, v in b.trace["flow_core"].items() if k != "supply"}
    s.expect(
        "maximal at cap 2: min_cut present",
        sorted(sol["min_cut"]),
        sorted(c.ids([0, 3])),
    )
    s.check(
        "cumulative solve certified",
        not certify.check_solve(c.g, list(c.g.capacity), [0], [3], 2, sol, tier="a"),
    )


@fixture(
    "global_capped_underfill_control",
    "global",
    ["2.2", "27"],
    "Attained value at the policy level: Global place(2) on one 3-unit edge places "
    "2 (min_cut None: a residual path is left; requested with bottleneck outputs), "
    "the cumulative solve passes the certificate; a mutated producer returning 1 "
    "under cap 2 is rejected; place(+inf) then places the last unit and min_cut is "
    "[0].",
)
def global_underfill(s: PSnap):
    c = Ctx(s, "one_edge", (2, [(0, 1, 3, 7)]))
    fg = c.ledger()
    pol = P.FlowPolicy(fg, 0, 1, glob())
    r = s.place("place(2)", pol, 2, outputs=P.Outputs(bottleneck=True))
    s.expect("capped", (r.placed, r.termination, r.min_cut), (2, "complete", None))
    payload = {k: v for k, v in r.trace["flow_core"].items() if k != "supply"}
    s.check(
        "certified", not certify.check_solve(c.g, [3], [0], [1], 2, payload, tier="a")
    )
    bad = dict(payload, placed=1, edge_flows=[1], cost_distribution=[[7, 1]])
    s.check(
        "underfill rejected",
        bool(certify.check_solve(c.g, [3], [0], [1], 2, bad, tier="a")),
    )
    r = s.place("place(+inf)", pol, INF, outputs=P.Outputs(bottleneck=True))
    s.expect("maximal", (r.placed, r.termination, r.min_cut), (1, "exhausted", [0]))


@fixture(
    "global_zero_growth_support_change",
    "global",
    ["9.4", "27", "20"],
    "Review-16 R1: edges 0->1 (cost 1, cap 1), 0->1 (cost 3, cap 2), 1->2 (cap 2). "
    "Owner X holds the cheap link and one unit of 1->2; Global place(1) takes the "
    "costly link (cost 4). X leaves; owner Y takes the costly link and the rest of "
    "1->2. Global place(1): the view admits no growth (placed 0, exhausted) but the "
    "re-solve moves the unit onto the released cheap link (cost 2): routes_changed "
    "true at placed 0; Y untouched; the next place(1) changes nothing (flags false).",
)
def global_support(s: PSnap):
    c = Ctx(s, "release", (3, [(0, 1, 1, 1), (0, 1, 2, 3), (1, 2, 2, 1)]))
    fg = c.ledger()
    x, _ = c.manual(fg, 0, 2, [0, 2], 1, "X")
    pol = P.FlowPolicy(fg, 0, 2, glob())
    r = s.place("Global place(1)", pol, 1)
    s.expect("costly", (r.placed, c.own(fg, pol)), (1, [0, 1, 1]))
    x.remove()
    y, _ = c.manual(fg, 0, 2, [1, 2], 1, "Y")
    r = s.place("Global place(1) at zero growth", pol, 1)
    s.expect(
        "moved",
        (r.placed, r.termination, r.routes_changed, c.own(fg, pol)),
        (0, "exhausted", True, [1, 0, 1]),
    )
    s.expect("Y untouched", c.own(fg, y), [0, 1, 1])
    r = s.place("Global place(1) unchanged", pol, 1)
    s.expect(
        "unchanged", (r.placed, r.routes_changed, r.weights_changed), (0, False, False)
    )


@fixture(
    "global_epochs",
    "global",
    ["8", "27", "9.1"],
    "Epochs for Global (review-16, revision 20): sources 0 (cost 3) and 1 (cost 1), "
    "each directly connected to the destination with capacity 2. place(0) creates "
    "nothing; place(1) opens the cohort with records for both sources (ids [0, 1]) "
    "and supply [0, 1] (no round robin: the controller picks the cheaper source); "
    "remove(); place(1) again supplies [0, 1] under ids [2, 3]. place(+inf) gives 4; "
    "a finite request and place(0) after +inf are rejected; remove(); place(1) is "
    "accepted.",
)
def global_epochs(s: PSnap):
    c = Ctx(s, "two_costs", two_direct((2, 2), (3, 1)))
    fg = c.ledger()
    pol = P.FlowPolicy(fg, [0, 1], 2, glob())
    r = s.place("place(0)", pol, 0)
    s.expect("no creation", (r.termination, fg.next_flow_id()), ("complete", 0))
    r = s.place("place(1)", pol, 1)
    s.expect(
        "epoch 1",
        (
            [x.flow_id[1] for x in pol.records()],
            [x.delivered for x in pol.records()],
            r.routes_changed,
        ),
        ([0, 1], [0, 1], True),
    )
    pol.remove()
    r = s.place("place(1) after remove", pol, 1)
    s.expect(
        "epoch 2",
        (
            [x.flow_id[1] for x in pol.records()],
            [x.delivered for x in pol.records()],
            r.routes_changed,
        ),
        ([2, 3], [0, 1], True),
    )
    r = s.place("place(+inf)", pol, INF)
    s.expect(
        "+inf", (r.placed, pol.placed_demand(), r.termination), (3, 4, "exhausted")
    )
    s.raises("finite after +inf", pol.place, 1)
    s.raises("place(0) after +inf", pol.place, 0)
    pol.remove()
    r = s.place("place(1) in a new epoch", pol, 1)
    s.expect("accepted", r.placed, 1)


@fixture(
    "epochs_place_remove_place",
    "lifetime",
    ["8", "9.3", "24"],
    "Section 8's epoch law for CostOnly over two sources: place(1) offers [1, 0] "
    "and creates both records (ids [0, 1], source 1 at zero delivery); remove(); "
    "place(1) offers [1, 0] again (requested reset) with ids [2, 3]; "
    "place(+inf); remove(); place(1) is accepted; remove() on an empty cohort "
    "resets the epoch state.",
)
def epochs(s: PSnap):
    c = Ctx(s, "two_direct", two_direct((5, 5)))
    fg = c.ledger()
    pol = P.FlowPolicy(fg, [0, 1], 2, co(MF))
    r = s.place("place(1)", pol, 1)
    s.expect(
        "epoch 1",
        (
            r.trace["offers"],
            [x.flow_id[1] for x in pol.records()],
            [x.delivered for x in pol.records()],
        ),
        ({0: 1, 1: 0}, [0, 1], [1, 0]),
    )
    pol.remove()
    r = s.place("place(1) after remove", pol, 1)
    s.expect(
        "epoch 2",
        (r.trace["offers"], [x.flow_id[1] for x in pol.records()]),
        ({0: 1, 1: 0}, [2, 3]),
    )
    s.place("place(+inf)", pol, INF)
    pol.remove()
    r = s.place("place(1) after +inf and remove", pol, 1)
    s.expect("accepted", (r.placed, r.trace["offers"]), (1, {0: 1, 1: 0}))
    pol.remove()
    pol.remove()
    s.expect(
        "reset state", (pol.requested, pol.infinite, pol.records()), (0, False, [])
    )


@fixture(
    "collapse_projection_every_priority",
    "selection",
    ["16", "9.1", "10"],
    "Selection is a projection at creation on the full static support (review-15 "
    "C4, review-16): equal-cost parallel links of capacities [3, 2], creation "
    "residual [0, 2] (another owner holds link 0). A CostOnly EB policy with "
    "multi_edge false keeps link 1 under PreferHigherPriority with Residual or "
    "Headroom and link 0 otherwise, for both multipath settings (12 cases).",
)
def collapse_projection(s: PSnap):
    c = Ctx(s, "parallel_3_2", parallel([3, 2]))
    for mp in (False, True):
        for tie in (search.DETERMINISTIC, search.PREFER):
            for prio in search.PRIORITIES:
                fg = c.ledger()
                c.manual(fg, 0, 1, [0], 3)
                sel = search.Selection(
                    multipath=mp, multi_edge=False, tie_break=tie, priority=prio
                )
                pol = P.FlowPolicy(fg, 0, 1, co(EB, sel=sel))
                s.place(f"mp={mp} {tie} {prio}", pol, 1, fg=fg)
                want = [1] if tie == search.PREFER and prio != "Capacity" else [0]
                s.expect(f"mp={mp} {tie} {prio}", c.routes_of(pol), [want])


@fixture(
    "eb_growth_equal_route_stop",
    "growth",
    ["8b", "9.2"],
    "Revision 24: with multipath routes growth never adds a flow on a DAG "
    "structurally equal to an existing flow's route. Parallel links of capacity 1 and 5, EB "
    "CapacityAware multipath, flow_count 1, no max: the first flow (both links) is "
    "capped at 2 by round robin over the 1-unit link; growth adds a flow on the "
    "5-unit link alone (cohort [2, 2], gain 2); the next candidate is that same "
    "route and growth stops: placed 4, exhausted (before revision 24 each further "
    "equal-route flow raised the bound by a rounding remainder).",
)
def eb_growth_equal_route(s: PSnap):
    c = Ctx(s, "parallel_1_5", parallel([1, 5]))
    fg = c.ledger()
    pol = P.FlowPolicy(fg, 0, 1, ca(EB, fc=1, sel=search.Selection()))
    r = s.place("request 100", pol, 100)
    s.expect(
        "result", (r.placed, r.termination, r.routes_changed), (4, "exhausted", True)
    )
    s.expect(
        "routes and shares", (c.routes_of(pol), c.shares(pol)), ([[0, 1], [1]], [2, 2])
    )
    last = r.trace["sources"][0]["stage4"][-1]
    s.expect("growth stop reason", last["stop"], "candidate equals a cohort route")


@fixture(
    "capacityaware_blocked_inf_remove_place",
    "lifetime",
    ["8", "9.4", "13"],
    "A CapacityAware policy whose static roster is non-empty but whose residual "
    "routes are blocked by another owner: place(+inf) creates nothing (no_route, no "
    "id consumed) yet opens the infinite state, so place(1) is rejected; remove() on "
    "the empty cohort resets the epoch state and place(1) is accepted (no_route while "
    "blocked); after the owner leaves, place(1) creates the cohort with the next id.",
)
def capacityaware_blocked(s: PSnap):
    c = Ctx(s, "one_link", parallel([3]))
    fg = c.ledger()
    other, _ = c.manual(fg, 0, 1, [0], 3, "owner blocks the only link")
    pol = P.FlowPolicy(fg, 0, 1, ca(MF))
    s.expect("static roster", pol.roster, [0])
    nxt = fg.next_flow_id()
    r = s.place("place(+inf)", pol, INF)
    s.expect(
        "+inf creates nothing",
        (r.placed, r.termination, r.routes_changed, fg.next_flow_id()),
        (0, "no_route", False, nxt),
    )
    s.raises("place(1) after +inf in the same epoch", pol.place, 1)
    pol.remove()
    s.expect("epoch reset", (pol.requested, pol.infinite), (0, False))
    r = s.place("place(1) after remove", pol, 1)
    s.expect(
        "accepted, still blocked",
        (r.placed, r.termination, fg.next_flow_id()),
        (0, "no_route", nxt),
    )
    other.remove()
    r = s.place("place(1) after release", pol, 1)
    s.expect(
        "created",
        (r.placed, r.termination, r.routes_changed, pol.flow_ids()[0][1]),
        (1, "complete", True, nxt),
    )


@fixture(
    "eb_growth_single_path_repeat",
    "growth",
    ["8b", "9.2"],
    "Revision 24 (owner review): a repeated single path is a legitimate growth "
    "candidate. Parallel links of capacity 1 and 2, EB CapacityAware single-path "
    "tunnels (Deterministic, Collapse keeps the lowest id), flow_count 1, request 3: "
    "the first tunnel takes the 1-link (1), growth adds a tunnel on the 2-link "
    "(cohort 2, shares [1, 1]), then a third tunnel repeating the 2-link path (cohort "
    "3, shares [1, 1, 1]): complete with three tunnels of 1. Stopping on the repeated "
    "path would have ended at 2 (exhausted).",
)
def eb_growth_single_path(s: PSnap):
    c = Ctx(s, "parallel_1_2", parallel([1, 2]))
    fg = c.ledger()
    pol = P.FlowPolicy(
        fg, 0, 1, ca(EB, fc=1, sel=tunnel("Residual", search.DETERMINISTIC))
    )
    r = s.place("request 3", pol, 3)
    s.expect(
        "result", (r.placed, r.termination, r.routes_changed), (3, "complete", True)
    )
    s.expect(
        "three tunnels of 1",
        (c.routes_of(pol), c.shares(pol)),
        ([[0], [1], [1]], [1, 1, 1]),
    )
