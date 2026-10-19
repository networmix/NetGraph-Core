"""Prototype 2 fixtures in integer units (APPROACH revision 24, sections 8, 9 and 16).

The fixtures of the grid prototype regenerated in units (the old names mapped in
``FIXTURE_MAP``), plus the review witnesses of review 13 to review 17; the deletion of T0 removed every
``CapacityAware`` fixture and part with that routing mode (README "T1 replaces
CapacityAware"). Graphs are
stated in input coordinates; bundles and expected per-edge arrays are input
positions translated through the inverse of ``input_index()``. Expected values are
derived in the fixture text or bodies (round robin, highest averages, hand-traced
transactions), never read back from the engine. "Manual" owners are ``Global``
policies whose support is one bundle (a controller restricted to it; step 3 converted the
former single-bundle ``Pinned`` MaxFlow owners by this rule: the same minimum-cost flow on the
same DAG). CostOnly Weighted configurations with the selection priority Residual take the
residual-snapshot weight basis (step 3: weights no longer come from the priority).

Graph shapes (input positions):
* ``diamond(s0, t0, s1, t1)``: 0 -[0]-> 1 -[1]-> 3 and 0 -[2]-> 2 -[3]-> 3; route 0
  = {0, 1}, route 1 = {2, 3};
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
from units import INF

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
    "eb_1_plus_2_vs_3": "eb_1_plus_2_vs_3",
    "eb_unequal_1_10": "eb_unequal_1_10",
    "four_calls_EB": "four_calls_EB",
    "four_calls_Lossy": "four_calls_Lossy",
    "four_calls_P": "four_calls_MaxFlow",
    "four_calls_W": "four_calls_W",
    "inf_every_kernel": "inf_every_kernel",
    "lossy_diamond_10_then_1": "lossy_diamond_10_then_1",
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


def parallel(caps, costs=None):
    costs = costs or [1] * len(caps)
    return 2, [(0, 1, c, w) for c, w in zip(caps, costs, strict=True)]


def co(kernel, adm=LOSSLESS, sel=None):
    sel = sel or search.Selection()
    # By rule (step 3): Weighted with the Residual priority read its weights from the residual
    # at preparation, which is the residual-snapshot weight basis.
    weights = "ResidualSnapshot" if kernel == W and sel.priority == "Residual" else None
    return P.FlowPolicyConfig(P.CostOnly(), kernel, adm, sel, weights)


def glob():
    # Step 2: Global takes neither a placement nor a selection.
    return P.FlowPolicyConfig(P.Global(), None, LOSSLESS, None)


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

    def manual(self, fg, src, dst, positions, amount, label="manual owner"):
        """A fixed-structure owner: Global with one bundle as its support."""
        cfg = P.FlowPolicyConfig(P.Global(self.ids(positions)), None, LOSSLESS, None)
        pol = P.FlowPolicy(fg, src, dst, cfg)
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


# ====================================================================== four calls


def _four_calls(s: PSnap, kernel, adm, expected):
    for shape, caps in (("equal", (5, 5)), ("unequal", (1, 10))):
        c = Ctx(s, f"diamond_{shape}", diamond(caps[0], caps[0], caps[1], caps[1]))
        for mode in ("CostOnly",):
            fg = c.ledger()
            cfg = co(kernel, adm)
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
    "CostOnly flow on the diamond: 1, 1, 1, 1.",
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
        (a.placed, a.termination, a.remaining, a.call_drops),
        (6, "complete", 4, [(tail, 4)]),
    )
    b = s.place("request 1", pol, 1)
    s.expect(
        "second",
        (b.placed, b.termination, b.remaining, b.call_drops),
        (0, "complete", 1, [(tail, 1)]),
    )
    rec = pol.records()[0]
    s.expect("record", (rec.offered, rec.delivered, rec.drops), (11, 6, {tail: 5}))


# ====================================================================== +inf


@fixture(
    "inf_every_kernel",
    "infinite",
    ["15", "8a", "8b"],
    "+inf for every kernel: CostOnly MaxFlow, Weighted and EB on 5/5 place 10 "
    "(exhausted, remaining +inf); Lossy injects every source member to its residual: 11 "
    "on the 1/10 diamond; on 10/10 over 1/10 delivered 11 with 20 offered and 9 dropped.",
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
    "eb_1_plus_2_vs_3",
    "units",
    ["16.4", "16.5"],
    "Three ample parallel links: EB place(1); place(2) equals place(3) exactly "
    "([1, 1, 1]) for one CostOnly flow.",
)
def eb_additivity(s: PSnap):
    c = Ctx(s, "parallel_10", parallel([10, 10, 10]))
    for mode in ("CostOnly",):
        finals = []
        for calls in ([1, 2], [3]):
            fg = c.ledger()
            cfg = co(EB)
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
    "Overlapping endpoint sets are rejected at construction (rule 13). Zero capacity "
    "is not a missing route (the contract's 3.5): CostOnly over zero-capacity "
    "links has the source in its roster and creates its zero-delivery record on "
    "the static support (exhausted: MaxFlow carries 0 there).",
)
def admission(s: PSnap):
    c = Ctx(s, "diamond_5", diamond(5, 5, 5, 5))
    fg = c.ledger()
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


# ====================================================================== flags


@fixture(
    "flag_zero_delivery_creation",
    "flags",
    ["9.4", "8"],
    "Revision-20 witness: a first place(1) on a static route whose unit another "
    "owner holds creates a zero-delivery record: routes_changed true, termination "
    "exhausted, the record stays (CostOnly); a second zero-delivery call keeps the same "
    "ids and sets no flag (revision 21).",
)
def flag_zero_delivery(s: PSnap):
    c = Ctx(s, "one_unit", parallel([1]))
    for mode in ("CostOnly",):
        fg = c.ledger()
        c.manual(fg, 0, 1, [0], 1)
        cfg = co(EB)
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
    "Lossy with pooled sources (replaces the zero-cost forward-versus-reverse witness; "
    "step 2, contract r8 5.1 and 9.1): sources 0 and 1 reach 3 through X (0->X 4, 1->X "
    "4, X->3 3) and source 1 also through Y (1->Y 4, Y->3 4, cost 2 each). Request 8 "
    "offers [4, 4]; source 1's routes are its own shortest walks (1->X->3 only, cost 2 < "
    "4); router X holds 8 and sends it all on X->3, attributed [4, 4]; the link carries "
    "3 = H(3, [4, 4]) = [2, 1]: source 0 delivers 2 and drops 2, source 1 delivers 1 "
    "and drops 3. Delivered 3, complete.",
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
        (2, [0, 0, 2, 0, 0]),
    )
    s.expect(
        "source 1",
        (recs[1].delivered, dict_to_input(c.g, recs[1].drops)),
        (1, [0, 0, 3, 0, 0]),
    )


@fixture(
    "sources_attribution",
    "sources",
    ["24", "26", "12", "9.3"],
    "Sources A (A->T 1, cost 1) and B (B->M 10, M->T 2, cost 1 each). Lossy "
    "request 20 (offers [10, 10]): A delivers 1 and drops 9 on A->T; B delivers 2 "
    "and drops 8 on M->T; placed 3, complete. Lossless (step 2): MaxFlow is a controller "
    "with the optimized origin, 1 + 2 = 3; Weighted and EB have the equal origin, the "
    "largest total whose round-robin shares fit: 2 ([1, 1]; 3 would give A 2). All "
    "exhausted. Route costs per record: 1 and 2.",
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
    for kernel, want in ((MF, (3, [1, 2])), (W, (2, [1, 1])), (EB, (2, [1, 1]))):
        fg = c.ledger()
        pol = P.FlowPolicy(fg, [0, 1], 3, co(kernel))
        r = s.place(f"{kernel} 20", pol, 20, fg=fg)
        s.expect(
            f"{kernel}",
            (r.placed, r.termination, [x.delivered for x in pol.records()]),
            (want[0], "exhausted", want[1]),
        )


@fixture(
    "sources_eb_per_source",
    "sources",
    ["26", "9.3"],
    "EqualBalanced with the equal origin (step 2, contract r8 3.4 and 5.1): A (A->T 4, "
    "cost 1) and B (B->T 10, cost 2). Lossless admits the largest total whose round-robin "
    "shares fit: 8 ([4, 4]; 9 would give A 5). +inf also 8. (Independent per-source "
    "placement, the rule before step 2, gave 4 + 10 = 14.)",
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
            (8, term, [4, 4]),
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
    "+inf on the attribution graph: MaxFlow 3 (optimized origin), Weighted 2 and EB 2 "
    "(equal origin: [1, 1]) (exhausted); Lossy injects every source member to its "
    "residual (A 1, B 10): delivered 3, offered 11, drops 8, exhausted.",
)
def sources_inf(s: PSnap):
    c = Ctx(s, "attribution", (4, [(0, 3, 1, 1), (1, 2, 10, 1), (2, 3, 2, 1)]))
    for kernel, want in ((MF, 3), (W, 2), (EB, 2)):
        fg = c.ledger()
        r = s.place(
            f"{kernel} +inf", P.FlowPolicy(fg, [0, 1], 3, co(kernel)), INF, fg=fg
        )
        s.expect(f"{kernel}", (r.placed, r.termination), (want, "exhausted"))
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
    "Endpoint-set validation: an empty source set and overlapping sets are "
    "rejected; an excluded source is skipped from the roster "
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
    "Shortest routes versus the whole graph: 0->2 (cost 2), 1->2 (cost 1), 2->3 "
    "(cost 1) and 0->3 (cost 5), unit capacities. CostOnly MaxFlow (step 2) is one "
    "controller on the union of the sources' shortest walks, 0->2->3 (cost 3) and "
    "1->2->3: the shared edge 2->3 limits +inf to 1 (0->3 is on no source's shortest "
    "walk). Global's unique optimum (cost 7 against 8) routes source 0 over 0->3 and "
    "source 1 over 1->2->3: 2.",
)
def sources_vs_global(s: PSnap):
    c = Ctx(
        s, "shared_edge", (4, [(0, 2, 1, 2), (1, 2, 1, 1), (2, 3, 1, 1), (0, 3, 1, 5)])
    )
    fg = c.ledger()
    r = s.place(
        "CostOnly MaxFlow +inf", P.FlowPolicy(fg, [0, 1], 3, co(MF)), INF, fg=fg
    )
    s.expect("union of the shortest walks", r.placed, 1)
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
    "Cumulative requests (16.7; step 2): sources A, B, C with direct links of 10. "
    "CostOnly MaxFlow is a controller with the optimized origin that re-solves the "
    "cumulative value from zero: place(1); place(1) supplies [1, 0, 0] then [2, 0, 0] "
    "(the root's arc to the lowest NodeId first, its path takes both units), identical "
    "to place(2). Lossy follows the cumulative round-robin split: two calls of 1 over "
    "two sources offer [1, 0] then [0, 1].",
)
def sources_split(s: PSnap):
    c = Ctx(s, "three_direct", (4, [(0, 3, 10, 1), (1, 3, 10, 1), (2, 3, 10, 1)]))
    fg = c.ledger()
    pol = P.FlowPolicy(fg, [0, 1, 2], 3, co(MF))
    a = s.place("call 1", pol, 1)
    b = s.place("call 2", pol, 1)
    s.expect(
        "supplies",
        (a.trace["flow_core"]["supply"], b.trace["flow_core"]["supply"]),
        ({0: 1, 1: 0, 2: 0}, {0: 2, 1: 0, 2: 0}),
    )
    fg2 = c.ledger()
    s.place("one call of 2", P.FlowPolicy(fg2, [0, 1, 2], 3, co(MF)), 2, fg=fg2)
    s.expect("1 + 1 equals 2", c.alloc(fg), c.alloc(fg2))
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
    "Review-13 source budget witness under step 2's origins: sources with capacities "
    "[1, 9], request 4. MaxFlow (the optimized origin) supplies [1, 3]: placed 4, "
    "complete. Weighted and EB (the equal origin) admit the largest total whose "
    "round-robin shares fit: 2 ([1, 1]; 3 would give the 1-link source 2), remaining 2, "
    "exhausted.",
)
def sources_1_9(s: PSnap):
    c = Ctx(s, "two_direct_1_9", two_direct((1, 9)))
    for kernel, want in (
        (MF, ([1, 3], 4, 0, "complete")),
        (W, ([1, 1], 2, 2, "exhausted")),
        (EB, ([1, 1], 2, 2, "exhausted")),
    ):
        fg = c.ledger()
        pol = P.FlowPolicy(fg, [0, 1], 2, co(kernel))
        r = s.place(f"{kernel} request 4", pol, 4, fg=fg)
        s.expect(
            f"{kernel}",
            (
                [x.delivered for x in pol.records()],
                r.placed,
                r.remaining,
                r.termination,
            ),
            want,
        )


@fixture(
    "sources_partition_1_1_vs_2_0",
    "sources",
    ["24", "16.4", "9.2"],
    "Source-set allocation of CostOnly MaxFlow (review-14 R3; step 2): 0->2 and 1->2 "
    "(2 each) share 2->3 (2). The controller re-solves the cumulative value from zero "
    "with one virtual source, so place(2); place(2) and place(4) both give deliveries "
    "[2, 0] (the root's arc to source 0 comes first): partition-independent (before step "
    "2 the independent per-source placement gave [1, 1] and [2, 0]).",
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
    s.expect("partition independence", out, [[2, 0], [2, 0]])


# ====================================================================== termination


@fixture(
    "termination_complete_and_no_route",
    "termination",
    ["9.4", "14"],
    "The terminations of the 14.2 table without a reduction: a Lossy finite offer "
    "consumed with drops is complete; v = 0 is complete; an empty static roster (no "
    "walk from the source; a zero-capacity walk would count) is no_route.",
)
def termination(s: PSnap):
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
    s.expect(
        "empty roster",
        (r.termination, r.by_source),
        ("no_route", [(0, None, 0, 0, 0, 0)]),
    )
    s.place("empty roster +inf", p0, INF)
    s.raises("finite after +inf with an empty roster", p0.place, 1)


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
    "2 (min_cut None: a residual path is left; requested with the min_cut output), "
    "the cumulative solve passes the certificate; a mutated producer returning 1 "
    "under cap 2 is rejected; place(+inf) then places the last unit and min_cut is "
    "[0].",
)
def global_underfill(s: PSnap):
    c = Ctx(s, "one_edge", (2, [(0, 1, 3, 7)]))
    fg = c.ledger()
    pol = P.FlowPolicy(fg, 0, 1, glob())
    r = s.place("place(2)", pol, 2, outputs=P.Outputs(bottleneck=True, min_cut=True))
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
    r = s.place(
        "place(+inf)", pol, INF, outputs=P.Outputs(bottleneck=True, min_cut=True)
    )
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
    "Section 8's epoch law for CostOnly MaxFlow over two sources: place(1) supplies "
    "[1, 0] and creates both records (ids [0, 1], source 1 at zero delivery); remove(); "
    "place(1) supplies [1, 0] again (requested reset) with ids [2, 3]; "
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
            r.trace["flow_core"]["supply"],
            [x.flow_id[1] for x in pol.records()],
            [x.delivered for x in pol.records()],
        ),
        ({0: 1, 1: 0}, [0, 1], [1, 0]),
    )
    pol.remove()
    r = s.place("place(1) after remove", pol, 1)
    s.expect(
        "epoch 2",
        (r.trace["flow_core"]["supply"], [x.flow_id[1] for x in pol.records()]),
        ({0: 1, 1: 0}, [2, 3]),
    )
    s.place("place(+inf)", pol, INF)
    pol.remove()
    r = s.place("place(1) after +inf and remove", pol, 1)
    s.expect("accepted", (r.placed, r.trace["flow_core"]["supply"]), (1, {0: 1, 1: 0}))
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
