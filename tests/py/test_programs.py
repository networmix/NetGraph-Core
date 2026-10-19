"""Forwarding programs: what every model exports, Global with a support, the weight
bases, ForwardingProgram.from_routes, Pinned over programs, the lossless bound, the loss rule
inside one program (cycles included), replay under a failure mask, stale forwarding with black
holes, the operations evaluate, replace and forward with their epochs and flags, the per-entry
records and the by_tunnel rows. Every expected value is derived by hand from the contract's rules;
none is read from a build. The independent
reference of the property tests is program_oracle.py."""

from __future__ import annotations

import pytest
from helpers import INF, exact, full_ledger_state, sweep_python_allocation
from private_hooks import needs_hooks

import netgraph_core as ngc

LOSSY = ngc.Admission.LOSSY
LOSSLESS = ngc.Admission.LOSSLESS
EB = ngc.FlowPlacement.EQUAL_BALANCED
W = ngc.FlowPlacement.WEIGHTED
MF = ngc.FlowPlacement.MAX_FLOW
UP = ngc.TunnelStatus.UP
FAILED = ngc.TunnelStatus.FAILED
BLACKHOLED = ngc.TunnelStatus.BLACKHOLED
T = ngc.Termination


# ---------------------------------------------------------------- builders


def graph(n, rows):
    """A Graph from (u, v, capacity[, cost]) rows written in EdgeId order (checked)."""
    g = ngc.Graph(
        n,
        [r[0] for r in rows],
        [r[1] for r in rows],
        [r[2] for r in rows],
        [r[3] if len(r) > 3 else 1 for r in rows],
        1.0,
    )
    assert [int(i) for i in g.input_index()] == list(range(len(rows)))
    return g


def path(g, edges):
    return ngc.Path(g, [[e] for e in edges])


def dag(g, edges):
    """The PredDAG of a set of EdgeIds."""
    per = [[] for _ in range(g.num_nodes())]
    for e in sorted(edges):
        per[int(g.edge_dst()[e])].append(e)
    offsets, parents, via = [0], [], []
    for v in range(g.num_nodes()):
        for e in per[v]:
            parents.append(int(g.edge_src()[e]))
            via.append(e)
        offsets.append(len(via))
    return ngc.PredDAG.from_csr(g, offsets, parents, via)


def tunnels(g, sources, destinations, specs, **options):
    """A supplied program of tunnels: specs are (head, [edges], weight)."""
    routes = [ngc.Route(h, path(g, es), weight=w) for h, es, w in specs]
    return ngc.ForwardingProgram.from_routes(
        g, sources, destinations, routes, **options
    )


def pinned(g, program, admission=LOSSLESS, elig=None, black_holes=None, ledger=None):
    cfg = ngc.FlowPolicyConfig(
        routing=ngc.Pinned(program, black_holes=black_holes), admission=admission
    )
    return ngc.FlowPolicy(
        ledger or ngc.FlowLedger(g),
        list(program.sources),
        list(program.destinations),
        cfg,
        elig,
    )


def cost_only(placement=EB, admission=LOSSLESS, **fields):
    return ngc.FlowPolicyConfig(
        routing=ngc.CostOnly(), placement=placement, admission=admission, **fields
    )


def global_(origin=None, support=None):
    return ngc.FlowPolicyConfig(routing=ngc.Global(support=support), origin=origin)


def rows(result):
    """The by_tunnel rows as (index, head, status, assigned, delivered, dropped)."""
    return [
        (t.index, t.head, t.status, t.assigned, t.delivered, t.dropped)
        for t in result.by_tunnel
    ]


def listing(program):
    """A program's tunnels as (head, edges, weight)."""
    return [(t.head, list(t.edges), t.weight) for t in program.tunnels]


def par2(a, b):
    return graph(2, [(0, 1, a), (0, 1, b)])


ONE = [(0, 1, 10)]
B1 = [(0, 3, 2), (1, 2, 10), (2, 3, 10, 2)]
B2 = [(0, 1, 4), (0, 2, 2), (1, 2, 10)]
ES = [(0, 1, 10), (0, 2, 10), (1, 3, 2), (2, 3, 10, 2)]
G1U = [(0, 1, 2), (0, 2, 10), (1, 3, 2), (2, 3, 10)]
PINLOSS = [(0, 1, 10), (0, 2, 10), (1, 3, 2), (2, 3, 10)]
DIAM = [(0, 1, 100), (0, 2, 100), (1, 3, 100), (2, 3, 100)]
POOL2 = [(0, 2, 10), (1, 2, 10), (2, 3, 10), (2, 3, 10)]
REL = [(0, 2, 10), (1, 2, 10), (2, 3, 10)]
TUN = [(0, 1, 10), (0, 2, 10), (1, 3, 10), (2, 3, 10)]
MB = [(0, 1, 10), (0, 2, 10), (0, 3, 10), (1, 4, 10), (2, 4, 10), (3, 4, 10)]
SURV = [(0, 3, 3), (1, 3, 1), (2, 3, 5)]
RB = [(0, 2, 5), (1, 2, 5), (1, 2, 5, 2)]
INC = [(0, 2, 1), (1, 2, 3)]
SHARE = [(0, 1, 4), (1, 2, 10), (1, 2, 10)]
HIER = [(0, 2, 1000), (1, 2, 1000), (2, 3, 1)]


def cyc(x, y):
    """CYC(x, y): tunnel 0 = [e6, e0, e1, e3, e5] meets X (e0) then Y (e3); tunnel 1 = [e7, e3,
    e4, e0, e2] meets Y then X."""
    return graph(
        6,
        [
            (0, 1, x),
            (1, 2, 1000),
            (1, 5, 1000),
            (2, 3, y),
            (3, 0, 1000),
            (3, 5, 1000),
            (4, 0, 1000),
            (4, 2, 1000),
        ],
    )


CYC_T0 = [6, 0, 1, 3, 5]
CYC_T1 = [7, 3, 4, 0, 2]


# ---------------------------------------------------------------- P: exported programs


def test_global_export():
    # B1, Global optimized and equal, place(12).
    g = graph(4, B1)
    p = ngc.FlowPolicy(ngc.FlowLedger(g), [0, 1], 3, global_())
    r = p.place(12)
    prog = r.program
    assert listing(prog) == [(0, [0], 2), (1, [1, 2], 10)]
    assert [t.index for t in prog.tunnels] == [0, 1]
    assert [t.local_index for t in prog.tunnels] == [0, 0]
    assert [t.destination for t in prog.tunnels] == [3, 3]
    assert [t.route_cost for t in prog.tunnels] == [1, 3]
    assert prog.source_weights == ((0, 2), (1, 10))
    assert prog.sources == (0, 1) and prog.destinations == (3,)
    assert [rec.entry for rec in p.records()] == [None, None]
    assert rows(r) == [(0, 0, UP, 2, 2, 0), (1, 1, UP, 10, 10, 0)]
    assert p.program.tunnels == prog.tunnels or listing(p.program) == listing(prog)
    eq = ngc.FlowPolicy(ngc.FlowLedger(g), [0, 1], 3, global_(ngc.Origin.EQUAL)).place(
        12
    )
    assert eq.placed == 4
    assert listing(eq.program) == [(0, [0], 2), (1, [1, 2], 2)]
    assert eq.program.source_weights is None


def test_single_source_decomposed():
    # One source is decomposed too; CostOnly MAX_FLOW.
    g = graph(4, ES)
    r = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 3, global_()).place(12)
    assert listing(r.program) == [(0, [0, 2], 2), (0, [1, 3], 10)]
    assert [t.local_index for t in r.program.tunnels] == [0, 1]
    assert [t.route_cost for t in r.program.tunnels] == [2, 3]
    assert r.total_cost == 34
    assert r.program.source_weights == ((0, 12),)
    mf = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 3, cost_only(MF)).place(12)
    assert mf.placed == 2
    assert listing(mf.program) == [(0, [0, 2], 2)]
    assert mf.program.source_weights == ((0, 2),)


def test_cost_only_split_program():
    # A CostOnly split program lists no tunnel.
    g = graph(4, POOL2)
    r = ngc.FlowPolicy(ngc.FlowLedger(g), [0, 1], 3, cost_only()).place(2)
    assert r.program.tunnels == ()
    assert r.program.source_weights is None
    assert r.program.sources == (0, 1) and r.program.destinations == (3,)
    assert r.by_tunnel == ()


def test_zero_supply_kept_and_lifecycle():
    # Every roster source has a source weight, zero included.
    g = graph(3, [(0, 2, 5), (1, 2, 5)])
    r = ngc.FlowPolicy(ngc.FlowLedger(g), [0, 1], 2, global_()).place(1)
    assert listing(r.program) == [(0, [0], 1)]
    assert r.program.source_weights == ((0, 1), (1, 0))
    # None before the first positive call and after remove().
    p = ngc.FlowPolicy(ngc.FlowLedger(g), [0, 1], 2, global_())
    assert p.program is None
    z = p.place(0)
    assert z.program is None and z.by_tunnel == ()
    p.place(3)
    assert p.program is not None
    p.remove()
    assert p.program is None and p.place(0).program is None


def test_replay_reproduces_the_solve_and_keeps_the_mask():
    # B1's optimized program forwarded by Pinned reproduces the solve.
    g = graph(4, B1)
    prog = ngc.FlowPolicy(ngc.FlowLedger(g), [0, 1], 3, global_()).place(12).program
    p = pinned(g, prog)
    r = p.place(12, ngc.OutputOptions(edge_flows=True))
    assert (r.placed, r.termination) == (12, T.COMPLETE)
    assert r.edge_flows.tolist() == [2, 10, 10]
    assert pinned(g, prog).place(INF).placed == 12
    # The policy's program is the configured one.
    assert listing(p.program) == listing(prog) and listing(r.program) == listing(prog)
    # The producer's mask is saved: no replay brings source 1 back.
    masked = ngc.FlowPolicy(
        ngc.FlowLedger(g), [0, 1], 3, global_(), g.eligibility(excluded_nodes=[1])
    ).place(12)
    assert masked.placed == 2 and listing(masked.program) == [(0, [0], 2)]
    lossy = pinned(g, masked.program, LOSSY).replace(12)
    assert lossy.placed == 2 and [row.offered for row in lossy.by_source] == [12, 0]


def test_empty_roster_program():
    # An empty roster installs an empty program.
    g = graph(2, ONE)
    p = ngc.FlowPolicy(
        ngc.FlowLedger(g), 0, 1, cost_only(), g.eligibility(excluded_nodes=[0])
    )
    r = p.place(1)
    assert r.termination == T.NO_ROUTE
    assert r.program is not None and r.program.tunnels == ()
    assert r.program.sources == (0,) and not r.routes_changed and not r.weights_changed


# ---------------------------------------------------------------- G: Global with a support


def test_global_support():
    # Only the trimmed support carries traffic.
    g = graph(4, ES)
    r = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 3, global_(support=dag(g, [1, 3]))).place(
        12, ngc.OutputOptions(min_cut=True)
    )
    assert (r.placed, r.total_cost) == (10, 30)
    assert r.min_cut is None and r.cost_distribution is None
    assert listing(r.program) == [(0, [1, 3], 10)]
    assert ngc.FlowPolicy(ngc.FlowLedger(g), 0, 3, global_()).place(12).placed == 12
    r2 = ngc.FlowPolicy(
        ngc.FlowLedger(g), 0, 3, global_(support=dag(g, [0, 1, 3]))
    ).place(12, ngc.OutputOptions(edge_flows=True))
    assert r2.placed == 10 and r2.edge_flows.tolist() == [0, 10, 0, 10]
    # The mask removes the support: an empty roster.
    r3 = ngc.FlowPolicy(
        ngc.FlowLedger(g),
        0,
        3,
        global_(support=dag(g, [1, 3])),
        g.eligibility(excluded_edges=[3]),
    ).place(12)
    assert (r3.termination, r3.placed, r3.program.tunnels) == (T.NO_ROUTE, 0, ())
    # The equal origin on a support.
    b2 = graph(3, B2)
    r4 = ngc.FlowPolicy(
        ngc.FlowLedger(b2), [0, 1], 2, global_(ngc.Origin.EQUAL, dag(b2, [1, 2]))
    ).place(12)
    assert (r4.placed, r4.total_cost) == (4, 4)


def test_global_support_wide_cost_and_ownership():
    # total_cost above 2^64.
    g = graph(2, [(0, 1, 1 << 40, 1 << 30)])
    r = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 1, global_(support=dag(g, [0]))).place(
        1 << 40
    )
    assert r.total_cost == 1 << 70
    # A support of another graph.
    other = graph(2, [(0, 1, 1)])
    with pytest.raises(ValueError, match=exact("support belongs to another graph")):
        ngc.FlowPolicy(ngc.FlowLedger(g), 0, 1, global_(support=dag(other, [0])))


# ---------------------------------------------------------------- W: weight bases


def test_explicit_weights():
    # Configured (explicit) weights.
    g = par2(10, 10)
    cfg = cost_only(W, LOSSY, weights=ngc.WeightBasis.EXPLICIT, edge_weights=[1, 2])
    r = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 1, cfg).place(
        12, ngc.OutputOptions(edge_flows=True)
    )
    assert r.edge_flows.tolist() == [4, 8]
    cap = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 1, cost_only(W, LOSSY)).place(
        12, ngc.OutputOptions(edge_flows=True)
    )
    assert cap.edge_flows.tolist() == [6, 6]
    ll = cost_only(W, weights=ngc.WeightBasis.EXPLICIT, edge_weights=[1, 2])
    assert ngc.FlowPolicy(ngc.FlowLedger(g), 0, 1, ll).place(INF).placed == 15
    assert ngc.capacity(g, 0, 1, ngc.QueryOptions(ll)).placed == 15
    # A positive weight on an edge of capacity 0 becomes 0.
    z = par2(0, 10)
    cfg = cost_only(W, LOSSY, weights=ngc.WeightBasis.EXPLICIT, edge_weights=[5, 1])
    r = ngc.FlowPolicy(ngc.FlowLedger(z), 0, 1, cfg).place(7)
    assert (r.placed, r.call_drops) == (7, [])
    # A zero weight prunes.
    cfg = cost_only(W, weights=ngc.WeightBasis.EXPLICIT, edge_weights=[0, 3])
    assert ngc.FlowPolicy(ngc.FlowLedger(g), 0, 1, cfg).place(12).placed == 10
    cfg = cost_only(W, LOSSY, weights=ngc.WeightBasis.EXPLICIT, edge_weights=[0, 3])
    r = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 1, cfg).place(12)
    assert (r.placed, r.call_drops) == (10, [(1, 2)])
    # Every weight 0.
    for adm, loss in ((LOSSY, 5), (LOSSLESS, 0)):
        cfg = cost_only(W, adm, weights=ngc.WeightBasis.EXPLICIT, edge_weights=[0, 0])
        r = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 1, cfg).place(5)
        assert (r.termination, r.placed, r.source_loss) == (T.NO_ACTIVE_WEIGHT, 0, loss)


def test_snapshot_weights():
    # Snapshot weights: the priority is no weight, the snapshot is frozen.
    g = par2(10, 10)

    def policy(admission, **fields):
        ledger = ngc.FlowLedger(g, [2, 10])
        return ngc.FlowPolicy(ledger, 0, 1, cost_only(W, admission, **fields))

    assert policy(LOSSLESS).place(12).placed == 4
    snap = ngc.WeightBasis.RESIDUAL_SNAPSHOT
    assert policy(LOSSLESS, weights=snap).place(12).placed == 12
    prio = ngc.PathSelection(priority=ngc.EdgePriority.RESIDUAL)
    assert policy(LOSSLESS, selection=prio).place(12).placed == 4
    p = policy(LOSSY, weights=snap)
    r1 = p.place(12, ngc.OutputOptions(edge_flows=True))
    assert r1.edge_flows.tolist() == [2, 10]
    r2 = p.place(6)
    assert (r2.placed, r2.call_drops, r2.call_source_loss) == (0, [(0, 1), (1, 5)], 0)
    # Weighted with the Headroom priority is accepted (CostOnly).
    sel = ngc.PathSelection(
        ngc.EdgeTieBreak.PREFER_HIGHER_PRIORITY, ngc.EdgePriority.HEADROOM, False
    )
    assert policy(LOSSLESS, selection=sel).place(12).placed == 10


def test_weight_texts():
    g = par2(10, 10)
    ledger = ngc.FlowLedger(g)
    cases = [
        (
            cost_only(weights=ngc.WeightBasis.CAPACITY),
            "weights require CostOnly placement=WEIGHTED",
        ),
        (cost_only(edge_weights=[1, 1]), "edge_weights require weights=EXPLICIT"),
        (
            cost_only(W, weights=ngc.WeightBasis.EXPLICIT),
            "weights=EXPLICIT requires edge_weights",
        ),
        (
            cost_only(W, weights=ngc.WeightBasis.CAPACITY, edge_weights=[1, 1]),
            "edge_weights require weights=EXPLICIT",
        ),
        (
            ngc.FlowPolicyConfig(
                routing=ngc.Global(), weights=ngc.WeightBasis.CAPACITY
            ),
            "Global does not accept weights",
        ),
        (
            ngc.FlowPolicyConfig(
                routing=ngc.CapacityAware(ngc.Tunnels(1)),
                weights=ngc.WeightBasis.CAPACITY,
            ),
            "T1 does not accept weights",
        ),
        (
            cost_only(W, weights=ngc.WeightBasis.EXPLICIT, edge_weights=[1, 2, 3]),
            "edge_weights: expected 2 entries",
        ),
        (
            cost_only(W, weights=ngc.WeightBasis.EXPLICIT, edge_weights=[1, 1 << 62]),
            "edge_weights[1] must be an integer in [0, 2^62)",
        ),
        (
            cost_only(W, weights=ngc.WeightBasis.EXPLICIT, edge_weights=[-1, 1]),
            "edge_weights[0] must be an integer in [0, 2^62)",
        ),
    ]
    for cfg, text in cases:
        with pytest.raises(ValueError, match=exact(text)):
            ngc.FlowPolicy(ledger, 0, 1, cfg)
    prog = tunnels(g, [0], [1], [(0, [0], 1)])
    with pytest.raises(ValueError, match=exact("Pinned does not accept weights")):
        ngc.FlowPolicy(
            ledger,
            0,
            1,
            ngc.FlowPolicyConfig(
                routing=ngc.Pinned(prog), weights=ngc.WeightBasis.CAPACITY
            ),
        )
    with pytest.raises(ValueError):
        cost_only(W, weights=ngc.WeightBasis.EXPLICIT, edge_weights=[True, 1])
    cfg = cost_only(W, weights=ngc.WeightBasis.EXPLICIT, edge_weights=[1, 2])
    assert cfg.edge_weights == (1, 2)


# ---------------------------------------------------------------- F: from_routes


def test_from_routes_texts():
    g = par2(4, 8)
    p0 = path(g, [0])
    d01 = dag(g, [0, 1])
    fr = ngc.ForwardingProgram.from_routes

    def raises(text, *args, **kwargs):
        with pytest.raises(ValueError, match=exact(text)):
            fr(g, *args, **kwargs)

    raises("sources: the node set is empty", [], [1], [ngc.Route(0, p0)])
    raises(
        "node 0 is in both endpoint sets: sources and destinations must not overlap",
        [0],
        [0],
        [ngc.Route(0, p0)],
    )
    raises("routes must not be empty", [0], [1], [])
    raises(
        "placement applies only to PredDAG routes",
        [0],
        [1],
        [ngc.Route(0, p0)],
        placement=EB,
    )
    raises(
        "weights applies only to PredDAG routes",
        [0],
        [1],
        [ngc.Route(0, p0)],
        weights=ngc.WeightBasis.CAPACITY,
    )
    raises(
        "edge_weights applies only to PredDAG routes",
        [0],
        [1],
        [ngc.Route(0, p0)],
        edge_weights=[1, 1],
    )
    raises(
        "residual applies only to PredDAG routes",
        [0],
        [1],
        [ngc.Route(0, p0)],
        residual=[1, 1],
    )
    dagr = [ngc.Route(0, d01)]
    raises(
        "from_routes placement must be EQUAL_BALANCED or WEIGHTED",
        [0],
        [1],
        dagr,
        placement=MF,
    )
    raises(
        "weights require placement=WEIGHTED",
        [0],
        [1],
        dagr,
        weights=ngc.WeightBasis.CAPACITY,
    )
    raises(
        "weights=EXPLICIT requires edge_weights",
        [0],
        [1],
        dagr,
        placement=W,
        weights=ngc.WeightBasis.EXPLICIT,
    )
    raises(
        "edge_weights require weights=EXPLICIT",
        [0],
        [1],
        dagr,
        placement=W,
        edge_weights=[1, 1],
    )
    raises(
        "weights=RESIDUAL_SNAPSHOT requires residual",
        [0],
        [1],
        dagr,
        placement=W,
        weights=ngc.WeightBasis.RESIDUAL_SNAPSHOT,
    )
    raises(
        "residual requires weights=RESIDUAL_SNAPSHOT",
        [0],
        [1],
        dagr,
        placement=W,
        residual=[4, 8],
    )
    raises(
        "edge_weights: expected 2 entries",
        [0],
        [1],
        dagr,
        placement=W,
        weights=ngc.WeightBasis.EXPLICIT,
        edge_weights=[1],
    )
    raises(
        "edge_weights[1] must be an integer in [0, 2^62)",
        [0],
        [1],
        dagr,
        placement=W,
        weights=ngc.WeightBasis.EXPLICIT,
        edge_weights=[1, 1 << 62],
    )
    raises(
        "residual[0] = 5 is outside [0, capacity_units[0] = 4]",
        [0],
        [1],
        dagr,
        placement=W,
        weights=ngc.WeightBasis.RESIDUAL_SNAPSHOT,
        residual=[5, 8],
    )
    other = par2(4, 8)
    raises(
        "routes[0]: bundle belongs to another graph",
        [0],
        [1],
        [ngc.Route(0, path(other, [0]))],
    )
    raises("routes[0]: head 1 is not a source", [0], [1], [ngc.Route(1, p0)])
    hop2 = ngc.Path(g, [[0, 1]])
    raises("routes[0]: a tunnel has one edge per hop", [0], [1], [ngc.Route(0, hop2)])
    d = graph(4, DIAM)
    fd = ngc.ForwardingProgram.from_routes
    with pytest.raises(
        ValueError, match=exact("routes[0]: the path starts at 0, not at head 1")
    ):
        fd(d, [0, 1], [3], [ngc.Route(1, path(d, [0, 2]))])
    with pytest.raises(
        ValueError, match=exact("routes[0]: the path ends at 1, not at a destination")
    ):
        fd(d, [0], [3], [ngc.Route(0, path(d, [0]))])
    with pytest.raises(
        ValueError,
        match=exact("routes[0]: the path meets destination 1 before its end"),
    ):
        fd(d, [0], [1, 3], [ngc.Route(0, path(d, [0, 2]))])
    with pytest.raises(
        ValueError, match=exact("routes[0]: no walk from head 0 to a destination")
    ):
        fd(d, [0], [3], [ngc.Route(0, dag(d, [2]))])
    for weight in (1 << 62, -1):
        raises(
            "routes[0]: weight must be an integer in [0, 2^62)",
            [0],
            [1],
            [ngc.Route(0, p0, weight=weight)],
        )
    with pytest.raises(ValueError):
        ngc.Route(0, p0, weight=True)
    # The routes in list order, each route's checks in order.
    raises(
        "routes[0]: weight must be an integer in [0, 2^62)",
        [0],
        [1],
        [ngc.Route(0, p0, weight=1 << 62), ngc.Route(5, p0)],
    )


def test_from_routes_entries():
    g = par2(4, 8)
    prog = ngc.ForwardingProgram.from_routes(
        g, [0], [1], [ngc.Route(0, path(g, [0]), weight=0), ngc.Route(0, path(g, [1]))]
    )
    assert [(t.index, t.weight) for t in prog.tunnels] == [(0, 0), (1, 1)]
    mixed = ngc.ForwardingProgram.from_routes(
        g,
        [0],
        [1],
        [
            ngc.Route(0, path(g, [1])),
            ngc.Route(0, dag(g, [0, 1])),
            ngc.Route(0, path(g, [0])),
        ],
    )
    assert [(t.index, t.local_index) for t in mixed.tunnels] == [(0, 0), (2, 2)]
    d = graph(4, DIAM)
    norm = ngc.ForwardingProgram.from_routes(
        d, [1, 0, 1], [3], [ngc.Route(0, path(d, [0, 2]))]
    )
    assert norm.sources == (0, 1)
    with pytest.raises(TypeError):
        ngc.ForwardingProgram.from_routes(g, [0], [1], [path(g, [0])])


# ---------------------------------------------------------------- C: Pinned's configuration


def test_pinned_texts():
    g = par2(4, 8)
    prog = tunnels(g, [0], [1], [(0, [0], 1)])
    ledger = ngc.FlowLedger(g)
    fields = [
        ({"placement": EB}, "Pinned does not accept placement"),
        ({"origin": ngc.Origin.EQUAL}, "Pinned does not accept origin"),
        ({"selection": ngc.PathSelection()}, "Pinned does not accept selection"),
        ({"weights": ngc.WeightBasis.CAPACITY}, "Pinned does not accept weights"),
        ({"edge_weights": [1, 1]}, "Pinned does not accept edge_weights"),
    ]
    for kw, text in fields:
        cfg = ngc.FlowPolicyConfig(routing=ngc.Pinned(prog), **kw)
        with pytest.raises(ValueError, match=exact(text)):
            ngc.FlowPolicy(ledger, 0, 1, cfg)
    three = graph(3, [(0, 1, 1), (0, 2, 1)])
    p3 = tunnels(three, [0], [1], [(0, [0], 1)])
    with pytest.raises(
        ValueError, match=exact("program endpoints do not match policy endpoints")
    ):
        ngc.FlowPolicy(
            ngc.FlowLedger(three),
            0,
            [1, 2],
            ngc.FlowPolicyConfig(routing=ngc.Pinned(p3)),
        )
    foreign = tunnels(par2(4, 8), [0], [1], [(0, [0], 1)])
    for kw in ({}, {"placement": EB}):
        cfg = ngc.FlowPolicyConfig(routing=ngc.Pinned(foreign), **kw)
        with pytest.raises(ValueError, match=exact("program belongs to another graph")):
            ngc.FlowPolicy(ledger, 0, 1, cfg)
    other = par2(4, 8)
    cfg = ngc.FlowPolicyConfig(
        routing=ngc.Pinned(prog, black_holes=other.eligibility(excluded_edges=[0]))
    )
    with pytest.raises(ValueError, match=exact("black_holes belongs to another graph")):
        ngc.FlowPolicy(ledger, 0, 1, cfg)
    with pytest.raises(
        ValueError, match=exact("foreign Eligibility: built for another graph")
    ):
        ngc.FlowPolicy(
            ledger,
            0,
            1,
            ngc.FlowPolicyConfig(routing=ngc.Pinned(prog)),
            other.eligibility(excluded_edges=[0]),
        )
    with pytest.raises(ValueError, match=exact("sources: the node set is empty")):
        ngc.FlowPolicy(ledger, [], 1, ngc.FlowPolicyConfig(routing=ngc.Pinned(foreign)))
    with pytest.raises(TypeError):
        ngc.Pinned(bundles=[dag(g, [0])])
    with pytest.raises(TypeError):
        ngc.Pinned([dag(g, [0])])
    with pytest.raises(TypeError):
        ngc.Pinned(None)
    # Several sources and destinations.
    d = graph(4, REL)
    many = tunnels(d, [0, 1], [3], [(0, [0, 2], 1), (1, [1, 2], 1)])
    assert pinned(d, many).place(INF).placed == 10


# ---------------------------------------------------------------- L: lossless bounds


def bound(g, specs, sources=(0,), destinations=(1,)):
    return (
        pinned(g, tunnels(g, list(sources), list(destinations), specs))
        .place(INF)
        .placed
    )


def test_lossless_bound_examples():
    # Duplicate routes, deduplicated routes, a static path filled, a zero first or last link, a
    # support's Pinned bound, weighted tunnels.
    assert bound(par2(4, 8), [(0, [0], 1), (0, [1], 1), (0, [1], 1)]) == 12
    assert bound(par2(4, 8), [(0, [0], 1), (0, [1], 1)]) == 8
    assert bound(par2(6, 4), [(0, [0], 1), (0, [1], 1)]) == 9
    assert bound(par2(0, 4), [(0, [0], 1), (0, [1], 1)]) == 0
    assert bound(par2(4, 0), [(0, [0], 1), (0, [1], 1)]) == 1
    assert bound(graph(4, PINLOSS), [(0, [0, 2], 1), (0, [1, 3], 1)], (0,), (3,)) == 4
    assert bound(par2(10, 10), [(0, [0], 1), (0, [1], 2)]) == 15


def test_duplicates_and_zero_weight_records():
    # Three records, entries 0, 1, 2; a weight-0 tunnel's record exists.
    g = par2(4, 8)
    p = pinned(g, tunnels(g, [0], [1], [(0, [0], 1), (0, [1], 1), (0, [1], 1)]))
    r = p.place(12)
    recs = p.records()
    assert [rec.entry for rec in recs] == [0, 1, 2]
    assert [rec.delivered for rec in recs] == [4, 4, 4]
    assert len({rec.flow_id for rec in recs}) == 3
    assert rows(r) == [(i, 0, UP, 4, 4, 0) for i in range(3)]
    z = pinned(g, tunnels(g, [0], [1], [(0, [0], 0), (0, [1], 1)]))
    r = z.place(INF)
    assert r.placed == 8
    assert [(rec.entry, rec.delivered) for rec in z.records()] == [(0, 0), (1, 8)]


def test_increments():
    # Increments, Lossless and Lossy.
    g = graph(3, INC)
    prog = tunnels(g, [0, 1], [2], [(0, [0], 1), (1, [1], 1)])
    p = pinned(g, prog)
    assert [p.place(3).placed, p.place(1).placed] == [2, 0]
    q = pinned(g, prog, LOSSY)
    r1 = q.place(3)
    assert (r1.placed, r1.call_drops) == (2, [(0, 1)])
    assert q.place(1).placed == 1
    # A Lossless increment's refused units.
    one = graph(2, ONE)
    ledger = ngc.FlowLedger(one)
    other = ngc.FlowPolicy(ledger, 0, 1, cost_only())
    other.place(8)
    p = pinned(one, tunnels(one, [0], [1], [(0, [0], 1)]), ledger=ledger)
    assert p.place(5).placed == 2
    other.remove()
    other.place(5)
    assert p.place(1).placed == 1


def test_infinite_and_mixed():
    g = par2(6, 4)
    prog = tunnels(g, [0], [1], [(0, [0], 1), (0, [1], 1)])
    r = pinned(g, prog).place(INF)
    assert (r.placed, r.termination, r.remaining) == (9, T.EXHAUSTED, INF)
    with pytest.raises(
        ValueError, match=exact("LOSSY forwarding requires a finite offer")
    ):
        pinned(g, prog, LOSSY).place(INF)
    h = par2(4, 8)
    mixed = ngc.ForwardingProgram.from_routes(
        h, [0], [1], [ngc.Route(0, path(h, [0])), ngc.Route(0, dag(h, [1]))]
    )
    p = pinned(h, mixed)
    assert p.place(INF).placed == 8
    assert [(rec.entry, rec.delivered) for rec in p.records()] == [(0, 4), (1, 4)]


def test_split_contexts():
    # Each identical split route splits in its own context, also pooled.
    g = par2(1, 1)
    d = dag(g, [0, 1])
    prog = ngc.ForwardingProgram.from_routes(
        g, [0], [1], [ngc.Route(0, d), ngc.Route(0, d)]
    )
    assert pinned(g, prog).place(INF).placed == 1
    r = pinned(g, prog, LOSSY).replace(2, outputs=ngc.OutputOptions(edge_flows=True))
    assert r.call_drops == [(0, 1)]
    pooled = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 1, cost_only())
    r = pooled.place(INF, ngc.OutputOptions(edge_flows=True))
    assert (r.placed, r.edge_flows.tolist()) == (2, [1, 1])


# ---------------------------------------------------------------- Y: losses inside one program


def lossy_rows(g, sources, destinations, specs, offer, **kw):
    r = pinned(g, tunnels(g, sources, destinations, specs), LOSSY, **kw).replace(offer)
    return r, [(t.assigned, t.delivered, t.dropped) for t in r.by_tunnel]


def test_loss_rule():
    one9 = graph(2, [(0, 1, 9)])
    r, got = lossy_rows(one9, [0], [1], [(0, [0], 2), (0, [0], 1)], 12)
    assert got == [(8, 6, 2), (4, 3, 1)]  # two entries merge on one link
    one1 = graph(2, [(0, 1, 1)])
    _, got = lossy_rows(one1, [0], [1], [(0, [0], 1), (0, [0], 1)], 2)
    assert got == [(1, 1, 0), (1, 0, 1)]  # a tie of the split keys
    one3 = graph(2, [(0, 1, 3)])
    _, got = lossy_rows(one3, [0], [1], [(0, [0], 5), (0, [0], 1), (0, [0], 1)], 7)
    assert got == [(5, 3, 2), (1, 0, 1), (1, 0, 1)]  # the rounding of the split
    hier = graph(4, HIER)
    _, got = lossy_rows(
        hier, [0, 1], [3], [(0, [0, 2], 1), (0, [0, 2], 1), (1, [1, 2], 1)], 4
    )
    assert [d for _, d, _ in got] == [1, 0, 0]  # the split hierarchy
    zero = graph(2, [(0, 1, 0)])
    r, got = lossy_rows(zero, [0], [1], [(0, [0], 1)], 2)
    assert got == [(2, 0, 2)] and r.call_drops == [(0, 2)]  # a link of service 0
    r, got = lossy_rows(
        graph(4, PINLOSS), [0], [3], [(0, [0, 2], 1), (0, [1, 3], 1)], 12
    )
    assert got == [(6, 2, 4), (6, 6, 0)] and r.call_drops == [(2, 4)]


def test_subdivision():
    # Subdividing a link before or after: a gate of 1, then three own links of 1.
    before = graph(
        6, [(0, 1, 1), (1, 2, 1), (1, 3, 1), (1, 4, 1), (2, 5, 1), (3, 5, 1), (4, 5, 1)]
    )
    _, got = lossy_rows(
        before,
        [0],
        [5],
        [(0, [0, 1, 4], 1), (0, [0, 2, 5], 1), (0, [0, 3, 6], 1)],
        3,
    )
    assert [d for _, d, _ in got] == [1, 0, 0]
    after = graph(
        7,
        [
            (0, 1, 1),
            (1, 2, 1),
            (1, 3, 1),
            (1, 4, 1),
            (2, 6, 1),
            (3, 5, 1),
            (4, 5, 1),
            (6, 5, 1),
        ],
    )
    _, got = lossy_rows(
        after,
        [0],
        [5],
        [(0, [0, 1, 4, 7], 1), (0, [0, 2, 5], 1), (0, [0, 3, 6], 1)],
        3,
    )
    assert [d for _, d, _ in got] == [1, 0, 0]


def test_cycles():
    # A compressed cycle: no link can drop.
    g = cyc(2, 2)
    _, got = lossy_rows(g, [4], [5], [(4, CYC_T0, 1), (4, CYC_T1, 1)], 2)
    assert [d for _, d, _ in got] == [1, 1]
    # The rounds from zero settle at round 3.
    g = cyc(1, 1)
    r, got = lossy_rows(g, [4], [5], [(4, CYC_T0, 1), (4, CYC_T1, 1)], 2)
    assert got == [(1, 1, 0), (1, 0, 1)] and r.call_drops == [(3, 1)]
    # Period two.
    g = cyc(1, 2)
    ledger = ngc.FlowLedger(g)
    p = pinned(
        g, tunnels(g, [4], [5], [(4, CYC_T0, 1), (4, CYC_T1, 2)]), LOSSY, ledger=ledger
    )
    before = (ledger.residual().tolist(), len(p.records()))
    with pytest.raises(
        ValueError, match=exact("lossy forwarding iteration entered a period-2 cycle")
    ):
        p.replace(3)
    assert (
        ledger.residual().tolist(),
        len(p.records()),
    ) == before and p.program is None
    # The round limit (64 * 3 rounds for three barriers).
    lim = graph(
        8,
        [
            (0, 1, 1),
            (1, 4, 1000),
            (2, 3, 2),
            (3, 0, 1000),
            (3, 7, 1000),
            (4, 5, 1),
            (5, 2, 1000),
            (5, 7, 1000),
            (6, 0, 1000),
            (6, 2, 1000),
        ],
    )
    prog = tunnels(
        lim, [6], [7], [(6, [8, 0, 1, 5, 6, 2, 4], 2), (6, [9, 2, 3, 0, 1, 5, 7], 3)]
    )
    with pytest.raises(
        ValueError,
        match=exact(
            "lossy forwarding iteration did not reach a fixed point within the round limit"
        ),
    ):
        pinned(lim, prog, LOSSY).replace(3)


def test_cycle_subdivision():
    # CYCS: Y split in two in series (e3: 2 -> 6, e8: 6 -> 3).
    def cycs(x, y, z):
        return graph(
            7,
            [
                (0, 1, x),
                (1, 2, 1000),
                (1, 5, 1000),
                (2, 6, y),
                (3, 0, 1000),
                (3, 5, 1000),
                (4, 0, 1000),
                (4, 2, 1000),
                (6, 3, z),
            ],
        )

    t0, t1 = [6, 0, 1, 3, 8, 5], [7, 3, 8, 4, 0, 2]
    r, got = lossy_rows(cycs(1, 1, 1), [4], [5], [(4, t0, 1), (4, t1, 1)], 2)
    assert got == [(1, 1, 0), (1, 0, 1)] and r.call_drops == [(3, 1)]
    g = cycs(1, 2, 2)
    with pytest.raises(
        ValueError, match=exact("lossy forwarding iteration entered a period-2 cycle")
    ):
        pinned(g, tunnels(g, [4], [5], [(4, t0, 1), (4, t1, 2)]), LOSSY).replace(3)


# ---------------------------------------------------------------- R: replay under a mask


def test_replay():
    # A Global program replayed with its first path failed.
    g = graph(4, G1U)
    prog = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 3, global_()).replace(12).program
    assert listing(prog) == [(0, [0, 2], 2), (0, [1, 3], 10)]
    fail = g.eligibility(excluded_edges=[2])
    p = pinned(g, prog, LOSSY, fail)
    r = p.replace(12)
    assert (r.placed, rows(r)) == (10, [(0, 0, FAILED, 0, 0, 0), (1, 0, UP, 12, 10, 2)])
    recs = p.records()
    assert recs[0].route_cost is None and recs[0].bundle.via_edges().tolist() == []
    assert pinned(g, prog, LOSSLESS, fail).place(INF).placed == 10
    zero = ngc.FlowLedger(g, [2, 0, 2, 10])
    r = pinned(g, prog, LOSSY, ledger=zero).replace(12)
    assert r.placed == 2 and [t.status for t in r.by_tunnel] == [UP, UP]
    # The surviving sources' captured weights.
    s = graph(4, SURV)
    surv = ngc.FlowPolicy(ngc.FlowLedger(s), [0, 1, 2], 3, global_()).place(9).program
    assert surv.source_weights == ((0, 3), (1, 1), (2, 5))
    q = pinned(s, surv, LOSSY, s.eligibility(excluded_nodes=[2]))
    r = q.replace(16)
    assert r.placed == 4 and [row.offered for row in r.by_source] == [12, 4, 0]
    assert (
        pinned(s, surv, LOSSLESS, s.eligibility(excluded_nodes=[2])).place(INF).placed
        == 4
    )


def test_replay_basis():
    # The topology basis keeps head 1; the supplied basis drops it.
    g = graph(3, RB)
    prog = ngc.FlowPolicy(ngc.FlowLedger(g), [0, 1], 2, global_()).place(10).program
    assert listing(prog) == [(0, [0], 5), (1, [1], 5)]
    mask = g.eligibility(excluded_edges=[1])
    r = pinned(g, prog, LOSSY, mask).replace(10)
    assert (r.placed, r.source_loss) == (5, 5)
    assert pinned(g, prog, LOSSLESS, mask).place(INF).placed == 1
    sup = tunnels(g, [0, 1], [2], [(0, [0], 1), (1, [1], 1)])
    r = pinned(g, sup, LOSSY, mask).replace(10)
    assert (r.placed, r.call_drops) == (5, [(0, 5)])
    assert pinned(g, sup, LOSSLESS, mask).place(INF).placed == 5
    # Every remaining weight 0.
    z = graph(3, [(0, 2, 5), (1, 2, 5)])
    zp = ngc.FlowPolicy(ngc.FlowLedger(z), [0, 1], 2, global_()).place(1).program
    m = z.eligibility(excluded_nodes=[0])
    r = pinned(z, zp, LOSSY, m).replace(4)
    assert (r.termination, r.source_loss, r.placed) == (T.NO_ACTIVE_WEIGHT, 4, 0)
    r = pinned(z, zp, LOSSLESS, m).replace(4)
    assert (r.termination, r.placed) == (T.NO_ACTIVE_WEIGHT, 0)


def test_replay_split_routes():
    # A stale remote tail repaired; a stale tunnel node repaired; frozen weights.
    g = graph(4, DIAM)
    prog = (
        ngc.FlowPolicy(ngc.FlowLedger(g), 0, 3, cost_only(admission=LOSSY))
        .place(1)
        .program
    )
    r = pinned(g, prog, LOSSY, g.eligibility(excluded_edges=[3])).replace(100)
    assert r.placed == 100
    t = graph(4, TUN)
    tp = tunnels(t, [0], [3], [(0, [0, 2], 1), (0, [1, 3], 1)])
    m = t.eligibility(excluded_nodes=[2])
    r = pinned(t, tp, LOSSY, m).replace(10)
    assert [s.status for s in r.by_tunnel] == [UP, FAILED] and r.placed == 10
    assert [s.assigned for s in r.by_tunnel] == [10, 0]
    assert pinned(t, tp, LOSSLESS, m).place(INF).placed == 10
    z = par2(0, 10)
    wp = ngc.FlowPolicy(ngc.FlowLedger(z), 0, 1, cost_only(W, LOSSY)).place(1).program
    fresh = par2(10, 10)
    r = pinned(z, wp, LOSSY, ledger=ngc.FlowLedger(z)).replace(5)
    assert (r.placed, r.call_drops) == (5, [])
    del fresh


# ---------------------------------------------------------------- S: stale forwarding


def test_stale_examples():
    g = graph(4, DIAM)
    p = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 3, cost_only(admission=LOSSY))
    p.replace(100)
    installed = p.program
    r = p.forward(100, black_holes=g.eligibility(excluded_edges=[3]))
    assert (r.placed, r.call_drops) == (50, [(3, 50)])  # a stale remote tail
    assert p.program is installed and not r.routes_changed and not r.weights_changed
    assert p.replace(100).placed == 100
    w = graph(2, [(0, 1, 10), (0, 1, 40)])
    q = ngc.FlowPolicy(ngc.FlowLedger(w), 0, 1, cost_only(W, LOSSY))
    q.replace(40)
    r = q.forward(40, black_holes=w.eligibility(excluded_edges=[0]))
    assert (r.placed, r.call_drops) == (32, [(0, 8)])  # a stale WCMP member
    wz = graph(2, [(0, 1, 0), (0, 1, 40)])
    q = ngc.FlowPolicy(ngc.FlowLedger(wz), 0, 1, cost_only(W, LOSSY))
    q.replace(1)
    assert q.forward(40).placed == 40  # a stale WCMP member of capacity 0
    rel = graph(4, REL)
    prog = tunnels(rel, [0, 1], [3], [(0, [0, 2], 1), (1, [1, 2], 1)])
    r = pinned(rel, prog, LOSSY).replace(20)
    assert [t.delivered for t in r.by_tunnel] == [5, 5]  # the relief's original program
    r = pinned(
        rel, prog, LOSSY, black_holes=rel.eligibility(excluded_edges=[0])
    ).replace(20)
    assert rows(r) == [
        (0, 0, BLACKHOLED, 10, 0, 10),
        (1, 1, UP, 10, 10, 0),
    ]  # no relief from a black hole


def test_stale_tunnels_and_pools():
    t = graph(4, TUN)
    tp = tunnels(t, [0], [3], [(0, [0, 2], 1), (0, [1, 3], 1)])
    bh = t.eligibility(excluded_nodes=[2])
    r = pinned(t, tp, LOSSY, black_holes=bh).replace(10)
    assert rows(r) == [(0, 0, UP, 5, 5, 0), (1, 0, BLACKHOLED, 5, 0, 5)]
    assert r.call_drops == [(1, 5)]
    assert pinned(t, tp, LOSSLESS, black_holes=bh).replace(10).placed == 1
    pool = graph(4, POOL2)
    p = ngc.FlowPolicy(ngc.FlowLedger(pool), [0, 1], 3, cost_only(admission=LOSSY))
    p.replace(0)
    r = p.forward(2, outputs=ngc.OutputOptions(edge_flows=True))
    assert (r.placed, r.edge_flows.tolist()) == (2, [1, 1, 1, 1])
    r = p.forward(2, black_holes=pool.eligibility(excluded_edges=[0, 3]))
    assert (r.placed, r.call_drops) == (1, [(0, 1)])  # a stale pooled remainder
    zero = graph(2, [(0, 1, 0)])
    q = ngc.FlowPolicy(ngc.FlowLedger(zero), 0, 1, cost_only(W, LOSSY))
    q.replace(0)
    r = q.forward(2)
    assert (r.source_loss, r.termination) == (2, T.NO_ACTIVE_WEIGHT)
    one = graph(2, ONE)
    q = ngc.FlowPolicy(ngc.FlowLedger(one), 0, 1, cost_only(admission=LOSSY))
    q.replace(0)
    r = q.forward(10, black_holes=one.eligibility(excluded_nodes=[0]))
    assert (r.source_loss, r.call_drops, r.placed) == (10, [], 0)
    op = tunnels(one, [0], [1], [(0, [0], 1)])
    r = pinned(one, op, LOSSY, black_holes=one.eligibility(excluded_nodes=[0])).replace(
        10
    )
    assert rows(r) == [(0, 0, BLACKHOLED, 0, 0, 0)] and r.source_loss == 10


def test_stale_with_mask_and_global():
    m = graph(5, MB)
    mp = tunnels(m, [0], [4], [(0, [0, 3], 1), (0, [1, 4], 1), (0, [2, 5], 1)])
    r = pinned(
        m,
        mp,
        LOSSY,
        m.eligibility(excluded_edges=[3]),
        m.eligibility(excluded_edges=[3, 4]),
    ).replace(10)
    assert rows(r) == [
        (0, 0, FAILED, 0, 0, 0),
        (1, 0, BLACKHOLED, 5, 0, 5),
        (2, 0, UP, 5, 5, 0),
    ]
    g = graph(4, G1U)
    prog = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 3, global_()).replace(12).program
    bh = g.eligibility(excluded_edges=[2])
    r = pinned(g, prog, LOSSY, black_holes=bh).replace(12)
    assert rows(r) == [(0, 0, BLACKHOLED, 2, 0, 2), (1, 0, UP, 10, 10, 0)]
    assert pinned(g, prog, LOSSLESS, black_holes=bh).replace(12).placed == 4
    allowed = g.eligibility()
    assert rows(pinned(g, prog, LOSSY, black_holes=allowed).replace(12)) == rows(
        pinned(g, prog, LOSSY).replace(12)
    )


# ---------------------------------------------------------------- O: operations


def test_evaluate_publishes_nothing():
    g = par2(10, 10)
    ledger = ngc.FlowLedger(g)
    p = ngc.FlowPolicy(ledger, 0, 1, cost_only())
    p.replace(6)

    def state():
        return ledger.residual().tolist(), [
            (x.flow_id, x.delivered, x.edges) for x in p.records()
        ]

    before = state()
    r = p.evaluate(12)
    assert (r.placed, r.delivery_delta) == (12, 6) and r.program is not None
    assert state() == before
    assert p.place(0).requested == 6
    gl = ngc.FlowPolicy(ngc.FlowLedger(graph(4, B1)), [0, 1], 3, global_())
    gl.place(2)
    out = ngc.OutputOptions(edge_flows=True, residual=True, min_cut=True)
    r = gl.evaluate(12, outputs=out)
    assert r.edge_flows.tolist() == [2, 10, 10] and r.residual.tolist() == [0, 0, 0]
    assert r.min_cut == [0, 1] and gl.placed_demand() == 2


def test_replace_and_epochs():
    g = par2(10, 10)
    ledger = ngc.FlowLedger(g)
    p = ngc.FlowPolicy(ledger, 0, 1, cost_only())
    p.place(4)
    old = p.records()[0].flow_id
    r = p.replace(6)
    assert (r.requested, r.admitted) == (6, 6)
    with pytest.raises(ValueError):
        ledger.record(old)
    assert p.records()[0].flow_id != old
    r = p.replace(INF)
    assert (r.placed, r.remaining, r.termination) == (20, INF, T.EXHAUSTED)
    q = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 1, cost_only())
    q.replace(10)
    r = q.replace(4)
    assert (r.placed, r.delivery_delta) == (4, -6)
    ids = [rec.flow_id for rec in q.records()]
    r = q.place(3)
    assert (r.placed, r.requested, r.admitted) == (3, 7, 7)
    assert [rec.flow_id for rec in q.records()] == ids
    lossy = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 1, cost_only(admission=LOSSY))
    lossy.replace(6)
    r = lossy.place(3)
    assert (r.placed, r.offered, r.requested) == (3, 9, 9)


def test_forward_keeps_the_program_and_flags():
    g = par2(10, 10)
    ledger = ngc.FlowLedger(g, [2, 10])
    a = ngc.FlowPolicy(
        ledger, 0, 1, cost_only(W, weights=ngc.WeightBasis.RESIDUAL_SNAPSHOT)
    )
    assert a.replace(6).placed == 6
    b = pinned(g, tunnels(g, [0], [1], [(0, [1], 1)]), ledger=ledger)
    b.place(5)
    r = a.forward(12)
    assert (r.placed, r.routes_changed, r.weights_changed) == (6, False, False)
    r = a.replace(12)
    assert (r.placed, r.routes_changed, r.weights_changed) == (7, False, True)
    r = a.replace(12)
    assert (r.placed, r.routes_changed, r.weights_changed) == (7, False, False)


def test_forward_texts_and_order():
    g = par2(10, 10)
    p = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 1, cost_only(admission=LOSSY))
    with pytest.raises(
        ValueError, match=exact("forward requires an installed program")
    ):
        p.forward(-1)
    p.place(1)
    with pytest.raises(
        ValueError,
        match=exact("offer must be an integer in [0, 2^62) or permitted +inf"),
    ):
        p.forward(-1)
    with pytest.raises(
        ValueError, match=exact("LOSSY forwarding requires a finite offer")
    ):
        p.forward(INF)
    other = par2(1, 1)
    with pytest.raises(ValueError, match=exact("black_holes belongs to another graph")):
        p.forward(1, black_holes=other.eligibility(excluded_edges=[0]))
    for call in (p.replace, p.evaluate):
        with pytest.raises(
            ValueError, match=exact("LOSSY forwarding requires a finite offer")
        ):
            call(INF)
    p.remove()
    with pytest.raises(
        ValueError, match=exact("forward requires an installed program")
    ):
        p.forward(1)
    ca = ngc.FlowPolicy(
        ngc.FlowLedger(g),
        0,
        1,
        ngc.FlowPolicyConfig(routing=ngc.CapacityAware(ngc.Tunnels(1))),
    )
    with pytest.raises(
        ValueError, match=exact("forward requires an installed program")
    ):
        ca.forward(1)
    pin = pinned(g, tunnels(g, [0], [1], [(0, [0], 1)]))
    with pytest.raises(ValueError, match=exact("Pinned forwards with replace")):
        pin.forward(INF, black_holes=other.eligibility(excluded_edges=[0]))


def test_zero_replace_and_first_flags():
    g = graph(4, POOL2)
    p = ngc.FlowPolicy(ngc.FlowLedger(g), [0, 1], 3, cost_only())
    r = p.replace(0)
    assert r.termination == T.COMPLETE and r.routes_changed and r.weights_changed
    assert [rec.delivered for rec in p.records()] == [0, 0]
    q = ngc.FlowPolicy(ngc.FlowLedger(g), [0, 1], 3, cost_only())
    r = q.place(2)
    assert r.routes_changed and r.weights_changed
    r = q.place(2)
    assert not r.routes_changed and not r.weights_changed


def test_controller_flags():
    # A controller's increment changes the weights; an offer-only change sets no flag.
    g = graph(2, [(0, 1, 20)])
    p = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 1, global_())
    r1 = p.place(5)
    assert (r1.routes_changed, r1.weights_changed) == (True, True)
    r2 = p.place(5)
    assert (r2.routes_changed, r2.weights_changed) == (False, True)
    r3 = p.replace(10)
    assert (r3.routes_changed, r3.weights_changed) == (False, False)
    b1 = graph(4, B1)
    q = ngc.FlowPolicy(ngc.FlowLedger(b1), [0, 1], 3, global_())
    assert listing(q.place(2).program) == [(0, [0], 2)]
    r = q.place(10)
    assert listing(r.program) == [(0, [0], 2), (1, [1, 2], 10)] and r.routes_changed
    pin = pinned(g, tunnels(g, [0], [1], [(0, [0], 1)]))
    flags = [
        (r.routes_changed, r.weights_changed) for r in map(pin.replace, (5, 10, 10))
    ]
    assert flags[1:] == [(False, False), (False, False)]
    cost = ngc.FlowPolicy(
        ngc.FlowLedger(par2(10, 10)), 0, 1, cost_only(admission=LOSSY)
    )
    cost.place(1)
    assert [cost.forward(x).routes_changed for x in (5, 7)] == [False, False]
    gu = graph(4, G1U)
    c = ngc.FlowPolicy(ngc.FlowLedger(gu), 0, 3, global_())
    c.replace(12)
    assert c.forward(12).placed == 12
    r = c.forward(20)
    assert (r.placed, r.termination) == (12, T.EXHAUSTED)


# ---------------------------------------------------------------- D, Z: records and rows


def test_pooled_records_in_pinned():
    # A CostOnly program in Pinned keeps the attribution of 5.1.
    g = graph(4, POOL2)
    prog = (
        ngc.FlowPolicy(ngc.FlowLedger(g), [0, 1], 3, cost_only(admission=LOSSY))
        .place(2)
        .program
    )
    p = pinned(g, prog, LOSSY)
    p.replace(2)
    assert [(rec.entry, rec.edges) for rec in p.records()] == [
        (0, [(0, 1), (2, 1)]),
        (1, [(1, 1), (3, 1)]),
    ]


def test_wide_rows():
    # TunnelFlow.assigned and dropped above 2^64.
    g = graph(2, [(0, 1, 0)])
    p = pinned(g, tunnels(g, [0], [1], [(0, [0], 1)]), LOSSY)
    big = (1 << 62) - 1
    for _ in range(5):
        r = p.place(big)
    total = 5 * big
    assert total > 1 << 64
    assert rows(r) == [(0, 0, UP, total, 0, total)]
    assert r.offered == total and r.drops == [(0, total)] and r.source_loss == 0
    rec = p.records()[0]
    assert rec.offered == total and rec.drops == [(0, total)]


# ---------------------------------------------------------------- E: evaluators


def test_query_options_offered():
    g = par2(6, 4)
    prog = tunnels(g, [0], [1], [(0, [0], 1), (0, [1], 1)])
    cfg = ngc.FlowPolicyConfig(routing=ngc.Pinned(prog))
    assert ngc.QueryOptions(cfg, offered=5).offered == 5
    assert ngc.QueryOptions(cfg).offered is None
    assert ngc.QueryOptions(cfg, offered=INF).offered == INF
    r = ngc.capacity(g, 0, 1, ngc.QueryOptions(cfg))
    assert (r.placed, r.termination) == (9, T.EXHAUSTED)
    pl = graph(4, PINLOSS)
    lp = tunnels(pl, [0], [3], [(0, [0, 2], 1), (0, [1, 3], 1)])
    lcfg = ngc.FlowPolicyConfig(routing=ngc.Pinned(lp), admission=LOSSY)
    assert ngc.capacity(pl, 0, 3, ngc.QueryOptions(lcfg, offered=12)).placed == 8
    texts = [
        (ngc.QueryOptions(lcfg), "Pinned LOSSY capacity requires finite offered"),
        (
            ngc.QueryOptions(lcfg, offered=INF),
            "Pinned LOSSY capacity requires finite offered",
        ),
        (
            ngc.QueryOptions(cost_only(admission=LOSSY), offered=5),
            "CostOnly LOSSY capacity uses first-hop injection",
        ),
        (
            ngc.QueryOptions(cost_only(), offered=5),
            "CostOnly LOSSLESS capacity does not accept offered",
        ),
        (
            ngc.QueryOptions(global_(), offered=5),
            "Global LOSSLESS capacity does not accept offered",
        ),
        (
            ngc.QueryOptions(cfg, offered=5),
            "Pinned LOSSLESS capacity does not accept offered",
        ),
        (
            ngc.QueryOptions(
                ngc.FlowPolicyConfig(routing=ngc.CapacityAware(ngc.Tunnels(1))),
                design_demand=4,
                offered=5,
            ),
            "T1 LOSSLESS capacity does not accept offered",
        ),
    ]
    graphs = {id(lcfg): pl}
    for opts, text in texts:
        gg = graphs.get(
            id(opts.config),
            pl
            if opts.config.routing.__class__ is ngc.Pinned
            and opts.config.admission == LOSSY
            else g,
        )
        with pytest.raises(ValueError, match=exact(text)):
            ngc.capacity(gg, 0, 3 if gg is pl else 1, opts)
    with pytest.raises(
        ValueError,
        match=exact("offered must be an integer in [0, 2^62) or permitted +inf"),
    ):
        ngc.QueryOptions(cfg, offered=-1)
    with pytest.raises(
        ValueError, match=exact("batch_capacity does not accept Pinned")
    ):
        ngc.batch_capacity(g, [([0], [1])], ngc.QueryOptions(cfg))


def test_pinned_sensitivity():
    pl = graph(4, PINLOSS)
    lp = tunnels(pl, [0], [3], [(0, [0, 2], 1), (0, [1, 3], 1)])
    cfg = ngc.FlowPolicyConfig(routing=ngc.Pinned(lp))
    s = ngc.sensitivity(pl, 0, 3, ngc.QueryOptions(cfg), workers=1)
    assert s.baseline.placed == 4
    assert [(e.edge, e.loss) for e in s.entries] == [(0, -6), (1, 2), (2, -6), (3, 2)]
    lcfg = ngc.FlowPolicyConfig(routing=ngc.Pinned(lp), admission=LOSSY)
    s = ngc.sensitivity(pl, 0, 3, ngc.QueryOptions(lcfg, offered=12), workers=1)
    assert s.baseline.placed == 8
    assert [(e.edge, e.loss) for e in s.entries] == [(0, -2), (1, 6), (2, -2), (3, 6)]
    s2 = ngc.sensitivity(pl, 0, 3, ngc.QueryOptions(lcfg, offered=12), workers=2)
    assert [(e.edge, e.loss) for e in s2.entries] == [
        (e.edge, e.loss) for e in s.entries
    ]
    es = graph(4, ES)
    sup = ngc.sensitivity(
        es, 0, 3, ngc.QueryOptions(global_(support=dag(es, [1, 3]))), workers=1
    )
    assert sup.baseline.placed == 10
    assert [(e.edge, e.loss) for e in sup.entries] == [(0, 0), (1, 10), (2, 0), (3, 10)]


# ---------------------------------------------------------------- properties (program_oracle.py)


def random_dag(rng, n):
    """A random DAG on n nodes (edges from lower to higher node, some parallel, capacities with
    zeros, costs 1 to 4), its rows sorted into EdgeId order."""
    rows = []
    for u in range(n):
        for v in range(u + 1, n):
            for _ in range(rng.choice((0, 1, 1, 2))):
                rows.append((u, v, rng.choice((0, 3, 5, 7, 10, 12)), rng.randint(1, 4)))
    rows.sort(key=lambda r: (r[0], r[3], r[1]))
    return rows


def simple_paths(rows, src, dst):
    """Every path (EdgeIds) from src to dst over the DAG rows."""
    out, stack = [], [(src, [])]
    while stack:
        u, p = stack.pop()
        if u == dst:
            out.append(p)
            continue
        for e, r in enumerate(rows):
            if r[0] == u:
                stack.append((r[1], [*p, e]))
    return sorted(out)


def random_tunnel_program(rng):
    """A random graph and a supplied tunnel program over it (sources {0, 1} or {0}, destination
    the last node), with weights 0 to 3; None when no source has a path. A quarter are the
    crossing gadget: head 0's tunnels 0, 3, 5 and 1, 4, 6 cross at node 2, so a
    mask can fail both while their union still reaches the destination (1, 3, 5)."""
    if rng.random() < 0.25:
        caps = [rng.choice((1, 3, 10)) for _ in range(7)]
        ends = [(0, 2), (0, 2), (1, 5), (2, 3), (2, 4), (3, 5), (4, 5)]
        rows = [(u, v, c, 1) for (u, v), c in zip(ends, caps, strict=True)]
        specs = [
            (0, [0, 3, 5], rng.randint(0, 3)),
            (0, [1, 4, 6], rng.randint(1, 3)),
            (1, [2], rng.randint(0, 3)),
        ]
        return 6, rows, [0, 1], specs
    n = rng.randint(4, 6)
    rows = random_dag(rng, n)
    if (
        rng.random() < 0.5
    ):  # two heads into one tight link: ties between heads at a split
        rows = [(0, 2, 10, 1), (1, 2, 10, 1)] + [r for r in rows if r[0] >= 2]
        rows += [(2, n - 1, rng.choice((1, 2, 3)), 1)]
        rows.sort(key=lambda r: (r[0], r[3], r[1]))
    sources = [0, 1] if rng.random() < 0.7 else [0]
    specs = []
    for s in sources:
        candidates = simple_paths(rows, s, n - 1)
        for p in rng.sample(
            candidates,
            min(len(candidates), rng.randint(0, 3) if s else rng.randint(1, 5)),
        ):
            specs.append((s, p, rng.randint(0, 3)))
    if not specs:
        return None
    rng.shuffle(specs)
    return n, rows, sources, specs


def test_tunnel_programs_against_the_oracle():
    # Random tunnel programs, with and without a failure mask: the by_tunnel amounts (8.1), the
    # LOSSY rows and edge loads (9.1, 9.2, nested split at each link), conservation, and the
    # LOSSLESS bound with its maximality (8.3).
    import random

    import program_oracle as orc

    rng = random.Random(20261018)
    checked = masked_checked = kept_without_entry = 0
    for _ in range(150):
        made = random_tunnel_program(rng)
        if made is None:
            continue
        n, rows_, sources, specs = made
        g = graph(n, rows_)
        caps = [r[2] for r in rows_]
        prog = tunnels(g, sources, [n - 1], specs)
        excluded = (
            []
            if rng.random() < 0.5
            else rng.sample(range(len(rows_)), min(len(rows_), rng.randint(1, 3)))
        )
        elig = g.eligibility(excluded_edges=excluded) if excluded else None
        alive = [not (set(es) & set(excluded)) for _, es, _ in specs]
        roster = orc.replay_roster(specs, rows_, excluded, [n - 1])
        kept_without_entry += any(
            not any(a for (h, _, _), a in zip(specs, alive, strict=True) if h == x)
            for x in roster
        )
        offer = rng.randint(0, 40)
        amounts, lost = orc.entry_amounts(offer, specs, alive, roster=roster)
        r = pinned(g, prog, LOSSY, elig=elig).replace(
            offer, outputs=ngc.OutputOptions(edge_flows=True)
        )
        assert [t.index for t in r.by_tunnel] == list(range(len(specs)))  # list order
        assert [t.assigned for t in r.by_tunnel] == amounts
        assert [t.status for t in r.by_tunnel] == [UP if a else FAILED for a in alive]
        want_rows, want_edges = orc.lossy(specs, amounts, caps)
        assert [(t.delivered, t.dropped) for t in r.by_tunnel] == want_rows
        assert r.edge_flows.tolist() == want_edges
        assert r.source_loss == lost
        assert (
            r.offered
            == offer
            == r.delivered + sum(d for _, d in r.drops) + r.source_loss
        )
        b = orc.lossless_bound(specs, caps, alive, roster=roster)
        lossless = pinned(g, prog, LOSSLESS, elig=elig)
        got = lossless.replace(INF, outputs=ngc.OutputOptions(edge_flows=True))
        assert got.placed == b
        assert all(f <= c for f, c in zip(got.edge_flows.tolist(), caps, strict=True))
        if b < 40:
            assert lossless.replace(b + 1).placed == b  # b + 1 does not fit
        checked += 1
        masked_checked += bool(excluded)
    # Not vacuous: masks, and heads the replay keeps without a surviving tunnel.
    assert checked > 100 and masked_checked > 30 and kept_without_entry > 0, (
        kept_without_entry
    )


def test_export_replay_property():
    # Over random graphs: a Global program (both origins) forwarded by Pinned on a fresh
    # ledger places the same amount on the same edges; its rows split each head's share over
    # its tunnels exactly by the tunnel weights.
    import random

    import program_oracle as orc

    rng = random.Random(7)
    checked = 0
    for _ in range(80):
        n = rng.randint(4, 6)
        rows_ = random_dag(rng, n)
        g = graph(n, rows_)
        sources = [0, 1] if rng.random() < 0.5 else [0]
        origin = rng.choice((ngc.Origin.OPTIMIZED, ngc.Origin.EQUAL))
        offer = rng.randint(1, 30)
        r = ngc.FlowPolicy(ngc.FlowLedger(g), sources, n - 1, global_(origin)).place(
            offer, ngc.OutputOptions(edge_flows=True)
        )
        if r.placed == 0:
            continue
        prog = r.program
        replay = pinned(g, prog).place(r.placed, ngc.OutputOptions(edge_flows=True))
        assert (replay.placed, replay.termination) == (r.placed, T.COMPLETE)
        assert replay.edge_flows.tolist() == r.edge_flows.tolist()
        specs = [(t.head, list(t.edges), t.weight) for t in prog.tunnels]
        amounts, lost = orc.entry_amounts(
            r.placed, specs, source_weights=prog.source_weights
        )
        assert lost == 0 and amounts == [t.weight for t in prog.tunnels]
        assert [t.assigned for t in replay.by_tunnel] == amounts
        checked += 1
    assert checked > 50


def test_global_support_optimality_certificate():
    # Global(support) over random supports: the value is the support's maximum flow (an
    # independent Edmonds-Karp), the flow uses only the support, fits and is conserved, and the
    # optimized origin leaves no negative-cost residual cycle (a minimum-cost maximum flow).
    import random

    import program_oracle as orc

    rng = random.Random(11)
    checked = 0
    for _ in range(80):
        n = rng.randint(4, 6)
        rows_ = random_dag(rng, n)
        if not rows_:
            continue
        g = graph(n, rows_)
        support = sorted(rng.sample(range(len(rows_)), rng.randint(1, len(rows_))))
        sources = [0, 1] if rng.random() < 0.5 else [0]
        r = ngc.FlowPolicy(
            ngc.FlowLedger(g), sources, n - 1, global_(support=dag(g, support))
        ).place(INF, ngc.OutputOptions(edge_flows=True))
        flow = r.edge_flows.tolist()
        arcs = [(rows_[e][0], rows_[e][1], rows_[e][2]) for e in support]
        value, _ = orc.max_flow(n, arcs, sources, [n - 1])
        assert r.placed == value
        assert all(f == 0 for e, f in enumerate(flow) if e not in support)
        assert all(f <= rows_[e][2] for e, f in enumerate(flow))
        for v in range(n):
            if v in sources or v == n - 1:
                continue
            assert sum(f for e, f in enumerate(flow) if rows_[e][1] == v) == sum(
                f for e, f in enumerate(flow) if rows_[e][0] == v
            )
        maximal, min_cost = orc.residual_certificates(
            n,
            [(*rows_[e][:3], rows_[e][3]) for e in support],
            [flow[e] for e in support],
            sources,
            [n - 1],
        )
        assert maximal and min_cost
        assert r.total_cost == sum(flow[e] * rows_[e][3] for e in support)
        checked += 1
    assert checked > 60


def test_mapping():
    # today's Pinned uses mapped onto programs (values derived by hand).
    split = ngc.FlowPlacement.WEIGHTED
    # One EqualBalanced DAG route over DIAM: Lossless 200; Lossy place(100) delivers 100.
    g = graph(4, DIAM)
    one = ngc.ForwardingProgram.from_routes(
        g, 0, 3, [ngc.Route(0, dag(g, [0, 1, 2, 3]))]
    )
    assert pinned(g, one).place(INF).placed == 200
    assert pinned(g, one, LOSSY).place(100).placed == 100
    # Weighted with capacity weights (the default basis): 200; with the residual snapshot
    # (100, 50, 100, 50), node 0 splits by H(x, (100, 50)): H(150) = (100, 50) fits, H(151) =
    # (101, 50) does not.
    w = ngc.ForwardingProgram.from_routes(
        g, 0, 3, [ngc.Route(0, dag(g, [0, 1, 2, 3]))], placement=split
    )
    assert pinned(g, w).place(INF).placed == 200
    snap = ngc.ForwardingProgram.from_routes(
        g,
        0,
        3,
        [ngc.Route(0, dag(g, [0, 1, 2, 3]))],
        placement=split,
        weights=ngc.WeightBasis.RESIDUAL_SNAPSHOT,
        residual=[100, 50, 100, 50],
    )
    assert pinned(g, snap).place(INF).placed == 150
    # Two one-edge DAG routes on PAR2(6, 4): EqualBalanced and Weighted both split by round
    # robin over the entries: 9; Global(support) over both edges: 10.
    g = par2(6, 4)
    for options in ({}, {"placement": split}):
        two = ngc.ForwardingProgram.from_routes(
            g, 0, 1, [ngc.Route(0, dag(g, [0])), ngc.Route(0, dag(g, [1]))], **options
        )
        assert pinned(g, two).place(INF).placed == 9
    assert (
        ngc.FlowPolicy(ngc.FlowLedger(g), 0, 1, global_(support=dag(g, [0, 1])))
        .place(INF)
        .placed
        == 10
    )
    # Two routes sharing e0 (4 units) on SHARE, Lossy place(10): RR(10, 2) = (5, 5) arrive at
    # e0, which carries H(4, (5, 5)) = (2, 2): records 2 and 2.
    g = graph(3, SHARE)
    shared = ngc.ForwardingProgram.from_routes(
        g, 0, 2, [ngc.Route(0, dag(g, [0, 1])), ngc.Route(0, dag(g, [0, 2]))]
    )
    p = pinned(g, shared, LOSSY)
    assert p.place(10).placed == 4
    assert [rec.delivered for rec in p.records()] == [2, 2]
    # One MaxFlow bundle becomes Global(support): ES over all four links, 12 at cost 34, one
    # record per source with the used edges and no entry.
    g = graph(4, ES)
    p = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 3, global_(support=dag(g, [0, 1, 2, 3])))
    r = p.place(INF)
    assert (r.placed, r.total_cost) == (12, 34)
    assert [(rec.entry, list(rec.bundle.via_edges())) for rec in p.records()] == [
        (None, [0, 1, 2, 3])
    ]
    # Lossy +inf is refused.
    g = par2(6, 4)
    two = tunnels(g, [0], [1], [(0, [0], 1), (0, [1], 1)])
    with pytest.raises(
        ValueError, match=exact("LOSSY forwarding requires a finite offer")
    ):
        pinned(g, two, LOSSY).place(INF)
    # Several sources and destinations are accepted.
    g = graph(4, B1)
    both = tunnels(g, [0, 1], [3], [(0, [0], 1), (1, [1, 2], 1)])
    assert pinned(g, both).place(4).placed == 4
    # Duplicate routes are distinct entries: 12 over three records.
    g = par2(4, 8)
    dup = tunnels(g, [0], [1], [(0, [0], 1), (0, [1], 1), (0, [1], 1)])
    p = pinned(g, dup)
    assert p.place(INF).placed == 12 and [rec.entry for rec in p.records()] == [0, 1, 2]
    # An empty route list.
    with pytest.raises(ValueError, match=exact("routes must not be empty")):
        ngc.ForwardingProgram.from_routes(g, 0, 1, [])
    # Pinned takes a program: no bundles keyword, a list is not a program, None is refused.
    with pytest.raises(TypeError):
        ngc.Pinned(bundles=[dag(g, [0])])
    with pytest.raises(TypeError):
        ngc.Pinned([dag(g, [0])])
    with pytest.raises(
        TypeError, match=exact("Pinned program must be a ForwardingProgram")
    ):
        ngc.Pinned(None)


# ---------------------------------------------------------------- facts of the retired fixtures


def test_pinned_facts_of_the_retired_reference_fixtures():
    # The Pinned facts that the removed reference fixtures held, derived again from the rules.
    # pinned_bundle_rules: a split route's record reports its cheapest walk as route_cost: the
    # DAG {0 -> 1 (cost 1), 1 -> 2 (cost 2), 0 -> 2 (cost 4)} has walks of cost 3 and 4.
    g = graph(3, [(0, 1, 5, 1), (0, 2, 5, 4), (1, 2, 5, 2)])
    prog = ngc.ForwardingProgram.from_routes(
        g, [0], [2], [ngc.Route(0, dag(g, [0, 1, 2]))]
    )
    p = pinned(g, prog)
    p.place(1)
    assert [(rec.entry, rec.route_cost) for rec in p.records()] == [(0, 3)]
    # eb_1_plus_2_vs_3: Lossless Pinned is additive: place(1); place(2) re-solves 3 from zero,
    # RR(3, 3) = (1, 1, 1), as place(3) does.
    g = graph(2, [(0, 1, 10), (0, 1, 10), (0, 1, 10)])
    prog = tunnels(g, [0], [1], [(0, [e], 1) for e in range(3)])
    flows = ngc.OutputOptions(edge_flows=True)
    split = pinned(g, prog)
    split.place(1)
    two = split.place(2, flows)
    one = pinned(g, prog).place(3, flows)
    assert two.edge_flows.tolist() == one.edge_flows.tolist() == [1, 1, 1]
    assert (two.placed, one.placed) == (2, 3)
    # admission_zero_and_no_route: every tunnel failed: the head reaches no destination over its
    # own routes, the replay roster is empty, a positive request ends NO_ROUTE.
    g = graph(4, TUN)
    prog = tunnels(g, [0], [3], [(0, [0, 2], 1), (0, [1, 3], 1)])
    r = pinned(g, prog, elig=g.eligibility(excluded_edges=[0, 1])).place(5)
    assert (r.placed, r.termination) == (0, T.NO_ROUTE)
    assert [t.status for t in r.by_tunnel] == [FAILED, FAILED]
    # flag_zero_delivery_creation: a first place(1) on a tunnel whose only unit another owner
    # holds creates the record at zero delivery and installs the program (both flags); a second
    # zero-delivery call keeps the id and sets no flag (an offer-only call).
    g = graph(2, [(0, 1, 1)])
    ledger = ngc.FlowLedger(g)
    ngc.FlowPolicy(ledger, 0, 1, cost_only()).place(1)
    p = pinned(g, tunnels(g, [0], [1], [(0, [0], 1)]), ledger=ledger)
    r = p.place(1)
    ids = [rec.flow_id for rec in p.records()]
    assert (r.placed, r.termination, r.routes_changed, r.weights_changed) == (
        0,
        T.EXHAUSTED,
        True,
        True,
    )
    assert [rec.delivered for rec in p.records()] == [0]
    r = p.place(1)
    assert (r.placed, r.routes_changed, r.weights_changed) == (0, False, False)
    assert [rec.flow_id for rec in p.records()] == ids


# ---------------------------------------------------------------- Python allocation failures


@needs_hooks
@pytest.mark.parametrize(
    "case",
    ["costonly replace", "costonly forward", "global evaluate", "pinned replace"],
)
def test_operations_are_all_or_nothing(case):
    # Every Python object of evaluate, replace and forward exists before anything is published:
    # the k-th Python allocation failing raises with the ledger and the policy as they were (the
    # installed program and the earlier call's records included), and a retry returns what a
    # call without failure returns.
    g = graph(4, PINLOSS)
    prog = tunnels(g, [0], [3], [(0, [0, 2], 1), (0, [1, 3], 1)])
    configs = {
        "costonly replace": (cost_only(admission=LOSSY), "replace"),
        "costonly forward": (cost_only(admission=LOSSY), "forward"),
        "global evaluate": (global_(), "evaluate"),
        "pinned replace": (
            ngc.FlowPolicyConfig(routing=ngc.Pinned(prog), admission=LOSSY),
            "replace",
        ),
    }
    cfg, method = configs[case]

    def make():
        ledger = ngc.FlowLedger(g)
        policy = ngc.FlowPolicy(ledger, 0, 3, cfg)
        policy.replace(2)

        def snapshot():
            return [
                full_ledger_state(ledger),
                repr(policy._state()),
                listing(policy.program),
            ]

        return getattr(policy, method), (6,), snapshot

    assert sweep_python_allocation(make) >= 3  # not vacuous


@needs_hooks
def test_from_routes_and_getters_survive_python_allocation_failures():
    # from_routes and the getters of shared values make their objects through the checked path:
    # a failed Python allocation raises, and the call that completes returns the clean value.
    g = graph(2, [(0, 1, 4), (0, 1, 8)])
    routes = [ngc.Route(0, path(g, [0])), ngc.Route(0, dag(g, [0, 1]), weight=2)]
    prog = ngc.ForwardingProgram.from_routes(g, [0], [1], routes)
    # A Global result: its program has no Python object yet, so the getter makes one.
    result = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 1, global_()).place(3)
    calls = [
        (ngc.ForwardingProgram.from_routes, (g, [0], [1], routes)),
        (lambda: result.program, ()),
        (lambda: result.by_tunnel, ()),
        (lambda: prog.tunnels, ()),
    ]
    for fn, args in calls:
        failures = sweep_python_allocation(lambda fn=fn, args=args: (fn, args, list))
        assert failures >= 1, fn


# ---------------------------------------------------------------- wide weights, zero offers, failed entries


@pytest.mark.parametrize("w", [1 << 61, (1 << 62) - 1])
@pytest.mark.parametrize("admission", [LOSSLESS, LOSSY])
@pytest.mark.parametrize("kind", ["tunnels", "explicit"])
def test_large_equal_weights_deliver_exactly_eight(kind, admission, w):
    # Eight equal legal weights sum to 2^64 or more; equal weights give RR(8, 8): one unit
    # on each of eight links of 100.
    g = graph(2, [(0, 1, 100)] * 8)
    if kind == "tunnels":
        p = pinned(g, tunnels(g, [0], [1], [(0, [i], w) for i in range(8)]), admission)
    else:
        p = ngc.FlowPolicy(
            ngc.FlowLedger(g),
            0,
            1,
            cost_only(
                W, admission, weights=ngc.WeightBasis.EXPLICIT, edge_weights=[w] * 8
            ),
        )
    r = p.replace(8, outputs=ngc.OutputOptions(edge_flows=True))
    assert r.placed == r.delivered == r.admitted == 8
    assert r.edge_flows.tolist() == [1] * 8
    assert r.source_loss == r.call_source_loss == 0


@pytest.mark.parametrize("admission", [LOSSLESS, LOSSY])
def test_large_unequal_weights_split_by_highest_averages(admission):
    # Weights 2^62 - 1, 2^61, 2^61 + 1, 1 (W = 2^63 + 1) on four links of 100, offer 9:
    # lower quotas floor(9 w / W) = (4, 2, 2, 0); the ninth unit goes to the largest w / (a + 1),
    # (2^62 - 1) / 5 against (2^61 + 1) / 3: H(9, w) = (5, 2, 2, 0), the oracle's too; explicit
    # WCMP and tunnels alike.
    import program_oracle as orc

    weights = [(1 << 62) - 1, 1 << 61, (1 << 61) + 1, 1]
    want = orc.ha(9, weights)
    assert want == [5, 2, 2, 0]
    g = graph(2, [(0, 1, 100)] * 4)
    for p in (
        pinned(
            g,
            tunnels(g, [0], [1], [(0, [i], w) for i, w in enumerate(weights)]),
            admission,
        ),
        ngc.FlowPolicy(
            ngc.FlowLedger(g),
            0,
            1,
            cost_only(
                W, admission, weights=ngc.WeightBasis.EXPLICIT, edge_weights=weights
            ),
        ),
    ):
        r = p.replace(9, outputs=ngc.OutputOptions(edge_flows=True))
        assert r.placed == 9 and r.edge_flows.tolist() == want and r.source_loss == 0


@pytest.mark.parametrize("method", ["replace", "evaluate", "forward"])
@pytest.mark.parametrize("case", ["no-route", "zero-weight"])
@pytest.mark.parametrize("admission", [LOSSLESS, LOSSY])
def test_full_zero_offer_is_complete(method, case, admission):
    # A full zero offer is COMPLETE (14.2's first row), after its preparation: an empty
    # roster and a program without a positive weight alike.
    g = graph(2, [(0, 1, 10)])
    cfg = cost_only(W, admission, weights=ngc.WeightBasis.EXPLICIT, edge_weights=[0])
    elig = g.eligibility(excluded_edges=[0]) if case == "no-route" else None
    p = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 1, cfg, elig)
    if method == "forward":
        p.replace(1)
    r = getattr(p, method)(0)
    assert r.termination == T.COMPLETE and r.placed == 0
    assert r.program is not None  # prepared


@pytest.mark.parametrize("method", ["replace", "evaluate"])
@pytest.mark.parametrize("admission", [LOSSLESS, LOSSY])
def test_full_zero_offer_is_complete_for_pinned(method, admission):
    # The same for Pinned: a program whose head is masked out, and one whose weights are all 0.
    g = graph(2, [(0, 1, 10)])
    masked = pinned(
        g,
        tunnels(g, [0], [1], [(0, [0], 1)]),
        admission,
        elig=g.eligibility(excluded_edges=[0]),
    )
    zero = pinned(g, tunnels(g, [0], [1], [(0, [0], 0)]), admission)
    for p in (masked, zero):
        assert getattr(p, method)(0).termination == T.COMPLETE


@pytest.mark.parametrize("admission", [LOSSLESS, LOSSY])
@pytest.mark.parametrize("offer", [0, 5])
def test_all_failed_entries_still_have_records(admission, offer):
    # One record per program entry also when every entry fails: zero, an empty bundle,
    # route cost None; replace retires and recreates them, an increment keeps them, evaluate
    # publishes nothing.
    g = graph(2, [(0, 1, 10)])
    p = pinned(
        g,
        tunnels(g, [0], [1], [(0, [0], 1)]),
        admission,
        elig=g.eligibility(excluded_edges=[0]),
    )
    assert p.evaluate(offer).placed == 0 and p.records() == []
    r = p.replace(offer)
    assert r.termination == (T.COMPLETE if offer == 0 else T.NO_ROUTE)
    rows = lambda: [  # noqa: E731
        (
            rec.entry,
            rec.delivered,
            rec.edges,
            rec.bundle.via_edges().tolist(),
            rec.route_cost,
        )
        for rec in p.records()
    ]
    assert rows() == [(0, 0, [], [], None)]
    first = [rec.flow_id for rec in p.records()]
    p.replace(offer)
    second = [rec.flow_id for rec in p.records()]
    assert rows() == [(0, 0, [], [], None)] and second != first  # retired and recreated
    r = p.place(3)
    assert (
        r.placed == 0 and [rec.flow_id for rec in p.records()] == second
    )  # an increment keeps it
    if admission == LOSSY:
        assert r.call_source_loss == 3


def test_oracle_witnesses():
    # The reference's replay roster is the saved basis, whatever survives, and the offer
    # that captured weights of 0 do not inject is source loss; each against the library, with the
    # old rule as a negative control (it must disagree on the same input).
    import program_oracle as orc

    rows_ = [
        (0, 2, 10),
        (0, 2, 10),
        (1, 5, 10),
        (2, 3, 10),
        (2, 4, 10),
        (3, 5, 10),
        (4, 5, 10),
    ]
    g = graph(6, rows_)
    specs = [(0, [0, 3, 5], 1), (0, [1, 4, 6], 1), (1, [2], 1)]
    excluded = [0, 6]
    alive = [False, False, True]
    r = pinned(
        g,
        tunnels(g, [0, 1], [5], specs),
        LOSSY,
        elig=g.eligibility(excluded_edges=excluded),
    ).replace(2)
    native = ([t.assigned for t in r.by_tunnel], r.source_loss)
    # Head 0's tunnels both fail, but its basis keeps 1, 3, 5: it stays in the equal pool and
    # loses its unit; head 1 carries the other.
    roster = orc.replay_roster(specs, rows_, excluded, [5])
    assert roster == [0, 1]
    assert orc.entry_amounts(2, specs, alive, roster=roster) == native == ([0, 0, 1], 1)
    # The old rule (the heads with a surviving tunnel) gives head 1 both units.
    assert orc.entry_amounts(2, specs, alive, roster=[1]) == ([0, 0, 2], 0) != native
    # Captured weights all 0 (place(1)'s Global program, whose source 1 weighs 0, with source 0
    # masked): the whole offer is lost.
    g = graph(3, [(0, 2, 5), (1, 2, 5)])
    prog = ngc.FlowPolicy(ngc.FlowLedger(g), [0, 1], 2, global_()).place(1).program
    r = pinned(g, prog, LOSSY, elig=g.eligibility(excluded_nodes=[0])).replace(4)
    gspecs = [(t.head, list(t.edges), t.weight) for t in prog.tunnels]
    galive = [t.head != 0 for t in prog.tunnels]
    got = orc.entry_amounts(
        4, gspecs, galive, source_weights=prog.source_weights, roster=[1]
    )
    assert got == ([0] * len(gspecs), 4) and r.source_loss == 4
    assert orc.entry_amounts(5, [(0, [0], 1)], source_weights=[(0, 0)]) == ([0], 5)


def test_global_record_bundle_is_its_share_inside_a_support():
    # A Global record's bundle is the set of edges that carry its share, also inside a supplied
    # support (semantics, "Records of the controllers"); it is the support itself only when the
    # share uses every support edge. Two paths to 3: 0-1-3 at cost 1 + 1 and 0-2-3 at cost 2 + 2,
    # capacity 10 each; the support is all four edges.
    g = graph(4, [(0, 1, 10, 1), (0, 2, 10, 2), (1, 3, 10, 1), (2, 3, 10, 2)])
    support = dag(g, [0, 1, 2, 3])
    p = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 3, global_(support=support))
    # 3 units: the cheaper path carries them all (cost 2 per unit).
    assert p.place(3).placed == 3
    [rec] = p.records()
    assert [(int(e), int(x)) for e, x in rec.edges] == [(0, 3), (2, 3)]
    assert rec.bundle.via_edges().tolist() == [0, 2]
    assert rec.route_cost == 2
    # 10 more: the cheap path fills at 10 and the other carries 3: every support edge, the same
    # value as the support.
    assert p.place(10).placed == 10
    [rec] = p.records()
    assert rec.bundle.via_edges().tolist() == [0, 1, 2, 3]
    assert rec.bundle == support
    assert rec.route_cost == 2


def test_change_flags_of_empty_programs_and_place_zero():
    # Semantics, change flags: the first nonempty program sets both flags; an empty program over
    # none (an empty roster, replace(0) included) sets neither; place(0) prepares and installs
    # nothing, while the full zero offer replace(0) installs the program and follows the rows.
    g = graph(3, [(0, 1, 10), (2, 1, 10)])  # source 0 reaches 1; nothing reaches 2
    empty = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 2, cost_only())
    r = empty.replace(0)
    assert r.program is not None and r.program.tunnels == ()
    assert (r.routes_changed, r.weights_changed) == (False, False)
    p = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 1, cost_only())
    r = p.place(0)
    assert r.program is None and p.program is None
    assert (r.routes_changed, r.weights_changed) == (False, False)
    r = p.replace(0)
    assert r.program is not None and p.program is not None
    assert (r.routes_changed, r.weights_changed) == (True, True)
