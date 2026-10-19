"""A valid client of the installed package: no type checker diagnostic, and it runs.

dev/stubs/run.sh type-checks it with Pyright against an installed wheel and an installation
built from the sdist (zero diagnostics expected), then runs it. assert_type pins the declared
types of the generated stub: a wrong dtype, a lost Optional or a widened result fails here.
"""

from fractions import Fraction
from typing import assert_type

import numpy as np
import numpy.typing as npt

import netgraph_core as ng

# The graph: integers as lists, tuples or NumPy arrays of any integer dtype.
g = ng.Graph(
    4,
    np.array([0, 0, 1, 2], dtype=np.int32),
    (1, 2, 3, 3),
    [10, 5, ng.to_units(2.5, 1.0), 4],
    np.array([1, 2, 1, 1], dtype=np.uint8),
    1.0,
)
assert_type(g.capacity_units(), npt.NDArray[np.int64])
assert_type(g.cost(), npt.NDArray[np.int64])
assert_type(g.edge_src(), npt.NDArray[np.int32])
assert_type(g.row_offsets(), npt.NDArray[np.int32])
assert_type(g.num_edges(), int)
first_capacity: int = int(g.capacity_units()[0])

elig = g.eligibility(excluded_nodes=[np.int64(2)])
assert_type(elig.edge_allowed(), npt.NDArray[np.bool_])
both = ng.Eligibility.intersect(elig, g.eligibility())
assert_type(both, ng.Eligibility)

# Units: a scalar gives an int, an array an int64 array; and back.
assert_type(ng.to_units(1.5, 0.5), int)
assert_type(ng.to_units([1.5, 2.0], 0.5), npt.NDArray[np.int64])
assert_type(ng.as_float(3, 0.5), float)
assert_type(ng.as_float(np.array([3, 4]), 0.5), npt.NDArray[np.float64])

# The search: optional results are checked before use.
sel = ng.PathSelection(ng.EdgeTieBreak.DETERMINISTIC, ng.EdgePriority.CAPACITY)
full = ng.spf(g, 0, ng.SpfOptions(ng.Direction.FORWARD, sel))
costs = full.costs
assert_type(costs, npt.NDArray[np.int64] | None)
if costs is not None:
    reached = int(costs[3])
near = ng.spf(g, [0], ng.SpfOptions(ng.Direction.FORWARD, sel, targets=3, elig=elig))
cost = near.target_cost
if cost is not None:
    assert_type(cost, int)
dag = near.dag
assert_type(dag.parents(), npt.NDArray[np.int32])
paths = ng.resolve_to_paths(g, full.dag, 0, 3, max_paths=None)
assert_type(paths, list[ng.Path])
for c, p in ng.ksp(g, 0, 3, ng.KspOptions(2, sel, max_path_cost=ng.Absolute(9))):
    assert_type(c, int)
    hops: list[list[int]] = p.hops()
same = ng.PredDAG.from_csr(
    g, dag.parent_offsets(), dag.parents(), np.asarray(dag.via_edges())
)
# An integer array serves as the sequence of hops (one row per hop): the walk 0 -> 1 -> 3.
walk = ng.Path(g, np.array([[0], [2]], dtype=np.int32))
assert_type(walk.cost(), int)

# The ledger and a policy.
ledger = ng.FlowLedger(g, residual=None)
policy = ng.FlowPolicy(
    ledger,
    0,
    [3],
    ng.FlowPolicyConfig(
        routing=ng.CapacityAware(ng.Tunnels(2), max_path_cost=ng.Stretch(1.5)),
        admission=ng.Admission.LOSSLESS,
    ),
)
result = policy.replace(
    10,
    design_demand=8,
    outputs=ng.OutputOptions(residual=True, reservation_edge_flows=True),
)
assert_type(result, ng.PlacementResult)
plan = result.plan
if plan is not None:
    assert_type(plan.design_demand, int | float)
    assert_type(plan.head_design, tuple[tuple[int, int], ...])
assert_type(
    result.program.tunnels[0].edges if result.program else None, tuple[int, ...] | None
)
# T2: the controller's tunnels; the plan's exact fields (a Fraction, ints above 2^64).
te16 = ng.FlowPolicy(
    ng.FlowLedger(g),
    0,
    3,
    ng.FlowPolicyConfig(
        routing=ng.Global(tunnels=ng.Tunnels(2, up_to=True), prefer_headroom=True),
        admission=ng.Admission.LOSSY,
    ),
)
t2 = te16.replace(4, design_demand=4).plan
if t2 is not None:
    assert_type(t2.fluid_share, Fraction | None)
    assert_type(t2.live_heads, tuple[int, ...] | None)
    assert_type(t2.route_cost, int | None)
    assert_type(t2.headroom_penalty, int | None)
    assert_type(t2.forwarding_bound, int | None)
    assert_type(t2.insufficient_resolution, bool | None)
assert_type(ng.Global(tunnels=ng.Tunnels(16)).tunnels, ng.Tunnels | None)
assert_type(ng.Global(tunnels=ng.Tunnels(16)).prefer_headroom, bool)
assert_type(result.program.tunnels[0].planned if result.program else 0, int)
held = result.reservation_edge_flows
if held is not None:
    assert_type(held, npt.NDArray[np.int64])
assert_type(ledger.reserved(), npt.NDArray[np.int64])
assert_type(result.remaining, int | float)
residual = result.residual
if residual is not None:
    assert_type(residual, npt.NDArray[np.int64])
for record in policy.records():
    assert_type(record.edges, list[tuple[int, int]])
    assert_type(record.route_cost, int | None)
    assert_type(ledger.record(record.flow_id), ng.FlowRecord)
policy.remove()

# The result's accounting: exact ints, None where a field does not apply, per-source rows.
assert_type(result.requested, int | float)
assert_type(result.offered, int | None)
assert_type(result.unadmitted, int | None)
assert_type(result.total_cost, int | None)
assert_type(result.call_drops, list[tuple[int, int]])
for row in result.by_source:
    assert_type(row, ng.SourceFlow)
    assert_type(row.admitted, int)
    assert_type(row.offered, int | None)

# Forwarding programs: built from routes (a Path or a PredDAG), exported by a policy, forwarded
# by Pinned under a mask and with black holes; the full-offer calls return results too.
route = ng.Route(0, ng.Path(g, [[0], [2]]), weight=np.int64(2))
program = ng.ForwardingProgram.from_routes(
    g, 0, [3], [ng.Route(0, full.dag), route], placement=ng.FlowPlacement.WEIGHTED
)
assert_type(program, ng.ForwardingProgram)
assert_type(program.tunnels, tuple[ng.Tunnel, ...])
assert_type(program.source_weights, tuple[tuple[int, int], ...] | None)
for tunnel in program.tunnels:
    assert_type(tunnel.edges, tuple[int, ...] | None)
pinned_cfg = ng.FlowPolicyConfig(
    routing=ng.Pinned(program, black_holes=elig), admission=ng.Admission.LOSSY
)
replay = ng.FlowPolicy(ng.FlowLedger(g), 0, 3, pinned_cfg, g.eligibility())
replaced = replay.replace(np.int64(4), outputs=None)
assert_type(replaced, ng.PlacementResult)
for tunnel_row in replaced.by_tunnel:
    assert_type(tunnel_row, ng.TunnelFlow)
    assert_type(tunnel_row.assigned, int)
    assert_type(tunnel_row.status, ng.TunnelStatus)
exported = ng.FlowPolicy(
    ng.FlowLedger(g), 0, 3, ng.FlowPolicyConfig(routing=ng.Global())
)
assert_type(exported.evaluate(4).program, ng.ForwardingProgram | None)
exported.replace(4)
assert_type(exported.forward(float("inf"), black_holes=None), ng.PlacementResult)
weighted = ng.FlowPolicyConfig(
    routing=ng.CostOnly(),
    placement=ng.FlowPlacement.WEIGHTED,
    weights=ng.WeightBasis.EXPLICIT,
    edge_weights=np.array([1, 2, 1, 1], dtype=np.int64),
)
assert_type(weighted.edge_weights, tuple[int, ...] | None)
lossy_opts = ng.QueryOptions(pinned_cfg, offered=4)
assert_type(lossy_opts.offered, int | float | None)
assert_type(ng.capacity(g, 0, 3, lossy_opts).by_tunnel, tuple[ng.TunnelFlow, ...])

# The evaluators, with optional arguments left out or None. Only CostOnly reads a cache; the
# controllers' origin defaults to OPTIMIZED.
cfg = ng.FlowPolicyConfig(
    routing=ng.CostOnly(),
    placement=ng.FlowPlacement.MAX_FLOW,
    origin=ng.Origin.EQUAL,
)
assert_type(cfg.origin, ng.Origin | None)
cache = ng.DagCache(g)
opts = ng.QueryOptions(cfg, elig=None, cache=cache)
assert_type(opts.elig, ng.Eligibility | None)
placed: int = ng.capacity(g, 0, 3, opts).placed
batch = ng.batch_capacity(g, [(0, 3), ([0, 1], np.array([3]))], opts, workers=1)
assert_type(batch, list[ng.PlacementResult])
# Integer arrays serve as the pairs too: one row (sources, destinations) per pair, or one pair.
from_rows = ng.batch_capacity(g, np.array([[0, 3], [1, 3]], dtype=np.int64), opts)
from_pair = ng.batch_capacity(g, [np.array([0, 3], dtype=np.int64)], opts)
assert_type(from_rows, list[ng.PlacementResult])
assert_type(from_pair, list[ng.PlacementResult])
global_cfg = ng.FlowPolicyConfig(routing=ng.Global())
sens = ng.sensitivity(g, 0, 3, ng.QueryOptions(global_cfg))
for entry in sens.entries:
    assert_type(entry.loss, int)
assert_type(cache.stats().hits, int)
routing = global_cfg.routing
if isinstance(routing, ng.Global):
    print("Global routing:", routing)
print(placed, len(batch), len(paths), ng.UNREACHABLE, ng.__version__)
