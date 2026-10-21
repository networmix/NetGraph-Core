"""A valid client of the installed package: no type checker diagnostic, and it runs.

dev/stubs/run.sh type-checks it with Pyright against an installed wheel and an installation
built from the sdist (zero diagnostics expected), then runs it. assert_type pins the declared
types of the generated stub: a wrong dtype, a lost Optional or a widened result fails here.
"""

from fractions import Fraction
from typing import SupportsIndex, assert_type

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
assert_type(g.edge_capacity_units(), npt.NDArray[np.int64])
assert_type(g.edge_costs(), npt.NDArray[np.int64])
assert_type(g.edge_sources(), npt.NDArray[np.int32])
assert_type(g.outgoing_offsets(), npt.NDArray[np.int32])
assert_type(g.num_edges(), int)
first_capacity: int = int(g.edge_capacity_units()[0])

eligibility = g.eligibility(excluded_nodes=[np.int64(2)])
assert_type(eligibility.edge_allowed(), npt.NDArray[np.bool_])
both = ng.Eligibility.intersect(eligibility, g.eligibility())
assert_type(both, ng.Eligibility)

# Units: scalars and Python sequences have known shapes; ndarray rank needs narrowing.
assert_type(ng.to_units(1.5, 0.5), int)
assert_type(ng.to_units(np.float32(0.1), np.int64(1)), int)
assert_type(ng.to_units(3, np.float32(0.5)), int)
assert_type(ng.to_units([1.5, 2.0], 0.5), npt.NDArray[np.int64])
assert_type(ng.to_units((np.int64(3), np.float32(2)), 0.5), npt.NDArray[np.int64])
assert_type(ng.to_units(np.array(1.5), 0.5), int | npt.NDArray[np.int64])
assert_type(ng.to_units(np.array([1.5, 2.0]), 0.5), int | npt.NDArray[np.int64])
converted = ng.to_units(np.array([1.5, 2.0]), 0.5)
if isinstance(converted, np.ndarray):
    assert_type(converted, npt.NDArray[np.int64])
    assert converted.tolist() == [3, 4]
assert type(ng.to_units(np.array(1.5), 0.5)) is int
assert_type(ng.as_float(3, 0.5), float)
assert_type(ng.as_float(np.int64(3), 0.5), float)
assert_type(ng.as_float([3, 4], 0.5), npt.NDArray[np.float64])
assert_type(ng.as_float(np.array([3, 4]), 0.5), npt.NDArray[np.float64])


class IndexValue:
    def __index__(self) -> int:
        return 3


def physical_index(x: SupportsIndex) -> float | npt.NDArray[np.float64]:
    assert_type(ng.as_float(x, 0.5), float | npt.NDArray[np.float64])
    return ng.as_float(x, 0.5)


assert physical_index(IndexValue()) == 1.5


# A structural integer protocol may also be a sequence. Sequence classification wins at runtime.
class IndexSequence(list[int]):
    def __index__(self) -> int:
        return 3


assert_type(ng.as_float(IndexSequence([3, 4]), 0.5), npt.NDArray[np.float64])
assert isinstance(physical_index(IndexSequence([3, 4])), np.ndarray)


# Integer subclasses keep the scalar overload even when they expose sequence methods.
class IntegerSequence(int):
    def __len__(self) -> int:
        return 2

    def __getitem__(self, index: int) -> int:
        if index in (0, 1):
            return (3, 4)[index]
        raise IndexError(index)


scalar_from_integer = ng.as_float(IntegerSequence(7), 0.5)
assert_type(scalar_from_integer, float)
assert type(scalar_from_integer) is float
assert scalar_from_integer == 3.5
assert scalar_from_integer.hex() == (3.5).hex()

# The search: optional results are checked before use.
full = ng.spf(
    g, 0, ng.SpfOptions(ng.Direction.FORWARD, tie_break=ng.RouteTieBreak.FEWER_HOPS)
)
costs = full.costs
assert_type(costs, npt.NDArray[np.int64] | None)
if costs is not None:
    reached = int(costs[3])
near = ng.spf(
    g,
    [0],
    ng.SpfOptions(
        ng.Direction.FORWARD, targets=3, multipath=False, eligibility=eligibility
    ),
)
cost = near.target_cost
if cost is not None:
    assert_type(cost, int)
dag = near.dag
assert_type(dag.parents(), npt.NDArray[np.int32])
paths = ng.resolve_to_paths(g, full.dag, 0, 3, max_paths=None)
assert_type(paths, list[ng.Path])
for c, p in ng.ksp(
    g, 0, 3, ng.KspOptions(2, parallel_edges=False, cost_ceiling=ng.AbsoluteCost(9))
):
    assert_type(c, int)
    hops: list[list[int]] = p.hops()
same = ng.PredecessorDAG.from_csr(
    g, dag.parent_offsets(), dag.parents(), np.asarray(dag.via_edges())
)
# An integer array serves as the sequence of hops (one row per hop): the walk 0 -> 1 -> 3.
walk = ng.Path(g, np.array([[0], [2]], dtype=np.int32))
assert_type(walk.cost(), int)

# The ledger and a policy.
ledger = ng.FlowLedger(g, service_baseline=None)
policy = ng.FlowPolicy(
    ledger,
    0,
    [3],
    ng.FlowPolicyConfig(
        routing=ng.DistributedTE(ng.TunnelCount(2), cost_ceiling=ng.CostStretch(1.5)),
        admission=ng.Admission.LOSSLESS,
    ),
)
result = policy.replace(
    10,
    design_demand=8,
    outputs=ng.OutputOptions(residual=True, edge_reservations=True),
)
assert_type(result, ng.PlacementResult)
plan = result.plan
if plan is not None:
    assert_type(plan.design_demand, int | float)
    assert_type(plan.design_by_head, tuple[tuple[int, int], ...])
assert_type(
    result.program.tunnels[0].edges if result.program else None, tuple[int, ...] | None
)
# T2: the controller's tunnels; the plan's exact fields (a Fraction, ints above 2^64).
te16 = ng.FlowPolicy(
    ng.FlowLedger(g),
    0,
    3,
    ng.FlowPolicyConfig(
        routing=ng.CentralizedTE(
            source_allocation=ng.SourceAllocation.EQUAL,
            tunnel_count=ng.TunnelCount(2, up_to=True),
            prefer_headroom=True,
        ),
        admission=ng.Admission.LOSSY,
    ),
)
t2 = te16.replace(4, design_demand=4).plan
if t2 is not None:
    assert_type(t2.max_fractional_bandwidth_per_tunnel, Fraction | None)
    assert_type(t2.heads_with_capacity, tuple[int, ...] | None)
    assert_type(t2.total_tunnel_cost, int | None)
    assert_type(t2.headroom_penalty, int | None)
    assert_type(t2.max_lossless_offer, int | None)
    assert_type(t2.insufficient_resolution, bool | None)
t2 = ng.CentralizedTE(
    source_allocation=ng.SourceAllocation.EQUAL, tunnel_count=ng.TunnelCount(16)
)
assert_type(
    t2.tunnel_count,
    ng.TunnelCount | None,
)
assert_type(t2.prefer_headroom, bool)
assert_type(t2.source_allocation, ng.SourceAllocation)
assert_type(result.program.tunnels[0].planned_bandwidth if result.program else 0, int)
held = result.edge_reservations
if held is not None:
    assert_type(held, npt.NDArray[np.int64])
assert_type(ledger.reserved(), npt.NDArray[np.int64])
assert_type(result.call_undelivered, int | float)
residual = result.residual
if residual is not None:
    assert_type(residual, npt.NDArray[np.int64])
for record in policy.records():
    assert_type(record.edge_flows, list[tuple[int, int]])
    assert_type(record.min_route_cost, int | None)
    assert_type(ledger.record(record.flow_id), ng.FlowRecord)
policy.remove()

# The result's accounting: exact ints, None where a field does not apply, per-source rows.
assert_type(result.requested, int | float)
assert_type(result.offered, int | None)
assert_type(result.unadmitted, int | None)
assert_type(result.total_cost, int | None)
assert_type(result.call_edge_drops, list[tuple[int, int]])
for row in result.by_source:
    assert_type(row, ng.SourceFlow)
    assert_type(row.admitted, int)
    assert_type(row.offered, int | None)

# Forwarding programs: built from routes (a Path or a PredecessorDAG), exported by a policy, forwarded
# by Pinned under a mask and with black holes; the full-offer calls return results too.
route = ng.RouteDefinition(0, ng.Path(g, [[0], [2]]), weight=np.int64(2))
program = ng.ForwardingProgram.from_routes(
    g,
    0,
    [3],
    [ng.RouteDefinition(0, full.dag), route],
    allocation=ng.RouteAllocation.WEIGHTED,
)
assert_type(program, ng.ForwardingProgram)
pooled = ng.ForwardingProgram.from_dag(
    g, [0], np.array([3]), full.dag, allocation=ng.RouteAllocation.EQUAL
)
assert_type(pooled, ng.ForwardingProgram)
assert_type(program.tunnels, tuple[ng.Tunnel, ...])
assert_type(program.source_weights, tuple[tuple[int, int], ...] | None)
for tunnel in program.tunnels:
    assert_type(tunnel.edges, tuple[int, ...] | None)
pinned_cfg = ng.FlowPolicyConfig(
    routing=ng.Pinned(program, black_holes=eligibility), admission=ng.Admission.LOSSY
)
replay = ng.FlowPolicy(ng.FlowLedger(g), 0, 3, pinned_cfg, g.eligibility())
replaced = replay.replace(np.int64(4), outputs=None)
assert_type(replaced, ng.PlacementResult)
for tunnel_row in replaced.by_tunnel:
    assert_type(tunnel_row, ng.TunnelFlow)
    assert_type(tunnel_row.assigned, int)
    assert_type(tunnel_row.status, ng.TunnelStatus)
exported = ng.FlowPolicy(
    ng.FlowLedger(g),
    0,
    3,
    ng.FlowPolicyConfig(
        routing=ng.CentralizedTE(source_allocation=ng.SourceAllocation.OPTIMIZED)
    ),
)
assert_type(exported.evaluate(4).program, ng.ForwardingProgram | None)
exported.replace(4)
assert_type(exported.forward(float("inf"), black_holes=None), ng.PlacementResult)
wcmp = ng.ShortestPaths(
    allocation=ng.RouteAllocation.WEIGHTED,
    weight_basis=ng.WeightBasis.EXPLICIT,
    edge_weights=np.array([1, 2, 1, 1], dtype=np.int64),
)
assert_type(wcmp.edge_weights, tuple[int, ...] | None)
assert_type(ng.FlowPolicyConfig(routing=wcmp).admission, ng.Admission)
lossy_opts = ng.EvaluationOptions(pinned_cfg, offer=4)
assert_type(lossy_opts.offer, int | float | None)
assert_type(
    ng.evaluate_demand(g, 0, 3, lossy_opts).by_tunnel, tuple[ng.TunnelFlow, ...]
)

# The evaluators, with optional arguments left out or None. Only ShortestPaths reads a cache;
# CentralizedTE names its source allocation.
controller = ng.ShortestPaths(
    allocation=ng.RouteAllocation.OPTIMIZED, source_allocation=ng.SourceAllocation.EQUAL
)
assert_type(controller.source_allocation, ng.SourceAllocation | None)
cfg = ng.FlowPolicyConfig(routing=controller)
cache = ng.ShortestPathCache(g)
opts = ng.EvaluationOptions(cfg, eligibility=None, cache=cache)
assert_type(opts.eligibility, ng.Eligibility | None)
call_delivered: int = ng.evaluate_demand(g, 0, 3, opts).call_delivered
batch = ng.evaluate_pairs(g, [(0, 3), ([0, 1], np.array([3]))], opts, workers=1)
assert_type(batch, list[ng.PlacementResult])
# Integer arrays serve as the pairs too: one row (sources, destinations) per pair, or one pair.
from_rows = ng.evaluate_pairs(g, np.array([[0, 3], [1, 3]], dtype=np.int64), opts)
from_pair = ng.evaluate_pairs(g, [np.array([0, 3], dtype=np.int64)], opts)
assert_type(from_rows, list[ng.PlacementResult])
assert_type(from_pair, list[ng.PlacementResult])
global_cfg = ng.FlowPolicyConfig(
    routing=ng.CentralizedTE(source_allocation=ng.SourceAllocation.OPTIMIZED)
)
sens = ng.evaluate_edge_failures(g, 0, 3, ng.EvaluationOptions(global_cfg))
for entry in sens.entries:
    assert_type(entry.delivery_reduction, int)
assert_type(cache.stats().hits, int)
routing = global_cfg.routing
if isinstance(routing, ng.CentralizedTE):
    print("CentralizedTE routing:", routing)
print(call_delivered, len(batch), len(paths), ng.UNREACHABLE, ng.__version__)

# CostStretch takes a Fraction, a float, a NumPy floating value or an integer; factor is a Fraction.
stretches = [
    ng.CostStretch(Fraction(6, 5)),
    ng.CostStretch(1.2),
    ng.CostStretch(np.float32(1.2)),
    ng.CostStretch(2),
]
assert_type(stretches[0].factor, Fraction)

# Several policies in one operation: the actions keep their offers as given, the list returns one
# result per action.
listed = ng.FlowLedger(g)
lossy_cfg = ng.FlowPolicyConfig(
    routing=ng.ShortestPaths(), admission=ng.Admission.LOSSY
)
first_policy = ng.FlowPolicy(listed, 0, 3, lossy_cfg)
second_policy = ng.FlowPolicy(listed, 1, 3, lossy_cfg)
actions = [ng.Replace(first_policy, 2), ng.Replace(second_policy, np.int64(2))]
assert_type(actions[0].offer, int | float)
assert_type(actions[0].design_demand, int | float | None)
assert_type(ng.Forward(first_policy, 1).black_holes, ng.Eligibility | None)
together = ng.apply_actions(listed, actions, mode=ng.BatchMode.SIMULTANEOUS)
assert_type(together, list[ng.PlacementResult])
in_order = ng.apply_actions(listed, (ng.Place(first_policy, 1),))
assert_type(in_order[0].call_delivered, int)

# Units' floating branches describe allowed positive infinity, including NumPy float32.
units_inf = np.float32(np.inf)
assert ng.EvaluationOptions(
    ng.FlowPolicyConfig(routing=ng.ShortestPaths()), offer=units_inf
).offer == float("inf")
