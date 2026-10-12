"""A valid client of the installed package: no type checker diagnostic, and it runs.

dev/stubs/run.sh type-checks it with Pyright against an installed wheel and an installation
built from the sdist (zero diagnostics expected), then runs it. assert_type pins the declared
types of the generated stub: a wrong dtype, a lost Optional or a widened result fails here.
"""

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
    pinned = ng.Pinned([p.to_dag()])
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
        routing=ng.CapacityAware(flow_count=2, max_path_cost=ng.Stretch(1.5)),
        placement=ng.FlowPlacement.EQUAL_BALANCED,
        admission=ng.Admission.LOSSLESS,
        selection=sel,
    ),
)
result = policy.place(float("inf"), ng.OutputOptions(residual=True))
assert_type(result, ng.PlacementResult)
assert_type(result.remaining, int | float)
residual = result.residual
if residual is not None:
    assert_type(residual, npt.NDArray[np.int64])
for record in policy.records():
    assert_type(record.edges, list[tuple[int, int]])
    assert_type(record.route_cost, int | None)
    assert_type(ledger.record(record.flow_id), ng.FlowRecord)
policy.remove()

# The evaluators, with optional arguments left out or None.
cfg = ng.FlowPolicyConfig(
    routing=ng.Global(),
    placement=ng.FlowPlacement.MAX_FLOW,
    admission=ng.Admission.LOSSLESS,
    selection=sel,
)
cache = ng.DagCache(g)
opts = ng.CapacityOptions(cfg, elig=None, cache=cache)
assert_type(opts.elig, ng.Eligibility | None)
placed: int = ng.capacity(g, 0, 3, opts).placed
batch = ng.batch_capacity(g, [(0, 3), ([0, 1], np.array([3]))], opts, workers=1)
assert_type(batch, list[ng.PlacementResult])
# Integer arrays serve as the pairs too: one row (sources, destinations) per pair, or one pair.
from_rows = ng.batch_capacity(g, np.array([[0, 3], [1, 3]], dtype=np.int64), opts)
from_pair = ng.batch_capacity(g, [np.array([0, 3], dtype=np.int64)], opts)
assert_type(from_rows, list[ng.PlacementResult])
assert_type(from_pair, list[ng.PlacementResult])
sens = ng.sensitivity(g, 0, 3, ng.SensitivityOptions(cfg))
for entry in sens.entries:
    assert_type(entry.loss, int)
assert_type(cache.stats().hits, int)
routing = cfg.routing
if isinstance(routing, ng.Global):
    print("Global routing:", routing)
print(placed, len(batch), len(paths), ng.UNREACHABLE, ng.__version__)
