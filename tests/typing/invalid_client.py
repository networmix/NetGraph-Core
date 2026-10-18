"""A deliberately invalid client: every marked line must get exactly its expected diagnostic.

dev/stubs/run.sh type-checks it with Pyright against an installed wheel and an installation
built from the sdist. Each line ending in "# expect: RULE[, RULE...]" must produce one
diagnostic of each listed Pyright rule, and no other line may produce any. It is never run:
each marked line is a mistake the compiled module rejects at run time (or, for the Optional
and dtype lines, a mistake in the caller's own code) that a type checker must catch before.
"""

import numpy as np
import numpy.typing as npt

import netgraph_core as ng

g = ng.Graph(2, [0], [1], [10], [1], 1.0)
sel = ng.PathSelection(ng.EdgeTieBreak.DETERMINISTIC, ng.EdgePriority.CAPACITY)
cfg = ng.FlowPolicyConfig(routing=ng.Global())

# The twelve result and factory classes have no constructor ("No constructor defined!").
ng.Eligibility()  # expect: reportCallIssue
ng.PredDAG()  # expect: reportCallIssue
ng.FlowRecord()  # expect: reportCallIssue
ng.PlacementResult()  # expect: reportCallIssue
ng.SourceFlow()  # expect: reportCallIssue
ng.DagCacheStats()  # expect: reportCallIssue
ng.SensitivityEntry()  # expect: reportCallIssue
ng.SensitivityResult()  # expect: reportCallIssue
ng.SpfResult()  # expect: reportCallIssue
ng.ForwardingProgram()  # expect: reportCallIssue
ng.Tunnel()  # expect: reportCallIssue
ng.TunnelFlow()  # expect: reportCallIssue

# dtype mistakes: the capacities are an int64 array of ints. (NumPy's own stubs type an element
# read arr[i] as Any, so "text: str = g.capacity_units()[0]" cannot be caught by any stub; an
# element read through .item() and the array itself carry the dtype.)
floats: npt.NDArray[np.float64] = g.capacity_units()  # expect: reportAssignmentType
text: str = g.capacity_units().item(0)  # expect: reportAssignmentType

# Wrong argument types (the module raises ValueError or TypeError).
ng.Graph("wrong", [0], [1], [10], [1], 1.0)  # expect: reportArgumentType
ng.Graph(2, [0], [1], [1.5], [1], 1.0)  # expect: reportArgumentType
ng.FlowLedger(g, residual=np.zeros(1, dtype=np.float64))  # expect: reportArgumentType
ng.spf(g, 0, ng.KspOptions(1, sel))  # expect: reportArgumentType
ng.capacity(g, 0, 1, cfg)  # expect: reportArgumentType
guard = ng.DiminishingReturns(2, 0.5)
ng.CapacityAware(max_path_cost=guard)  # expect: reportArgumentType
ng.to_units("1.5", 1.0)  # expect: reportCallIssue, reportArgumentType
# FlowPolicyConfig takes keywords only.
ng.FlowPolicyConfig(ng.Global())  # expect: reportCallIssue
dag = ng.Path(g, [[0]]).to_dag()
ng.Pinned([dag])  # expect: reportArgumentType
# Integer arrays serve as hops and pairs; float arrays do not.
ng.Path(g, np.zeros((1, 1), dtype=np.float64))  # expect: reportArgumentType
float_pairs = np.zeros((1, 2), dtype=np.float64)
ng.batch_capacity(g, float_pairs, ng.QueryOptions(cfg))  # expect: reportArgumentType

# Missing Optional checks: these results may be None.
result = ng.spf(g, 0, ng.SpfOptions(ng.Direction.FORWARD, sel, targets=1))
cost: int = result.target_cost  # expect: reportAssignmentType
first = result.costs[0]  # expect: reportOptionalSubscript
flows = ng.capacity(g, 0, 1, ng.QueryOptions(cfg)).edge_flows
total = flows.sum()  # expect: reportOptionalMemberAccess
best = ng.capacity(g, 0, 1, ng.QueryOptions(cfg))
cost_total: int = best.total_cost  # expect: reportAssignmentType
weights = best.program.source_weights  # expect: reportOptionalMemberAccess
