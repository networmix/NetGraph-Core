"""A deliberately invalid client: every marked line must get exactly its expected diagnostic.

dev/stubs/run.sh type-checks it with Pyright against an installed wheel and an installation
built from the sdist. Each line ending in "# expect: RULE[, RULE...]" must produce one
diagnostic of each listed Pyright rule, and no other line may produce any. It is never run:
each marked line is a mistake the compiled module rejects at run time (or, for the Optional
and dtype lines, a mistake in the caller's own code) that a type checker must catch before.
"""

import typing
from decimal import Decimal
from fractions import Fraction

import numpy as np
import numpy.typing as npt

import netgraph_core as ng

g = ng.Graph(2, [0], [1], [10], [1], 1.0)
opt = ng.SourceAllocation.OPTIMIZED
cfg = ng.FlowPolicyConfig(routing=ng.CentralizedTE(source_allocation=opt))

# The twelve result and factory classes have no constructor ("No constructor defined!").
ng.Eligibility()  # expect: reportCallIssue
ng.PredecessorDAG()  # expect: reportCallIssue
ng.FlowRecord()  # expect: reportCallIssue
ng.PlacementResult()  # expect: reportCallIssue
ng.SourceFlow()  # expect: reportCallIssue
ng.ShortestPathCacheStats()  # expect: reportCallIssue
ng.SensitivityEntry()  # expect: reportCallIssue
ng.SensitivityResult()  # expect: reportCallIssue
ng.SpfResult()  # expect: reportCallIssue
ng.ForwardingProgram()  # expect: reportCallIssue
ng.Tunnel()  # expect: reportCallIssue
ng.TunnelFlow()  # expect: reportCallIssue

# dtype mistakes: the capacities are an int64 array of ints. (NumPy's own stubs type an element
# read arr[i] as Any, so "text: str = g.edge_capacity_units()[0]" cannot be caught by any stub; an
# element read through .item() and the array itself carry the dtype.)
# fmt: off
floats: npt.NDArray[np.float64] = g.edge_capacity_units()  # expect: reportAssignmentType
# fmt: on
text: str = g.edge_capacity_units().item(0)  # expect: reportAssignmentType

# Wrong argument types (the module raises ValueError or TypeError).
ng.Graph("wrong", [0], [1], [10], [1], 1.0)  # expect: reportArgumentType
ng.Graph(2, [0], [1], [1.5], [1], 1.0)  # expect: reportArgumentType
# fmt: off
ng.FlowLedger(g, service_baseline=np.zeros(1, dtype=np.float64))  # expect: reportArgumentType
# fmt: on
ng.spf(g, 0, ng.KspOptions(1))  # expect: reportArgumentType
ng.SpfOptions(ng.Direction.FORWARD, multipath=None)  # expect: reportArgumentType
ng.CostStretch("1.2")  # expect: reportArgumentType
ng.evaluate_demand(g, 0, 1, cfg)  # expect: reportArgumentType
ng.EvaluationOptions(cfg, offer="4")  # expect: reportArgumentType
guard = ng.TunnelCount(2)
ng.DistributedTE(guard, cost_ceiling=guard)  # expect: reportArgumentType
ng.CentralizedTE(source_allocation=opt, tunnel_count=2)  # expect: reportArgumentType
# fmt: off
ng.CentralizedTE(source_allocation=opt, prefer_headroom="yes")  # expect: reportArgumentType
# fmt: on
# CentralizedTE names its source allocation; the allocations and weights are the models' fields.
ng.CentralizedTE()  # expect: reportCallIssue
sp = ng.ShortestPaths()
ng.FlowPolicyConfig(routing=sp, placement=None)  # expect: reportCallIssue
ng.ShortestPaths(allocation="EQUAL")  # expect: reportArgumentType
ng.to_units("1.5", 1.0)  # expect: reportCallIssue, reportArgumentType
# FlowPolicyConfig takes keywords only.
ng.FlowPolicyConfig(ng.CentralizedTE(source_allocation=opt))  # expect: reportCallIssue
dag = ng.Path(g, [[0]]).to_dag()
ng.Pinned([dag])  # expect: reportArgumentType
# Integer arrays serve as hops and pairs; float arrays do not.
ng.Path(g, np.zeros((1, 1), dtype=np.float64))  # expect: reportArgumentType
float_pairs = np.zeros((1, 2), dtype=np.float64)
# fmt: off
ng.evaluate_pairs(g, float_pairs, ng.EvaluationOptions(cfg))  # expect: reportArgumentType
# fmt: on

# Missing Optional checks: these results may be None.
result = ng.spf(g, 0, ng.SpfOptions(ng.Direction.FORWARD, targets=1))
cost: int = result.target_cost  # expect: reportAssignmentType
first = result.costs[0]  # expect: reportOptionalSubscript
flows = ng.evaluate_demand(g, 0, 1, ng.EvaluationOptions(cfg)).edge_flows
total = flows.sum()  # expect: reportOptionalMemberAccess
best = ng.evaluate_demand(g, 0, 1, ng.EvaluationOptions(cfg))
cost_total: int = best.total_cost  # expect: reportAssignmentType
weights = best.program.source_weights  # expect: reportOptionalMemberAccess

# A list takes actions, its mode by keyword only, and an action takes a number.
lister = ng.FlowLedger(g)
listed_policy = ng.FlowPolicy(lister, 0, 1, cfg)
ng.apply_actions(lister, [listed_policy])  # expect: reportArgumentType
ng.apply_actions(lister, [], ng.BatchMode.ORDERED)  # expect: reportCallIssue
ng.Replace(listed_policy, "1")  # expect: reportArgumentType


class PolicySubclass(ng.FlowPolicy):  # expect: reportGeneralTypeIssues
    pass


# Unknown ndarray rank cannot promise an array, and a general index protocol may be a sequence.
a0 = np.array(1.5)
c: npt.NDArray[np.int64] = ng.to_units(a0, 0.5)  # expect: reportAssignmentType


def require_float(x: typing.SupportsIndex) -> float:
    return ng.as_float(x, 0.5)  # expect: reportReturnType


ng.to_units(1.5, Fraction(1, 2))  # expect: reportCallIssue, reportArgumentType
ng.as_float(3, Fraction(1, 2))  # expect: reportCallIssue, reportArgumentType


class IndexResolution:
    def __index__(self) -> int:
        return 1


class FloatResolution:
    def __float__(self) -> float:
        return 1.0


ng.Graph(2, [0], [1], [10], [1], Fraction(1))  # expect: reportArgumentType
ng.Graph(2, [0], [1], [10], [1], Decimal(1))  # expect: reportArgumentType
ng.Graph(2, [0], [1], [10], [1], IndexResolution())  # expect: reportArgumentType
ng.Graph(2, [0], [1], [10], [1], FloatResolution())  # expect: reportArgumentType


class BothResolution(IndexResolution, FloatResolution):
    pass


ng.Graph(2, [0], [1], [10], [1], BothResolution())  # expect: reportArgumentType
ng.Graph(2, [0], [1], [10], [1], np.array(1.0))  # expect: reportArgumentType

ng.CentralizedTE(
    source_allocation=ng.SourceAllocation.WEIGHTED,
    source_weights=[1.5],  # expect: reportArgumentType
)
