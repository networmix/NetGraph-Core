"""NetGraph-Core: graph, flow and placement models on exact integer units.

Every public name is re-exported from the compiled ``_netgraph_core`` module, so users
import everything from here. Type checkers read the generated stub ``__init__.pyi`` next
to this file (``py.typed`` marks the package as typed).

- The graph on integer units: ``Graph``, ``Eligibility``, ``to_units``, ``as_float``.
- Route structures: ``PredecessorDAG``, ``Path``.
- Route choice and placement settings: ``RouteTieBreak``, ``RouteAllocation``, ``Admission``,
  ``SourceAllocation``, ``WeightBasis``.
- The ledger of placed traffic: ``FlowLedger`` and its ``FlowRecord`` snapshots.
- Forwarding programs: ``ForwardingProgram`` (with its ``Tunnel`` entries,
  ``ForwardingProgram.from_routes`` over ``RouteDefinition`` values and
  ``ForwardingProgram.from_dag``), ``TunnelStatus``.
- The policy engine: ``FlowPolicy``, ``FlowPolicyConfig``, the routings ``ShortestPaths``,
  ``DistributedTE`` (T1, with ``TunnelCount``), ``CentralizedTE`` (with ``tunnel_count``: T2), ``Pinned``, and
  ``PlacementResult`` (with its ``SourceFlow`` and ``TunnelFlow`` rows and its ``TunnelPlan``),
  ``OutputOptions``, ``Termination``.
- Several policies in one operation: ``apply_actions`` over ``Place``, ``Replace`` and ``Forward``
  actions, in a ``BatchMode`` (ORDERED or SIMULTANEOUS).
- The shortest-path cache: ``ShortestPathCache`` (building and frozen modes), ``ShortestPathCacheMiss``.
- The evaluators: ``evaluate_demand``, ``evaluate_pairs``, ``evaluate_edge_failures``.
- The searches: ``spf`` (``Direction``, ``SpfOptions``, ``SpfResult``, ``UNREACHABLE``),
  ``ksp`` (``KspOptions``), ``resolve_to_paths``.

Every public call either completes or raises and leaves the library's state as it was, with two
limits: (1) pybind11 does not check every allocation when it constructs a bound object from
Python, parses the arguments of a call that passes keywords, or converts an attribute or a
returned value into a new Python object outside the library's checked result construction;
running out of memory exactly there may crash the process, with the library's state unchanged.
(2) A ``ShortestPathCache`` passed to the call: a building one keeps the complete entries and
statistics of the lookups that succeeded before the call failed; a frozen one counts one miss
when the call raised ``ShortestPathCacheMiss``. docs/semantics.md states every rule exactly.
"""

from importlib.metadata import version as _version

from _netgraph_core import (
    UNREACHABLE,
    AbsoluteCost,
    Admission,
    BatchMode,
    CentralizedTE,
    CostStretch,
    Direction,
    DistributedTE,
    Eligibility,
    EvaluationOptions,
    FlowLedger,
    FlowPolicy,
    FlowPolicyConfig,
    FlowRecord,
    Forward,
    ForwardingProgram,
    Graph,
    KspOptions,
    OutputOptions,
    Path,
    Pinned,
    Place,
    PlacementResult,
    PredecessorDAG,
    Replace,
    RouteAllocation,
    RouteDefinition,
    RouteTieBreak,
    SensitivityEntry,
    SensitivityResult,
    ShortestPathCache,
    ShortestPathCacheMiss,
    ShortestPathCacheStats,
    ShortestPaths,
    SourceAllocation,
    SourceFlow,
    SpfOptions,
    SpfResult,
    Termination,
    Tunnel,
    TunnelCount,
    TunnelFlow,
    TunnelPlan,
    TunnelStatus,
    WeightBasis,
    apply_actions,
    as_float,
    evaluate_demand,
    evaluate_edge_failures,
    evaluate_pairs,
    ksp,
    resolve_to_paths,
    spf,
    to_units,
)

__version__ = _version("netgraph-core")

__all__ = [
    "__version__",
    "Graph",
    "Eligibility",
    "to_units",
    "as_float",
    "PredecessorDAG",
    "Path",
    "RouteTieBreak",
    "RouteAllocation",
    "Admission",
    "SourceAllocation",
    "WeightBasis",
    "FlowLedger",
    "FlowRecord",
    "TunnelStatus",
    "Tunnel",
    "RouteDefinition",
    "ForwardingProgram",
    "AbsoluteCost",
    "CostStretch",
    "TunnelCount",
    "ShortestPaths",
    "DistributedTE",
    "CentralizedTE",
    "Pinned",
    "FlowPolicyConfig",
    "OutputOptions",
    "Termination",
    "PlacementResult",
    "SourceFlow",
    "TunnelFlow",
    "TunnelPlan",
    "FlowPolicy",
    "Place",
    "Replace",
    "Forward",
    "BatchMode",
    "apply_actions",
    "ShortestPathCache",
    "ShortestPathCacheMiss",
    "ShortestPathCacheStats",
    "EvaluationOptions",
    "SensitivityEntry",
    "SensitivityResult",
    "evaluate_demand",
    "evaluate_pairs",
    "evaluate_edge_failures",
    "Direction",
    "SpfOptions",
    "SpfResult",
    "UNREACHABLE",
    "spf",
    "KspOptions",
    "ksp",
    "resolve_to_paths",
]
