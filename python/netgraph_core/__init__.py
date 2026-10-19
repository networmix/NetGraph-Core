"""NetGraph-Core: graph, flow and placement models on exact integer units.

Every public name is re-exported from the compiled ``_netgraph_core`` module, so users
import everything from here. Type checkers read the generated stub ``__init__.pyi`` next
to this file (``py.typed`` marks the package as typed).

- The graph on integer units: ``Graph``, ``Eligibility``, ``to_units``, ``as_float``.
- Route structures: ``PredecessorDAG``, ``Path``.
- Route choice and placement settings: ``RouteTieBreak``, ``RouteAllocation``, ``Admission``,
  ``SourceAllocation``, ``WeightBasis``.
- The ledger of placed traffic: ``FlowLedger`` and its ``FlowRecord`` snapshots.
- Forwarding programs: ``ForwardingProgram`` (with its ``Tunnel`` entries and
  ``ForwardingProgram.from_routes`` over ``RouteDefinition`` values), ``TunnelStatus``.
- The policy engine: ``FlowPolicy``, ``FlowPolicyConfig``, the routings ``ShortestPaths``,
  ``DistributedTE`` (T1, with ``TunnelCount``), ``CentralizedTE`` (with ``tunnel_count``: T2), ``Pinned``, and
  ``PlacementResult`` (with its ``SourceFlow`` and ``TunnelFlow`` rows and its ``TunnelPlan``),
  ``OutputOptions``, ``Termination``.
- The shortest-path cache: ``ShortestPathCache`` (building and frozen modes), ``ShortestPathCacheMiss``.
- The evaluators: ``evaluate_demand``, ``evaluate_pairs``, ``evaluate_edge_failures``.
- The searches: ``spf`` (``Direction``, ``SpfOptions``, ``SpfResult``, ``UNREACHABLE``),
  ``ksp`` (``KspOptions``), ``resolve_to_paths``.

``_netgraph_core`` also holds a private entry of the flow core, never re-exported here.
"""

from importlib.metadata import version as _version

from _netgraph_core import (
    UNREACHABLE,
    AbsoluteCost,
    Admission,
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
    ForwardingProgram,
    Graph,
    KspOptions,
    OutputOptions,
    Path,
    Pinned,
    PlacementResult,
    PredecessorDAG,
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
