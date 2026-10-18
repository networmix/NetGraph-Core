"""NetGraph-Core: graph, flow and placement models on exact integer units.

Every public name is re-exported from the compiled ``_netgraph_core`` module, so users
import everything from here. Type checkers read the generated stub ``__init__.pyi`` next
to this file (``py.typed`` marks the package as typed).

- The graph on integer units: ``Graph``, ``Eligibility``, ``to_units``, ``as_float``.
- Route structures: ``PredDAG``, ``Path``.
- Route selection and placement settings: ``PathSelection``, ``EdgeTieBreak``,
  ``EdgePriority``, ``FlowPlacement``, ``Admission``, ``Origin``.
- The ledger of placed traffic: ``FlowLedger`` and its ``FlowRecord`` snapshots.
- The policy engine: ``FlowPolicy``, ``FlowPolicyConfig``, the routings ``CostOnly``,
  ``CapacityAware``, ``Global``, ``Pinned``, and ``PlacementResult`` (with its ``SourceFlow``
  rows), ``OutputOptions``, ``Termination``.
- The shortest-path cache: ``DagCache`` (building and frozen modes), ``DagCacheMiss``.
- The evaluators: ``capacity``, ``batch_capacity``, ``sensitivity``.
- The searches: ``spf`` (``Direction``, ``SpfOptions``, ``SpfResult``, ``UNREACHABLE``),
  ``ksp`` (``KspOptions``), ``resolve_to_paths``.

``_netgraph_core`` also holds a private entry of the flow core, never re-exported here.
"""

from importlib.metadata import version as _version

from _netgraph_core import (
    MAX_FLOWS_PER_CALL,
    UNREACHABLE,
    Absolute,
    Admission,
    CapacityAware,
    CostOnly,
    DagCache,
    DagCacheMiss,
    DagCacheStats,
    DiminishingReturns,
    Direction,
    EdgePriority,
    EdgeTieBreak,
    Eligibility,
    FlowLedger,
    FlowPlacement,
    FlowPolicy,
    FlowPolicyConfig,
    FlowRecord,
    Global,
    Graph,
    KspOptions,
    Origin,
    OutputOptions,
    Path,
    PathSelection,
    Pinned,
    PlacementResult,
    PredDAG,
    QueryOptions,
    SensitivityEntry,
    SensitivityResult,
    SourceFlow,
    SpfOptions,
    SpfResult,
    Stretch,
    Termination,
    as_float,
    batch_capacity,
    capacity,
    ksp,
    resolve_to_paths,
    sensitivity,
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
    "PredDAG",
    "Path",
    "EdgeTieBreak",
    "EdgePriority",
    "PathSelection",
    "FlowPlacement",
    "Admission",
    "Origin",
    "FlowLedger",
    "FlowRecord",
    "Absolute",
    "Stretch",
    "DiminishingReturns",
    "CostOnly",
    "CapacityAware",
    "Global",
    "Pinned",
    "FlowPolicyConfig",
    "OutputOptions",
    "Termination",
    "PlacementResult",
    "SourceFlow",
    "FlowPolicy",
    "MAX_FLOWS_PER_CALL",
    "DagCache",
    "DagCacheMiss",
    "DagCacheStats",
    "QueryOptions",
    "SensitivityEntry",
    "SensitivityResult",
    "capacity",
    "batch_capacity",
    "sensitivity",
    "Direction",
    "SpfOptions",
    "SpfResult",
    "UNREACHABLE",
    "spf",
    "KspOptions",
    "ksp",
    "resolve_to_paths",
]
