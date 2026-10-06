"""Runtime package surface for :mod:`netgraph_core`.

All public symbols are re-exported from the compiled ``_netgraph_core`` module so
downstream users import everything from a single place. Typing shims supplied in
``_docs.py`` keep editors happy without affecting runtime behavior.

The package carries the graph on integer units (``Graph``, ``Eligibility``
with ``intersect``, ``to_units``, ``as_float``), the route structures (``PredDAG``,
``Path``), the selection and placement axes (``PathSelection``, ``EdgeTieBreak``,
``EdgePriority``, ``FlowPlacement``, ``Admission``), the integer ledger (``FlowLedger``
with its ``FlowRecord`` snapshots), the policy engine (``FlowPolicy`` with
``FlowPolicyConfig`` over the routing variants ``CostOnly``, ``CapacityAware``,
``Global``, ``Pinned``; ``PlacementResult``, ``OutputOptions``, ``Termination``), the
``DagCache`` (building mode) and the evaluators (``capacity``, ``batch_capacity``,
``sensitivity``). The public search API (shortest paths, k shortest paths, path
resolution) and the frozen cache are not part of this build; ``_netgraph_core`` also holds the private entries
of the flow core and the internal search, which are never re-exported here.
"""

# pyright: reportMissingImports=false

from importlib.metadata import version
from typing import TYPE_CHECKING

from _netgraph_core import (
    MAX_FLOWS_PER_CALL,
    Absolute,
    Admission,
    CapacityAware,
    CapacityOptions,
    CostOnly,
    DagCache,
    DagCacheStats,
    DiminishingReturns,
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
    OutputOptions,
    Path,
    PathSelection,
    Pinned,
    PlacementResult,
    PredDAG,
    SensitivityEntry,
    SensitivityOptions,
    SensitivityResult,
    SensitivityStatus,
    Stretch,
    Termination,
    as_float,
    batch_capacity,
    capacity,
    sensitivity,
    to_units,
)

__version__ = version("netgraph-core")

# Provide richer type information for editors/type-checkers without affecting runtime.
# Every name here must also appear in the runtime import above; the aliases only
# change what a type checker sees. Keep this list in sync with _docs.py, and keep
# _docs.py in sync with the bindings -- a wrong stub is worse than a missing one,
# because it turns "unchecked" into "confidently wrong".
if TYPE_CHECKING:  # pragma: no cover - typing-only
    from ._docs import (  # noqa: I001
        MAX_FLOWS_PER_CALL as MAX_FLOWS_PER_CALL,
        Absolute as Absolute,
        Admission as Admission,
        CapacityAware as CapacityAware,
        CapacityOptions as CapacityOptions,
        CostOnly as CostOnly,
        DagCache as DagCache,
        DagCacheStats as DagCacheStats,
        DiminishingReturns as DiminishingReturns,
        EdgePriority as EdgePriority,
        EdgeTieBreak as EdgeTieBreak,
        Eligibility as Eligibility,
        FlowLedger as FlowLedger,
        FlowPlacement as FlowPlacement,
        FlowPolicy as FlowPolicy,
        FlowPolicyConfig as FlowPolicyConfig,
        FlowRecord as FlowRecord,
        Global as Global,
        Graph as Graph,
        OutputOptions as OutputOptions,
        Path as Path,
        PathSelection as PathSelection,
        Pinned as Pinned,
        PlacementResult as PlacementResult,
        PredDAG as PredDAG,
        SensitivityEntry as SensitivityEntry,
        SensitivityOptions as SensitivityOptions,
        SensitivityResult as SensitivityResult,
        SensitivityStatus as SensitivityStatus,
        Stretch as Stretch,
        Termination as Termination,
        as_float as as_float,
        batch_capacity as batch_capacity,
        capacity as capacity,
        sensitivity as sensitivity,
        to_units as to_units,
    )

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
    "FlowPolicy",
    "MAX_FLOWS_PER_CALL",
    "DagCache",
    "DagCacheStats",
    "CapacityOptions",
    "SensitivityOptions",
    "SensitivityStatus",
    "SensitivityEntry",
    "SensitivityResult",
    "capacity",
    "batch_capacity",
    "sensitivity",
]
