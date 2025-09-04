# Re-export types from extension
# pyright: reportMissingImports=false
from _netgraph_core import (
    CostBucket,
    CostDistribution,
    EdgeSelection,
    EdgeTieBreak,
    FlowGraph,
    FlowIndex,
    FlowPlacement,
    FlowPolicy,
    FlowState,
    FlowSummary,
    MinCut,
    PathAlg,
    PredDAG,
    StrictMultiDiGraph,
    batch_max_flow,
    ksp,
    max_flow,
    spf,
)

from ._version import __version__

# Provide richer type information for editors/type-checkers without affecting runtime.
try:  # pragma: no cover - typing-only import
    from typing import TYPE_CHECKING as _TYPE_CHECKING

    if _TYPE_CHECKING:  # noqa: SIM108
        from ._docs import (  # noqa: I001
            CostBucket as CostBucket,
            CostDistribution as CostDistribution,
            EdgeSelection as EdgeSelection,
            EdgeTieBreak as EdgeTieBreak,
            FlowPlacement as FlowPlacement,
            FlowSummary as FlowSummary,
            MinCut as MinCut,
            PredDAG as PredDAG,
        )
except Exception:
    # Safe fallback if _docs.py changes; runtime bindings above remain authoritative.
    pass

__all__ = [
    "__version__",
    "StrictMultiDiGraph",
    "FlowGraph",
    "spf",
    "ksp",
    "max_flow",
    "batch_max_flow",
    "EdgeSelection",
    "EdgeTieBreak",
    "PathAlg",
    "FlowPlacement",
    "PredDAG",
    "FlowIndex",
    "FlowState",
    "FlowPolicy",
    "MinCut",
    "CostBucket",
    "CostDistribution",
    "FlowSummary",
]
