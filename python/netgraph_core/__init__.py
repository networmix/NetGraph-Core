# Re-export types from extension
# pyright: reportMissingImports=false
from _netgraph_core import (
    Algorithms,
    Backend,
    EdgeSelection,
    EdgeTieBreak,
    FlowGraph,
    FlowIndex,
    FlowPlacement,
    FlowPolicy,
    FlowPolicyConfig,
    FlowState,
    FlowSummary,
    Graph,
    MinCut,
    PathAlg,
    PredDAG,
    StrictMultiDiGraph,
)

from ._version import __version__

# Provide richer type information for editors/type-checkers without affecting runtime.
try:  # pragma: no cover - typing-only import
    from typing import TYPE_CHECKING as _TYPE_CHECKING

    if _TYPE_CHECKING:  # noqa: SIM108
        from ._docs import (  # noqa: I001
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
    "Graph",
    "EdgeSelection",
    "EdgeTieBreak",
    "PathAlg",
    "FlowPlacement",
    "FlowPolicy",
    "FlowPolicyConfig",
    "PredDAG",
    "FlowIndex",
    "FlowState",
    "MinCut",
    "FlowSummary",
    "Backend",
    "Algorithms",
]
