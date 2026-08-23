"""Runtime package surface for :mod:`netgraph_core`.

All public symbols are re-exported from the compiled ``_netgraph_core`` module so
downstream users import everything from a single place. Typing shims supplied in
``_docs.py`` keep editors happy without affecting runtime behavior.
"""

# pyright: reportMissingImports=false

from importlib.metadata import version
from typing import TYPE_CHECKING

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
    profiling_dump,
    profiling_enabled,
    profiling_reset,
)

__version__ = version("netgraph-core")

# Provide richer type information for editors/type-checkers without affecting runtime.
# Every name here must also appear in the runtime import above; the aliases only
# change what a type checker sees. Keep this list in sync with _docs.py, and keep
# _docs.py in sync with the bindings -- a wrong stub is worse than a missing one,
# because it turns "unchecked" into "confidently wrong".
if TYPE_CHECKING:  # pragma: no cover - typing-only
    from ._docs import (  # noqa: I001
        Algorithms as Algorithms,
        Backend as Backend,
        EdgeSelection as EdgeSelection,
        EdgeTieBreak as EdgeTieBreak,
        FlowGraph as FlowGraph,
        FlowIndex as FlowIndex,
        FlowPlacement as FlowPlacement,
        FlowPolicy as FlowPolicy,
        FlowPolicyConfig as FlowPolicyConfig,
        FlowState as FlowState,
        FlowSummary as FlowSummary,
        Graph as Graph,
        MinCut as MinCut,
        PathAlg as PathAlg,
        PredDAG as PredDAG,
        StrictMultiDiGraph as StrictMultiDiGraph,
        profiling_dump as profiling_dump,
        profiling_enabled as profiling_enabled,
        profiling_reset as profiling_reset,
    )

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
    "profiling_enabled",
    "profiling_dump",
    "profiling_reset",
]
