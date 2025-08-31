# Re-export types from extension
# pyright: reportMissingImports=false
from _netgraph_core import (
    CostBucket,
    CostDistribution,
    EdgeSelect,
    FlowPlacement,
    FlowSummary,
    MinCut,
    PredDAG,
    StrictMultiDiGraph,
    batch_max_flow,
    calc_max_flow,
    ksp,
    max_flow,
    spf,
    spf_residual,
)

from ._version import __version__

# Provide richer type information for editors/type-checkers without affecting runtime.
try:  # pragma: no cover - typing-only import
    from typing import TYPE_CHECKING as _TYPE_CHECKING

    if _TYPE_CHECKING:  # noqa: SIM108
        from ._docs import (  # noqa: I001
            CostBucket as CostBucket,
            CostDistribution as CostDistribution,
            EdgeSelect as EdgeSelect,
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
    "spf",
    "ksp",
    "calc_max_flow",
    "max_flow",
    "spf_residual",
    "batch_max_flow",
    "EdgeSelect",
    "FlowPlacement",
    "PredDAG",
    "MinCut",
    "CostBucket",
    "CostDistribution",
    "FlowSummary",
]
