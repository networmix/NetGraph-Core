"""Runtime package surface for :mod:`netgraph_core`.

All public symbols are re-exported from the compiled ``_netgraph_core`` module so
downstream users import everything from a single place. Typing shims supplied in
``_docs.py`` keep editors happy without affecting runtime behavior.

The streamlined Core is being rebuilt in sub-PRs (``dev/design/2026-10-streamline``):
this build carries the graph on integer units (``StrictMultiDiGraph``,
``Eligibility`` with ``intersect``, ``to_units``, ``as_float``), the route structures
(``PredDAG``, ``Path``), the selection and placement axes (``PathSelection``,
``EdgeTieBreak``, ``EdgePriority``, ``FlowPlacement``, ``Admission``) and the integer
ledger (``FlowGraph`` with its ``FlowRecord`` snapshots), plus the private entries
of the flow core, the kernels and the internal search in ``_netgraph_core``; the
policies, the public searches and the evaluators return in the following sub-PRs.
"""

# pyright: reportMissingImports=false

from importlib.metadata import version
from typing import TYPE_CHECKING

from _netgraph_core import (
    Admission,
    EdgePriority,
    EdgeTieBreak,
    Eligibility,
    FlowGraph,
    FlowPlacement,
    FlowRecord,
    Path,
    PathSelection,
    PredDAG,
    StrictMultiDiGraph,
    as_float,
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
        Admission as Admission,
        EdgePriority as EdgePriority,
        EdgeTieBreak as EdgeTieBreak,
        Eligibility as Eligibility,
        FlowGraph as FlowGraph,
        FlowPlacement as FlowPlacement,
        FlowRecord as FlowRecord,
        Path as Path,
        PathSelection as PathSelection,
        PredDAG as PredDAG,
        StrictMultiDiGraph as StrictMultiDiGraph,
        as_float as as_float,
        to_units as to_units,
    )

__all__ = [
    "__version__",
    "StrictMultiDiGraph",
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
    "FlowGraph",
    "FlowRecord",
]
