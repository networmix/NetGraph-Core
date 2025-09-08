"""Typing and docstrings for netgraph_core public API.

This module is intentionally light; runtime implementations live in the
compiled extension. The goal is to provide type hints and help text.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import TYPE_CHECKING, Optional

if TYPE_CHECKING:  # only for typing; runtime comes from extension
    import numpy as np  # type: ignore[reportMissingImports]
    # Use forward references (strings) below; avoid importing from the package
    # to prevent circular imports during type checking.


class EdgeTieBreak(Enum):
    DETERMINISTIC = 1
    PREFER_HIGHER_RESIDUAL = 2


@dataclass
class EdgeSelection:
    multi_edge: bool = True
    require_capacity: bool = False
    tie_break: EdgeTieBreak = EdgeTieBreak.DETERMINISTIC


class FlowPlacement(Enum):
    """How to place flow across equal-cost predecessors during augmentation."""

    PROPORTIONAL = 1
    EQUAL_BALANCED = 2


class PredDAG:
    """Compact predecessor DAG representation.

    Arrays are int32; offsets has length N+1.
    """

    parent_offsets: np.ndarray
    parents: np.ndarray
    via_edges: np.ndarray

    def resolve_to_paths(
        self,
        src: int,
        dst: int,
        *,
        split_parallel_edges: bool = False,
        max_paths: Optional[int] = None,
    ) -> list[tuple[tuple[int, tuple[int, ...]], ...]]: ...


class PathAlg(Enum):
    SPF = 1


class FlowGraph:
    """Opaque graph handle provided by the runtime extension (typing stub only)."""

    ...


class StrictMultiDiGraph:
    """Opaque graph structure provided by the runtime extension (typing stub only)."""

    ...


class FlowIndex:
    src: int
    dst: int
    flowClass: int
    flowId: int


class FlowPolicy:
    def __init__(
        self,
        path_alg: PathAlg = PathAlg.SPF,
        flow_placement: FlowPlacement = FlowPlacement.PROPORTIONAL,
        selection: EdgeSelection = ...,  # default provided by runtime binding
        min_flow_count: int = 1,
        max_flow_count: Optional[int] = None,
        max_path_cost: Optional[int] = None,
        max_path_cost_factor: Optional[float] = None,
        shortest_path: bool = False,
        reoptimize_flows_on_each_placement: bool = False,
        max_no_progress_iterations: int = 100,
        max_total_iterations: int = 10000,
        diminishing_returns_enabled: bool = True,
        diminishing_returns_window: int = 8,
        diminishing_returns_epsilon_frac: float = 1e-3,
    ) -> None: ...

    def flow_count(self) -> int: ...
    def placed_demand(self) -> float: ...
    def place_demand(
        self,
        flow_graph: "FlowGraph",
        src: int,
        dst: int,
        flowClass: int,
        volume: float,
        target_per_flow: Optional[float] = None,
        min_flow: Optional[float] = None,
    ) -> tuple[float, float]: ...

    def rebalance_demand(
        self,
        flow_graph: "FlowGraph",
        src: int,
        dst: int,
        flowClass: int,
        target: float,
    ) -> tuple[float, float]: ...

    def remove_demand(self, flow_graph: "FlowGraph") -> None: ...

    @property
    def flows(self) -> dict[tuple[int, int, int, int], tuple[int, int, int, float]]: ...


@dataclass(frozen=True)
class Path:
    nodes: np.ndarray
    edges: np.ndarray
    cost: float


class MinCut:
    edges: list[int]


@dataclass(frozen=True)
class FlowSummary:
    total_flow: float
    min_cut: MinCut
    costs: "np.ndarray"  # int64[K]
    flows: "np.ndarray"  # float64[K]
    edge_flows: "np.ndarray"
    residual_capacity: "np.ndarray"
    reachable_nodes: "np.ndarray"


class Algorithms:
    def spf(
        self,
        graph: "FlowGraph",
        src: int,
        *,
        dst: Optional[int] = None,
        selection: Optional[EdgeSelection] = None,
        residual: Optional["np.ndarray"] = None,
        node_mask: Optional["np.ndarray"] = None,
        edge_mask: Optional["np.ndarray"] = None,
    ): ...
    def ksp(
        self,
        graph: "FlowGraph",
        src: int,
        dst: int,
        *,
        k: int,
        max_cost_factor: Optional[float] = None,
        unique: bool = True,
        node_mask: Optional["np.ndarray"] = None,
        edge_mask: Optional["np.ndarray"] = None,
    ): ...
    def max_flow(
        self,
        graph: "FlowGraph",
        src: int,
        dst: int,
        *,
        flow_placement: FlowPlacement = FlowPlacement.PROPORTIONAL,
        shortest_path: bool = False,
        with_edge_flows: bool = False,
        with_reachable: bool = False,
        with_residuals: bool = False,
        node_mask: Optional["np.ndarray"] = None,
        edge_mask: Optional["np.ndarray"] = None,
    ): ...
    def batch_max_flow(
        self,
        graph: "FlowGraph",
        pairs: "np.ndarray",
        *,
        node_masks: Optional[list["np.ndarray"]] = None,
        edge_masks: Optional[list["np.ndarray"]] = None,
        flow_placement: FlowPlacement = FlowPlacement.PROPORTIONAL,
        shortest_path: bool = False,
        with_edge_flows: bool = False,
        with_reachable: bool = False,
        with_residuals: bool = False,
    ) -> list[FlowSummary]: ...
