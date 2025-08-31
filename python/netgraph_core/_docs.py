"""Typing and docstrings for netgraph_core public API.

This module is intentionally light; runtime implementations live in the
compiled extension. The goal is to provide type hints and help text.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import TYPE_CHECKING, Optional

if TYPE_CHECKING:  # only for typing; runtime comes from extension
    from . import StrictMultiDiGraph  # type: ignore
import numpy as np


class EdgeSelect(Enum):
    """Edge selection policies for shortest-path traversal.

    Mirrors netgraph.core.types.EdgeSelect (C++ enum).
    """

    ALL_MIN_COST = 1
    SINGLE_MIN_COST = 2
    ALL_MIN_COST_WITH_CAP_REMAINING = 3
    ALL_ANY_COST_WITH_CAP_REMAINING = 4
    SINGLE_MIN_COST_WITH_CAP_REMAINING = 5
    SINGLE_MIN_COST_WITH_CAP_REMAINING_LOAD_FACTORED = 6
    USER_DEFINED = 99


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


@dataclass(frozen=True)
class Path:
    nodes: np.ndarray
    edges: np.ndarray
    cost: float


def spf(
    g: "StrictMultiDiGraph",
    src: int,
    dst: Optional[int] = None,
    *,
    edge_select: EdgeSelect = EdgeSelect.ALL_MIN_COST,
    multipath: bool = True,
    node_mask: Optional[np.ndarray] = None,
    edge_mask: Optional[np.ndarray] = None,
    eps: float = 1e-10,
):
    """Shortest paths wrapper.

    Returns: (dist: float64[N], pred: PredDAG)
    """
    ...


def spf_residual(
    g: "StrictMultiDiGraph",
    src: int,
    dst: Optional[int] = None,
    *,
    edge_select: EdgeSelect = EdgeSelect.ALL_MIN_COST_WITH_CAP_REMAINING,
    multipath: bool = True,
    eps: float = 1e-12,
    residual: np.ndarray | list[float] = ...,
    node_mask: Optional[np.ndarray] = None,
    edge_mask: Optional[np.ndarray] = None,
):
    """Residual-aware shortest paths wrapper.

    Uses per-edge residual capacities to filter/select edges.
    Returns: (dist: float64[N], pred: PredDAG)
    """
    ...


def ksp(
    g: "StrictMultiDiGraph",
    src: int,
    dst: int,
    *,
    k: int,
    max_cost_factor: Optional[float] = None,
    unique: bool = True,
    node_mask: Optional[np.ndarray] = None,
    edge_mask: Optional[np.ndarray] = None,
    eps: float = 1e-10,
): ...


def calc_max_flow(
    g: "StrictMultiDiGraph",
    src: int,
    dst: int,
    *,
    flow_placement: FlowPlacement = FlowPlacement.PROPORTIONAL,
    shortest_path: bool = False,
    eps: float = 1e-10,
    with_edge_flows: bool = False,
    node_mask: Optional[np.ndarray] = None,
    edge_mask: Optional[np.ndarray] = None,
): ...


@dataclass(frozen=True)
class CostBucket:
    cost: float
    share: float


@dataclass(frozen=True)
class CostDistribution:
    buckets: list[CostBucket]


@dataclass(frozen=True)
class MinCut:
    edges: list[int]


@dataclass(frozen=True)
class FlowSummary:
    total_flow: float
    min_cut: MinCut
    cost_distribution: CostDistribution
    edge_flows: list[float]


def batch_max_flow(
    g: "StrictMultiDiGraph",
    pairs: np.ndarray,
    *,
    flow_placement: FlowPlacement = FlowPlacement.PROPORTIONAL,
    shortest_path: bool = False,
    eps: float = 1e-10,
    with_edge_flows: bool = False,
    threads: Optional[int] = None,
    seed: Optional[int] = None,
    node_masks: Optional[list[np.ndarray]] = None,
    edge_masks: Optional[list[np.ndarray]] = None,
) -> list[FlowSummary]: ...
