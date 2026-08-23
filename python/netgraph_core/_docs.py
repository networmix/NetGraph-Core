"""Typing and docstrings for netgraph_core public API.

This module is intentionally light; runtime implementations live in the
compiled extension. The goal is to provide type hints and help text.

TODO: replace this hand-written stub with generated `.pyi` (e.g. via
pybind11-stubgen) to prevent drift from the runtime bindings.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, ClassVar, Optional

if TYPE_CHECKING:  # only for typing; runtime comes from extension
    import numpy as np  # type: ignore[reportMissingImports]
    # Use forward references (strings) below; avoid importing from the package
    # to prevent circular imports during type checking.


class EdgeTieBreak:
    """Tie-break rule among equal-cost parallel edges (pybind11 enum)."""

    DETERMINISTIC: ClassVar[EdgeTieBreak]
    PREFER_HIGHER_RESIDUAL: ClassVar[EdgeTieBreak]
    __members__: ClassVar[dict[str, EdgeTieBreak]]

    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...


class EdgeSelection:
    """Edge selection policy. The constructor is keyword-only."""

    multi_edge: bool
    require_capacity: bool
    tie_break: EdgeTieBreak

    def __init__(
        self,
        *,
        multi_edge: bool = True,
        require_capacity: bool = False,
        tie_break: EdgeTieBreak = ...,
    ) -> None: ...


class FlowPlacement:
    """How to place flow across equal-cost predecessors during augmentation.

    PROPORTIONAL (WCMP-like): Distributes flow proportionally to available capacity.
        May be used iteratively (e.g., for max-flow).

    EQUAL_BALANCED (ECMP): Single-pass admission on a fixed shortest-path DAG (Dijkstra).
        Computes one global scale so no edge is oversubscribed under equal
        per-edge splits, places once, and stops. Re-invoking on updated
        residuals changes the next-hop set (progressive traffic-engineering behavior).
        ECMP = Equal-Cost Multi-Path; WCMP = Weighted-Cost Multi-Path.
    """

    PROPORTIONAL: ClassVar[FlowPlacement]
    EQUAL_BALANCED: ClassVar[FlowPlacement]
    __members__: ClassVar[dict[str, FlowPlacement]]

    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...


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


class PathAlg:
    SPF: ClassVar[PathAlg]
    __members__: ClassVar[dict[str, PathAlg]]

    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...


class Backend:
    """Backend for algorithm execution (typing stub only)."""

    @staticmethod
    def cpu() -> "Backend":
        """Create a CPU backend."""
        ...


class Graph:
    """Opaque graph handle provided by the runtime extension (typing stub only)."""

    ...


class FlowGraph:
    """Flow graph for multi-commodity flow tracking.

    Views return read-only NumPy arrays referencing internal C++ memory.
    The FlowGraph must remain alive while views are in use.
    """

    def __init__(self, graph: "StrictMultiDiGraph") -> None: ...

    def capacity_view(self) -> "np.ndarray":
        """Return read-only float64 view of edge capacities."""
        ...

    def residual_view(self) -> "np.ndarray":
        """Return read-only float64 view of residual capacities."""
        ...

    def edge_flow_view(self) -> "np.ndarray":
        """Return read-only float64 view of placed flows."""
        ...

    @property
    def graph(self) -> "StrictMultiDiGraph": ...

    def place(
        self,
        index: "FlowIndex",
        src: int,
        dst: int,
        dag: "PredDAG",
        amount: float,
        flow_placement: FlowPlacement = ...,
    ) -> float: ...

    def remove(self, index: "FlowIndex") -> None: ...
    def remove_by_class(self, cls: int) -> None: ...
    def reset(self) -> None: ...
    def get_flow_edges(self, index: "FlowIndex") -> list[tuple[int, float]]: ...
    def get_flow_path(self, index: "FlowIndex") -> list[int]: ...


class FlowState:
    """Flow state tracker for a single-commodity flow problem.

    Tracks per-edge residual capacity and flow placement for a graph.

    **Lifetime & Memory Safety:**
    - Views (capacity_view, residual_view, edge_flow_view) return read-only
      NumPy arrays that directly reference internal C++ memory (zero-copy).
    - The FlowState must remain alive while these views are in use.
    - Arrays passed as 'residual' to __init__ or reset() are **copied** internally.
    - Masks passed to compute_min_cut() are **copied** before releasing GIL (thread-safe).

    **Array Requirements:**
    - All arrays must be 1-D, C-contiguous, and have correct dtype (float64 for
      residual, bool for masks).
    - Mask lengths must match num_nodes/num_edges exactly or TypeError is raised.
    """

    def __init__(
        self, graph: "StrictMultiDiGraph", residual: "Optional[np.ndarray]" = None
    ) -> None:
        """Initialize FlowState with a graph and optional residual capacities.

        Args:
            graph: The network graph (kept alive automatically)
            residual: Optional 1-D float64 array of initial residual capacities.
                     If not provided, uses graph capacities. **Array is copied.**

        Raises:
            TypeError: If residual has wrong dtype, ndim, or length.
        """
        ...

    def reset(self, residual: "Optional[np.ndarray]" = None) -> None:
        """Reset flow state to initial capacities.

        Args:
            residual: Optional 1-D float64 array. **Array is copied.**
        """
        ...

    def capacity_view(self) -> "np.ndarray":
        """Return read-only float64 view of edge capacities (zero-copy)."""
        ...

    def residual_view(self) -> "np.ndarray":
        """Return read-only float64 view of residual capacities (zero-copy)."""
        ...

    def edge_flow_view(self) -> "np.ndarray":
        """Return read-only float64 view of placed flows (zero-copy)."""
        ...

    def place_on_dag(
        self,
        src: int,
        dst: int,
        dag: "PredDAG",
        requested_flow: float = float("inf"),
        flow_placement: FlowPlacement = ...,
    ) -> float:
        """Place flow along a predecessor DAG.

        EqualBalanced is **single-pass ECMP admission** on the provided DAG:
        we compute one global scale so no edge is oversubscribed under equal per-edge
        splits, apply it once, and return. Re-invoking on updated residuals changes
        the next-hop set (progressive behavior).

        Returns:
            Amount of flow actually placed (may be less than requested).
        """
        ...

    def place_max_flow(
        self,
        src: int,
        dst: int,
        flow_placement: FlowPlacement = ...,
        shortest_path: bool = False,
        require_capacity: bool = True,
        *,
        node_mask: "Optional[np.ndarray]" = None,
        edge_mask: "Optional[np.ndarray]" = None,
    ) -> float:
        """Place maximum flow from src to dst.

        Args:
            require_capacity: Whether to require edges to have capacity.
                - True (default): Routes adapt to residuals (SDN/TE, progressive fill).
                - False: Routes based on costs only (IP/IGP, fixed routing).
            node_mask: Optional 1-D bool array (True=allowed). Copied for thread safety.
            edge_mask: Optional 1-D bool array (True=allowed). Copied for thread safety.

        For *IP-style ECMP* max-flow, use: require_capacity=False, shortest_path=True,
        flow_placement=EQUAL_BALANCED.

        Caution: with require_capacity=True + EQUAL_BALANCED, this method iterates shortest-path +
        placement (progressive fill), which differs from single-pass ECMP admission.

        Returns:
            Total flow placed.
        """
        ...

    def compute_min_cut(
        self,
        src: int,
        *,
        node_mask: "Optional[np.ndarray]" = None,
        edge_mask: "Optional[np.ndarray]" = None,
    ) -> "MinCut":
        """Compute minimum cut from source based on current residual.

        Args:
            src: Source node for reachability analysis
            node_mask: Optional 1-D bool array of length num_nodes. **Copied for thread safety.**
            edge_mask: Optional 1-D bool array of length num_edges. **Copied for thread safety.**

        Returns:
            MinCut object with edges in the cut set.

        Raises:
            TypeError: If masks have wrong dtype, ndim, or length.
        """
        ...


class StrictMultiDiGraph:
    """Immutable directed multigraph with CSR and reverse-CSR adjacency.

    Edges are reordered by (cost, src, dst) at construction, so EdgeId values are
    positions in that canonical order rather than the caller's input order. Use
    `ext_edge_ids` to map results back to caller-side identifiers.

    Every `*_view()` method returns a fresh, writable COPY of the internal buffer,
    so mutating the result cannot violate the graph's immutability. (This differs
    from `FlowState`/`FlowGraph` views, which are live read-only views.)
    """

    @staticmethod
    def from_arrays(
        num_nodes: int,
        src: "np.ndarray",
        dst: "np.ndarray",
        capacity: "np.ndarray",
        cost: "np.ndarray",
        ext_edge_ids: "np.ndarray",
    ) -> "StrictMultiDiGraph":
        """Build a graph from parallel edge arrays (all must be 1-D, C-contiguous).

        Args:
            num_nodes: Node count; every src/dst must lie in [0, num_nodes).
            src: int32[E] source node ids.
            dst: int32[E] destination node ids.
            capacity: float64[E] capacities, each >= 0.
            cost: int64[E] costs, each >= 0. Their TOTAL must stay below 2**62,
                since SPF accumulates path costs as int64.
            ext_edge_ids: int64[E] caller-side ids, permuted alongside the edges.

        Raises:
            TypeError: On a wrong dtype, a non-contiguous array, or a wrong shape.
            ValueError: On a negative capacity/cost, a node id out of range, a
                length mismatch, or a total cost at or above 2**62.
        """
        ...

    def num_nodes(self) -> int: ...
    def num_edges(self) -> int: ...
    def capacity_view(self) -> "np.ndarray":
        """Copy of per-edge capacities, float64[E]."""
        ...

    def cost_view(self) -> "np.ndarray":
        """Copy of per-edge costs, int64[E]."""
        ...

    def edge_src_view(self) -> "np.ndarray":
        """Copy of per-edge source node ids, int32[E]."""
        ...

    def edge_dst_view(self) -> "np.ndarray":
        """Copy of per-edge destination node ids, int32[E]."""
        ...

    def ext_edge_ids_view(self) -> "np.ndarray":
        """Copy of caller-supplied external edge ids, int64[E]."""
        ...

    def row_offsets_view(self) -> "np.ndarray":
        """CSR row offsets over outgoing edges, int32[num_nodes + 1]."""
        ...

    def col_indices_view(self) -> "np.ndarray":
        """CSR neighbour node ids for outgoing edges, int32[E]."""
        ...

    def adj_edge_index_view(self) -> "np.ndarray":
        """EdgeId for each CSR outgoing entry, int32[E]."""
        ...

    def in_row_offsets_view(self) -> "np.ndarray":
        """Reverse-CSR row offsets over incoming edges, int32[num_nodes + 1]."""
        ...

    def in_col_indices_view(self) -> "np.ndarray":
        """Reverse-CSR predecessor node ids, int32[E]."""
        ...

    def in_adj_edge_index_view(self) -> "np.ndarray":
        """EdgeId for each reverse-CSR entry, int32[E]."""
        ...


class FlowIndex:
    src: int
    dst: int
    flowClass: int
    flowId: int


class FlowPolicyConfig:
    """Configuration for FlowPolicy behavior.

    Key Parameters:
        multipath: Controls whether individual flows split across multiple equal-cost paths.
            - True (default): Hash-based ECMP - each flow uses a DAG with ALL equal-cost paths.
              Flow volume is split across these paths according to flow_placement strategy.
              This models router ECMP behavior where packets are hashed across next-hops.
            - False: Tunnel-based ECMP - each flow uses a SINGLE path (one tunnel/LSP).
              Multiple flows can share the same path. This models MPLS LSP semantics where
              each LSP follows one specific path, and ECMP means balancing ACROSS LSPs.
    """

    path_alg: PathAlg
    flow_placement: FlowPlacement
    selection: EdgeSelection
    require_capacity: bool
    multipath: bool
    min_flow_count: int
    max_flow_count: Optional[int]
    max_path_cost: Optional[int]
    max_path_cost_factor: Optional[float]
    shortest_path: bool
    reoptimize_flows_on_each_placement: bool
    max_no_progress_iterations: int
    max_total_iterations: int
    diminishing_returns_enabled: bool
    diminishing_returns_window: int
    diminishing_returns_epsilon_frac: float


class FlowPolicy:
    """Flow policy for demand placement.

    When static_paths is empty the policy may refresh the DAG per round using
    residual-aware shortest paths. This progressively prunes saturated next-hops
    (traffic-engineering style) and differs from one-shot ECMP admission.

    Args:
        algorithms: Algorithms instance (kept alive by FlowPolicy)
        graph: Graph handle (kept alive by FlowPolicy)
        config: FlowPolicyConfig with all policy parameters
        node_mask: Optional 1-D bool array (True=allowed). **Copied for thread safety.**
        edge_mask: Optional 1-D bool array (True=allowed). **Copied for thread safety.**

    Raises:
        TypeError: If masks have wrong dtype, ndim, or length.
    """

    def __init__(
        self,
        algorithms: "Algorithms",
        graph: "Graph",
        config: "FlowPolicyConfig",
        *,
        node_mask: "Optional[np.ndarray]" = None,
        edge_mask: "Optional[np.ndarray]" = None,
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


class MinCut:
    edges: "np.ndarray"  # int32[K], copy of the cut edge ids


class FlowSummary:
    total_flow: float
    min_cut: MinCut
    costs: "np.ndarray"  # int64[K]
    flows: "np.ndarray"  # float64[K]
    edge_flows: "np.ndarray"
    residual_capacity: "np.ndarray"
    reachable_nodes: "np.ndarray"


class Algorithms:
    """Core graph algorithms with thread-safe bindings.

    **Thread Safety & Array Handling:**
    - All input arrays (residual, node_mask, edge_mask) are **copied** before
      releasing the GIL, ensuring thread safety even if Python code mutates
      the arrays concurrently.
    - All arrays must be 1-D, C-contiguous, with correct dtype.
    - Mask lengths must match num_nodes/num_edges exactly or TypeError is raised.

    **Lifetime Management:**
    - Algorithms instances must be kept alive while Graph handles derived from
      them are in use (automatic via Python reference counting in normal usage).
    """

    def __init__(self, backend: "Backend") -> None: ...

    def build_graph(self, graph: "StrictMultiDiGraph") -> "Graph":
        """Build a graph handle from StrictMultiDiGraph.

        The graph object is kept alive automatically.
        """
        ...

    def build_graph_from_arrays(
        self,
        num_nodes: int,
        src: "np.ndarray",
        dst: "np.ndarray",
        capacity: "np.ndarray",
        cost: "np.ndarray",
        ext_edge_ids: "np.ndarray",
    ) -> "Graph":
        """Build graph directly from arrays (graph is owned by the handle)."""
        ...

    def spf(
        self,
        graph: "Graph",
        src: int,
        dst: Optional[int] = None,
        *,
        selection: Optional[EdgeSelection] = None,
        residual: Optional["np.ndarray"] = None,
        node_mask: Optional["np.ndarray"] = None,
        edge_mask: Optional["np.ndarray"] = None,
        multipath: bool = True,
        dtype: str = "float64",
    ) -> tuple["np.ndarray", "PredDAG"]:
        """Shortest path first algorithm.

        Args:
            graph: Graph handle
            src: Source node
            dst: Optional destination (if None, compute from source to all)
            selection: Edge selection policy
            residual: Optional 1-D float64 array of residuals. **Copied for thread safety.**
            node_mask: Optional 1-D bool mask (length num_nodes). **Copied for thread safety.**
                       True = node allowed, False = node excluded.
                       **If source node is masked (False), returns empty DAG with all distances at infinity.**
            edge_mask: Optional 1-D bool mask (length num_edges). **Copied for thread safety.**
                       True = edge allowed, False = edge excluded.
            multipath: Whether to track multiple equal-cost paths
            dtype: "float64" (inf for unreachable) or "int64" (max for unreachable)

        Returns:
            (distances, predecessor_dag)

        Note:
            When the source node is masked out (node_mask[src] == False), the algorithm
            immediately returns an empty predecessor DAG with all distances set to infinity,
            as no traversal can begin from an excluded source.

        Raises:
            TypeError: If arrays have wrong dtype, ndim, or length
                (residual must have length num_edges; masks must match
                num_nodes / num_edges).
            ValueError: If src/dst out of range.
        """
        ...

    def ksp(
        self,
        graph: "Graph",
        src: int,
        dst: int,
        *,
        k: int,
        max_cost_factor: Optional[float] = None,
        unique: bool = True,
        node_mask: Optional["np.ndarray"] = None,
        edge_mask: Optional["np.ndarray"] = None,
        dtype: str = "float64",
    ) -> list[tuple["np.ndarray", "PredDAG"]]:
        """K shortest paths algorithm.

        Args:
            k: Number of paths to find
            max_cost_factor: Optional maximum cost factor relative to shortest
            node_mask: Optional 1-D bool mask. **Copied for thread safety.**
            edge_mask: Optional 1-D bool mask. **Copied for thread safety.**

        Returns:
            List of (distances, dag) tuples, up to k paths.

        Note:
            Tie-breaking is deterministic by edge ID. EdgeSelection policy is not used.
        """
        ...

    def max_flow(
        self,
        graph: "Graph",
        src: int,
        dst: int,
        *,
        flow_placement: FlowPlacement = ...,
        shortest_path: bool = False,
        require_capacity: bool = True,
        with_edge_flows: bool = False,
        with_reachable: bool = False,
        with_residuals: bool = False,
        node_mask: Optional["np.ndarray"] = None,
        edge_mask: Optional["np.ndarray"] = None,
    ) -> tuple[float, "FlowSummary"]:
        """Maximum flow algorithm.

        Behavior depends on require_capacity parameter:

        - require_capacity=True (default): "True max-flow" - require edges to have capacity,
          exclude saturated links (SDN/TE behavior). Uses all available paths.

        - require_capacity=False: "IP ECMP" - routes based on costs only (IP/IGP behavior).
          Routes computed once from topology, ignoring capacity. Models OSPF/IS-IS routing.

        Common configurations:

        1. True max-flow (SDN/TE):
           require_capacity=True, shortest_path=False

        2. IP ECMP max-flow:
           require_capacity=False, shortest_path=True, flow_placement=EQUAL_BALANCED

        3. IP WCMP max-flow:
           require_capacity=False, shortest_path=True, flow_placement=PROPORTIONAL

        Args:
            require_capacity: Whether to require edges to have capacity.
            shortest_path: Whether to restrict to lowest-cost tier only.
            flow_placement: How to split flow across parallel edges (ECMP vs WCMP).
            node_mask: Optional 1-D bool mask. **Copied for thread safety.**
            edge_mask: Optional 1-D bool mask. **Copied for thread safety.**

        Returns:
            (total_flow, flow_summary)
        """
        ...

    def batch_max_flow(
        self,
        graph: "Graph",
        pairs: "np.ndarray",
        *,
        node_masks: Optional[list["np.ndarray"]] = None,
        edge_masks: Optional[list["np.ndarray"]] = None,
        flow_placement: FlowPlacement = ...,
        shortest_path: bool = False,
        require_capacity: bool = True,
        with_edge_flows: bool = False,
        with_reachable: bool = False,
        with_residuals: bool = False,
    ) -> list[FlowSummary]:
        """Batch maximum flow for multiple (src,dst) pairs.

        Args:
            pairs: int32 array of shape [B, 2] with (src, dst) pairs
            node_masks: Optional list of B bool masks. **Each copied for thread safety.**
            edge_masks: Optional list of B bool masks. **Each copied for thread safety.**
            require_capacity: If True, exclude saturated edges (SDN/TE). If False, route by cost only (IP/IGP).

        Returns:
            List of FlowSummary objects, one per pair.
        """
        ...

    def sensitivity_analysis(
        self,
        graph: "Graph",
        src: int,
        dst: int,
        *,
        flow_placement: FlowPlacement = ...,
        shortest_path: bool = False,
        require_capacity: bool = True,
        node_mask: Optional["np.ndarray"] = None,
        edge_mask: Optional["np.ndarray"] = None,
    ) -> list[tuple[int, float]]:
        """Sensitivity analysis to identify critical edges that constrain flow.

        Computes baseline flow, then tests removing each saturated edge to
        measure how much the total flow would be reduced.

        The `shortest_path` parameter controls the routing semantics:

        - shortest_path=False (default): Full max-flow analysis (SDN/TE mode).
          Identifies edges critical for achieving maximum possible flow.

        - shortest_path=True: Shortest-path-only analysis (IP/IGP mode).
          Identifies edges critical for flow under ECMP routing. Edges on
          unused longer paths are not reported as critical.

        Args:
            graph: Graph handle
            src: Source node
            dst: Destination node
            flow_placement: Flow placement strategy
            shortest_path: If True, use single-pass shortest-path flow (IP/IGP).
                          If False, use full iterative max-flow (SDN/TE).
            require_capacity: If True, exclude saturated edges from routing.
            node_mask: Optional 1-D bool mask. **Copied for thread safety.**
            edge_mask: Optional 1-D bool mask. **Copied for thread safety.**

        Returns:
            List of (edge_id, flow_delta) tuples. Each edge_id is an edge
            whose removal would reduce total flow by flow_delta.
        """
        ...
