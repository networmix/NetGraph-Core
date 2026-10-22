# Reviewed overrides of the generated stub python/netgraph_core/__init__.pyi (applied by
# dev/stubs/generate.py; regenerate with bash dev/stubs/run.sh generate).
#
# pybind11-stubgen can only write what the compiled signatures show. The bindings take most
# integer arguments as Python objects (to validate exact, unnarrowed values themselves), return
# NumPy arrays without a declared dtype and optional values as plain objects, so the raw stub
# says typing.Any, numpy.ndarray or a bare list there. Each override below restores the domain
# the binding code accepts and the type it returns, read from bindings/python/*.cpp and checked
# on the built module.
#
# Rules (enforced by the generator): an override changes annotations only (its parameter names,
# kinds and defaults are the generated ones); "# reason:" says why; "# generated:" quotes the
# exact generated line it replaces, and the generation fails when that line changes, so a
# changed binding is reviewed again. Overloads replace one generated declaration. The one
# addition is the marker of a class the module refuses to construct ("No constructor defined"):
# __init__(self, _no_public_constructor: typing.Never, /), so a type checker rejects Cls().
#
# The domains, written out in full (a stub has no private aliases):
# - an integer: typing.SupportsIndex, a Python int or a NumPy integer scalar (any object with
#   __index__; a bool is rejected at run time, which a type cannot express);
# - integers: collections.abc.Sequence[typing.SupportsIndex] | numpy.typing.NDArray[numpy.integer],
#   one-dimensional (the rank and the values are checked at run time);
# - a node set: one integer or integers;
# - a physical number or resolution: int | float | numpy.integer | numpy.floating.
# Arrays that the module returns carry their exact dtype; their shape (one-dimensional) and
# whether they are read-only views are in the docstrings.
#
# The overrides name the classes of the generated stub, which this file does not declare.
# ruff: noqa: F821

import collections.abc
import fractions
import typing

import numpy
import numpy.typing

class Eligibility:
    # reason: no constructor: built by Graph.eligibility and Eligibility.intersect
    # generated: (none)
    def __init__(self, _no_public_constructor: typing.Never, /) -> None: ...
    # reason: a read-only view of the bool flags, one per edge
    # generated: def edge_allowed(self) -> numpy.ndarray:
    def edge_allowed(self) -> numpy.typing.NDArray[numpy.bool_]: ...
    # reason: a read-only view of the bool flags, one per node
    # generated: def node_allowed(self) -> numpy.ndarray:
    def node_allowed(self) -> numpy.typing.NDArray[numpy.bool_]: ...

class Graph:
    # reason: num_nodes is an integer, the four edge arrays are integers (exact values checked by the binding)
    # generated: def __init__(self, num_nodes: typing.Any, edge_sources: typing.Any, edge_destinations: typing.Any, edge_capacity_units: typing.Any, edge_costs: typing.Any, resolution: int | float | numpy.integer | numpy.floating) -> None:
    def __init__(
        self,
        num_nodes: typing.SupportsIndex,
        edge_sources: collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer],
        edge_destinations: collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer],
        edge_capacity_units: collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer],
        edge_costs: collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer],
        resolution: int | float | numpy.integer | numpy.floating,
    ) -> None: ...
    # reason: a read-only int64 view
    # generated: def edge_capacity_units(self) -> numpy.ndarray:
    def edge_capacity_units(self) -> numpy.typing.NDArray[numpy.int64]: ...
    # reason: a read-only int64 view
    # generated: def edge_costs(self) -> numpy.ndarray:
    def edge_costs(self) -> numpy.typing.NDArray[numpy.int64]: ...
    # reason: a read-only int32 view
    # generated: def edge_destinations(self) -> numpy.ndarray:
    def edge_destinations(self) -> numpy.typing.NDArray[numpy.int32]: ...
    # reason: a read-only int32 view
    # generated: def edge_sources(self) -> numpy.ndarray:
    def edge_sources(self) -> numpy.typing.NDArray[numpy.int32]: ...
    # reason: the excluded node ids and EdgeIds are integers
    # generated: def eligibility(self, excluded_nodes: typing.Any = (), excluded_edges: typing.Any = ()) -> Eligibility:
    def eligibility(
        self,
        excluded_nodes: collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer] = (),
        excluded_edges: collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer] = (),
    ) -> Eligibility: ...
    # reason: a read-only int32 view
    # generated: def incoming_edges(self) -> numpy.ndarray:
    def incoming_edges(self) -> numpy.typing.NDArray[numpy.int32]: ...
    # reason: a read-only int32 view
    # generated: def incoming_offsets(self) -> numpy.ndarray:
    def incoming_offsets(self) -> numpy.typing.NDArray[numpy.int32]: ...
    # reason: a read-only int32 view
    # generated: def input_positions(self) -> numpy.ndarray:
    def input_positions(self) -> numpy.typing.NDArray[numpy.int32]: ...
    # reason: a read-only int32 view
    # generated: def outgoing_offsets(self) -> numpy.ndarray:
    def outgoing_offsets(self) -> numpy.typing.NDArray[numpy.int32]: ...

class PredecessorDAG:
    # reason: no constructor: built by PredecessorDAG.from_csr, Path.to_dag and the searches
    # generated: (none)
    def __init__(self, _no_public_constructor: typing.Never, /) -> None: ...
    # reason: parent_offsets, parents and via_edges are integers
    # generated: def from_csr(graph: Graph, parent_offsets: typing.Any, parents: typing.Any, via_edges: typing.Any) -> PredecessorDAG:
    @staticmethod
    def from_csr(
        graph: Graph,
        parent_offsets: collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer],
        parents: collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer],
        via_edges: collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer],
    ) -> PredecessorDAG: ...
    # reason: a read-only int32 view
    # generated: def parents(self) -> numpy.ndarray:
    def parents(self) -> numpy.typing.NDArray[numpy.int32]: ...
    # reason: a read-only int32 view
    # generated: def via_edges(self) -> numpy.ndarray:
    def via_edges(self) -> numpy.typing.NDArray[numpy.int32]: ...

class Path:
    # reason: hops is a sequence of hops, each one integers (EdgeIds); a two-dimensional integer
    # array is such a sequence too (one row per hop)
    # generated: def __init__(self, graph: Graph, hops: typing.Any) -> None:
    def __init__(
        self,
        graph: Graph,
        hops: collections.abc.Sequence[
            collections.abc.Sequence[typing.SupportsIndex]
            | numpy.typing.NDArray[numpy.integer]
        ]
        | numpy.typing.NDArray[numpy.integer],
    ) -> None: ...

class FlowRecord:
    # reason: no constructor: snapshots come from FlowLedger.records/record and FlowPolicy.records/record
    # generated: (none)
    def __init__(self, _no_public_constructor: typing.Never, /) -> None: ...
    # reason: a list of (EdgeId, units) pairs; units are exact Python ints
    # generated: def edge_drops(self) -> list:
    @property
    def edge_drops(self) -> list[tuple[int, int]]: ...
    # reason: a list of (EdgeId, units) pairs
    # generated: def edge_flows(self) -> list:
    @property
    def edge_flows(self) -> list[tuple[int, int]]: ...
    # reason: a list of (EdgeId, units) pairs
    # generated: def edge_reservations(self) -> list:
    @property
    def edge_reservations(self) -> list[tuple[int, int]]: ...
    # reason: an exact Python int for a Lossy flow, None otherwise
    # generated: def offered(self) -> typing.Any:
    @property
    def offered(self) -> int | None: ...
    # reason: an int, None when the route is empty
    # generated: def min_route_cost(self) -> typing.Any:
    @property
    def min_route_cost(self) -> int | None: ...

class FlowLedger:
    # reason: the service and reservation baselines are integers, or None for the capacities
    # generated: def __init__(self, graph: Graph, service_baseline: typing.Any = None, *, reservation_baseline: typing.Any = None) -> None:
    def __init__(
        self,
        graph: Graph,
        service_baseline: collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer]
        | None = None,
        *,
        reservation_baseline: collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer]
        | None = None,
    ) -> None: ...
    # reason: sources is a node set; eligibility may be None (all allowed)
    # generated: def saturated_boundary_edges(self, sources: typing.Any, eligibility: Eligibility = None) -> list[int]:
    def saturated_boundary_edges(
        self,
        sources: typing.SupportsIndex
        | collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer],
        eligibility: Eligibility | None = None,
    ) -> list[int]: ...
    # reason: a read-only int64 view
    # generated: def edge_capacity_units(self) -> numpy.ndarray:
    def edge_capacity_units(self) -> numpy.typing.NDArray[numpy.int64]: ...
    # reason: a FlowId is a tuple (ledger uid, id) of integers
    # generated: def record(self, flow_id: typing.Any) -> FlowRecord:
    def record(
        self, flow_id: tuple[typing.SupportsIndex, typing.SupportsIndex]
    ) -> FlowRecord: ...
    # reason: the new baselines are integers, or None for the capacities
    # generated: def reset(self, service_baseline: typing.Any = None, *, reservation_baseline: typing.Any = None) -> None:
    def reset(
        self,
        service_baseline: collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer]
        | None = None,
        *,
        reservation_baseline: collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer]
        | None = None,
    ) -> None: ...

class AbsoluteCost:
    # reason: cost is an integer
    # generated: def __init__(self, cost: typing.Any) -> None:
    def __init__(self, cost: typing.SupportsIndex) -> None: ...

class CostStretch:
    # reason: factor is a Fraction, a float, a NumPy floating scalar or an integer (a bool is rejected at run time)
    # generated: def __init__(self, factor: typing.Any) -> None:
    def __init__(
        self,
        factor: fractions.Fraction | float | numpy.floating | typing.SupportsIndex,
    ) -> None: ...
    # reason: the exact factor, a fractions.Fraction built on each access
    # generated: def factor(self) -> typing.Any:
    @property
    def factor(self) -> fractions.Fraction: ...

class TunnelCount:
    # reason: per_head is an integer
    # generated: def __init__(self, per_head: typing.Any, *, up_to: bool = False) -> None:
    def __init__(
        self, per_head: typing.SupportsIndex, *, up_to: bool = False
    ) -> None: ...

class DistributedTE:
    # reason: the ceiling is AbsoluteCost, CostStretch or None; the search minimum is an integer
    # generated: def __init__(self, tunnel_count: TunnelCount, *, cost_ceiling: typing.Any = None, min_unreserved_bandwidth: typing.Any = 0, tie_break: RouteTieBreak = ..., parallel_edges: bool = True) -> None:
    def __init__(
        self,
        tunnel_count: TunnelCount,
        *,
        cost_ceiling: AbsoluteCost | CostStretch | None = None,
        min_unreserved_bandwidth: typing.SupportsIndex = 0,
        tie_break: RouteTieBreak = ...,
        parallel_edges: bool = True,
    ) -> None: ...
    # reason: the ceiling given at construction: AbsoluteCost, CostStretch or None
    # generated: def cost_ceiling(self) -> typing.Any:
    @property
    def cost_ceiling(self) -> AbsoluteCost | CostStretch | None: ...

class Tunnel:
    # reason: no constructor: a tunnel of ForwardingProgram.tunnels
    # generated: (none)
    def __init__(self, _no_public_constructor: typing.Never, /) -> None: ...
    # reason: the path's EdgeIds in path order, None for a tunnel that is down
    # generated: def edges(self) -> typing.Any:
    @property
    def edges(self) -> tuple[int, ...] | None: ...

class RouteDefinition:
    # reason: head and weight are integers (weight checked by from_routes); route a Path or a PredecessorDAG (TypeError otherwise)
    # generated: def __init__(self, head: typing.Any, route: typing.Any, *, weight: typing.Any = 1) -> None:
    def __init__(
        self,
        head: typing.SupportsIndex,
        route: Path | PredecessorDAG,
        *,
        weight: typing.SupportsIndex = 1,
    ) -> None: ...

class ForwardingProgram:
    # reason: no constructor: made by from_routes or from_dag and returned by the policies
    # generated: (none)
    def __init__(self, _no_public_constructor: typing.Never, /) -> None: ...
    # reason: sources and destinations are node sets; routes a sequence of RouteDefinition; edge_weights and residual integers or None
    # generated: def from_routes(graph: Graph, sources: typing.Any, destinations: typing.Any, routes: typing.Any, *, allocation: RouteAllocation | None = None, weight_basis: WeightBasis | None = None, edge_weights: typing.Any = None, residual: typing.Any = None) -> ForwardingProgram:
    @staticmethod
    def from_routes(
        graph: Graph,
        sources: typing.SupportsIndex
        | collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer],
        destinations: typing.SupportsIndex
        | collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer],
        routes: collections.abc.Sequence[RouteDefinition],
        *,
        allocation: RouteAllocation | None = None,
        weight_basis: WeightBasis | None = None,
        edge_weights: collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer]
        | None = None,
        residual: collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer]
        | None = None,
    ) -> ForwardingProgram: ...
    # reason: sources and destinations are node sets; dag a PredecessorDAG; edge_weights and residual integers or None
    # generated: def from_dag(graph: Graph, sources: typing.Any, destinations: typing.Any, dag: typing.Any, *, allocation: RouteAllocation | None = None, weight_basis: WeightBasis | None = None, edge_weights: typing.Any = None, residual: typing.Any = None) -> ForwardingProgram:
    @staticmethod
    def from_dag(
        graph: Graph,
        sources: typing.SupportsIndex
        | collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer],
        destinations: typing.SupportsIndex
        | collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer],
        dag: PredecessorDAG,
        *,
        allocation: RouteAllocation | None = None,
        weight_basis: WeightBasis | None = None,
        edge_weights: collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer]
        | None = None,
        residual: collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer]
        | None = None,
    ) -> ForwardingProgram: ...
    # reason: the destination set, ascending
    # generated: def destinations(self) -> tuple:
    @property
    def destinations(self) -> tuple[int, ...]: ...
    # reason: (source, weight) pairs, or None for equal sources
    # generated: def source_weights(self) -> typing.Any:
    @property
    def source_weights(self) -> tuple[tuple[int, int], ...] | None: ...
    # reason: the source set, ascending
    # generated: def sources(self) -> tuple:
    @property
    def sources(self) -> tuple[int, ...]: ...
    # reason: the program's tunnels, ascending index
    # generated: def tunnels(self) -> tuple:
    @property
    def tunnels(self) -> tuple[Tunnel, ...]: ...

class ShortestPaths:
    # reason: edge_weights one integer per edge (a sequence or an integer array) or None
    # generated: def __init__(self, *, allocation: RouteAllocation = ..., source_allocation: SourceAllocation | None = None, weight_basis: WeightBasis | None = None, edge_weights: typing.Any = None, tie_break: RouteTieBreak = ..., multipath: bool = True, parallel_edges: bool = True) -> None:
    def __init__(
        self,
        *,
        allocation: RouteAllocation = ...,
        source_allocation: SourceAllocation | None = None,
        weight_basis: WeightBasis | None = None,
        edge_weights: collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer]
        | None = None,
        tie_break: RouteTieBreak = ...,
        multipath: bool = True,
        parallel_edges: bool = True,
    ) -> None: ...
    # reason: the explicit weights, one int per edge, or None
    # generated: def edge_weights(self) -> typing.Any:
    @property
    def edge_weights(self) -> tuple[int, ...] | None: ...

class CentralizedTE:
    # reason: node-indexed source weights; typed support and tunnel count
    # generated: def __init__(self, *, source_allocation: SourceAllocation, source_weights: typing.Any = None, support: typing.Any = None, tunnel_count: typing.Any = None, prefer_headroom: bool = False) -> None:
    def __init__(
        self,
        *,
        source_allocation: SourceAllocation,
        source_weights: collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer]
        | None = None,
        support: PredecessorDAG | None = None,
        tunnel_count: TunnelCount | None = None,
        prefer_headroom: bool = False,
    ) -> None: ...
    # reason: owned immutable node-indexed weights
    # generated: def source_weights(self) -> typing.Any:
    @property
    def source_weights(self) -> tuple[int, ...] | None: ...

class Pinned:
    # reason: black_holes may be None
    # generated: def __init__(self, program: ForwardingProgram, *, black_holes: Eligibility = None) -> None:
    def __init__(
        self, program: ForwardingProgram, *, black_holes: Eligibility | None = None
    ) -> None: ...

class FlowPolicyConfig:
    # reason: routing is one of the four routing modes (TypeError otherwise)
    # generated: def __init__(self, *, routing: typing.Any, admission: Admission = ...) -> None:
    def __init__(
        self,
        *,
        routing: ShortestPaths | DistributedTE | CentralizedTE | Pinned,
        admission: Admission = ...,
    ) -> None: ...
    # reason: the routing mode given at construction
    # generated: def routing(self) -> typing.Any:
    @property
    def routing(self) -> ShortestPaths | DistributedTE | CentralizedTE | Pinned: ...

class SourceFlow:
    # reason: no constructor: one row of PlacementResult.by_source
    # generated: (none)
    def __init__(self, _no_public_constructor: typing.Never, /) -> None: ...
    # reason: an exact int (a 128-bit counter)
    # generated: def admitted(self) -> typing.Any:
    @property
    def admitted(self) -> int: ...
    # reason: an exact int (a 128-bit counter)
    # generated: def total_edge_drops(self) -> typing.Any:
    @property
    def total_edge_drops(self) -> int: ...
    # reason: an exact int when LOSSY, None when LOSSLESS
    # generated: def offered(self) -> typing.Any:
    @property
    def offered(self) -> int | None: ...
    # reason: an exact int (a 128-bit counter)
    # generated: def source_drops(self) -> typing.Any:
    @property
    def source_drops(self) -> int: ...

class TunnelFlow:
    # reason: no constructor: one row of PlacementResult.by_tunnel
    # generated: (none)
    def __init__(self, _no_public_constructor: typing.Never, /) -> None: ...
    # reason: an exact int (a 128-bit counter)
    # generated: def assigned(self) -> typing.Any:
    @property
    def assigned(self) -> int: ...
    # reason: an exact int (a 128-bit counter)
    # generated: def total_edge_drops(self) -> typing.Any:
    @property
    def total_edge_drops(self) -> int: ...

class PlacementResult:
    # reason: no constructor: returned by FlowPolicy.place and the evaluators
    # generated: (none)
    def __init__(self, _no_public_constructor: typing.Never, /) -> None: ...
    # reason: an exact int (a 128-bit counter)
    # generated: def admitted(self) -> typing.Any:
    @property
    def admitted(self) -> int: ...
    # reason: one SourceFlow per input source
    # generated: def by_source(self) -> tuple:
    @property
    def by_source(self) -> tuple[SourceFlow, ...]: ...
    # reason: one TunnelFlow per program tunnel
    # generated: def by_tunnel(self) -> tuple:
    @property
    def by_tunnel(self) -> tuple[TunnelFlow, ...]: ...
    # reason: a list of (EdgeId, units) pairs
    # generated: def call_edge_drops(self) -> list:
    @property
    def call_edge_drops(self) -> list[tuple[int, int]]: ...
    # reason: an exact int when LOSSY, None when LOSSLESS
    # generated: def offered(self) -> typing.Any:
    @property
    def offered(self) -> int | None: ...
    # reason: an exact int, or float('inf') after an infinite request
    # generated: def requested(self) -> typing.Any:
    @property
    def requested(self) -> int | float: ...
    # reason: an exact int (a 128-bit counter)
    # generated: def source_drops(self) -> typing.Any:
    @property
    def source_drops(self) -> int: ...
    # reason: an exact int for the controllers, None otherwise
    # generated: def total_cost(self) -> typing.Any:
    @property
    def total_cost(self) -> int | None: ...
    # reason: an exact int, None after an unlimited LOSSLESS request
    # generated: def unadmitted(self) -> typing.Any:
    @property
    def unadmitted(self) -> int | None: ...
    # reason: a list of (EdgeId, units) pairs; units are exact Python ints
    # generated: def edge_drops(self) -> list:
    @property
    def edge_drops(self) -> list[tuple[int, int]]: ...
    # reason: an int64 array when requested, else None
    # generated: def edge_flows(self) -> typing.Any:
    @property
    def edge_flows(self) -> numpy.typing.NDArray[numpy.int64] | None: ...
    # reason: an int, or float('inf') after an infinite request
    # generated: def call_undelivered(self) -> typing.Any:
    @property
    def call_undelivered(self) -> int | float: ...
    # reason: an int64 array when requested, else None
    # generated: def residual(self) -> typing.Any:
    @property
    def residual(self) -> numpy.typing.NDArray[numpy.int64] | None: ...
    # reason: an int64 array when requested, else None
    # generated: def edge_reservations(self) -> typing.Any:
    @property
    def edge_reservations(self) -> numpy.typing.NDArray[numpy.int64] | None: ...

class TunnelPlan:
    # reason: no constructor: the plan of a T1 or T2 result
    # generated: (none)
    def __init__(self, _no_public_constructor: typing.Never, /) -> None: ...
    # reason: an int, or float('inf') for an unlimited design
    # generated: def design_demand(self) -> typing.Any:
    @property
    def design_demand(self) -> int | float: ...
    # reason: an exact Fraction (T2), or None for T1
    # generated: def max_fractional_bandwidth_per_tunnel(self) -> typing.Any:
    @property
    def max_fractional_bandwidth_per_tunnel(self) -> fractions.Fraction | None: ...
    # reason: (head, share) pairs
    # generated: def design_by_head(self) -> tuple:
    @property
    def design_by_head(self) -> tuple[tuple[int, int], ...]: ...
    # reason: an exact int, possibly above 2^64 (T2 with prefer_headroom), or None
    # generated: def headroom_penalty(self) -> typing.Any:
    @property
    def headroom_penalty(self) -> int | None: ...
    # reason: the live heads, ascending (T2), or None for T1
    # generated: def heads_with_capacity(self) -> typing.Any:
    @property
    def heads_with_capacity(self) -> tuple[int, ...] | None: ...
    # reason: an exact int, possibly above 2^64 (T2), or None for T1
    # generated: def total_tunnel_cost(self) -> typing.Any:
    @property
    def total_tunnel_cost(self) -> int | None: ...

class ShortestPathCacheStats:
    # reason: no constructor: returned by ShortestPathCache.stats
    # generated: (none)
    def __init__(self, _no_public_constructor: typing.Never, /) -> None: ...

class FlowPolicy:
    # reason: sources and destinations are node sets; eligibility and cache may be None
    # generated: def __init__(self, ledger: FlowLedger, sources: typing.Any, destinations: typing.Any, config: FlowPolicyConfig, eligibility: Eligibility = None, cache: ShortestPathCache = None) -> None:
    def __init__(
        self,
        ledger: FlowLedger,
        sources: typing.SupportsIndex
        | collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer],
        destinations: typing.SupportsIndex
        | collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer],
        config: FlowPolicyConfig,
        eligibility: Eligibility | None = None,
        cache: ShortestPathCache | None = None,
    ) -> None: ...
    # reason: offer is integer units or float('inf') (a NumPy floating +inf too); outputs may be None; the result is a PlacementResult
    # generated: def place(self, offer: typing.Any, outputs: typing.Any = None) -> typing.Any:
    def place(
        self,
        offer: typing.SupportsIndex | float | numpy.floating,
        outputs: OutputOptions | None = None,
    ) -> PlacementResult: ...
    # reason: offer and design_demand are integer units or float('inf') (design_demand may be None); outputs may be None; the result is a PlacementResult
    # generated: def evaluate(self, offer: typing.Any, *, design_demand: typing.Any = None, outputs: typing.Any = None) -> typing.Any:
    def evaluate(
        self,
        offer: typing.SupportsIndex | float | numpy.floating,
        *,
        design_demand: typing.SupportsIndex | float | numpy.floating | None = None,
        outputs: OutputOptions | None = None,
    ) -> PlacementResult: ...
    # reason: offer is integer units or float('inf'); black_holes and outputs may be None; the result is a PlacementResult
    # generated: def forward(self, offer: typing.Any, *, black_holes: typing.Any = None, outputs: typing.Any = None) -> typing.Any:
    def forward(
        self,
        offer: typing.SupportsIndex | float | numpy.floating,
        *,
        black_holes: Eligibility | None = None,
        outputs: OutputOptions | None = None,
    ) -> PlacementResult: ...
    # reason: offer and design_demand are integer units or float('inf') (design_demand may be None); outputs may be None; the result is a PlacementResult
    # generated: def replace(self, offer: typing.Any, *, design_demand: typing.Any = None, outputs: typing.Any = None) -> typing.Any:
    def replace(
        self,
        offer: typing.SupportsIndex | float | numpy.floating,
        *,
        design_demand: typing.SupportsIndex | float | numpy.floating | None = None,
        outputs: OutputOptions | None = None,
    ) -> PlacementResult: ...
    # reason: a FlowId is a tuple (ledger uid, id) of integers
    # generated: def record(self, flow_id: typing.Any) -> FlowRecord:
    def record(
        self, flow_id: tuple[typing.SupportsIndex, typing.SupportsIndex]
    ) -> FlowRecord: ...

class Place:
    # reason: policy is a FlowPolicy; offer is integer units or a float (a NumPy floating value too), kept as given
    # generated: def __init__(self, policy: typing.Any, offer: typing.Any) -> None:
    def __init__(
        self, policy: FlowPolicy, offer: typing.SupportsIndex | float | numpy.floating
    ) -> None: ...
    # reason: the offer as given, an int or a float
    # generated: def offer(self) -> typing.Any:
    @property
    def offer(self) -> int | float: ...
    # reason: the FlowPolicy given
    # generated: def policy(self) -> typing.Any:
    @property
    def policy(self) -> FlowPolicy: ...

class Replace:
    # reason: policy is a FlowPolicy; offer and design_demand are integer units or a float (a NumPy floating value too), kept as given; design_demand may be None
    # generated: def __init__(self, policy: typing.Any, offer: typing.Any, *, design_demand: typing.Any = None) -> None:
    def __init__(
        self,
        policy: FlowPolicy,
        offer: typing.SupportsIndex | float | numpy.floating,
        *,
        design_demand: typing.SupportsIndex | float | numpy.floating | None = None,
    ) -> None: ...
    # reason: the design as given, an int or a float, or None
    # generated: def design_demand(self) -> typing.Any:
    @property
    def design_demand(self) -> int | float | None: ...
    # reason: the offer as given, an int or a float
    # generated: def offer(self) -> typing.Any:
    @property
    def offer(self) -> int | float: ...
    # reason: the FlowPolicy given
    # generated: def policy(self) -> typing.Any:
    @property
    def policy(self) -> FlowPolicy: ...

class Forward:
    # reason: policy is a FlowPolicy; offer is integer units or a float (a NumPy floating value too), kept as given; black_holes may be None
    # generated: def __init__(self, policy: typing.Any, offer: typing.Any, *, black_holes: typing.Any = None) -> None:
    def __init__(
        self,
        policy: FlowPolicy,
        offer: typing.SupportsIndex | float | numpy.floating,
        *,
        black_holes: Eligibility | None = None,
    ) -> None: ...
    # reason: the Eligibility given, or None
    # generated: def black_holes(self) -> typing.Any:
    @property
    def black_holes(self) -> Eligibility | None: ...
    # reason: the offer as given, an int or a float
    # generated: def offer(self) -> typing.Any:
    @property
    def offer(self) -> int | float: ...
    # reason: the FlowPolicy given
    # generated: def policy(self) -> typing.Any:
    @property
    def policy(self) -> FlowPolicy: ...

class EvaluationOptions:
    # reason: eligibility, cache, outputs, design_demand and offer may be None; design_demand and offer are integer units or a float, kept as given
    # generated: def __init__(self, config: FlowPolicyConfig, *, eligibility: Eligibility = None, cache: ShortestPathCache = None, outputs: typing.Any = None, design_demand: typing.Any = None, offer: typing.Any = None) -> None:
    def __init__(
        self,
        config: FlowPolicyConfig,
        *,
        eligibility: Eligibility | None = None,
        cache: ShortestPathCache | None = None,
        outputs: OutputOptions | None = None,
        design_demand: typing.SupportsIndex | float | numpy.floating | None = None,
        offer: typing.SupportsIndex | float | numpy.floating | None = None,
    ) -> None: ...
    # reason: the design as given, an int or a float, or None
    # generated: def design_demand(self) -> typing.Any:
    @property
    def design_demand(self) -> int | float | None: ...
    # reason: the offer as given, an int or a float, or None
    # generated: def offer(self) -> typing.Any:
    @property
    def offer(self) -> int | float | None: ...

class SensitivityEntry:
    # reason: no constructor: the entries of a SensitivityResult
    # generated: (none)
    def __init__(self, _no_public_constructor: typing.Never, /) -> None: ...

class SensitivityResult:
    # reason: no constructor: returned by evaluate_edge_failures
    # generated: (none)
    def __init__(self, _no_public_constructor: typing.Never, /) -> None: ...

class KspOptions:
    # reason: k is an integer; the ceiling is AbsoluteCost, CostStretch or None; residual is integers or None; eligibility may be None (each validated by ksp)
    # generated: def __init__(self, k: typing.Any, *, tie_break: RouteTieBreak = ..., parallel_edges: bool = True, cost_ceiling: typing.Any = None, residual: typing.Any = None, eligibility: Eligibility = None) -> None:
    def __init__(
        self,
        k: typing.SupportsIndex,
        *,
        tie_break: RouteTieBreak = ...,
        parallel_edges: bool = True,
        cost_ceiling: AbsoluteCost | CostStretch | None = None,
        residual: collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer]
        | None = None,
        eligibility: Eligibility | None = None,
    ) -> None: ...
    # reason: the k given at construction (an integer)
    # generated: def k(self) -> typing.Any:
    @property
    def k(self) -> typing.SupportsIndex: ...
    # reason: the ceiling given at construction
    # generated: def cost_ceiling(self) -> typing.Any:
    @property
    def cost_ceiling(self) -> AbsoluteCost | CostStretch | None: ...
    # reason: the residual given at construction (the same object)
    # generated: def residual(self) -> typing.Any:
    @property
    def residual(
        self,
    ) -> (
        collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer]
        | None
    ): ...

class SpfOptions:
    # reason: targets is a node set or None; residual is integers or None; min_residual is an integer; eligibility may be None (each validated by spf)
    # generated: def __init__(self, direction: Direction, *, targets: typing.Any = None, tie_break: RouteTieBreak = ..., multipath: bool = True, parallel_edges: bool = True, residual: typing.Any = None, min_residual: typing.Any = 1, eligibility: Eligibility = None) -> None:
    def __init__(
        self,
        direction: Direction,
        *,
        targets: typing.SupportsIndex
        | collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer]
        | None = None,
        tie_break: RouteTieBreak = ...,
        multipath: bool = True,
        parallel_edges: bool = True,
        residual: collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer]
        | None = None,
        min_residual: typing.SupportsIndex = 1,
        eligibility: Eligibility | None = None,
    ) -> None: ...
    # reason: the min_residual given at construction (an integer)
    # generated: def min_residual(self) -> typing.Any:
    @property
    def min_residual(self) -> typing.SupportsIndex: ...
    # reason: the residual given at construction (the same object)
    # generated: def residual(self) -> typing.Any:
    @property
    def residual(
        self,
    ) -> (
        collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer]
        | None
    ): ...
    # reason: the targets given at construction (the same object)
    # generated: def targets(self) -> typing.Any:
    @property
    def targets(
        self,
    ) -> (
        typing.SupportsIndex
        | collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer]
        | None
    ): ...

class SpfResult:
    # reason: no constructor: returned by spf
    # generated: (none)
    def __init__(self, _no_public_constructor: typing.Never, /) -> None: ...
    # reason: a read-only int64 array for a full search, None with targets
    # generated: def costs(self) -> typing.Any:
    @property
    def costs(self) -> numpy.typing.NDArray[numpy.int64] | None: ...
    # reason: always a PredecessorDAG
    # generated: def dag(self) -> typing.Any:
    @property
    def dag(self) -> PredecessorDAG: ...
    # reason: the nearest targets' cost, None when no target is reached and for a full search
    # generated: def target_cost(self) -> typing.Any:
    @property
    def target_cost(self) -> int | None: ...

# reason: vectors give float64 arrays; int and NumPy integer scalars give float; another index protocol may also be a sequence
# generated: def as_float(x: typing.Any, resolution: int | float | numpy.integer | numpy.floating) -> typing.Any:
@typing.overload
def as_float(
    x: collections.abc.Sequence[typing.SupportsIndex]
    | numpy.typing.NDArray[numpy.integer],
    resolution: int | float | numpy.integer | numpy.floating,
) -> numpy.typing.NDArray[numpy.float64]: ...
@typing.overload
def as_float(
    x: int | numpy.integer, resolution: int | float | numpy.integer | numpy.floating
) -> float: ...
@typing.overload
def as_float(
    x: typing.SupportsIndex, resolution: int | float | numpy.integer | numpy.floating
) -> float | numpy.typing.NDArray[numpy.float64]: ...

# reason: pairs of node sets, each pair a length-2 sequence (sources, destinations), an integer array included (so pairs may be an integer array of shape (k, 2), or (k, 2, m) for sets of m nodes); options is EvaluationOptions; per_pair_eligibility has one Eligibility or None per pair; workers is an integer or None; one PlacementResult per pair
# generated: def evaluate_pairs(graph: Graph, pairs: typing.Any, options: typing.Any, per_pair_eligibility: typing.Any = None, workers: typing.Any = None) -> list:
def evaluate_pairs(
    graph: Graph,
    pairs: collections.abc.Sequence[
        collections.abc.Sequence[
            typing.SupportsIndex
            | collections.abc.Sequence[typing.SupportsIndex]
            | numpy.typing.NDArray[numpy.integer]
        ]
        | numpy.typing.NDArray[numpy.integer]
    ]
    | numpy.typing.NDArray[numpy.integer],
    options: EvaluationOptions,
    per_pair_eligibility: collections.abc.Sequence[Eligibility | None] | None = None,
    workers: typing.SupportsIndex | None = None,
) -> list[PlacementResult]: ...

# reason: ledger is a FlowLedger; actions a sequence of Place, Replace and Forward; mode a BatchMode; outputs may be None (TypeError otherwise)
# generated: def apply_actions(ledger: typing.Any, actions: typing.Any, *, mode: typing.Any = ..., outputs: typing.Any = None) -> list[PlacementResult]:
def apply_actions(
    ledger: FlowLedger,
    actions: collections.abc.Sequence[Place | Replace | Forward],
    *,
    mode: BatchMode = ...,
    outputs: OutputOptions | None = None,
) -> list[PlacementResult]: ...

# reason: sources and destinations are node sets; options is EvaluationOptions (TypeError otherwise)
# generated: def evaluate_demand(graph: Graph, sources: typing.Any, destinations: typing.Any, options: typing.Any) -> typing.Any:
def evaluate_demand(
    graph: Graph,
    sources: typing.SupportsIndex
    | collections.abc.Sequence[typing.SupportsIndex]
    | numpy.typing.NDArray[numpy.integer],
    destinations: typing.SupportsIndex
    | collections.abc.Sequence[typing.SupportsIndex]
    | numpy.typing.NDArray[numpy.integer],
    options: EvaluationOptions,
) -> PlacementResult: ...

# reason: source and destination are integers; options is KspOptions; a list of (cost, Path)
# generated: def ksp(graph: Graph, source: typing.Any, destination: typing.Any, options: typing.Any) -> list:
def ksp(
    graph: Graph,
    source: typing.SupportsIndex,
    destination: typing.SupportsIndex,
    options: KspOptions,
) -> list[tuple[int, Path]]: ...

# reason: source, destination and max_paths are integers (max_paths may be None); a list of Path
# generated: def resolve_to_paths(graph: Graph, dag: PredecessorDAG, source: typing.Any, destination: typing.Any, split_parallel_edges: bool = False, max_paths: typing.Any = None) -> list:
def resolve_to_paths(
    graph: Graph,
    dag: PredecessorDAG,
    source: typing.SupportsIndex,
    destination: typing.SupportsIndex,
    split_parallel_edges: bool = False,
    max_paths: typing.SupportsIndex | None = None,
) -> list[Path]: ...

# reason: sources and destinations are node sets; options is EvaluationOptions; workers is an integer or None
# generated: def evaluate_edge_failures(graph: Graph, sources: typing.Any, destinations: typing.Any, options: typing.Any, workers: typing.Any = None) -> typing.Any:
def evaluate_edge_failures(
    graph: Graph,
    sources: typing.SupportsIndex
    | collections.abc.Sequence[typing.SupportsIndex]
    | numpy.typing.NDArray[numpy.integer],
    destinations: typing.SupportsIndex
    | collections.abc.Sequence[typing.SupportsIndex]
    | numpy.typing.NDArray[numpy.integer],
    options: EvaluationOptions,
    workers: typing.SupportsIndex | None = None,
) -> SensitivityResult: ...

# reason: roots is a node set; options is SpfOptions (TypeError otherwise)
# generated: def spf(graph: Graph, roots: typing.Any, options: typing.Any) -> typing.Any:
def spf(
    graph: Graph,
    roots: typing.SupportsIndex
    | collections.abc.Sequence[typing.SupportsIndex]
    | numpy.typing.NDArray[numpy.integer],
    options: SpfOptions,
) -> SpfResult: ...

# reason: scalars give int, Python sequences int64 arrays; unknown-rank arrays may be 0-d (int) or 1-d (array)
# generated: def to_units(x: typing.Any, resolution: int | float | numpy.integer | numpy.floating) -> typing.Any:
@typing.overload
def to_units(
    x: int | float | numpy.floating | numpy.integer,
    resolution: int | float | numpy.integer | numpy.floating,
) -> int: ...
@typing.overload
def to_units(
    x: collections.abc.Sequence[int | float | numpy.floating | numpy.integer],
    resolution: int | float | numpy.integer | numpy.floating,
) -> numpy.typing.NDArray[numpy.int64]: ...
@typing.overload
def to_units(
    x: numpy.typing.NDArray[numpy.floating | numpy.integer],
    resolution: int | float | numpy.integer | numpy.floating,
) -> int | numpy.typing.NDArray[numpy.int64]: ...
