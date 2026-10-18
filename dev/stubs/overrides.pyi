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
# - a number: float | numpy.floating | typing.SupportsIndex.
# Arrays that the module returns carry their exact dtype; their shape (one-dimensional) and
# whether they are read-only views are in the docstrings.
#
# The overrides name the classes of the generated stub, which this file does not declare.
# ruff: noqa: F821

import collections.abc
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
    # generated: def __init__(self, num_nodes: typing.Any, src: typing.Any, dst: typing.Any, capacity_units: typing.Any, cost: typing.Any, resolution: typing.SupportsFloat | typing.SupportsIndex) -> None:
    def __init__(
        self,
        num_nodes: typing.SupportsIndex,
        src: collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer],
        dst: collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer],
        capacity_units: collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer],
        cost: collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer],
        resolution: typing.SupportsFloat | typing.SupportsIndex,
    ) -> None: ...
    # reason: a read-only int64 view
    # generated: def capacity_units(self) -> numpy.ndarray:
    def capacity_units(self) -> numpy.typing.NDArray[numpy.int64]: ...
    # reason: a read-only int64 view
    # generated: def cost(self) -> numpy.ndarray:
    def cost(self) -> numpy.typing.NDArray[numpy.int64]: ...
    # reason: a read-only int32 view
    # generated: def edge_dst(self) -> numpy.ndarray:
    def edge_dst(self) -> numpy.typing.NDArray[numpy.int32]: ...
    # reason: a read-only int32 view
    # generated: def edge_src(self) -> numpy.ndarray:
    def edge_src(self) -> numpy.typing.NDArray[numpy.int32]: ...
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
    # generated: def in_adj_edge_index(self) -> numpy.ndarray:
    def in_adj_edge_index(self) -> numpy.typing.NDArray[numpy.int32]: ...
    # reason: a read-only int32 view
    # generated: def in_row_offsets(self) -> numpy.ndarray:
    def in_row_offsets(self) -> numpy.typing.NDArray[numpy.int32]: ...
    # reason: a read-only int32 view
    # generated: def input_index(self) -> numpy.ndarray:
    def input_index(self) -> numpy.typing.NDArray[numpy.int32]: ...
    # reason: a read-only int32 view
    # generated: def row_offsets(self) -> numpy.ndarray:
    def row_offsets(self) -> numpy.typing.NDArray[numpy.int32]: ...

class PredDAG:
    # reason: no constructor: built by PredDAG.from_csr, Path.to_dag and the searches
    # generated: (none)
    def __init__(self, _no_public_constructor: typing.Never, /) -> None: ...
    # reason: offsets, parents and via are integers
    # generated: def from_csr(graph: Graph, offsets: typing.Any, parents: typing.Any, via: typing.Any) -> PredDAG:
    @staticmethod
    def from_csr(
        graph: Graph,
        offsets: collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer],
        parents: collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer],
        via: collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer],
    ) -> PredDAG: ...
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
    # generated: def drops(self) -> list:
    @property
    def drops(self) -> list[tuple[int, int]]: ...
    # reason: a list of (EdgeId, units) pairs
    # generated: def edges(self) -> list:
    @property
    def edges(self) -> list[tuple[int, int]]: ...
    # reason: an exact Python int for a Lossy flow, None otherwise
    # generated: def offered(self) -> typing.Any:
    @property
    def offered(self) -> int | None: ...
    # reason: an int, None when the route is empty
    # generated: def route_cost(self) -> typing.Any:
    @property
    def route_cost(self) -> int | None: ...

class FlowLedger:
    # reason: the baseline residual is integers, or None for the capacities
    # generated: def __init__(self, graph: Graph, residual: typing.Any = None) -> None:
    def __init__(
        self,
        graph: Graph,
        residual: collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer]
        | None = None,
    ) -> None: ...
    # reason: sources is a node set; elig may be None (all allowed)
    # generated: def bottleneck(self, sources: typing.Any, elig: Eligibility = None) -> list[int]:
    def bottleneck(
        self,
        sources: typing.SupportsIndex
        | collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer],
        elig: Eligibility | None = None,
    ) -> list[int]: ...
    # reason: a read-only int64 view
    # generated: def capacity_units(self) -> numpy.ndarray:
    def capacity_units(self) -> numpy.typing.NDArray[numpy.int64]: ...
    # reason: a FlowId is a tuple (ledger uid, id) of integers
    # generated: def record(self, flow_id: typing.Any) -> FlowRecord:
    def record(
        self, flow_id: tuple[typing.SupportsIndex, typing.SupportsIndex]
    ) -> FlowRecord: ...
    # reason: the new baseline is integers, or None for the capacities
    # generated: def reset(self, residual: typing.Any = None) -> None:
    def reset(
        self,
        residual: collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer]
        | None = None,
    ) -> None: ...

class Absolute:
    # reason: cost is an integer
    # generated: def __init__(self, cost: typing.Any) -> None:
    def __init__(self, cost: typing.SupportsIndex) -> None: ...

class Stretch:
    # reason: factor is a number (a float, a NumPy floating scalar or an integer)
    # generated: def __init__(self, factor: typing.Any) -> None:
    def __init__(
        self, factor: float | numpy.floating | typing.SupportsIndex
    ) -> None: ...

class DiminishingReturns:
    # reason: window is an integer, fraction a number
    # generated: def __init__(self, window: typing.Any, fraction: typing.Any) -> None:
    def __init__(
        self,
        window: typing.SupportsIndex,
        fraction: float | numpy.floating | typing.SupportsIndex,
    ) -> None: ...

class CapacityAware:
    # reason: the counts and the bandwidth are integers (max_flow_count may be None); the ceiling is Absolute, Stretch or None; diminishing is DiminishingReturns or None
    # generated: def __init__(self, *, flow_count: typing.Any = 1, max_flow_count: typing.Any = None, max_path_cost: typing.Any = None, min_flow_bandwidth: typing.Any = 0, reoptimize: bool = False, diminishing: typing.Any = None) -> None:
    def __init__(
        self,
        *,
        flow_count: typing.SupportsIndex = 1,
        max_flow_count: typing.SupportsIndex | None = None,
        max_path_cost: Absolute | Stretch | None = None,
        min_flow_bandwidth: typing.SupportsIndex = 0,
        reoptimize: bool = False,
        diminishing: DiminishingReturns | None = None,
    ) -> None: ...
    # reason: the ceiling given at construction: Absolute, Stretch or None
    # generated: def max_path_cost(self) -> typing.Any:
    @property
    def max_path_cost(self) -> Absolute | Stretch | None: ...

class Tunnel:
    # reason: no constructor: a tunnel of ForwardingProgram.tunnels
    # generated: (none)
    def __init__(self, _no_public_constructor: typing.Never, /) -> None: ...
    # reason: the path's EdgeIds in path order
    # generated: def edges(self) -> tuple:
    @property
    def edges(self) -> tuple[int, ...]: ...

class Route:
    # reason: head and weight are integers (weight checked by from_routes); bundle a Path or a PredDAG (TypeError otherwise)
    # generated: def __init__(self, head: typing.Any, bundle: typing.Any, *, weight: typing.Any = 1) -> None:
    def __init__(
        self,
        head: typing.SupportsIndex,
        bundle: Path | PredDAG,
        *,
        weight: typing.SupportsIndex = 1,
    ) -> None: ...

class ForwardingProgram:
    # reason: no constructor: made by from_routes and returned by the policies
    # generated: (none)
    def __init__(self, _no_public_constructor: typing.Never, /) -> None: ...
    # reason: sources and destinations are node sets; routes a sequence of Route; edge_weights and residual integers or None
    # generated: def from_routes(graph: Graph, sources: typing.Any, destinations: typing.Any, routes: typing.Any, *, placement: FlowPlacement | None = None, weights: WeightBasis | None = None, edge_weights: typing.Any = None, residual: typing.Any = None) -> ForwardingProgram:
    @staticmethod
    def from_routes(
        graph: Graph,
        sources: typing.SupportsIndex
        | collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer],
        destinations: typing.SupportsIndex
        | collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer],
        routes: collections.abc.Sequence[Route],
        *,
        placement: FlowPlacement | None = None,
        weights: WeightBasis | None = None,
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

class Global:
    # reason: support is a PredDAG or None
    # generated: def __init__(self, *, support: typing.Any = None) -> None:
    def __init__(self, *, support: PredDAG | None = None) -> None: ...

class Pinned:
    # reason: black_holes may be None
    # generated: def __init__(self, program: ForwardingProgram, *, black_holes: Eligibility = None) -> None:
    def __init__(
        self, program: ForwardingProgram, *, black_holes: Eligibility | None = None
    ) -> None: ...

class FlowPolicyConfig:
    # reason: routing is one of the four routing modes (TypeError otherwise)
    # generated: def __init__(self, *, routing: typing.Any, placement: FlowPlacement | None = None, admission: Admission = ..., selection: PathSelection | None = None, origin: Origin | None = None, weights: WeightBasis | None = None, edge_weights: typing.Any = None) -> None:
    def __init__(
        self,
        *,
        routing: CostOnly | CapacityAware | Global | Pinned,
        placement: FlowPlacement | None = None,
        admission: Admission = ...,
        selection: PathSelection | None = None,
        origin: Origin | None = None,
        weights: WeightBasis | None = None,
        edge_weights: collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer]
        | None = None,
    ) -> None: ...
    # reason: the explicit weights, one int per edge, or None
    # generated: def edge_weights(self) -> typing.Any:
    @property
    def edge_weights(self) -> tuple[int, ...] | None: ...
    # reason: the routing mode given at construction
    # generated: def routing(self) -> typing.Any:
    @property
    def routing(self) -> CostOnly | CapacityAware | Global | Pinned: ...

class SourceFlow:
    # reason: no constructor: one row of PlacementResult.by_source
    # generated: (none)
    def __init__(self, _no_public_constructor: typing.Never, /) -> None: ...
    # reason: an exact int (a 128-bit counter)
    # generated: def admitted(self) -> typing.Any:
    @property
    def admitted(self) -> int: ...
    # reason: an exact int (a 128-bit counter)
    # generated: def dropped(self) -> typing.Any:
    @property
    def dropped(self) -> int: ...
    # reason: an exact int when LOSSY, None when LOSSLESS
    # generated: def offered(self) -> typing.Any:
    @property
    def offered(self) -> int | None: ...
    # reason: an exact int (a 128-bit counter)
    # generated: def source_loss(self) -> typing.Any:
    @property
    def source_loss(self) -> int: ...

class TunnelFlow:
    # reason: no constructor: one row of PlacementResult.by_tunnel
    # generated: (none)
    def __init__(self, _no_public_constructor: typing.Never, /) -> None: ...
    # reason: an exact int (a 128-bit counter)
    # generated: def assigned(self) -> typing.Any:
    @property
    def assigned(self) -> int: ...
    # reason: an exact int (a 128-bit counter)
    # generated: def dropped(self) -> typing.Any:
    @property
    def dropped(self) -> int: ...

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
    # generated: def call_drops(self) -> list:
    @property
    def call_drops(self) -> list[tuple[int, int]]: ...
    # reason: an exact int when LOSSY, None when LOSSLESS
    # generated: def offered(self) -> typing.Any:
    @property
    def offered(self) -> int | None: ...
    # reason: an exact int, or float('inf') after an infinite request
    # generated: def requested(self) -> typing.Any:
    @property
    def requested(self) -> int | float: ...
    # reason: an exact int (a 128-bit counter)
    # generated: def source_loss(self) -> typing.Any:
    @property
    def source_loss(self) -> int: ...
    # reason: an exact int for the controllers, None otherwise
    # generated: def total_cost(self) -> typing.Any:
    @property
    def total_cost(self) -> int | None: ...
    # reason: an exact int, None after an unlimited LOSSLESS request
    # generated: def unadmitted(self) -> typing.Any:
    @property
    def unadmitted(self) -> int | None: ...
    # reason: a list of (EdgeId, units) pairs; units are exact Python ints
    # generated: def drops(self) -> list:
    @property
    def drops(self) -> list[tuple[int, int]]: ...
    # reason: an int64 array when requested, else None
    # generated: def edge_flows(self) -> typing.Any:
    @property
    def edge_flows(self) -> numpy.typing.NDArray[numpy.int64] | None: ...
    # reason: an int, or float('inf') after an infinite request
    # generated: def remaining(self) -> typing.Any:
    @property
    def remaining(self) -> int | float: ...
    # reason: an int64 array when requested, else None
    # generated: def residual(self) -> typing.Any:
    @property
    def residual(self) -> numpy.typing.NDArray[numpy.int64] | None: ...

class DagCacheStats:
    # reason: no constructor: returned by DagCache.stats
    # generated: (none)
    def __init__(self, _no_public_constructor: typing.Never, /) -> None: ...

class FlowPolicy:
    # reason: sources and destinations are node sets; elig and cache may be None
    # generated: def __init__(self, ledger: FlowLedger, sources: typing.Any, destinations: typing.Any, config: FlowPolicyConfig, elig: Eligibility = None, cache: DagCache = None) -> None:
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
        elig: Eligibility | None = None,
        cache: DagCache | None = None,
    ) -> None: ...
    # reason: volume is integer units or float('inf') (a NumPy floating +inf too); outputs may be None; the result is a PlacementResult
    # generated: def place(self, volume: typing.Any, outputs: typing.Any = None) -> typing.Any:
    def place(
        self,
        volume: typing.SupportsIndex | float | numpy.floating,
        outputs: OutputOptions | None = None,
    ) -> PlacementResult: ...
    # reason: offer is integer units or float('inf'); outputs may be None; the result is a PlacementResult
    # generated: def evaluate(self, offer: typing.Any, *, outputs: typing.Any = None) -> typing.Any:
    def evaluate(
        self,
        offer: typing.SupportsIndex | float | numpy.floating,
        *,
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
    # reason: offer is integer units or float('inf'); outputs may be None; the result is a PlacementResult
    # generated: def replace(self, offer: typing.Any, *, outputs: typing.Any = None) -> typing.Any:
    def replace(
        self,
        offer: typing.SupportsIndex | float | numpy.floating,
        *,
        outputs: OutputOptions | None = None,
    ) -> PlacementResult: ...
    # reason: a FlowId is a tuple (ledger uid, id) of integers
    # generated: def record(self, flow_id: typing.Any) -> FlowRecord:
    def record(
        self, flow_id: tuple[typing.SupportsIndex, typing.SupportsIndex]
    ) -> FlowRecord: ...

class QueryOptions:
    # reason: elig, cache, outputs and offered may be None; offered is integer units or float('inf')
    # generated: def __init__(self, config: FlowPolicyConfig, *, elig: Eligibility = None, cache: DagCache = None, outputs: typing.Any = None, offered: typing.Any = None) -> None:
    def __init__(
        self,
        config: FlowPolicyConfig,
        *,
        elig: Eligibility | None = None,
        cache: DagCache | None = None,
        outputs: OutputOptions | None = None,
        offered: typing.SupportsIndex | float | numpy.floating | None = None,
    ) -> None: ...
    # reason: the offer given, an int or float('inf'), or None
    # generated: def offered(self) -> typing.Any:
    @property
    def offered(self) -> int | float | None: ...
    # reason: None when no cache was given
    # generated: def cache(self) -> DagCache:
    @property
    def cache(self) -> DagCache | None: ...
    # reason: None when no Eligibility was given
    # generated: def elig(self) -> Eligibility:
    @property
    def elig(self) -> Eligibility | None: ...

class SensitivityEntry:
    # reason: no constructor: the entries of a SensitivityResult
    # generated: (none)
    def __init__(self, _no_public_constructor: typing.Never, /) -> None: ...

class SensitivityResult:
    # reason: no constructor: returned by sensitivity
    # generated: (none)
    def __init__(self, _no_public_constructor: typing.Never, /) -> None: ...

class KspOptions:
    # reason: k is an integer; the ceiling is Absolute, Stretch or None; residual is integers or None; elig may be None (each validated by ksp)
    # generated: def __init__(self, k: typing.Any, selection: PathSelection, *, max_path_cost: typing.Any = None, residual: typing.Any = None, elig: Eligibility = None) -> None:
    def __init__(
        self,
        k: typing.SupportsIndex,
        selection: PathSelection,
        *,
        max_path_cost: Absolute | Stretch | None = None,
        residual: collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer]
        | None = None,
        elig: Eligibility | None = None,
    ) -> None: ...
    # reason: None when no Eligibility was given
    # generated: def elig(self) -> Eligibility:
    @property
    def elig(self) -> Eligibility | None: ...
    # reason: the k given at construction (an integer)
    # generated: def k(self) -> typing.Any:
    @property
    def k(self) -> typing.SupportsIndex: ...
    # reason: the ceiling given at construction
    # generated: def max_path_cost(self) -> typing.Any:
    @property
    def max_path_cost(self) -> Absolute | Stretch | None: ...
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
    # reason: targets is a node set or None; residual is integers or None; min_residual is an integer; elig may be None (each validated by spf)
    # generated: def __init__(self, direction: Direction, selection: PathSelection, *, targets: typing.Any = None, residual: typing.Any = None, min_residual: typing.Any = 1, elig: Eligibility = None) -> None:
    def __init__(
        self,
        direction: Direction,
        selection: PathSelection,
        *,
        targets: typing.SupportsIndex
        | collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer]
        | None = None,
        residual: collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer]
        | None = None,
        min_residual: typing.SupportsIndex = 1,
        elig: Eligibility | None = None,
    ) -> None: ...
    # reason: None when no Eligibility was given
    # generated: def elig(self) -> Eligibility:
    @property
    def elig(self) -> Eligibility | None: ...
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
    # reason: always a PredDAG
    # generated: def dag(self) -> typing.Any:
    @property
    def dag(self) -> PredDAG: ...
    # reason: the nearest targets' cost, None when no target is reached and for a full search
    # generated: def target_cost(self) -> typing.Any:
    @property
    def target_cost(self) -> int | None: ...

# reason: one integer gives a float, integers give a float64 array (the array overload first: a NumPy array also has __index__)
# generated: def as_float(x: typing.Any, resolution: typing.SupportsFloat | typing.SupportsIndex) -> typing.Any:
@typing.overload
def as_float(
    x: collections.abc.Sequence[typing.SupportsIndex]
    | numpy.typing.NDArray[numpy.integer],
    resolution: typing.SupportsFloat | typing.SupportsIndex,
) -> numpy.typing.NDArray[numpy.float64]: ...
@typing.overload
def as_float(
    x: typing.SupportsIndex, resolution: typing.SupportsFloat | typing.SupportsIndex
) -> float: ...

# reason: pairs of node sets, each pair a length-2 sequence (sources, destinations), an integer array included (so pairs may be an integer array of shape (k, 2), or (k, 2, m) for sets of m nodes); opts is QueryOptions; per_pair_elig has one Eligibility or None per pair; workers is an integer or None; one PlacementResult per pair
# generated: def batch_capacity(graph: Graph, pairs: typing.Any, opts: typing.Any, per_pair_elig: typing.Any = None, workers: typing.Any = None) -> typing.Any:
def batch_capacity(
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
    opts: QueryOptions,
    per_pair_elig: collections.abc.Sequence[Eligibility | None] | None = None,
    workers: typing.SupportsIndex | None = None,
) -> list[PlacementResult]: ...

# reason: sources and destinations are node sets; opts is QueryOptions (TypeError otherwise)
# generated: def capacity(graph: Graph, sources: typing.Any, destinations: typing.Any, opts: typing.Any) -> typing.Any:
def capacity(
    graph: Graph,
    sources: typing.SupportsIndex
    | collections.abc.Sequence[typing.SupportsIndex]
    | numpy.typing.NDArray[numpy.integer],
    destinations: typing.SupportsIndex
    | collections.abc.Sequence[typing.SupportsIndex]
    | numpy.typing.NDArray[numpy.integer],
    opts: QueryOptions,
) -> PlacementResult: ...

# reason: src and dst are integers; opts is KspOptions; a list of (cost, Path)
# generated: def ksp(graph: Graph, src: typing.Any, dst: typing.Any, opts: typing.Any) -> list:
def ksp(
    graph: Graph, src: typing.SupportsIndex, dst: typing.SupportsIndex, opts: KspOptions
) -> list[tuple[int, Path]]: ...

# reason: src, dst and max_paths are integers (max_paths may be None); a list of Path
# generated: def resolve_to_paths(graph: Graph, dag: PredDAG, src: typing.Any, dst: typing.Any, split_parallel_edges: bool = False, max_paths: typing.Any = None) -> list:
def resolve_to_paths(
    graph: Graph,
    dag: PredDAG,
    src: typing.SupportsIndex,
    dst: typing.SupportsIndex,
    split_parallel_edges: bool = False,
    max_paths: typing.SupportsIndex | None = None,
) -> list[Path]: ...

# reason: sources and destinations are node sets; opts is QueryOptions; workers is an integer or None
# generated: def sensitivity(graph: Graph, sources: typing.Any, destinations: typing.Any, opts: typing.Any, workers: typing.Any = None) -> typing.Any:
def sensitivity(
    graph: Graph,
    sources: typing.SupportsIndex
    | collections.abc.Sequence[typing.SupportsIndex]
    | numpy.typing.NDArray[numpy.integer],
    destinations: typing.SupportsIndex
    | collections.abc.Sequence[typing.SupportsIndex]
    | numpy.typing.NDArray[numpy.integer],
    opts: QueryOptions,
    workers: typing.SupportsIndex | None = None,
) -> SensitivityResult: ...

# reason: roots is a node set; opts is SpfOptions (TypeError otherwise)
# generated: def spf(graph: Graph, roots: typing.Any, opts: typing.Any) -> typing.Any:
def spf(
    graph: Graph,
    roots: typing.SupportsIndex
    | collections.abc.Sequence[typing.SupportsIndex]
    | numpy.typing.NDArray[numpy.integer],
    opts: SpfOptions,
) -> SpfResult: ...

# reason: a scalar number gives an int, a sequence or a 1-D array of numbers an int64 array (a 0-d array gives an int, which the array overload cannot say)
# generated: def to_units(x: typing.Any, resolution: typing.SupportsFloat | typing.SupportsIndex) -> typing.Any:
@typing.overload
def to_units(
    x: float | numpy.floating | numpy.integer,
    resolution: typing.SupportsFloat | typing.SupportsIndex,
) -> int: ...
@typing.overload
def to_units(
    x: collections.abc.Sequence[float | numpy.floating | numpy.integer]
    | numpy.typing.NDArray[numpy.floating | numpy.integer],
    resolution: typing.SupportsFloat | typing.SupportsIndex,
) -> numpy.typing.NDArray[numpy.int64]: ...
