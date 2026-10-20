# Generated from the compiled module; do not edit.
"""NetGraph-Core: graph, flow and placement models on exact integer units.

Every public name is re-exported from the compiled ``_netgraph_core`` module, so users
import everything from here. Type checkers read the generated stub ``__init__.pyi`` next
to this file (``py.typed`` marks the package as typed).

- The graph on integer units: ``Graph``, ``Eligibility``, ``to_units``, ``as_float``.
- Route structures: ``PredecessorDAG``, ``Path``.
- Route choice and placement settings: ``RouteTieBreak``, ``RouteAllocation``, ``Admission``,
  ``SourceAllocation``, ``WeightBasis``.
- The ledger of placed traffic: ``FlowLedger`` and its ``FlowRecord`` snapshots.
- Forwarding programs: ``ForwardingProgram`` (with its ``Tunnel`` entries and
  ``ForwardingProgram.from_routes`` over ``RouteDefinition`` values), ``TunnelStatus``.
- The policy engine: ``FlowPolicy``, ``FlowPolicyConfig``, the routings ``ShortestPaths``,
  ``DistributedTE`` (T1, with ``TunnelCount``), ``CentralizedTE`` (with ``tunnel_count``: T2), ``Pinned``, and
  ``PlacementResult`` (with its ``SourceFlow`` and ``TunnelFlow`` rows and its ``TunnelPlan``),
  ``OutputOptions``, ``Termination``.
- Several policies in one operation: ``apply_actions`` over ``Place``, ``Replace`` and ``Forward``
  actions, in a ``BatchMode`` (ORDERED or SIMULTANEOUS).
- The shortest-path cache: ``ShortestPathCache`` (building and frozen modes), ``ShortestPathCacheMiss``.
- The evaluators: ``evaluate_demand``, ``evaluate_pairs``, ``evaluate_edge_failures``.
- The searches: ``spf`` (``Direction``, ``SpfOptions``, ``SpfResult``, ``UNREACHABLE``),
  ``ksp`` (``KspOptions``), ``resolve_to_paths``.

``_netgraph_core`` also holds a private entry of the flow core, never re-exported here."""

from __future__ import annotations

import collections.abc
import fractions
import typing

import numpy
import numpy.typing

__all__: list[str] = [
    "__version__",
    "Graph",
    "Eligibility",
    "to_units",
    "as_float",
    "PredecessorDAG",
    "Path",
    "RouteTieBreak",
    "RouteAllocation",
    "Admission",
    "SourceAllocation",
    "WeightBasis",
    "FlowLedger",
    "FlowRecord",
    "TunnelStatus",
    "Tunnel",
    "RouteDefinition",
    "ForwardingProgram",
    "AbsoluteCost",
    "CostStretch",
    "TunnelCount",
    "ShortestPaths",
    "DistributedTE",
    "CentralizedTE",
    "Pinned",
    "FlowPolicyConfig",
    "OutputOptions",
    "Termination",
    "PlacementResult",
    "SourceFlow",
    "TunnelFlow",
    "TunnelPlan",
    "FlowPolicy",
    "Place",
    "Replace",
    "Forward",
    "BatchMode",
    "apply_actions",
    "ShortestPathCache",
    "ShortestPathCacheMiss",
    "ShortestPathCacheStats",
    "EvaluationOptions",
    "SensitivityEntry",
    "SensitivityResult",
    "evaluate_demand",
    "evaluate_pairs",
    "evaluate_edge_failures",
    "Direction",
    "SpfOptions",
    "SpfResult",
    "UNREACHABLE",
    "spf",
    "KspOptions",
    "ksp",
    "resolve_to_paths",
]
__version__: str

class RouteAllocation:
    """
    How traffic spreads over the routes (ShortestPaths's kernel; the interior split of supplied PredecessorDAG
    routes, EQUAL or WEIGHTED).
    OPTIMIZED: the minimum-cost flow core over the routes (a controller restricted to them), re-solved
    over the free capacity plus the policy's own allocation (it may move its own earlier placement,
    never another's). WEIGHTED (WCMP, weighted-cost multipath): every node splits in proportion to
    weights fixed when the routes are prepared (highest averages (D'Hondt): assign each next unit
    to the largest weight / (assigned + 1), ties to the lower EdgeId; the weights come from
    WeightBasis). EQUAL (ECMP, equal-cost
    multipath): every node splits equally, by round robin.

    Members:

      EQUAL

      WEIGHTED

      OPTIMIZED
    """

    EQUAL: typing.ClassVar[RouteAllocation]
    OPTIMIZED: typing.ClassVar[RouteAllocation]
    WEIGHTED: typing.ClassVar[RouteAllocation]
    __members__: typing.ClassVar[dict[str, RouteAllocation]]

    @typing.overload
    def __eq__(self, other: RouteAllocation, /) -> bool: ...
    @typing.overload
    def __eq__(self, other: typing.Any, /) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __init__(self, value: typing.SupportsInt | typing.SupportsIndex) -> None: ...
    def __int__(self) -> int: ...
    @typing.overload
    def __ne__(self, other: RouteAllocation, /) -> bool: ...
    @typing.overload
    def __ne__(self, other: typing.Any, /) -> bool: ...
    def __repr__(self) -> str: ...
    def __setstate__(
        self, state: typing.SupportsInt | typing.SupportsIndex, /
    ) -> None: ...
    def __str__(self) -> str: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class Admission:
    """
    LOSSLESS admits the largest total the routes carry with no loss at all; LOSSY sends the whole
    finite offer hop by hop and records where it is lost (every model except the optimizing ones:
    ShortestPaths OPTIMIZED and CentralizedTE without tunnel_count).

    Members:

      LOSSLESS

      LOSSY
    """

    LOSSLESS: typing.ClassVar[Admission]
    LOSSY: typing.ClassVar[Admission]
    __members__: typing.ClassVar[dict[str, Admission]]

    @typing.overload
    def __eq__(self, other: Admission, /) -> bool: ...
    @typing.overload
    def __eq__(self, other: typing.Any, /) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __init__(self, value: typing.SupportsInt | typing.SupportsIndex) -> None: ...
    def __int__(self) -> int: ...
    @typing.overload
    def __ne__(self, other: Admission, /) -> bool: ...
    @typing.overload
    def __ne__(self, other: typing.Any, /) -> bool: ...
    def __repr__(self) -> str: ...
    def __setstate__(
        self, state: typing.SupportsInt | typing.SupportsIndex, /
    ) -> None: ...
    def __str__(self) -> str: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class SourceAllocation:
    """
    How a demand's traffic enters at its sources, separately from how it is routed. EQUAL: round
    robin over the roster (ascending NodeId, the first sources one unit more); a source whose share
    cannot go anywhere keeps it and loses it (LOSSLESS then admits only totals that give it
    nothing). OPTIMIZED: the controller chooses the injections (one virtual source).

    Members:

      EQUAL

      OPTIMIZED
    """

    EQUAL: typing.ClassVar[SourceAllocation]
    OPTIMIZED: typing.ClassVar[SourceAllocation]
    __members__: typing.ClassVar[dict[str, SourceAllocation]]

    @typing.overload
    def __eq__(self, other: SourceAllocation, /) -> bool: ...
    @typing.overload
    def __eq__(self, other: typing.Any, /) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __init__(self, value: typing.SupportsInt | typing.SupportsIndex) -> None: ...
    def __int__(self) -> int: ...
    @typing.overload
    def __ne__(self, other: SourceAllocation, /) -> bool: ...
    @typing.overload
    def __ne__(self, other: typing.Any, /) -> bool: ...
    def __repr__(self) -> str: ...
    def __setstate__(
        self, state: typing.SupportsInt | typing.SupportsIndex, /
    ) -> None: ...
    def __str__(self) -> str: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class WeightBasis:
    """
    Where WEIGHTED member weights come from (ShortestPaths WEIGHTED and supplied PredecessorDAG routes).
    CAPACITY (the default): each member's static capacity. EXPLICIT: edge_weights, one integer per
    edge (a positive weight on an edge of capacity 0 becomes 0). RESIDUAL_SNAPSHOT: each member's free
    capacity read once when the routes are prepared, against one common view, and frozen for the
    epoch. A member of weight 0 is pruned with the branches it leaves without a destination, inside
    the route only. The tie rule never chooses weights.

    Members:

      CAPACITY

      EXPLICIT

      RESIDUAL_SNAPSHOT
    """

    CAPACITY: typing.ClassVar[WeightBasis]
    EXPLICIT: typing.ClassVar[WeightBasis]
    RESIDUAL_SNAPSHOT: typing.ClassVar[WeightBasis]
    __members__: typing.ClassVar[dict[str, WeightBasis]]

    @typing.overload
    def __eq__(self, other: WeightBasis, /) -> bool: ...
    @typing.overload
    def __eq__(self, other: typing.Any, /) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __init__(self, value: typing.SupportsInt | typing.SupportsIndex) -> None: ...
    def __int__(self) -> int: ...
    @typing.overload
    def __ne__(self, other: WeightBasis, /) -> bool: ...
    @typing.overload
    def __ne__(self, other: typing.Any, /) -> bool: ...
    def __repr__(self) -> str: ...
    def __setstate__(
        self, state: typing.SupportsInt | typing.SupportsIndex, /
    ) -> None: ...
    def __str__(self) -> str: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class RouteTieBreak:
    """
    How a route choice breaks a cost tie; cost always comes first. FEWER_HOPS: then fewer hops,
    then the smaller EdgeId sequence. HIGHER_CAPACITY / HIGHER_AVAILABLE_CAPACITY /
    HIGHER_HEADROOM: then the larger path bottleneck of static capacity / availability /
    availability over static capacity (exact; 0 at static capacity 0), then fewer hops, then the
    EdgeId sequence. Availability: unreserved bandwidth for DistributedTE, the service view for
    ShortestPaths, the residual for spf and ksp, static capacity without one. The same order
    keeps one of each (u, v) pair's parallel edges when parallel_edges is False. It never
    supplies weights and never enables a controller objective.

    Members:

      FEWER_HOPS

      HIGHER_CAPACITY

      HIGHER_AVAILABLE_CAPACITY

      HIGHER_HEADROOM
    """

    FEWER_HOPS: typing.ClassVar[RouteTieBreak]
    HIGHER_AVAILABLE_CAPACITY: typing.ClassVar[RouteTieBreak]
    HIGHER_CAPACITY: typing.ClassVar[RouteTieBreak]
    HIGHER_HEADROOM: typing.ClassVar[RouteTieBreak]
    __members__: typing.ClassVar[dict[str, RouteTieBreak]]

    @typing.overload
    def __eq__(self, other: RouteTieBreak, /) -> bool: ...
    @typing.overload
    def __eq__(self, other: typing.Any, /) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __init__(self, value: typing.SupportsInt | typing.SupportsIndex) -> None: ...
    def __int__(self) -> int: ...
    @typing.overload
    def __ne__(self, other: RouteTieBreak, /) -> bool: ...
    @typing.overload
    def __ne__(self, other: typing.Any, /) -> bool: ...
    def __repr__(self) -> str: ...
    def __setstate__(
        self, state: typing.SupportsInt | typing.SupportsIndex, /
    ) -> None: ...
    def __str__(self) -> str: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class Direction:
    """
    The orientation of a search: FORWARD from the roots along the edges, or REVERSE toward the roots
    against the edges (the DAG it returns keeps the edges' own direction, so its walks lead to the
    roots).

    Members:

      FORWARD

      REVERSE
    """

    FORWARD: typing.ClassVar[Direction]
    REVERSE: typing.ClassVar[Direction]
    __members__: typing.ClassVar[dict[str, Direction]]

    @typing.overload
    def __eq__(self, other: Direction, /) -> bool: ...
    @typing.overload
    def __eq__(self, other: typing.Any, /) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __init__(self, value: typing.SupportsInt | typing.SupportsIndex) -> None: ...
    def __int__(self) -> int: ...
    @typing.overload
    def __ne__(self, other: Direction, /) -> bool: ...
    @typing.overload
    def __ne__(self, other: typing.Any, /) -> bool: ...
    def __repr__(self) -> str: ...
    def __setstate__(
        self, state: typing.SupportsInt | typing.SupportsIndex, /
    ) -> None: ...
    def __str__(self) -> str: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class Eligibility:
    """
    Which nodes and edges of one graph may be used: one flag per node and per edge, fixed when built
    (it does not depend on traffic). Built by Graph.eligibility(excluded_nodes, excluded_edges) and
    Eligibility.intersect.
    """

    def __init__(self, _no_public_constructor: typing.Never, /) -> None: ...
    @staticmethod
    def intersect(a: Eligibility, b: Eligibility) -> Eligibility:
        """
        The intersection of two values of one graph (ValueError for values of different graphs):
        an operand itself when the result equals it (a first, then b), otherwise a new value with a
        fresh uid; equal intersections built separately get distinct uids.
        """

    def edge_allowed(self) -> numpy.typing.NDArray[numpy.bool_]:
        """
        bool[E] read-only view (not excluded and both endpoints allowed)
        """

    def graph_uid(self) -> int:
        """
        uid of the graph this value belongs to
        """

    def node_allowed(self) -> numpy.typing.NDArray[numpy.bool_]:
        """
        bool[N] read-only view
        """

    def uid(self) -> int:
        """
        0 for the all-allowed value, a fresh uid otherwise
        """

class Graph:
    """
    Immutable directed multigraph on integer units (parallel edges and self-loops allowed). Edges are
    given in input order and numbered by sorting on (source, cost, destination), ties kept in input
    order: that number is the public EdgeId. input_positions()[e] is the input position of edge e.
    Per-edge values use EdgeId order; incoming_edges() instead lists EdgeIds in incoming-slot order.
    Costs are positive integers; capacities are int64 units of the explicit resolution q (a power of
    two in [2^-32, 2^32]); physical capacities enter through to_units.
    """

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
        resolution: typing.SupportsFloat | typing.SupportsIndex,
    ) -> None:
        """
        Validates, in order, and reports the first defect: the resolution (a power of two in
        [2^-32, 2^32]; a bool is a TypeError), the unit values (every capacity in [0, 2^62)), the array
        shapes, zero or negative costs, the remaining bounds (endpoints in range, then the capacity
        and cost sums below 2^62); a value of a wrong type (a bool or a float for an integer, a
        non-sequence, an array of another rank or dtype) raises TypeError, a value outside its domain
        ValueError. Integers are read without narrowing, so the error does not depend on how the
        input is represented.
        """

    def edge_capacity_units(self) -> numpy.typing.NDArray[numpy.int64]:
        """
        int64[E] capacities in units, id order (read-only view)
        """

    def edge_cost_gcd(self) -> int:
        """
        gcd of the edge costs (0 without edges)
        """

    def edge_costs(self) -> numpy.typing.NDArray[numpy.int64]:
        """
        int64[E] costs, id order (read-only view)
        """

    def edge_destinations(self) -> numpy.typing.NDArray[numpy.int32]:
        """
        int32[E] destination node per edge (read-only view)
        """

    def edge_sources(self) -> numpy.typing.NDArray[numpy.int32]:
        """
        int32[E] source node per edge (read-only view)
        """

    def eligibility(
        self,
        excluded_nodes: collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer] = (),
        excluded_edges: collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer] = (),
    ) -> Eligibility:
        """
        Static eligibility from excluded node ids and EdgeIds (duplicates and unsorted lists
        accepted, out-of-range ids rejected); an edge is allowed when it is not excluded and
        both endpoints are allowed. Empty lists mean all allowed (uid 0).
        """

    def incoming_edges(self) -> numpy.typing.NDArray[numpy.int32]:
        """
        int32[E]: the in-edges of v in (cost, source, input) order, as EdgeIds
        """

    def incoming_offsets(self) -> numpy.typing.NDArray[numpy.int32]:
        """
        int32[N + 1]: the in-edges of v are incoming_edges()[incoming_offsets()[v] : incoming_offsets()[v + 1]]
        """

    def input_positions(self) -> numpy.typing.NDArray[numpy.int32]:
        """
        int32[E]: EdgeId -> input position (read-only view); build the inverse yourself if needed
        """

    def max_edge_cost(self) -> int:
        """
        largest edge cost (0 without edges)
        """

    def num_edges(self) -> int:
        """
        E: EdgeIds are 0 .. E - 1
        """

    def num_nodes(self) -> int:
        """
        N: node ids are 0 .. N - 1
        """

    def outgoing_offsets(self) -> numpy.typing.NDArray[numpy.int32]:
        """
        int32[N + 1]: the out-edges of u are the ids [outgoing_offsets()[u], outgoing_offsets()[u + 1])
        """

    def resolution(self) -> float:
        """
        the resolution q (units of capacity)
        """

    def uid(self) -> int:
        """
        process-unique id (never reused)
        """

class PredecessorDAG:
    """
    Immutable route structure of one graph, a predecessor DAG (directed, no cycle): for each node v,
    the edges that reach v on the route and the node each comes from (its parent). The entries of v
    are parents()[k], via_edges()[k] for k in [parent_offsets()[v], parent_offsets()[v + 1]), in
    strictly ascending EdgeId, so two DAGs of the same graph with the same entries have equal arrays
    and compare equal; DAGs of different graphs never compare equal.
    parents() and via_edges() are read-only views that keep the DAG alive; parent_offsets() is a new
    read-only array on each call (the DAG stores only the nodes that have entries). Built by
    PredecessorDAG.from_csr, Path.to_dag(), the searches and the policies (FlowRecord.route).
    """

    def __init__(self, _no_public_constructor: typing.Never, /) -> None: ...
    @staticmethod
    def from_csr(
        graph: Graph,
        parent_offsets: collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer],
        parents: collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer],
        via_edges: collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer],
    ) -> PredecessorDAG:
        """
        Validates, in order, before any traversal (ValueError on the first defect): parent_offsets of
        length N + 1, parent_offsets[0] == 0, monotone, parent_offsets[N] == len(parents) ==
        len(via_edges); then every parent and via edge in range, every via edge joining (parent, child),
        no duplicate (child, via edge), acyclic; within a node the entries in ascending via EdgeId.
        """

    def __eq__(self, arg0: typing.Any) -> bool: ...
    def __hash__(self) -> int: ...
    def __repr__(self) -> str: ...
    def graph_uid(self) -> int:
        """
        uid of the graph this DAG belongs to
        """

    def parent_offsets(self) -> numpy.typing.NDArray[numpy.int32]:
        """
        int32[N + 1]: the entries of node v are [parent_offsets[v], parent_offsets[v + 1]); a new read-only array per call
        """

    def parents(self) -> numpy.typing.NDArray[numpy.int32]:
        """
        int32 parent NodeId per entry (read-only view)
        """

    def via_edges(self) -> numpy.typing.NDArray[numpy.int32]:
        """
        int32 via EdgeId per entry (read-only view)
        """

class Path:
    """
    Immutable route of one graph, hop by hop: each hop lists one or more parallel edges (members)
    from the same u to the same v, in the caller's order. cost = the sum over hops of the cost of
    each hop's cheapest member.
    """

    def __eq__(self, arg0: typing.Any) -> bool: ...
    def __hash__(self) -> int: ...
    def __init__(
        self,
        graph: Graph,
        hops: collections.abc.Sequence[
            collections.abc.Sequence[typing.SupportsIndex]
            | numpy.typing.NDArray[numpy.integer]
        ]
        | numpy.typing.NDArray[numpy.integer],
    ) -> None:
        """
        Validates, in order (ValueError on the first defect): every member id in range, no duplicate
        member, a non-empty hop list, every hop non-empty, the members of a hop sharing (u, v),
        contiguous hops, no repeated node.
        """

    def __repr__(self) -> str: ...
    def cost(self) -> int:
        """
        sum over hops of the minimum member cost
        """

    def destination(self) -> int:
        """
        the last node of the walk
        """

    def graph_uid(self) -> int:
        """
        uid of the graph this path belongs to
        """

    def hops(self) -> list[list[int]]:
        """
        the hops as lists of EdgeIds (the caller's member order)
        """

    def source(self) -> int:
        """
        the first node of the walk
        """

    def to_dag(self) -> PredecessorDAG:
        """
        the equivalent PredecessorDAG (every member of every hop)
        """

class FlowRecord:
    """
    One flow of a ledger, as a snapshot: flow_id (ledger uid, id), policy_id (the policy's id),
    sources, destinations, admission, entry_id (a T1, T2 or Pinned policy's program entry, else None),
    route (PredecessorDAG: the flow's route), min_route_cost (its minimum walk cost; None when empty),
    offered (Lossy only, else None), delivered, edge_flows [(EdgeId, units)] (carried load) and
    edge_drops [(EdgeId, units)] (cumulative, ascending EdgeId), edge_reservations [(EdgeId, units)]
    (T1's strict holds, ascending EdgeId). offered and the edge_drops units are exact Python ints (they
    may exceed 64 bits).
    """

    def __init__(self, _no_public_constructor: typing.Never, /) -> None: ...
    def __repr__(self) -> str: ...
    @property
    def admission(self) -> Admission:
        """
        LOSSLESS or LOSSY
        """

    @property
    def delivered(self) -> int:
        """
        the units the flow delivers to its destinations (its total)
        """

    @property
    def destinations(self) -> list[int]:
        """
        the flow's destination node ids, ascending
        """

    @property
    def edge_drops(self) -> list[tuple[int, int]]:
        """
        Lossy: the units dropped per edge [(EdgeId, units)], summed over calls, positive entries in
        ascending EdgeId (exact ints); empty for a Lossless flow
        """

    @property
    def edge_flows(self) -> list[tuple[int, int]]:
        """
        the flow's allocation [(EdgeId, units)], positive entries in ascending EdgeId
        """

    @property
    def edge_reservations(self) -> list[tuple[int, int]]:
        """
        T1: the strict holds of an up tunnel [(EdgeId, units)], its request on every link of its path,
        ascending EdgeId; empty for every other record (holds never police traffic)
        """

    @property
    def entry_id(self) -> int | None:
        """
        the program entry the record carries (T1, T2 and Pinned: the entry identifier); None for the
        records kept per roster source (ShortestPaths, CentralizedTE)
        """

    @property
    def flow_id(self) -> tuple[int, int]:
        """
        (ledger uid, id): ids are assigned in creation order and never reused in the ledger
        """

    @property
    def min_route_cost(self) -> int | None:
        """
        the cost of the cheapest walk of the route from the sources to the destinations; None when
        the route has no such walk
        """

    @property
    def offered(self) -> int | None:
        """
        Lossy: the units offered to the flow, summed over calls (an exact int, may exceed 64 bits);
        None for a Lossless flow
        """

    @property
    def policy_id(self) -> int:
        """
        the ledger-local id of the policy that created the flow
        """

    @property
    def route(self) -> PredecessorDAG:
        """
        the flow's route (an immutable PredecessorDAG shared with the ledger)
        """

    @property
    def sources(self) -> list[int]:
        """
        the flow's source node ids, ascending
        """

class FlowLedger:
    """
    The integer ledger of one graph, in two projections that are never added: carried traffic (the
    service baseline; carried() is the sum of every record's edge_flows and residual() = baseline -
    carried()) and the strict holds of T1 tunnels (the reservation baseline; reserved() is the sum of
    every record's edge_reservations and unreserved() = reservation baseline - reserved()), all int64[E]
    snapshots. Traffic may use bandwidth that is held: holds never police. The baselines are replaced
    only by reset(); the records are the flows placed on it by policies (the ledger has no public
    placement, removal or reservation). reset() and copy() are rejected while a policy is bound (a
    FlowPolicy built on this ledger still exists).
    """

    def __copy__(self) -> FlowLedger: ...
    def __deepcopy__(self, memo: typing.Any) -> FlowLedger: ...
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
    ) -> None:
        """
        service_baseline and reservation_baseline: the baselines as integer units, one entry per edge in
        [0, the edge's capacity] (a physical baseline goes through to_units first; service_baseline
        checked first); None: the capacities. A wrong type (not a 1-D sequence of integers, a non-integer
        entry) raises TypeError, a wrong length or an entry outside its range ValueError.
        """

    def carried(self) -> numpy.typing.NDArray[numpy.int64]:
        """
        int64[E] snapshot: the sum of every record's edge_flows
        """

    def copy(self) -> FlowLedger:
        """
        a new ledger (new uid) with the same baselines, records and holds, ids kept under the new uid;
        ValueError while a policy is bound
        """

    def edge_capacity_units(self) -> numpy.typing.NDArray[numpy.int64]:
        """
        int64[E] the graph's capacities in units (read-only view)
        """

    def graph(self) -> Graph:
        """
        the ledger's graph
        """

    def record(
        self, flow_id: tuple[typing.SupportsIndex, typing.SupportsIndex]
    ) -> FlowRecord:
        """
        one record by FlowId (ValueError for a FlowId of another ledger, a retired or an unknown id)
        """

    def records(self) -> list[FlowRecord]:
        """
        live records as snapshots in FlowId order: those of bound policies and orphans (records left by a
        policy destroyed without remove())
        """

    def reserved(self) -> numpy.typing.NDArray[numpy.int64]:
        """
        int64[E] snapshot: strict holds per edge (T1), the sum of every record's edge_reservations
        """

    def reset(
        self,
        service_baseline: collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer]
        | None = None,
        *,
        reservation_baseline: collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer]
        | None = None,
    ) -> None:
        """
        New baselines (None: the capacities); orphaned records with their allocations and holds are
        cleared and their ids retired (ids are never reused). ValueError while a policy is bound.
        """

    def residual(self) -> numpy.typing.NDArray[numpy.int64]:
        """
        int64[E] snapshot: service_baseline - carried()
        """

    def saturated_boundary_edges(
        self,
        sources: typing.SupportsIndex
        | collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer],
        eligibility: Eligibility | None = None,
    ) -> list[int]:
        """
        Saturation report: the eligible edges with capacity >= 1 unit that leave the set of nodes the
        allowed sources can reach (following eligible edges forward with residual >= 1 unit), in
        ascending EdgeId. It ignores the destinations, so it is not a source-destination cut.
        """

    def uid(self) -> int:
        """
        process-unique ledger uid (the first half of every FlowId)
        """

    def unreserved(self) -> numpy.typing.NDArray[numpy.int64]:
        """
        int64[E] snapshot: reservation baseline minus strict holds, per edge
        """

class TunnelStatus:
    """
    A tunnel's state. UP; INSUFFICIENT_BANDWIDTH and COST_CEILING_EXCEEDED: T1's setup left it down (no path
    passes the bandwidth filter, or the cheapest exceeds the cost ceiling): no path, weight 0, and it
    stays down in every forwarding of its program; and in one forwarding only (a by_tunnel row)
    FAILED: its path touches a node or an edge the forwarding policy's Eligibility excludes, so it
    gets nothing and its head splits over its other entries (replay), and BLACKHOLED: up, but its
    path reaches a black hole of stale forwarding, where its traffic is dropped.

    Members:

      UP

      INSUFFICIENT_BANDWIDTH

      COST_CEILING_EXCEEDED

      FAILED

      BLACKHOLED
    """

    BLACKHOLED: typing.ClassVar[TunnelStatus]
    COST_CEILING_EXCEEDED: typing.ClassVar[TunnelStatus]
    FAILED: typing.ClassVar[TunnelStatus]
    INSUFFICIENT_BANDWIDTH: typing.ClassVar[TunnelStatus]
    UP: typing.ClassVar[TunnelStatus]
    __members__: typing.ClassVar[dict[str, TunnelStatus]]

    @typing.overload
    def __eq__(self, other: TunnelStatus, /) -> bool: ...
    @typing.overload
    def __eq__(self, other: typing.Any, /) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __init__(self, value: typing.SupportsInt | typing.SupportsIndex) -> None: ...
    def __int__(self) -> int: ...
    @typing.overload
    def __ne__(self, other: TunnelStatus, /) -> bool: ...
    @typing.overload
    def __ne__(self, other: typing.Any, /) -> bool: ...
    def __repr__(self) -> str: ...
    def __setstate__(
        self, state: typing.SupportsInt | typing.SupportsIndex, /
    ) -> None: ...
    def __str__(self) -> str: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class Tunnel:
    """
    One tunnel (a path entry) of a ForwardingProgram, read-only. entry_id: the entry identifier (T1 and
    T2: j*H + i for round j and roster position i, T2's dead heads leaving holes; CentralizedTE and ShortestPaths
    OPTIMIZED: the path's position in the decomposition order; supplied routes: the route's position in
    the list), never a position in program.tunnels; head; local_index: its position among the head's entries of both kinds;
    destination: where its path ends; edges: its EdgeIds in path order; path_cost: the sum of their
    costs (all three None for a tunnel that is down); status: UP, INSUFFICIENT_BANDWIDTH or
    COST_CEILING_EXCEEDED; requested_bandwidth and reserved_bandwidth: T1's strict bandwidth (its request,
    and what it holds on every link of its path: the request when UP, else 0; 0 for the other models);
    planned_bandwidth: T2's planned bandwidth
    u (metadata, never held; 0 for the other models); weight: its integer share among the head's
    entries (T1: 1 up, 0 down; T2: 1; a controller path: its solved amount).
    """

    def __init__(self, _no_public_constructor: typing.Never, /) -> None: ...
    def __repr__(self) -> str: ...
    @property
    def destination(self) -> int | None: ...
    @property
    def edges(self) -> tuple[int, ...] | None:
        """
        the path's EdgeIds in path order; None for a tunnel that is down
        """

    @property
    def entry_id(self) -> int: ...
    @property
    def head(self) -> int: ...
    @property
    def local_index(self) -> int: ...
    @property
    def path_cost(self) -> int | None: ...
    @property
    def planned_bandwidth(self) -> int: ...
    @property
    def requested_bandwidth(self) -> int: ...
    @property
    def reserved_bandwidth(self) -> int: ...
    @property
    def status(self) -> TunnelStatus: ...
    @property
    def weight(self) -> int: ...

class RouteDefinition:
    """
    A supplied route for ForwardingProgram.from_routes: a Path with one edge per hop (a tunnel) or a
    PredecessorDAG (a split route) starting at `head`; weight is its integer share among the head's routes
    (0 allowed: it carries nothing). from_routes validates it.
    """

    def __init__(
        self,
        head: typing.SupportsIndex,
        route: Path | PredecessorDAG,
        *,
        weight: typing.SupportsIndex = 1,
    ) -> None:
        """
        RouteDefinition(head, route, *, weight=1): head a NodeId, route a Path or a PredecessorDAG (TypeError
        otherwise), weight an integer (checked in [0, 2^62) by from_routes).
        """

    def __repr__(self) -> str: ...
    @property
    def head(self) -> int:
        """
        the route's head (a source of the program)
        """

    @property
    def route(self) -> Path | PredecessorDAG:
        """
        the Path (a tunnel) or the PredecessorDAG (a split route)
        """

    @property
    def weight(self) -> int:
        """
        the route's integer share among the head's routes
        """

class ForwardingProgram:
    """
    An immutable value: one demand's routes and integer weights. sources and destinations (the normalized
    sets); tunnels: the path entries in ascending index (split routes are not listed); source_weights: the
    captured source supplies ((source, weight) for every roster source, ascending) of the optimized
    controllers, None for equal sources (round robin over the roster). Each head splits its share by
    highest averages over its entries' weights; identical routes are never merged. ShortestPaths and CentralizedTE
    produce one (FlowPolicy.program, PlacementResult.program), and so do T1 and T2 (their tunnels);
    from_routes builds one from given routes; Pinned forwards one. No public constructor.
    """

    def __init__(self, _no_public_constructor: typing.Never, /) -> None: ...
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
    ) -> ForwardingProgram:
        """
        One shared hop-by-hop program on a supplied PredecessorDAG of the graph -> ForwardingProgram: the sources
        allocated equally and pooled at every shared router exactly as ShortestPaths EQUAL/WEIGHTED (one
        split context), the DAG trimmed to the walks from each source to its first destinations. ValueError for
        the first failure, in this order: the endpoint sets; `dag belongs to another graph`; allocation OPTIMIZED
        (`from_dag allocation must be EQUAL or WEIGHTED`); weight_basis only with allocation=WEIGHTED, EXPLICIT
        with edge_weights, RESIDUAL_SNAPSHOT with residual; edge_weights' length and entries in [0, 2^62) and
        residual's (one entry per edge in [0, capacity]). A dag of another type raises TypeError. The roster is
        the sources that reach a destination in the DAG (possibly none): one split route per roster source,
        identified by roster position, no tunnel; a replay rebuilds the roster on the DAG.
        """

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
    ) -> ForwardingProgram:
        """
        Validate and freeze supplied routes -> ForwardingProgram (ValueError for the first failure, in this
        order): the endpoint sets; `routes must not be empty`; when no route is a PredecessorDAG, any option given
        (`<allocation, weight_basis, edge_weights or residual> applies only to PredecessorDAG routes`); allocation
        OPTIMIZED (`from_routes allocation must be EQUAL or WEIGHTED`); weight_basis only with allocation=WEIGHTED,
        EXPLICIT with edge_weights, RESIDUAL_SNAPSHOT with residual; edge_weights' length and entries in [0,
        2^62) and residual's (one entry per edge in [0, capacity]); then each route in list order: its
        route's graph, its head a source, a Path with one edge per hop that starts at the head, ends at a
        destination and meets none before, a PredecessorDAG with a walk from the head to a destination, its weight
        in [0, 2^62). Order and duplicates are kept; a route of weight 0 is validated in full. allocation
        (None: EQUAL) is the interior split of the PredecessorDAG routes; each split route is its own
        split context (a router splits only that route's share). The roster is the sources with a route.
        """

    @property
    def destinations(self) -> tuple[int, ...]:
        """
        the destination set, ascending
        """

    @property
    def source_weights(self) -> tuple[tuple[int, int], ...] | None:
        """
        the captured source weights ((source, weight), ascending), or None for equal sources
        """

    @property
    def sources(self) -> tuple[int, ...]:
        """
        the source set, ascending
        """

    @property
    def tunnels(self) -> tuple[Tunnel, ...]:
        """
        the program's tunnels (Tunnel), ascending index
        """

class AbsoluteCost:
    """
    Cost ceiling: an absolute, inclusive path-cost budget (cost >= 0).
    """

    def __eq__(self, arg0: typing.Any) -> bool: ...
    def __hash__(self) -> int: ...
    def __init__(self, cost: typing.SupportsIndex) -> None:
        """
        AbsoluteCost(cost): a 64-bit integer budget in cost units (checked >= 0 where it is used)
        """

    def __repr__(self) -> str: ...
    @property
    def cost(self) -> int:
        """
        the budget: a path qualifies when its cost is at most this
        """

class CostStretch:
    """
    Cost ceiling: a path qualifies when cost * d <= n * base, exactly, where n/d is the factor as
    one reduced fraction and base the static shortest cost under the Eligibility.
    """

    def __eq__(self, arg0: typing.Any) -> bool: ...
    def __hash__(self) -> int: ...
    def __init__(
        self, factor: fractions.Fraction | float | numpy.floating | typing.SupportsIndex
    ) -> None:
        """
        CostStretch(factor): in this order, a bool (or numpy.bool_) raises TypeError; a fractions.Fraction
        is taken as it is; a Python float or any NumPy floating value is Fraction(str(f)), the shortest
        decimal its own type prints (1.2, np.float64(1.2) and np.float32(1.2) are all 6/5), and a NaN or
        an infinity raises ValueError; an integer, NumPy integer or any __index__ object is n/1 exactly;
        anything else raises TypeError. Then a reduced numerator outside [-2^63, 2^63) or denominator
        outside [1, 2^63) raises ValueError. A factor below 1 is valid here and refused where the
        ceiling is used.
        """

    def __repr__(self) -> str: ...
    @property
    def factor(self) -> fractions.Fraction:
        """
        the factor as a fractions.Fraction (built on each access)
        """

class TunnelCount:
    """
    The tunnel count of T1 (DistributedTE) and T2, per head: per_head in [1, 10000] (checked when a
    policy or an evaluator uses it). up_to=True (T2 only; T1 rejects it): one common count k in
    1..per_head for every head: the smallest common count whose program carries the live heads'
    design shares; if none does, the largest live-share score, then the smallest count. It does not
    promise whole-demand lossless delivery.
    """

    def __eq__(self, arg0: typing.Any) -> bool: ...
    def __hash__(self) -> int: ...
    def __init__(self, per_head: typing.SupportsIndex, *, up_to: bool = False) -> None:
        """
        TunnelCount(per_head, *, up_to=False): a 64-bit integer per_head (a bool, a float or None raises
        TypeError and an integer outside 64 bits ValueError here; [1, 10000] is checked where it is
        used) and a bool
        """

    def __repr__(self) -> str: ...
    @property
    def per_head(self) -> int:
        """
        the tunnels per head (the largest count for up_to)
        """

    @property
    def up_to(self) -> bool:
        """
        True: the smallest count up to per_head whose program carries the live heads' design shares (else the
        best live-share score, then the smallest); T2 only
        """

class ShortestPaths:
    """
    Routing like IGP forwarding (OSPF, IS-IS routers) on each roster source's static shortest-path
    DAG (read from the ShortestPathCache; the route choice is applied at the epoch's first positive call),
    fixed for the epoch; one record per roster source. EQUAL (ECMP, the default) and WEIGHTED
    (WCMP) pool the sources (source allocation EQUAL): where their routes meet, a router splits
    their traffic once. OPTIMIZED is a controller on the union of those routes (source allocation
    OPTIMIZED by default, or EQUAL). An edge of capacity 0 stays in the DAG: EQUAL still sends it
    its share, which it cannot carry, and WEIGHTED gives it weight 0.
    """

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
    ) -> None:
        """
        ShortestPaths(*, allocation=EQUAL, source_allocation=None, weight_basis=None, edge_weights=None,
        tie_break=FEWER_HOPS, multipath=True, parallel_edges=True): ShortestPaths() is ECMP.
        source_allocation None resolves to EQUAL for EQUAL and WEIGHTED (the only allowed value) and to
        OPTIMIZED for OPTIMIZED; weight_basis (WEIGHTED only; None: CAPACITY) and edge_weights (one
        integer per edge, given exactly when weight_basis=EXPLICIT; one array shared by every policy
        built from this model). The combination is checked when a policy or an evaluator uses the model
        (ValueError, in the order of the model matrix). Here, a value of the wrong type (a multipath or parallel_edges that
        is not a bool, an enumeration of another class, an edge_weights that is not a 1-D sequence of
        integers or has a bool or a float entry) raises TypeError, and an edge_weights entry beyond
        int64 ValueError.
        """

    def __repr__(self) -> str: ...
    @property
    def allocation(self) -> RouteAllocation:
        """
        the RouteAllocation: EQUAL, WEIGHTED or OPTIMIZED
        """

    @property
    def edge_weights(self) -> tuple[int, ...] | None:
        """
        the explicit weights (one int per edge, a tuple), or None
        """

    @property
    def multipath(self) -> bool:
        """
        True: every equal-cost route; False: one canonical path per roster source toward the destination set
        """

    @property
    def parallel_edges(self) -> bool:
        """
        True: every eligible parallel edge; False: one per (u, v) pair by the tie rule
        """

    @property
    def source_allocation(self) -> SourceAllocation | None:
        """
        the SourceAllocation, or None (a policy's config: resolved to the allocation's default)
        """

    @property
    def tie_break(self) -> RouteTieBreak:
        """
        the RouteTieBreak of equal-cost choices
        """

    @property
    def weight_basis(self) -> WeightBasis | None:
        """
        the WeightBasis, or None (a policy's config: CAPACITY resolved for WEIGHTED)
        """

class DistributedTE:
    """
    T1: every head sets up its own strict tunnels by CSPF, no controller (like RSVP-TE head ends).
    The design volume D is split equally over the roster heads and each head's share equally over its
    tunnel_count.per_head tunnels; in round order (round j, heads ascending: tunnel index j*H + i) each tunnel
    keeps the allowed links whose unreserved capacity (the ledger's, minus the holds of every other
    owner and of this setup's earlier tunnels) is at least max(request, min_unreserved_bandwidth) (a
    request of 0 with a minimum of 0 has no filter) and takes the shortest path in the tie_break's
    single-path order; it comes up holding exactly its request on every link (UP), or stays down
    holding nothing: INSUFFICIENT_BANDWIDTH (no path passes the filter) or COST_CEILING_EXCEEDED (the cheapest exceeds
    cost_ceiling: AbsoluteCost is inclusive; CostStretch(f) is f times the smallest static shortest cost over
    all roster heads). Traffic enters equally at the heads and splits equally over each head's up
    tunnels; a head with no tunnel up loses its share. Holds never police traffic.
    """

    def __init__(
        self,
        tunnel_count: TunnelCount,
        *,
        cost_ceiling: AbsoluteCost | CostStretch | None = None,
        min_unreserved_bandwidth: typing.SupportsIndex = 0,
        tie_break: RouteTieBreak = ...,
        parallel_edges: bool = True,
    ) -> None:
        """
        DistributedTE(tunnel_count, *, cost_ceiling=None, min_unreserved_bandwidth=0, tie_break=FEWER_HOPS,
        parallel_edges=True): the fields are checked when a policy is built with them (the count in
        [1, 10000] and exact, the ceiling, min_unreserved_bandwidth in [0, 2^62) units). tie_break and
        parallel_edges are the CSPF choice, scored on unreserved bandwidth. Here, a
        min_unreserved_bandwidth outside 64 bits raises ValueError; one of the wrong type, a tunnel_count
        that is not a TunnelCount, a ceiling of the wrong class, a parallel_edges that is not a bool and an
        unknown keyword raise TypeError.
        """

    def __repr__(self) -> str: ...
    @property
    def cost_ceiling(self) -> AbsoluteCost | CostStretch | None:
        """
        the inclusive route cost ceiling (AbsoluteCost or CostStretch); None: no ceiling
        """

    @property
    def min_unreserved_bandwidth(self) -> int:
        """
        units in [0, 2^62): a search minimum only (each attempt keeps the links with at least
        max(request, this) unreserved units), never a reservation
        """

    @property
    def parallel_edges(self) -> bool:
        """
        True: every eligible parallel edge; False: one per (u, v) pair by the tie rule
        """

    @property
    def tie_break(self) -> RouteTieBreak:
        """
        the RouteTieBreak of the CSPF choice
        """

    @property
    def tunnel_count(self) -> TunnelCount:
        """
        the tunnel count per head (a TunnelCount)
        """

class CentralizedTE:
    """
    Routing: an ideal central controller that sees the whole network (with `tunnel_count`: T2, below).
    Every positive call re-solves the minimum-cost flow of what the policy already delivers plus this
    call's request (delivered() + offer; what an earlier call could not place is not requested
    again) over the eligible graph, or exactly over `support` (a PredecessorDAG of the policy's graph, trimmed to the
    walks from the sources to the allowed destinations under the Eligibility; every other edge
    carries nothing), on the free capacity plus its own allocation (it may move its own traffic,
    never another's). It takes no route choice. Its source_allocation is required: OPTIMIZED, the
    controller chooses how much enters at each source; EQUAL: every roster source sends its
    round-robin share of the largest total whose shares all fit. Its program is the solve cut into
    paths per source, each a tunnel weighted by its amount. A positive call on a network beyond the
    flow core's int32 numbering raises ValueError (see FlowPolicy.place).
    T2, CentralizedTE(source_allocation=..., tunnel_count=TunnelCount(N)) (no route choice either; with
    support, every step uses only the support's edges allowed by the policy's mask): a controller (like
    a PCE) lays out, for every live head (one reaching a destination over links with free capacity
    >= 1; the other heads are dead and lose their share), N single-path tunnels with one common
    whole-unit planned bandwidth: first as large as possible (U), then for a finite design no larger
    than needed (u), then the cheapest routes, then with prefer_headroom=True the routes that leave
    the most spare capacity, then one fixed unit-by-unit tie. TunnelCount(L, up_to=True): the smallest
    common count k <= L whose program carries the live heads' design shares (else the largest
    live-share score, then the smallest count); no promise of whole-demand lossless delivery.
    Traffic splits equally over the heads and over each head's tunnels; planned bandwidth is never
    held.
    """

    def __init__(
        self,
        *,
        source_allocation: SourceAllocation,
        support: PredecessorDAG | None = None,
        tunnel_count: TunnelCount | None = None,
        prefer_headroom: bool = False,
    ) -> None:
        """
        CentralizedTE(*, source_allocation, support=None, tunnel_count=None, prefer_headroom=False):
        source_allocation is required (EQUAL or OPTIMIZED; None or a missing argument raises TypeError);
        support a PredecessorDAG of the policy's graph (checked when a policy is built with it: `support
        belongs to another graph`) or None for every eligible edge; tunnel_count a TunnelCount (T2;
        the count in [1, 10000] is checked when a policy is built with it) or None; prefer_headroom (T2
        only; True without tunnel_count raises `prefer_headroom requires tunnel_count` when a policy is
        built with it): the headroom penalty as T2's objective after route cost
        """

    def __repr__(self) -> str: ...
    @property
    def prefer_headroom(self) -> bool:
        """
        True: T2 minimizes the headroom penalty after route cost (T2 only)
        """

    @property
    def source_allocation(self) -> SourceAllocation:
        """
        the SourceAllocation: EQUAL or OPTIMIZED
        """

    @property
    def support(self) -> PredecessorDAG | None:
        """
        the supplied support (a PredecessorDAG), or None
        """

    @property
    def tunnel_count(self) -> TunnelCount | None:
        """
        T2's tunnel count per head (a TunnelCount), or None
        """

class Pinned:
    """
    Routing: forwarding over a given program, computing nothing (the program's endpoints must be
    the policy's). The policy's Eligibility repairs failures: a tunnel touching an excluded element is
    FAILED, a split route loses its excluded members, and the roster is rebuilt from the program's
    basis (a head that no longer reaches a destination leaves the equal pool). black_holes: the nodes
    and edges this Eligibility excludes are silent failures; the program is forwarded unchanged and
    traffic reaching them is dropped there (stale forwarding; None: no black holes). One record per
    program entry.
    """

    def __init__(
        self, program: ForwardingProgram, *, black_holes: Eligibility | None = None
    ) -> None:
        """
        Pinned(program, *, black_holes=None): a ForwardingProgram and an Eligibility of the same graph (both
        kept by reference; checked when a policy is built with them)
        """

    def __repr__(self) -> str: ...
    @property
    def black_holes(self) -> Eligibility | None:
        """
        the black holes (an Eligibility), or None
        """

    @property
    def program(self) -> ForwardingProgram:
        """
        the forwarded program (shared, immutable)
        """

class FlowPolicyConfig:
    """
    One demand's settings: the routing model, which holds every option it reads (ShortestPaths,
    DistributedTE (T1), CentralizedTE, CentralizedTE(tunnel_count=...) (T2) or Pinned), and the admission
    (LOSSLESS by default; LOSSY for ShortestPaths EQUAL and WEIGHTED, T1, T2 and Pinned). The
    combination is checked when a policy or an evaluator uses the configuration (ValueError), in
    this order: the enumeration values; the model matrix (route and source allocation, admission,
    the ineffective tie rule, the pooled single-path rule, prefer_headroom, the weight combinations,
    and last a ShortestPathCache on a model that does not read it); the tunnel count, an up_to count
    for DistributedTE, the ceilings, min_unreserved_bandwidth and edge_weights' length and entries;
    the tunnel limit (roster heads times the tunnel count at most 2^30 - 1); a Pinned program's
    endpoints.
    """

    def __init__(
        self,
        *,
        routing: ShortestPaths | DistributedTE | CentralizedTE | Pinned,
        admission: Admission = ...,
    ) -> None: ...
    def __repr__(self) -> str: ...
    @property
    def admission(self) -> Admission:
        """
        LOSSLESS or LOSSY
        """

    @property
    def routing(self) -> ShortestPaths | DistributedTE | CentralizedTE | Pinned:
        """
        the routing model: ShortestPaths, DistributedTE, CentralizedTE or Pinned
        """

class OutputOptions:
    """
    Optional outputs of a placement (all false by default; absent outputs are None):
    edge_flows (the policy's own allocation per edge after the call), residual and
    saturated_boundary_edges (the ledger's after the call), edge_reservations (the policy's own
    strict holds per edge after the call; zeros for a model without holds), and min_cut (CentralizedTE with
    the optimized source allocation, when its solve reached the maximum flow; None for every other
    model, source allocation or call).
    """

    def __init__(
        self,
        *,
        edge_flows: bool = False,
        residual: bool = False,
        saturated_boundary_edges: bool = False,
        edge_reservations: bool = False,
        min_cut: bool = False,
    ) -> None: ...
    @property
    def edge_flows(self) -> bool:
        """
        request PlacementResult.edge_flows
        """

    @property
    def edge_reservations(self) -> bool:
        """
        request PlacementResult.edge_reservations
        """

    @property
    def min_cut(self) -> bool:
        """
        request PlacementResult.min_cut
        """

    @property
    def residual(self) -> bool:
        """
        request PlacementResult.residual
        """

    @property
    def saturated_boundary_edges(self) -> bool:
        """
        request PlacementResult.saturated_boundary_edges
        """

class Termination:
    """
    How a call ended, the first that applies: COMPLETE for a zero request; NO_ROUTE for a positive
    request with an empty roster (no eligible source statically reaches an allowed destination);
    NO_ACTIVE_WEIGHT when the program has an entry with a path (an up tunnel, a controller path or
    a split route) but no head can use one (no positive weight: WEIGHTED over edges of capacity 0
    only, routes of weight 0, every up tunnel failed); a program whose tunnels are all down is not
    this case; a LOSSY finite request: COMPLETE (delivery and losses account for the whole offer,
    which does not mean nothing was lost); otherwise COMPLETE when the whole request was delivered,
    else EXHAUSTED (an infinite request is EXHAUSTED). A controller whose maximum is 0 over existing
    routes is EXHAUSTED, not NO_ROUTE.

    Members:

      COMPLETE

      NO_ROUTE

      EXHAUSTED

      NO_ACTIVE_WEIGHT
    """

    COMPLETE: typing.ClassVar[Termination]
    EXHAUSTED: typing.ClassVar[Termination]
    NO_ACTIVE_WEIGHT: typing.ClassVar[Termination]
    NO_ROUTE: typing.ClassVar[Termination]
    __members__: typing.ClassVar[dict[str, Termination]]

    @typing.overload
    def __eq__(self, other: Termination, /) -> bool: ...
    @typing.overload
    def __eq__(self, other: typing.Any, /) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __init__(self, value: typing.SupportsInt | typing.SupportsIndex) -> None: ...
    def __int__(self) -> int: ...
    @typing.overload
    def __ne__(self, other: Termination, /) -> bool: ...
    @typing.overload
    def __ne__(self, other: typing.Any, /) -> bool: ...
    def __repr__(self) -> str: ...
    def __setstate__(
        self, state: typing.SupportsInt | typing.SupportsIndex, /
    ) -> None: ...
    def __str__(self) -> str: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class SourceFlow:
    """
    One normalized input source's forwarding in the epoch (read-only; zero rows included: a
    masked or unreachable source): source, offered (LOSSY: its injections; None when
    LOSSLESS), admitted (LOSSLESS: what it delivers; LOSSY: offered), delivered, source_drops
    (offered units lost at the source for lack of a usable route) and total_edge_drops (its losses
    on links). LOSSY: offered = delivered + source_drops + total_edge_drops. Exact ints.
    """

    def __init__(self, _no_public_constructor: typing.Never, /) -> None: ...
    def __repr__(self) -> str: ...
    @property
    def admitted(self) -> int:
        """
        LOSSLESS: the units it delivers; LOSSY: offered
        """

    @property
    def delivered(self) -> int:
        """
        the units it delivers to the destinations
        """

    @property
    def offered(self) -> int | None:
        """
        LOSSY: the units this source injected; None when LOSSLESS
        """

    @property
    def source(self) -> int:
        """
        the source node id
        """

    @property
    def source_drops(self) -> int:
        """
        LOSSY: its offered units lost at the source (no usable route); 0 when LOSSLESS
        """

    @property
    def total_edge_drops(self) -> int:
        """
        LOSSY: its units dropped on links; 0 when LOSSLESS
        """

class TunnelFlow:
    """
    One program tunnel's forwarding in the epoch (read-only): entry_id (its entry identifier) and
    head, status (TunnelStatus: UP; INSUFFICIENT_BANDWIDTH or COST_CEILING_EXCEEDED for a tunnel
    T1's setup left down; FAILED, BLACKHOLED in this forwarding), assigned (its share of the
    injections: LOSSY the shares sent to it, LOSSLESS what it carries), delivered and
    total_edge_drops (its losses on links). Exact ints (assigned and total_edge_drops may exceed 64
    bits).
    """

    def __init__(self, _no_public_constructor: typing.Never, /) -> None: ...
    def __repr__(self) -> str: ...
    @property
    def assigned(self) -> int:
        """
        its share of the injections in the epoch (LOSSY: sent to it; LOSSLESS: what it carries)
        """

    @property
    def delivered(self) -> int:
        """
        the units it delivers
        """

    @property
    def entry_id(self) -> int:
        """
        the tunnel's entry identifier
        """

    @property
    def head(self) -> int:
        """
        the tunnel's head
        """

    @property
    def status(self) -> TunnelStatus:
        """
        UP, a down tunnel's setup status (INSUFFICIENT_BANDWIDTH, COST_CEILING_EXCEEDED), FAILED or
        BLACKHOLED in this forwarding
        """

    @property
    def total_edge_drops(self) -> int:
        """
        its units dropped on links
        """

class TunnelPlan:
    """
    The tunnel models' setup summary (read-only): tunnels_per_head (T2 up to: the chosen count),
    design_demand (the design volume D in units), design_by_head ((head, RR(D, H, i)) for every
    roster head in ascending NodeId) and total_reserved_bandwidth (T1's strict holds, counted once
    per up tunnel: the sum of their requests; 0 for T2). T2 only (None for T1): heads_with_capacity,
    insufficient_resolution (live heads exist but no whole unit per tunnel fits, U == 0: the
    program has no tunnels), max_fractional_bandwidth_per_tunnel (the exact fluid optimum s*, a
    Fraction), max_bandwidth_per_tunnel (U), planned_bandwidth_per_tunnel (u),
    total_planned_bandwidth (K*u), total_tunnel_cost (the route-cost sum, can exceed 2^64),
    headroom_penalty (with prefer_headroom only) and max_lossless_offer (the program's lossless
    maximum; a dead head's share blocks it). Exact ints.
    """

    def __init__(self, _no_public_constructor: typing.Never, /) -> None: ...
    def __repr__(self) -> str: ...
    @property
    def design_by_head(self) -> tuple[tuple[int, int], ...]:
        """
        ((head, share), ...): every roster head's share RR(D, H, i), ascending NodeId
        """

    @property
    def design_demand(self) -> int | float:
        """
        the design volume D in units; float('inf') when unlimited
        """

    @property
    def headroom_penalty(self) -> int | None:
        """
        T2 with prefer_headroom: the sum over the links' tunnel slots j of floor(2^32 (B - V + j*u) / B)
        (B: static capacity; V: free capacity); None without prefer_headroom and for T1
        """

    @property
    def heads_with_capacity(self) -> tuple[int, ...] | None:
        """
        T2: the heads that reach a destination over links with free capacity >= 1 (ascending); None for T1
        """

    @property
    def insufficient_resolution(self) -> bool | None:
        """
        T2: True when live heads exist but no whole unit per tunnel fits (max_bandwidth_per_tunnel 0,
        no tunnels)
        """

    @property
    def max_bandwidth_per_tunnel(self) -> int | None:
        """
        T2: U, the largest common whole-unit size
        """

    @property
    def max_fractional_bandwidth_per_tunnel(self) -> fractions.Fraction | None:
        """
        T2: the exact fluid optimum s* (a fractions.Fraction; 0 without a live head); None for T1
        """

    @property
    def max_lossless_offer(self) -> int | None:
        """
        T2: F, the program's lossless maximum on the call's view (a dead head's share blocks it)
        """

    @property
    def planned_bandwidth_per_tunnel(self) -> int | None:
        """
        T2: u, every tunnel's planned bandwidth
        """

    @property
    def total_planned_bandwidth(self) -> int | None:
        """
        T2: K*u over the K tunnels
        """

    @property
    def total_reserved_bandwidth(self) -> int:
        """
        the strict holds of the up tunnels, each tunnel's request counted once
        """

    @property
    def total_tunnel_cost(self) -> int | None:
        """
        T2: the sum of cost * tunnel count over the links (exact; can exceed 2^64)
        """

    @property
    def tunnels_per_head(self) -> int:
        """
        the tunnels per head (T2 up to: the chosen count)
        """

class PlacementResult:
    """
    The result of one call. This call: call_delivered (units delivered by the call; replace, forward
    and evaluate: the full new delivery), call_undelivered (request - call_delivered; float('inf') for
    an infinite request), delivery_delta (new minus old delivery, may be negative), call_edge_drops
    [(EdgeId, units)], call_source_drops, termination, routes_changed, weights_changed and
    reservations_changed. The
    demand's state after the call (totals of the epoch): program (the installed ForwardingProgram;
    evaluate: the computed one; None before the first preparation), plan (a TunnelPlan for T1 and T2;
    None for every other model), requested (float('inf') after an infinite request), offered (LOSSY: the traffic sent; None when LOSSLESS), admitted, delivered, unadmitted
    (LOSSLESS: requested - delivered, None after an infinite request; LOSSY: 0), source_drops,
    edge_drops [(EdgeId, units)], by_source (one SourceFlow per normalized input source, ascending
    NodeId), by_tunnel (one TunnelFlow per program tunnel, ascending entry_id; an alternative view, never
    added to by_source), total_cost (controllers: CentralizedTE and ShortestPaths OPTIMIZED; else None)
    and marginal_cost_curve (unrestricted CentralizedTE with the optimized source allocation only).
    LOSSY: offered = source_drops + sum(edge_drops) + delivered. The rows sum to the totals, except for
    LOSSY with an empty roster, where the whole offer is the demand's source drops and every row is
    zero. Optional outputs, None unless requested: edge_flows, residual and edge_reservations
    (int64[E]), saturated_boundary_edges, min_cut (they describe the same
    state; evaluate: the evaluated one). Amounts that accumulate are exact ints (they may exceed 64 bits).
    """

    def __init__(self, _no_public_constructor: typing.Never, /) -> None: ...
    def __repr__(self) -> str: ...
    @property
    def admitted(self) -> int:
        """
        LOSSLESS: the total carried; LOSSY: offered
        """

    @property
    def by_source(self) -> tuple[SourceFlow, ...]:
        """
        one SourceFlow per normalized input source, ascending NodeId (zero rows included)
        """

    @property
    def by_tunnel(self) -> tuple[TunnelFlow, ...]:
        """
        one TunnelFlow per program tunnel, ascending index (empty without tunnels)
        """

    @property
    def call_delivered(self) -> int:
        """
        the units placed (delivered) by this call
        """

    @property
    def call_edge_drops(self) -> list[tuple[int, int]]:
        """
        this call's drops on links [(EdgeId, units)], ascending EdgeId (empty unless LOSSY)
        """

    @property
    def call_source_drops(self) -> int:
        """
        this call's offered units lost at sources without a usable route (LOSSY)
        """

    @property
    def call_undelivered(self) -> int | float:
        """
        the request minus call_delivered, in units (an amount, not a status); float('inf') for an infinite
        request
        """

    @property
    def delivered(self) -> int:
        """
        the units the demand delivers
        """

    @property
    def delivery_delta(self) -> int:
        """
        the delivery after the call minus before it (call_delivered for a placement; may be negative
        for a full-offer call)
        """

    @property
    def edge_drops(self) -> list[tuple[int, int]]:
        """
        the epoch's drops on links [(EdgeId, units)], ascending EdgeId (empty unless LOSSY)
        """

    @property
    def edge_flows(self) -> numpy.typing.NDArray[numpy.int64] | None:
        """
        if requested: the policy's own allocation per edge after the call (int64[E], EdgeId order),
        else None
        """

    @property
    def edge_reservations(self) -> numpy.typing.NDArray[numpy.int64] | None:
        """
        if requested: the policy's own strict holds per edge after the call (int64[E], EdgeId order; zeros
        for a model without holds), else None
        """

    @property
    def marginal_cost_curve(self) -> list[tuple[int, int]] | None:
        """
        unrestricted CentralizedTE with the optimized source allocation: the marginal cost curve [(cost, units)] of the policy's
        cumulative flow (a zero call reports it unchanged; None before the epoch's first positive
        call); in ascending cost, delivering `units` more raised the minimum total flow cost by `cost`
        per unit (more delivery can reroute earlier flow, so these are not the costs of the final
        paths). None for the other models and the equal source allocation
        """

    @property
    def min_cut(self) -> list[int] | None:
        """
        if requested (OutputOptions(min_cut=True)), unrestricted CentralizedTE with the optimized source allocation, when the solve
        reached the maximum flow: the minimum cut closest to the sources, every allowed edge leaving
        the nodes the roster sources reach in the residual network (free capacity plus the policy's
        own allocation), an edge with nothing free (capacity 0, or held by others) included,
        ascending EdgeId; [] for an empty roster. None otherwise (another model, the equal source allocation, a
        support, a call below the maximum, a zero call, forward)
        """

    @property
    def offered(self) -> int | None:
        """
        LOSSY: the traffic sent in the epoch (the first-hop injection of an infinite request); None
        when LOSSLESS
        """

    @property
    def plan(self) -> TunnelPlan | None:
        """
        T1 and T2: the TunnelPlan of the program after the call; None for every other model (Pinned over
        a tunnel program too)
        """

    @property
    def program(self) -> ForwardingProgram | None:
        """
        the installed ForwardingProgram after the call (evaluate: the computed one); None before the first
        preparation
        """

    @property
    def requested(self) -> int | float:
        """
        the sum of the epoch's requests; float('inf') after an infinite one
        """

    @property
    def reservations_changed(self) -> bool:
        """
        True when a tunnel's requested or reserved amount differs from the installed program's, by
        index (an absent tunnel holds 0): only T1 sets it; forward and Pinned never do
        """

    @property
    def residual(self) -> numpy.typing.NDArray[numpy.int64] | None:
        """
        if requested: the ledger's residual after the call (int64[E], EdgeId order), else None
        """

    @property
    def routes_changed(self) -> bool:
        """
        True when the call's program differs from the installed one in its ordered heads, paths, split
        tables or duplicates (no program counts as the empty program; two programs without an entry
        are equal). An offer-only change (forward, Pinned, ShortestPaths increments) and an identical
        replace set nothing
        """

    @property
    def saturated_boundary_edges(self) -> list[int] | None:
        """
        if requested: the ledger's saturation report after the call
        (FlowLedger.saturated_boundary_edges of the policy's sources and Eligibility), else None
        """

    @property
    def source_drops(self) -> int:
        """
        LOSSY: the offered units lost at sources (no usable route, or an empty roster); 0 when LOSSLESS
        """

    @property
    def termination(self) -> Termination:
        """
        the call's outcome (Termination)
        """

    @property
    def total_cost(self) -> int | None:
        """
        controllers (CentralizedTE, ShortestPaths OPTIMIZED): the exact sum of cost * flow over the policy's
        allocation; None for the other models (T1, T2, Pinned)
        """

    @property
    def unadmitted(self) -> int | None:
        """
        LOSSLESS: requested - delivered, None after an infinite request; LOSSY: 0
        """

    @property
    def weights_changed(self) -> bool:
        """
        True when the call's program differs from the installed one in its injection (source weights)
        or any entry or split-table weight (a controller's re-solve that changes a path's amount; an
        absent entry weighs 0)
        """

class ShortestPathCacheStats:
    """
    ShortestPathCache statistics (atomic per field, not a consistent snapshot): hits, misses,
    entries, bytes (the entries' retained payload), live_values (values alive anywhere).
    """

    def __init__(self, _no_public_constructor: typing.Never, /) -> None: ...
    @property
    def bytes(self) -> int:
        """
        the retained payload of the entries (their DAGs' stored arrays)
        """

    @property
    def entries(self) -> int:
        """
        the entries held
        """

    @property
    def hits(self) -> int:
        """
        lookups answered from an entry (they persist across clear())
        """

    @property
    def live_values(self) -> int:
        """
        the cached values alive anywhere, entries included
        """

    @property
    def misses(self) -> int:
        """
        lookups without an entry: a building cache searched, a frozen one raised ShortestPathCacheMiss
        """

class ShortestPathCacheMiss(Exception):
    """
    A lookup in a frozen ShortestPathCache whose key is absent: a dedicated exception (not a
    RuntimeError or ValueError); the evaluators propagate it unchanged.
    """

class ShortestPathCache:
    """
    A cache of static shortest-path DAGs (from costs and the Eligibility only) for one graph, keyed by
    (the root set, the direction, the Eligibility uid). It stores the full DAG and applies each
    caller's route choice afterwards, so the route choice is never cached. ShortestPaths policies and the
    evaluators read it.
    Building mode (a new cache): one writer, no concurrent readers; a miss searches and inserts;
    clear() empties it (the hit and miss counters persist), and values already handed out stay valid.
    Frozen mode (freeze()): read-only and shareable by every worker; a miss counts one miss and raises
    ShortestPathCacheMiss; clear() raises ValueError; thaw() returns to building once the caller has joined
    every reader.
    When a public call fails: a building cache passed to it keeps the complete entries and statistics
    of the lookups that succeeded before the failure (never a partial entry; a kept entry never
    changes a result); a frozen cache's statistics stay unchanged, except one miss when the call
    raised ShortestPathCacheMiss. A frozen cache counts its hits only when the call succeeds, Python result
    included.
    """

    def __init__(self, graph: Graph) -> None:
        """
        An empty cache in building mode, bound to graph (by its uid).
        """

    def clear(self) -> None:
        """
        empty a building cache; ValueError when frozen
        """

    def freeze(self) -> None:
        """
        enter frozen mode (no effect when frozen)
        """

    def graph_uid(self) -> int:
        """
        uid of the graph this cache is bound to
        """

    def stats(self) -> ShortestPathCacheStats:
        """
        the current ShortestPathCacheStats (each field read atomically)
        """

    def thaw(self) -> None:
        """
        return to building mode (no effect when building); join every reader first
        """

class FlowPolicy:
    """
    The placement engine for one demand, bound to its ledger at construction. Every call (place,
    evaluate, replace, forward) is all-or-nothing: a failed call raises and leaves the ledger and the
    policy as they were (with the two limits in the module docstring). remove() takes every flow of the
    policy back out, clears its program and starts a new epoch (the request total starts again from 0);
    replace and forward start one too. A policy dropped without remove()
    leaves its records in the ledger as orphans: readable, never changed again, cleared by
    FlowLedger.reset(). sources and destinations are one NodeId or a sequence; the two sets must not
    overlap. ValueError for invalid endpoints, configurations and requests and for objects of another
    graph, ledger or policy; TypeError for a ledger (None included), config, eligibility, cache or
    outputs of the wrong class, a node id or a request of a wrong type (a bool, a float node id, a
    string) and an unknown keyword.
    """

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
    ) -> None:
        """
        Bind a policy to ledger. eligibility (None: everything allowed) and cache (None: none; only ShortestPaths
        routing reads it, another model rejects it) belong to the ledger's graph; the policy keeps the
        ledger, eligibility and cache alive (a Pinned policy shares its program and black holes). The roster
        (the sources that statically reach an allowed destination, over allowed edges whatever their
        capacity, the tunnel models' heads; Pinned: the program's heads that reach an allowed destination over its basis
        under eligibility) is fixed here. Everything is validated before the policy binds, in this order: the endpoints, the
        graph of eligibility, the program, the support, the cache and the black holes (`foreign Eligibility:
        built for another graph`, `program belongs to another graph`, `support belongs to another graph`,
        `foreign ShortestPathCache: bound to another graph`, `black_holes belongs to another graph`), the
        configuration (FlowPolicyConfig's order, the pooled single-path rule and T2's several-sources rule
        included), the tunnel limit (`tunnel limit exceeded:
        <H> heads x <N> tunnels > 1073741823`) and a Pinned program's endpoints (`program endpoints do not
        match policy endpoints`).
        """

    def delivered(self) -> int:
        """
        the sum of delivered over the policy's live records
        """

    def evaluate(
        self,
        offer: typing.SupportsIndex | float | numpy.floating,
        *,
        design_demand: typing.SupportsIndex | float | numpy.floating | None = None,
        outputs: OutputOptions | None = None,
    ) -> PlacementResult:
        """
        A full, independent computation and forwarding of offer (int units or float('inf')) from the
        policy's own released view (the residual plus its own carried load; T1 sets up on the unreserved
        capacity with its own holds released; T2 lays out on the released residual), as replace would do it -> PlacementResult (the computed
        program and plan, the flags replace would set, the outputs of the evaluated state); publishes
        nothing, creates no record, books no hold and leaves the epoch as it was. design_demand as for
        replace. LOSSY requires a finite offer.
        """

    def forward(
        self,
        offer: typing.SupportsIndex | float | numpy.floating,
        *,
        black_holes: Eligibility | None = None,
        outputs: OutputOptions | None = None,
    ) -> PlacementResult:
        """
        Forward a full new offer over the installed program, with no route computation and no setup, on the
        released view -> PlacementResult; a new epoch like replace, the program, its plan and T1's holds
        kept (no flag; T1's records are recreated with the same edge_reservations). LOSSLESS admits the
        program's lossless bound up to the offer; LOSSY sends the offer. black_holes (an Eligibility of the
        graph): the nodes and edges it excludes are black holes for this call (stale forwarding: traffic
        reaching one is dropped there; a black-holed source loses its share). In this order: Pinned raises
        `Pinned forwards with replace`; no installed program `forward requires an installed program`;
        the offer; LOSSY +inf `LOSSY forwarding requires a finite offer`; black_holes
        of another graph `black_holes belongs to another graph`.
        """

    def place(
        self,
        offer: typing.SupportsIndex | float | numpy.floating,
        outputs: OutputOptions | None = None,
    ) -> PlacementResult:
        """
        Place offer (int units in [0, 2^62), or float('inf')) all or nothing -> PlacementResult. T1 and T2
        raise `<DistributedTE or CentralizedTE with tunnel_count> does not support incremental placement`
        (also for 0).
        Requests add up within an epoch (from construction, the last remove(), replace() or forward());
        after float('inf') any finite request is rejected, 0 included. A request of 0 places nothing and
        reports the current state. ShortestPaths prepares its program at the epoch's first positive request;
        LOSSLESS admits the largest total in [admitted, admitted + offer] that the program carries with no
        loss (a refused unit is not offered again); LOSSY forwards the new round-robin (or captured-weight)
        injections, each head splitting its increment afresh (earlier drops are not offered again); a
        controller re-solves admitted + offer. A LOSSY float('inf') is ShortestPaths's first-hop experiment;
        Pinned raises `LOSSY forwarding requires a finite offer`. outputs (OutputOptions; None: none)
        requests edge_flows, residual, saturated_boundary_edges, edge_reservations and min_cut. Everything is
        validated before any change; a
        failed call raises with the ledger and the policy unchanged. A solve of the flow core (CentralizedTE
        or ShortestPaths OPTIMIZED) on a network beyond its int32 numbering raises ValueError:
        more than INT32_MAX - 2 nodes, or more than 2^30 - 1 edges, sources and destinations together. A
        LOSSY forwarding whose losses wait on each other in a cycle raises ValueError when its rounds
        alternate between two states or reach the round limit, RuntimeError when the result's check fails.
        """

    def record(
        self, flow_id: tuple[typing.SupportsIndex, typing.SupportsIndex]
    ) -> FlowRecord:
        """
        one of the policy's live records by FlowId (ValueError for a FlowId of another ledger or policy, a
        retired or an unknown id)
        """

    def records(self) -> list[FlowRecord]:
        """
        the policy's live records as snapshots in FlowId order
        """

    def remove(self) -> None:
        """
        revert every flow and hold of the policy, retire the ids, clear the program, open a new epoch
        """

    def replace(
        self,
        offer: typing.SupportsIndex | float | numpy.floating,
        *,
        design_demand: typing.SupportsIndex | float | numpy.floating | None = None,
        outputs: OutputOptions | None = None,
    ) -> PlacementResult:
        """
        Recompute from the policy's own released view and forward the full offer, published atomically ->
        PlacementResult: the old records and holds are retired (their ids never live again), new ones
        created, and a new epoch starts at this offer (delivery may go down: delivery_delta). ShortestPaths
        prepares its routes again; a controller re-solves; T1 sets its tunnels up again for the design and
        books their holds; T2 lays its tunnels out again for the design (planned bandwidth, never held);
        Pinned forwards its program. design_demand (T1 and T2: int units, or float('inf'); None: equal to
        the offer): T1 needs a finite design (`DistributedTE requires a finite design demand`), T2 takes an unlimited
        one (its largest common size); every other model raises `<model> does not accept design_demand`.
        In this order: the offer, LOSSY's finite offer, then the design demand (the model, its domain,
        finiteness).
        """

    @property
    def config(self) -> FlowPolicyConfig:
        """
        the configuration as validated at construction, ShortestPaths' absent
        source_allocation and weight_basis resolved to the allocation's defaults
        """

    @property
    def program(self) -> ForwardingProgram | None:
        """
        the installed ForwardingProgram; None before the first preparation and after remove()
        """

class BatchMode:
    """
    How apply_actions runs its list. ORDERED: each action as its scalar call, in list order, on the
    state the earlier actions left (earlier actions are served first). SIMULTANEOUS: every action forwarded at
    the same time on one view, losses shared on links across the demands.

    Members:

      ORDERED

      SIMULTANEOUS
    """

    ORDERED: typing.ClassVar[BatchMode]
    SIMULTANEOUS: typing.ClassVar[BatchMode]
    __members__: typing.ClassVar[dict[str, BatchMode]]

    @typing.overload
    def __eq__(self, other: BatchMode, /) -> bool: ...
    @typing.overload
    def __eq__(self, other: typing.Any, /) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __init__(self, value: typing.SupportsInt | typing.SupportsIndex) -> None: ...
    def __int__(self) -> int: ...
    @typing.overload
    def __ne__(self, other: BatchMode, /) -> bool: ...
    @typing.overload
    def __ne__(self, other: typing.Any, /) -> bool: ...
    def __repr__(self) -> str: ...
    def __setstate__(
        self, state: typing.SupportsInt | typing.SupportsIndex, /
    ) -> None: ...
    def __str__(self) -> str: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class Place:
    """
    An incremental place(offer) of policy in an apply_actions list (ORDERED only). An immutable value:
    policy (a FlowPolicy, kept alive) and offer (int units or a float, kept as given; another type
    raises TypeError here, a value outside the domain raises the call's ValueError at apply_actions).
    """

    def __init__(
        self, policy: FlowPolicy, offer: typing.SupportsIndex | float | numpy.floating
    ) -> None: ...
    @property
    def offer(self) -> int | float:
        """
        the offer as given (int units or a float)
        """

    @property
    def policy(self) -> FlowPolicy:
        """
        the FlowPolicy
        """

class Replace:
    """
    A replace(offer, design_demand) of policy in an apply_actions list; SIMULTANEOUS takes it only for
    ShortestPaths EQUAL or WEIGHTED and for Pinned. An immutable value: policy (kept alive), offer and
    design_demand (int units or a float, kept as given; None: absent; another type raises TypeError
    here, a value outside the domain raises the call's ValueError at apply_actions).
    """

    def __init__(
        self,
        policy: FlowPolicy,
        offer: typing.SupportsIndex | float | numpy.floating,
        *,
        design_demand: typing.SupportsIndex | float | numpy.floating | None = None,
    ) -> None: ...
    @property
    def design_demand(self) -> int | float | None:
        """
        the design demand as given (int units or a float); None: absent
        """

    @property
    def offer(self) -> int | float:
        """
        the offer as given (int units or a float)
        """

    @property
    def policy(self) -> FlowPolicy:
        """
        the FlowPolicy
        """

class Forward:
    """
    A forward(offer, black_holes) of policy over its installed program in an apply_actions list. An
    immutable value: policy and black_holes (an Eligibility of the policy's graph or None) kept alive,
    offer (int units or a float, kept as given; another type raises TypeError here, a value outside
    the domain raises the call's ValueError at apply_actions).
    """

    def __init__(
        self,
        policy: FlowPolicy,
        offer: typing.SupportsIndex | float | numpy.floating,
        *,
        black_holes: Eligibility | None = None,
    ) -> None: ...
    @property
    def black_holes(self) -> Eligibility | None:
        """
        the black holes (an Eligibility); None: none
        """

    @property
    def offer(self) -> int | float:
        """
        the offer as given (int units or a float)
        """

    @property
    def policy(self) -> FlowPolicy:
        """
        the FlowPolicy
        """

class EvaluationOptions:
    """
    The options of every evaluator (evaluate_demand, evaluate_pairs, evaluate_edge_failures): config
    (the configuration), eligibility (the static Eligibility; None: all allowed), cache (the
    ShortestPathCache; None: private per worker), outputs (the optional outputs of each result;
    evaluate_edge_failures: of its baseline), design_demand and
    offer (int units or a float, kept as given; another type raises TypeError). offer: a finite offer
    for every model (None: the uncapped LOSSLESS question, or for ShortestPaths EQUAL/WEIGHTED with
    LOSSY the first-hop experiment; required for DistributedTE, Pinned and, without a finite design,
    CentralizedTE with tunnel_count under LOSSY). design_demand: DistributedTE always, required and
    finite; CentralizedTE with tunnel_count: None means the offer when it is finite, else unlimited;
    every other model leaves it None. The evaluators check both after the configuration, in this order:
    an explicit offer (`offer must be an integer in [0, 2^62)`, +inf included), the design
    (`<model> does not accept design_demand`, `design_demand must be an integer in [0, 2^62) or
    permitted +inf`, `DistributedTE requires a finite design demand`), DistributedTE's missing design
    (`DistributedTE evaluation requires design_demand`), the missing LOSSY offer (`<DistributedTE or
    CentralizedTE or Pinned> LOSSY evaluation requires a finite offer`).
    """

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
    @property
    def cache(self) -> ShortestPathCache | None:
        """
        the supplied ShortestPathCache; None: a private cache per worker
        """

    @property
    def config(self) -> FlowPolicyConfig:
        """
        the FlowPolicyConfig of every evaluated policy (evaluate_edge_failures: the baseline and every deletion)
        """

    @property
    def design_demand(self) -> int | float | None:
        """
        the tunnel models' design volume as given (int units or a float); None: absent
        """

    @property
    def eligibility(self) -> Eligibility | None:
        """
        the static Eligibility; None: everything allowed
        """

    @property
    def offer(self) -> int | float | None:
        """
        the offer as given (int units or a float); None: absent
        """

    @property
    def outputs(self) -> OutputOptions:
        """
        the OutputOptions of every result (evaluate_edge_failures: of the baseline)
        """

class SensitivityEntry:
    """
    One candidate edge: edge and delivery_reduction (int units: the baseline's delivery minus the
    delivery without the edge; signed).
    """

    def __init__(self, _no_public_constructor: typing.Never, /) -> None: ...
    @property
    def delivery_reduction(self) -> int: ...
    @property
    def edge(self) -> int: ...

class SensitivityResult:
    """
    baseline (the PlacementResult of evaluate_demand()) and entries (ascending EdgeId; an
    absent edge has reduction 0).
    """

    def __init__(self, _no_public_constructor: typing.Never, /) -> None: ...
    @property
    def baseline(self) -> PlacementResult: ...
    @property
    def entries(self) -> list[SensitivityEntry]: ...

class KspOptions:
    """
    The options of ksp: k (>= 1); tie_break (RouteTieBreak) orders paths of equal cost;
    parallel_edges = True lists each parallel edge as a separate path, False keeps one per (u, v)
    pair by the tie rule (no multipath field: k sets the number of paths); cost_ceiling (None,
    AbsoluteCost(c) with c >= 0, or CostStretch(n/d) with n/d >= 1 times the static shortest
    cost under the eligibility, compared exactly, residual ignored; inclusive); residual (int64 units per
    EdgeId, each in [0, edge_capacity_units]: edges with less than 1 unit are skipped, and the values
    feed the availability; it is the only capacity filter, and without it an edge of capacity 0 is
    searched like any other); eligibility (None: everything allowed). A parallel_edges that is not a
    bool raises TypeError here; the other fields are validated by ksp, in the order its docstring
    gives.
    """

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
    @property
    def cost_ceiling(self) -> AbsoluteCost | CostStretch | None:
        """
        the inclusive cost ceiling as given: None, AbsoluteCost or CostStretch (validated by ksp)
        """

    @property
    def eligibility(self) -> Eligibility | None:
        """
        the static Eligibility; None: everything allowed
        """

    @property
    def k(self) -> typing.SupportsIndex:
        """
        the number of paths as given (validated by ksp: an integer >= 1)
        """

    @property
    def parallel_edges(self) -> bool:
        """
        True: each parallel edge a separate path; False: one per (u, v) pair by the tie rule
        """

    @property
    def residual(
        self,
    ) -> (
        collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer]
        | None
    ):
        """
        the residual units per EdgeId as given (validated by ksp); None: no capacity filter;
        priorities read static capacity
        """

    @property
    def tie_break(self) -> RouteTieBreak:
        """
        the RouteTieBreak of paths of equal cost
        """

class SpfOptions:
    """
    The options of spf: the direction; targets (one NodeId or a sequence; None: a full search)
    stops the search at the nearest targets (anycast: any nearest target will do); tie_break
    (RouteTieBreak) orders equal-cost choices; multipath (True: every equal-cost walk as a DAG,
    False: the one canonical path, targets required); parallel_edges (True: every eligible parallel
    edge, False: one per (u, v) pair by the tie rule); a HIGHER_* tie_break with multipath and
    parallel_edges both True is rejected by spf; residual (int64 units per EdgeId, each in [0,
    edge_capacity_units]) makes the search skip edges with fewer than min_residual residual units
    (an integer >= 1; without a residual it has no effect) and supplies the availability: it is the
    only capacity filter, and without it an edge of capacity 0 is searched like any other;
    eligibility (None: everything allowed). A direction outside its enumeration raises ValueError
    here, a multipath or parallel_edges that is not a bool TypeError; the other fields are validated
    by spf, in the order its docstring gives.
    """

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
    @property
    def direction(self) -> Direction:
        """
        FORWARD or REVERSE
        """

    @property
    def eligibility(self) -> Eligibility | None:
        """
        the static Eligibility; None: everything allowed
        """

    @property
    def min_residual(self) -> typing.SupportsIndex:
        """
        the gating threshold in units as given (validated by spf: an integer >= 1)
        """

    @property
    def multipath(self) -> bool:
        """
        True: every equal-cost walk (a DAG); False: the canonical single path
        """

    @property
    def parallel_edges(self) -> bool:
        """
        True: every eligible parallel edge; False: one per (u, v) pair by the tie rule
        """

    @property
    def residual(
        self,
    ) -> (
        collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer]
        | None
    ):
        """
        the residual units per EdgeId as given (validated by spf); None: no capacity filter;
        priorities read static capacity
        """

    @property
    def targets(
        self,
    ) -> (
        typing.SupportsIndex
        | collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer]
        | None
    ):
        """
        the targets as given (validated by spf); None: a full search
        """

    @property
    def tie_break(self) -> RouteTieBreak:
        """
        the RouteTieBreak of equal-cost choices
        """

class SpfResult:
    """
    The result of spf: costs (full searches: a read-only int64 array, UNREACHABLE for the nodes not
    reached; None with targets), dag (PredecessorDAG of the tight edges, the edges on some shortest walk,
    one per (u, v) pair with parallel_edges = False, each in its own direction; with targets only the
    walks between the roots and a nearest target, for REVERSE from the target to a root; for
    multipath = False the canonical single path) and
    target_cost (with targets: the nearest targets' cost, None when no target is reached; None for a
    full search).
    """

    def __init__(self, _no_public_constructor: typing.Never, /) -> None: ...
    def __repr__(self) -> str: ...
    @property
    def costs(self) -> numpy.typing.NDArray[numpy.int64] | None:
        """
        full search: each node's distance from the roots for FORWARD, or to the roots for REVERSE
        (read-only int64[N], UNREACHABLE for nodes not reached); None with targets
        """

    @property
    def dag(self) -> PredecessorDAG:
        """
        the PredecessorDAG of the result, in the edges' original orientation
        """

    @property
    def target_cost(self) -> int | None:
        """
        with targets: the cost of the nearest targets, None when no target is reached; None for a
        full search
        """

def apply_actions(
    ledger: FlowLedger,
    actions: collections.abc.Sequence[Place | Replace | Forward],
    *,
    mode: BatchMode = ...,
    outputs: OutputOptions | None = None,
) -> list[PlacementResult]:
    """
    Run the actions of several policies of ledger (a sequence of Place, Replace and Forward) as one
    transaction -> list[PlacementResult], one per action in input order, each the scalar call's fields
    for its policy. The whole list is validated first, then computed, then published once (one exchange
    of every policy's records, then their epochs): a failed call raises with the ledger and every
    policy unchanged (with the two limits in the module docstring). A policy appears at most once, and
    every policy belongs to ledger.
    mode ORDERED: each action runs in list order on a private working state exactly as its scalar call
    would on the state the earlier actions left (their carried load and holds); new records take
    FlowIds in list order. Earlier actions are served first: a reservation gives no protection. mode SIMULTANEOUS:
    Replace of ShortestPaths EQUAL or WEIGHTED (routes prepared on the view) and of Pinned, Forward of
    any policy with an installed program, every one LOSSY with a finite offer, forwarded at the same
    time on one view: the residual plus every participant's old carried load (outside owners and orphans
    stay; holds are never added). At a link that cannot carry its arrivals the loss is shared by highest
    averages over the demands (ties by policy creation order, never by list position), then over each
    demand's sources, then its tunnels; each participant's black holes are its own. Permuting the list
    changes nothing. outputs (OutputOptions; None: none) as for the scalar calls: residual and
    saturated_boundary_edges read the state after that action (ORDERED) or after the whole list
    (SIMULTANEOUS, one residual for every entry, each entry's boundary from its own sources and mask).
    Checked in this order: the arguments' types (TypeError: `ledger must be a FlowLedger`, `actions must
    be a sequence of Place, Replace or Forward (got <type>)`, `actions[<i>] must be a Place, Replace or
    Forward (got <type>)`, `mode must be a BatchMode (got <type>)`, `outputs must be an OutputOptions
    (got <type>)`); mode's value (`mode: <v> is not a value of BatchMode`); an empty list returns [];
    then each action in list order, all of its checks before the next action's: `list policy belongs
    to another ledger`, `duplicate policy in list`, SIMULTANEOUS's own rules (`SIMULTANEOUS does not
    accept Place`, `SIMULTANEOUS Replace requires ShortestPaths EQUAL or WEIGHTED, or Pinned`,
    `SIMULTANEOUS requires admission=LOSSY`, `SIMULTANEOUS requires finite offers`), then the scalar
    call's checks in its order. Only then is anything computed; an error of the computation (a
    loss-iteration error, a solver or forwarding size limit, the record identity limit counted for the
    whole list) publishes nothing.
    """

@typing.overload
def as_float(
    x: collections.abc.Sequence[typing.SupportsIndex]
    | numpy.typing.NDArray[numpy.integer],
    resolution: typing.SupportsFloat | typing.SupportsIndex,
) -> numpy.typing.NDArray[numpy.float64]:
    """
    Units to physical value(s): x * q, correctly rounded (exact below 2^53 units). An int, of
    any size, gives a float, an integer array a float64 array; a result beyond the finite
    doubles raises ValueError. x is an integer position (any object with __index__, never a bool;
    anything else is a TypeError); q must be a power of two in [2^-32, 2^32] (a float, a NumPy float
    or an integer).
    """

@typing.overload
def as_float(
    x: typing.SupportsIndex, resolution: typing.SupportsFloat | typing.SupportsIndex
) -> float: ...
def evaluate_demand(
    graph: Graph,
    sources: typing.SupportsIndex
    | collections.abc.Sequence[typing.SupportsIndex]
    | numpy.typing.NDArray[numpy.integer],
    destinations: typing.SupportsIndex
    | collections.abc.Sequence[typing.SupportsIndex]
    | numpy.typing.NDArray[numpy.integer],
    options: EvaluationOptions,
) -> PlacementResult:
    """
    The evaluation of a fresh FlowPolicy with options.config on a fresh ledger of graph at options.offer
    (None: +inf; DistributedTE and CentralizedTE with tunnel_count: evaluate(offer, design_demand=...))
    -> PlacementResult. LOSSLESS delivers up to a finite offer; LOSSY accounts for exactly that offer
    (ShortestPaths EQUAL/WEIGHTED LOSSY without an offer: the first-hop experiment). The options are checked as
    the EvaluationOptions docstring states.
    The exact maximum flow is CentralizedTE(source_allocation=OPTIMIZED); with source_allocation=EQUAL, the largest
    total whose round-robin source shares are all carried. A supplied frozen cache raises ShortestPathCacheMiss
    for an absent key, unchanged.
    A failed call raises; a supplied ShortestPathCache then behaves as the ShortestPathCache docstring states.
    """

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
) -> SensitivityResult:
    """
    The baseline evaluate_demand() and one SensitivityEntry per candidate edge. Candidates: with unrestricted
    CentralizedTE and the optimized source allocation, the edges that carry flow in the baseline (removing any
    other edge keeps the baseline feasible, at any offer); otherwise every allowed edge, whatever its capacity (an edge of capacity 0
    can carry a route, so deleting it can change the result, and with the equal source allocation a deletion that
    removes a source from the roster can raise the total: a negative reduction). delivery_reduction =
    baseline.call_delivered minus the call_delivered of evaluate_demand() without the edge, in signed
    units, never thresholded. Two
    shortcuts give exactly that value: unrestricted optimized CentralizedTE whose baseline is a maximum flow
    (no offer, or an offer above what it delivers) repairs a copy of the baseline's solve instead of solving
    again (at an offer the baseline delivers, each deletion is solved at the offer), and ShortestPaths
    gives reduction 0 without placing to an edge outside the static shortest-path support reachable from
    the roster sources (no route can use it). CentralizedTE with a support solves every deletion afresh; Pinned
    replays the same program under each deletion. workers as evaluate_pairs
    (NGRAPH_CORE_SENSITIVITY_THREADS). Only ShortestPaths reads the cache, and each placed deletion uses a
    new Eligibility uid, so a frozen cache raises ShortestPathCacheMiss there.
    A failed deletion fails the call (the lowest failing candidate's exception, no partial result); a
    supplied ShortestPathCache then behaves as the ShortestPathCache docstring states.
    """

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
) -> list[PlacementResult]:
    """
    evaluate_demand() for each (sources, destinations) pair, results in pair order. Pair i uses
    options.eligibility intersected with per_pair_eligibility[i] (None entries allowed; a given sequence has one entry per pair,
    also when empty). workers >= 1; None: NGRAPH_CORE_BATCH_THREADS or 1, and 1 when a cache is
    supplied. Everything is validated before a worker starts, the arguments first, then each pair
    completely in pair order, so the error of the lowest failing pair is raised. Pinned is rejected
    (`evaluate_pairs does not accept Pinned`: a program has fixed endpoints);
    an empty batch returns [].
    Caches: a supplied building cache needs workers == 1; a frozen cache is shared by every worker
    (pass workers explicitly), and its ShortestPathCacheMiss reaches the caller unchanged (the lowest failing
    pair's). An intersection that differs from both operands gets a new Eligibility uid, which a
    frozen cache cannot hold: pass the combined per-pair values with eligibility None, or a common value
    that contains them.
    A failed call raises; a supplied ShortestPathCache then behaves as the ShortestPathCache docstring states.
    """

def ksp(
    graph: Graph,
    source: typing.SupportsIndex,
    destination: typing.SupportsIndex,
    options: KspOptions,
) -> list[tuple[int, Path]]:
    """
    The k shortest simple source -> destination paths over the eligible edges (only those with residual >= 1 when
    a residual is given; an edge of capacity 0 is searched like any other) -> [(cost, Path)], each hop
    a single edge. Order: lowest cost; then, under a HIGHER_* tie_break only, the highest bottleneck
    (the smallest tie value along the path); then fewest hops; then the smallest EdgeId
    sequence. The result is the first k paths of that order
    within the inclusive cost ceiling. Validated in this order before any work (ValueError): source and
    destination in range, source == destination, k >= 1, the ceiling, the residual (length, then every
    entry), a foreign Eligibility. An excluded endpoint, no static route under a CostStretch ceiling, or
    no path gives [].
    """

def resolve_to_paths(
    graph: Graph,
    dag: PredecessorDAG,
    source: typing.SupportsIndex,
    destination: typing.SupportsIndex,
    split_parallel_edges: bool = False,
    max_paths: typing.SupportsIndex | None = None,
) -> list[Path]:
    """
    The source -> destination walks of a DAG -> [Path], in increasing (cost, hops, edge-id sequence):
    a hop costs its cheapest member and is keyed by its lowest member id. split_parallel_edges = False
    keeps every parallel member of a (u, v) pair in one hop (members ascending); True gives one path
    per member combination, sorted by that path's own (cost, hops, edge-id sequence). max_paths
    (>= 1) truncates; None lists every walk.
    ValueError: a DAG of another graph, source or destination out of range, source == destination,
    max_paths < 1. An endpoint outside the DAG, or no walk, gives [].
    """

def spf(
    graph: Graph,
    roots: typing.SupportsIndex
    | collections.abc.Sequence[typing.SupportsIndex]
    | numpy.typing.NDArray[numpy.integer],
    options: SpfOptions,
) -> SpfResult:
    """
    Shortest paths from the roots for FORWARD, or to the roots for REVERSE (one NodeId or a sequence;
    every root at distance 0), over the eligible edges (only those with residual >= min_residual when
    a residual is given; an edge of capacity 0 is searched like any other) -> SpfResult. Without
    targets every reachable node gets its distance: costs and the full DAG of tight edges (the edges
    on some shortest walk; one per (u, v) pair with parallel_edges = False). With targets the search
    stops at the nearest targets' cost:
    target_cost and the DAG of the walks between the roots and a nearest target (for REVERSE, from
    the target to a root, in the edges' own direction);
    with multipath = False, instead the canonical single path, the walk minimising (cost, -bottleneck
    under a HIGHER_* tie_break, hops, EdgeId sequence). Validated in this order before any work
    (ValueError): the roots (non-empty, ids in range), the targets (the same), a node in both sets,
    multipath = False without targets, a HIGHER_* tie_break with multipath and parallel_edges both
    True ('spf tie_break requires multipath=False or parallel_edges=False'), the residual (length,
    then every entry in [0, edge_capacity_units]),
    min_residual >= 1, a foreign Eligibility. Excluded roots and targets are skipped; all roots
    excluded gives an empty result.
    """

@typing.overload
def to_units(
    x: float | numpy.floating | numpy.integer,
    resolution: typing.SupportsFloat | typing.SupportsIndex,
) -> int:
    """
    Physical value(s) to units of the resolution q, rounding half up: an int for a scalar, an
    int64 array for an array. NaN, infinity, negative values, an integer too large for NumPy and
    results at or above 2^62 units raise ValueError; a bool (also as an entry), a string or an array
    of another dtype or rank, or any type but int, float and the NumPy numbers (a Decimal, a
    Fraction, an object known by __float__, __index__ or __array__) TypeError; an int or float
    subclass is read by its value. q must be a power of two in [2^-32, 2^32]: a float, a NumPy
    float or an integer (any other type is a TypeError; an integer beyond the float range a
    ValueError).
    """

@typing.overload
def to_units(
    x: collections.abc.Sequence[float | numpy.floating | numpy.integer]
    | numpy.typing.NDArray[numpy.floating | numpy.integer],
    resolution: typing.SupportsFloat | typing.SupportsIndex,
) -> numpy.typing.NDArray[numpy.int64]: ...

UNREACHABLE: int = 9223372036854775807
