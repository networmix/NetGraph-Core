# NetGraph-Core semantics

This document says exactly what NetGraph-Core computes: how amounts are counted, how edges are
numbered, which route a search returns, how each routing model spreads traffic, what a call reports
and when it fails, how the shortest-path cache behaves, and what the evaluators mean. It
is the full reference behind the API docstrings. For installation and a first example, see the
[README](../README.md#installation).

Each section starts with a short explanation in plain words and then gives the exact rules. When
the two seem to differ, the exact rules are the authority. Names are the Python spellings
(`import netgraph_core as ngc`); [the last section](#11-c-names) lists the C++ spellings. Names that
start with an underscore are private test and verification entries, not part of the API.

Contents:

1. [Units: the one numeric domain](#1-units-the-one-numeric-domain)
2. [The graph](#2-the-graph)
3. [Route values: PredecessorDAG and Path](#3-route-values-predecessordag-and-path)
4. [Search](#4-search)
5. [Route allocations and admission](#5-route-allocations-and-admission)
6. [The ledger](#6-the-ledger)
7. [Policies](#7-policies)
8. [Failures, threads and lifetimes](#8-failures-threads-and-lifetimes)
9. [The shortest-path cache](#9-the-shortest-path-cache)
10. [Evaluators](#10-evaluators)
11. [C++ names](#11-c-names)

## Terms

- **Unit, resolution `q`.** Traffic is counted in whole units; one unit is `q` of a physical
  quantity (for example `q = 2**-12` Gbit/s). The graph fixes `q` when it is built.
- **Amount.** A whole number of units: a capacity, a residual, a request, a delivery or a drop.
- **Residual.** The capacity of an edge that is still free, in units.
- **EdgeId.** The public number of an edge. It is not the edge's position in the input arrays
  (see [section 2](#2-the-graph)).
- **Tight edge.** An edge `u -> v` that lies on a shortest walk of a search. Searching from the
  roots, the distance of `u` plus the edge's cost equals the distance of `v`; searching toward
  them (reverse), where a distance is the distance to the roots, the distance of `v` plus the
  cost equals the distance of `u`.
- **Shortest-path DAG.** The set of all tight edges: every shortest route at once. A DAG (directed
  acyclic graph) is a directed graph with no directed cycle.
- **Parallel members.** Several edges between the same two nodes in the same direction.
- **Canonical path.** The one walk a search returns when a single path is asked for, chosen by
  the fixed order of [section 4](#the-canonical-single-path).
- **Bottleneck of a walk.** The smallest priority value (residual, capacity or headroom) among the
  walk's edges.
- **A flow** is one routed amount on one route, recorded in the ledger as a `FlowRecord`.
- **Policy.** The routing model of one demand, the traffic from a set of sources to a set of
  destinations (a `FlowPolicy`): it computes that demand's routes, forwards its traffic over them
  and removes it again.
- **Roster.** The sources of a policy that take part in placement, fixed when the policy is built.
- **Source allocation.** How a demand's total enters at its roster sources: in equal round-robin
  shares (`SourceAllocation.EQUAL`), highest-averages shares of declared weights (`WEIGHTED`),
  or as a controller chooses (`OPTIMIZED`).
- **Round robin `RR(x, k, i)`.** `x` units over `k` members: `x // k` each, and one more for the
  first `x % k` members.
- **Highest averages `H(x, w)`.** `x` units over members with integer weights `w`, given one at a
  time to the member with the largest `w[i] / (assigned[i] + 1)`, ties to the earlier member; a
  zero weight gets nothing.
- **Program.** A `ForwardingProgram`: the routes and integer weights of one demand and the rule
  that injects its traffic at the sources. `ShortestPaths`, `CentralizedTE`, T1 (`DistributedTE`)
  and T2 produce one; `Pinned` forwards a given one without computing anything
  ([section 7](#forwarding-programs)).
- **Entry.** One route of a program: a **tunnel** (one path, one edge per hop) or a **split
  route** (a DAG whose routers split by fixed rules). Each entry has an integer identifier, its
  `entry_id`, and a weight: its share among its head's entries.
- **Head.** A roster source, as the start of its entries.
- **T1.** `DistributedTE(TunnelCount(N))`: every head sets up its own N strict-bandwidth tunnels by
  CSPF (constrained shortest path first: a shortest path among the links with enough unreserved
  bandwidth), with no controller ([section 7](#t1-tunnels-set-up-by-their-heads)).
- **T2.** `CentralizedTE(source_allocation=..., tunnel_count=TunnelCount(N))`: a central
  controller lays out N single-path tunnels for every live head with one common whole-unit planned
  bandwidth ([section 7](#t2-the-controllers-equal-tunnels)).
- **Hold.** Bandwidth an up T1 tunnel reserves on every link of its path. Holds are the ledger's
  second projection (`reserved`, `unreserved`), kept apart from the carried traffic (`carried`,
  `residual`); they never police or queue traffic.
- **Design demand.** The volume T1 and T2 size their tunnels for (`design_demand`), separately
  from the offer they forward. T2's planned bandwidth is metadata, never a hold.
- **Split context.** The traffic a split route's routers split: only that route's own share. A
  `ShortestPaths` program and a `from_dag` program are one context pooled over their sources;
  every split route given to `from_routes` is its own.
- **Black hole.** A node or edge that fails silently: forwarding goes on unchanged and traffic
  that reaches it disappears there (stale forwarding).
- **Controller.** A routing that solves a flow problem for the whole demand: `CentralizedTE`
  without a tunnel count, and `ShortestPaths` with `allocation=OPTIMIZED`.
- **Offer.** The traffic a call asks a policy to carry (`place(offer)`, `replace(offer)`, ...).
- **Epoch.** The period from a policy's construction, a `remove()`, a `replace()` or a `forward()`
  to the next of these; the offers of `place` add up within it.
- **Lossless / Lossy.** Admission control: Lossless places only what fits end to end; Lossy
  forwards hop by hop and drops what does not fit on each edge.
- **Static.** Depending only on the graph, the costs and the Eligibility, not on traffic already
  placed (a static walk, the static shortest cost).
- **Route choice.** The tie rule (`RouteTieBreak`), `multipath` and `parallel_edges` of a model or
  a search ([section 4](#the-route-choice)): cost always comes first, and the tie rule only
  orders equal-cost choices.
- **Member graph.** The edges a search may use: eligible (and, when a residual is supplied, at
  least the search's residual threshold), whatever their capacity, and one member per pair with
  `parallel_edges=False`.
- **Route allocation.** The rule that spreads one flow over its route (`RouteAllocation`):
  `OPTIMIZED`, `WEIGHTED` or `EQUAL`.
- **Route of a record.** The edges one record may use, as a `PredecessorDAG` (`FlowRecord.route`).
- **View.** The capacity a solve may use: the free capacity plus the policy's own allocation.
- **Minimum cut.** A set of edges whose removal separates the sources from the destinations,
  with the smallest total capacity; that capacity equals the maximum flow.
- **Marginal cost curve.** `(cost, units)` pairs in ascending cost: delivering `units` more units
  raises the minimum total cost of the flow by `cost` per unit. Delivering more can reroute flow
  already placed, so the pairs are not the costs of the paths in the final allocation.

The examples below run in sequence and share these helpers. Each output block is the printed
output of the example above it.

```python
import netgraph_core as ngc

def graph(num_nodes, edges, q=1.0):
    """edges: (source, destination, capacity units, cost) tuples in input order."""
    sources, destinations, capacities, costs = zip(*edges)
    return ngc.Graph(num_nodes, list(sources), list(destinations), list(capacities), list(costs), q)

def config(routing, admission=ngc.Admission.LOSSLESS):
    """A configuration: the routing model and the admission."""
    return ngc.FlowPolicyConfig(routing=routing, admission=admission)
```

## 1. Units: the one numeric domain

**In plain words.** Every amount of traffic is a whole number of units, and the graph says how
big a unit is. Physical numbers enter through `to_units` and leave through `as_float`; nothing
in between is a float, so totals add up exactly and nothing is lost to rounding. Costs are a
separate kind of number: positive integers that are never converted.

**Exact rules.**

- **Resolution.** A graph's resolution `q` is a power of two in `[2^-32, 2^32]`, given
  explicitly when the graph is built; there is no default. The resolution is the precision of
  every result: there is no finer grid. A model that needs finer splits uses a smaller `q`.
- **Amounts.** Capacities, residuals, allocations, requests, deliveries, offers and drops are
  integers in units of `q`. In Python they are `int` (exact, unbounded); per-edge arrays are
  NumPy `int64[E]`. Cumulative counters that are sums over many calls (a Lossy record's
  `offered` and its per-edge `edge_drops`) may exceed 64 bits; Python returns them as exact `int`.
  Cumulative counters must remain below `2^127`. A call that would reach that limit raises
  `ValueError('<counter> limit exceeded')` and publishes nothing; the cumulative request is
  checked before computation.
- **Internal limits.** A forwarding's network numbers its positions in 32 bits: its forwarded
  program entries, occurrences (every link visit of every entry: a tunnel's path length, a split
  route's member links), split-route routers, LOSSY events (the links that can drop plus the
  routers) and their dependency arcs, checked in that order, are each at most `2^31 - 1`; a call
  beyond raises `ValueError('forwarding <counter> limit exceeded: <n> > 2147483647')` (`<counter>`
  one of `entry`, `occurrence`, `router`, `event`, `dependency`; `<n>` the exact count). The
  entry, occurrence and router limits are checked after validation and before the network is
  allocated; the event and dependency limits, known only after the loads of LOSSY's first
  forward (each demand's traffic stopped at its own black holes), before the
  event structures are allocated. Either refusal changes nothing. A list of several demands
  ([several demands at once](#several-demands-at-once)) counts its one joint network, every
  participant together. Identities are
  never reused: a ledger issues at most `2^63 - 1` policy identifiers and `2^63 - 1` FlowId
  counters, and the process at most `2^63 - 1` identities shared by graphs, ledgers and masks (0,
  the all-allowed mask's, is never issued); the next one raises `RuntimeError('policy identity
  limit exceeded')`, `RuntimeError('record identity limit exceeded')` (a publication's new records
  counted together, a list's included) or `RuntimeError('object identity limit exceeded')` before anything changes.
- **Domain.** Every edge capacity, and the sum of all capacities, is below `2^62` units. Edge
  costs are integers of at least 1 and their sum is below `2^62`, so no path cost reaches the
  sentinel `UNREACHABLE` (the largest `int64`). An offer is an integer in `[0, 2^62)` units or
  positive infinity (`float('inf')`, "as much as possible"). Everything outside the domain is a
  `ValueError` raised before anything changes.
- **`to_units(x, q)`** reads each physical input as its binary64 value `b`: Python floats and
  NumPy float64 keep their exact value, float32 widens exactly, and integers become the nearest
  double, ties to even (`2^53 + 1` becomes `2^53`, `2^53 + 3` becomes `2^53 + 4`). It returns
  `floor(b/q + 1/2)`, evaluated exactly on rational numbers, not by floating addition of `1/2`.
  A scalar or 0-d NumPy array gives an `int`; a Python sequence or 1-d NumPy array gives an
  `int64` array. Negative zero gives 0. NaN, infinity, negative values and results at or above
  `2^62` units raise `ValueError`. This is the only entry for physical numbers.
- **Binary values and precision.** A float is never read through its printed decimal; conversion
  cannot recover a source literal. `2^53` is the limit of consecutive integers among doubles,
  not a size below which decimal fractions are exact. A decimal is exactly a double only if
  its reduced denominator is a power of two and it fits binary64; even `0.1` is rounded on
  parsing. A decimal and its parsed double give the same units exactly when both lie in the
  same cell `[(n - 1/2)*q, (n + 1/2)*q)`. No magnitude or unit count guarantees this: a parse
  error can cross or land on a boundary at an ordinary size. At `q = 2^-21`, `5501195.1796`
  denotes the exact binary tie `11536842473288.5` units, giving `11536842473289`; the typed
  decimal would be `11536842473288.4992` units. At `q = 2^-18`, the computed `400 + 2^-19`
  is the tie `104857600.5` and gives `104857601`, although its printed decimal
  `400.00000190734863` is just below the tie. `np.float32(0.1)` is exactly
  `13421773/134217728`, giving `429496736` units at `q = 2^-32`; decimal `0.1` would give
  `429496730`. `CostStretch` deliberately uses the printed decimal instead: a stretch is a
  ratio a person writes (`1.2` means `6/5`), compared exactly with ratios of integer costs.
  A physical amount is the measured or computed double the caller holds.
- **Ties and exact scaling.** The ties are exactly the accepted `b = (2*n + 1)*q/2` for integers
  `0 <= n < 2^52`; every tie gives `n + 1`, never the even neighbour. At `q = 1`, `0.5`
  gives 1, `2.5` gives 3, `3.5` gives 4, and the last tie `2^52 - 1/2` gives `2^52`.
  For `q = 2^k`, scaling `y = b*2^-k` is exact while normal; only below the normal range
  (`y < 2^-1022`) may it lose bits, where the answer is 0 either way. Overflow is refused by
  the `2^62` bound. Below `2^52`, `y - floor(y)` is exact and compared with `1/2`; at and
  above `2^52` every double is integral, so `y` is the result. These double operations thus
  introduce no error beyond the specified quantization.
- **`as_float(u, q)`** returns the correctly rounded binary64 value of the exact product `u*q`
  (nearest, ties to even). Here `u` may be any signed, unbounded integer, unlike graph and
  flow amounts. An integer scalar gives a `float`, an integer vector a `float64` array.
  An int or NumPy integer, a subclass included, is never read as a sequence.
  The result is finite exactly when `abs(u*q) < 2^1024 - 2^970`; otherwise it raises
  `ValueError('as_float: <u> units times the resolution exceed the largest finite double')`.
  Every `abs(u) <= 2^53`, including the endpoint, converts exactly. Beyond `2^53`, an integer
  of binary length `L` converts exactly iff it is divisible by `2^(L - 53)`; `2^53 + 1`
  gives `2^53*q`, and `2^53 + 3` gives `(2^53 + 4)*q`. Costs are never converted.
- **Round trip.** `to_units(as_float(u, q), q) == u` for every `u` in `[0, 2^53]` at every
  resolution, and for larger `u` in `[0, 2^62)` exactly when `u` is a double. Otherwise the
  result is the units of the rounded value. The first failure, `2^53 + 1`, comes back as
  `2^53`; `2^62 - 1` becomes `2^62*q`, which `to_units` refuses.
- **Error and sums.** Every successful `to_units` result `u` is itself a double, so `as_float`
  adds no error: `e = as_float(to_units(b, q), q) - b` satisfies `-q/2 < e <= q/2`, with
  equality at the upper end exactly at ties. This bound concerns `b`. Against an argument
  `a` changed by promotion (a large integer or a wider NumPy float where available), the
  error is `e + (b - a)`: `2^53 + 1` at `q = 1` has error 0 against `b`, but -1 against
  the original integer. There is no exact-integer path for physical inputs.
  Rounding `N` values separately gives total error in `(-N*q/2, N*q/2]`. Ties contribute
  `+q/2` each; other errors go either way, so a sum has no promised sign. At `q = 1`,
  `[0.5, 0.5]` gives `[1, 1]` (2 units for a total of 1), while `[0.25, 0.25]` gives
  `[0, 0]` (0 units for a total of 0.5). When non-tie fractional parts are symmetric about
  `1/2`, their mean error cancels and mean error per value is `(q/2)*P(tie)`; a continuous
  uniform fractional part has none. The Core assumes no distribution. Round a demand once,
  then split whole units to preserve its total; rounding `1.0` once gives 1 unit.
  Choose one `q` before graph construction and use the constructed graph's `resolution()` afterwards. Choosing
  `q` and reporting rounding belong to the caller (NetGraph's units report); Core converts
  and validates, without a selector or rounding report.
- **Routability.** An edge with less than 1 unit of capacity carries nothing, but it stays in
  routing: searches, routes and rosters treat it like any other edge, and only a supplied
  residual (below the search's threshold) or an Eligibility removes an edge from a search. An
  equal split still sends such an edge its share, which it cannot carry
  ([section 5](#5-route-allocations-and-admission)). A physical offer below `q / 2` converts to
  0 units. An offer of 0 units is accepted after validation and does nothing (its result is
  `COMPLETE`); it is not a removal: `remove()` is.
- **Splitting units.** When a node forwards `I` units over `k` member edges:
  - **round robin** (`EQUAL`): every member gets `I // k` units and the first `I % k`
    members, in ascending EdgeId, get one more;
  - **highest averages** (`WEIGHTED`, the D'Hondt rule): the units are handed out one at a time,
    each to the member with the largest `w / (a + 1)` (its weight `w` over its current count
    `a` plus one), ties to the lower EdgeId. Equal weights give round robin.

  Both rules are exact and never lose a unit: a share that rounds down to zero leaves its unit to
  a sibling. Both are monotone: when `I` grows, no member's share shrinks. For one split into `k`
  members, each member's absolute rounding error relative to that split's actual integer input and
  effective weights is below `k` units (`k * q` physically), and below one unit for equal positive
  weights. This is a local bound, not an end-to-end delivery bound, and highest averages is not
  the closest proportional rounding: weights `[100, 1, 1, 1]` and 50 units give `[50, 0, 0, 0]`
  (total absolute error `300/103` against the exact shares `[5000/103, 50/103, 50/103, 50/103]`,
  where `[49, 1, 0, 0]` would have `200/103`). With all weights zero, no unit is assigned and a
  proportional share is undefined.
- **Conservation.** At every node, inflow plus supply equals outflow plus delivery plus drops,
  exactly. The ledger's totals are exactly the sums of its records.

```python
print(ngc.to_units([1.5, 2.5, 0.4, 0.5], 1.0).tolist())   # half up
print(ngc.to_units(10.0, 2.0**-12), ngc.as_float(40960, 2.0**-12))
```

```text
[2, 3, 0, 1]
40960 10.0
```

## 2. The graph

**In plain words.** A `Graph` is built once from arrays and never changes. Every edge gets a
public number, its EdgeId, by sorting the edges by source, cost and destination; only edges
that agree on all three keep their input order among themselves. A map gives back each edge's
input position. All per-edge arrays the library reads or returns are in EdgeId order, and every
tie rule of the library ("the lowest EdgeId") refers to this order. A caller who keeps per-edge
values in input order (explicit weights, a residual, exclusions) maps them through
`input_positions()` first.

**Exact rules.**

- **Construction.** `Graph(num_nodes, edge_sources, edge_destinations, edge_capacity_units,
  edge_costs, resolution)` builds an immutable directed multigraph: parallel edges and self-loops
  are allowed. The arrays are in input order (the caller's). Capacities are already in units; use
  `to_units` for physical values.
- **Validation order.** When an input has several defects, the first in this order is reported (a
  wrong type as a `TypeError`, a value outside its domain as a `ValueError`): the resolution (a
  bool is a `TypeError`); the capacity values (each an integer in `[0, 2^62)`; floats and booleans
  are rejected with `TypeError`, not rounded); the array shapes; a zero or negative cost; the remaining
  bounds (node ids in range, then the capacity and cost sums below `2^62`, both added up edge by
  edge in input order: the first edge at which a running sum reaches `2^62` reports that sum, the
  capacity sum before the cost sum at the same edge). Integers are read
  without narrowing, so the reported error does not depend on how the same input is represented.
- **Size.** Node ids and EdgeIds are `int32`: a graph has at most `INT32_MAX` (`2^31 - 1`) nodes
  and at most as many edges. The constructor accepts every graph within these bounds, and the
  searches and route values work on all of them. The flow core (`CentralizedTE` routing, the
  `OPTIMIZED` allocation and the evaluators that use them) numbers a solve's network in `int32` as
  well: it adds a virtual
  source and a virtual sink, and gives every edge and every virtual arc (one per source, one per
  destination) two slots, one to add flow and one to cancel it. A solve therefore raises
  `ValueError`, before the solver sizes anything and with nothing changed, on a network of more
  than `INT32_MAX - 2` nodes, or with more than `2^30 - 1` edges, sources and destinations
  together. A `CentralizedTE` solve counts the graph's nodes, the eligible edges, the roster
  sources and the eligible destinations; a `ShortestPaths` `OPTIMIZED` solve counts the nodes of
  its route's edges and endpoints, those edges, its sources and its destinations. Only graphs of
  billions of nodes or about a billion edges reach these limits.
- **EdgeId.** The EdgeId of an edge is its rank in the stable order `(source, cost, destination,
  input position)`. `input_positions()[e]` is the input position of edge `e`; the inverse map is
  the caller's to build (one pass over `E`).
- **EdgeId order is the coordinate of every per-edge value.** Residual arrays, explicit weights,
  Eligibility exclusions, DAG and Path edges, flows, records, `saturated_boundary_edges`, `min_cut`
  and the entries of `evaluate_edge_failures` all use EdgeIds. A caller that keeps per-edge data in input order translates it first: an
  in-range id in the wrong coordinate is silently a different edge, and no check can catch that.
- **Layout guarantees.** The out-edges of node `u` are exactly the ids
  `outgoing_offsets()[u]` to `outgoing_offsets()[u + 1] - 1`, in `(cost, destination, input)`
  order, and `edge_destinations()[e]` is the head of `e`. The in-edges of node `v` are
  `incoming_edges()[incoming_offsets()[v]:incoming_offsets()[v + 1]]`, as EdgeIds in
  `(cost, source, input)` order. Both offset arrays have `N + 1` entries, starting at 0 and ending
  at `E`, also for an empty graph.
- **Accessors.** `edge_capacity_units()`, `edge_costs()` (`int64[E]`), `edge_sources()`,
  `edge_destinations()`, `input_positions()`, `outgoing_offsets()`, `incoming_offsets()` and
  `incoming_edges()` (`int32`) are read-only NumPy views that keep the graph alive;
  `max_edge_cost()` and `edge_cost_gcd()` are 0 for a graph without edges; `num_nodes()`,
  `num_edges()` and `resolution()` read the sizes; `uid()` is a process-unique number, never
  reused.
- **Ids belong to one graph.** Building a graph from different arrays, or with another release of
  the library, may renumber the edges. Every value bound to a graph (Eligibility, PredecessorDAG,
  Path, program, ledger, cache) is rejected with another graph.

```python
# Input order: edge 0 = 1->2, edge 1 = 0->2 (cost 5), edge 2 = 0->1.
g = ngc.Graph(3, [1, 0, 0], [2, 2, 1], [7, 8, 9], [1, 5, 1], 1.0)
print("sources", g.edge_sources().tolist(), "destinations", g.edge_destinations().tolist(), "costs", g.edge_costs().tolist())
print("input_positions", g.input_positions().tolist(), "outgoing_offsets", g.outgoing_offsets().tolist())
```

```text
sources [0, 0, 1] destinations [1, 2, 2] costs [1, 5, 1]
input_positions [2, 1, 0] outgoing_offsets [0, 2, 3, 3]
```

### Eligibility

**In plain words.** An `Eligibility` says which nodes and edges of one graph may be used, for
example to model a failure or a maintenance window. It is built from lists of exclusions and
never changes.

**Exact rules.**

- `graph.eligibility(excluded_nodes=(), excluded_edges=())` builds the value. Duplicates and
  unsorted lists are accepted; out-of-range ids raise `ValueError`. An edge is allowed when it is
  not excluded and both its endpoints are allowed. `node_allowed()` and `edge_allowed()` are
  read-only boolean views.
- Empty lists give the all-allowed value, whose `uid()` is 0. Every other value has a fresh uid;
  equal lists built twice give two different uids. Wherever an Eligibility is optional, `None`
  means everything is allowed.
- `Eligibility.intersect(a, b)` allows what both allow. When the result equals one operand, that
  operand is returned (same uid); otherwise the result is a new value with a new uid. Both
  operands must belong to the same graph.
- Uids matter for the [DAG cache](#9-the-shortest-path-cache): its key contains the Eligibility uid, so
  build each state once and reuse the value.
- An excluded endpoint takes no part in a search or a policy: an excluded root or target in a
  set is skipped, and an excluded single endpoint gives an empty result.

## 3. Route values: PredecessorDAG and Path

**In plain words.** A route is either a DAG of predecessor edges (the shape a search returns:
possibly many equal-cost walks) or a `Path` (one explicit walk, each hop possibly a bundle of
parallel edges). Both are immutable, bound to their graph and checked once when built.

**Exact rules.**

- **`PredecessorDAG`.** For node `v`, its entries are `parents()[k]` (the parent node) and
  `via_edges()[k]` (the edge from that parent) for `k` from `parent_offsets()[v]` to
  `parent_offsets()[v + 1] - 1`, in strictly ascending EdgeId. This is the one canonical form:
  two DAGs of the same graph are equal exactly when their entries are equal; DAGs of different
  graphs are never equal. `parents()` and `via_edges()` are read-only views; `parent_offsets()`
  is built on demand.
- **`PredecessorDAG.from_csr(graph, parent_offsets, parents, via_edges)`** validates, before
  reading any structure: `parent_offsets` has `N + 1` entries, starts at 0, never decreases and
  ends at `len(parents) == len(via_edges)`; then every id is in range, every via edge joins its
  parent to its
  node, no `(node, via)` pair repeats, the DAG has no cycle and every node's entries are in
  ascending via EdgeId. A wrong type raises `TypeError`; structural and range defects raise
  `ValueError`, in the order above.
- **`Path(graph, hops)`.** `hops` is a list of hops, each a non-empty list of EdgeIds that share
  one `(u, v)` pair (parallel members). Validation, in order: every member id in range, no
  member twice, at least one hop, no empty hop, the members of a hop share `(u, v)`, consecutive
  hops connect, no node is visited twice. The members keep the caller's order. `cost()` is the
  sum over the hops of each hop's cheapest member; `source()`, `destination()` and `hops()` read
  the walk;
  `to_dag()` gives the equivalent `PredecessorDAG`. A single-member walk is `[[e1], [e2], ...]`.
- **`resolve_to_paths(graph, dag, source, destination, split_parallel_edges=False,
  max_paths=None)`** lists the `source -> destination` walks of a DAG as `Path` values in
  increasing `(cost, hop count, edge-id
  sequence)`, where a hop costs its cheapest member and is represented in the key by its lowest
  member id. With `split_parallel_edges=False` the parallel members of a pair stay together in
  one hop (ascending); with `True` every member combination is a separate path, and the paths are
  sorted by each path's own `(cost, hop count, edge-id sequence)`. `max_paths` (at least 1) keeps
  the first paths of that order. `source == destination`, a DAG of another graph or
  `max_paths < 1` raise `ValueError`; an endpoint outside the DAG, or no walk, gives `[]`.

## 4. Search

**In plain words.** `spf` finds shortest paths from a set of roots, forward or in reverse,
either to every node or only to the nearest of a set of targets. It returns the shortest-path
DAG, which contains every shortest route at once, or, on request, one canonical path chosen by a
fixed rule. `ksp` lists the K shortest simple paths in a fixed order. All searches see every
eligible edge, whatever its capacity: an edge of capacity 0 is a route like any other. Only a
supplied residual filters edges by capacity.

### The route choice

The searches and the models that choose routes (`ShortestPaths`, T1's CSPF) share one rule set.
Cost always comes first; the rest only orders equal-cost choices.

- **`tie_break: RouteTieBreak`**, one of four values: `FEWER_HOPS` (the default: cost, then fewer
  hops, then the smaller EdgeId sequence), `HIGHER_CAPACITY`, `HIGHER_AVAILABLE_CAPACITY` and
  `HIGHER_HEADROOM` (cost, then the larger bottleneck of the named edge score, then fewer hops,
  then the smaller EdgeId sequence).
- **Edge scores.** `HIGHER_CAPACITY` reads the static capacity; `HIGHER_AVAILABLE_CAPACITY` the
  availability of the current search; `HIGHER_HEADROOM` the availability divided by the static
  capacity, compared exactly as a fraction (an edge of capacity 0 scores 0). The availability is
  the supplied `residual` for `spf` and `ksp`, the ledger's residual for `ShortestPaths`, the
  unreserved view for T1's CSPF, and the static capacity where no view is supplied (then
  `HIGHER_HEADROOM` reads 1 on an edge of positive capacity and 0 on an edge of capacity 0).
- **`multipath: bool`** (`ShortestPaths` and `spf` only; default `True`): `True` keeps every
  equal-cost route (a DAG), `False` the one canonical path. `None` is not a value. T1 always
  builds single-path tunnels; `ksp` always returns individual paths.
- **`parallel_edges: bool`** (default `True`): `True` keeps every eligible parallel member of a
  `(u, v)` pair; `False` keeps one per pair by the member rule.
- **The member rule** (`parallel_edges=False`). Of the eligible parallel members of a pair, keep
  the one with the lowest cost; among those, under a `HIGHER_*` rule, the one with the highest
  score; then the lowest EdgeId. This choice happens before any route is searched, so it can
  change which walks exist.
- **No silent no-op.** With `multipath=True` and `parallel_edges=True` a tie rule has nothing to
  choose, so a `HIGHER_*` value there is rejected: `spf tie_break requires multipath=False or
  parallel_edges=False` (and `ShortestPaths tie_break requires multipath=False or
  parallel_edges=False` for a model). `FEWER_HOPS` is always valid and never prunes the all-routes
  DAG by hop count. T1 and `ksp` always order paths and need no such rule.
- **The tie rule never supplies weights** (`ShortestPaths` `WEIGHTED` reads a `WeightBasis`) and
  never enables a controller objective (T2's headroom objective is its own option).
- **Static eligibility** (always): the Eligibility. Capacity is never a condition: an edge of
  capacity 0 is searched like any other.
- **Residual gating** (only when a `residual` array is supplied): `residual` is an `int64` array
  in EdgeId order with every entry in `[0, edge_capacity_units[e]]`. It excludes edges whose
  residual is below the threshold (`min_residual` units for `spf`, 1 unit for `ksp`) and supplies
  the availability of the tie scores. It is read during the call only and must not change during
  it.

### Shortest paths

`spf(graph, roots, SpfOptions(direction, *, targets=None, tie_break=FEWER_HOPS, multipath=True,
parallel_edges=True, residual=None, min_residual=1, eligibility=None))` returns an `SpfResult` with `costs`, `dag` and `target_cost`.

- **Roots and targets** are one node id or a sequence. Every root starts at distance 0.
  Duplicates are removed. Validated in this order before any work (`ValueError`): the roots
  (non-empty, ids in range), the targets (the same), a node in both sets, `multipath=False`
  without targets (`multipath = false requires targets (the single path joins a root to a
  target)`), the ineffective tie rule above, the residual (length, then every entry),
  `min_residual` (an integer of at least 1), a foreign Eligibility. Excluded roots and targets are
  skipped; when every root is excluded the result is empty, not an error.
- **Direction.** `FORWARD` searches from the roots along the edges. `REVERSE` searches toward the
  roots against the edges, but the DAG it returns keeps the edges' own orientation: its walks
  run from a target toward a root. A forward and a reverse search with the roles of roots and
  targets exchanged describe the same shortest walks.
- **Without targets** every reachable node gets its final distance: `costs` is a read-only `int64[N]` array
  with `UNREACHABLE` for the nodes not reached, `dag` holds every tight edge, `target_cost` is
  `None`.
- **With targets** the search stops at the cost of the nearest targets (anycast: the nearest of
  a set). `target_cost` is that cost (`None` when no target is reached), `dag` keeps only the
  walks that end in a target at that cost, and `costs` is `None`.
- **The DAG is unique.** With `multipath=True`, costs are positive, so the member graph, the
  roots, the direction and the targets (or their absence) determine one DAG. The member graph
  changes through `parallel_edges`, the member rule, static eligibility and residual gating; with
  the member graph fixed, the tie rule, the capacity and residual values and the C++ queue choice
  do not change the DAG. A forward search with targets and a reverse search with the
  roots and the targets exchanged return the same DAG. Without targets the two directions can
  return different DAGs, also when the reverse search is rooted at a node the forward search
  reaches: each keeps the edges tight in its own direction.

### The canonical single path

With `multipath=False` (targets required), `spf` returns the one walk from a root to a target
that is first in this order:

1. the lowest **cost**;
2. under a `HIGHER_*` tie rule only: the highest **bottleneck**, the smallest edge score along
   the walk (`FEWER_HOPS` skips this step);
3. the fewest **hops**;
4. the lexicographically smallest **sequence of EdgeIds** along the walk.

The order is total, so the walk is unique. The same walk is returned by a forward search, by a
reverse search with the roles exchanged, and by `ksp` with `k=1`. Note the order of steps 3 and
4: a walk with fewer hops wins even when its first EdgeId is larger. These steps compare whole
walks, never their first edge alone. The member rule above is a separate, earlier choice: with
`parallel_edges=False` each pair is first reduced to one member, and the canonical path is then
chosen among the walks that remain.

```python
F = ngc.Direction.FORWARD
# Two equal-cost walks 0->1->3 (capacity 1) and 0->2->3 (capacity 5).
g = graph(4, [(0, 1, 1, 1), (0, 2, 5, 1), (1, 3, 1, 1), (2, 3, 5, 1)])
for tie in (ngc.RouteTieBreak.FEWER_HOPS, ngc.RouteTieBreak.HIGHER_CAPACITY):
    r = ngc.spf(g, 0, ngc.SpfOptions(F, targets=3, tie_break=tie, multipath=False))
    print(tie.name, r.target_cost, r.dag.via_edges().tolist())

# Equal cost 2: the direct edge (EdgeId 1) against the two-hop walk [0, 2].
g = graph(3, [(0, 1, 1, 1), (1, 2, 1, 1), (0, 2, 1, 2)])
r = ngc.spf(g, 0, ngc.SpfOptions(F, targets=2, multipath=False))
print("fewest hops:", [p.hops() for p in ngc.resolve_to_paths(g, r.dag, 0, 2)])
full = ngc.spf(g, 0, ngc.SpfOptions(F))
print("costs", full.costs.tolist(), "dag", full.dag.via_edges().tolist())
```

```text
FEWER_HOPS 2 [0, 2]
HIGHER_CAPACITY 2 [1, 3]
fewest hops: [[[1]]]
costs [0, 1, 2] dag [0, 1, 2]
```

### K shortest paths

`ksp(graph, source, destination, KspOptions(k, *, tie_break=FEWER_HOPS, parallel_edges=True,
cost_ceiling=None, residual=None, eligibility=None))`
returns a list of `(cost, Path)` pairs; every hop of every Path has one member.

- **Order.** The simple `source -> destination` paths of the member graph in increasing `(cost,
  -bottleneck under a HIGHER_* rule, hop count, edge-id sequence)`: the canonical order above
  applied to whole paths. The result is the first `k` paths of that order that are within the
  ceiling.
- **Parallel members.** `parallel_edges=True` makes every eligible member a separate path;
  `parallel_edges=False` keeps one member per pair by the member rule. `ksp` has no `multipath`.
- **Ceiling.** `cost_ceiling` is `None`, `AbsoluteCost(c)` (cost at most `c`; `AbsoluteCost(0)` is
  a valid budget) or `CostStretch(f)` (cost at most `f` times the static shortest cost of the
  member graph under the Eligibility, residual ignored, edges of capacity 0 included). Both are
  inclusive and compared exactly ([CostStretch](#cost-ceilings-and-the-exact-stretch)). With a
  stretch, no static route gives `[]`.
- **Residual gating** at 1 unit when a residual is supplied.
- Validated in this order before any work (`ValueError` unless stated): `source` and
  `destination` in range, `source == destination` (`node <v> is in both endpoint sets: source and
  destination must not overlap`), `k >= 1` (`k must be >= 1 (got <k>)`), the ceiling (a value of
  another type: `TypeError("cost_ceiling must be AbsoluteCost, CostStretch or None")`; then
  `cost_ceiling AbsoluteCost cost must be >= 0 (got <c>)` or `cost_ceiling CostStretch factor must
  be >= 1 (got <n>/<d>)`), the residual (length, then every entry), a foreign Eligibility.
  `KspOptions` keeps the ceiling as given; `ksp` checks it in this place. An excluded endpoint or
  no path gives `[]`.

Internally, each pending alternative remembers the spur index that generated it. When that
path is accepted, Lawler's partition needs searches only at that index and later ones: earlier
prefix regions already have their pending minimum. The generating index is included in the
expansion, and all preceding edges still contribute their cost, exact bottleneck and excluded
vertices. The index never enters the path key. Each spur returns the minimum whole prefixed
path under that same key, so this partition preserves the complete ordered result of full Yen.

```python
g = graph(3, [(0, 1, 3, 1), (0, 1, 5, 1), (1, 2, 9, 1)])   # parallel members 0 and 1
higher = ngc.RouteTieBreak.HIGHER_CAPACITY
for options in (ngc.KspOptions(5), ngc.KspOptions(5, tie_break=higher),
                ngc.KspOptions(5, tie_break=higher, parallel_edges=False)):
    print([(c, p.hops()) for c, p in ngc.ksp(g, 0, 2, options)])
r = ngc.spf(g, 0, ngc.SpfOptions(F, targets=2))
print([p.hops() for p in ngc.resolve_to_paths(g, r.dag, 0, 2)],
      [p.hops() for p in ngc.resolve_to_paths(g, r.dag, 0, 2, split_parallel_edges=True)])
```

```text
[(2, [[0], [2]]), (2, [[1], [2]])]
[(2, [[1], [2]]), (2, [[0], [2]])]
[(2, [[1], [2]])]
[[[0, 1], [2]]] [[[0], [2]], [[1], [2]]]
```

### Cost ceilings and the exact stretch

`AbsoluteCost(cost)` is an inclusive cost bound, an integer (`AbsoluteCost cost must be an integer
(got <type>)`, `TypeError`); a negative cost constructs and is refused where it is used.
`CostStretch(factor)` is one exact reduced fraction `n/d` with `d >= 1`, converted in this order:

1. a boolean (Python `bool` or `numpy.bool_`) raises `TypeError`;
2. a `fractions.Fraction` is taken as it is;
3. a float (a Python `float` or any NumPy floating value) is `Fraction(str(f))`, the shortest
   decimal that value's own type prints: `CostStretch(1.2)`, `CostStretch(numpy.float32(1.2))` and
   `CostStretch(Fraction(6, 5))` are all exactly 6/5, never the binary value of the float. A NaN
   or an infinity raises `ValueError` (`CostStretch factor must be finite (got <value>)`);
4. an integer (an `int`, a NumPy integer, any object with `__index__`) is `n/1`;
5. anything else raises `TypeError` (`CostStretch factor must be a Fraction, a float or an integer
   (got <type>)`).

Then the reduced components are checked: `CostStretch fraction numerator must be in [-2^63, 2^63)`,
then `CostStretch fraction denominator must be in [1, 2^63)` (`Fraction(2**64, 2**62)` reduces to
4/1 and is accepted). `factor` returns the `Fraction`; equality, `hash` and `repr` follow it. A
factor below 1 is a valid value, refused where the ceiling is used. A path of cost `c` passes
`CostStretch(n/d)` against a baseline `c0` exactly when `c * d <= n * c0`, both products computed
exactly in 128 bits. For `ksp` the baseline is the static shortest cost; for T1 it is one baseline
for the whole pool of heads ([T1](#t1-tunnels-set-up-by-their-heads)).

```python
from fractions import Fraction
print(ngc.CostStretch(1.2).factor, ngc.CostStretch(1.2) == ngc.CostStretch(Fraction(6, 5)),
      ngc.CostStretch(Fraction(2**64, 2**62)))
# Two links 0 -> 1 of cost 5 and 6: a stretch of 1.2 admits the cost-6 link at equality, 11/10 does not.
g = graph(2, [(0, 1, 1, 5), (0, 1, 1, 6)])
for factor in (1.2, Fraction(11, 10)):
    paths = ngc.ksp(g, 0, 1, ngc.KspOptions(2, cost_ceiling=ngc.CostStretch(factor)))
    print(factor, [cost for cost, _ in paths])
```

```text
6/5 True CostStretch(Fraction(4, 1))
1.2 [5, 6]
11/10 [5]
```

## 5. Route allocations and admission

**In plain words.** A route allocation decides how one demand spreads over its routes. `OPTIMIZED`
behaves like a controller that pushes as much as the routes can carry, cheapest walks first.
`WEIGHTED` is a router splitting traffic at every node by fixed weights (WCMP). `EQUAL` is a router
splitting equally at every node (ECMP). Admission is a separate choice: `LOSSLESS` places only
what fits end to end; `LOSSY` forwards hop by hop, as routers do, and drops on each edge what does
not fit.

**Exact rules.**

- **Split graph.** A demand uses the part of its route (a DAG or a path) that reaches its
  destination set through the edges its allocation may use (its members). `OPTIMIZED` uses the
  edges with capacity in its view (the free capacity plus the policy's own allocation); `EQUAL`
  uses every edge of the route, an edge of capacity 0 included: a router splits over every next
  hop it has, so that edge receives its round-robin share and carries none of it (Lossless can then
  admit nothing on that branch; Lossy drops the share on that edge); `WEIGHTED` uses the edges of
  positive weight (below) and drops a branch left without such an edge, inside the existing route
  only (it never takes a longer bypass). Lossy keeps its allocation's membership.
- **`OPTIMIZED`** solves the minimum-cost flow of the demand's cumulative amount on its route,
  over the free capacity plus the policy's own current allocation. It may move the policy's own
  earlier allocation (that is what lets it reach the maximum); it never moves another policy's.
  It fills the cheapest walks of a mixed-cost route first. Per-edge changes between calls can be
  negative. With `ShortestPaths` it solves the whole demand at once, over the union of its
  sources' routes ([section 7](#routing-models)). `OPTIMIZED` is Lossless only: a controller never
  drops. A solve beyond the flow core's size limit raises `ValueError` ([section 2](#2-the-graph),
  Size).
- **`WEIGHTED`** splits at every node by highest averages over integer weights frozen when the
  route is prepared. The **weight basis** (`WeightBasis`) says where they come from: `CAPACITY`
  (the default) the member's static capacity (of the outgoing link, not of the branch behind it);
  `EXPLICIT` one integer per edge in EdgeId order given as `edge_weights`, each in `[0, 2^62)`,
  validated in full (a positive weight on an edge of capacity 0 becomes 0); `RESIDUAL_SNAPSHOT` the
  member's free capacity, read once when the program is prepared (a policy: its view at that
  call's entry, the ledger's residual plus, for `evaluate` and `replace`, the policy's own released
  traffic; `from_routes` and `from_dag`: the `residual` they are given). The tie rule never
  chooses weights. Weights are never renormalised; a member whose residual later shrinks limits
  the demand. Equal positive weights give exactly round robin: `H(x, [w] * k) = RR(x, k)`.
- **`EQUAL`** splits at every node by round robin over its members, which are frozen when the
  program is prepared.
- **Lossless filling** (`WEIGHTED`, `EQUAL`). On a fixed program, a demand's delivered amount `T`
  only grows: a call raises it to the smallest of `T` plus the call's offer and `T_max`, the
  largest amount whose split fits every member edge (found by exact binary search). At `T_max` at
  least one member edge is exactly full; there is no safety margin. Because the split rules are
  monotone, raising `T` never takes units away from an edge.
- **Lossy forwarding** (`ShortestPaths` with `WEIGHTED` or `EQUAL`, T1, T2 and `Pinned`). Each
  call's offer is split at every node by the allocation's rule; each edge carries the smaller of
  what reaches it and its residual, and the rest is dropped on that edge (where several routes of
  one program, or of several demands forwarded together, meet on an edge, [the loss
  rule](#losses-on-shared-edges) divides what it carries). The delivery is what reaches the destination set. A share that has no edge to go on (a
  `WEIGHTED` node whose members all have weight 0) is lost where it is: at a source it is a source
  drop (`source_drops`). Records keep `offered`, `delivered` and per-edge `edge_drops`. Lossy calls
  are independent: nothing is offered again and two calls are not equivalent to one. An infinite
  Lossy offer is refused (`LOSSY forwarding requires a finite offer`) except for `ShortestPaths`
  `place(float('inf'))`, the **first-hop experiment**: every source injects on each of its first
  edges exactly that edge's free capacity, at once, and forwarding runs as usual.

On a diamond with paths of capacity 10 and 2, equal-share forwarding stops at 5 units: round
robin gives the odd unit to the lower EdgeId (3 and 2), and a sixth unit would put 3 units on the
2-unit path. The weighted and optimized allocations carry 12:

```python
g = graph(4, [(0, 1, 10, 1), (1, 3, 10, 1), (0, 2, 2, 1), (2, 3, 2, 1)])
allocations = ngc.RouteAllocation.EQUAL, ngc.RouteAllocation.WEIGHTED, ngc.RouteAllocation.OPTIMIZED
for allocation in allocations:
    ledger = ngc.FlowLedger(g)
    r = ngc.FlowPolicy(ledger, 0, 3, config(ngc.ShortestPaths(allocation=allocation))).place(100)
    print(allocation.name, r.call_delivered, r.call_undelivered, r.termination.name, ledger.carried().tolist())

ledger = ngc.FlowLedger(g)
pol = ngc.FlowPolicy(ledger, 0, 3, config(ngc.ShortestPaths(), admission=ngc.Admission.LOSSY))
r = pol.place(10)
print("LOSSY", r.call_delivered, r.call_undelivered, r.termination.name, "drops", r.call_edge_drops)
rec = pol.records()[0]
print("record", rec.offered, rec.delivered, rec.edge_drops)
```

```text
EQUAL 5 95 EXHAUSTED [3, 2, 3, 2]
WEIGHTED 12 88 EXHAUSTED [10, 2, 10, 2]
OPTIMIZED 12 88 EXHAUSTED [10, 2, 10, 2]
LOSSY 7 3 COMPLETE drops [(1, 3)]
record 10 7 [(1, 3)]
```

The Lossy call is `COMPLETE` although 3 units were dropped: the whole offer of 10 was forwarded,
and what did not fit is accounted as drops, not as a shortfall ([section 7](#termination)).

An edge of capacity 0 on the only shortest route is still that route. Equal splitting sends it
the whole offer: Lossless places nothing, Lossy drops all of it there. Weighted splitting gives
it weight 0 and has no other edge in its route, so nothing can be sent (`NO_ACTIVE_WEIGHT`) and
a Lossy offer is lost at the source; neither takes the costlier detour. A controller
(`CentralizedTE`) uses the detour:

```python
# Edge 0: 0 -> 1 of capacity 0 (cost 1); the detour 0 -> 2 -> 1 (5 units, cost 2).
g = graph(3, [(0, 1, 0, 1), (0, 2, 5, 1), (2, 1, 5, 1)])
r = ngc.spf(g, 0, ngc.SpfOptions(F, targets=1))
print("spf", r.target_cost, r.dag.via_edges().tolist())
for allocation in (ngc.RouteAllocation.EQUAL, ngc.RouteAllocation.WEIGHTED):
    for admission in (ngc.Admission.LOSSLESS, ngc.Admission.LOSSY):
        cfg = config(ngc.ShortestPaths(allocation=allocation), admission=admission)
        r = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 1, cfg).place(5)
        print(allocation.name, admission.name, r.call_delivered, r.call_edge_drops, r.call_source_drops,
              r.termination.name)
cfg = ngc.FlowPolicyConfig(routing=ngc.CentralizedTE(source_allocation=ngc.SourceAllocation.OPTIMIZED))
print("CentralizedTE", ngc.FlowPolicy(ngc.FlowLedger(g), 0, 1, cfg).place(5).call_delivered)
```

```text
spf 1 [0]
EQUAL LOSSLESS 0 [] 0 EXHAUSTED
EQUAL LOSSY 0 [(0, 5)] 0 COMPLETE
WEIGHTED LOSSLESS 0 [] 0 NO_ACTIVE_WEIGHT
WEIGHTED LOSSY 0 [] 5 NO_ACTIVE_WEIGHT
CentralizedTE 5
```

Weights: two parallel links of 10 and an offer of 12. Capacity weights split it 6 and 6, explicit
weights 1 and 2 give 4 and 8. With only 2 units free on the first link (a service baseline of 2
and 10), equal splitting admits 4, while weights snapshot from the free capacity (2 and 10) admit
all 12:

```python
g = graph(2, [(0, 1, 10, 1), (0, 1, 10, 1)])
W = ngc.RouteAllocation.WEIGHTED
for model in (ngc.ShortestPaths(allocation=W),
              ngc.ShortestPaths(allocation=W, weight_basis=ngc.WeightBasis.EXPLICIT, edge_weights=[1, 2])):
    ledger = ngc.FlowLedger(g)
    pol = ngc.FlowPolicy(ledger, 0, 1, config(model))
    print(pol.config.routing.weight_basis.name, pol.place(12).call_delivered, ledger.carried().tolist())
snapshot = ngc.ShortestPaths(allocation=W, weight_basis=ngc.WeightBasis.RESIDUAL_SNAPSHOT)
for model in (ngc.ShortestPaths(), snapshot):
    print(model.allocation.name, ngc.FlowPolicy(ngc.FlowLedger(g, [2, 10]), 0, 1, config(model)).place(12).call_delivered)
```

```text
CAPACITY 12 [6, 6]
EXPLICIT 12 [4, 8]
EQUAL 4
WEIGHTED 12
```

## 6. The ledger

**In plain words.** A `FlowLedger` holds what is placed on one graph, in two projections that are
never added together. The **carried traffic** is a service baseline (the starting residual) and
one record per flow; its residual is always the baseline minus everything the records carry. The
**strict holds** are a reservation baseline and the bandwidth that up T1 tunnels reserve on their
paths; its unreserved capacity is that baseline minus every hold. A hold never polices traffic:
traffic may use bandwidth that is held, and a tunnel may carry more than it holds. Policies
place, hold and remove; the ledger has no public placement, reservation or removal of its own.
Its own `reset()`, rejected while a policy is bound, replaces both baselines and clears the
records.

**Exact rules.**

- **`FlowLedger(graph, service_baseline=None, *, reservation_baseline=None)`.** The service
  baseline B and the reservation baseline R0 are each an integer array in EdgeId order with every
  entry in `[0, edge_capacity_units[e]]` (a physical baseline goes through `to_units` first;
  anything above the capacity is rejected); `None` means the capacities. `service_baseline` is
  checked first, then `reservation_baseline`, with the same messages naming each
  (`service_baseline has length <n>, expected <E> (one entry per edge)`, then
  `service_baseline[<e>] = <v> is outside [0, edge_capacity_units[<e>] = <c>]`). The baselines
  never change except through `reset`. A ledger that should model load from outside the library
  starts from reduced baselines.
- **Snapshots.** `carried()` L is the exact sum of every record's `edge_flows` and `residual()` is
  `B - L`; `reserved()` S is the exact sum of every record's `edge_reservations` and
  `unreserved()` is `R0 - S`; all four are `int64[E]` copies. Each projection keeps its own bound
  on every edge (`0 <= L <= B`, `0 <= S <= R0`); L and S are never added, and `max(L, S)` is no
  occupancy. T1's setup reads `unreserved` (other owners' holds included); every forwarding, every
  `CentralizedTE` solve and T2's layout read `residual` (a T2 plan books no hold). Traffic without
  a hold can use bandwidth that T1 holds: there is no priority queueing, preemption or policing.
  `edge_capacity_units()` is the graph's view and `graph()` returns the graph.
- **Records.** `records()` lists the live records (of bound policies and orphans) as snapshots
  in FlowId order; `record(flow_id)` reads one. A FlowId is the pair `(ledger uid, id)`; ids are
  never reused. A FlowId of another ledger, or a retired or unknown id, raises `ValueError`.
- **`FlowRecord`** fields: `flow_id`; `policy_id` (the ledger-local creation-order identifier of
  the policy, kept after it becomes an orphan); `sources` and `destinations` (node sets);
  `admission`; `entry_id` (a T1, T2 or `Pinned` record: the identifier of the program entry it
  carries; `None` for the records kept per roster source); `route` (the record's route as a
  `PredecessorDAG`); `min_route_cost` (the cheapest walk of the route, cached; `None` for an empty
  route); `offered` (Lossy records only, else `None`); `delivered`; `edge_flows`
  (`[(EdgeId, units)]`, ascending, the carried traffic), `edge_drops` (`[(EdgeId, units)]`,
  cumulative, ascending) and `edge_reservations` (`[(EdgeId, units)]`, ascending: an up T1
  tunnel's request on every edge of its path; `[]` for every other record). A record has no
  allocation field: the policy's configuration says how its traffic was computed and the program
  how it is forwarded.
- **`saturated_boundary_edges(sources, eligibility=None)`** is a saturation report of the carried
  traffic (holds play no part): the edges (eligible, at least 1 unit of capacity) that leave the
  set of nodes reachable from the allowed sources over edges with at least 1 unit of residual, in
  ascending EdgeId. It follows edges forward only. It ignores the destinations, so it is not a
  source-destination cut, and it can be empty while a route exists.
- **Lifetimes.** The graph outlives the ledger and the ledger outlives its policies (Python keeps
  them alive for you).

| Ledger state | Reached by | `reset(service_baseline=None, *, reservation_baseline=None)` | `copy()` |
| --- | --- | --- | --- |
| a policy is bound | building a policy on it | `ValueError` | `ValueError` |
| no policy bound | no policy built yet, or every policy destroyed | allowed: orphaned records with their allocations and holds are cleared, their ids retired, the new baselines set (`None`: the capacities) | a new ledger with a new uid and the same baselines, records, holds and ids |

A policy destroyed without `remove()` leaves its records, allocations and holds in the ledger as
**orphans**: readable, never changed again, cleared only by `reset`.

The two projections of one link of 10, where a T1 tunnel holds 6 and
carries 4; a second demand on the same ledger places the residual 6, whatever is held:

```python
g = graph(2, [(0, 1, 10, 1)])
ledger = ngc.FlowLedger(g)
t1 = ngc.FlowPolicy(ledger, 0, 1, ngc.FlowPolicyConfig(routing=ngc.DistributedTE(ngc.TunnelCount(1))))
t1.replace(4, design_demand=6)
print(ledger.carried().tolist(), ledger.residual().tolist(),
      ledger.reserved().tolist(), ledger.unreserved().tolist())
print(ngc.FlowPolicy(ledger, 0, 1, config(ngc.ShortestPaths())).place(10).call_delivered)
```

```text
[4] [6] [6] [4]
6
```

## 7. Policies

**In plain words.** A `FlowPolicy` is the routing model of one demand between a set of sources
and a set of destinations. `place(offer)` asks it to carry more; it routes and fills its flows,
commits everything at once, and reports how much was delivered, what stopped it, where traffic
was lost and how much each source carried (T1 and T2, whose tunnels are sized for a design
volume, take whole offers only). `evaluate`, `replace` and `forward` ask for a whole new offer
instead (only computed, published, or forwarded over the routes already installed). `remove()`
takes all of its traffic back out. A demand with several sources is one demand, not several: how
its traffic enters at the sources is its **source allocation**, and where the routes of its
sources meet, their traffic is added and split once. Several policies of one ledger can act in one
operation, `apply_actions`: one after another, or all at the same time with their losses shared
([several demands at once](#several-demands-at-once)).

### Configuration

    FlowPolicyConfig(*, routing, admission=Admission.LOSSLESS)
    ShortestPaths(*, allocation=RouteAllocation.EQUAL, source_allocation=None, weight_basis=None,
                  edge_weights=None, tie_break=RouteTieBreak.FEWER_HOPS, multipath=True,
                  parallel_edges=True)
    DistributedTE(tunnel_count, *, cost_ceiling=None, min_unreserved_bandwidth=0,
                  tie_break=RouteTieBreak.FEWER_HOPS, parallel_edges=True)          # T1
    CentralizedTE(*, source_allocation, source_weights=None, support=None, tunnel_count=None,
                  prefer_headroom=False)                                          # T2 with tunnel_count
    Pinned(program, *, black_holes=None)
    TunnelCount(per_head, *, up_to=False)

`FlowPolicyConfig` holds only the routing model and the admission; every other option belongs to
the model that reads it. Every value takes keywords only after its required arguments, and every
field is read-only (a model is an immutable value; its arrays are shared by reference, so policies
built from one model read one `edge_weights` array).

| Routing | Models | Fields |
| --- | --- | --- |
| `ShortestPaths()` | IGP forwarding (as OSPF or IS-IS routers do): static shortest paths, never refreshed. `ShortestPaths()` is ECMP; `allocation=WEIGHTED` is WCMP; `allocation=OPTIMIZED` is a controller restricted to those paths | `allocation`, `source_allocation`, `weight_basis`, `edge_weights` and the route choice (`tie_break`, `multipath`, `parallel_edges`) |
| `DistributedTE(TunnelCount(N))` (T1) | distributed strict-bandwidth tunnels, as RSVP-TE head ends set them up: every head runs CSPF for its own tunnels, no controller ([T1](#t1-tunnels-set-up-by-their-heads)) | `tunnel_count`; `cost_ceiling`: `None`, `AbsoluteCost(c)` or `CostStretch(f)`; `min_unreserved_bandwidth`: a search minimum in units; `tie_break`, `parallel_edges` |
| `CentralizedTE(source_allocation=...)` | an ideal controller that sees the whole network, or only the edges of `support`; with `tunnel_count`, T2: the controller lays out equal tunnels for every head ([T2](#t2-the-controllers-equal-tunnels)) | `source_allocation` (required); `source_weights`: one integer per graph node, required exactly for `WEIGHTED` without tunnels; `support`: a `PredecessorDAG` of the graph, or `None` for every eligible edge; `tunnel_count`: a `TunnelCount(N)` or `TunnelCount(L, up_to=True)`, or `None`; `prefer_headroom`: T2's headroom objective after route cost |
| `Pinned(program)` | forwarding over a given program, computing nothing | `program`: a `ForwardingProgram`; `black_holes`: an `Eligibility` whose excluded nodes and edges are black holes, or `None` |

**Two stages of validation.** A constructor converts its arguments: an unknown keyword (also every
field that belongs to another model, such as `allocation` on `FlowPolicyConfig` or
`source_allocation` on `DistributedTE`), a missing required argument (`routing`, `tunnel_count`,
`CentralizedTE`'s `source_allocation`, `program`), `None` where none is allowed, and a value of the
wrong type raise `TypeError` before any policy exists. A bool field (and a bool search option)
takes only a `bool` or a `numpy.bool_`, known by its type (`multipath=None`, `multipath=1` and an
object with `__bool__` are `TypeError`). An integer field of a configuration value (`TunnelCount`
`per_head`, `AbsoluteCost` `cost`, `min_unreserved_bandwidth`, a `RouteDefinition`'s `head` and
`weight`) given a value of another type raises `TypeError` (`TunnelCount per_head must be an
integer (got bool)`), and so does an `edge_weights` that is not a 1-D sequence of integers (a
scalar, a string, a 2-D array, a bool or float array, a bool or float entry:
`edge_weights[0] must be an integer (got bool)`); an integer outside 64 bits raises
`ValueError` (`TunnelCount per_head = <v> is out of range [...]`), and an `edge_weights` entry
beyond 64 bits `ValueError('edge_weights[<i>] must be an integer in [0, 2^62)')`. Every other
check of a well-typed configuration raises `ValueError` when a policy or an evaluator uses it, in
the order below, before any work, a zero offer included.

**The model matrix.** The options belong to each row's model; any other well-typed combination
raises `ValueError`:

| Model | route allocation | source allocation: default / allowed | `admission` | route choice | weights |
| --- | --- | --- | --- | --- | --- |
| `ShortestPaths()` | `EQUAL` (ECMP, the default) | `EQUAL` / `EQUAL` | `LOSSLESS`, `LOSSY` | `tie_break`, `multipath`, `parallel_edges` | none |
| `ShortestPaths(allocation=WEIGHTED)` | `WEIGHTED` (WCMP) | `EQUAL` / `EQUAL` | `LOSSLESS`, `LOSSY` | the same | `weight_basis`: `CAPACITY` (default), `EXPLICIT`, `RESIDUAL_SNAPSHOT` |
| `ShortestPaths(allocation=OPTIMIZED)` | `OPTIMIZED` | `OPTIMIZED` / `EQUAL`, `OPTIMIZED` | `LOSSLESS` | the same | none |
| `CentralizedTE(source_allocation=...)`, optionally `support` | optimal | required / EQUAL, OPTIMIZED, WEIGHTED | `LOSSLESS` | none | `source_weights` exactly with WEIGHTED |
| T1: `DistributedTE(TunnelCount(N))` | tunnels | `EQUAL`, fixed | `LOSSLESS`, `LOSSY` | `tie_break`, `parallel_edges` | none |
| T2: `CentralizedTE(source_allocation=..., tunnel_count=...)`, optionally `support` | tunnels | required; `EQUAL`; `OPTIMIZED` only with one source | `LOSSLESS`, `LOSSY` | none (`prefer_headroom`) | none |
| `Pinned(program)` | the program's | the program's | `LOSSLESS`, `LOSSY` | none | none |

`ShortestPaths`' `source_allocation=None` resolves to the allocation's default and its
`weight_basis=None` to `CAPACITY` for `WEIGHTED`; a policy's `config` is the configuration it was
built with, these two fields resolved (the model value itself never changes). The texts of the
matrix, in the order they are checked:

1. the enumeration values (`admission`, then the model's `allocation`, `source_allocation`,
   `weight_basis`, `tie_break`), each `<field>: <v> is not a value of <Enum>` (reachable with a
   value built from an integer);
2. the route allocation and the source allocation: first `source_allocation=WEIGHTED requires
   CentralizedTE without tunnel_count`; then `ShortestPaths with EQUAL or WEIGHTED requires
   source_allocation=EQUAL`; for T2, `CentralizedTE with tunnel_count and several sources requires
   source_allocation=EQUAL` (`OPTIMIZED` with more than one input source, counted after removing
   repeats and before any Eligibility; with one source it is the same calculation as `EQUAL`);
3. the admission: `ShortestPaths OPTIMIZED requires admission=LOSSLESS`; `CentralizedTE without
   tunnel_count requires admission=LOSSLESS`;
4. the route choice: `ShortestPaths tie_break requires multipath=False or parallel_edges=False` (a
   `HIGHER_*` rule with every route and every parallel edge kept, also for `OPTIMIZED`, whose route
   the rule selects); then the **pooled single-path rule**: `ShortestPaths with multiple input
   sources requires tie_break=FEWER_HOPS when multipath=False` (more than one input source, counted
   after removing repeats and before any Eligibility: a bottleneck choice could differ per source
   at a router the sources share);
5. `prefer_headroom requires tunnel_count`;
6. the weights: `weight_basis requires allocation=WEIGHTED`; `weight_basis=EXPLICIT requires
   edge_weights`; `edge_weights require weight_basis=EXPLICIT`; for CentralizedTE,
   `source_allocation=WEIGHTED requires source_weights` and `source_weights require source_allocation=WEIGHTED`;
7. last, a cache given to a model that does not read one: `<DistributedTE or CentralizedTE or
   Pinned> does not use ShortestPathCache` (only `ShortestPaths` reads a `ShortestPathCache`).

Then each model's own fields: the tunnel count in `[1, 10000]` (`tunnel count must be in [1,
10000]`; `TunnelCount(per_head)` itself takes any 64-bit integer); for T1 an "up to" count
(`DistributedTE requires an exact tunnel count`), then the ceiling (`cost_ceiling AbsoluteCost
cost must be >= 0 (got <c>)`, `cost_ceiling CostStretch factor must be >= 1 (got <n>/<d>)`), then
`min_unreserved_bandwidth must be an integer in [0, 2^62)`; `edge_weights`' length and entries
(`edge_weights: expected <E> entries`, `edge_weights[<i>] must be an integer in [0, 2^62)`).

`FlowPolicy(ledger, sources, destinations, config, eligibility=None, cache=None)` binds the policy
to its ledger. `sources` and `destinations` are one node id or a sequence; the two sets must not
overlap. Construction validates everything before the policy binds, and when several defects
meet, the first in this order is reported: the endpoint sets (empty, out of range, the lowest node
in both); the graph of each object the policy is given: the `Eligibility` (`foreign Eligibility:
built for another graph`), the program (`program belongs to another graph`), the support
(`support belongs to another graph`, with or without a tunnel count, an empty support included),
the `ShortestPathCache` (`foreign ShortestPathCache: bound to another graph`), the black holes
(`black_holes belongs to another graph`); the matrix and the models' fields above; the tunnel limit
(T1, T2: the roster's heads times the tunnel count, for "up to" the largest count, is at most
2^30 - 1: `tunnel limit exceeded: <H> heads x <N> tunnels > 1073741823`; no tunnel is made at
construction); `Pinned`'s endpoints (`program endpoints do not match policy endpoints`; any number
of each).

`source_weights` uses the same conversion as `edge_weights`: first the entire shape (rank and
nested sequences), then every integer type, then signed-int64 storage overflow. Thus `[1.5,[1]]`
reports the second entry's shape, and `[2**64,1.5]` reports the second entry's type. Ordinary
sequences such as range and tuples are accepted; unsigned arrays do not bypass overflow checks.
Policy/evaluator validation checks the graph node count (`source_weights: expected <N> entries`),
then every entry in NodeId order (`source_weights[<i>] must be an integer in [0, 2^62)`) after the
model/cache checks, even at zero offer. Each individual weight is bounded; their sum may exceed
64 bits. Cumulative Pinned LOSSY increments use exact wide apportionment through the existing
counter domain below 2^127, preserving declared proportions after many offers.

### The source allocation

The source allocation says how a demand's total enters at its roster sources, separately from how
it is routed:

- **`EQUAL`**: a total `x` enters as the round robin `RR(x, k, i) = x // k + (1 if i < x % k else
  0)` at roster member `i` (`k` members, ascending node id): the first members get one extra
  unit. A source whose share cannot go anywhere keeps it and loses it; the share is never moved
  to another source. Lossless therefore admits only totals that give such a source nothing, and
  Lossy loses its share.
- **`OPTIMIZED`** (controllers only): the controller chooses how much enters at each source (one
  virtual source joined to every roster source), whatever maximises the delivery and then
  minimises its cost.

- **`WEIGHTED`** (`CentralizedTE` without `tunnel_count` only): each roster source injects its
  highest-averages share `H(x, w)` of the total. `source_weights` is required: one integer per graph
  node in `[0, 2^62)`, copied once and exposed as an immutable tuple. Every entry is validated,
  including unused and masked nodes, but only roster weights participate. Zero weights receive
  nothing; all-zero roster weights give no injection, with no equal fallback. Other source
  allocations refuse this vector; ShortestPaths and every tunnel configuration refuse WEIGHTED.
`DistributedTE` always allocates equally and `Pinned` takes its program's injection; neither
accepts the keyword.

### The source roster

The roster is fixed when the policy is built. For `ShortestPaths`, T1, T2 and `CentralizedTE` it
is the allowed sources that have a static walk to an allowed destination over allowed edges,
whatever their capacity, in ascending node id: a source whose only routes cross an edge of
capacity 0 is in the roster, receives its share and cannot place it (a failure is expressed with
an Eligibility, which removes the source); `CentralizedTE` with a support counts only walks inside
the support. For `Pinned` it is the program's heads that still reach a destination under the
policy's Eligibility ([replay](#replay-under-a-failure)). The roster decides which sources receive shares; for the
computing models it also decides which records exist (T1 and T2: its heads' entries), while `Pinned`
keeps one record per configured program entry even when no head survives the replay (such a
record stays at zero). `by_source`
reports every input source, a source outside the roster with a zero row. A positive offer with
an empty roster ends `NO_ROUTE`: a computing model creates nothing, and `Pinned` creates or keeps its
per-entry records at zero.

### Routing models

- **`ShortestPaths` with `EQUAL` or `WEIGHTED`** forwards like IGP routers. At the epoch's
  first positive call it prepares a **program**, the routes and splits used for the whole epoch:
  each roster source's route is its static shortest-path DAG toward the destination set (with
  `multipath=False`, its canonical path), read from the [shortest-path cache](#9-the-shortest-path-cache),
  with the route choice applied on the ledger's residual at that call's entry (the
  `HIGHER_AVAILABLE_CAPACITY` and `HIGHER_HEADROOM` scores read that residual); `WEIGHTED` freezes
  its weights from its weight basis there ([section 5](#5-route-allocations-and-admission)). The
  program is built by the one pooled DAG builder that
  [`ForwardingProgram.from_dag`](#supplied-dags-from_dag) uses, from the union of these selected
  routes: at a router several sources pass, their traffic for the demand is added (with the
  router's own injection when it is a source) and split once over its members, by round robin
  (`EQUAL`) or by highest averages over the weights (`WEIGHTED`). Every source that passes a
  router has the same continuation there, so the union lends no source a link outside its own
  route. A roster of one source with `multipath=True` and a static weight basis (not
  `RESIDUAL_SNAPSHOT`) takes a one-pass preparation instead, which selects the route and builds its
  split in the same pass; its plan is the one the pooled builder makes of that route, so the
  program is the same. The program is never refreshed or recomputed within the epoch. It is the
  policy's `program`: one split context pooled over the roster, one split route per roster source
  (entry identifier = roster position, weight 1), equal sources; it lists no tunnel. With
  `FEWER_HOPS` or `HIGHER_CAPACITY` and `EQUAL`, `CAPACITY` or `EXPLICIT` weights, a preparation
  reads nothing that changes while the policy lives, so the policy keeps the program it first
  published and uses it again at each later epoch's preparation: routes, order, costs and splits
  equal a fresh preparation. Each policy keeps only its own program, until it is destroyed;
  nothing is added to the cache or its statistics. `evaluate` uses a kept program but never keeps
  one. `HIGHER_AVAILABLE_CAPACITY`, `HIGHER_HEADROOM` and `RESIDUAL_SNAPSHOT` prepare afresh at
  every epoch.
- **`ShortestPaths` with `OPTIMIZED`** is a controller restricted to shortest paths: it solves like
  `CentralizedTE` (below), but only over the union of the roster sources' routes prepared as
  above, edges of capacity 0 included (so a zero-capacity shortest route is never bypassed).
- **`DistributedTE` (T1)** sets up its tunnels for a design volume on every `replace` and
  `evaluate` and forwards them as its program; it has no incremental `place`
  ([T1](#t1-tunnels-set-up-by-their-heads)).
- **`CentralizedTE(tunnel_count=...)` (T2)** lays out its tunnels for a design volume on every
  `replace` and `evaluate` and forwards them as its program; it has no incremental `place`
  ([T2](#t2-the-controllers-equal-tunnels)). The controller rules below are `CentralizedTE`'s
  without a tunnel count.
- **`CentralizedTE`** re-solves, on every positive call, the minimum-cost flow of the amount the
  policy already delivers plus this call's offer (`delivered()` plus the offer; what an earlier
  call could not place is not offered again) from its roster to its whole destination set over
  every eligible edge, on the free capacity plus its own allocation (it may move its own traffic,
  never another owner's): the largest delivery up to that amount, and among those the cheapest.
  Hence `place(a)` followed by `place(b)` ends in the same allocation as `place(a + b)` on an
  unchanged network. With `EQUAL` source allocation the virtual source's arc to roster source `i`
  carries exactly `RR(T, k, i)`: the delivered total `T` is the largest whose equal shares all fit
  (found by exact binary search over the totals; a total fits only when the solve reaches it), and
  that total is solved from zero, so the result does not depend on the search. WEIGHTED uses
  `H(T, w)` instead of round robin in the same exact search: every root arc must be saturated,
  then minimum cost is chosen. Each share is nondecreasing in T, so feasible totals form a prefix.
  Its maximum is not an unrestricted minimum-cut certificate: `min_cut` and
  `marginal_cost_curve` are None. A
  positive call on a network beyond the flow core's size limit raises `ValueError`
  ([section 2](#2-the-graph), Size).
- **`CentralizedTE(support=dag)`** is the same controller restricted to the edges of `dag`,
  trimmed to the walks from the sources to the allowed destinations inside it (an edge that leads
  nowhere inside the support is never used; a mask removes its edges as usual). It reports
  `total_cost` but neither `marginal_cost_curve` nor `min_cut`, and `evaluate_edge_failures`
  solves each deletion afresh.
- **The controllers' programs.** After every positive call a controller decomposes its solved
  flow into paths, each roster source's share in a canonical order (for each source in node id
  order, walk the smallest positive EdgeId to the first destination, take the bottleneck or the
  remaining supply, subtract, repeat), one tunnel per path with its solved amount as weight
  (identifiers in decomposition order), also with a single source. With the `OPTIMIZED` source
  allocation the program keeps the solved supplies as **captured source weights**
  (`source_weights`: every roster source, a zero supply included). With `WEIGHTED` it instead
  keeps each roster source's declared weight, zeros included; with `EQUAL` it injects by
  round robin. Forwarded on the same view, the program reproduces the solve exactly.
- **Records of the controllers** (`CentralizedTE`, `ShortestPaths` with `OPTIMIZED`): one record
  per roster source, holding that source's share of the solved flow (its edges and its delivery),
  extracted in a canonical order. A `CentralizedTE` record's route is the set of edges that carry
  the share, rebuilt on every positive call under the same FlowId; a source with no share has an
  empty route and `min_route_cost` `None`. A `ShortestPaths` `OPTIMIZED` record's route is its
  source's prepared route, and the share stays inside it.
- **`Pinned`** forwards its program as given, computing nothing: injection, entries, weights and
  split tables are the program's, replayed under the policy's Eligibility and with its black
  holes ([forwarding programs](#forwarding-programs)). It keeps one record per program entry
  (`entry_id` is the entry's identifier), created at the epoch's first positive call, a record of
  weight 0 or of a `FAILED` tunnel included (that one with an empty route and `min_route_cost`
  `None`). A positive offer with an empty replay roster ends `NO_ROUTE`.

### Several sources and destinations

- **Destinations** (`ShortestPaths`): a route goes to the nearest destinations, and one flow may
  deliver to several of them at that cost (a T1 tunnel is one path, to the allowed destination its
  CSPF reaches first). A record's `delivered` is its total; the
  delivery to one destination is the sum of the record's `edge_flows` entering it.
- **`ShortestPaths` split routing pools its sources.** Each record (one per roster source) holds the
  part of the program's traffic that is that source's. At every router the members' amounts are
  given back to the sources member by member, in member order, by highest averages over what each
  source still holds there; on an edge that drops traffic (Lossy), what the edge carries is
  divided by highest averages over what each source brought to it. Highest averages `H(x, w)`
  gives `x` units one at a time to the member with the largest `w[i] / (assigned[i] + 1)`, ties to
  the earlier member; a zero weight gets nothing. So the records sum exactly to the totals, and a
  source's record uses only edges of its own route.
- **The controllers** (`CentralizedTE`, `ShortestPaths` with `OPTIMIZED`) optimise the whole
  roster and destination set together; `CentralizedTE` may use a costlier destination to carry
  more.
- **T1 and T2 heads** share the design and the offer by round robin over the roster in ascending
  node id; a head's share goes only to its own tunnels, and what one head cannot carry is lost at
  that head, never moved to another.
- **`Pinned`** follows its program: its injection over the replay roster, then each head's entries.

Two sources with one unit each meet at node 2, whose two members lead to the destination. Pooled,
node 2 splits the two units once, one per member, and gives one member to each source; split
separately, each source would send its unit to the first member:

```python
# Edges 0: 0 -> 2, 1: 1 -> 2, then the two members 2: 2 -> 3 and 3: 2 -> 3.
g = graph(4, [(0, 2, 5, 1), (1, 2, 5, 1), (2, 3, 1, 1), (2, 3, 1, 1)])
pol = ngc.FlowPolicy(ngc.FlowLedger(g), [0, 1], 3, config(ngc.ShortestPaths()))
r = pol.place(2)
print(r.call_delivered, r.termination.name, [(row.source, row.delivered) for row in r.by_source])
print([(rec.sources, rec.edge_flows) for rec in pol.records()])
```

```text
2 COMPLETE [(0, 1), (1, 1)]
[([0], [(0, 1), (2, 1)]), ([1], [(1, 1), (3, 1)])]
```

### Equal source allocation and a source that cannot pass anything

A roster source whose every route crosses an edge with no free capacity stays in the roster: it
receives its share and loses it. Lossless admits only totals that give it nothing, which is the
number of roster sources before it, because the first remainder units go to the lowest node ids:
0 when it is the first source. This holds in every model. Lossy loses its share where its program
puts it: on the edge that limits it (`edge_drops`), or at the source when it has no member to send
it on (`source_drops`). An Eligibility that excludes the source or the edge removes it from the
roster instead.

```python
# Source 0 reaches node 2 only over the 0-unit edge 0; source 1 over the 10-unit edge 1.
g = graph(3, [(0, 2, 0, 1), (1, 2, 10, 1)])
equal = ngc.FlowPolicyConfig(routing=ngc.CentralizedTE(source_allocation=ngc.SourceAllocation.EQUAL))
optimized = ngc.FlowPolicyConfig(routing=ngc.CentralizedTE(source_allocation=ngc.SourceAllocation.OPTIMIZED))
for name, cfg in (("ShortestPaths", config(ngc.ShortestPaths())),
                  ("CentralizedTE EQUAL", equal),
                  ("CentralizedTE OPTIMIZED", optimized)):
    r = ngc.FlowPolicy(ngc.FlowLedger(g), [0, 1], 2, cfg).place(10)
    print(name, r.call_delivered, r.termination.name)
lossy = config(ngc.ShortestPaths(), admission=ngc.Admission.LOSSY)
r = ngc.FlowPolicy(ngc.FlowLedger(g), [0, 1], 2, lossy).place(10)
print("LOSSY", r.call_delivered, r.call_edge_drops, [(row.source, row.total_edge_drops) for row in r.by_source])
r = ngc.FlowPolicy(ngc.FlowLedger(g), [0, 1], 2, config(ngc.ShortestPaths()),
                   g.eligibility(excluded_nodes=[0])).place(10)
print("source 0 excluded", r.call_delivered)
```

```text
ShortestPaths 0 EXHAUSTED
CentralizedTE EQUAL 0 EXHAUSTED
CentralizedTE OPTIMIZED 10 COMPLETE
LOSSY 5 [(0, 5)] [(0, 5), (1, 0)]
source 0 excluded 10
```

### Forwarding programs

**In plain words.** A program is what a routing decided, written down: which routes each source
uses and in which integer proportions. Every computing model produces one (`policy.program`,
`result.program`). `Pinned(program)` forwards a program without computing anything: to replay a
solution on another ledger, to study a failure (the routers drop the failed routes and split over
the survivors) or to study the moment before anyone reacts (the failed elements swallow what
reaches them). `ForwardingProgram.from_routes` builds a program from routes you give, each its own
entry; `ForwardingProgram.from_dag` builds one shared hop-by-hop program on a DAG you give, its
sources pooled exactly as `ShortestPaths` pools them.

**Exact rules.**

- **The value.** A `ForwardingProgram` is immutable and shared by reference: `sources` and
  `destinations` (the normalised input sets), `tunnels` (the tunnel entries in ascending
  `entry_id`; split routes are not listed) and `source_weights` (`None` for equal sources). A
  `Tunnel` has `entry_id`, `head`, `local_index` (its position among its head's entries of both
  kinds), `destination`, `edges` (in path order), `path_cost` (these three `None` for a tunnel
  T1's setup left down), `status` (`UP`, or a down tunnel's `INSUFFICIENT_BANDWIDTH` or
  `COST_CEILING_EXCEEDED`), `requested_bandwidth` and `reserved_bandwidth` (T1's request and what
  it holds on every edge of its path: the request when up, else 0; 0 for the other producers),
  `planned_bandwidth` (T2's planned bandwidth u, never held; 0 for the other producers) and
  `weight`. It has no public constructor.
- **What a program shows.** A program can be replayed without exporting its tables. The public
  value does not expose per-router effective weights, the split tables or a reachable-source list
  (its `sources` are the normalised input sources); a split program may have no tunnels while it
  carries traffic. Records and `OutputOptions` give the carried and dropped amounts per edge.
- **Injection.** A total `x` enters at the roster (the heads with an entry; for T1 and T2 the
  policy's roster, a head without an up tunnel or without any tunnel included) either by round
  robin (**equal sources**) or by highest averages over the captured source weights; each head
  then splits its share by highest averages over its entries' weights (equal weights give round
  robin). A head whose entries all have weight 0 loses its share at the source. Identical routes
  are never merged: each entry takes its own share. `ShortestPaths` split routing, `from_dag`, T1,
  T2, equal controllers and supplied routes use equal sources; optimized controllers use their
  solved supplies.
- **Entry identifiers.** Every entry has an integer `entry_id` that identifies it within one
  program (not across recomputations): T1 and T2 tunnels are numbered `j * H + i` (round `j`,
  roster position `i`, `H` heads; T2's dead heads leave holes); controller paths are numbered in
  decomposition order; a `ShortestPaths` split program and a `from_dag` program number their routes
  by roster position; supplied routes keep their position in the `routes` list (a split route uses
  its position too, so tunnel identifiers can skip it). Records (`FlowRecord.entry_id`) and
  `by_tunnel` rows (`TunnelFlow.entry_id`) use the same identifier.
- **Split contexts.** A split route given to `from_routes` is its own context: its routers split
  only that route's share, even where another split route reaches the same router with the same
  members. A `ShortestPaths` program and a `from_dag` program pool all their sources in one
  context.
- **`ForwardingProgram.from_routes(graph, sources, destinations, routes, *, allocation=None,
  weight_basis=None, edge_weights=None, residual=None)`** with `RouteDefinition(head, route, *,
  weight=1)` (`route` a `Path` with one edge per hop, a tunnel, or a `PredecessorDAG`, a split
  route). Order and duplicates are kept; the roster is the sources with a route. `allocation` is
  the interior split of the `PredecessorDAG` routes, `EQUAL` (`None`) or `WEIGHTED`, whose weights
  `weight_basis`, `edge_weights` and `residual` choose as in
  [section 5](#5-route-allocations-and-admission) (`residual` is the free capacity
  `RESIDUAL_SNAPSHOT` reads; it has no default). A split route keeps the edges of its DAG on the
  walks from its head to a destination. The first failure, in this order, raises `ValueError`:
  the endpoint sets; `routes must not be empty`; when no route is a `PredecessorDAG`, any option
  given (`<allocation, weight_basis, edge_weights or residual> applies only to PredecessorDAG
  routes`); `OPTIMIZED` (`from_routes allocation must be EQUAL or WEIGHTED`); the weights'
  combination (`weight_basis requires allocation=WEIGHTED`, `weight_basis=EXPLICIT requires
  edge_weights`, `edge_weights require weight_basis=EXPLICIT`, `weight_basis=RESIDUAL_SNAPSHOT
  requires residual`, `residual requires weight_basis=RESIDUAL_SNAPSHOT`); `edge_weights`' length
  and entries, `residual`'s (the ledger's baseline messages); then each route in list order:
  `routes[<i>]: route belongs to another graph`, `routes[<i>]: head <h> is not a source`, for a
  `Path` `routes[<i>]: a tunnel has one edge per hop`, `routes[<i>]: the path starts at <v>, not at
  head <h>`, `routes[<i>]: the path ends at <v>, not at a destination`, `routes[<i>]: the path
  meets destination <v> before its end`, for a `PredecessorDAG` `routes[<i>]: no walk from head
  <h> to a destination`, then `routes[<i>]: weight must be an integer in [0, 2^62)`. A route of
  weight 0 is validated in full and carries nothing.
- **The lossless bound.** A Lossless call forwards each entry's share without drops, adds the
  loads per edge and admits the largest total that leaves no share without a route and puts no
  edge above its free capacity (exact bisection: the totals that fit form a prefix).

#### Supplied DAGs: from_dag

`ForwardingProgram.from_dag(graph, sources, destinations, dag, *, allocation=None,
weight_basis=None, edge_weights=None, residual=None)` builds one shared hop-by-hop program on a
`PredecessorDAG` of the graph, trimmed to the walks from each source to its first destinations.
It is the mechanism `ShortestPaths` `EQUAL`/`WEIGHTED` uses for its own routes (one builder):
sources are allocated equally and pooled at every shared router with exactly the same splitting
and attribution, in one split context. `allocation` is `EQUAL` (`None`) or `WEIGHTED`; the weights
and the snapshot follow [section 5](#5-route-allocations-and-admission), and `RESIDUAL_SNAPSHOT`
needs an explicit `residual`. The roster is the sources that reach a destination in the DAG,
ignoring capacity and weights; an empty roster is allowed (an empty program). There is one split
route per roster source, identified by roster position, and no tunnel; a `Pinned` policy keeps one
record per route, attributed to its source. A replay rebuilds the roster on the DAG and keeps the
pooled context; stale forwarding keeps the tables and the roster. The first failure, in this
order, raises `ValueError`: the endpoint sets; `dag belongs to another graph`; the allocation
(`from_dag allocation must be EQUAL or WEIGHTED`) and the weights' combination as for
`from_routes`; `edge_weights`' and `residual`'s lengths and entries. A `dag` of another type is a
`TypeError` (`dag must be a PredecessorDAG (got <type>)`).

Two sources whose routes meet at node 2, which has two links of 1 to the destination. `from_dag`
pools the two arrivals and splits them once, one per link; two separate split routes of
`from_routes` each send their unit to the first link:

```python
g = graph(4, [(0, 2, 10, 1), (1, 2, 10, 1), (2, 3, 1, 1), (2, 3, 1, 1)])
dag = ngc.PredecessorDAG.from_csr(g, [0, 0, 0, 2, 4], [0, 1, 2, 2], [0, 1, 2, 3])
pooled = ngc.ForwardingProgram.from_dag(g, [0, 1], 3, dag)
separate = ngc.ForwardingProgram.from_routes(g, [0, 1], 3, [ngc.RouteDefinition(0, dag), ngc.RouteDefinition(1, dag)])
for name, prog in (("from_dag", pooled), ("from_routes", separate)):
    cfg = ngc.FlowPolicyConfig(routing=ngc.Pinned(prog))
    print(name, prog.tunnels, ngc.FlowPolicy(ngc.FlowLedger(g), [0, 1], 3, cfg).place(float("inf")).call_delivered)
sp = ngc.FlowPolicy(ngc.FlowLedger(g), [0, 1], 3, config(ngc.ShortestPaths())).place(float("inf"))
print("ShortestPaths", sp.call_delivered)
```

```text
from_dag () 2
from_routes () 1
ShortestPaths 2
```

#### Replay under a failure

A `Pinned` policy's Eligibility is the failure being studied (local repair):

- A tunnel that touches an excluded node or edge is `FAILED` and gets nothing; replay never
  computes a replacement path. A tunnel T1's setup left down stays down with its status: replay
  never brings one up.
- A split route loses its excluded members and the branches left without a destination; members
  of capacity 0 stay (`EQUAL` still sends them their share).
- The roster is rebuilt from the program's **basis**, ignoring capacities and weights: for an
  unrestricted `CentralizedTE`, T1 and T2 without a support the topology under the producing
  policy's exclusions; for `ShortestPaths`, `from_dag`, a supported `CentralizedTE` (with or
  without a tunnel count) and supplied routes the original routes or support; each together with
  the policy's exclusions. A head that no longer reaches a destination leaves the pool (equal sources: its
  share moves to the others); a head that still reaches one but has no surviving entry keeps its
  share and loses it. Captured source weights keep the remaining heads' saved weights; when all of
  them are 0, a positive offer is lost at the sources. The call ends NO_ACTIVE_WEIGHT when the
  program has an entry with a path. A program without one (an all-zero WEIGHTED program) ends
  COMPLETE under LOSSY, with the offer in demand-level `source_drops` and every per-source row
  zero; under LOSSLESS it admits nothing and ends EXHAUSTED.
- An edge whose capacity is now 0 is not a failure: its traffic is lost on that edge, the tunnel
  stays `UP`. Weights stay as the program froze them.

#### Stale forwarding: black holes

`Pinned(program, black_holes=e)` and `forward(offer, black_holes=e)` forward the program unchanged
while the nodes and edges `e` excludes swallow what reaches them (`None`, or an all-allowed value,
means none). Nothing is pruned, re-split or reweighted: every router keeps its members and
weights and applies the program's whole-unit rule to the traffic that actually arrives. A
black-holed edge carries nothing (all that arrives is dropped there); a black-holed node acts the
same on every edge entering it; a black-holed source loses its own share as a source drop. Lossy
drops what reaches a black hole where it arrives: it takes no part in sharing an edge (the loss
rule below) and makes no edge wait, and the traffic behind it can carry nothing. Lossless admits
only totals that put nothing on a black hole (a black-holed edge serves 0). A tunnel that is up but whose path reaches a black hole is `BLACKHOLED` in
`by_tunnel`, its traffic dropped at the first black hole on its path. The policy's Eligibility
applies first: an element both excluded and black-holed is excluded. A down tunnel keeps its
setup status. A program forwarded by `Pinned` books no hold, whatever its producer.

```python
# Two parallel edges of 4 and 8: the duplicate tunnel takes its own share.
g = graph(2, [(0, 1, 4, 1), (0, 1, 8, 1)])
def tunnel(g, head, edges, weight=1):
    return ngc.RouteDefinition(head, ngc.Path(g, [[e] for e in edges]), weight=weight)
def pinned(program, admission=ngc.Admission.LOSSLESS, black_holes=None):
    return ngc.FlowPolicyConfig(routing=ngc.Pinned(program, black_holes=black_holes),
                                admission=admission)
for routes in ([tunnel(g, 0, [0]), tunnel(g, 0, [1]), tunnel(g, 0, [1])],
               [tunnel(g, 0, [0]), tunnel(g, 0, [1])]):
    prog = ngc.ForwardingProgram.from_routes(g, 0, 1, routes)
    r = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 1, pinned(prog)).place(float("inf"))
    print(r.call_delivered, [(t.entry_id, t.assigned) for t in r.by_tunnel])

# Two tunnels 0 -> 1 -> 3 and 0 -> 2 -> 3, every edge 10: node 2 fails.
g = graph(4, [(0, 1, 10, 1), (0, 2, 10, 1), (1, 3, 10, 1), (2, 3, 10, 1)])
prog = ngc.ForwardingProgram.from_routes(g, 0, 3, [tunnel(g, 0, [0, 2]), tunnel(g, 0, [1, 3])])
down = g.eligibility(excluded_nodes=[2])
lossy = pinned(prog, ngc.Admission.LOSSY)
r = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 3, lossy, down).replace(10)
print("repair", r.call_delivered, [(t.status.name, t.assigned) for t in r.by_tunnel])
r = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 3, pinned(prog, ngc.Admission.LOSSY, down)).replace(10)
print("stale", r.call_delivered, r.call_edge_drops, [(t.status.name, t.assigned) for t in r.by_tunnel])
print("stale lossless", ngc.FlowPolicy(ngc.FlowLedger(g), 0, 3, pinned(prog, black_holes=down))
      .place(float("inf")).call_delivered)
```

```text
12 [(0, 4), (1, 4), (2, 4)]
8 [(0, 4), (1, 4)]
repair 10 [('UP', 10), ('FAILED', 0)]
stale 5 [(1, 5)] [('UP', 5), ('BLACKHOLED', 5)]
stale lossless 1
```

The repaired head sends all 10 over its surviving tunnel. Stale, node 0 still splits 5 and 5 and
the second tunnel's 5 are dropped on edge 1, the edge entering node 2; Lossless admits 1 unit (2
units would put one on the black hole).

A program of a controller replays its solve; under a mask its failed tunnel's share moves to the
head's surviving tunnel:

```python
g = graph(4, [(0, 1, 2, 1), (0, 2, 10, 1), (1, 3, 2, 1), (2, 3, 10, 1)])
solved = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 3, ngc.FlowPolicyConfig(routing=ngc.CentralizedTE(source_allocation=ngc.SourceAllocation.OPTIMIZED)))
prog = solved.replace(12).program
print([(t.entry_id, list(t.edges), t.weight) for t in prog.tunnels], prog.source_weights)
r = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 3, pinned(prog, ngc.Admission.LOSSY),
                   g.eligibility(excluded_edges=[2])).replace(12)
print(r.call_delivered, r.call_edge_drops, [(t.status.name, t.assigned, t.delivered) for t in r.by_tunnel])
```

```text
[(0, [0, 2], 2), (1, [1, 3], 10)] ((0, 12),)
10 [(1, 2)] [('FAILED', 0, 0), ('UP', 12, 10)]
```

The first tunnel is `FAILED`; the head's 12 go to its other tunnel by its weight 10 alone, which
carries 10 and drops 2 on edge 1.

#### Failure studies

A failure has three states, all built from these primitives: **stale forwarding** (no router has
reacted yet: the black holes above), **local repair** (the routers drop the failed routes and split
over the survivors: the mask of a `Pinned` policy) and **converged routing** (a new computation on
the failed network: a new policy with that mask). Save the baseline program before the failure.
`evaluate_edge_failures` recomputes or repairs under each deletion
([section 10](#10-evaluators)); it never turns a deletion into a black hole, so a stale study of
every link is a loop over `Pinned(base, black_holes=...)`. A failure never changes the graph or
another policy's mask; on a shared ledger, the old owner is removed or replaced before its
replacement is modelled, since a new policy is a different demand that does not release the old
one's traffic or holds.

### Losses on shared edges

Where traffic meets on an edge that cannot carry it all (Lossy), the edge carries `T = min(free
capacity, sum of arrivals)` and divides it by highest averages, level by level: first over the
demands by what each demand's traffic brings (one program is one demand; in a SIMULTANEOUS list
each policy is a demand, ties to the one created first), then within a demand over its heads, then
within a head over its entries by their arrivals. Each part continues with its share; the rest is
dropped on that edge. A pooled `ShortestPaths` demand (several sources whose routes meet) splits
its total once at each router, and each member's amount is shared among the sources by highest
averages over what each source still brings. Every demand is evaluated at once: forwarded first
with no edge shared and its traffic dropped only at its own black holes, an edge whose load fits
can never drop and waits for nothing; every other shared edge waits for all its possible arrivals,
and the edges are evaluated in dependency order (the order of independent edges does not change
the result). Traffic that can carry nothing in that first forward (a share of 0, a split member
whose whole-unit share is 0, traffic behind its own demand's black hole) takes no part: it loads
no edge there, makes no edge wait and never keeps two edges in series from passing traffic through
as one, so adding a demand that can carry nothing never changes another demand's result, error or
round count. Splitting an edge into two in series with the same capacity changes nothing (a black
hole on it becomes one on its first part, or on both). Traffic that reaches its own demand's black
hole is dropped there and takes no part in the sharing.

When edges wait on each other in a cycle (one tunnel meets edge X before Y, another Y before X,
and both can drop; the tunnels may belong to different demands), only that group is recalculated, round by round, from nothing carried, until
a round gives back exactly what went in; the answer is checked independently against the rule
before anything is published. A group that alternates between two states raises
`ValueError("lossy forwarding iteration entered a period-2 cycle")`; one that does not settle
within 64 rounds per waiting edge raises `ValueError("lossy forwarding iteration did not reach a
fixed point within the round limit")`; a failed check, an implementation fault, raises
`RuntimeError("lossy forwarding fixed-point verification failed")`. Each error publishes nothing
and does not claim that no answer exists.

```python
# Two tunnels of one head share the second edge of 2 on tunnel 0 only.
g = graph(4, [(0, 1, 10, 1), (0, 2, 10, 1), (1, 3, 2, 1), (2, 3, 10, 1)])
prog = ngc.ForwardingProgram.from_routes(g, 0, 3, [tunnel(g, 0, [0, 2]), tunnel(g, 0, [1, 3])])
r = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 3, pinned(prog, ngc.Admission.LOSSY)).place(12)
print(r.call_delivered, r.call_edge_drops, [(t.assigned, t.delivered, t.total_edge_drops) for t in r.by_tunnel])

# One edge of 9 shared by two tunnels of weights 2 and 1: the offer splits 8 and 4, the edge 6 and 3.
g = graph(2, [(0, 1, 9, 1)])
prog = ngc.ForwardingProgram.from_routes(g, 0, 1, [tunnel(g, 0, [0], 2), tunnel(g, 0, [0], 1)])
r = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 1, pinned(prog, ngc.Admission.LOSSY)).place(12)
print(r.call_delivered, [(t.assigned, t.delivered, t.total_edge_drops) for t in r.by_tunnel])

# Two tunnels cross edges X (e0, 1 unit) and Y (e3, 2 units) in opposite orders, weights 1 and 2.
g = graph(6, [(0, 1, 1, 1), (1, 2, 1000, 1), (1, 5, 1000, 1), (2, 3, 2, 1), (3, 0, 1000, 1),
              (3, 5, 1000, 1), (4, 0, 1000, 1), (4, 2, 1000, 1)])
prog = ngc.ForwardingProgram.from_routes(
    g, 4, 5, [tunnel(g, 4, [6, 0, 1, 3, 5], 1), tunnel(g, 4, [7, 3, 4, 0, 2], 2)])
try:
    ngc.FlowPolicy(ngc.FlowLedger(g), 4, 5, pinned(prog, ngc.Admission.LOSSY)).place(3)
except ValueError as error:
    print(error)
```

```text
8 [(2, 4)] [(6, 2, 4), (6, 6, 0)]
9 [(8, 6, 2), (4, 3, 1)]
lossy forwarding iteration entered a period-2 cycle
```

In the first program the second edge of 2 on tunnel 0 drops 4 of its 6. In the cycle, the rounds
from zero alternate between two states, so the call raises and publishes nothing.

### T1: tunnels set up by their heads

**In plain words.** `DistributedTE(TunnelCount(N))` models distributed RSVP-TE: there is no
controller, and each head (a roster source) sets up its own N tunnels. A design volume D is split
equally over the heads and each head's share equally over its N tunnels; each tunnel asks for
that bandwidth. One tunnel at a time, in rounds over the heads, the head keeps the links with
enough unreserved bandwidth and takes the shortest path: the tunnel comes up and holds exactly its
request on every link of its path, or it stays down and holds nothing. Real traffic then enters
equally at the heads and splits equally over each head's tunnels that are up; a head with no
tunnel up loses its share, and no one moves it elsewhere. Holds are the ledger's second
projection ([section 6](#6-the-ledger)): they decide which tunnels come up, never what traffic
gets through.

**Exact rules.**

- **Sizing.** `H` roster heads in ascending node id, `N = tunnel_count.per_head` (an "up to" count
  is refused). Head `i`'s design share is `RR(D, H, i)` (`plan.design_by_head`), its tunnel `j`'s
  request `RR(share, N, j)`: the tunnels of one head differ by at most one unit. Tunnels are
  numbered in round order (round `j = 0..N-1`, the heads ascending within a round): entry
  identifier `j * H + i`, local index `j`; this is also the setup order. An empty roster gives an
  empty program.
- **Small shares and zero requests.** When a head's share is below N units, some requests are 0.
  All N tunnels are still set up: a tunnel with request 0 comes up holding nothing, like an
  RSVP-TE tunnel with no bandwidth; its search has no bandwidth filter (unless
  `min_unreserved_bandwidth` is positive), so it may cross an edge of capacity 0, where the
  traffic later sent on it is lost on that edge. A design of 0 is the extreme case: every request
  is 0 and an offer is forwarded over the tunnels.
- **The setup view U** of an edge is the ledger's unreserved capacity with this policy's own holds
  released (never above the capacity). Every attempt runs, also after failures, on U minus the
  requests of the tunnels already up in this setup.
- **One attempt (CSPF)** with request `b` keeps the allowed edges with `U >= max(b,
  min_unreserved_bandwidth)` (a request of 0 with a minimum of 0 keeps every allowed edge, one of
  capacity 0 included), then with `parallel_edges=False` one member per pair by the member rule
  (the filter first, then the member rule, at every attempt), then takes the canonical single path
  of [section 4](#the-canonical-single-path) on static costs toward the allowed destinations: the
  lowest cost; under a `HIGHER_*` rule the larger bottleneck score; fewer hops; the smaller EdgeId
  sequence. The scores read the setup view: `HIGHER_AVAILABLE_CAPACITY` reads U,
  `HIGHER_HEADROOM` U divided by the capacity (exactly; 0 on an edge of capacity 0),
  `HIGHER_CAPACITY` the capacity.
- **The ceiling.** `AbsoluteCost(c)` admits a path of cost at most `c`. `CostStretch(n/d)` admits
  a cost `c` exactly when `c * d <= n * c0`, `c0` the smallest static shortest cost over all
  roster heads under the policy's Eligibility, edges of capacity 0 included, before any bandwidth
  filter: one baseline for the whole pool, not per head (a head whose own shortest cost is 3 can
  fail `CostStretch(2)` because another head has a zero-capacity shortest route of cost 1), fixed
  when the policy is built.
- **Statuses.** `UP`: the tunnel holds `b` on every edge of its path (its record's
  `edge_reservations`) and has weight 1. `INSUFFICIENT_BANDWIDTH`: no path passes the filter.
  `COST_CEILING_EXCEEDED`: the cheapest filtered path exceeds the ceiling. A down tunnel keeps its
  request, has no destination, path or cost (`None`), has weight 0 and holds nothing. There is no
  partial hold, retry, growth or count reduction; `min_unreserved_bandwidth` filters searches and
  is never held. A tunnel count is a number of setup attempts per head, not a number of tunnels
  guaranteed up, and two tunnels may take the same links.
- **Forwarding** is the program's: a total `x` enters as `RR(x, H, i)` at head `i`, then round
  robin over that head's up tunnels in local order. A head with no up tunnel keeps its share and
  loses it (Lossy: a source drop; Lossless admits only totals that give it nothing,
  [equal source allocation](#equal-source-allocation-and-a-source-that-cannot-pass-anything)); the
  demand is never renormalised over the heads that have tunnels. Changing only the offer never
  changes the setup. When every tunnel is up and U is at most the residual on every edge, D itself
  is carried losslessly. Holds do not police: a tunnel carries what forwarding gives it within the
  residual, more than its hold when the edges have room, and its traffic meets every other
  demand's on the residual. T1 reports no `total_cost`, `min_cut` or `marginal_cost_curve`.
- **Operations.** No `place` (`DistributedTE does not support incremental placement`, `place(0)`
  too). `replace(offer, design_demand=D)` releases the policy's own holds and traffic, sets up
  again and forwards; `evaluate` computes the same and publishes nothing (no record, no hold);
  `forward(offer)` forwards the installed program and keeps its tunnels and holds (the new records
  carry the same `edge_reservations`). `design_demand=None` means `D = offer`; an infinite offer
  needs a finite design (`DistributedTE requires a finite design demand`). To study a change of
  the offer alone, keep D explicit: `replace(offer + 1)` without it designs again for the new
  offer. A full zero offer such as `replace(0, design_demand=D)` still sets up and books the
  holds. `remove()` releases both projections. T1 is not monotone in D: more design demand can
  bring more tunnels up or fewer.
- **Records and rows.** One record per program entry in identifier order (`entry_id`), a down
  tunnel's included (zero, an empty route, `min_route_cost` `None`). `by_tunnel` shows a down
  tunnel's setup status and `UP`, or `BLACKHOLED` in stale forwarding; a `Pinned` replay of a T1
  program shows `FAILED` for an up tunnel its mask touches, keeps a down tunnel down and holds
  nothing.
- **The plan.** `result.plan` is a `TunnelPlan` of the program after the call:
  `tunnels_per_head` (N), `design_demand` (D), `design_by_head` (`(head, RR(D, H, i))` for every
  roster head) and `total_reserved_bandwidth` (the up tunnels' requests, each counted once).
  T2 also returns a `TunnelPlan`, with the fields described below; all other routings, including
  `Pinned` over a T1 program, return `None`.

A T1 call: head 0's only edge has 2 units, head 1's has 10, two tunnels per head:

```python
g = graph(3, [(0, 2, 2, 1), (1, 2, 10, 1)])
ledger = ngc.FlowLedger(g)
lossy_t1 = ngc.FlowPolicyConfig(routing=ngc.DistributedTE(ngc.TunnelCount(2)),
                                admission=ngc.Admission.LOSSY)
pol = ngc.FlowPolicy(ledger, [0, 1], 2, lossy_t1)
r = pol.replace(12, design_demand=12)
print([t.status.name for t in r.program.tunnels])
print([(t.entry_id, t.head, t.requested_bandwidth, t.reserved_bandwidth, t.edges) for t in r.program.tunnels])
print(r.call_delivered, r.call_source_drops, r.plan.total_reserved_bandwidth, r.plan.design_by_head)
print(ledger.reserved().tolist(), ledger.unreserved().tolist(), ledger.residual().tolist())
r = pol.forward(8)
print(r.call_delivered, r.call_source_drops, r.reservations_changed, ledger.reserved().tolist())
```

```text
['INSUFFICIENT_BANDWIDTH', 'UP', 'INSUFFICIENT_BANDWIDTH', 'UP']
[(0, 0, 3, 0, None), (1, 1, 3, 3, (1,)), (2, 0, 3, 0, None), (3, 1, 3, 3, (1,))]
6 6 6 ((0, 6), (1, 6))
[0, 6] [2, 4] [2, 4]
4 4 False [0, 6]
```

Each head's share is 6, so every tunnel asks for 3. Head 0's edge has 2 units unreserved: both of
its tunnels are down. Head 1's two tunnels come up on edge 1 and hold 6 of its 10. The offer of 12
enters as 6 and 6: head 0 loses its 6, head 1's tunnels carry 3 and 3. `forward(8)` keeps the
tunnels and the holds; head 0 still loses its share, now 4.

### T2: the controller's equal tunnels

**In plain words.** `CentralizedTE(source_allocation=..., tunnel_count=TunnelCount(N))` models a
central controller (a PCE or an SDN controller) that programs, for every head, the same number of
single-path tunnels, all planned with one common whole-unit bandwidth. It first makes that common
bandwidth as large as possible, then (for a finite design volume) no larger than needed, then
picks the cheapest routes, then (with `prefer_headroom=True`) the routes that leave the most spare
capacity, and finally breaks any remaining tie by one fixed unit-by-unit procedure. Real traffic
then splits equally over the heads and over each head's tunnels; whole-unit remainders may use
room beyond the plan. A planned bandwidth is metadata: the ledger holds nothing for it.
`TunnelCount(L, up_to=True)` gives every head one common count `k <= L`: the smallest whose
program carries the live heads' design shares; if none does, the largest live-share score, then
the smallest count. It does not promise whole-demand lossless delivery. A head that cannot reach a destination over free capacity is dead: it keeps
its share and loses it; a weak but live head limits the common size of every tunnel, by design
(heads that must not limit each other are separate demands). A `support` DAG restricts every step
to its edges.

**Exact rules.**

- **Heads.** The roster (`H` heads in ascending node id) is topological, as for T1. A head is
  **live** if it reaches an allowed destination over allowed links with free capacity `V >= 1`,
  `V` the residual with this policy's own traffic released (`plan.heads_with_capacity`); the
  others are **dead**: they keep their roster position and their share, get no tunnel and lose
  the share. `a` live heads, `k` tunnels each (`N`, or the chosen count), `K = a * k` tunnels.
- **On a support.** `CentralizedTE(source_allocation=..., support=dag, tunnel_count=...)` is T2
  with the allowed links restricted to the edges of `dag`, intersected with the policy's
  Eligibility. Every step (roster, live heads, sizing, count choice, both objectives, the
  decomposition, the program's replay basis) uses this one allowed set, and the amounts, plans,
  routes and flags equal those of the same configuration without a support and with the DAG's
  edge set composed into the Eligibility. Parallel members and EdgeIds are kept; the DAG is not
  read as separate routes.
- **Sizing.** Tunnels of a common size `s > 0` exist exactly when `K` path units flow from a
  virtual root (an arc of capacity `k` to every live head) to the destinations with link
  capacities `min(K, floor(V / s))` (paths stop at the first destination).
  `max_bandwidth_per_tunnel` U is the largest integer size that passes (0 when none does);
  `max_fractional_bandwidth_per_tunnel` is the exact fluid optimum `s*` (a `Fraction` `V[e] / j`
  with `j <= K`; `U = floor(s*)`).
- **Working size.** For a finite design D, `planned_bandwidth_per_tunnel` u is `min(U, max(1, max
  over live heads of ceil(RR(D, H, i) / k)))`, the shares taken over the whole roster, dead heads
  included; an unlimited design uses `u = U`. Every tunnel is planned at u;
  `total_planned_bandwidth` is `K * u`. Once U >= 1, every live head keeps all k tunnels even when
  its share is below k units, including a design of 0 (the tunnels are planned at 1 and carry what
  forwarding gives them).
- **Routes.** At u, among the count vectors `x` (the tunnels on each link) that route `k` tunnels
  per live head with `x[e] <= min(K, floor(V[e] / u))`: first the smallest route-cost sum `J =
  sum of cost[e] * x[e]` (`total_tunnel_cost`; idle tunnels count); then, with
  `prefer_headroom=True`, the smallest `headroom_penalty`: each link with static capacity `B > 0`
  charges its `j`-th tunnel `floor(2^32 * (B - V + j * u) / B)` (the occupied fraction of the
  link, rounded down to 2^-32); then the canonical augmentation: `K` single units from zero, each
  on the residual path from the root to the sink of the smallest `(J, Phi)` price, then the fewest
  arcs, then the smallest sequence of arc keys (`(0, EdgeId, direction)` for links, `(1, head,
  direction)` for root arcs, `(2, destination, direction)` for sink arcs; forward 0, reverse 1); a
  reverse arc undoes an earlier unit. Without `prefer_headroom` the second price is 0. `J` and
  `Phi` are exact Python integers and can exceed 64 bits.
- **The program.** `x` is cut, live head by live head in ascending node id, into unit walks along
  the smallest positive EdgeId to a destination; each head's `k` walks are sorted by EdgeId
  sequence (duplicates are separate tunnels). The tunnels are in round order: local index `j`,
  then roster position `i`, with the entry identifier `j * H + i` (`H` counts dead heads too, so
  a dead head leaves holes). Every tunnel is `UP` with weight 1, `planned_bandwidth = u` and
  `requested_bandwidth = reserved_bandwidth = 0`; a plan books no hold.
- **Up to L.** For each `k` in `1..L` the program at `U(k)` scores `G(min(D, F'(k)))`: `G(z)` is
  the part of a round-robin total `z` that the live heads receive, `F'(k)` the largest total whose
  live shares fit on that program (dead shares ignored); an unlimited design scores `G(F'(k))`.
  The smallest `k` whose score is the live part of D wins; if none, the largest score, then the
  smallest `k`; the chosen `k` is then laid out at its own u (`plan.tunnels_per_head`). Delivery
  is not monotone in `k`, so only two proved bounds skip a `k`: with `M = H * k` roster positions
  the score lies between `G(min(D, M * U(k)))` and `G(min(D, M * (U(k) + 1) - 1))`, and `U(k)`
  never grows with `k` (once it is 0, every later count scores 0). "Up to L" finds the smallest
  common count that carries the live heads' design shares (or the best score); it does not
  promise lossless delivery of the whole demand and does not optimise an arbitrary offer.
- **Count outcomes.** An empty roster gives an empty program (`NO_ROUTE` for a positive offer).
  No live head, or live heads with `U = 0` (`insufficient_resolution` True,
  `max_bandwidth_per_tunnel` 0, `max_fractional_bandwidth_per_tunnel` positive), give a program
  **without tunnels that keeps the roster**: `policy.program` is set, `forward` works and every
  head's share is lost at its head (Lossy `COMPLETE`, Lossless `EXHAUSTED` for a positive offer).
  Otherwise every live head has `k` tunnels. No partial program is made for one head; check the
  plan before assuming `len(tunnels) == N * H`.
- **Forwarding** is the program's: a total `x` enters as `RR(x, H, i)` over the whole roster,
  then round robin over the head's `k` tunnels in local order; a dead head's share is lost.
  `max_lossless_offer` F is the program's lossless maximum on the call's view (a dead head's share
  blocks every total that gives it a unit: F is then 0 or a few remainder units). With `u = U` and
  every head live, `0 <= F - K * U < K`.
- **Configuration.** `source_allocation` is stated: `EQUAL` with any number of heads; `OPTIMIZED`
  only with one input source, the same calculation (`CentralizedTE with tunnel_count and several
  sources requires source_allocation=EQUAL`, counted before any Eligibility). `prefer_headroom`
  (`False` by default) switches the headroom objective on; T2 has no route choice of its own:
  it always keeps every parallel link, uses single-path tunnels and the canonical final tie.
  Admission `LOSSLESS` or `LOSSY`. No path-cost ceiling. The tunnel count and the tunnel limit as
  for T1 (`L` counts for "up to").
- **Operations.** No `place` (`CentralizedTE with tunnel_count does not support incremental
  placement`, `place(0)` too). `replace(offer, design_demand=D)` releases the policy's own
  traffic, lays the tunnels out again for D and forwards; `design_demand=None` means `D = offer`,
  and `float('inf')` (or an infinite Lossless offer) is an unlimited design. `evaluate` computes
  the same and publishes nothing; `forward(offer)` forwards the installed program and keeps its
  plan. A planned size change over the same routes sets no flag. T2 reports no `total_cost`,
  `min_cut` or `marginal_cost_curve`.
- **Records and rows.** One record per tunnel in identifier order (a dead head has no entry and no
  record).
- **The plan.** `tunnels_per_head`, `design_demand`, `design_by_head` (empty for an unlimited
  design) and `total_reserved_bandwidth` (0), then T2's own fields: `heads_with_capacity`,
  `insufficient_resolution`, `max_fractional_bandwidth_per_tunnel`, `max_bandwidth_per_tunnel`,
  `planned_bandwidth_per_tunnel`, `total_planned_bandwidth`, `total_tunnel_cost`,
  `headroom_penalty` (`None` without `prefer_headroom`) and `max_lossless_offer`. With no live head
  every amount is 0 (`tunnels_per_head` is N, or 1 for "up to"). A plan is metadata of the
  computation, neither the installed program nor its delivery.

A T2 call: two branches from node 0 to node 3 of 2 and 10 (at q = 2^-12: 8192 and 40960 units),
16 tunnels, `prefer_headroom=True`, then the same program replayed with the 2-unit branch's second
link down:

```python
q = 2.0 ** -12
g = ngc.Graph(4, edge_sources=[0, 0, 1, 2], edge_destinations=[1, 2, 3, 3], edge_capacity_units=ngc.to_units([2, 10, 2, 10], q),
              edge_costs=[1, 1, 1, 1], resolution=q)
te16 = ngc.FlowPolicyConfig(routing=ngc.CentralizedTE(source_allocation=ngc.SourceAllocation.EQUAL,
                                                     tunnel_count=ngc.TunnelCount(16), prefer_headroom=True),
                            admission=ngc.Admission.LOSSY)
pol = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 3, te16)
r = pol.replace(49152)
p = r.plan
print(p.max_bandwidth_per_tunnel, p.planned_bandwidth_per_tunnel, p.max_fractional_bandwidth_per_tunnel, p.max_lossless_offer, r.call_delivered)
print(p.total_planned_bandwidth, p.total_tunnel_cost, p.headroom_penalty, p.heads_with_capacity, p.insufficient_resolution)
print([t.edges for t in r.program.tunnels][:3], [t.planned_bandwidth for t in r.program.tunnels][:3])
fail = g.eligibility(excluded_edges=[2])
replay = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 3,
                        ngc.FlowPolicyConfig(routing=ngc.Pinned(pol.program), admission=ngc.Admission.LOSSY),
                        eligibility=fail)
s = replay.replace(49152)
print([t.status.name for t in s.by_tunnel][:3], s.call_delivered)
```

```text
2925 2925 20480/7 46812 47104
46800 32 73610035200 (0,) False
[(0, 2), (0, 2), (1, 3)] [2925, 2925, 2925]
['FAILED', 'FAILED', 'UP'] 40960
```

The fluid optimum is 20480/7 (a size that fills both branches exactly), so the common whole-unit
size is 2925: 2 tunnels fit the 8192-unit branch and 14 the 40960-unit one. The plan carries 46800;
forwarding the whole offer delivers 47104 because remainder units use the room left on the wide
branch. In the replay the two tunnels on the failed branch are `FAILED` and the head splits its
whole share over the other 14.

A dead head, a head too weak for one unit per tunnel, and "up to":

```python
lossy_t2 = lambda tunnel_count: ngc.FlowPolicyConfig(
    routing=ngc.CentralizedTE(source_allocation=ngc.SourceAllocation.EQUAL, tunnel_count=tunnel_count),
    admission=ngc.Admission.LOSSY)
dead = graph(3, [(0, 2, 0, 1), (1, 2, 10, 1)])
r = ngc.FlowPolicy(ngc.FlowLedger(dead), [0, 1], 2, lossy_t2(ngc.TunnelCount(2))).replace(10)
print(r.plan.heads_with_capacity, [t.entry_id for t in r.program.tunnels], r.plan.planned_bandwidth_per_tunnel, r.call_delivered, r.call_source_drops)
tiny = graph(3, [(0, 2, 1, 1), (1, 2, 100, 1)])
r = ngc.FlowPolicy(ngc.FlowLedger(tiny), [0, 1], 2, lossy_t2(ngc.TunnelCount(16))).replace(200)
print(r.plan.insufficient_resolution, r.plan.max_fractional_bandwidth_per_tunnel, r.program.tunnels, r.call_source_drops)
r = ngc.FlowPolicy(ngc.FlowLedger(tiny), [0, 1], 2, lossy_t2(ngc.TunnelCount(3, up_to=True))).replace(200)
print(r.plan.tunnels_per_head, [t.edges for t in r.program.tunnels], r.call_delivered)
```

```text
(1,) [1, 3] 3 5 5
True 1/16 () 200
1 [(0,), (1,)] 101
```

Head 0's only link has no free capacity: it is dead, keeps its share of 5 and loses it, and its
positions 0 and 2 are holes. Behind links of 1 and 100, 16 tunnels per head fit no whole unit
(`s*` = 1/16): the program has no tunnels and the whole offer is lost at the heads. Up to 3, one
tunnel per head fits a unit; it carries 1 of head 0's 100 and all of head 1's.

The headroom objective on two equal-cost parallel links of 10, two tunnels for a design of 2:

```python
par = graph(2, [(0, 1, 10, 1), (0, 1, 10, 1)])
for on in (False, True):
    cfg = ngc.FlowPolicyConfig(routing=ngc.CentralizedTE(source_allocation=ngc.SourceAllocation.EQUAL,
                                                         tunnel_count=ngc.TunnelCount(2), prefer_headroom=on),
                               admission=ngc.Admission.LOSSY)
    r = ngc.FlowPolicy(ngc.FlowLedger(par), 0, 1, cfg).replace(12, design_demand=2)
    print(on, [t.edges for t in r.program.tunnels], r.plan.headroom_penalty, r.plan.max_lossless_offer, r.call_delivered)
```

```text
False [(0,), (0,)] None 10 10
True [(0,), (1,)] 858993458 20 12
```

Both plans size the tunnels at 1 unit (U = 10). Without the objective the cost tie goes to the
smallest keys: both tunnels on link 0, which then bounds forwarding at 10. With it, a second slot
on link 0 would pay `floor(2^32 * 2 / 10)` against `floor(2^32 / 10)` on link 1, so the tunnels
spread and the program forwards 20.

T2 on a support: three parallel links from node 0 to node 1 of (capacity, cost) (6, 1), (6, 1) and
(100, 2); the support holds the first two. The cheap pair carries one tunnel each and the third
link is never used; with the second link masked, both tunnels share the first:

```python
sup = graph(2, [(0, 1, 6, 1), (0, 1, 6, 1), (0, 1, 100, 2)])
support = ngc.PredecessorDAG.from_csr(sup, [0, 0, 2], [0, 0], [0, 1])
on_support = ngc.CentralizedTE(source_allocation=ngc.SourceAllocation.EQUAL, support=support,
                               tunnel_count=ngc.TunnelCount(2))
for mask in (None, sup.eligibility(excluded_edges=[1])):
    pol = ngc.FlowPolicy(ngc.FlowLedger(sup), 0, 1, ngc.FlowPolicyConfig(routing=on_support), mask)
    r = pol.replace(float("inf"))
    print(r.plan.max_bandwidth_per_tunnel, [t.edges for t in r.program.tunnels], r.plan.max_lossless_offer)
```

```text
6 [(0,), (1,)] 12
3 [(0,), (0,)] 6
```

### The placement call

`place(offer, outputs=None)` takes an offer in units (an `int`) or `float('inf')`.

1. **Validation** happens before anything changes. Within an epoch, a finite offer after an
   infinite one is rejected (`finite placement after infinite placement requires a new epoch`, the
   cumulative request would be undefined), `place(0)` included.
2. **An offer of 0** is `COMPLETE` and does nothing: it creates no flow and prepares nothing.
   Before any program exists (a new policy, or after `remove()`) its result is all zeros, every
   flag false; with a program installed it reports the demand's current state with every flag
   false. It is not a removal: `remove()` is.
3. **`ShortestPaths` split routing and `Pinned`.** `ShortestPaths` prepares the program at the
   epoch's first positive call (one record per roster source, also at zero delivery); `Pinned`
   creates its records, one per entry, there. Lossless then admits the largest total `a` between
   what it already delivers and that plus the offer whose forwarding of `RR(a)` from zero fits the
   free capacity plus the policy's own allocation with no drop and gives no share to a source
   without a member (exact binary search: the totals that fit form a prefix); a unit refused by
   an earlier call is not offered again (the next call targets what is delivered plus its own
   offer). Lossy injects at the sources this call's increase of the injection,
   `RR(offered + v) - RR(offered)` (or the difference of the highest averages over captured source
   weights), and forwards it hop by hop on the free capacity: each edge carries the smaller of
   what reaches it and its residual and drops the rest; old drops are never offered again. An
   infinite Lossy offer is rejected for `Pinned` (`LOSSY forwarding requires a finite offer`); for
   `ShortestPaths` it is the first-hop experiment of [section 5](#5-route-allocations-and-admission).
4. **The controllers** re-solve as described above.
5. **T1 and T2** have no `place`: every call raises `DistributedTE does not support incremental
   placement` or `CentralizedTE with tunnel_count does not support incremental placement`, before
   anything else (an offer of 0 too); their calls are `replace`, `evaluate` and `forward`.
6. **Commit.** All changes to the ledger and the policy are committed together at the end, or
   none are ([section 8](#8-failures-threads-and-lifetimes)).

What the call guarantees:

- `call_delivered` is the increase in delivery made by this call, and `0 <= call_delivered <=
  offer`. For Lossy it is this call's delivery. `call_undelivered = offer - call_delivered`
  (`float('inf')` for an infinite offer; for Lossy it includes what was lost); it is never a
  backlog to retry.
- A `place` call never reduces `delivered()`, the sum of `delivered` over the policy's records.
- The first positive offer creates the policy's records even when none of them can deliver
  anything yet. Such a zero-delivery record keeps its FlowId for the rest of the epoch and can fill
  on a later call.
- **Additivity.** With the rest of the network unchanged between the calls, `place(a); place(b)`
  equals `place(a + b)` exactly for Lossless `ShortestPaths` split routing, the controllers (value,
  allocation and cost curve) and Lossless `Pinned`. Lossy calls are not additive: each forwards its
  own injection on the free capacity the earlier calls left, and nothing is offered again.

### Epochs

`remove()` takes every flow of the policy out of the ledger (allocations reverted, records gone,
ids retired), clears the installed program and starts a new epoch: the cumulative request, the
round-robin position of the source split and the finite-or-infinite state are reset. The next
positive offer prepares again (a `ShortestPaths` policy with a static route choice uses its kept,
equal program; see Routing models) and creates new flows with new ids; ids are never reused.
Until the next positive `place` or a `replace`, `program` reads empty, also when the policy keeps
a program. A repeated `remove()` does nothing. A policy destroyed without `remove()` leaves its traffic and
holds in the ledger as orphans until `FlowLedger.reset()`; nothing is redistributed.
`replace` and `forward` also start a new epoch (below).

```python
g = graph(4, [(0, 1, 10, 1), (1, 3, 10, 1), (0, 2, 2, 1), (2, 3, 2, 1)])
ledger = ngc.FlowLedger(g)
pol = ngc.FlowPolicy(ledger, [0, 1], 3, config(ngc.ShortestPaths()))
for v in (1, 1, 0):
    r = pol.place(v)
    print(v, r.call_delivered, r.routes_changed, [(x.flow_id[1], x.delivered) for x in pol.records()])
pol.remove()
r = pol.place(1)
print(1, r.call_delivered, r.routes_changed, [(x.flow_id[1], x.delivered) for x in pol.records()])
```

```text
1 1 True [(0, 1), (1, 0)]
1 1 False [(0, 1), (1, 1)]
0 0 False [(0, 1), (1, 1)]
1 1 True [(2, 1), (3, 0)]
```

The first unit goes to source 0 and the second to source 1 (round robin of the cumulative
offer); the first call creates both records, source 1's at zero delivery. After `remove()` the
split starts again at source 0, under new ids.

### Operations: evaluate, replace and forward

Besides `place`, a policy has three full-offer calls. Each starts from the policy's own **released
view**: the ledger's residual with the policy's own traffic given back.

| Call | What it does | Routings |
| --- | --- | --- |
| `evaluate(offer, *, design_demand=None, outputs=None)` | computes and forwards `offer` from the released view as a fresh epoch would (T1: sets up on its unreserved view with its own holds released), and reports it; publishes nothing, creates no record and books no hold | all |
| `replace(offer, *, design_demand=None, outputs=None)` | the same computation, then published at once: the old records and holds are retired, new ones created, the new program installed, a new epoch begun with `requested = offer`; delivery may go down | all |
| `forward(offer, *, black_holes=None, outputs=None)` | forwards a full new offer over the installed program with no route computation and no setup (weights stay frozen, a controller does not re-solve, T1 keeps its tunnels and holds, T2 its tunnels and plan), on the released view, then publishes like `replace`; `black_holes` forwards it stale for this call | every computing model with an installed program; `Pinned` uses `replace` |

- `call_delivered` of a full-offer call is its whole new delivery and `delivery_delta` the new delivery
  minus the old (it may be negative). A later `place` continues the epoch a `replace` or a
  `forward` began (its records keep their ids).
- With other owners fixed, identical full inputs give identical results, except fresh record ids.
- `forward`'s checks, in order: `Pinned forwards with replace`; the offer; a Lossy `+inf`
  (`LOSSY forwarding requires a finite offer`, as for `evaluate` and `replace` and Lossy `Pinned`);
  the black holes' type and graph (`black_holes belongs to another graph`); then, under the ledger
  gate, `forward requires an installed program` (also after `remove()`).
- **`design_demand`** belongs to T1 and T2 alone (every other routing raises `<model> does not
  accept design_demand`, the model `ShortestPaths`, `CentralizedTE` or `Pinned`): an integer in
  `[0, 2^62)` or `float('inf')` (else `design_demand must be an integer in [0, 2^62) or permitted
  +inf`); `None` means equal to the offer; T1 needs it finite (`DistributedTE requires a finite
  design demand`, also for an infinite offer without a design); T2 takes an unlimited one (its
  largest common size). It never accumulates. `evaluate` and `replace` check the offer (`offer
  must be an integer in [0, 2^62) or permitted +inf`), then a Lossy `+inf`, then the design (the
  routing, its domain, finiteness).
- A program reused elsewhere, or read again under a different mask, goes through a new `Pinned`
  policy; a program forwarded by `Pinned` never books holds, whatever its producer (read a
  result's `edge_reservations` or the ledger for actual holds).

```python
g = graph(2, [(0, 1, 10, 1), (0, 1, 10, 1)])
ledger = ngc.FlowLedger(g)
pol = ngc.FlowPolicy(ledger, 0, 1, config(ngc.ShortestPaths()))
r = pol.replace(6)
print(r.call_delivered, [x.flow_id[1] for x in pol.records()])
r = pol.evaluate(12)
print(r.call_delivered, r.delivery_delta, ledger.carried().tolist())
r = pol.replace(4)
print(r.call_delivered, r.delivery_delta, [x.flow_id[1] for x in pol.records()])
r = pol.forward(7)
print(r.call_delivered, r.routes_changed, r.weights_changed, ledger.carried().tolist())
```

```text
6 [0]
12 6 [3, 3]
4 -2 [1]
7 False False [4, 3]
```

### Termination

The call's `termination` (of `place` and of the full-offer calls alike) is the first row that
applies (an invalid call raises instead):

| Termination | When |
| --- | --- |
| `COMPLETE` | the offer is 0 |
| `NO_ROUTE` | the offer is positive and the roster is empty (`Pinned`: the replay roster) |
| `NO_ACTIVE_WEIGHT` | the program has an entry with a path (an up tunnel, a controller path or a split route), but nothing can be sent anywhere: every `WEIGHTED` member has weight 0, or no head with a positive injection weight has a surviving entry of positive weight (routes of weight 0, every up tunnel failed). A program whose tunnels are all down (T1) or that has no tunnels (T2 with every head dead, or U = 0) is not this case: it goes on to the rows below, every share lost at its head |
| `COMPLETE` | Lossy with a finite offer: the whole offer was forwarded, and delivery and losses account for all of it (it does not mean nothing was lost) |
| `COMPLETE` | Lossless: the whole offer was admitted |
| `EXHAUSTED` | otherwise: Lossless below the offer, any infinite offer (the first-hop experiment included), a controller whose maximum is 0 over existing routes |

`call_undelivered` is an amount, never a status: a Lossy call can be `COMPLETE` with a positive
`call_undelivered`. The live WEIGHTED controller has one additional outcome: positive offer,
nonempty roster and all roster weights zero gives NO_ACTIVE_WEIGHT, zero delivery and all finite
traffic unadmitted. It installs an empty program with captured zero weights. Pinned replay follows
its fixed-program rule: LOSSY COMPLETE with the whole offer lost at the sources, LOSSLESS
EXHAUSTED. This is the one producer/replay termination difference; empty roster is still NO_ROUTE
and zero offer COMPLETE.

```python
def show(r):
    print(r.call_delivered, r.termination.name, [(row.source, row.delivered) for row in r.by_source])

# Every tunnel down: five edges, the design 14 asks 7 per tunnel and no path has 7.
g = graph(4, [(0, 1, 6, 1), (0, 2, 6, 3), (1, 2, 4, 1), (1, 3, 6, 3), (2, 3, 6, 1)])
for admission in (ngc.Admission.LOSSY, ngc.Admission.LOSSLESS):
    all_down = ngc.FlowPolicyConfig(routing=ngc.DistributedTE(ngc.TunnelCount(2)), admission=admission)
    r = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 3, all_down).replace(14, design_demand=14)
    print(admission.name, [t.status.name for t in r.by_tunnel], r.call_source_drops, end=" ")
    show(r)

# One edge for source 1 to node 4.
g = graph(5, [(0, 2, 1, 1), (0, 3, 1, 1), (2, 4, 1, 1), (3, 4, 1, 1), (1, 4, 1, 1)])

# Weights from the residual where another demand holds the only edge: no positive weight.
ledger = ngc.FlowLedger(g)
ngc.FlowPolicy(ledger, 1, 4, config(ngc.ShortestPaths())).place(1)
weighted = ngc.FlowPolicyConfig(routing=ngc.ShortestPaths(allocation=ngc.RouteAllocation.WEIGHTED,
                                                          weight_basis=ngc.WeightBasis.RESIDUAL_SNAPSHOT))
show(ngc.FlowPolicy(ledger, 1, 4, weighted).place(1))
```

```text
LOSSY ['INSUFFICIENT_BANDWIDTH', 'INSUFFICIENT_BANDWIDTH'] 14 0 COMPLETE [(0, 0)]
LOSSLESS ['INSUFFICIENT_BANDWIDTH', 'INSUFFICIENT_BANDWIDTH'] 0 0 EXHAUSTED [(0, 0)]
0 NO_ACTIVE_WEIGHT [(1, 0)]
```

Both of T1's tunnels are down, so the program has no entry with a path: Lossy forwards the whole
offer and loses it at the head (`COMPLETE`), Lossless admits nothing (`EXHAUSTED`). The Weighted
route has an entry with a path whose only member has weight 0: `NO_ACTIVE_WEIGHT`.

### Flags

`routes_changed`, `weights_changed` and `reservations_changed` start false in every call and
record only changes that were committed; a call that raises publishes nothing (`evaluate` reports
what `replace` would set). They compare the installed program before and after the call (no
program counts as an empty one): `routes_changed` the ordered heads, each entry's kind and route
(a tunnel's path, a split route's members), duplicates included; `weights_changed` the injection
(equal or the captured weights, by value) and every entry and split-table weight, an entry absent
from one program counting as weight 0; `reservations_changed` every tunnel's `requested_bandwidth`
and `reserved_bandwidth` by identifier, an absent tunnel counting as 0 (only T1 sets it: T2's planned
bandwidth is no reservation, and a T2 layout over the same routes at another planned size sets
no flag). Hence:

| Committed change | `routes_changed` | `weights_changed` |
| --- | --- | --- |
| the epoch's first nonempty program (`ShortestPaths` and the controllers' first positive call, a `replace`, `replace(0)` included) | true | true when some entry has a positive weight (a T1 program whose tunnels are all down: false) |
| a T1 setup whose paths or up and down tunnels differ from the installed one's (requests changed too: `reservations_changed`) | true | true when a tunnel came up or went down |
| an empty program over no program or an empty one (a call with an empty roster, `replace(0)` included) | false | false |
| a controller's re-solve that moves a source's share to other edges, also with `call_delivered == 0` | true | true unless the entries keep their number and amounts |
| a controller's re-solve that changes only a path's amount (`CentralizedTE` `place(5)` then `place(5)` on one path) | false | true |
| an offer-only call over an unchanged program: `forward`, `Pinned`, `ShortestPaths` increments, an identical `replace` (new ids) | false | false |
| nothing of the above: `place(0)` (an offer of 0, which prepares and installs nothing, unlike a full zero offer such as `replace(0)`, which follows the rows above), rejected candidates only | unchanged | unchanged |

`reservations_changed` is true exactly when some tunnel's requested or reserved amount changed by
identifier: a T1 `replace` (or `evaluate`) whose design or setup differs, even over unchanged
paths; it is never set by `forward`, by `Pinned` or by the other routings. The flags describe
routes, weights and holds only: delivered and dropped amounts can change with all three false.

```python
# A first offer on a route whose only unit another demand holds.
g = graph(5, [(0, 2, 1, 1), (0, 3, 1, 1), (2, 4, 1, 1), (3, 4, 1, 1), (1, 4, 1, 1)])
ledger = ngc.FlowLedger(g)
ngc.FlowPolicy(ledger, 1, 4, config(ngc.ShortestPaths())).place(1)
pol = ngc.FlowPolicy(ledger, 1, 4, config(ngc.ShortestPaths()))
r = pol.place(1)
print(r.call_delivered, r.termination.name, r.routes_changed, r.weights_changed,
      [x.delivered for x in pol.records()])
```

```text
0 EXHAUSTED True True [0]
```

### The result

`place`, `evaluate`, `replace` and `forward` (and every evaluator) return a `PlacementResult`. Its
fields describe either this call or the demand's state after it (totals over the epoch); a
per-call loss field carries the `call_` prefix, a state field does not, so the two are never
confused.

- **This call**: `call_delivered`, `call_undelivered`, `termination`, `routes_changed`,
  `weights_changed` and `reservations_changed` as above; `delivery_delta`, the change of the
  delivery (for `place` it equals `call_delivered`); `call_edge_drops`, this call's drops per edge
  as `[(EdgeId, units)]`; `call_source_drops`, the units this call lost at sources.
- **The demand's state**: `program` (the installed `ForwardingProgram`, `evaluate`'s computed one;
  `None` before the epoch's first preparation and after `remove()`), `plan` (T1 and T2: the
  `TunnelPlan` of that program, [T1](#t1-tunnels-set-up-by-their-heads),
  [T2](#t2-the-controllers-equal-tunnels); `None` for every other routing, `Pinned` included),
  `requested` (the sum of the offers, `float('inf')` after an infinite one), `offered`,
  `admitted`, `delivered`, `unadmitted`, `source_drops`, `edge_drops` (per edge, as
  `[(EdgeId, units)]`), `by_source`, `by_tunnel` and `total_cost`. Every amount is an exact Python
  integer and may exceed 64 bits.
- **Accounting.** Lossless: `offered` is `None`, `source_drops` 0, no edge drops, `admitted =
  delivered`, and `requested = unadmitted + delivered` for a finite offer (`unadmitted` is `None`
  after an infinite one). Lossy: `offered` is the traffic actually sent, `admitted = offered`,
  `unadmitted = 0`, and `offered = source_drops + sum(edge_drops) + delivered`; `offered =
  requested` except in the first-hop experiment, where `requested` is infinite and `offered` the
  finite injection. `delivered + sum(edge_drops)` is what entered the network. `edge_flows` are
  carried loads per edge: summing them counts every hop, not the delivery.
- **`by_source`** is a tuple of `SourceFlow` rows, one per input source in ascending node id, a
  source outside the roster included (a zero row): `source`, `offered` (`None` when Lossless),
  `admitted`, `delivered`, `source_drops` and `total_edge_drops` (its share of `edge_drops`). The
  rows sum to the totals, with one exception: a Lossy positive offer with an empty roster reports
  `offered = requested` and all of it as `source_drops`, lost before entering, while every row
  stays zero.
- **`by_tunnel`** is a tuple of `TunnelFlow` rows, one per program tunnel in ascending `entry_id`
  (empty without tunnels): `entry_id`, `head`, `status` (`UP`; `INSUFFICIENT_BANDWIDTH` or
  `COST_CEILING_EXCEEDED` for a tunnel T1's setup left down; `FAILED` or `BLACKHOLED` in this
  forwarding), `assigned` (its injected share), `delivered` and `total_edge_drops`. It is another
  view of the same traffic as `by_source`, never added to it; a head's source drops are not
  attributed to one of its down tunnels.
- **`total_cost`**: the controllers only, the exact `sum(cost[e] * flow[e])` of the policy's
  allocation; `None` otherwise (T1, T2 and `Pinned` included).
- **`marginal_cost_curve`**: unrestricted `CentralizedTE` with `OPTIMIZED` source allocation
  only, the marginal cost curve of the policy's cumulative flow (an offer of 0 reports it
  unchanged; `None` before the epoch's first positive call). `None` for every other model,
  `EQUAL` source allocation included: the marginal costs of an equal flow need not increase. For
  other models, the sum of `(min_route_cost, delivered)` over the records is a histogram of
  delivery by each record's cheapest walk cost: the cost every unit travels at for `ShortestPaths`
  split routing, T1 and T2, whose routes have one cost, but not for a split route whose walks
  differ in cost.
- **Optional outputs**, `None` unless requested through `OutputOptions(edge_flows=False,
  residual=False, saturated_boundary_edges=False, edge_reservations=False, min_cut=False)`:
  `edge_flows` is the policy's own allocation per edge after the call (the ledger's total is
  `carried()`); `residual` and `saturated_boundary_edges` are the ledger's after the call;
  `edge_reservations` is the policy's own holds per edge after the call (zeros for a routing
  without holds; the ledger's total is `reserved()`); `min_cut` is described below.
- **`saturated_boundary_edges` and `min_cut` are different things.** `saturated_boundary_edges` is
  the saturation report of [section 6](#6-the-ledger), available in every model. `min_cut`
  (requested by `OutputOptions(min_cut=True)`) exists only for unrestricted `CentralizedTE` with
  `OPTIMIZED` source allocation: the edges of a minimum cut of the policy's view, the one on the
  sources' side: the allowed edges that leave the nodes the sources can still reach once the flow
  is maximal (counting the option to undo flow on an edge), in ascending EdgeId, an edge with
  nothing free in the view (capacity 0, or held by other policies) included: it adds 0 to the
  cut's value and is listed because every allowed edge leaving that set is listed; an edge the
  policy's Eligibility excludes is never listed; `[]` for an empty roster. It is present only when
  the solve reached the maximum on its view; when an offer below the maximum stopped it, on an
  offer of 0 (nothing is solved), and in every other model, `EQUAL` sources and a supplied support
  included, `min_cut` is `None`.

```python
g = graph(4, [(0, 1, 10, 1), (1, 3, 10, 1), (0, 2, 5, 2), (2, 3, 5, 2)])
ledger = ngc.FlowLedger(g)
pol = ngc.FlowPolicy(ledger, 0, 3, ngc.FlowPolicyConfig(routing=ngc.CentralizedTE(source_allocation=ngc.SourceAllocation.OPTIMIZED)))
want = ngc.OutputOptions(saturated_boundary_edges=True, min_cut=True)
r = pol.place(12, want)
print(r.call_delivered, r.marginal_cost_curve, r.min_cut, r.saturated_boundary_edges, r.termination.name)
r = pol.place(float("inf"), want)
print(r.call_delivered, r.marginal_cost_curve, r.min_cut, r.saturated_boundary_edges, r.termination.name)
print(r.requested, r.admitted, r.unadmitted, r.total_cost)
```

```text
12 [(2, 10), (4, 2)] None [0] COMPLETE
3 [(2, 10), (4, 5)] [0, 1] [0, 1] EXHAUSTED
inf 15 None 40
```

After 12 units the cheap path (cost 2) is full and the solve stopped below the maximum, so there
is no cut; `saturated_boundary_edges` already reports the full edge 0. The second call adds 3 units at cost 4
and reaches the maximum of 15, at a total cost of 10 x 2 + 5 x 4 = 40.

### Several demands at once

    apply_actions(ledger, actions, *, mode=BatchMode.ORDERED, outputs=None) -> list[PlacementResult]
    Place(policy, offer)
    Replace(policy, offer, *, design_demand=None)
    Forward(policy, offer, *, black_holes=None)

`apply_actions` runs the actions of several policies of one ledger as one operation and returns
one `PlacementResult` per action, in input order, each with the fields of that policy's scalar
call (`place`, `replace` or `forward`). The whole list is validated first and published once: if
anything fails, the ledger and every policy are as they were. A policy appears at most once, and
every policy belongs to `ledger`. An action is an immutable value that keeps its policy (and its
black holes) alive; its offer and design demand are kept as given (an int or a float: another type
raises `TypeError` when the action is built), and the call checks their domain.

- **`BatchMode.ORDERED`** runs each action in list order exactly as its scalar call would, on the
  state the earlier actions left (their traffic and their holds); the policies not reached yet keep
  their old state. Normal outcomes (exhaustion, source drops, a down tunnel) are kept. Earlier
  actions are served first: on one edge of 10, an ECMP replace of 10 and then a T1 tunnel's replace of 10
  deliver 10 and 0, though the tunnel holds its 10, because a hold reserves setup bandwidth, not
  service. New records take FlowIds in list order.
- **`BatchMode.SIMULTANEOUS`** forwards every action at the same time on fixed routes, the losses
  shared across the demands by [the loss rule](#losses-on-shared-edges). It takes `Replace` of
  `ShortestPaths` `EQUAL` or `WEIGHTED` (its routes prepared on the list's view) and of `Pinned`,
  and `Forward` of any policy with an installed program; every participant is LOSSY with a finite
  offer. Controllers compute nothing here and T1 sets nothing up: prepare them first (scalar calls
  or an ORDERED list), then forward them; a T1 owner keeps its holds. The view is the residual
  plus every participant's old traffic, released together: owners outside the list and orphans
  stay as fixed load, and holds are never added. Repeating a list with the same outside load gives
  the same state, and permuting it changes nothing. Each participant's black holes are its own:
  what reaches one is dropped there and takes no part in sharing that edge, while the others share
  it as usual. Removing a participant later frees only its own traffic; the others do not
  rebalance until they are forwarded again.
- **Outputs.** `outputs` applies to every entry. `residual` and `saturated_boundary_edges` read the
  ledger after that action (ORDERED) or after the whole list (SIMULTANEOUS: one residual for every
  entry, each entry's boundary from its own sources and Eligibility).
- **Validation, in this order:** the arguments' types (`TypeError`: `ledger must be a FlowLedger`,
  `actions must be a sequence of Place, Replace or Forward (got <type>)`, `actions[<i>] must be a
  Place, Replace or Forward (got <type>)`, `mode must be a BatchMode (got <type>)`, `outputs must be
  an OutputOptions (got <type>)`); the mode's value (`mode: <v> is not a value of BatchMode`); an
  empty list checks ledger usability under the gate before returning `[]`. Each action's
  argument checks run in list order: `list policy belongs to another ledger`, `duplicate policy in list`, the SIMULTANEOUS
  rules (`SIMULTANEOUS does not accept Place`, `SIMULTANEOUS Replace requires ShortestPaths EQUAL or
  WEIGHTED, or Pinned`, `SIMULTANEOUS requires admission=LOSSY`, `SIMULTANEOUS requires finite
  offers`), then the scalar call's argument checks. After all actions pass, check ledger usability
  and every action's mutable state under the gate, in list order. Only then compute anything. An
  error of the computation (a loss cycle that does not settle, a solver or forwarding size limit,
  an allocation failure) publishes nothing.

```python
# Two ECMP demands of 10 on one edge of 10, at the same time and one after the other.
g = graph(2, [(0, 1, 10, 1)])
ecmp = ngc.FlowPolicyConfig(routing=ngc.ShortestPaths(), admission=ngc.Admission.LOSSY)
for mode in (ngc.BatchMode.SIMULTANEOUS, ngc.BatchMode.ORDERED):
    led = ngc.FlowLedger(g)
    a = ngc.FlowPolicy(led, 0, 1, ecmp)
    b = ngc.FlowPolicy(led, 0, 1, ecmp)
    res = ngc.apply_actions(led, [ngc.Replace(a, 10), ngc.Replace(b, 10)], mode=mode)
    print(mode.name, [r.call_delivered for r in res], [r.call_edge_drops for r in res])

# Edges 0: 0 -> 2 of 5, 1: 1 -> 2 of 10, 2: 2 -> 3 of 10; a from 0 and b from 1, 5 each.
g = graph(4, [(0, 2, 5, 1), (1, 2, 10, 1), (2, 3, 10, 1)])
led = ngc.FlowLedger(g)
a, b = (ngc.FlowPolicy(led, s, 3, ecmp) for s in (0, 1))
out = ngc.OutputOptions(residual=True, saturated_boundary_edges=True)
for r in ngc.apply_actions(led, [ngc.Replace(a, 5), ngc.Replace(b, 5)],
                           mode=ngc.BatchMode.SIMULTANEOUS, outputs=out):
    print(r.call_delivered, r.residual.tolist(), r.saturated_boundary_edges)
```

```text
SIMULTANEOUS [5, 5] [[(0, 5)], [(0, 5)]]
ORDERED [10, 0] [[], [(0, 10)]]
5 [0, 5, 0] [0]
5 [0, 5, 0] [2]
```

At the same time the edge shares its 10 between the two demands; one after the other the first
takes it all. In the second list both entries report the residual after the whole list, and each
its own boundary: a's source reaches nothing beyond its full edge 0, b's reaches node 2, whose exit
is full.

## 8. Failures, threads and lifetimes

**In plain words.** A call either does everything or changes nothing: if it fails, for any
reason, every ledger, policy and cache is as it was and you get one exception. There are two
stated exceptions to this, described below. Immutable values can be shared freely between
threads. Calls on one ledger serialize; different ledgers run in parallel. Cache leases serialize
building calls and admit frozen readers together.

**Exact rules.**

- **All or nothing.** A failure anywhere in a public call (validation, memory allocation, a
  worker thread, the construction of the Python result) leaves every ledger, policy, cache and
  internal buffer as it was before the call and raises one exception in the caller's thread.
  Worker threads are always joined first; an exception never escapes a worker. When several pairs
  or candidates of an evaluator fail, the exception of the lowest index is raised.
- **Limit 1: pybind11's own object machinery.** Constructing a bound object from Python (any of
  the library's classes), parsing the arguments of a call that passes keywords, and converting the
  result of a method into a new Python object outside the library's checked construction (for
  example `PredecessorDAG.from_csr`, `Graph.eligibility`) use pybind11's own allocation paths,
  which do not all check for failure. A memory failure at exactly that point may crash the process
  instead of raising. Library state is never left half-changed by it. The results of `spf`, `ksp`,
  `resolve_to_paths`, `ForwardingProgram.from_routes` and `from_dag`, `FlowPolicy.place`,
  `evaluate`, `replace` and `forward`, `apply_actions`, `evaluate_demand`, `evaluate_pairs` and
  `evaluate_edge_failures` are built through a checked path and stay inside the guarantee, and so
  is every read of a field or property and every `repr` of the library's values (the configuration
  values and their enumerations, `FlowPolicy.config`, `program`, `plan`, `by_source`, `by_tunnel`,
  `tunnels`, `termination`, a tunnel's `status`, a record's `route`, `SensitivityResult.baseline`
  and `entries`, among others), and `records()` and `record()` of a policy and a ledger,
  `FlowLedger.graph` and `Eligibility.intersect`.
- **Limit 2: a `ShortestPathCache` passed to the call.** A building cache is a memo: when a call fails
  after some cache lookups succeeded (a policy constructor failing after its lookup, an evaluator
  failing midway), the cache keeps those complete entries and their statistics. It never holds a
  partial entry, every value it handed out stays valid, and a kept entry never changes any later
  result. A failed lookup counts one miss. A frozen cache never changes its entries; a failed
  call changes none of its statistics, except one miss when the call raised `ShortestPathCacheMiss`.
- **Threads and the GIL.** Each ledger has one native gate protecting its projections and every
  bound policy, including reads, lists, removal and policy destruction. Calls on one ledger take
  effect one at a time; different ledgers can compute in parallel. Native work releases the GIL,
  including policy calls, searches, factories and array conversions. Immutable reads keep the GIL
  and take no gate: graph, mask, route, program, configuration and result fields, `policy.config`,
  cache `graph_uid()` and `stats()`. Each statistics field is read atomically, without promising a
  consistent snapshot of all fields.
- **Call phases.** Convert and copy inputs, pin dependencies and check every argument and immutable
  object in the first pass's order. Python conversions and protocols keep the GIL; linear checks of
  a just-converted argument may keep it too. Native checks that grow with the graph release it,
  before any lease or gate: the tunnel-limit roster search in policy construction and the evaluators
  (once per pair in a batch), and `from_routes`' walk of each supplied PredecessorDAG. `evaluate_pairs`'
  workers and the other evaluators' main evaluation reuse the checked rosters; each deletion of
  `evaluate_edge_failures` searches under its own mask. Publishing calls also allocate and register their fixed Python
  results here: a policy, a `PlacementResult`, a list of n results, or a `SensitivityResult`.
  Release the GIL, acquire at most one cache lease, then the ledger gate. Check mutable state,
  compute privately, join workers and publish without allocation or failure. Frozen hits publish
  last, inside the lease. Release every gate and lease before reacquiring the GIL to return or
  translate an exception. Reads, searches and factories wrap their owned native snapshots after
  unlocking. Python conversions and protocols execute with no native lock held.
- **Validation passes.** In the documented validation order, all argument and immutable checks
  precede all mutable-state checks. Thus an invalid offer or foreign black holes precede
  `forward requires an installed program`; a later invalid list action precedes an earlier
  action's epoch error. Installed programs, finite `place` calls after infinity, requested limits,
  cache mode, ledger usability and identity availability are checked under the gate or lease,
  in their documented order, action by action. No list action computes before both passes finish.
- **Free-threaded Python.** The cp314t build supports the same concurrency rules. It detaches its
  Python thread state before waiting for a gate or lease. Step-one conversion still copies every
  borrowed array; changing an argument concurrently with its conversion remains a caller race.
- **Lifetimes.** The graph outlives the ledger, the ledger outlives its policies; Python keeps
  the graph, ledger, Eligibility and cache alive while a policy needs them.
- **Immutable values.** The configuration values (`FlowPolicyConfig`, the four models,
  `TunnelCount`, `AbsoluteCost`, `CostStretch`, `OutputOptions`, `EvaluationOptions`,
  `RouteDefinition`), the search options (`SpfOptions`, `KspOptions`) and the actions of a list
  (`Place`, `Replace`, `Forward`) have read-only fields: no option object changes after it is
  built. A policy validates its configuration once and keeps that snapshot (`policy.config`); an
  evaluator validates and copies its options at call entry, and `apply_actions` its actions.
- **Exception classes.** One rule for the class of an input error, at every public entry
  (constructors, search options, call arguments, evaluator options): a Python type the position
  never takes is a `TypeError` (a bool, Python's or `numpy.bool_`, at any numeric position; a float
  where only an integer is taken, `k=1.0` and `k=float('inf')` included; a string; `None` where none
  is allowed; an object of another class; for a vector a non-sequence, a string, a nested entry, a
  NumPy array of another rank or of a non-integer dtype, or an entry of a wrong type); an accepted
  type whose value is outside its domain is a `ValueError` (out of range, beyond 64 bits, a wrong
  length, NaN, `-inf`, a finite float where floats are taken only for `+inf`: an offer, a design
  demand). The text is the same either way, and a text that covers both causes (`offer must be an
  integer in [0, 2^62) or permitted +inf`, `a FlowId is a tuple (ledger uid, id)`, `each pair must
  be (sources, destinations)`) takes the class of its cause. A NumPy integer is an integer, a NumPy
  floating value a float; an int at or above `2^63` is an integer (a `ValueError` where its range is
  bounded); a `fractions.Fraction` is taken only by `CostStretch`. `ValueError` also for an invalid
  combination, a foreign object (a value of another graph, ledger or policy) or a call in the wrong
  state;
  `TypeError` for a field or option that does not exist (for example
  `ShortestPaths(tunnel_count=ngc.TunnelCount(2))`), a missing required argument, `None` where
  none is allowed, an argument of the wrong type (a bool field given anything but a `bool` or
  `numpy.bool_`; an integer field or `edge_weights` of a configuration value given a non-integer), or an argument
  that cannot be converted to its parameter's type. The non-integer number positions (a
  `resolution` of `Graph`, `to_units` and `as_float`; `to_units`' values) take exactly an `int`
  (not a bool), a `float`, a NumPy integer or a NumPy floating scalar, an `int` or `float` subclass
  by its value (never through its own `__int__`, `__float__` or `__array__`), and for `to_units`
  an array of an integer or floating dtype or a sequence of those scalars: a `Decimal`, a
  `Fraction`, an object known by `__float__`, `__index__` or `__array__`, and a 0-d array as a
  resolution are a `TypeError`. Every integer position, `as_float`'s amounts included, takes any
  object with `__index__` (never a bool), read exactly. An integer beyond the float range
  as a resolution is a value outside its domain (`ValueError: resolution must be a power of two in
  [2^-32, 2^32]`); `to_units` refuses a bool entry beside numbers (`TypeError`, which NumPy would
  read as a number; an all-bool sequence keeps its dtype text), and an integer too large for NumPy
  beside numbers is a value (`ValueError`). Every converted value passes on unchanged any exception raised by the value's own conversion (its `__float__` or
  `__index__`), whatever its class. `ShortestPathCacheMiss` for a miss in a frozen cache (its own
  class, neither a `ValueError` nor a `RuntimeError`); `MemoryError` for an allocation failure
  (everything rolled back); `RuntimeError` for an internal invariant failure.

```python
g = graph(4, [(0, 1, 10, 1), (1, 3, 10, 1), (0, 2, 2, 1), (2, 3, 2, 1)])
ledger = ngc.FlowLedger(g)
pol = ngc.FlowPolicy(ledger, 0, 3, config(ngc.ShortestPaths()))
print(pol.place(float("inf")).call_delivered)
try:
    pol.place(1)          # a finite offer after +inf in the same epoch
except ValueError:
    print("ValueError; unchanged:", pol.delivered(), ledger.carried().tolist())
```

```text
5
ValueError; unchanged: 5 [3, 2, 3, 2]
```

## 9. The shortest-path cache

**In plain words.** Many placements need the same static shortest-path DAG toward the same
destinations. A `ShortestPathCache` remembers those DAGs for one graph. It stores the full DAG and
applies each caller's route choice afterwards, so one entry serves every route choice. A cache is either
building (one writer adds entries) or frozen (read-only, shared by many threads, and a miss is an
error).

**Exact rules.**

- **Binding and key.** `ShortestPathCache(graph)` is bound to one graph (`graph_uid()`). The key is
  `(the root set as a sorted set, the direction, the Eligibility uid)`. The value is the full
  static shortest-path DAG: the DAG of every equal-cost walk (`multipath=True`,
  `parallel_edges=True`) with no residual, edges of capacity 0 included. Values are
  immutable and stay valid after `clear()` and after the cache is destroyed.
- **The route choice is never cached.** The single path and the member rule are applied when a
  program is prepared, on the cached DAG, reading the scores of `HIGHER_AVAILABLE_CAPACITY` and
  `HIGHER_HEADROOM` from the residual at that moment and `HIGHER_CAPACITY` from static capacities.
  Because the cached DAG keeps every equal-cost parallel member, a member that a static choice
  would have dropped is still available when its residual makes it the right one.
- **Who reads it.** `ShortestPaths` policy construction reads the cache (keyed by its destination
  set, reverse direction, and Eligibility uid) and retains the support DAG. Later scalar calls and
  lists use that retained support without a cache lease. The evaluators read their supplied cache when they run a `ShortestPaths`
  configuration. The other models do not read it and reject a supplied cache (`ValueError`, for
  example `CentralizedTE does not use ShortestPathCache`), as the last check of the model matrix.

| Operation | Building | Frozen |
| --- | --- | --- |
| lookup hit | returns the value | returns the value; safe from many threads |
| lookup miss | searches, inserts, returns | raises `ShortestPathCacheMiss` |
| `freeze()` | becomes frozen | no effect |
| `clear()` | empties the cache; hit and miss counters persist | `ValueError` |
| `thaw()` | no effect | waits for active calls, then becomes building |
| destruction | callers keep it alive during calls | callers keep it alive during calls; returned values stay valid |

- **Leases.** Each building call owns one exclusive lease for its whole computation; frozen
  calls own shared leases. `freeze()`, `thaw()` and `clear()` take exclusive admission and wait
  for active calls. Admission is FIFO: a queued lifecycle call cannot be overtaken by later
  builders or readers. Evaluator workers borrow their outer call's lease and never acquire it
  again. The cache lease precedes any ledger gate.
- **Statistics.** `stats()` returns `hits`, `misses`, `entries`, `bytes` (the entries' retained
  arrays, without container overhead) and `live_values` (values alive anywhere, cache entries
  included). Each field is atomic; together they are not one consistent snapshot. An invalid
  query (foreign graph, bad root) counts nothing. A frozen cache's statistics follow the public
  call: its hits count only when the call succeeds, its Python result included; a failed call
  changes no statistic, except one miss for a call that raises `ShortestPathCacheMiss`.
- **Workers.** A supplied building cache is accepted only with `workers == 1` (`ValueError`
  otherwise). A frozen cache may be shared by every worker. Without a supplied cache, each
  worker uses its own private building cache.
- **Uids and frozen caches.** A per-pair Eligibility that the evaluator must intersect with a
  common one gets a new uid when the result differs from both, and that key cannot be in a
  frozen cache. To reuse prefilled entries, pass already-combined per-pair values with
  `eligibility=None`, or a common value that contains them (so the intersection returns the
  prefilled operand).
- **Scope.** Use one cache per series of analyses on one graph, and drop it at the end of the
  series.

```python
g = graph(5, [(0, 1, 10, 1), (1, 3, 10, 1), (0, 2, 2, 1), (2, 3, 2, 1), (3, 4, 9, 1)])
cache = ngc.ShortestPathCache(g)
opts = ngc.EvaluationOptions(config(ngc.ShortestPaths()), cache=cache)
rs = ngc.evaluate_pairs(g, [(0, 3), (1, 3)], opts)          # building: workers defaults to 1
print([r.call_delivered for r in rs], cache.stats().hits, cache.stats().misses)
cache.freeze()
rs = ngc.evaluate_pairs(g, [(0, 3), (1, 3)], opts, workers=2)  # frozen: shared by 2 workers
print([r.call_delivered for r in rs], cache.stats().hits, cache.stats().misses)
try:
    ngc.evaluate_demand(g, 0, 4, opts)                # destination 4: not in the cache
except ngc.ShortestPathCacheMiss:
    print("ShortestPathCacheMiss; misses =", cache.stats().misses)
```

```text
[5, 10] 1 1
[5, 10] 3 1
ShortestPathCacheMiss; misses = 2
```

## 10. Evaluators

**In plain words.** The evaluators answer questions about a model without touching your ledgers.
`evaluate_demand` asks how much one model carries between two node sets on a fresh, empty ledger,
either as much as possible or at a given offer. `evaluate_pairs` asks that for many independent
pairs. `evaluate_edge_failures` asks, for each link, how much the answer drops when that link is
removed, at the same offer. A fresh ledger means these answers can differ from placing into a
ledger that other demands already occupy: to ask about a loaded network, use `FlowPolicy.evaluate`
on a ledger with that load.

**Exact rules.**

- **`EvaluationOptions(config, *, eligibility=None, cache=None, outputs=None, design_demand=None,
  offer=None)`** is the one options type of the three evaluators: the configuration, the static
  Eligibility (`None`: everything allowed), the `ShortestPathCache` (`None`: a private cache per
  worker), the optional outputs of each result (`evaluate_edge_failures`: of its baseline), the
  tunnel models' design (`design_demand`) and the `offer`. The design and the offer are kept as
  given (an `int` or a float; another type, a `bool` or a string, raises `TypeError` at
  construction); the evaluator checks them.
- **A finite offer for every model.** `offer`, when given, is a finite integer in `[0, 2^62)` for
  every model under each admission it supports, in all three evaluators: the offer of every pair,
  of the baseline and of every failed network. Each evaluation computes the same amounts as a
  fresh policy's `evaluate(offer, design_demand=...)`: Lossless delivers up to the offer, Lossy
  accounts for exactly that offer. An absent offer is the uncapped Lossless question or, only for
  `ShortestPaths` `EQUAL`/`WEIGHTED` with Lossy, the first-hop experiment.

  | Model | Lossless | Lossy |
  | --- | --- | --- |
  | `ShortestPaths` `EQUAL` / `WEIGHTED` | `offer` `None`: the uncapped bound; finite: delivery up to it | `offer` finite: that offer; `None`: the first-hop experiment |
  | `ShortestPaths` `OPTIMIZED`, `CentralizedTE` without a tunnel count | `offer` `None`: the uncapped optimum; finite: the optimum up to it | refused by the matrix |
  | T1 | `design_demand` required and finite; `offer` `None`: the lossless bound of those tunnels; finite: up to it | `design_demand` required and finite; `offer` required |
  | T2 | `design_demand` `None`: the offer if it is finite, else unlimited; `offer` `None`: the bound at that design | `offer` finite; `None`: the design if finite, else an error |
  | `Pinned` | `offer` `None`: the fixed-program bound; finite: up to it | `offer` required |

- **Validation order** (before any worker starts): the policy's construction checks (the
  configuration, [section 7](#configuration)); then an explicitly given offer (`offer must be an
  integer in [0, 2^62)`, an explicit `float('inf')` included); then the design: only T1 and T2
  take one (`<model> does not accept design_demand`, checked before its range), in `[0, 2^62)` or
  `float('inf')` (`design_demand must be an integer in [0, 2^62) or permitted +inf`), T1's finite
  (`DistributedTE requires a finite design demand`); then T1's missing design (`DistributedTE
  evaluation requires design_demand`); then the missing Lossy offer (`DistributedTE LOSSY
  evaluation requires a finite offer`, `Pinned LOSSY evaluation requires a finite offer`, and
  `CentralizedTE LOSSY evaluation requires a finite offer` for T2 when no finite design gives one).
  The model in these texts is the class name, T2 included.
- **`evaluate_demand(graph, sources, destinations, options)`** evaluates a `FlowPolicy` with
  `options.config` on a fresh ledger of the graph and returns its `PlacementResult`; it commits
  nothing anywhere. The exact maximum flow is the configuration
  `CentralizedTE(source_allocation=OPTIMIZED)`: its result carries the marginal cost curve and,
  with `OutputOptions(min_cut=True)`, the minimum cut. With `EQUAL` source allocation it is the
  largest total whose equal shares all fit.
- **`evaluate_pairs(graph, pairs, options, per_pair_eligibility=None, workers=None)`** runs
  `evaluate_demand` for each `(sources, destinations)` pair and returns the results in pair
  order. Nothing is shared between pairs (no reservations, no capacity): the results are
  independent capacities, not allocations that coexist; demands that share capacity are policies
  on one ledger. Pair `i` uses `options.eligibility` intersected with `per_pair_eligibility[i]`
  (entries may be `None`; a given list has one entry per pair). All checks finish before evaluation workers start: first the arguments (the worker count, the graph of `options.eligibility` and
  of the cache, `Pinned`, which is rejected because a program fixes its endpoints (`evaluate_pairs does not accept
  Pinned`), the configuration, the offer and design rules above, then the length of `per_pair_eligibility`), then each
  pair completely (its endpoint sets, their overlap, its Eligibility's graph, the pooled single-path rule of
  [section 7](#configuration), T2's several-sources rule and, for the tunnel models, the tunnel limit for its sources)
  in pair order, so the error of the lowest failing pair is the one raised. A supplied building cache's
  `workers == 1` rule reads the cache's mode and comes last, under the cache's lease. An empty batch returns `[]`.
- **`evaluate_edge_failures(graph, sources, destinations, options, workers=None)`** returns a
  `SensitivityResult`: `baseline` (the result of `evaluate_demand`, with the requested outputs)
  and `entries`, one `SensitivityEntry(edge, delivery_reduction)` per candidate edge, in ascending
  EdgeId.
  - **Candidates.** With unrestricted `CentralizedTE` and `OPTIMIZED` source allocation: the edges
    that carry flow in the baseline (removing an edge without flow keeps the baseline feasible, so
    the delivery stays, at any offer). With every other configuration, `EQUAL` source allocation
    included: every allowed edge,
    whatever its capacity (an edge of capacity 0 can carry a route and keep a source in the
    roster, so deleting it can change the result either way). An edge without an entry has
    reduction 0.
  - **`delivery_reduction`** is `baseline.delivered` minus the delivery of `evaluate_demand` with
    that edge also excluded, at the same offer and design, in units, signed and never rounded or
    thresholded. It is not the traffic the deleted edge carried. A negative value means the model
    carries more without the edge (for example, equal splitting no longer sends a share down a
    thin path, or a source behind it no longer holds the equal allocation back). Unrestricted
    `CentralizedTE` with `OPTIMIZED` source allocation, when its baseline is a maximum flow (no
    offer, or an offer above what it delivers), obtains the same delivery reduction as a fresh solve through
    exact residual certificates or the first repair from a private baseline copy. A positive-flow
    saturated edge across residual strongly connected components loses exactly its flow; other
    edges lose their flow minus the amount the first repair can reroute. The intermediate amounts
    may violate conservation and are never published or used as a new baseline. Each actual repair
    resets independently. These shortcuts preserve the complete baseline and the candidate entries/order.
    At an offer the baseline delivers (the offer binds, perhaps below the maximum), and with
    `EQUAL` source allocation or a support, each deletion is solved afresh. T1 sets up again for the same design
    under each deletion; T2 redoes its sizing, count choice and objectives, on the same
    supplied support when there is one. `Pinned` replays the same program under each deletion (a
    repaired fixed-route question, not a new optimum and not stale forwarding): deleting an edge of
    a tunnel fails that tunnel. `ShortestPaths` gives reduction 0 without recomputation to an edge
    outside the static shortest-path support reachable from the roster sources (no route can use
    it); these are equivalent forms of the same definition.
  - **Cache.** Each deletion uses a new Eligibility uid, so a `ShortestPaths` deletion that has to
    be placed raises `ShortestPathCacheMiss` through a frozen cache. The other models reject a
    supplied cache. Private caches are cleared between candidates.
- **Workers.** An explicit `workers` (an integer of at least 1) is always honoured; the
  evaluator uses at most one thread per pair or candidate. With `workers=None` the default is 1
  when a cache is supplied, and otherwise comes from the environment:

  | Variable | Default worker count of | Accepted values |
  | --- | --- | --- |
  | `NGRAPH_CORE_BATCH_THREADS` | `evaluate_pairs` | an integer in `[1, 4096]`; unset or anything else means 1 |
  | `NGRAPH_CORE_SENSITIVITY_THREADS` | `evaluate_edge_failures` | the same |

  The variables only supply defaults; they never override an explicit `workers`. The C++
  functions take the worker count as an explicit argument.
- **Failures** follow [section 8](#8-failures-threads-and-lifetimes): a failed evaluator call
  raises (also a solve beyond the flow core's size limit of [section 2](#2-the-graph), with
  `CentralizedTE` or `OPTIMIZED`), and a failed deletion fails the whole call (the exception of the
  lowest failing candidate; no entry is returned); a supplied building cache keeps the entries of
  the lookups that succeeded before the failure; a supplied frozen cache's statistics are
  unchanged, except one miss when the call raises `ShortestPathCacheMiss`.

```python
g = graph(4, [(0, 1, 10, 1), (1, 3, 10, 1), (0, 2, 2, 1), (2, 3, 2, 1)])
for name, cfg in (("ShortestPaths", config(ngc.ShortestPaths())),
                  ("CentralizedTE", ngc.FlowPolicyConfig(routing=ngc.CentralizedTE(source_allocation=ngc.SourceAllocation.OPTIMIZED)))):
    s = ngc.evaluate_edge_failures(g, 0, 3, ngc.EvaluationOptions(cfg))
    print(name, s.baseline.call_delivered, [(e.edge, e.delivery_reduction) for e in s.entries])
```

```text
ShortestPaths 5 [(0, 3), (1, -5), (2, 3), (3, -5)]
CentralizedTE 12 [(0, 10), (1, 2), (2, 10), (3, 2)]
```

Under equal splitting the thin path (edges 1 and 3) holds the whole diamond back to 5 units;
without it the wide path alone carries 10, a reduction of -5.

A failure study at a fixed offer, against the first-hop experiment: two parallel links of 10, ECMP,
Lossy. With an offer of 4 every network delivers 4; without an offer the experiment injects each
remaining link's capacity, 20 before a deletion and 10 after, a different question:

```python
g = graph(2, [(0, 1, 10, 1), (0, 1, 10, 1)])
lossy = config(ngc.ShortestPaths(), admission=ngc.Admission.LOSSY)
for offer in (4, None):
    s = ngc.evaluate_edge_failures(g, 0, 1, ngc.EvaluationOptions(lossy, offer=offer))
    print(offer, s.baseline.call_delivered, [e.delivery_reduction for e in s.entries])
try:
    ngc.evaluate_demand(g, 0, 1, ngc.EvaluationOptions(lossy, offer=float("inf")))
except ValueError as error:
    print(error)
```

```text
4 4 [0, 0]
None 20 [10, 10]
offer must be an integer in [0, 2^62)
```

### Limits

- **One route allocation per program.** A program applies one route allocation at every router.
  Weighted ingress above routers that should split equally is expressible exactly when every
  member of those downstream equal splits has positive static capacity: give each member of such a
  split the same positive `EXPLICIT` weight, since `H(x, [w] * k) = RR(x, k)` for any positive `w`.
  A downstream equal split with a zero-capacity member differs (`WEIGHTED` prunes that member,
  `EQUAL` keeps it) and is not expressible.
- **Certificates and exports.** The public optimality certificates are `min_cut` and
  `marginal_cost_curve` of unrestricted optimized `CentralizedTE`. A program has no serialized
  form for moving it between processes, and its per-router tables are not exported.
- **Not modelled:** traffic classes and priority queues, RSVP setup and holding priorities and
  preemption, bandwidth policing, automatically disjoint or SRLG-diverse tunnels, a stretch
  baseline per head inside one pooled T1 demand, growth of a tunnel count over time, and a joint
  optimum or joint bound over several demands. Per-pair offers, designs or configurations are
  separate `evaluate_demand` calls, and a study of the three failure states is a loop of the
  caller's.

## 11. C++ names

The C++ library (namespace `netgraph::core`) has the same semantics. Its public headers are
self-contained:

| Header | Declares |
| --- | --- |
| `netgraph/core/types.hpp` | `NodeId`, `EdgeId` (`int32`), the cost and amount aliases (`int64`), `kDomainLimit` (`2^62`), `kUnreachable` |
| `netgraph/core/graph.hpp` | `Graph`, `to_units`, `as_float` |
| `netgraph/core/eligibility.hpp` | `Eligibility` (`allows_node`, `allows_edge`) |
| `netgraph/core/route_structures.hpp` | `PredecessorDAG`, `Path` |
| `netgraph/core/options.hpp` | `Direction`, `AbsoluteCost`, `CostStretch`, `RouteTieBreak`, `RouteAllocation`, `WeightBasis`, `Admission`, `SourceAllocation` |
| `netgraph/core/forwarding_program.hpp` | `ForwardingProgram` (`from_routes`, `from_dag`), `Tunnel`, `TunnelStatus`, `RouteDefinition`, `RouteOptions` (the options of `from_routes` and `from_dag`) |
| `netgraph/core/search.hpp` | `spf`, `ksp`, `resolve_to_paths`, `SpfOptions`, `SpfResult`, `KspOptions`, `SpfQueue` |
| `netgraph/core/flow_record.hpp` | `FlowId`, `Counter128` (the 128-bit cumulative counter), `FlowRecord` |
| `netgraph/core/flow_ledger.hpp` | `FlowLedger` |
| `netgraph/core/flow_policy.hpp` | `FlowPolicy`, its configuration (`FlowPolicyConfig`, `ShortestPaths`, `DistributedTE`, `CentralizedTE`, `Pinned`, `TunnelCount`, `OutputOptions`) and result types (`PlacementResult`, `SourceFlow`, `TunnelFlow`, `TunnelPlan`), `Termination`; the lists: `apply_actions`, `BatchMode`, `Place`, `Replace`, `Forward`, `Action` |
| `netgraph/core/shortest_path_cache.hpp` | `ShortestPathCache`, `ShortestPathCacheMiss` |
| `netgraph/core/evaluators.hpp` | `evaluate_demand`, `evaluate_pairs`, `evaluate_edge_failures`, `EvaluationOptions`, `SensitivityEntry`, `SensitivityResult` |

Spelling differences: Python enum members and constants are upper case
(`RouteTieBreak.FEWER_HOPS`, `UNREACHABLE`), C++ enumerators are the same words in PascalCase
(`RouteTieBreak::FewerHops`, `RouteAllocation::Optimized`) and C++ constants carry a `k` prefix
(`kUnreachable`). Where Python returns `None` for an unreachable cost (`target_cost`,
`min_route_cost`), C++ returns `kUnreachable`; an infinite offer is an empty `std::optional`, and
its `call_undelivered` and `requested` are empty too. An optional field that Python reports as
`None` (`ShortestPaths`' `source_allocation` and `weight_basis` before resolution, `offered`,
`unadmitted`, `total_cost`, `marginal_cost_curve`, `plan` and T2's fields of a `TunnelPlan`,
`CentralizedTE`'s `tunnel_count`, `FlowRecord.entry_id`, a down tunnel's `destination`, `edges`
and `path_cost`, the optional outputs) is an empty `std::optional` in C++;
`max_fractional_bandwidth_per_tunnel` is a reduced `(numerator, denominator)` pair,
`total_tunnel_cost` and `headroom_penalty` are `Counter128`; `CostStretch` is constructed from a
numerator and a denominator (reduced on construction, `numerator()` and `denominator()`); a
`design_demand` or an evaluator's `offer` is a `std::optional<std::optional<Amount>>`, absent for
`None` and an empty inner value for `float('inf')`; and a `None` program, support or black holes
a null `std::shared_ptr` (`edge_weights` and `source_weights` are shared arrays, `std::shared_ptr<const
std::vector<Amount>>`). `CentralizedTE` has no default constructor: it is built from its
`SourceAllocation`. The epoch's amounts are `Counter128` values. A Python `ValueError` is a
`std::invalid_argument` in C++. `spf` takes a trailing C++-only `SpfQueue` argument that chooses
the search's internal queue; the results never depend on it. `apply_actions` takes a
`std::span<const Action>`, where `Action` is a `std::variant<Place, Replace, Forward>` whose
values borrow their policy (a `std::reference_wrapper<FlowPolicy>`) and black holes (a `const
Eligibility*`, null for none), which must outlive the call.
