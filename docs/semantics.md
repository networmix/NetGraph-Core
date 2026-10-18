# NetGraph-Core semantics

This document says exactly what NetGraph-Core computes: how amounts are counted, how edges are
numbered, which route a search returns, how a placement spreads traffic, what a placement call
reports and when it fails, how the shortest-path cache behaves, and what the evaluators mean. It
is the full reference behind the API docstrings. For installation and a first example, see the
[README](../README.md#installation).

Each section starts with a short explanation in plain words and then gives the exact rules. When
the two seem to differ, the exact rules are the authority. Names are the Python spellings
(`import netgraph_core as ngc`); [the last section](#11-c-names) lists the C++ spellings. Names that
start with an underscore are private test and verification entries, not part of the API.

Contents:

1. [Units: the one numeric domain](#1-units-the-one-numeric-domain)
2. [The graph](#2-the-graph)
3. [Route values: PredDAG and Path](#3-route-values-preddag-and-path)
4. [Search](#4-search)
5. [Placement kernels and admission](#5-placement-kernels-and-admission)
6. [The ledger](#6-the-ledger)
7. [Policies](#7-policies)
8. [Failures, threads and lifetimes](#8-failures-threads-and-lifetimes)
9. [The DAG cache](#9-the-dag-cache)
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
  destinations (a `FlowPolicy`): it creates, fills, grows and removes that demand's flows.
- **Cohort.** The flows a policy manages for one source; growth adds a flow to it
  (`CapacityAware`).
- **Roster.** The sources of a policy that take part in placement, fixed when the policy is built.
- **Origin.** How a demand's total enters at its roster sources: in equal round-robin shares
  (`EQUAL`) or as a controller chooses (`OPTIMIZED`).
- **Round robin `RR(x, k, i)`.** `x` units over `k` members: `x // k` each, and one more for the
  first `x % k` members.
- **Highest averages `H(x, w)`.** `x` units over members with integer weights `w`, given one at a
  time to the member with the largest `w[i] / (assigned[i] + 1)`, ties to the earlier member; a
  zero weight gets nothing.
- **Program.** The routes and splits a `CostOnly` policy prepares once per epoch and forwards
  over.
- **Controller.** A routing that solves a flow problem for the whole demand: `Global`, and
  `CostOnly` with `MAX_FLOW`.
- **Offer.** The part of a request that a call asks one source to place.
- **Epoch.** The period from a policy's construction, or a `remove()`, to the next `remove()`;
  requests add up within it.
- **Lossless / Lossy.** Admission control: Lossless places only what fits end to end; Lossy
  forwards hop by hop and drops what does not fit on each edge.
- **Static.** Depending only on the graph, the costs and the Eligibility, not on traffic already
  placed (a static walk, the static shortest cost).
- **Member graph.** The edges a search may use: eligible (and, when a residual is supplied, at
  least the search's residual threshold), whatever their capacity, and one member per pair with
  `multi_edge=False`.
- **Kernel.** The rule that spreads one flow over its route: `MAX_FLOW`, `WEIGHTED` or
  `EQUAL_BALANCED`.
- **Bundle.** The route of one flow, as a `PredDAG`.
- **View.** The capacity a solve may use: the free capacity plus the flow's (or the `Global`
  policy's) own allocation.
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
    """edges: (src, dst, capacity_units, cost) tuples in input order."""
    src, dst, cap, cost = zip(*edges)
    return ngc.Graph(num_nodes, list(src), list(dst), list(cap), list(cost), q)

DET = ngc.PathSelection(ngc.EdgeTieBreak.DETERMINISTIC, ngc.EdgePriority.CAPACITY)
PREFER = ngc.PathSelection(ngc.EdgeTieBreak.PREFER_HIGHER_PRIORITY, ngc.EdgePriority.CAPACITY)

def config(routing, placement=ngc.FlowPlacement.EQUAL_BALANCED,
           admission=ngc.Admission.LOSSLESS, selection=DET):
    """A configuration with a kernel (every routing but Global, which takes none)."""
    return ngc.FlowPolicyConfig(routing=routing, placement=placement,
                                admission=admission, selection=selection)
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
  `offered` and its per-edge `drops`) may exceed 64 bits; Python returns them as exact `int`.
- **Domain.** Every edge capacity, and the sum of all capacities, is below `2^62` units. Edge
  costs are integers of at least 1 and their sum is below `2^62`, so no path cost reaches the
  sentinel `UNREACHABLE` (the largest `int64`). A request is an integer in `[0, 2^62)` units or
  positive infinity (`float('inf')`, "as much as possible"). Everything outside the domain is a
  `ValueError` raised before anything changes.
- **`to_units(x, q)`** converts physical values to units, rounding half up: a scalar gives an
  `int`, an array an `int64` array. NaN, infinity, negative values and results at or above
  `2^62` units raise `ValueError`. It is the only way physical numbers enter the library.
- **`as_float(x, q)`** converts units to physical values (`x * q`, correctly rounded, exact
  below `2^53` units): an `int` gives a `float`, an integer array a `float64` array. A result
  beyond the finite doubles raises `ValueError`. Costs are never converted.
- **Routability.** An edge with less than 1 unit of capacity carries nothing, but it stays in
  routing: searches, routes and rosters treat it like any other edge, and only a supplied
  residual (below the search's threshold) or an Eligibility removes an edge from a search. An
  equal split still sends such an edge its share, which it cannot carry
  ([section 5](#5-placement-kernels-and-admission)). A physical request below `q / 2` converts to
  0 units. A request of 0 units is accepted after
  validation and does nothing (its result is `COMPLETE`).
- **Splitting units.** When a node forwards `I` units over `k` member edges:
  - **round robin** (`EQUAL_BALANCED`): every member gets `I // k` units and the first `I % k`
    members, in ascending EdgeId, get one more;
  - **highest averages** (`WEIGHTED`, the D'Hondt rule): the units are handed out one at a time,
    each to the member with the largest `w / (a + 1)` (its weight `w` over its current count
    `a` plus one), ties to the lower EdgeId. Equal weights give round robin.

  Both rules are exact and never lose a unit: a share that rounds down to zero leaves its unit to
  a sibling. Both are monotone: when `I` grows, no member's share shrinks.
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
tie rule of the library ("the lowest EdgeId") refers to this order.

**Exact rules.**

- **Construction.** `Graph(num_nodes, src, dst, capacity_units, cost, resolution)` builds an
  immutable directed multigraph: parallel edges and self-loops are allowed. The arrays are in
  input order (the caller's). Capacities are already in units; use `to_units` for physical
  values.
- **Validation order.** When an input has several defects, the first in this order is reported as a
  `ValueError`: the resolution; the capacity values (each an integer in `[0, 2^62)`; floats and
  booleans are rejected, not rounded); the array shapes; a zero or negative cost; the remaining
  bounds (node ids in range, then the capacity and cost sums below `2^62`, both added up edge by
  edge in input order: the first edge at which a running sum reaches `2^62` reports that sum, the
  capacity sum before the cost sum at the same edge). Integers are read
  without narrowing, so the reported error does not depend on how the same input is represented.
- **Size.** Node ids and EdgeIds are `int32`: a graph has at most `INT32_MAX` (`2^31 - 1`) nodes
  and at most as many edges. The constructor accepts every graph within these bounds, and the
  searches and route values work on all of them. The flow core (`Global` routing, the `MAX_FLOW` kernel and the
  evaluators that use them) numbers a solve's network in `int32` as well: it adds a virtual
  source and a virtual sink, and gives every edge and every virtual arc (one per source, one per
  destination) two slots, one to add flow and one to cancel it. A solve therefore raises
  `ValueError`, before the solver sizes anything and with nothing changed, on a network of more
  than `INT32_MAX - 2` nodes, or with more than `2^30 - 1` edges, sources and destinations
  together. A `Global` solve counts the graph's nodes, the eligible edges, the roster sources and
  the eligible destinations; a `MAX_FLOW` solve counts the nodes of the flow's member edges and
  endpoints, its member edges, its source and its destinations. Only graphs of billions of nodes
  or about a billion edges reach these limits.
- **EdgeId.** The EdgeId of an edge is its rank in the stable order `(src, cost, dst, input
  position)`. `input_index()[e]` is the input position of edge `e`; the inverse map is the
  caller's to build (one pass over `E`).
- **EdgeId order is the coordinate of every per-edge value.** Residual arrays, Eligibility
  exclusions, DAG and Path edges, flows, records, `bottleneck`, `min_cut` and sensitivity entries
  all use EdgeIds. A caller that keeps per-edge data in input order translates it first: an
  in-range id in the wrong coordinate is silently a different edge, and no check can catch that.
- **Layout guarantees.** The out-edges of node `u` are exactly the ids
  `row_offsets()[u]` to `row_offsets()[u + 1] - 1`, in `(cost, dst, input)` order, and
  `edge_dst()[e]` is the head of `e`. The in-edges of node `v` are
  `in_adj_edge_index()[in_row_offsets()[v]:in_row_offsets()[v + 1]]`, as EdgeIds in
  `(cost, src, input)` order. Both offset arrays have `N + 1` entries, starting at 0 and ending
  at `E`, also for an empty graph.
- **Accessors.** `capacity_units()`, `cost()` (`int64[E]`), `edge_src()`, `edge_dst()`,
  `input_index()`, `row_offsets()`, `in_row_offsets()` and `in_adj_edge_index()` (`int32`) are
  read-only NumPy views that keep the graph alive; `max_cost()` and `cost_gcd()` are 0 for a
  graph without edges; `uid()` is a process-unique number, never reused.
- **Ids belong to one graph.** Building a graph from different arrays, or with another release of
  the library, may renumber the edges. Every value bound to a graph (Eligibility, PredDAG,
  Path, ledger, cache) is rejected with another graph.

```python
# Input order: edge 0 = 1->2, edge 1 = 0->2 (cost 5), edge 2 = 0->1.
g = ngc.Graph(3, [1, 0, 0], [2, 2, 1], [7, 8, 9], [1, 5, 1], 1.0)
print("src", g.edge_src().tolist(), "dst", g.edge_dst().tolist(), "cost", g.cost().tolist())
print("input_index", g.input_index().tolist(), "row_offsets", g.row_offsets().tolist())
```

```text
src [0, 0, 1] dst [1, 2, 2] cost [1, 5, 1]
input_index [2, 1, 0] row_offsets [0, 2, 3, 3]
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
- Uids matter for the [DAG cache](#9-the-dag-cache): its key contains the Eligibility uid, so
  build each state once and reuse the value.
- An excluded endpoint takes no part in a search or a policy: an excluded root or target in a
  set is skipped, and an excluded single endpoint gives an empty result.

## 3. Route values: PredDAG and Path

**In plain words.** A route is either a DAG of predecessor edges (the shape a search returns:
possibly many equal-cost walks) or a `Path` (one explicit walk, each hop possibly a bundle of
parallel edges). Both are immutable, bound to their graph and checked once when built.

**Exact rules.**

- **`PredDAG`.** For node `v`, its entries are `parents()[k]` (the parent node) and
  `via_edges()[k]` (the edge from that parent) for `k` from `parent_offsets()[v]` to
  `parent_offsets()[v + 1] - 1`, in strictly ascending EdgeId. This is the one canonical form:
  two DAGs of the same graph are equal exactly when their entries are equal; DAGs of different
  graphs are never equal. `parents()` and `via_edges()` are read-only views; `parent_offsets()`
  is built on demand.
- **`PredDAG.from_csr(graph, offsets, parents, via)`** validates, before reading any structure:
  `offsets` has `N + 1` entries, starts at 0, never decreases and ends at
  `len(parents) == len(via)`; then every id is in range, every via edge joins its parent to its
  node, no `(node, via)` pair repeats, the DAG has no cycle and every node's entries are in
  ascending via EdgeId. Any failure is a `ValueError`.
- **`Path(graph, hops)`.** `hops` is a list of hops, each a non-empty list of EdgeIds that share
  one `(u, v)` pair (parallel members). Validation, in order: every member id in range, no
  member twice, at least one hop, no empty hop, the members of a hop share `(u, v)`, consecutive
  hops connect, no node is visited twice. The members keep the caller's order. `cost()` is the
  sum over the hops of each hop's cheapest member; `src()`, `dst()` and `hops()` read the walk;
  `to_dag()` gives the equivalent `PredDAG`. A single-member walk is `[[e1], [e2], ...]`.
- **`resolve_to_paths(graph, dag, src, dst, split_parallel_edges=False, max_paths=None)`** lists
  the `src -> dst` walks of a DAG as `Path` values in increasing `(cost, hop count, edge-id
  sequence)`, where a hop costs its cheapest member and is represented in the key by its lowest
  member id. With `split_parallel_edges=False` the parallel members of a pair stay together in
  one hop (ascending); with `True` every member combination is a separate path, and the paths are
  sorted by each path's own `(cost, hop count, edge-id sequence)`. `max_paths` (at least 1) keeps
  the first paths of that order. `src == dst`, a DAG of
  another graph or `max_paths < 1` raise `ValueError`; an endpoint outside the DAG, or no walk,
  gives `[]`.

## 4. Search

**In plain words.** `spf` finds shortest paths from a set of roots, forward or in reverse,
either to every node or only to the nearest of a set of targets. It returns the shortest-path
DAG, which contains every shortest route at once, or, on request, one canonical path chosen by a
fixed rule. `ksp` lists the K shortest simple paths in a fixed order. All searches see every
eligible edge, whatever its capacity: an edge of capacity 0 is a route like any other. Only a
supplied residual filters edges by capacity.

### Options shared by the searches

- **`PathSelection(tie_break, priority, multipath=True, multi_edge=True)`.** `tie_break` and
  `priority` have no default.
  - `multipath=True` returns the DAG of every equal-cost walk; `False` the one canonical path.
  - `multi_edge=True` keeps every eligible parallel member of a `(u, v)` pair; `False` keeps
    one per pair (the member rule below).
  - `tie_break`: `DETERMINISTIC` never compares priority values; ties fall to the lowest
    EdgeId. `PREFER_HIGHER_PRIORITY` prefers the higher priority value first, then the lowest
    EdgeId.
  - `priority`: the per-edge value that `PREFER_HIGHER_PRIORITY` reads (and that the Weighted
    kernel uses for membership and weights): `RESIDUAL` (residual units), `CAPACITY` (static
    capacity units) or `HEADROOM` (residual divided by capacity, compared exactly as a
    fraction; 0 on an edge of capacity 0). A search without a residual array reads the static
    capacity as the residual: `RESIDUAL` and `CAPACITY` read the capacity, `HEADROOM` reads 1 on
    an edge of positive capacity and 0 on an edge of capacity 0.
- **The member rule** (`multi_edge=False`). Of the eligible parallel members of a pair, keep the
  one with the lowest cost; among those, under `PREFER_HIGHER_PRIORITY`, the one with the
  highest priority value; then the lowest EdgeId. This choice happens before any route is
  searched, so it can change which walks exist.
- **Static eligibility** (always): the Eligibility. Capacity is never a condition.
- **Residual gating** (only when a `residual` array is supplied): `residual` is an `int64` array
  in EdgeId order with every entry in `[0, capacity_units[e]]`. It excludes edges whose residual
  is below the threshold (`min_residual` units for `spf`, 1 unit for `ksp`) and supplies the
  residual values of the priority. It is read during the call only and must not change during
  it.

### Shortest paths

`spf(graph, roots, SpfOptions(direction, selection, targets=None, residual=None,
min_residual=1, elig=None))` returns an `SpfResult` with `costs`, `dag` and `target_cost`.

- **Roots and targets** are one node id or a sequence. Every root starts at distance 0.
  Duplicates are removed. Validated in this order before any work (`ValueError`): the roots
  (non-empty, ids in range), the targets (the same), a node in both sets, `multipath=False`
  without targets, the residual (length, then every entry), `min_residual` (an integer of at
  least 1), a foreign Eligibility. Excluded roots and targets are skipped; when every root is
  excluded the result is empty, not an error.
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
  changes through `multi_edge`, the member rule, static eligibility and residual gating; with the
  member graph fixed, the tie mode, the priority, the capacity and residual values and the C++
  queue choice do not change the DAG. A forward search with targets and a reverse search with the
  roots and the targets exchanged return the same DAG. Without targets the two directions can
  return different DAGs, also when the reverse search is rooted at a node the forward search
  reaches: each keeps the edges tight in its own direction.

### The canonical single path

With `multipath=False` (targets required), `spf` returns the one walk from a root to a target
that is first in this order:

1. the lowest **cost**;
2. under `PREFER_HIGHER_PRIORITY` only: the highest **bottleneck**, the smallest priority value
   along the walk (`DETERMINISTIC` skips this step);
3. the fewest **hops**;
4. the lexicographically smallest **sequence of EdgeIds** along the walk.

The order is total, so the walk is unique. The same walk is returned by a forward search, by a
reverse search with the roles exchanged, and by `ksp` with `k=1`. Note the order of steps 3 and
4: a walk with fewer hops wins even when its first EdgeId is larger. The member rule above is a
separate, earlier choice: with `multi_edge=False` each pair is first reduced to one member, and
the canonical path is then chosen among the walks that remain.

```python
F = ngc.Direction.FORWARD
# Two equal-cost walks 0->1->3 (capacity 1) and 0->2->3 (capacity 5).
g = graph(4, [(0, 1, 1, 1), (0, 2, 5, 1), (1, 3, 1, 1), (2, 3, 5, 1)])
for sel in (DET, PREFER):
    single = ngc.PathSelection(sel.tie_break, sel.priority, multipath=False)
    r = ngc.spf(g, 0, ngc.SpfOptions(F, single, targets=3))
    print(sel.tie_break.name, r.target_cost, r.dag.via_edges().tolist())

# Equal cost 2: the direct edge (EdgeId 1) against the two-hop walk [0, 2].
g = graph(3, [(0, 1, 1, 1), (1, 2, 1, 1), (0, 2, 1, 2)])
single = ngc.PathSelection(ngc.EdgeTieBreak.DETERMINISTIC, ngc.EdgePriority.CAPACITY,
                           multipath=False)
r = ngc.spf(g, 0, ngc.SpfOptions(F, single, targets=2))
print("fewest hops:", [p.hops() for p in ngc.resolve_to_paths(g, r.dag, 0, 2)])
full = ngc.spf(g, 0, ngc.SpfOptions(F, DET))
print("costs", full.costs.tolist(), "dag", full.dag.via_edges().tolist())
```

```text
DETERMINISTIC 2 [0, 2]
PREFER_HIGHER_PRIORITY 2 [1, 3]
fewest hops: [[[1]]]
costs [0, 1, 2] dag [0, 1, 2]
```

### K shortest paths

`ksp(graph, src, dst, KspOptions(k, selection, max_path_cost=None, residual=None, elig=None))`
returns a list of `(cost, Path)` pairs; every hop of every Path has one member.

- **Order.** The simple `src -> dst` paths of the member graph in increasing `(cost,
  -bottleneck under PREFER_HIGHER_PRIORITY, hop count, edge-id sequence)`: the canonical order
  above applied to whole paths. The result is the first `k` paths of that order that are within
  the ceiling.
- **Parallel members.** `multi_edge=True` makes every eligible member a separate path;
  `multi_edge=False` keeps one member per pair by the member rule. `multipath` is
  ignored.
- **Ceiling.** `max_path_cost` is `None`, `Absolute(c)` (cost at most `c`, `c >= 0`; `Absolute(0)`
  is a valid budget) or `Stretch(f)` (cost at most `f` times the static shortest cost of the
  member graph under the Eligibility, residual ignored, edges of capacity 0 included; `f` finite
  and at least 1). Both are
  inclusive and compared exactly: `f` is read as the exact binary fraction a float is. With
  `Stretch`, no static route gives `[]`.
- **Residual gating** at 1 unit when a residual is supplied.
- Validated in this order before any work (`ValueError`): `src` and `dst` in range, `src ==
  dst`, `k >= 1`, the ceiling, the residual (length, then every entry), a foreign Eligibility. An
  excluded endpoint or no path gives `[]`.

```python
g = graph(3, [(0, 1, 3, 1), (0, 1, 5, 1), (1, 2, 9, 1)])   # parallel members 0 and 1
collapse = ngc.PathSelection(ngc.EdgeTieBreak.PREFER_HIGHER_PRIORITY,
                             ngc.EdgePriority.CAPACITY, multi_edge=False)
for sel in (DET, PREFER, collapse):
    print([(c, p.hops()) for c, p in ngc.ksp(g, 0, 2, ngc.KspOptions(5, sel))])
r = ngc.spf(g, 0, ngc.SpfOptions(F, DET, targets=2))
print([p.hops() for p in ngc.resolve_to_paths(g, r.dag, 0, 2)],
      [p.hops() for p in ngc.resolve_to_paths(g, r.dag, 0, 2, split_parallel_edges=True)])
```

```text
[(2, [[0], [2]]), (2, [[1], [2]])]
[(2, [[1], [2]]), (2, [[0], [2]])]
[(2, [[1], [2]])]
[[[0, 1], [2]]] [[[0], [2]], [[1], [2]]]
```

## 5. Placement kernels and admission

**In plain words.** A kernel decides how one flow spreads over its route. `MAX_FLOW` behaves like
a controller that pushes as much as the route can carry, cheapest walks first. `WEIGHTED` is a
router splitting traffic at every node by fixed weights (WCMP). `EQUAL_BALANCED` is a router
splitting equally at every node (ECMP), or equal-share tunnels. Admission is a separate choice:
`LOSSLESS` places only what fits end to end; `LOSSY` forwards hop by hop, as routers do, and
drops on each edge what does not fit.

**Exact rules.**

- **Split graph.** A flow uses the part of its route (DAG or bundle) that reaches its destination
  set through the edges its kernel may use (its member edges). Which edges: `MAX_FLOW` uses the
  edges with capacity in its view (the free capacity plus the flow's own allocation);
  `EQUAL_BALANCED` uses every edge of the route, an edge of capacity 0 included: a router splits
  over every next hop it has, so that edge receives its round-robin share and carries none of it
  (Lossless can then admit nothing on that branch; Lossy drops the share on that edge);
  `WEIGHTED` uses the edges of positive weight, by its `priority` the edges with positive residual
  when the flow was created (`RESIDUAL`) or with at least 1 unit of static capacity
  (`CAPACITY`), and drops a branch left without such an edge (never taking a route outside its
  own). `WEIGHTED` with `HEADROOM` is rejected. Lossy keeps its kernel's membership.
- **`MAX_FLOW`** solves the minimum-cost flow of the flow's cumulative amount on its route, over
  the free capacity plus the flow's own current allocation. It may move the flow's own earlier
  allocation (that is what lets it reach the maximum); it never moves another flow's. It fills
  the cheapest walks of a mixed-cost bundle first. Per-edge changes between calls can be
  negative. With `CostOnly` routing it solves the whole demand at once, over the union of its
  sources' routes ([section 7](#routing-modes)). `MAX_FLOW` is Lossless only: a controller never
  drops. A flow beyond the flow core's
  size limit raises `ValueError` ([section 2](#2-the-graph), Size).
- **`WEIGHTED`** splits at every node by highest averages over weights frozen when the flow is
  created: the members' priority values at that moment (residual units under `RESIDUAL`, static
  capacity units under `CAPACITY`). Weights are never renormalised; a member whose residual later
  shrinks limits the flow.
- **`EQUAL_BALANCED`** splits at every node by round robin over its members, which are frozen
  when the flow is created.
- **Lossless filling** (`WEIGHTED`, `EQUAL_BALANCED`). On a fixed route and plan, a flow's delivered
  amount `T` only grows (an equal-share cohort that gains a tunnel re-splits its total; see [the
  placement call](#the-placement-call)): a call raises it to the smallest of `T` plus the call's
  budget and `T_max`, the largest amount whose allocation fits every member edge (found by exact
  binary search). At `T_max` at least one member edge is exactly full; there is no safety margin.
  Because the split rules are monotone, raising `T` never takes units away from an edge.
- **Lossy forwarding** (`WEIGHTED` or `EQUAL_BALANCED`, and only with `CostOnly` or `Pinned`
  routing). Each call's offer is split at every node by the kernel's rule; each edge carries the
  smaller of its share and its residual, and the rest is dropped on that edge. The delivery is
  what reaches the destination set. An infinite offer injects into each first edge exactly its
  residual. A share that has no edge to go on (a `WEIGHTED` node whose members all have weight 0)
  is lost where it is: at a source it is `source_loss`. Records keep `offered`, `delivered` and
  per-edge `drops`. Lossy calls are independent: nothing is re-offered and two calls are not
  equivalent to one.

On a diamond with paths of capacity 10 and 2, equal-share forwarding stops at 5 units: round
robin gives the odd unit to the lower EdgeId (3 and 2), and a sixth unit would put 3 units on the
2-unit path. The weighted and maximum-flow kernels carry 12:

```python
g = graph(4, [(0, 1, 10, 1), (1, 3, 10, 1), (0, 2, 2, 1), (2, 3, 2, 1)])
kernels = ngc.FlowPlacement.EQUAL_BALANCED, ngc.FlowPlacement.WEIGHTED, ngc.FlowPlacement.MAX_FLOW
for kernel in kernels:
    ledger = ngc.FlowLedger(g)
    r = ngc.FlowPolicy(ledger, 0, 3, config(ngc.CostOnly(), placement=kernel)).place(100)
    print(kernel.name, r.placed, r.remaining, r.termination.name, ledger.allocated().tolist())

ledger = ngc.FlowLedger(g)
pol = ngc.FlowPolicy(ledger, 0, 3, config(ngc.CostOnly(), admission=ngc.Admission.LOSSY))
r = pol.place(10)
print("LOSSY", r.placed, r.remaining, r.termination.name, "drops", r.call_drops)
rec = pol.records()[0]
print("record", rec.offered, rec.delivered, rec.drops)
```

```text
EQUAL_BALANCED 5 95 EXHAUSTED [3, 2, 3, 2]
WEIGHTED 12 88 EXHAUSTED [10, 2, 10, 2]
MAX_FLOW 12 88 EXHAUSTED [10, 2, 10, 2]
LOSSY 7 3 COMPLETE drops [(1, 3)]
record 10 7 [(1, 3)]
```

The Lossy call is `COMPLETE` although 3 units were dropped: the whole offer of 10 was forwarded,
and what did not fit is accounted as drops, not as a shortfall ([section 7](#outcomes)).

An edge of capacity 0 on the only shortest route is still that route. Equal splitting sends it
the whole offer: Lossless places nothing, Lossy drops all of it there. Weighted splitting gives
it weight 0 and has no other edge in its route, so nothing can be sent (`NO_ACTIVE_WEIGHT`) and
a Lossy offer is lost at the source; neither takes the costlier detour. A controller (`Global`)
uses the detour:

```python
# Edge 0: 0 -> 1 of capacity 0 (cost 1); the detour 0 -> 2 -> 1 (5 units, cost 2).
g = graph(3, [(0, 1, 0, 1), (0, 2, 5, 1), (2, 1, 5, 1)])
r = ngc.spf(g, 0, ngc.SpfOptions(F, DET, targets=1))
print("spf", r.target_cost, r.dag.via_edges().tolist())
for kernel in (ngc.FlowPlacement.EQUAL_BALANCED, ngc.FlowPlacement.WEIGHTED):
    for admission in (ngc.Admission.LOSSLESS, ngc.Admission.LOSSY):
        cfg = config(ngc.CostOnly(), placement=kernel, admission=admission)
        r = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 1, cfg).place(5)
        print(kernel.name, admission.name, r.placed, r.call_drops, r.call_source_loss,
              r.termination.name)
cfg = ngc.FlowPolicyConfig(routing=ngc.Global())
print("Global", ngc.FlowPolicy(ngc.FlowLedger(g), 0, 1, cfg).place(5).placed)
```

```text
spf 1 [0]
EQUAL_BALANCED LOSSLESS 0 [] 0 EXHAUSTED
EQUAL_BALANCED LOSSY 0 [(0, 5)] 0 COMPLETE
WEIGHTED LOSSLESS 0 [] 0 NO_ACTIVE_WEIGHT
WEIGHTED LOSSY 0 [] 5 NO_ACTIVE_WEIGHT
Global 5
```

## 6. The ledger

**In plain words.** A `FlowLedger` holds the traffic placed on one graph: a starting residual
(the baseline) and one record per flow. Policies place and remove its flows; it has no public
placement or removal of its own. Its own `reset()`, rejected while a policy is bound, replaces
the baseline and clears the records. Its residual is always the baseline minus everything its
records hold.

**Exact rules.**

- **`FlowLedger(graph, residual=None)`.** The baseline is an integer array in EdgeId order with
  every entry in `[0, capacity_units[e]]` (a physical baseline goes through `to_units` first;
  anything above the capacity is rejected); `None` means the capacities. The baseline never
  changes except through `reset`.
- **Snapshots.** `allocated()` is the exact sum of every record's edges and `residual()` is
  `baseline - allocated()`, both `int64[E]` copies. `capacity_units()` is the graph's view.
  `graph()` returns the graph.
- **Records.** `records()` lists the live records (of bound policies and orphans) as snapshots
  in FlowId order; `record(flow_id)` reads one. A FlowId is the pair `(ledger uid, id)`; ids are
  never reused. A FlowId of another ledger, or a retired or unknown id, raises `ValueError`.
- **`FlowRecord`** fields: `flow_id`; `owner` (the policy's id); `sources` and `destinations`
  (node sets); `admission`; `bundle` (the flow's route as a `PredDAG`);
  `route_cost` (the cheapest walk of the bundle, cached; `None` for an empty bundle); `offered`
  (Lossy records only, else `None`); `delivered`; `edges` (`[(EdgeId, units)]`, ascending) and
  `drops` (`[(EdgeId, units)]`, cumulative, ascending).
- **`bottleneck(sources, elig=None)`** is a saturation report: the edges (eligible, at least 1
  unit of capacity) that leave the set of nodes reachable from the allowed sources over edges
  with at least 1 unit of residual, in ascending EdgeId. It follows edges forward only. It ignores
  the destinations, so it is not a source-destination cut, and it can be empty while a route
  exists.
- **Lifetimes.** The graph outlives the ledger and the ledger outlives its policies (Python keeps
  them alive for you).

| Ledger state | Reached by | `reset(residual)` | `copy()` |
| --- | --- | --- | --- |
| a policy is bound | building a policy on it | `ValueError` | `ValueError` |
| no policy bound | no policy built yet, or every policy destroyed | allowed: orphaned records and their allocations are cleared, their ids retired, the new baseline set | a new ledger with a new uid and the same records and ids |

A policy destroyed without `remove()` leaves its records and allocations in the ledger as
**orphans**: readable, never changed again, cleared only by `reset`.

## 7. Policies

**In plain words.** A `FlowPolicy` is the routing model of one demand between a set of sources
and a set of destinations. `place(volume)` asks it to carry more; it routes, fills and possibly
grows its flows, commits everything at once, and reports how much was placed, what stopped it,
where traffic was lost and how much each source carried. `remove()` takes all of its traffic
back out. A demand with several sources is one demand, not several: how its traffic enters at
the sources is its **origin**, and where the routes of its sources meet, their traffic is added
and split once.

### Configuration

`FlowPolicyConfig(*, routing, placement=None, admission=Admission.LOSSLESS, selection=None,
origin=None)` takes keywords only. The routing is one of:

| Routing | Models | Fields (defaults) |
| --- | --- | --- |
| `CostOnly()` | IGP forwarding (as OSPF or IS-IS routers do): static shortest paths, never refreshed; with `MAX_FLOW`, a controller restricted to those paths | none |
| `CapacityAware(...)` | distributed CSPF tunnels (constrained shortest path first: each tunnel head picks a shortest path with enough free capacity) | `flow_count=1`, `max_flow_count=None`, `max_path_cost=None`, `min_flow_bandwidth=0`, `reoptimize=False`, `diminishing=None` |
| `Global()` | an ideal controller that sees the whole network | none |
| `Pinned(bundles)` | explicit routes given by the caller | a list of `PredDAG` (a `Path` converts with `to_dag()`); copied when the policy is built |

Each routing accepts exactly the combinations of this table; any other raises `ValueError` before
any work, a request of 0 included. A field marked "none" must be left `None`:

| Routing | `placement` | `origin` default / allowed | `admission` | `selection` |
| --- | --- | --- | --- | --- |
| `CostOnly` | `EQUAL_BALANCED` (ECMP) or `WEIGHTED` (WCMP) | `EQUAL` / `EQUAL` | `LOSSLESS`, `LOSSY` | all four fields |
| `CostOnly` | `MAX_FLOW` | `OPTIMIZED` / both | `LOSSLESS` | all four fields |
| `Global` | none | `OPTIMIZED` / both | `LOSSLESS` | none |
| `CapacityAware` | any kernel (required) | `EQUAL` / `EQUAL` | `LOSSLESS` | all four fields |
| `Pinned` | any kernel (required) | none | `LOSSLESS`, `LOSSY` (not with `MAX_FLOW`) | only `priority`, with `WEIGHTED` |

`selection=None` means `PathSelection()` (`DETERMINISTIC`, `CAPACITY`, `multipath=True`,
`multi_edge=True`) for the routings that read it, and `origin=None` means the routing's default.
A policy's `config()` returns the configuration with both filled in (`selection` stays `None` for
`Global` and `origin` for `Pinned`). Further rules:

- `WEIGHTED` with priority `HEADROOM` is rejected.
- `flow_count >= 1`; `max_flow_count` is `None` (no ceiling) or at least `flow_count`;
  `min_flow_bandwidth` is in [0, 2^62) units; `max_path_cost` is `None`, `Absolute(c)` with
  `c >= 0` or `Stretch(f)` with `f` finite and at least 1 (as for `ksp`); `diminishing` is `None`
  (off) or `DiminishingReturns(window, fraction)` with `window >= 1` and `fraction` strictly
  between 0 and 1.
- Only `CostOnly` reads a `DagCache`; any other routing given one is rejected.
- `Pinned`: one source and one destination; at least one bundle; no duplicate bundles; every
  bundle runs from the policy's source to its destination and has at least one walk.
- **The pooled single-path rule.** `CostOnly` with more than one source (counted after removing
  repeats, before any Eligibility) rejects `multipath=False` together with
  `PREFER_HIGHER_PRIORITY`: the bottleneck choice could then differ per source at a router the
  sources share.

`FlowPolicy(ledger, sources, destinations, config, elig=None, cache=None)` binds the policy to
its ledger. `sources` and `destinations` are one node id or a sequence; the two sets must not
overlap. Construction validates everything before the policy binds (a field that a routing does
not have is a `TypeError`; every other defect a `ValueError`), and when several defects meet, the
first in this order is reported: the endpoint sets (empty, out of range, the lowest node in both);
a foreign `Eligibility` or `DagCache`; enum values outside their enumeration; the table above,
its row (routing, then placement) and then its columns left to right (placement, origin,
admission, selection); the cache rule; `CapacityAware`'s own fields; the endpoint counts
(`Pinned`'s single endpoints, the pooled single-path rule); `Pinned`'s bundles.

### The origin

The origin says how a demand's total enters at its roster sources, separately from how it is
routed:

- **`EQUAL`**: a total `x` enters as the round robin `RR(x, k, i) = x // k + (1 if i < x % k else
  0)` at roster member `i` (`k` members, ascending node id): the first members get one extra
  unit. A source whose share cannot go anywhere keeps it and loses it; the share is never moved
  to another source. Lossless therefore admits only totals that give such a source nothing, and
  Lossy loses its share.
- **`OPTIMIZED`** (controllers only): the controller chooses how much enters at each source (one
  virtual source joined to every roster source), whatever maximises the delivery and then
  minimises its cost.

There is no weighted origin; sources with fixed or uneven amounts are separate demands.

### The source roster

The roster is fixed when the policy is built. For `CostOnly`, `CapacityAware` and `Global` it is
the allowed sources that have a static walk to an allowed destination over allowed edges, whatever
their capacity, in ascending node id: a source whose only routes cross an edge of capacity 0 is in
the roster, receives its share and cannot place it (a failure is expressed with an Eligibility,
which removes the source). For `Pinned` it is the single source when at least one bundle is
usable. The roster decides which sources receive shares and which records exist; `by_source`
reports every input source, a source outside the roster with a zero row. A positive request with
an empty roster creates nothing and ends `NO_ROUTE`.

### Routing modes

- **`CostOnly` with `EQUAL_BALANCED` or `WEIGHTED`** forwards like IGP routers. At the epoch's
  first positive call it prepares a **program**, the routes and splits used for the whole epoch:
  each roster source's route is its static shortest-path DAG toward the destination set (with
  `multipath=False`, its canonical path), read from the [DAG cache](#9-the-dag-cache), with the
  selection applied on the ledger's residual at that call's entry (priorities read from that
  residual for `RESIDUAL` and `HEADROOM`, from the static capacities for `CAPACITY`; `WEIGHTED`
  freezes its weights there). The program is the union of these routes: at a router several
  sources pass, their traffic for the demand is added (with the router's own injection when it
  is a source) and split once over its members, by round robin (`EQUAL_BALANCED`) or by highest
  averages over the weights (`WEIGHTED`). Every source that passes a router has the same
  continuation there, so the union lends no source a link outside its own route. The program is
  never refreshed, grown or reoptimised within the epoch.
- **`CostOnly` with `MAX_FLOW`** is a controller restricted to shortest paths: it solves like
  `Global` (below), but only over the union of the roster sources' routes prepared as above,
  edges of capacity 0 included (so a zero-capacity shortest route is never bypassed).
- **`CapacityAware`** routes by the canonical path or the DAG of [section 4](#4-search) on the
  current free capacity, excluding edges with less than `max(1, min_flow_bandwidth)` units free and
  routes above `max_path_cost` (`Stretch` is relative to the static shortest cost, computed when the
  policy is built). A source first creates `flow_count` flows, one search at a time, each search
  seeing the capacity the earlier ones reserved; with `reoptimize=True` it moves a retained flow to
  a newly found route only if the flow keeps its whole amount there; while the offer is not met it
  grows its cohort by one flow at a time, up to `max_flow_count`. Growth stops when the new route is
  a multipath DAG equal to an existing flow's route (it would add nothing); a repeated single path
  is a real candidate (tunnels of equal share fill a fatter path with more tunnels). `diminishing`
  stops growth once at least `window` growth steps were accepted in the call and the last `window`
  gains sum below `fraction` times the source's (finite) offer. A call creates at most
  `MAX_FLOWS_PER_CALL` (10,000) flows over all sources. This mode is greedy: route order matters and
  the result is not an upper bound.
- **`Global`** re-solves, on every positive call, the minimum-cost flow of the amount the
  policy already delivers plus this call's request (`placed_demand()` plus the request; what an
  earlier call could not place is not requested again) from its roster to its whole
  destination set over every eligible edge, on the free capacity plus its own allocation (it may
  move its own traffic, never another owner's): the largest delivery up to that amount, and among
  those the cheapest. Hence `place(a)` followed by `place(b)` ends in the same allocation as
  `place(a + b)` on an unchanged network. With `EQUAL` origin the virtual source's arc to roster
  source `i` carries exactly `RR(T, k, i)`: the delivered total `T` is the largest whose equal
  shares all fit (found by exact binary search over the totals; a total fits only when the solve
  reaches it), and that total is solved once more from zero, so the result does not depend on the
  search. A positive call on a network beyond the flow core's size limit raises `ValueError`
  ([section 2](#2-the-graph), Size).
- **Records of the controllers** (`Global`, `CostOnly` with `MAX_FLOW`): one record per roster
  source, holding that source's share of the solved flow (its edges and its delivery), extracted
  in a canonical order. A `Global` record's bundle is the set of edges that carry the share,
  rebuilt on every positive call under the same FlowId; a source with no share has an empty
  bundle and `route_cost` `None`. A `CostOnly` `MAX_FLOW` record's bundle is its source's
  prepared route, and the share stays inside it.
- **`Pinned`** creates one flow per usable bundle, in the order given. A bundle is usable when at
  least one of its walks survives the Eligibility, whatever the capacities of its edges; its
  `route_cost` is the cheapest such walk. Bundles that become equal after that pruning keep separate flows. When no bundle is
  usable, `usable_bundles()` is 0 and a positive request ends `NO_ROUTE`.

### Several sources and destinations

- **Destinations** (`CostOnly`, `CapacityAware`): a route goes to the nearest destinations, and
  one flow may deliver to several of them at that cost. A record's `delivered` is its total; the
  delivery to one destination is the sum of the record's `edges` entering it.
- **`CostOnly` split routing pools its sources.** Each record (one per roster source) holds the
  part of the program's traffic that is that source's. At every router the members' amounts are
  given back to the sources member by member, in member order, by highest averages over what each
  source still holds there; on an edge that drops traffic (Lossy), what the edge carries is
  divided by highest averages over what each source brought to it. Highest averages `H(x, w)`
  gives `x` units one at a time to the member with the largest `w[i] / (assigned[i] + 1)`, ties to
  the earlier member; a zero weight gets nothing. So the records sum exactly to the totals, and a
  source's record uses only edges of its own route.
- **The controllers** (`Global`, `CostOnly` with `MAX_FLOW`) optimise the whole roster and
  destination set together; `Global` may use a costlier destination to carry more.
- **`CapacityAware` sources** act as independent tunnel heads. The policy's cumulative request is
  split over the roster by round robin in ascending node id, and each source is offered the
  increase of its share in this call (an infinite request offers infinity to every source).
  Sources are placed one after another, in ascending node id, on the shared free capacity; what
  one source cannot place is not passed to another. A split request is therefore not equivalent
  to one large request.
- **`Pinned`** has a single source and a single destination.

Two sources with one unit each meet at node 2, whose two members lead to the destination. Pooled,
node 2 splits the two units once, one per member, and gives one member to each source; split
separately, each source would send its unit to the first member:

```python
# Edges 0: 0 -> 2, 1: 1 -> 2, then the two members 2: 2 -> 3 and 3: 2 -> 3.
g = graph(4, [(0, 2, 5, 1), (1, 2, 5, 1), (2, 3, 1, 1), (2, 3, 1, 1)])
pol = ngc.FlowPolicy(ngc.FlowLedger(g), [0, 1], 3, config(ngc.CostOnly()))
r = pol.place(2)
print(r.placed, r.termination.name, [(row.source, row.delivered) for row in r.by_source])
print([(rec.sources, rec.edges) for rec in pol.records()])
```

```text
2 COMPLETE [(0, 1), (1, 1)]
[([0], [(0, 1), (2, 1)]), ([1], [(1, 1), (3, 1)])]
```

### The equal origin and a source that cannot pass anything

A roster source whose every route crosses an edge with no free capacity stays in the roster: it
receives its share and loses it. Lossless admits only totals that give it nothing, which is the
number of roster sources before it, because the first remainder units go to the lowest node ids:
0 when it is the first source. Lossy loses its share where its program puts it: on the edge that
limits it (`drops`), or at the source when it has no member to send it on (`source_loss`). An
Eligibility that excludes the source or the edge removes it from the roster instead.

```python
# Source 0 reaches node 2 only over the 0-unit edge 0; source 1 over the 10-unit edge 1.
g = graph(3, [(0, 2, 0, 1), (1, 2, 10, 1)])
equal = ngc.FlowPolicyConfig(routing=ngc.Global(), origin=ngc.Origin.EQUAL)
for name, cfg in (("CostOnly", config(ngc.CostOnly())),
                  ("Global EQUAL", equal),
                  ("Global OPTIMIZED", ngc.FlowPolicyConfig(routing=ngc.Global()))):
    r = ngc.FlowPolicy(ngc.FlowLedger(g), [0, 1], 2, cfg).place(10)
    print(name, r.placed, r.termination.name)
lossy = config(ngc.CostOnly(), admission=ngc.Admission.LOSSY)
r = ngc.FlowPolicy(ngc.FlowLedger(g), [0, 1], 2, lossy).place(10)
print("LOSSY", r.placed, r.call_drops, [(row.source, row.dropped) for row in r.by_source])
r = ngc.FlowPolicy(ngc.FlowLedger(g), [0, 1], 2, config(ngc.CostOnly()),
                   g.eligibility(excluded_nodes=[0])).place(10)
print("source 0 excluded", r.placed)
```

```text
CostOnly 0 EXHAUSTED
Global EQUAL 0 EXHAUSTED
Global OPTIMIZED 10 COMPLETE
LOSSY 5 [(0, 5)] [(0, 5), (1, 0)]
source 0 excluded 10
```

### The placement call

`place(volume, outputs=None)` takes a request in units (an `int`) or `float('inf')`.

1. **Validation** happens before anything changes. Within an epoch, a finite request after an
   infinite one is rejected (the cumulative request would be undefined), `place(0)` included.
2. **A request of 0** is `COMPLETE` and does nothing: it creates no flow and prepares nothing. It
   reports the demand's state, with both flags false.
3. **`CostOnly` split routing** prepares the program at the epoch's first positive call (one
   record per roster source, also at zero delivery). Lossless then admits the largest total `a`
   between what it already delivers and that plus the request whose forwarding of `RR(a)` from
   zero fits the free capacity plus the policy's own allocation with no drop and gives no share to
   a source without a member (exact binary search: the totals that fit form a prefix). Lossy
   injects at the sources this call's increase of the round robin, `RR(offered + v) -
   RR(offered)`, and forwards it hop by hop on the free capacity: each edge carries the smaller of
   what reaches it and its residual and drops the rest. An infinite Lossy request is the
   **first-hop experiment**: each source injects on each of its first edges exactly that edge's
   residual, on top of the traffic other sources send through it, and forwarding runs as usual.
4. **The controllers** re-solve as described above.
5. **`CapacityAware` and `Pinned`** run, per source and on one shared working copy of the free
   capacity: creation (while the source has no flows), reoptimisation (`CapacityAware` with
   `reoptimize`), filling (each flow, or each equal-share cohort, raised toward its maximum),
   growth (`CapacityAware`).
6. **Commit.** All changes to the ledger and the policy are committed together at the end, or
   none are ([section 8](#8-failures-threads-and-lifetimes)).

What the call guarantees:

- `placed` is the increase in delivery made by this call, and `0 <= placed <= volume`. For Lossy
  it is this call's delivery. `remaining = volume - placed` (`float('inf')` for an infinite
  request; for Lossy it includes what was lost).
- A call never reduces `placed_demand()`, the sum of `delivered` over the policy's records. The
  total of an equal-share cohort never decreases; one tunnel's amount can decrease when growth
  re-splits a larger total over more tunnels (`[3, 3]` plus a new tunnel can become
  `[3, 2, 2]`).
- The first positive request creates the policy's records even when none of them can deliver
  anything yet (`CapacityAware` only creates a flow on a route it finds). Such a zero-delivery
  record keeps its FlowId for the rest of the epoch and can fill on a later call.
- **Additivity.** With the rest of the network unchanged between the calls, `place(a); place(b)`
  equals `place(a + b)` exactly for Lossless `CostOnly` split routing, the controllers (value,
  allocation and cost curve) and a single Lossless `Pinned` flow or cohort. Lossy calls are not
  additive: each forwards its own injection on the free capacity the earlier calls left, and
  nothing is offered again. `CapacityAware` sources are not additive.

### Epochs

`remove()` takes every flow of the policy out of the ledger (allocations reverted, records gone,
ids retired) and starts a new epoch: the cumulative request, the round-robin position of the
source split, the finite-or-infinite state and the prepared program are reset. The next positive
request prepares again and creates new flows with new ids; ids are never reused.

```python
g = graph(4, [(0, 1, 10, 1), (1, 3, 10, 1), (0, 2, 2, 1), (2, 3, 2, 1)])
ledger = ngc.FlowLedger(g)
pol = ngc.FlowPolicy(ledger, [0, 1], 3, config(ngc.CostOnly()))
for v in (1, 1, 0):
    r = pol.place(v)
    print(v, r.placed, r.routes_changed, [(x.flow_id[1], x.delivered) for x in pol.records()])
pol.remove()
r = pol.place(1)
print(1, r.placed, r.routes_changed, [(x.flow_id[1], x.delivered) for x in pol.records()])
```

```text
1 1 True [(0, 1), (1, 0)]
1 1 False [(0, 1), (1, 1)]
0 0 False [(0, 1), (1, 1)]
1 1 True [(2, 1), (3, 0)]
```

The first unit goes to source 0 and the second to source 1 (round robin of the cumulative
request); the first call creates both records, source 1's at zero delivery. After `remove()` the
split starts again at source 0, under new ids.

### Termination

The call's `termination` is the first row that applies (an invalid call raises instead):

| Termination | When |
| --- | --- |
| `COMPLETE` | the request is 0 |
| `NO_ROUTE` | the request is positive and the roster is empty (`Pinned`: no usable bundle) |
| `NO_ACTIVE_WEIGHT` | routes exist, but none has an edge it may use: every `WEIGHTED` member has weight 0, so the program sends nothing anywhere |
| `COMPLETE` | Lossy with a finite request: the whole request was forwarded, and delivery and losses account for all of it (it does not mean nothing was lost) |
| `COMPLETE` | Lossless: the whole request was admitted |
| `EXHAUSTED` | otherwise: Lossless below the request, any infinite request (the first-hop experiment included), a controller whose maximum is 0 over existing routes |

`CapacityAware` keeps its own reduction over its sources. Each roster source ends the call
served (its offer is 0 or was delivered) or with the first local outcome that applies:
`LIMIT_FLOWS` (creation or growth was stopped by `MAX_FLOWS_PER_CALL`, or by `max_flow_count` when
`max_flow_count > flow_count`), `LIMIT_DIMINISHING` (the `diminishing` guard stopped growth),
`NO_ROUTE` (the source has no flows after the call), else `EXHAUSTED` (every flow at its maximum and
growth found nothing that helps, a costlier helpful route may still exist; or a fixed-size cohort
is full). A request of 0 is `COMPLETE` and an empty roster `NO_ROUTE`; otherwise the call is
`COMPLETE` when every source is served, else the highest outcome among the unserved sources in the
order `LIMIT_FLOWS` > `LIMIT_DIMINISHING` > `EXHAUSTED` > `NO_ROUTE`. `LIMIT_FLOWS` is reported when
the ceiling is reached before the offer is met, even if no further route would have existed.

`remaining` is an amount, never a status: a Lossy call can be `COMPLETE` with a positive
`remaining`.

```python
# Source 0 has two disjoint 1-unit tunnels to node 4; source 1 one 1-unit edge.
g = graph(5, [(0, 2, 1, 1), (0, 3, 1, 1), (2, 4, 1, 1), (3, 4, 1, 1), (1, 4, 1, 1)])
tunnel = ngc.PathSelection(ngc.EdgeTieBreak.DETERMINISTIC, ngc.EdgePriority.RESIDUAL,
                           multipath=False)

def show(r):
    print(r.placed, r.termination.name, [(row.source, row.delivered) for row in r.by_source])

# Request 10: each source is offered 5. Source 0 stops at its ceiling of 2 tunnels.
pol = ngc.FlowPolicy(ngc.FlowLedger(g), [0, 1], 4,
                     config(ngc.CapacityAware(max_flow_count=2), ngc.FlowPlacement.MAX_FLOW,
                            selection=tunnel))
show(pol.place(10))

# Weighted splitting by residual where another demand holds the only edge: no positive weight.
ledger = ngc.FlowLedger(g)
ngc.FlowPolicy(ledger, 1, 4, config(ngc.CostOnly())).place(1)
weighted = config(ngc.CostOnly(), ngc.FlowPlacement.WEIGHTED,
                  selection=ngc.PathSelection(priority=ngc.EdgePriority.RESIDUAL))
show(ngc.FlowPolicy(ledger, 1, 4, weighted).place(1))
```

```text
3 LIMIT_FLOWS [(0, 2), (1, 1)]
0 NO_ACTIVE_WEIGHT [(1, 0)]
```

### Flags

`routes_changed` and `weights_changed` start false in every call and record only changes that
were committed. Once set in a call, a flag stays set; a call with no committed structural change
ends with both false, even when it places nothing.

| Committed change | `routes_changed` | `weights_changed` |
| --- | --- | --- |
| records are created, in any routing mode, also at zero delivery and in a new epoch after `remove()` (`CapacityAware` only with a route, possibly on a later call) | true | unchanged |
| a retained flow is moved to a new route (`reoptimize`) | true | true when a `WEIGHTED` flow's weights were captured again (also when the call places nothing) |
| growth adds a flow (also on a route equal to an existing flow's: routes belong to flows) | true | unchanged |
| a `Global` call changes the edges some source's share uses, also with `placed == 0` | true | a controller never sets it |
| nothing of the above: a request of 0, rejected candidates only, `Global` with every source's edges unchanged | unchanged | unchanged |

A candidate route equal to the flow's current route is not a change.

```python
# A first request on a route whose only unit another demand holds.
g = graph(5, [(0, 2, 1, 1), (0, 3, 1, 1), (2, 4, 1, 1), (3, 4, 1, 1), (1, 4, 1, 1)])
ledger = ngc.FlowLedger(g)
ngc.FlowPolicy(ledger, 1, 4, config(ngc.CostOnly())).place(1)
pol = ngc.FlowPolicy(ledger, 1, 4, config(ngc.CostOnly()))
r = pol.place(1)
print(r.placed, r.termination.name, r.routes_changed, [x.delivered for x in pol.records()])
```

```text
0 EXHAUSTED True [0]
```

### The result

`place` returns a `PlacementResult`. Its fields describe either this call or the demand's state
after it (totals over the epoch); a per-call loss field carries the `call_` prefix, a state field
does not, so the two are never confused.

- **This call**: `placed`, `remaining`, `termination`, `routes_changed` and `weights_changed` as
  above; `delivery_delta`, the change of the delivery (for `place` it equals `placed`);
  `call_drops`, this call's drops per edge as `[(EdgeId, units)]`; `call_source_loss`, the units
  this call lost at sources.
- **The demand's state**: `requested` (the sum of the requests, `float('inf')` after an infinite
  one), `offered`, `admitted`, `delivered`, `unadmitted`, `source_loss`, `drops` (per edge, as
  `[(EdgeId, units)]`), `by_source` and `total_cost`. Every amount is an exact Python integer and
  may exceed 64 bits.
- **Accounting.** Lossless: `offered` is `None`, `source_loss` 0, no drops, `admitted = delivered`,
  and `requested = unadmitted + delivered` for a finite request (`unadmitted` is `None` after an
  infinite one). Lossy: `offered` is the traffic actually sent, `admitted = offered`, `unadmitted =
  0`, and `offered = source_loss + sum(drops) + delivered`; `offered = requested` except in the
  first-hop experiment, where `requested` is infinite and `offered` the finite injection.
  `delivered + sum(drops)` is what entered the network.
- **`by_source`** is a tuple of `SourceFlow` rows, one per input source in ascending node id, a
  source outside the roster included (a zero row): `source`, `offered` (`None` when Lossless),
  `admitted`, `delivered`, `source_loss` and `dropped` (its share of `drops`). The rows sum to the
  totals, with one exception: a Lossy positive request with an empty roster reports `offered =
  requested` and all of it as `source_loss`, lost before entering, while every row stays zero.
- **`total_cost`**: the controllers only, the exact `sum(cost[e] * flow[e])` of the policy's
  allocation; `None` otherwise.
- **`cost_distribution`**: `Global` with `OPTIMIZED` origin only, the marginal cost curve of the
  policy's cumulative flow (a request of 0 reports it unchanged; `None` before the epoch's first
  positive call). `None` for every other mode, `EQUAL` origin included: the marginal costs of an
  equal flow need not increase. For other modes, the sum of `(route_cost, delivered)` over the
  records is a histogram of delivery by each record's cheapest walk cost: the cost every unit
  travels at for `CostOnly` split routing and `CapacityAware`, whose routes have one cost, but not
  for a `Pinned` bundle whose walks differ in cost.
- **Optional outputs**, `None` unless requested through `OutputOptions(edge_flows=False,
  residual=False, bottleneck=False, min_cut=False)`: `edge_flows` is the policy's own allocation
  per edge after the call (the ledger's total is `allocated()`); `residual` and `bottleneck` are
  the ledger's after the call; `min_cut` is described below.
- **`bottleneck` and `min_cut` are different things.** `bottleneck` is the saturation report of
  [section 6](#6-the-ledger), available in every mode. `min_cut` exists only for `Global` with
  `OPTIMIZED` origin: the edges of a minimum cut of the policy's view, the one on the sources' side:
  the allowed edges that leave the nodes the sources can still reach once the flow is maximal
  (counting the option to undo flow on an edge), in ascending EdgeId, an edge with nothing free in
  the view (capacity 0, or held by other policies) included: it adds 0 to the cut's value and is
  listed because every allowed edge leaving that set is listed; `[]` for an empty roster. It is
  present only when the solve reached the maximum on its view; when a request below the maximum
  stopped it, on a request of 0 (nothing is solved), and in every other mode, `min_cut` is `None`.

```python
g = graph(4, [(0, 1, 10, 1), (1, 3, 10, 1), (0, 2, 5, 2), (2, 3, 5, 2)])
ledger = ngc.FlowLedger(g)
pol = ngc.FlowPolicy(ledger, 0, 3, ngc.FlowPolicyConfig(routing=ngc.Global()))
want = ngc.OutputOptions(bottleneck=True, min_cut=True)
r = pol.place(12, want)
print(r.placed, r.cost_distribution, r.min_cut, r.bottleneck, r.termination.name)
r = pol.place(float("inf"), want)
print(r.placed, r.cost_distribution, r.min_cut, r.bottleneck, r.termination.name)
print(r.requested, r.admitted, r.unadmitted, r.total_cost)
```

```text
12 [(2, 10), (4, 2)] None [0] COMPLETE
3 [(2, 10), (4, 5)] [0, 1] [0, 1] EXHAUSTED
inf 15 None 40
```

After 12 units the cheap path (cost 2) is full and the solve stopped below the maximum, so there
is no cut; `bottleneck` already reports the full edge 0. The second call adds 3 units at cost 4
and reaches the maximum of 15, at a total cost of 10 x 2 + 5 x 4 = 40.

## 8. Failures, threads and lifetimes

**In plain words.** A call either does everything or changes nothing: if it fails, for any
reason, every ledger, policy and cache is as it was and you get one exception. There are two
stated exceptions to this, described below. Immutable values can be shared freely between
threads; ledgers, policies and building caches have one writer at a time.

**Exact rules.**

- **All or nothing.** A failure anywhere in a public call (validation, memory allocation, a
  worker thread, the construction of the Python result) leaves every ledger, policy, cache and
  internal buffer as it was before the call and raises one exception in the caller's thread.
  Worker threads are always joined first; an exception never escapes a worker. When several pairs
  or candidates of an evaluator fail, the exception of the lowest index is raised.
- **Limit 1: pybind11's own object machinery.** Constructing a bound object from Python (any of
  the library's classes) and converting a returned value into a new Python object or container
  outside the library's checked result construction (for example `PredDAG.from_csr`,
  `Graph.eligibility`, `SensitivityResult.baseline`, `FlowPolicy.records()`) use pybind11's own
  allocation paths, which do not all check for failure. A memory failure at exactly that point
  may crash the process instead of raising. Library state is never left half-changed by it. The
  results of `spf`, `ksp`, `resolve_to_paths`, `FlowPolicy.place`, `capacity`, `batch_capacity`
  and `sensitivity` are built through a checked path and stay inside the guarantee.
- **Limit 2: a `DagCache` passed to the call.** A building cache is a memo: when a call fails
  after some cache lookups succeeded (a policy constructor failing after its lookup, an evaluator
  failing midway), the cache keeps those complete entries and their statistics. It never holds a
  partial entry, every value it handed out stays valid, and a kept entry never changes any later
  result. A failed lookup counts one miss. A frozen cache never changes its entries; a failed
  call changes none of its statistics, except one miss when the call raised `DagCacheMiss`.
- **Threads.** `Graph`, `Eligibility`, `PredDAG`, `Path` and a frozen `DagCache` may be shared
  by any number of threads. A `FlowLedger`, a `FlowPolicy` and a building `DagCache` have one
  writer and no concurrent readers. An array lent to a call (a residual) must not change during
  the call.
- **Lifetimes.** The graph outlives the ledger, the ledger outlives its policies; Python keeps
  the graph, ledger, Eligibility and cache alive while a policy needs them.
- **Exception classes.** `ValueError` for a value outside its domain, an invalid combination, a
  foreign object (a value of another graph, ledger or policy) or a call in the wrong state;
  `TypeError` for a field or option that does not exist (for example `CostOnly(flow_count=2)`),
  an argument of the wrong type, or an argument that cannot be converted to its parameter's type:
  an integer beyond the float range where a float is expected (the `resolution` of `Graph`,
  `to_units` and `as_float`) is a `TypeError`. The fields that convert their own values (the
  factor of `Stretch`, the fraction of `DiminishingReturns`, the request of `place`) raise
  `ValueError` for such an integer and pass on unchanged any exception raised by the value's own
  conversion (its `__float__` or `__index__`), whatever its class. `DagCacheMiss` for a miss in
  a frozen cache (its own class, neither a `ValueError` nor a `RuntimeError`).

```python
g = graph(4, [(0, 1, 10, 1), (1, 3, 10, 1), (0, 2, 2, 1), (2, 3, 2, 1)])
ledger = ngc.FlowLedger(g)
pol = ngc.FlowPolicy(ledger, 0, 3, config(ngc.CostOnly()))
print(pol.place(float("inf")).placed)
try:
    pol.place(1)          # a finite request after +inf in the same epoch
except ValueError:
    print("ValueError; unchanged:", pol.placed_demand(), ledger.allocated().tolist())
```

```text
5
ValueError; unchanged: 5 [3, 2, 3, 2]
```

## 9. The DAG cache

**In plain words.** Many placements need the same static shortest-path DAG toward the same
destinations. A `DagCache` remembers those DAGs for one graph. It stores the full DAG and applies
each caller's selection afterwards, so one entry serves every selection. A cache is either
building (one writer adds entries) or frozen (read-only, shared by many threads, and a miss is an
error).

**Exact rules.**

- **Binding and key.** `DagCache(graph)` is bound to one graph (`graph_uid()`). The key is
  `(the root set as a sorted set, the direction, the Eligibility uid)`. The value is the full
  static shortest-path DAG: the DAG with `multipath=True`, `multi_edge=True` and no residual. Values are
  immutable and stay valid after `clear()` and after the cache is destroyed.
- **Selection is never cached.** The single path and the member rule are applied when a flow is
  created, on the cached DAG, reading priorities from the residual at that moment for
  `RESIDUAL` and `HEADROOM` and from static capacities for `CAPACITY`. Because the cached DAG keeps
  every equal-cost parallel member, a member that a static choice would have dropped is still
  available when its residual makes it the right one.
- **Who reads it.** `CostOnly` policies read the cache (keyed by their destination set, reverse
  direction, and Eligibility uid); so do the evaluators when they run a `CostOnly` configuration.
  The other routing modes do not read it and reject a supplied cache (`ValueError`, for example
  `Global does not use DagCache`).

| Operation | Building | Frozen |
| --- | --- | --- |
| lookup hit | returns the value | returns the value; safe from many threads |
| lookup miss | searches, inserts, returns | raises `DagCacheMiss` |
| `freeze()` | becomes frozen | no effect |
| `clear()` | empties the cache; hit and miss counters persist | `ValueError` |
| `thaw()` | no effect | becomes building; join every reader first |
| destruction | allowed | allowed after joining every reader; values already handed out stay valid |

- **Statistics.** `stats()` returns `hits`, `misses`, `entries`, `bytes` (the entries' retained
  arrays, without container overhead) and `live_values` (values alive anywhere, cache entries
  included). Each field is atomic; together they are not one consistent snapshot. An invalid
  query (foreign graph, bad root) counts nothing. A frozen cache's statistics follow the public
  call: its hits count only when the call succeeds, its Python result included; a failed call
  changes no statistic, except one miss for a call that raises `DagCacheMiss`.
- **Workers.** A supplied building cache is accepted only with `workers == 1` (`ValueError`
  otherwise). A frozen cache may be shared by every worker. Without a supplied cache, each
  worker uses its own private building cache.
- **Uids and frozen caches.** A per-pair Eligibility that the evaluator must intersect with a
  common one gets a new uid when the result differs from both, and that key cannot be in a
  frozen cache. To reuse prefilled entries, pass already-combined per-pair values with
  `elig=None`, or a common value that contains them (so the intersection returns the
  prefilled operand).
- **Scope.** Use one cache per series of analyses on one graph, and drop it at the end of the
  series.

```python
g = graph(5, [(0, 1, 10, 1), (1, 3, 10, 1), (0, 2, 2, 1), (2, 3, 2, 1), (3, 4, 9, 1)])
cache = ngc.DagCache(g)
opts = ngc.QueryOptions(config(ngc.CostOnly()), cache=cache)
rs = ngc.batch_capacity(g, [(0, 3), (1, 3)], opts)          # building: workers defaults to 1
print([r.placed for r in rs], cache.stats().hits, cache.stats().misses)
cache.freeze()
rs = ngc.batch_capacity(g, [(0, 3), (1, 3)], opts, workers=2)  # frozen: shared by 2 workers
print([r.placed for r in rs], cache.stats().hits, cache.stats().misses)
try:
    ngc.capacity(g, 0, 4, opts)                # destination 4: not in the cache
except ngc.DagCacheMiss:
    print("DagCacheMiss; misses =", cache.stats().misses)
```

```text
[5, 10] 1 1
[5, 10] 3 1
DagCacheMiss; misses = 2
```

## 10. Evaluators

**In plain words.** The evaluators answer questions about a model without touching your ledgers.
`capacity` asks how much one model can carry between two node sets by placing an unlimited
request on a fresh, empty ledger. `batch_capacity` asks that for many pairs. `sensitivity` asks,
for each link, how much the answer drops when that link is removed. A fresh ledger means these
answers can differ from placing into a ledger that other demands already occupy.

**Exact rules.**

- **`QueryOptions(config, *, elig=None, cache=None, outputs=None)`** is the one options type of
  the three evaluators: the configuration, the static Eligibility (`None`: everything allowed),
  the DagCache (`None`: a private cache per worker) and the optional outputs of each result
  (`sensitivity`: of its baseline).
- **`capacity(graph, sources, destinations, opts)`** is `place(float('inf'))` of a `FlowPolicy`
  with `opts.config` on a fresh ledger of the graph, and returns its `PlacementResult`. The exact maximum flow is the configuration
  `Global()` (its default origin is `OPTIMIZED`): its result carries the marginal cost curve and,
  with `OutputOptions(min_cut=True)`, the minimum cut. With `EQUAL` origin it is the largest total
  whose equal shares all fit.
- **`batch_capacity(graph, pairs, opts, per_pair_elig=None, workers=None)`** runs `capacity` for
  each `(sources, destinations)` pair and returns the results in pair order. Pair `i` uses
  `opts.elig` intersected with `per_pair_elig[i]` (entries may be `None`; a given list has one
  entry per pair). Everything is validated before any worker starts: first the arguments (the
  worker count, the graph of `opts.elig` and of the cache, the cache rule, `Pinned`, which is
  rejected because its bundles fix the endpoints, the configuration, then the length of
  `per_pair_elig`), then each pair completely (its endpoint sets, their overlap, its
  Eligibility's graph, the pooled single-path rule of [section 7](#configuration) for its
  sources) in pair order, so the error of the lowest failing pair is the one raised.
  An empty batch returns `[]`.
- **`sensitivity(graph, sources, destinations, opts, workers=None)`** returns a `SensitivityResult`: `baseline` (the result of
  `capacity`, with the requested outputs) and `entries`, one `SensitivityEntry(edge, loss)` per
  candidate edge, in ascending EdgeId.
  - **Candidates.** With `Global` and `OPTIMIZED` origin: the edges that carry flow in the
    baseline (removing an edge without flow cannot lower a maximum flow). With every other
    configuration, `EQUAL` origin included: every allowed edge, whatever its capacity (an edge of
    capacity 0 can carry a route and keep a source in the roster, so deleting it can change the
    result either way). An edge without an entry has loss 0.
  - **Loss** is `baseline.placed` minus the `placed` of `capacity` with that edge also excluded,
    in units, signed and never rounded or thresholded. A negative loss means the model carries
    more without the edge (for example, equal splitting no longer sends a share down a thin
    path, or a source behind it no longer holds the equal origin back). `Global` with `OPTIMIZED`
    origin computes each deletion by repairing a copy of the baseline's solve, which gives exactly
    the value of a fresh solve; with `EQUAL` origin each deletion is solved afresh. `CostOnly`
    gives loss 0 without recomputation to an edge outside the static shortest-path support
    reachable from the roster sources (no route can use it); these are equivalent forms of the
    same definition.
  - **Cache.** Each deletion uses a new Eligibility uid, so a `CostOnly` deletion that has to
    be placed raises `DagCacheMiss` through a frozen cache. The other routing modes reject a
    supplied cache. Private caches are cleared between candidates.
- **Workers.** An explicit `workers` (an integer of at least 1) is always honoured; the
  evaluator uses at most one thread per pair or candidate. With `workers=None` the default is 1
  when a cache is supplied, and otherwise comes from the environment:

  | Variable | Default worker count of | Accepted values |
  | --- | --- | --- |
  | `NGRAPH_CORE_BATCH_THREADS` | `batch_capacity` | an integer in `[1, 4096]`; unset or anything else means 1 |
  | `NGRAPH_CORE_SENSITIVITY_THREADS` | `sensitivity` | the same |

  The variables only supply defaults; they never override an explicit `workers`. The C++
  functions take the worker count as an explicit argument.
- **Failures** follow [section 8](#8-failures-threads-and-lifetimes): a failed evaluator call
  raises (also a solve beyond the flow core's size limit of [section 2](#2-the-graph), with
  `Global` or `MAX_FLOW`), and a failed deletion fails the whole `sensitivity` call (the
  exception of the lowest failing candidate; no entry is returned); a supplied building cache
  keeps the entries of the lookups that succeeded before the failure; a supplied frozen cache's
  statistics are unchanged, except one miss when the call raises `DagCacheMiss`.

```python
g = graph(4, [(0, 1, 10, 1), (1, 3, 10, 1), (0, 2, 2, 1), (2, 3, 2, 1)])
for name, cfg in (("CostOnly", config(ngc.CostOnly())),
                  ("Global", ngc.FlowPolicyConfig(routing=ngc.Global()))):
    s = ngc.sensitivity(g, 0, 3, ngc.QueryOptions(cfg))
    print(name, s.baseline.placed, [(e.edge, e.loss) for e in s.entries])
```

```text
CostOnly 5 [(0, 3), (1, -5), (2, 3), (3, -5)]
Global 12 [(0, 10), (1, 2), (2, 10), (3, 2)]
```

Under equal splitting the thin path (edges 1 and 3) holds the whole diamond back to 5 units;
without it the wide path alone carries 10, a loss of -5.

## 11. C++ names

The C++ library (namespace `netgraph::core`) has the same semantics. Its public headers are
self-contained:

| Header | Declares |
| --- | --- |
| `netgraph/core/types.hpp` | `NodeId`, `EdgeId` (`int32`), the cost and amount aliases (`int64`), `kDomainLimit` (`2^62`), `kUnreachable` |
| `netgraph/core/graph.hpp` | `Graph`, `to_units`, `as_float` |
| `netgraph/core/eligibility.hpp` | `Eligibility` |
| `netgraph/core/pred_dag.hpp` | `PredDAG`, `Path` |
| `netgraph/core/options.hpp` | `Direction`, `Absolute`, `Stretch`, `EdgeTieBreak`, `EdgePriority`, `PathSelection`, `FlowPlacement`, `Admission`, `Origin` |
| `netgraph/core/search.hpp` | `spf`, `ksp`, `resolve_to_paths`, `SpfOptions`, `SpfResult`, `KspOptions`, `SpfQueue` |
| `netgraph/core/flow_record.hpp` | `FlowId`, `Counter128` (the 128-bit cumulative counter), `FlowRecord` |
| `netgraph/core/flow_ledger.hpp` | `FlowLedger` |
| `netgraph/core/flow_policy.hpp` | `FlowPolicy`, its configuration and result types (`SourceFlow` among them), `Termination`, `kMaxFlowsPerCall` |
| `netgraph/core/dag_cache.hpp` | `DagCache`, `DagCacheMiss` |
| `netgraph/core/evaluators.hpp` | `capacity`, `batch_capacity`, `sensitivity` and their option and result types |

Spelling differences: Python enum members and constants are upper case
(`EdgeTieBreak.DETERMINISTIC`, `MAX_FLOWS_PER_CALL`, `UNREACHABLE`), C++ enumerators are capitalised
(`EdgeTieBreak::Deterministic`) and C++ constants carry a `k` prefix (`kMaxFlowsPerCall`,
`kUnreachable`). Where Python returns `None` for an unreachable cost (`target_cost`, `route_cost`),
C++ returns `kUnreachable`; an infinite request is an empty `std::optional` volume, and its
`remaining` and `requested` are empty too. An optional field that Python reports as `None`
(`placement`, `selection` and `origin` of a configuration, `offered`, `unadmitted`, `total_cost`,
`cost_distribution`, the optional outputs) is an empty `std::optional` in C++; the epoch's amounts
are `Counter128` values. A Python `ValueError` is a `std::invalid_argument` in C++. `spf` takes a
trailing C++-only `SpfQueue` argument that chooses the search's internal queue; the results never
depend on it.
