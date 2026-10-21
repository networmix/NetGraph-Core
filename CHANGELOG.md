# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

A rebuilt Core (breaking): one numeric domain (integer units), one flow core, one placement
engine with four routing models, forwarding programs and three evaluators. It replaces the 0.11
API.

### Added

- **Several demands at once**: `apply_actions(ledger, actions, *, mode=BatchMode.ORDERED,
  outputs=None)` runs `Place(policy, offer)`, `Replace(policy, offer, *, design_demand=None)` and
  `Forward(policy, offer, *, black_holes=None)` actions of several policies of one ledger as one
  transaction: validated whole in two passes (arguments and immutable objects for every action,
  then mutable state under the ledger gate, each pass in list order), then computed, then published
  once, with one `PlacementResult` per action in input order. The actions are immutable values
  that keep their policy and black holes alive and their offers as given. `BatchMode.ORDERED` runs
  each action as its scalar call on the state the earlier ones left (earlier actions are served
  first; a hold protects nothing); new FlowIds follow list order. `BatchMode.SIMULTANEOUS` takes
  `Replace` of `ShortestPaths` `EQUAL` or `WEIGHTED` and of `Pinned` and `Forward` of any installed
  program, every one LOSSY with a finite offer, and forwards them at the same time on one view (the
  residual plus every participant's old traffic; outside owners, orphans and holds stay): a full
  edge is shared by highest averages over the demands (ties by policy creation order, never list
  position), then within a demand over its heads, then its entries; a pooled demand's router
  attributes each member's amount to its sources; each participant's black holes are its own;
  loss cycles across demands are recalculated in rounds as within one program. Outputs read the
  state after each action (ORDERED) or after the whole list (SIMULTANEOUS). The forwarding size
  limits and the record identity limit count the whole list. C++: `apply_actions(FlowLedger&,
  std::span<const Action>, BatchMode, const OutputOptions&)` with `Action = std::variant<Place,
  Replace, Forward>`.

### Changed

- **Performance tooling**: A/B/A harness builds use 64-byte function alignment on both the
  reference and candidate to reduce timing variation due to code placement. Wheel build flags are unchanged.
- **Exact sensitivity**: unrestricted optimized `CentralizedTE` with a maximum baseline uses
  residual-component certificates, or stops a deletion repair when its
  scalar loss is known. The complete baseline, candidate order and reductions are preserved;
  every actual repair resets independently and publishes no intermediate allocation. Binding
  finite offers and other model domains retain their existing computation.
- **K shortest paths**: alternatives retain their generating spur index and skip searches
  in earlier prefix regions (Lawler’s partition), preserving complete path order, exact tie
  keys, masks and cost ceilings.
- **Concurrency**: native work releases the GIL. One gate per ledger protects reads and whole
  scalar/list publications; independent ledgers run in parallel. Policy construction and evaluators
  take one optional building-exclusive or frozen-shared cache lease, with FIFO lifecycle admission;
  workers borrow that lease. Native first-pass checks release the GIL before any lease or gate:
  the tunnel-limit roster search in policy construction and the evaluators (once per pair in a
  batch), and `from_routes`' walk of each supplied PredecessorDAG. `evaluate_pairs`' workers and
  the other evaluators' main evaluation reuse the checked rosters; each deletion of
  `evaluate_edge_failures` searches under its own mask. Immutable reads remain lock-free. cp314t detaches thread state
  before native waits. Publishing calls prepare Python results before native work; reads wrap
  owned snapshots after unlocking. `FlowPolicy` is final in Python, and failed construction
  leaves an unconstructed object that raises a clean `TypeError`.
- **Validation order**: argument and immutable checks run first, then mutable-state checks under
  the ledger gate or cache lease. Invalid offers, black holes and later list actions now precede
  missing-program/epoch errors; evaluator argument errors precede the building-cache worker rule.
  With a frozen empty cache, an invalid offer now precedes a cache miss in all three evaluators,
  and no miss is counted.
- **Cache clear**: its frozen-cache refusal now says to thaw it, then clear it; thaw waits for readers itself.
- **C++ signatures**: `FlowPolicy::program()`, `ShortestPathCache::freeze()` and `thaw()` can fail
  while acquiring a native lock and are no longer `noexcept`. Policy destruction remains
  `noexcept` and allocates nothing. A moved-from ledger retains its immutable graph owner so
  existing graph references remain valid; its public calls still reject the unusable ledger.

- **Graph**: `Graph(num_nodes, edge_sources, edge_destinations, edge_capacity_units, edge_costs,
  resolution)` (replacing 0.11's `StrictMultiDiGraph`) is the one constructor: capacities are
  `int64` units of an explicit resolution (a power of two in `[2^-32, 2^32]`), costs are positive
  integers, and the public `EdgeId` is the rank of an edge in `(source, cost, destination, input)`
  with `input_positions()` mapping ids to input positions. Validation reports the first defect in
  a fixed order (the resolution, the unit values, the shapes, zero or negative costs, the remaining
  bounds including the cost-sum guard), on the exact integers also beyond int64. Accessors
  (`edge_capacity_units()`, `edge_costs()`, `edge_sources()`, `edge_destinations()`,
  `outgoing_offsets()`, `incoming_offsets()`, `incoming_edges()`) return read-only NumPy views
  that keep their graph alive, also when empty; `max_edge_cost()` and `edge_cost_gcd()` read the
  costs.
  Construction numbers the edges by counting sort, in expected time linear in edges plus nodes
  for typical costs; above 65,536 distinct costs, or when the hash table that ranks the costs
  examines too many slots, it sorts all costs instead. Moving a graph no longer allocates.
  `Eligibility` is built from excluded node ids and EdgeIds, with `Eligibility.intersect(a, b)`.
  `to_units` (round half up) and `as_float` are the only conversions between physical values and
  units.
- **Route structures**: `PredecessorDAG.from_csr(graph, parent_offsets, parents, via_edges)` with
  the complete validation (shape, ranges, endpoints, duplicates, cycles, ascending via within a
  node); a DAG is stored over its member nodes and presents the CSR shape on demand.
  `Path(graph, hops)` with `source()`, `destination()`, `cost()` and `to_dag()`.
- **Route choice**: one tie rule, `RouteTieBreak` (`FEWER_HOPS`, the default, `HIGHER_CAPACITY`,
  `HIGHER_AVAILABLE_CAPACITY`, `HIGHER_HEADROOM`): cost first, then for a `HIGHER_*` value the
  larger bottleneck of its edge score (static capacity, the availability of the search, or their
  exact ratio), then fewer hops, then the smaller EdgeId sequence. `multipath: bool` (every
  equal-cost route, or the canonical path; `ShortestPaths` and `spf` only) and `parallel_edges:
  bool` (every parallel link, or one per pair) sit on the model or the search options. A
  `HIGHER_*` rule with every route and every parallel edge kept has nothing to choose and is
  rejected (`spf tie_break requires multipath=False or parallel_edges=False`, `ShortestPaths
  tie_break requires multipath=False or parallel_edges=False`). The tie rule never supplies
  weights.
- **Ledger**: `FlowLedger(graph, service_baseline=None, *, reservation_baseline=None)`
  (replacing 0.11's `FlowGraph`) is the integer ledger, in two projections that are never added:
  the carried traffic (`carried()`, `residual()`) and the strict holds of T1 tunnels
  (`reserved()`, `unreserved()`), each baseline validated within `[0, edge_capacity_units]`.
  Records as snapshots in FlowId order (`records()`, `record(flow_id)` with `FlowId = (ledger uid,
  id)`, ids never reused), `edge_capacity_units()`, `saturated_boundary_edges(sources,
  eligibility=None)` (the carried traffic only), `reset(service_baseline=None, *,
  reservation_baseline=None)` and `copy()` (both rejected while a policy is bound; `reset()` clears
  the records and holds a destroyed policy left behind; a copy keeps the holds). `FlowRecord`
  carries `flow_id`, `policy_id`, the endpoints as sets, `admission`, `entry_id` (the program entry
  a T1, T2 or Pinned record carries), `route` (a `PredecessorDAG`), `min_route_cost`, `offered`
  (Lossy only), `delivered`, `edge_flows`, `edge_drops` and `edge_reservations`; counters are
  128-bit (exact Python ints). Policies place, hold and remove; holds never police traffic.
- **Route allocations**: `RouteAllocation` (`EQUAL`, `WEIGHTED`, `OPTIMIZED`) and `Admission`
  (`LOSSLESS`, `LOSSY`). `EQUAL` splits equally by round robin from the lowest EdgeId; `WEIGHTED`
  splits by highest averages (D'Hondt, EdgeId ties) over integer weights frozen when the program
  is prepared, from a `WeightBasis` (`CAPACITY`, `EXPLICIT` with `edge_weights`,
  `RESIDUAL_SNAPSHOT`); `OPTIMIZED` re-solves the flow with the flow core on its route over the
  free capacity plus its own allocation; Lossy splits each offer per call and records the drops on
  each edge and what is lost at a source.
- **Configuration**: `FlowPolicyConfig(*, routing, admission=LOSSLESS)` holds the routing model
  and the admission; every other option belongs to the model that reads it:
  `ShortestPaths(*, allocation=EQUAL, source_allocation=None, weight_basis=None,
  edge_weights=None, tie_break=FEWER_HOPS, multipath=True, parallel_edges=True)` (ECMP by
  default; `WEIGHTED` is WCMP; `OPTIMIZED` a controller restricted to the shortest routes),
  `DistributedTE(tunnel_count, *, cost_ceiling=None, min_unreserved_bandwidth=0,
  tie_break=FEWER_HOPS, parallel_edges=True)` (T1), `CentralizedTE(*, source_allocation,
  support=None, tunnel_count=None, prefer_headroom=False)` (the source allocation is required;
  with `tunnel_count`, T2) and `Pinned(program, *, black_holes=None)`; `TunnelCount(per_head, *,
  up_to=False)`. `SourceAllocation` (`EQUAL`, `OPTIMIZED`) says how a demand's total enters at its
  sources. Every configuration value is immutable (read-only fields; a model's `edge_weights` is
  one shared array) and validated in two stages: construction raises `TypeError` for an unknown
  keyword (every field of another model included), a missing required argument, `None` where none
  is allowed, a bool field given anything but a `bool` or `numpy.bool_`, and an integer field
  given a non-integer (`TunnelCount per_head must be an integer (got bool)`); a policy or an
  evaluator raises `ValueError` for the rest, in a stated order: the endpoints, the ownership of
  the mask, program, support, cache and black holes, the enumeration values, the model matrix
  (route and source allocation, admission, the route-choice rules including the pooled
  single-path rule, `prefer_headroom requires tunnel_count`, the weight combinations, and last a
  cache on a model that does not read it), then the models' own fields, the tunnel limit and
  Pinned's endpoints. A policy's `config` is the configuration it was built with,
  `ShortestPaths`' absent source allocation and weight basis resolved.
- **Policies**: `FlowPolicy(ledger, sources, destinations, config, eligibility=None, cache=None)`
  is the one placement engine. `place(offer, outputs=None)` forwards `ShortestPaths`' pooled
  program or re-solves the controllers (`CentralizedTE` and `ShortestPaths` `OPTIMIZED`: the
  cumulative minimum-cost flow on the residual plus the policy's own allocation, under either
  source allocation); T1 and T2 have no `place`. `evaluate(offer, *, design_demand=None)` computes
  a whole new offer on the policy's released view and publishes nothing, `replace(offer, *,
  design_demand=None)` publishes it (new records, a new epoch; `design_demand` is the tunnel
  models' design, `None` meaning the offer), `forward(offer, *, black_holes=None)` sends a whole
  new offer over the installed program (no route computation; T1's holds and T2's plan kept), and
  `remove()`, `records()`, `record(flow_id)`, `delivered()`, `config` and `program` complete the
  policy. Every call returns `PlacementResult`: this call's `call_delivered`,
  `call_undelivered`, `delivery_delta`, `call_edge_drops`, `call_source_drops`, `termination`
  (`COMPLETE`, `NO_ROUTE`, `EXHAUSTED`, `NO_ACTIVE_WEIGHT`) and the flags `routes_changed`,
  `weights_changed` and `reservations_changed`; the demand's `requested`, `offered`, `admitted`,
  `delivered`, `unadmitted`, `source_drops`, `edge_drops`, `by_source` (`SourceFlow` rows),
  `by_tunnel` (`TunnelFlow` rows), `total_cost`, `marginal_cost_curve`, `program` and `plan` (the
  tunnel models' `TunnelPlan`); and `edge_flows`, `residual`, `saturated_boundary_edges`,
  `edge_reservations` and `min_cut` when requested through `OutputOptions`. An offer of 0 is
  validated like any other.
- **Pooled sources**: a demand with several sources is one demand. `ShortestPaths` with `EQUAL`
  or `WEIGHTED` prepares one program per epoch at its first positive call (each roster source's
  selected route on the residual at that call's entry; their union), splits the traffic of all
  sources once at every router they share, and gives each source its part by highest averages
  (records sum exactly to the totals); Lossless admits the largest equal total that fits with no
  drop, so a source that cannot pass anything holds the demand at the number of sources before it.
  `ShortestPaths` `OPTIMIZED` is one controller on the union of the sources' routes, as
  `CentralizedTE` is on the whole graph. With several sources, `multipath=False` requires
  `tie_break=FEWER_HOPS` (`ShortestPaths with multiple input sources requires tie_break=FEWER_HOPS
  when multipath=False`).
- **T1 tunnels**: `DistributedTE(TunnelCount(N))` models distributed RSVP-TE head ends, with no
  controller. A design volume D is split by round robin over the roster heads and each head's
  share over its N tunnels; in round order (round `j`, heads ascending: entry identifier
  `j * H + i`) each tunnel runs CSPF on the unreserved capacity (the ledger's, this policy's own
  holds released, minus this setup's earlier tunnels): the allowed links with at least
  `max(request, min_unreserved_bandwidth)` units, then the member rule with
  `parallel_edges=False`, then the canonical single path on static costs with the tie scores read
  on that view, then the cost ceiling (`AbsoluteCost` inclusive; `CostStretch` relative to the
  smallest static shortest cost over all roster heads, one pool baseline). A tunnel is `UP` and
  holds its request on every link of its path, or down (`INSUFFICIENT_BANDWIDTH`,
  `COST_CEILING_EXCEEDED`) with no path and weight 0. Traffic enters by round robin over the heads
  and then over each head's up tunnels; a head with no up tunnel loses its share.
- **Exact cost stretch**: `CostStretch(factor)` is one exact reduced fraction: a `Fraction` as it
  is, an integer as `n/1`, a float (Python or NumPy) as `Fraction(str(f))`, the decimal it prints,
  so `CostStretch(1.2)` is exactly 6/5; a bool or another type raises `TypeError`, a NaN or an
  infinity `ValueError`. `factor` returns the `Fraction`; equality, hash and `repr` follow it. A
  ceiling compares `cost * d <= n * c0` exactly in 128 bits. A factor below 1 is refused where the
  ceiling is used (`cost_ceiling CostStretch factor must be >= 1 (got <n>/<d>)`); a ceiling of
  another type is `TypeError("cost_ceiling must be AbsoluteCost, CostStretch or None")` in
  `DistributedTE` and `ksp` alike.
- **T2 tunnels**: `CentralizedTE(source_allocation=..., tunnel_count=TunnelCount(N))` models a
  central controller: every live head (one reaching a destination over links with free capacity
  of at least one unit; the others are dead, keep their roster position and lose their share) gets
  N single-path tunnels with one common whole-unit planned bandwidth. The largest common size U
  comes from quota flows, with the exact fluid optimum as a `Fraction`; a finite design lowers it
  to the largest live design share per tunnel; then, among the count vectors at that size, the
  smallest route-cost sum, then with `prefer_headroom=True` the smallest headroom penalty, then one
  canonical unit-by-unit augmentation; tunnels in round order with identifier `j * H + i` (dead
  heads leave holes), weight 1, `planned_bandwidth` u, never held. `TunnelCount(L, up_to=True)`
  chooses the smallest common count up to L whose program carries the live part of the design.
  No live head, or no whole unit per tunnel (`insufficient_resolution`), installs a program
  without tunnels that keeps the roster. `OPTIMIZED` source allocation only with one input source.
  With a `support`, T2 uses only the support's edges allowed by the policy's mask, in every step,
  and gives what the same configuration gives with those edges composed into the mask.
- **Forwarding programs**: `ForwardingProgram`, the routes and integer weights of one demand
  (tunnels `Tunnel`, split routes) with its injection (equal sources, or the captured source
  weights of an optimized controller). Every computing model produces one (`policy.program`,
  `result.program`); `ForwardingProgram.from_routes(graph, sources, destinations, routes, *,
  allocation, weight_basis, edge_weights, residual)` builds one from `RouteDefinition(head, route,
  *, weight=1)` (order and duplicates kept, every split route its own split context), and
  `ForwardingProgram.from_dag(graph, sources, destinations, dag, *, allocation, weight_basis,
  edge_weights, residual)` builds one shared hop-by-hop program on a supplied DAG, its sources
  pooled exactly as `ShortestPaths` pools them (one builder serves both). `Pinned(program)`
  forwards one: replayed under the policy's Eligibility (a tunnel touching an excluded element is
  `FAILED`, its head splits over the survivors) and, with `black_holes`, stale (the program
  unchanged, traffic reaching a black hole dropped there, `BLACKHOLED`). Lossy forwarding divides a
  shared edge by nested highest averages (demands, heads, then entries; what reaches a black hole
  takes no part), evaluates the program as a whole,
  and recalculates edges that wait on each other in a cycle in rounds (a period-2 or round-limit
  `ValueError`; an independent check before publication).
- **Failure guarantee**: every public call either completes or raises and changes nothing, also
  when memory runs out, a worker thread fails or the Python result cannot be built. Two limits:
  running out of memory inside pybind11's own object construction (outside the library's checked
  result construction) may crash the process; a `ShortestPathCache` passed to the call: a building
  cache keeps the complete entries and statistics of the lookups that succeeded before the
  failure, and a frozen cache keeps its entries and changes none of its statistics, except one
  miss when the call raised `ShortestPathCacheMiss`.
- **Search**: `spf(graph, roots, SpfOptions(direction, *, targets=None, tie_break=FEWER_HOPS,
  multipath=True, parallel_edges=True, residual=None, min_residual=1, eligibility=None))` returns
  `SpfResult`. It searches from a root set (every root at distance 0), `Direction.FORWARD` or
  `REVERSE`. Without targets it returns `costs` (`UNREACHABLE` for the nodes not reached) and the
  full DAG of tight edges; with targets it stops at the nearest targets' cost and returns
  `target_cost` and the DAG of the walks between the roots and a nearest target. With
  `multipath=False` it returns the canonical single path. A residual keeps only the edges with at
  least `min_residual` residual units. Every search runs on one traversal core, with a bucket queue
  when the largest edge cost is at most 65,535 times the costs' greatest common divisor, and a
  binary heap otherwise.
- **KSP and paths**: `ksp(graph, source, destination, KspOptions(k, *, tie_break=FEWER_HOPS,
  parallel_edges=True, cost_ceiling=None, residual=None, eligibility=None))` returns the K shortest
  simple paths as `(cost, Path)` pairs in the single-path order, truncated at `k` and at an
  inclusive ceiling. `resolve_to_paths(graph, dag, source, destination, split_parallel_edges=False,
  max_paths=None)` lists a DAG's walks in (cost, hops, edge ids) order.
- **ShortestPathCache**: `ShortestPathCache(graph)` is keyed by the root set, the direction and
  the Eligibility. Its values are the full static shortest-path DAGs (immutable, shared, valid
  after `clear()`), and `stats()` reports `hits`, `misses`, `entries`, `bytes` and `live_values`.
  - **Building mode:** one writer; a miss searches and inserts.
  - **Frozen mode:** `freeze()` makes the cache read-only and shareable by every worker of an
    evaluator. A frozen lookup miss raises `ShortestPathCacheMiss`, an exception the evaluators
    pass on unchanged. `clear()` raises when frozen, and `thaw()` returns to building mode after
    waiting for active calls. Lifecycle admission is FIFO, with no overtaking by later readers.
- **Evaluators**: one options type, `EvaluationOptions(config, *, eligibility=None, cache=None,
  outputs=None, design_demand=None, offer=None)`, the design and the offer kept as given and
  checked by the evaluator; `evaluate_demand(graph, sources, destinations, options)` (with
  `CentralizedTE(source_allocation=OPTIMIZED)` the exact maximum flow, its marginal cost curve and
  its minimum cut), `evaluate_pairs(graph, pairs, options, per_pair_eligibility=None,
  workers=None)` (pair order; the lowest failing index raises; Pinned is rejected) and
  `evaluate_edge_failures(graph, sources, destinations, options, workers=None)` (every candidate
  edge and its signed `delivery_reduction`; a failed deletion fails the call). A finite `offer` is
  accepted for every model under each admission it supports, in all three evaluators, and every
  pair, baseline and failed network answers that same offer; an absent offer is the uncapped
  Lossless question or the first-hop experiment of `ShortestPaths` Lossy. T1 always names a finite
  design; T2 takes the offer as its design when none is given. A supplied building cache needs
  `workers == 1`; a frozen cache is shared by every worker. With `workers=None` the count is 1 when
  a cache is supplied, otherwise `NGRAPH_CORE_BATCH_THREADS` and
  `NGRAPH_CORE_SENSITIVITY_THREADS` (an integer in [1, 4096]; 1 when unset or invalid).
- **Flow core**: one exact integer minimum-cost flow solver (successive shortest paths with
  potentials, then Dinic phases for each cost level; each edge is its own arc, visited in
  ascending EdgeId) serves `CentralizedTE`, the `OPTIMIZED` allocation, T2's sizing (every cost
  taken as 0: a maximum flow) and the evaluators (`evaluate_edge_failures` of an unrestricted
  optimized `CentralizedTE` deletes only the edges its baseline loads, and when that baseline is a
  maximum flow, with no offer or an offer above it, repairs the baseline's solution for each removed
  edge instead of solving again; at an offer that binds it solves each deletion at the offer); the
  private entry
  `_netgraph_core._solve_min_cost_flow` returns the value, the marginal cost curve, the edge flows,
  the minimum cut, the potentials and six operation counters. The solver numbers its network in
  `int32`, so a solve raises `ValueError`, before it sizes anything, on a network of more than
  `INT32_MAX - 2` nodes or with more than `2^30 - 1` edges, sources and destinations together.
- **Zero-capacity links stay in routing**: `spf`, `ksp`, the cached supports, the `ShortestPaths`
  routes, every roster and the stretch baselines keep edges below one unit of capacity; only a
  supplied residual (or an `Eligibility`) filters. Without a residual, `HIGHER_HEADROOM` reads 1 on
  an edge of positive capacity and 0 on an edge of capacity 0. `EQUAL` keeps a zero-capacity
  member and its round-robin share; `WEIGHTED` gives it weight 0 and drops a branch left without a
  member, without routing around it; `OPTIMIZED` keeps it and never bypasses it.
- **C++ headers** (under `netgraph/core/`): `types.hpp`, `graph.hpp` (`Graph`, `to_units`,
  `as_float`; it includes `eligibility.hpp`), `eligibility.hpp`, `options.hpp` (`Direction`,
  `AbsoluteCost`, `CostStretch`, `RouteTieBreak`, `RouteAllocation`, `WeightBasis`, `Admission`,
  `SourceAllocation`), `route_structures.hpp`, `forwarding_program.hpp` (`ForwardingProgram`,
  `Tunnel`, `TunnelStatus`, `RouteDefinition`, `RouteOptions`), `flow_record.hpp` (`FlowId`,
  `Counter128`, `FlowRecord`), `flow_ledger.hpp`, `flow_policy.hpp`, `shortest_path_cache.hpp`,
  `evaluators.hpp` and `search.hpp` (`spf`, `ksp`, `resolve_to_paths` and their options; the
  queue choice `SpfQueue { Auto, Heap }` is a C++-only trailing parameter of `spf`); each
  compiles alone and includes no `detail/` header.
- **Typing**: the package ships `__init__.pyi`, a type stub generated from the compiled module
  (pybind11-stubgen with reviewed overrides: exact integer and array types, NumPy result dtypes,
  optional values, and no constructor for the result classes the module refuses to construct),
  next to `py.typed`. The stub carries the binding docstrings. `bash dev/stubs/run.sh generate`
  regenerates it; `make stubs-check` checks it against a fresh generation, the wheel and sdist
  contents, and typed client programs.
- **Input error classes**: one rule at every public entry. A type the position never takes raises
  `TypeError`: a bool (Python's or `numpy.bool_`) at any numeric position, a float where only an
  integer is taken, a string, `None` where none is allowed, an object of another class; for a vector
  a non-sequence, a nested entry, an array of another rank or of a non-integer dtype. A value
  outside its domain raises `ValueError`: a range, a length, NaN, `-inf`, a finite float as an offer
  or a design demand. Texts are unchanged. Against this release's earlier development snapshots,
  the wrong-type errors of call arguments, search options, graph arrays and node ids, vectors,
  amounts, `None` for `graph`, `ledger` or `dag`, FlowIds, `evaluate_pairs`' pairs,
  `Eligibility.intersect`, `Path` hops and `to_units` become `TypeError` (were `ValueError`):
  `place(True)` and `SpfOptions(FORWARD, min_residual=1.0)` at the call raise `TypeError`, while
  `place(1.5)` and `min_residual=0` stay `ValueError`. The non-integer number positions (a
  `resolution` of `Graph`, `to_units` and `as_float`, `to_units`' values) take exactly ints, floats
  and NumPy numbers, a subclass by its value: a bool (a resolution read as 1.0 or 0.0 before), a
  `Decimal`, a `Fraction`, an object known by `__float__`, `__index__` or `__array__`, or a 0-d array
  as a resolution is refused with `TypeError`; an integer beyond the float range as a resolution is now a `ValueError` with the
  resolution's text (was pybind11's `TypeError`). `to_units` refuses a bool entry beside numbers
  (`[True, 1.0]`, accepted before) with `TypeError`. Integer positions, `as_float`'s amounts
  included, keep Python's integer protocol (`__index__`, never a bool).
- **Results that change** against this release's earlier development snapshots, on purpose: a
  zero-capacity link stays in routing, so a direct edge of capacity 0 and cost 1 beside a cost-2
  route makes the shortest cost 1 (was 2), and an equal LOSSLESS split that sends that edge its
  share delivers 0 (was 5); pooled sources whose routes meet at a router split there jointly, so
  two sources reaching two unit-capacity members deliver 2 with LOSSY (was 1); `CostStretch(1.2)`
  is the exact fraction 6/5, so a baseline cost of 5 admits a path of cost 6 at equality and the
  tunnel case that delivered 0 now delivers 1; a black-holed edge no longer waits in a loss cycle,
  so a cycle that never settles raises its round-limit error after 64 rounds per edge that does
  wait (a loaded black-holed edge inside such a cycle added 64 more), with every amount, every
  settling round and every period-2 error unchanged; traffic that can carry nothing (a tunnel,
  split-route member or pooled source whose share is 0, or traffic behind its own demand's black
  hole) no longer counts in the loads that decide which edges can drop, makes an edge wait or keeps
  two edges in series from passing traffic through as one, so adding a demand that can carry nothing
  no longer changes another demand's result, error or round count (it could turn a loss cycle's
  error into a result or a result into an error).

### Removed

- The 0.11 library and its bindings: `FlowState`, the previous `FlowGraph` and `FlowPolicy` and
  their configuration, the backend and algorithms facade (with its `Graph` handle), the search
  entry points (`shortest_paths`, `shortest_paths_to`, `spf_dag`, `spf_to`, `k_shortest_paths`
  with its `unique` flag, `PredDAG.resolve_to_paths`, `NGRAPH_CORE_SPF_QUEUE`), max-flow,
  batch max-flow and the previous sensitivity, the SPF memo and `StateStamp`, `kMinFlow`,
  `kEpsilon`, `constants.hpp`, and the profiling facility (the Python accessors, the C++ scopes
  and `NGRAPH_CORE_PROFILE`).
- The hand-written typing module `netgraph_core/_docs.py`, replaced by the generated stub.
- The growth model of the earlier capacity-aware routing (replaced by T1's tunnels): its fields
  `flow_count`, `max_flow_count`, `min_flow_bandwidth`, `reoptimize` and `diminishing`,
  `DiminishingReturns`, the per-call creation cap `MAX_FLOWS_PER_CALL` (C++ `kMaxFlowsPerCall`),
  the terminations `LIMIT_FLOWS` and `LIMIT_DIMINISHING` and the per-source reduction of outcomes.
- The sensitivity entry status (`SensitivityStatus`, `SensitivityEntry.status` and `.reason`):
  successful entries were always `EVALUATED`; a failed deletion fails the call.
- The names of this release's earlier development snapshots. Nothing old survives as an alias, a
  warning or an error that names a replacement; each name below is simply gone. Public names:

| Name | Replacement |
| --- | --- |
| `PredDAG` | `PredecessorDAG` |
| `EdgeTieBreak` | `RouteTieBreak` |
| `EdgePriority` | removed: the score is part of the `RouteTieBreak` value (`HIGHER_CAPACITY`, `HIGHER_AVAILABLE_CAPACITY`, `HIGHER_HEADROOM`) |
| `PathSelection` | removed: `tie_break`, `multipath`, `parallel_edges` on `ShortestPaths`, `DistributedTE`, `SpfOptions`, `KspOptions` |
| `SpfOptions` | same name; `selection` replaced by `tie_break`, `multipath`, `parallel_edges` |
| `KspOptions` | same name; `selection` replaced by `tie_break`, `parallel_edges` (no multipath) |
| `Absolute` | `AbsoluteCost` |
| `Stretch` | `CostStretch` |
| `DagCache` | `ShortestPathCache` |
| `DagCacheMiss` | `ShortestPathCacheMiss` |
| `DagCacheStats` | `ShortestPathCacheStats` |
| `FlowPlacement` | `RouteAllocation` |
| `Origin` | `SourceAllocation` |
| `Tunnels` | `TunnelCount` |
| `CostOnly` | `ShortestPaths` |
| `CapacityAware` | `DistributedTE` |
| `Global` | `CentralizedTE` |
| `FlowPolicyConfig` | same name; only `routing` and `admission` remain |
| `ForwardingProgram` | same name; new `from_dag` |
| `Route` | `RouteDefinition` |
| `QueryOptions` | `EvaluationOptions` |
| `CapacityOptions` | `EvaluationOptions` |
| `SensitivityOptions` | `EvaluationOptions` |
| `capacity` | `evaluate_demand` |
| `batch_capacity` | `evaluate_pairs` |
| `sensitivity` | `evaluate_edge_failures` |

  Members, parameters and enum values:

| Kind | Name | Replacement |
| --- | --- | --- |
| parameter | `batch_capacity:opts` | `evaluate_pairs:options` |
| parameter | `batch_capacity:per_pair_elig` | `evaluate_pairs:per_pair_eligibility` |
| parameter | `capacity:opts` | `evaluate_demand:options` |
| parameter | `CapacityAware.__init__:max_path_cost` | `DistributedTE.__init__:cost_ceiling` |
| parameter | `CapacityAware.__init__:tunnels` | `DistributedTE.__init__:tunnel_count` |
| property | `CapacityAware.max_path_cost` | `DistributedTE.cost_ceiling` |
| property | `CapacityAware.tunnels` | `DistributedTE.tunnel_count` |
| field | `EdgePriority.RESIDUAL` | removed; `RouteTieBreak.HIGHER_AVAILABLE_CAPACITY` |
| field | `EdgeTieBreak.DETERMINISTIC` | `RouteTieBreak.FEWER_HOPS` |
| field | `EdgeTieBreak.PREFER_HIGHER_PRIORITY` | removed; `RouteTieBreak.HIGHER_CAPACITY`, `HIGHER_AVAILABLE_CAPACITY` or `HIGHER_HEADROOM`, by the old priority |
| parameter | `FlowLedger.__init__:residual` | `FlowLedger.__init__:service_baseline` |
| parameter | `FlowLedger.__init__:unreserved` | `FlowLedger.__init__:reservation_baseline` |
| method | `FlowLedger.allocated` | `FlowLedger.carried` |
| method | `FlowLedger.bottleneck` | `FlowLedger.saturated_boundary_edges` |
| parameter | `FlowLedger.bottleneck:elig` | `FlowLedger.saturated_boundary_edges:eligibility` |
| method | `FlowLedger.capacity_units` | `FlowLedger.edge_capacity_units` |
| parameter | `FlowLedger.reset:residual` | `FlowLedger.reset:service_baseline` |
| parameter | `FlowLedger.reset:unreserved` | `FlowLedger.reset:reservation_baseline` |
| field | `FlowPlacement.EQUAL_BALANCED` | `RouteAllocation.EQUAL` |
| field | `FlowPlacement.MAX_FLOW` | `RouteAllocation.OPTIMIZED` |
| parameter | `FlowPolicy.__init__:elig` | `FlowPolicy.__init__:eligibility` |
| parameter | `FlowPolicy.place:volume` | `FlowPolicy.place:offer` |
| method | `FlowPolicy.placed_demand` | `FlowPolicy.delivered` |
| parameter | `FlowPolicyConfig.__init__:origin` | `ShortestPaths(source_allocation=...)`, `CentralizedTE(source_allocation=...)` (required there); none on `DistributedTE` or `Pinned` |
| parameter | `FlowPolicyConfig.__init__:placement` | `ShortestPaths(allocation=...)` (`from_routes` and `from_dag` have their own `allocation`) |
| parameter | `FlowPolicyConfig.__init__:weights` | `ShortestPaths(weight_basis=...)` |
| property | `FlowPolicyConfig.origin` | `ShortestPaths.source_allocation`, `CentralizedTE.source_allocation` |
| property | `FlowPolicyConfig.placement` | `ShortestPaths.allocation` |
| property | `FlowPolicyConfig.weights` | `ShortestPaths.weight_basis` |
| property | `FlowRecord.bundle` | `FlowRecord.route` |
| property | `FlowRecord.drops` | `FlowRecord.edge_drops` |
| property | `FlowRecord.edges` | `FlowRecord.edge_flows` |
| property | `FlowRecord.entry` | `FlowRecord.entry_id` |
| property | `FlowRecord.owner` | `FlowRecord.policy_id` |
| field | `FlowRecord.reserved_edges` | `FlowRecord.edge_reservations` |
| property | `FlowRecord.route_cost` | `FlowRecord.min_route_cost` |
| parameter | `ForwardingProgram.from_routes:placement` | `ForwardingProgram.from_routes:allocation` |
| parameter | `ForwardingProgram.from_routes:weights` | `ForwardingProgram.from_routes:weight_basis` |
| parameter | `Global.__init__:tunnels` | `CentralizedTE.__init__:tunnel_count` |
| property | `Global.tunnels` | `CentralizedTE.tunnel_count` |
| parameter | `Graph.__init__:capacity_units` | `Graph.__init__:edge_capacity_units` |
| parameter | `Graph.__init__:cost` | `Graph.__init__:edge_costs` |
| parameter | `Graph.__init__:dst` | `Graph.__init__:edge_destinations` |
| parameter | `Graph.__init__:src` | `Graph.__init__:edge_sources` |
| method | `Graph.capacity_units` | `Graph.edge_capacity_units` |
| method | `Graph.cost` | `Graph.edge_costs` |
| method | `Graph.cost_gcd` | `Graph.edge_cost_gcd` |
| method | `Graph.edge_dst` | `Graph.edge_destinations` |
| method | `Graph.edge_src` | `Graph.edge_sources` |
| method | `Graph.in_adj_edge_index` | `Graph.incoming_edges` |
| method | `Graph.in_row_offsets` | `Graph.incoming_offsets` |
| method | `Graph.input_index` | `Graph.input_positions` |
| method | `Graph.max_cost` | `Graph.max_edge_cost` |
| method | `Graph.row_offsets` | `Graph.outgoing_offsets` |
| parameter | `ksp:dst` | `ksp:destination` |
| parameter | `ksp:opts` | `ksp:options` |
| parameter | `ksp:src` | `ksp:source` |
| parameter | `KspOptions.__init__:elig` | `KspOptions.__init__:eligibility` |
| parameter | `KspOptions.__init__:max_path_cost` | `KspOptions.__init__:cost_ceiling` |
| property | `KspOptions.elig` | `KspOptions.eligibility` |
| property | `KspOptions.max_path_cost` | `KspOptions.cost_ceiling` |
| parameter | `OutputOptions.__init__:bottleneck` | `OutputOptions.__init__:saturated_boundary_edges` |
| parameter | `OutputOptions.__init__:reservation_edge_flows` | `OutputOptions.__init__:edge_reservations` |
| property | `OutputOptions.bottleneck` | `OutputOptions.saturated_boundary_edges` |
| property | `OutputOptions.reservation_edge_flows` | `OutputOptions.edge_reservations` |
| method | `Path.dst` | `Path.destination` |
| method | `Path.src` | `Path.source` |
| parameter | `PathSelection.__init__:multi_edge` | `parallel_edges=` on `ShortestPaths`, `DistributedTE`, `SpfOptions`, `KspOptions` |
| property | `PathSelection.multi_edge` | the model's or the search options' `parallel_edges` |
| property | `PlacementResult.bottleneck` | `PlacementResult.saturated_boundary_edges` |
| property | `PlacementResult.call_drops` | `PlacementResult.call_edge_drops` |
| property | `PlacementResult.call_source_loss` | `PlacementResult.call_source_drops` |
| property | `PlacementResult.cost_distribution` | `PlacementResult.marginal_cost_curve` |
| property | `PlacementResult.drops` | `PlacementResult.edge_drops` |
| property | `PlacementResult.placed` | `PlacementResult.call_delivered` |
| property | `PlacementResult.remaining` | `PlacementResult.call_undelivered` |
| field | `PlacementResult.reservation_edge_flows` | `PlacementResult.edge_reservations` |
| property | `PlacementResult.source_loss` | `PlacementResult.source_drops` |
| parameter | `PredDAG.from_csr:offsets` | `PredecessorDAG.from_csr:parent_offsets` |
| parameter | `PredDAG.from_csr:via` | `PredecessorDAG.from_csr:via_edges` |
| parameter | `QueryOptions.__init__:elig` | `EvaluationOptions.__init__:eligibility` |
| parameter | `QueryOptions.__init__:offered` | `EvaluationOptions.__init__:offer` |
| property | `QueryOptions.elig` | `EvaluationOptions.eligibility` |
| property | `QueryOptions.offered` | `EvaluationOptions.offer` |
| parameter | `resolve_to_paths:dst` | `resolve_to_paths:destination` |
| parameter | `resolve_to_paths:src` | `resolve_to_paths:source` |
| parameter | `Route.__init__:bundle` | `RouteDefinition.__init__:route` |
| property | `Route.bundle` | `RouteDefinition.route` |
| parameter | `sensitivity:opts` | `evaluate_edge_failures:options` |
| property | `SensitivityEntry.loss` | `SensitivityEntry.delivery_reduction` |
| property | `SourceFlow.dropped` | `SourceFlow.total_edge_drops` |
| property | `SourceFlow.source_loss` | `SourceFlow.source_drops` |
| parameter | `spf:opts` | `spf:options` |
| parameter | `SpfOptions.__init__:elig` | `SpfOptions.__init__:eligibility` |
| property | `SpfOptions.elig` | `SpfOptions.eligibility` |
| property | `Tunnel.index` | `Tunnel.entry_id` |
| field | `Tunnel.planned` | `Tunnel.planned_bandwidth` |
| field | `Tunnel.requested` | `Tunnel.requested_bandwidth` |
| field | `Tunnel.reserved` | `Tunnel.reserved_bandwidth` |
| property | `Tunnel.route_cost` | `Tunnel.path_cost` |
| property | `TunnelFlow.dropped` | `TunnelFlow.total_edge_drops` |
| property | `TunnelFlow.index` | `TunnelFlow.entry_id` |
| field | `TunnelPlan.common_share` | `TunnelPlan.max_bandwidth_per_tunnel` |
| field | `TunnelPlan.count` | `TunnelPlan.tunnels_per_head` |
| field | `TunnelPlan.fluid_share` | `TunnelPlan.max_fractional_bandwidth_per_tunnel` |
| field | `TunnelPlan.forwarding_bound` | `TunnelPlan.max_lossless_offer` |
| field | `TunnelPlan.head_design` | `TunnelPlan.design_by_head` |
| field | `TunnelPlan.live_heads` | `TunnelPlan.heads_with_capacity` |
| field | `TunnelPlan.planned_per_tunnel` | `TunnelPlan.planned_bandwidth_per_tunnel` |
| field | `TunnelPlan.planned_total` | `TunnelPlan.total_planned_bandwidth` |
| field | `TunnelPlan.reserved_total` | `TunnelPlan.total_reserved_bandwidth` |
| field | `TunnelPlan.route_cost` | `TunnelPlan.total_tunnel_cost` |
| parameter | `Tunnels.__init__:count` | `TunnelCount.__init__:per_head` |
| property | `Tunnels.count` | `TunnelCount.per_head` |
| field | `TunnelStatus.COST_CEILING` | `TunnelStatus.COST_CEILING_EXCEEDED` |
| field | `TunnelStatus.NO_BANDWIDTH` | `TunnelStatus.INSUFFICIENT_BANDWIDTH` |
| parameter | `FlowPolicyConfig.__init__:edge_weights` | `ShortestPaths(edge_weights=...)` |
| property | `FlowPolicyConfig.edge_weights` | `ShortestPaths.edge_weights` |
| parameter | `FlowPolicyConfig.__init__:selection` | removed; its fields became `tie_break`, `multipath`, `parallel_edges` of the model |
| property | `FlowPolicyConfig.selection` | removed |
| parameter | `PathSelection.__init__:tie_break`, `priority` | one `tie_break: RouteTieBreak` |
| parameter | `PathSelection.__init__:multipath` (bool or None) | `multipath: bool = True` on `ShortestPaths` and `SpfOptions` only |
| parameter | `SpfOptions.__init__:selection`, `KspOptions.__init__:selection` | removed; `tie_break`, `parallel_edges` (and `multipath` for `SpfOptions`) as keywords |
| parameter | `Global.__init__:support` + `tunnels` together (rejected before) | `CentralizedTE(source_allocation=..., support=dag, tunnel_count=...)` allowed |
| parameter | `Global.__init__` without an origin | `CentralizedTE(source_allocation=...)` required |
| property | `Stretch.factor` (float) | `CostStretch.factor` (a `fractions.Fraction`) |
| parameter | `Stretch.__init__:factor` float | a `Fraction`, an integer or a float read as `Fraction(str(f))` |
| parameter | `QueryOptions.__init__:offered` (some models only) | `EvaluationOptions(offer=...)`, finite, every model |
| property | `FlowRecord.placement` | removed: the configuration selects the computation, the program the forwarding |
| property | `PlacementResult.local_outcomes` | removed: `by_source` carries the per-source accounting (no per-source termination) |
| method | `FlowPolicy.usable_bundles` | removed: inspect the program (`FlowPolicy.program`) |
| property | `Pinned.bundles` | `Pinned.program`, built by `ForwardingProgram.from_routes` or `from_dag` |

  C++ declarations and headers (C++ enum values are the Python words in PascalCase):

| Kind | C++ name | Replacement |
| --- | --- | --- |
| file | `include/netgraph/core/dag_cache.hpp` | `include/netgraph/core/shortest_path_cache.hpp` |
| file | `include/netgraph/core/pred_dag.hpp` | `include/netgraph/core/route_structures.hpp` |
| file | `include/netgraph/core/program.hpp` | `include/netgraph/core/forwarding_program.hpp` |
| type | `Absolute` | `AbsoluteCost` |
| function | `batch_capacity` | `evaluate_pairs` |
| parameter | `batch_capacity:g` | `evaluate_pairs:graph` |
| function | `capacity` | `evaluate_demand` |
| parameter | `capacity:g` | `evaluate_demand:graph` |
| type | `CapacityAware` | `DistributedTE` |
| field | `CapacityAware::max_path_cost` | `DistributedTE::cost_ceiling` |
| type | `CostOnly` | `ShortestPaths` |
| type | `DagCache` | `ShortestPathCache` |
| constructor_or_destructor | `DagCache::DagCache` | `ShortestPathCache::ShortestPathCache` |
| parameter | `DagCache::lookup:elig` | `ShortestPathCache::lookup:eligibility` |
| parameter | `DagCache::lookup:g` | `ShortestPathCache::lookup:graph` |
| constructor_or_destructor | `DagCache::~DagCache` | `ShortestPathCache::~ShortestPathCache` |
| type | `DagCacheMiss` | `ShortestPathCacheMiss` |
| constructor_or_destructor | `DagCacheMiss::DagCacheMiss` | `ShortestPathCacheMiss::ShortestPathCacheMiss` |
| enum-value | `EdgePriority::Residual` | removed with `EdgePriority` |
| type | `EdgeTieBreak` | `RouteTieBreak` |
| enum-value | `EdgeTieBreak::Deterministic` | `RouteTieBreak::FewerHops` |
| enum-value | `EdgeTieBreak::PreferHigherPriority` | removed; `RouteTieBreak::HigherCapacity`, `HigherAvailableCapacity` or `HigherHeadroom` |
| function | `Eligibility::edge_ok` | `Eligibility::allows_edge` |
| function | `Eligibility::node_ok` | `Eligibility::allows_node` |
| function | `FlowLedger::allocated` | `FlowLedger::carried` |
| function | `FlowLedger::bottleneck` | `FlowLedger::saturated_boundary_edges` |
| parameter | `FlowLedger::bottleneck:elig` | `FlowLedger::saturated_boundary_edges:eligibility` |
| function | `FlowLedger::capacity_units` | `FlowLedger::edge_capacity_units` |
| parameter | `FlowLedger::FlowLedger:residual` | `FlowLedger::FlowLedger:service_baseline` |
| parameter | `FlowLedger::reset:residual` | `FlowLedger::reset:service_baseline` |
| type | `FlowPlacement` | `RouteAllocation` |
| enum-value | `FlowPlacement::EqualBalanced` | `RouteAllocation::Equal` |
| enum-value | `FlowPlacement::MaxFlow` | `RouteAllocation::Optimized` |
| parameter | `FlowPolicy::FlowPolicy:elig` | `FlowPolicy::FlowPolicy:eligibility` |
| parameter | `FlowPolicy::place:volume` | `FlowPolicy::place:offer` |
| function | `FlowPolicy::placed_demand` | `FlowPolicy::delivered` |
| field | `FlowPolicyConfig::origin` | `ShortestPaths::source_allocation`, `CentralizedTE::source_allocation` |
| field | `FlowPolicyConfig::placement` | `ShortestPaths::allocation` (the option lives in the model) |
| field | `FlowPolicyConfig::weights` | `ShortestPaths::weight_basis` |
| field | `FlowRecord::bundle` | `FlowRecord::route` |
| field | `FlowRecord::drops` | `FlowRecord::edge_drops` |
| field | `FlowRecord::edges` | `FlowRecord::edge_flows` |
| field | `FlowRecord::entry` | `FlowRecord::entry_id` |
| field | `FlowRecord::owner` | `FlowRecord::policy_id` |
| field | `FlowRecord::route_cost` | `FlowRecord::min_route_cost` |
| parameter | `ForwardingProgram::from_routes:g` | `ForwardingProgram::from_routes:graph` |
| type | `Global` | `CentralizedTE` |
| function | `Graph::capacity_units` | `Graph::edge_capacity_units` |
| function | `Graph::cost` | `Graph::edge_costs` |
| function | `Graph::cost_gcd` | `Graph::edge_cost_gcd` |
| function | `Graph::edge_dst` | `Graph::edge_destinations` |
| function | `Graph::edge_src` | `Graph::edge_sources` |
| parameter | `Graph::Graph:capacity_units` | `Graph::Graph:edge_capacity_units` |
| parameter | `Graph::Graph:cost` | `Graph::Graph:edge_costs` |
| parameter | `Graph::Graph:dst` | `Graph::Graph:edge_destinations` |
| parameter | `Graph::Graph:src` | `Graph::Graph:edge_sources` |
| function | `Graph::in_adj_edge_index` | `Graph::incoming_edges` |
| function | `Graph::in_row_offsets` | `Graph::incoming_offsets` |
| function | `Graph::input_index` | `Graph::input_positions` |
| function | `Graph::max_cost` | `Graph::max_edge_cost` |
| function | `Graph::row_offsets` | `Graph::outgoing_offsets` |
| parameter | `ksp:g` | `ksp:graph` |
| field | `KspOptions::elig` | `KspOptions::eligibility` |
| field | `KspOptions::max_path_cost` | `KspOptions::cost_ceiling` |
| type | `Origin` | `SourceAllocation` |
| field | `OutputOptions::bottleneck` | `OutputOptions::saturated_boundary_edges` |
| function | `Path::dst` | `Path::destination` |
| parameter | `Path::Path:g` | `Path::Path:graph` |
| function | `Path::src` | `Path::source` |
| field | `PathSelection::multi_edge` | `parallel_edges` of the model or the search options (`PathSelection` removed) |
| field | `PlacementResult::bottleneck` | `PlacementResult::saturated_boundary_edges` |
| field | `PlacementResult::call_drops` | `PlacementResult::call_edge_drops` |
| field | `PlacementResult::call_source_loss` | `PlacementResult::call_source_drops` |
| field | `PlacementResult::cost_distribution` | `PlacementResult::marginal_cost_curve` |
| field | `PlacementResult::drops` | `PlacementResult::edge_drops` |
| field | `PlacementResult::placed` | `PlacementResult::call_delivered` |
| field | `PlacementResult::remaining` | `PlacementResult::call_undelivered` |
| field | `PlacementResult::source_loss` | `PlacementResult::source_drops` |
| type | `PredDAG` | `PredecessorDAG` |
| parameter | `PredDAG::from_csr:g` | `PredecessorDAG::from_csr:graph` |
| type | `QueryOptions` | `EvaluationOptions` |
| field | `QueryOptions::elig` | `EvaluationOptions::eligibility` |
| field | `QueryOptions::offered` | `EvaluationOptions::offer` |
| parameter | `resolve_to_paths:g` | `resolve_to_paths:graph` |
| type | `Route` | `RouteDefinition` |
| field | `Route::bundle` | `RouteDefinition::route` |
| field | `RouteOptions::placement` | `RouteOptions::allocation` |
| field | `RouteOptions::weights` | `RouteOptions::weight_basis` |
| function | `sensitivity` | `evaluate_edge_failures` |
| parameter | `sensitivity:g` | `evaluate_edge_failures:graph` |
| field | `SensitivityEntry::loss` | `SensitivityEntry::delivery_reduction` |
| field | `SourceFlow::dropped` | `SourceFlow::total_edge_drops` |
| field | `SourceFlow::source_loss` | `SourceFlow::source_drops` |
| parameter | `spf:g` | `spf:graph` |
| field | `SpfOptions::elig` | `SpfOptions::eligibility` |
| type | `Stretch` | `CostStretch` |
| field | `Tunnel::index` | `Tunnel::entry_id` |
| field | `Tunnel::route_cost` | `Tunnel::path_cost` |
| field | `TunnelFlow::dropped` | `TunnelFlow::total_edge_drops` |
| field | `TunnelFlow::index` | `TunnelFlow::entry_id` |

### Fixed

- Fixed double frees and leaks when native allocation fails while constructing Python result
  wrappers for `spf`, `ksp` and `resolve_to_paths`; these failures now raise `MemoryError`.

## [0.11.0] - 2026-10-02

### Changed

- **Shortest Paths**: `shortest_paths` / `shortest_paths_to` (and therefore `spf`, `spf_to`, max-flow, `FlowPolicy` and KSP spur searches) now run Dijkstra on a monotone bucket queue whenever the graph's edge costs fit 65,536 buckets (`max_cost / gcd(costs) + 1`), which covers unit, uniform and small-integer cost graphs including zero-cost pseudo attachments; wider or int64-offset costs keep the binary heap. The active bucket is ordered by the same `(-bottleneck, node)` key as the heap, so the pop order and every output (distances, predecessor order, single-path capacity ties, early exit) are bit-identical; a differential test in the C++ suite enforces this. Scratch arrays are now a thread-local workspace reset per search (only the touched nodes when few were touched), which removes the O(N) initialisation that dominated small destination-limited searches on large graphs. Measured sequentially on an M4 Max: 1.4–1.8x on weighted random graphs, 1.1–1.3x on unit-cost graphs, 2.4–3.0x on near-destination queries; path graphs (frontier of one node) are about 1.3x slower in absolute terms of a few microseconds, and graphs that fall back to the heap (`W > 65536`) are within a few percent of the previous heap on full-source searches while their small destination-limited searches are still about 2.8x faster. `StrictMultiDiGraph` exposes `max_cost()` / `cost_gcd()`; the C++ functions take an optional trailing `SpfQueue` for tests, and `NGRAPH_CORE_SPF_QUEUE=heap|bucket` overrides the choice process-wide. The search also skips settled neighbours before scanning their parallel edges and sizes the predecessor output from live per-node counts, which makes dense unit-cost fabrics about 1.3x faster; a frontier of one node is held inline, which halves the path-graph penalty.
- **Max-Flow**: `calc_max_flow` and `FlowState::place_max_flow` with `require_capacity=false` computed the identical cost-only shortest-path DAG on every tier; it is now computed once per call. Results are unchanged.
- **Internals**: `FlowPolicy` no longer copies the path DAG when returning a bundle, and the `spf` / `spf_to` bindings move the DAG into the Python result instead of copying it.

## [0.10.0] - 2026-09-13

### Fixed

- **Flow Placement**: `EQUAL_BALANCED` built its split set from the DAG edges that still had residual, so re-placing on a DAG that contained a member filled by an earlier placement silently re-hashed over the remaining members. A forwarding table does not react to load: the split set is now every DAG edge with capacity, and a full member drives the single admission scale to 0, so nothing more is admitted on that DAG. Progression past a full member comes from recomputing the DAG with a residual-aware SPF, which `place_max_flow` and `FlowPolicy` already do, so every documented path is unchanged. The one observable change is `place_max_flow` with `require_capacity=false` and `shortest_path=false`, which filled a cost-only DAG in several passes and now places once, the single-pass answer that cost-only routing is documented to give.

### Removed

- **BREAKING**: `EQUAL_BALANCED_FIXED`, added in 0.9.0 as a separate mode with the behaviour above. `EQUAL_BALANCED` now has that behaviour, so the value is gone; `EQUAL_BALANCED_LOSSY` is now value 3.

## [0.9.0] - 2026-09-13

### Added

- **Shortest Paths**: `shortest_paths_to(dst, ...)` / `Algorithms.spf_to` runs Dijkstra over the in-adjacency and returns distances *to* one destination plus a forward-oriented `PredDAG` valid for placement from any node toward it (the union of the per-source shortest-path DAGs). `fanout_edges` adds edges regardless of cost, for an origin whose first hop is a traffic split rather than a routing decision, such as a pseudo source whose demand originates evenly at every attached real source; a fan-out edge must leave a node with no incoming DAG entry, which keeps the result acyclic.
- **Flow Placement**: two equal-balanced modes for hop-by-hop ECMP with a load-blind forwarding table. `EQUAL_BALANCED_FIXED` takes the split set from the DAG edges with *capacity* rather than residual and keeps the single global admission scale, so a member filled since the DAG was computed drives the scale to 0 (lossless hash-ECMP admission; the existing `EQUAL_BALANCED` drops such members from the split, which is progressive/TE behaviour). `EQUAL_BALANCED_LOSSY` splits over the same set without a scale: each edge carries `min(share, residual)`, the excess is dropped, deficits propagate downstream, and the placed amount is what reaches the destination. `FlowGraph.place_with_drops` returns the dropped volume per edge; `FlowState::place_on_dag` gained an optional `drops` collector. On a 100/10 parallel pair offered 100 units, `EQUAL_BALANCED` and `EQUAL_BALANCED_FIXED` admit 20 and `EQUAL_BALANCED_LOSSY` delivers 60. `FlowPolicy` treats the new modes as equal-balanced except that the lossy mode skips the equalizing rebalance, so with pinned routes each LSP carries what fits and `placed` is the delivered total.

### Changed

- **Flow Policy**: a cost-only policy (`require_capacity=false`) no longer passes the residual to SPF when an equal-balanced per-flow target is set. A residual forces capacity-aware edge selection, so such policies silently routed around saturated edges instead of following costs; they now route on cost alone, and the per-flow target still bounds how much each flow requests.

- **Flow Policy**: `get_path_bundle` now memoizes raw SPF results keyed by the exact inputs that can vary per call (src, dst, residual content, residual-awareness). EqualBalanced placement and rebalance rounds re-request bundles against residual state that repeats -- 94% of SPF calls in a measured place/rebalance cycle were exact input repeats, largely remove+place round-trips restoring identical bytes -- and those calls are now elided. Matching is exact (FlowGraph state stamp fast path, full residual `memcmp` content path), so all outputs are bit-identical; a corpus hash over max-flow, SPF, KSP and policy outputs is unchanged. Measured on place/rebalance churn: 31-70% faster; single EqualBalanced placement: 6-26% faster; Proportional mode skips the memo (it measured 0% repeat inputs) and is unaffected. `FlowGraph` gained an internal monotonic state stamp to support this; no public API change.

## [0.8.0] - 2026-08-24

### Added

- **Flow Policy**: `set_static_paths(src, dst, paths)` pins a demand to explicit path bundles (MPLS-style routing), finishing a feature that was ported but never bound and whose port bound every flow to the first bundle. One flow per usable bundle, each bound to its own bundle in supply order; bundles are validated against the graph and pruned against the policy's masks, and a bundle with no surviving `src->dst` walk is down (no reroute). `EqualBalanced` spreads over the up bundles only. With static paths the policy neither grows its flow set nor reoptimizes; `max_path_cost`/`-factor`, `min_flow_count` and `reoptimize_flows_on_each_placement` are inert.
- **Shortest Paths**: `PredDAG.from_edges(graph, edges)` builds a single-path `PredDAG` from a contiguous edge-id sequence (`make_path_dag` in C++) -- the missing constructor for operator-defined explicit paths. Extracting it also deduplicates the two identical path-to-DAG conversion blocks in `k_shortest_paths` (outputs unchanged).

### Fixed

- **Max-Flow**: `calc_max_flow` could return less than the true maximum, since the tier loop augments only along forward SPF DAGs and never cancels an earlier placement (the reported `min_cut` then contradicted `total_flow`). Added a residual completion phase with reverse arcs. **This is the default configuration** of `max_flow`/`batch_max_flow` (`Proportional` + `require_capacity=True` + `shortest_path=False`), so recorded expected values, regression baselines and stored results must be regenerated: `total_flow` increases where the previous value was suboptimal, and `min_cut`, `edge_flows`, `residual_capacity` and `reachable_nodes` change with it. `sensitivity_analysis`, which is built on `calc_max_flow`, is affected too.
- **Shortest Paths**: Zero-cost edges produced a cyclic predecessor DAG, making `EqualBalanced` placement return 0 flow and `resolve_to_paths` hang. Equal-cost predecessors are now recorded only while a node is unsettled, and `resolve_to_paths` guards against cycles in caller-supplied DAGs. Side effect on graphs containing zero-cost edges: among nodes mutually reachable at equal distance the settle order now decides which alternatives are kept, so a `PredDAG` can hold fewer equal-cost predecessors than in 0.7.2 and ECMP fan-out narrows accordingly. Graphs with strictly positive costs are unaffected (verified bit-identical).
- **K-Shortest Paths**: Spur enumeration materialized every equal-cost path, exponential in ECMP fan-out (a 67-node ECMP ladder needed ~4 s and ~3.5 GiB for `k=3`, and a 73-node one ~19 s; both are now instant in ~34 MB). Enumeration is now bounded by the number of candidates still acceptable; tie-breaking among equal-cost paths may differ, but path counts and costs are unchanged.
- **Flow Policy**: `place_demand` silently routed a second `(src, dst)` pair over the first pair's paths, and neither `place_demand` nor `rebalance_demand` checked that the supplied `FlowGraph` wrapped the policy's own graph -- with equal edge counts that selected paths on one topology and placed flow on another. Both now raise `invalid_argument` (`ValueError` from Python); use one policy per demand, or call `remove_demand()` first to retarget it.
- **Graph Construction**: `from_arrays` now rejects a total edge cost at or above 2^62, which overflows the int64 path arithmetic in SPF and silently corrupts results. The bound is on the sum of all edge costs -- a conservative upper bound on any path -- so it can reject a graph whose individual paths could not actually overflow.
- **Flow State**: `compute_min_cut` and max-flow reachability derived placed flow from `capacity - residual`, overstating it for a `FlowState` built with a custom `residual_init`. Both now use `edge_flow`.
- **K-Shortest Paths**: `k_shortest_paths` read `paths.back()` on an empty vector when `max_cost_factor < 1.0` put the cost ceiling below the shortest path, crashing for `k > 1` (the `k == 1` early return masked it). All `k` now return no paths for such a factor.
- **Python Bindings**: `batch_max_flow` rejected a valid int32 `pairs` array on Windows, because the dtype check compared buffer format strings and NumPy spells int32 as `NPY_LONG` on LLP64 but `NPY_INT` on LP64. It now compares dtype equivalence, and rejects non-contiguous `pairs` rather than reading them as if packed.
- **Build**: `NETGRAPH_CORE_SANITIZE` passed its flags as one quoted string, so sanitizer builds never compiled, and the test target was missing them entirely. `make sanitize-test` also no longer sets `detect_leaks=1` on macOS, where it aborts.
- **Build**: Coverage builds now use `-fprofile-update=atomic` on GCC; the default non-atomic counters corrupt under the threads used by `batch_max_flow` and `sensitivity_analysis`.

### Changed

- **Max-Flow**: `FlowSummary.costs`/`flows` is now a cost-weighted breakdown of the total flow rather than a list of path costs. Entries added by the new completion phase are *marginal* costs -- the augmenting path's forward edge costs minus the cost of the flow it cancels -- so an entry need not correspond to any traversable `src->dst` path. Consumers that treat these as path costs (for example a cost-distribution report) should be reviewed.
- **Flow Policy**: `max_path_cost` no longer implies a relative factor of 1.0. Previously `max_path_cost_factor` defaulted to 1.0 whenever either limit was set, so setting `max_path_cost` alone also rejected every path costing more than the best path found -- making the absolute limit unreachable. The two limits are now independent, and a policy with only `max_path_cost` set admits every path up to that cost. (Found while fixing undefined behavior in the same expression: with an unreachable destination, `best_path_cost_` is the `INT64_MAX` sentinel and multiplying it by the factor overflowed the cast back to `int64`.)
- **Shortest Paths**: Replaced nested per-node predecessor vectors with flat intrusive lists, removing the allocation churn that dominated the hot path (~3-4x faster; output unchanged).
- **Flow Placement**: `place_on_dag` now rebuilds its `(parent, child)` edge groups into a reused arena instead of nested per-node vectors, which was ~70% of its runtime (~2-2.4x faster; output unchanged).
- **Max-Flow**: Parallelized `batch_max_flow` across source/destination pairs using `std::async`; workers claim pairs from a shared counter, so a batch whose costs are unevenly distributed still parallelizes. Thread count from `NGRAPH_CORE_BATCH_THREADS` env or hardware concurrency; set it to `1` when calling from your own worker pool to avoid oversubscription.
- **BREAKING** **C++ API**: `FlowPolicy::set_static_paths` changed signature from `(std::vector<std::tuple<NodeId, NodeId, PredDAG, Cost>>)` to `(NodeId src, NodeId dst, std::vector<PredDAG> bundles)`, and its behavior changed: the previous implementation handed every flow the first matching bundle instead of one bundle per flow. It was never bound to Python, so only C++ callers are affected.
- **Flow Policy**: `place_demand` computes one SPF when seeding initial flows instead of repeating an identical one per flow; seeding cost no longer scales with `min_flow_count`.
- **Flow Policy**: The `EqualBalanced` rebalance recursion (`place_demand` -> `rebalance_demand` -> `place_demand`) is now an iterative loop with numerically identical results (returned leftover may differ by ~1 ulp from floating-point summation order); with many pinned bundles of heterogeneous capacity the recursion depth grew like `U * ln(imbalance/kMinFlow)`, a stack-overflow risk on worker threads.
- **Python Bindings**: Three input errors now raise where they previously did not, or raise a different type -- adjust `except` clauses accordingly. `Algorithms.spf` raises `TypeError` (was `ValueError`) for a wrong `residual` length, matching every other length check; `batch_max_flow` raises `ValueError` for out-of-range node ids that previously produced a zero-flow summary absorbed silently into the batch; and `Algorithms.ksp` validates `dtype` up front, where an invalid value previously went unnoticed whenever the query returned no paths.
- **Python Bindings**: A wrong-typed `graph`/`algorithms` argument now raises `TypeError` instead of an opaque `RuntimeError: Unable to cast ... to C++ type '?'`.
- **Python Typing**: Corrected `_docs.py` (pybind11 enums and classes are not `enum.Enum`/`dataclass`, `MinCut.edges` is an `int32` array, removed a `Path` class that never existed, and added the constructors for `FlowIndex`, `FlowPolicyConfig` and the three enums) and widened the type-checking re-exports from 6 names to all 16, so `Algorithms`, `FlowPolicy`, `StrictMultiDiGraph` and friends are no longer `Unknown` to type checkers despite the shipped `py.typed`. A missing stub constructor is now worse than no stub, so a test asserts every stub constructor matches its binding.
- **Tests**: Several assertions could never fail and now do: a self-loop assertion sat inside `except Exception: pass`, a thread-safety test ended in `assert True` while its workers called `pytest.fail()` from a thread where it cannot fail the test, and two tautologies (`assert size >= 0`) were replaced with the properties actually under test.
- **Docs**: Corrected README/CONTRIBUTING (required CMake is 3.23 not 3.15, `make py-test` does not exist, `make test` does not collect coverage), the `sensitivity_analysis` description in README and `backend.hpp` (it measures flow *lost* on edge removal, not gain from relaxing capacity), the zero-copy and GIL claims, and header contracts for `from_arrays`, `PredDAG`, `FlowSummary.costs`, `calc_max_flow` and `FlowPolicy`.

### Removed

- **Dead code**: `.gcovr.cfg` (never read -- gcovr's default config name has no leading dot, and `make cov` passes every flag explicitly), the unreachable `FlowGraph` include and `expect_flow_conservation()` helper in the C++ test utilities, unused test helpers/fixtures, `_USE_MATH_DEFINES` (no `M_*` macro is used), and unused `<optional>`/`<stack>`/`<unordered_map>`/`<unordered_set>`/`<cmath>` includes across `src/` and `include/`.
- **BREAKING** **Python Bindings**: `Algorithms.build_graph_from_arrays`. Nothing used it (no test, and not NetGraph), and the `Graph` handle it returned could not be passed to `FlowGraph` or `FlowPolicy`, so it only ever served the stateless algorithms while looking like a peer of `build_graph`. Build the graph with `StrictMultiDiGraph.from_arrays()` and pass it to `Algorithms.build_graph()`.
- **BREAKING** **C++ API**: The `build_graph(std::shared_ptr<const StrictMultiDiGraph>)` overload on `Algorithms` and `Backend`, which existed solely to serve that binding. It was a pure virtual, so an out-of-tree `Backend` implementation must delete its now-dangling `override` to compile. Callers wanting the handle to own the graph construct it directly: `GraphHandle{my_shared_ptr}`.
- **BREAKING** **Python Bindings**: The `Flow` class (bound from `FlowRecord`). It was reachable only as `_netgraph_core.Flow`, never exported from `netgraph_core`, and no binding constructed or returned one, so nothing could obtain an instance; the name also collided with the C++ alias `using Flow = double`.

## [0.7.2] - 2026-03-26

### Fixed

- Added missing `<stdexcept>` include for MSVC compatibility

## [0.7.1] - 2026-03-26

### Added

- Windows AMD64 wheel builds in release workflow
- Windows test coverage in CI

## [0.7.0] - 2026-03-26

### Changed

- Relicensed from BSD-3-Clause to MIT
- Parallelized sensitivity analysis across candidate edges using `std::async`; thread count controlled by `NGRAPH_CORE_SENSITIVITY_THREADS` env or hardware concurrency
- Thread-local profiling stats to remove global mutex from the hot `record()` path; public API unchanged

## [0.6.0] - 2026-02-26

### Changed

- **License**: Relicensed from GPL-3.0-or-later to BSD-3-Clause.

## [0.5.0] - 2026-02-18

### Changed

- **License**: Relicensed from AGPL-3.0 to GPL-3.0.

## [0.4.1] - 2026-02-17

### Added

- **Build**: Add Linux ARM64 (aarch64) wheel builds using native GitHub ARM runners in the release workflow.

## [0.4.0] - 2026-02-17

### Added

- **Python 3.14**: Add support for Python 3.14 in CI tests, wheel builds, and packaging.
- **Free Threading**: Declare free-threading compatibility (`py::mod_gil_not_used`) for Python 3.13t/3.14t builds. Build and publish free-threaded (`cp314t`) wheels.

### Changed

- **CI**: Upgrade cibuildwheel from v2 to v3.3.1 for Python 3.14 wheel builds.

### Fixed

- **Python Bindings**: Fix object leak on Python 3.14 caused by pybind11 not clearing managed dicts during deallocation. Replace `dynamic_attr` + `_graph_ref` pattern with `py::keep_alive` for FlowState, FlowGraph, and FlowPolicy.

## [0.3.6] - 2026-02-17

### Fixed

- **Python Bindings**: Fix build failure with pybind11 3.x typed tuples by adding explicit `-> py::tuple` return type on the `spf` lambda.

### Changed

- **CI**: Add CentOS Stream 9 build and test jobs.

## [0.3.5] - 2026-02-08

### Fixed

- **K-Shortest Paths**: Fix off-by-one in Yen's prefix comparison and PredDAG CSR fill order. Paths visiting nodes out of numerical order (e.g., `[0,3,2]`) previously produced malformed predecessor DAGs.

## [0.3.4] - 2026-02-02

### Changed

- **Build**: Simplified version handling.

## [0.3.3] - 2026-02-02

### Fixed

- **License**: Replaced incomplete AGPLv3 LICENSE file with complete official text from GNU.

## [0.3.2] - 2025-12-12

### Fixed

- **Profiling**: Fix ODR violation that caused empty stats when profiling was enabled. Moved `profiling_enabled()` and `ProfilingStats::instance()` definitions from inline header to `profiling.cpp` to ensure a single instance when static library is linked into the Python module.

## [0.3.1] - 2025-12-08

### Added

- **Profiling**: Runtime profiling infrastructure for C++ hot paths (`shortest_paths_core`, `place_demand`, `place_on_dag`).
  - Enable via `NGRAPH_CORE_PROFILE=1` environment variable.
  - Python API: `profiling_enabled()`, `profiling_dump()`, `profiling_reset()`.
  - Minimal overhead when disabled (single static bool check per instrumented scope).
  - ~2% overhead when enabled.

### Changed

- **Build**: Default optimizations: LTO, loop unrolling, `-fno-math-errno`. Add `make install-native` for CPU-specific builds.

## [0.3.0] - 2025-12-06

### Changed

- **BREAKING**: Minimum Python version raised to 3.11 (was 3.9)

## [0.2.3] - 2025-12-06

### Changed

- **Python bindings**: `StrictMultiDiGraph.from_arrays` now requires `ext_edge_ids` so callers always supply stable external edge identifiers.
- **FlowPolicy**: construction is now config-only (via `FlowPolicyConfig`), dropping the parameter-heavy constructor.

### Fixed

- **FlowGraph**: `get_flow_path` now filters only below-`kEpsilon` noise so paths are reconstructed even when per-edge allocations are smaller than `kMinFlow`.

## [0.2.2] - 2025-12-05

### Fixed

- **Flow Placement**: EqualBalanced placement now correctly returns 0 when the shortest path has no capacity with `require_capacity=False`. Previously, flow could be incorrectly reported on partial paths that didn't reach the destination.

## [0.2.1] - 2025-12-01

### Fixed

- **Shortest Paths**: In single-path mode, ties between equal-cost paths are now broken by preferring higher bottleneck capacity. Improves flow placement when multiple equal-cost paths exist with different capacities.
- **Flow Placement**: Use epsilon threshold in `place_on_dag()` to fix placement of very small flow fractions in large fanout networks.

### Changed

- **Build**: Use `uv` build frontend for wheel builds.
- **Build**: Drop 32-bit Linux (i686) wheels.

## [0.2.0] - 2025-11-25

### Added

- **Sensitivity Analysis**: Added `shortest_path` parameter to `sensitivity_analysis()`.
  - `shortest_path=False` (default): Uses full max-flow; reports all saturated edges across all cost tiers.
  - `shortest_path=True`: Uses single-tier shortest-path flow; reports only edges used under ECMP routing.
- Python type stub documentation for `Algorithms.sensitivity_analysis()`.

## [0.1.0] - 2025-11-23

### Added

- **Core Library**: Initial release of C++ implementation for graph algorithms and flow tracking.
- **Graph Structures**:
  - `StrictMultiDiGraph`: Immutable directed multigraph using CSR (Compressed Sparse Row) adjacency.
  - `FlowGraph`: Manages flow state, per-flow edge allocations, and residual capacities.
- **Algorithms**:
  - Shortest paths (Dijkstra variant returning a DAG for ECMP; supports node/edge masking and residual-aware tie-breaking).
  - K-Shortest paths (Yen's algorithm).
  - Max-flow (Successive Shortest Path with ECMP/WCMP placement; supports capacity-aware (TE) and cost-only (IP) routing modes).
  - Sensitivity analysis (identifies bottlenecks).
- **Flow Policy**:
  - **Modeling**: Unified configuration for IP routing (cost-based ECMP) and Traffic Engineering (capacity-aware TE).
  - **Placement**: `Proportional` (WCMP) and `EqualBalanced` (ECMP) strategies.
  - **Lifecycle**: Manages demand placement, static/dynamic path selection, and re-optimization.
  - **Constraints**: Enforces limits on path cost, stretch factor, and flow counts.
- **Python Bindings**:
  - Python 3.9+ support via pybind11.
  - NumPy integration using zero-copy views where applicable.
  - Releases GIL during long-running graph algorithms.
- **Testing**:
  - Python and C++ test suites.
