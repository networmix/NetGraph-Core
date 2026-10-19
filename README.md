# NetGraph-Core

C++ graph engine for network flow analysis, traffic engineering simulation, and capacity planning.

> **Unreleased, breaking.** This version replaces the 0.11 API (see the Unreleased section of
> [CHANGELOG.md](CHANGELOG.md)): the graph on integer units, the minimum-cost flow core, the route
> structures (`PredDAG`, `Path`), the integer ledger `FlowLedger` with the placement kernels, the
> policy engine `FlowPolicy`, the searches (`spf`, `ksp`, `resolve_to_paths`), the `DagCache` in
> building and frozen modes, and the evaluators `capacity`, `batch_capacity` and `sensitivity`.
> The previous API (0.11 and earlier) is in git history and on PyPI.

## Overview

- **Exact integer amounts**: a graph carries an explicit resolution `q` (a power of two in
  `[2^-32, 2^32]`); capacities and flows are `int64` units of `q`. Physical values enter only
  through `to_units` (round half up) and leave through `as_float`.
- **Positive costs**: edge costs are integers `>= 1`, so the shortest-path DAG (every edge on
  some shortest path) has no cycle and does not depend on search order or tie rules.
- **Zero capacity is not absence**: an edge of capacity 0 stays in routing, as a configured but
  empty link does: searches, routes and rosters use it, an equal split still sends it a share
  (which is lost), a weighted split gives it weight 0. Only an explicit residual or an
  `Eligibility` removes an edge from a search.
- **Stable edge ids**: edges are numbered by sorting on `(src, cost, dst)`, ties kept in input
  order; that number is the public `EdgeId`, and `input_index()` maps it back to the input position.
- **One flow core**: an exact integer minimum-cost flow solver (the most traffic, sent as cheaply
  as possible). It returns the flow value, the marginal cost curve (the extra minimum total cost
  per additional delivered unit, as `(cost, units)` segments), the flow on every edge, a minimum
  cut (the full edges that separate the sources from the destinations) and the node potentials
  (per-node prices that prove the flow is cheapest). It works one cost level at a time: shortest
  paths with potentials find the next cheapest cost, then Dinic phases fill every path of that
  cost.

## Core Features

The exact behaviour of every call (units and rounding, edge ids, which route a search returns,
how the kernels split traffic, outcomes and flags, failures, the cache and the evaluators) is
specified in [docs/semantics.md](docs/semantics.md).

### Graph

- **`Graph(num_nodes, src, dst, capacity_units, cost, resolution)`**: immutable
  directed multigraph that stores each node's out-edges and in-edges as arrays; parallel edges and
  self-loops are allowed. Validation reports the first defect in a fixed order (the resolution, the unit
  values, the shapes, zero or negative costs, the remaining bounds). Accessors (`capacity_units()`, `cost()`, `edge_src()`, `edge_dst()`,
  `input_index()`, `row_offsets()`, `in_row_offsets()`, `in_adj_edge_index()`) return read-only
  NumPy views that keep the graph alive.
- **`Eligibility`**: which nodes and edges of one graph may be used (immutable), built by
  `graph.eligibility(excluded_nodes, excluded_edges)`; `Eligibility.intersect(a, b)` returns an
  operand when the result equals it, otherwise a new value.

### Route structures

- **`PredDAG`**: the route structure a search returns, a predecessor DAG (directed, no cycle):
  for each node, the edges that reach it on the route and the node each comes from. `parents()`
  and `via_edges()` list these entries, in ascending EdgeId within a node; `parent_offsets()`
  (built on each call) says where each node's entries start.
  `PredDAG.from_csr(graph, offsets, parents, via)` validates the shape, the ranges and endpoints,
  duplicates, acyclicity and the canonical order before accepting a DAG.
- **`Path(graph, hops)`**: a route given hop by hop, where a hop may hold several parallel edges
  between the same two nodes;
  `cost()` is the sum over hops of the minimum member cost and `to_dag()` gives the equivalent
  `PredDAG`.
- **`PathSelection(tie_break=DETERMINISTIC, priority=CAPACITY, multipath=None,
  multi_edge=True)`** with `EdgeTieBreak` and `EdgePriority`: the selection axes the searches and
  policies read (`multipath=None` is the model's natural shape: the DAG of equal-cost walks for
  searches and CostOnly, one path for a tunnel).

### Ledger and kernels

- **`FlowLedger(graph, residual=None, *, unreserved=None)`**: the integer ledger of one graph,
  in two projections that are never added: the carried traffic (a service baseline `residual`;
  `allocated()` is the exact sum of the records' edges and `residual()` the baseline minus it)
  and the strict holds of T1 tunnels (a reservation baseline `unreserved`; `reserved()` is the
  sum of the records' holds and `unreserved()` the baseline minus it). Holds never police
  traffic. The records (`FlowRecord`: flow id, owner, endpoints, admission, program entry,
  bundle, route cost, offered, delivered, per-edge allocation, drops and holds) are the flows
  placed on it; `bottleneck(sources)` lists the full edges leaving the set of nodes the sources
  can still reach. `reset()` and `copy()` are rejected while a policy is bound to the ledger.
  Policies place, hold and remove; `reset(residual=None, *, unreserved=None)` replaces an
  unbound ledger's baselines and clears its records.
- **Kernels**, the rules that spread one flow over its route (`FlowPlacement`: `MAX_FLOW`,
  `WEIGHTED`, `EQUAL_BALANCED`; `Admission`: `LOSSLESS`, `LOSSY`): `MAX_FLOW` re-solves the flow
  with the flow core on its DAG, over the free capacity plus its own allocation; `WEIGHTED` splits
  at every node in proportion to integer weights frozen when the route is prepared, from a
  `WeightBasis` (`CAPACITY`, `EXPLICIT`, `RESIDUAL_SNAPSHOT`) (highest averages: whole-unit shares
  as close to proportional as possible); `EQUAL_BALANCED` splits equally (round
  robin); `LOSSY` forwards hop by hop and records what it drops.

### Policies

- **`FlowPolicy(ledger, sources, destinations, config, elig=None, cache=None)`**: the placement
  engine for one demand, bound to its ledger. `FlowPolicyConfig(*, routing, placement=None,
  admission=LOSSLESS, selection=None, origin=None, weights=None, edge_weights=None)` chooses the
  routing, its kernel (required for CostOnly, taken by no other routing), the admission, the path
  selection (default `PathSelection()`), the weight basis of CostOnly `WEIGHTED` and the origin: `EQUAL` (every source sends the same round-robin share) or
  `OPTIMIZED` (a controller chooses how much enters where; the default of the controllers). The
  routing is one of:
  - `CostOnly()`: static shortest paths, as routers forward with an IGP such as OSPF; the
    sources of one demand are pooled (where their routes meet, a router splits their traffic
    once); with `MAX_FLOW`, a controller restricted to those paths;
  - `CapacityAware(Tunnels(N), max_path_cost=None, min_unreserved_bandwidth=0)` (T1): every
    head sets up its own N strict-bandwidth tunnels, as RSVP-TE head ends do, with no controller:
    a design volume is split equally over the heads and their tunnels, each tunnel runs CSPF (the
    shortest path among the links with enough unreserved bandwidth, within an `Absolute` or
    `Stretch` cost ceiling) and either comes up holding its request on every link or stays down
    (`NO_BANDWIDTH`, `COST_CEILING`); traffic then enters equally at the heads and splits equally
    over their up tunnels;
  - `Global(support=None)`: an ideal controller: the minimum-cost flow from all the policy's
    sources to its destinations, on the free capacity plus its own allocation, optionally only
    over the edges of a support DAG;
  - `Pinned(program, black_holes=None)`: forwarding over a given `ForwardingProgram`, computing
    nothing; the policy's Eligibility fails the tunnels it touches (their heads split over the
    survivors) and `black_holes` makes failed elements swallow what reaches them (stale
    forwarding).
- **`ForwardingProgram`**: the routes and integer weights of one demand. CostOnly, Global and T1
  produce one (`policy.program`, `result.program`); `ForwardingProgram.from_routes(graph,
  sources, destinations, routes)` builds one from `Route(head, Path or PredDAG, weight=1)`
  entries.
- **`place(volume, outputs=None)`** adds `volume` to the policy's request as one all-or-nothing
  step (CostOnly forwards over the program it prepares once per epoch; the controllers re-solve
  what they already deliver plus this request; T1 has no `place`) and returns
  a `PlacementResult`: this call's `placed`, `remaining`, `delivery_delta`, `call_drops`,
  `call_source_loss`, `termination` and flags; the demand's `requested`, `offered`, `admitted`,
  `delivered`, `unadmitted`, `source_loss`, `drops`, per-source rows `by_source`, `total_cost`
  (controllers), `cost_distribution` (optimized Global), the `program`, T1's `plan` and
  `reservations_changed`, and per-tunnel rows `by_tunnel`; and, when requested through
  `OutputOptions`, `edge_flows`, `residual`, `bottleneck`, `reservation_edge_flows` and `min_cut`
  (optimized Global).
  `evaluate(offer, design_demand=None)` computes a whole new offer without publishing it,
  `replace(offer, design_demand=None)` publishes it (T1 sets its tunnels up for the design,
  `None` meaning the offer), `forward(offer, black_holes=None)` sends it over the installed
  program (T1 keeps its tunnels and holds); `remove()` reverts the policy's flows and holds; `records()`, `record(flow_id)`, `placed_demand()` and `program` read its
  state.
- **Failure guarantee**: every public call either completes or raises one exception and leaves
  every ledger, policy, cache and scratch buffer as it was, whatever fails (validation, memory, a
  worker thread, building the Python result). Two limits: pybind11 does not check every
  allocation when it constructs a bound object from Python or converts a returned value into a
  new Python object outside the library's checked result construction, so running out of memory
  exactly there may crash the process (Core state stays unchanged); and a `DagCache` passed to the
  call: a building cache keeps the complete entries and statistics of the lookups that succeeded
  before the failure, and a frozen cache keeps its entries and statistics unchanged, except that a
  call raising `DagCacheMiss` counts one miss
  ([docs/semantics.md](docs/semantics.md#8-failures-threads-and-lifetimes)).

### Search

`spf(graph, roots, SpfOptions(direction, selection, targets=None, residual=None,
min_residual=1, elig=None))` searches from a root set, forward or reverse.
- **Without targets:** it returns `costs` (`UNREACHABLE` for nodes not reached) and the full
  DAG of tight edges, the edges that lie on some shortest path (one per node pair with
  `multi_edge=False`).
- **With targets:** it stops at the nearest targets and returns `target_cost` and the DAG of
  the walks between the roots and the nearest targets (for a reverse search, from a target to
  a root).
- **With `multipath=False`:** it returns the canonical single path: the walk of lowest cost;
  among those, under `PREFER_HIGHER_PRIORITY` only, the one of highest bottleneck (the smallest
  priority value along the walk); then the fewest hops; then the lexicographically smallest
  sequence of EdgeIds ([docs/semantics.md](docs/semantics.md#the-canonical-single-path)).
- **With `multi_edge=False`** (any search): each node pair is first reduced to one parallel
  member (lowest cost, then under `PREFER_HIGHER_PRIORITY` the highest priority, then the lowest
  EdgeId).

`ksp(graph, src, dst, KspOptions(k, selection, max_path_cost=None, residual=None, elig=None))`
returns the K shortest simple paths as `(cost, Path)` pairs, in the order of
[docs/semantics.md](docs/semantics.md#k-shortest-paths) and up to
an inclusive cost ceiling. `resolve_to_paths(graph, dag, src, dst, split_parallel_edges=False,
max_paths=None)` lists a DAG's walks as `Path` values.

### DagCache

`DagCache(graph)` keeps the static shortest-path DAGs (from costs and the Eligibility, not
traffic) that policies build, keyed by the root set, the direction and the Eligibility. Values
are immutable and shared, and `stats()` reports hits, misses, entries, retained bytes and live
values.
- **Building mode (new caches):** one writer; a miss searches and inserts. An evaluator given a
  building cache runs on one worker (an explicit `workers` above 1 is rejected).
- **Frozen mode (`freeze()`):** read-only and shared by every worker of an evaluator (pass
  `workers` explicitly; with a cache the default is 1); a miss raises `DagCacheMiss`. `thaw()`
  returns to building mode once every reader has finished.

### Evaluators

All three take one options type, `QueryOptions(config, *, elig=None, cache=None,
outputs=None, design_demand=None, offered=None)`. `capacity(graph, sources, destinations, opts)` is `place(+inf)` of a policy on a fresh ledger
(a Lossy Pinned policy places the finite `offered` instead; T1 evaluates the tunnels it sets up
for the required finite `design_demand`, forwarding `+inf` or, Lossy, the finite `offered`); `Global()` gives the
exact maximum flow with its marginal cost curve, and its minimum cut with
`outputs=OutputOptions(min_cut=True)`. `batch_capacity(graph, pairs,
opts, per_pair_elig=None, workers=None)` evaluates pairs in pair order (the lowest failing
index raises; Pinned, whose program fixes its endpoints, is rejected). `sensitivity(graph, sources, destinations, opts, workers=None)` reports, for every candidate edge, the signed loss in placed units when the
edge is removed. With `workers=None` the worker count is 1 when a cache is supplied, otherwise
it comes from `NGRAPH_CORE_BATCH_THREADS` and `NGRAPH_CORE_SENSITIVITY_THREADS`.

### Flow core

One native integer minimum-cost flow core serves Global routing, the MAX_FLOW kernel and the
evaluators: the largest flow up to a cap from a source set to a destination set over the eligible
edges, sent at the lowest total cost, starting from zero flow. For sensitivity it repairs the
baseline's solution for each removed edge instead of solving again. Its private Python entry
`_netgraph_core._place_global_units(graph, view, sources, destinations, elig, cap)` returns the
payload `{placed, cost_distribution, edge_flows, min_cut, potentials, counters}` for the
flow checks of the test suites.

## Quick Start

```python
import netgraph_core as ngc

q = 2.0**-12
# Two parallel paths from node 0 to node 3 (physical capacities 10 and 5).
caps = ngc.to_units([10.0, 10.0, 5.0, 5.0], q)
graph = ngc.Graph(4, [0, 1, 0, 2], [1, 3, 2, 3], caps, [1, 1, 2, 2], q)

config = ngc.FlowPolicyConfig(routing=ngc.Global())
result = ngc.capacity(graph, 0, 3, ngc.QueryOptions(config))
print(ngc.as_float(result.placed, q))   # 15.0
print(result.cost_distribution)         # [(2, 40960), (4, 20480)]
print(result.total_cost)                # 163840

# T1: two tunnels sized for 12; the second finds no path with 6 units unreserved.
t1 = ngc.FlowPolicyConfig(routing=ngc.CapacityAware(ngc.Tunnels(2)), admission=ngc.Admission.LOSSY)
policy = ngc.FlowPolicy(ngc.FlowLedger(graph), 0, 3, t1)
r = policy.replace(ngc.to_units(12.0, q), design_demand=ngc.to_units(12.0, q))
print([t.status.name for t in r.program.tunnels])  # ['UP', 'NO_BANDWIDTH']
print(ngc.as_float(r.placed, q))                    # 10.0
print(ngc.as_float(r.plan.reserved_total, q))       # 6.0
```

## Installation

```bash
pip install netgraph-core
```

Or from source (the shipped configuration; `pip install -e .` is the editable development build,
which carries the private test hooks):

```bash
pip install .
```

The package is typed: it ships `py.typed` and a type stub generated from the compiled module.

### Build Optimizations

Default builds include LTO and loop unrolling. For local development:

```bash
make install-native   # CPU-specific optimizations (not portable)
```

## Repository Structure

```
src/                    # C++ implementation
include/netgraph/core/  # Public C++ headers
bindings/python/        # pybind11 bindings
python/netgraph_core/   # Python package
tests/cpp/              # C++ tests (googletest)
tests/py/               # Python tests (pytest)
dev/                    # Developer scripts (not shipped): checks, coverage,
dev/perf/               #   and the A/B/A timing tool
```

## Development

```bash
make dev           # Setup: venv, dependencies, pre-commit hooks
make check         # Rebuild the editable extension, then everything check-ci runs
make check-ci      # Lint, types, C++ tests (hooks on and off), Python tests (CI entrypoint)
make test          # Python tests
make cpp-test      # C++ tests (private hooks on)
make cpp-test-off  # C++ tests without the hooks (the shipped configuration)
make sanitize-test # C++ tests under ASan/UBSan
make tsan-test     # C++ tests under ThreadSanitizer
make cov           # Combined coverage report (C++ + Python; its own venv and builds)
make stubs-check   # The generated type stub: drift, package contents, typed clients
make format        # Apply ruff formatting (make hooks runs every pre-commit fixer)
```

`make check` and `make check-ci` only verify; they never edit a file.

Performance claims are timed with the A/B/A tool (old build, new build, then the old build again
on one machine, so drift shows), described in [dev/perf/README.md](dev/perf/README.md)
(`make perf-tool-test` tests the tool itself; CI runs it on one macOS job).

The editable development build (`make install`, `make dev`) carries the private test hooks
(CMake option `NETGRAPH_CORE_PRIVATE_HOOKS`; a wheel build does not); the test suites skip the
tests that need them when they are absent.

## Environment Variables

| Variable | Effect |
| --- | --- |
| `NGRAPH_CORE_BATCH_THREADS` | Default worker count of `batch_capacity` (an integer in [1, 4096]; 1 when unset or invalid; not read when a cache is supplied, which defaults to 1; an explicit `workers` wins). |
| `NGRAPH_CORE_SENSITIVITY_THREADS` | Default worker count of `sensitivity` (as above). |

## Requirements

- **C++:** C++20 compiler (GCC 10+, Clang 12+, MSVC 2019+)
- **Python:** 3.11+
- **Build:** CMake 3.23+, scikit-build-core
- **Dependencies:** pybind11, NumPy

## License

[MIT License](LICENSE)
