# NetGraph-Core

NetGraph-Core is a C++ engine with Python bindings that computes how traffic crosses a network:
the maximum flow between two sets of nodes, how ECMP or WCMP shortest-path routing, RSVP-TE
tunnels or an ideal controller would carry a demand, and how much a link failure costs. All
amounts are exact integers, so a result never depends on floating-point tolerance. It is a
library for developers who build network analysis tools. For a ready tool, with networks and
failure scenarios in YAML and reports, use [NetGraph](https://github.com/networmix/NetGraph),
which is built on this engine.

## Installation

```bash
pip install netgraph-core
```

The package needs Python 3.11 or later and NumPy, which pip installs with it. Where pip finds no
prebuilt wheel for your platform and Python version, it builds the package from source, as
`pip install .` does in a checkout. A source build also needs a C++20 compiler (GCC 10+,
Clang 12+, MSVC 2019+) and CMake 3.23+; pip fetches the Python build tools, scikit-build-core and
pybind11. `pip install -e .` is the editable development build, which also carries private test
hooks. The package is typed: it ships `py.typed` and a type stub generated from the compiled
module.

## Quick Start

The engine counts bandwidth in whole units. Each graph fixes the size of one unit, `q`, with its
`resolution` argument. Here `q` is 1 Gbit/s, so capacities and results read directly as Gbit/s
([Bandwidth units and q](#bandwidth-units-and-q) covers other values). Edge costs are positive
integers, and the cost of a route is the sum of its edge costs.

Node A reaches node C directly over 10 Gbit/s, or through B over two 20 Gbit/s edges. Every edge
has cost 1, so the direct route costs 1 and the route through B costs 2. Nodes are numbered from
0: A is 0, B is 1, C is 2.

```text
               B (1)
   20 Gbit/s  / \  20 Gbit/s
             /   \
        A (0)-----C (2)
           10 Gbit/s

        every edge: cost 1
```

How much can A send to C with shortest-path routing, and how much with an ideal controller?

```python
import netgraph_core as ngc

A, B, C = 0, 1, 2
graph = ngc.Graph(
    num_nodes=3,
    edge_sources=[A, A, B],
    edge_destinations=[C, B, C],
    edge_capacity_units=[10, 20, 20],
    edge_costs=[1, 1, 1],
    resolution=1.0,  # q: one unit is 1 Gbit/s
)

shortest_path = ngc.FlowPolicyConfig(routing=ngc.ShortestPaths())
optimized_te = ngc.FlowPolicyConfig(
    routing=ngc.CentralizedTE(source_allocation=ngc.SourceAllocation.OPTIMIZED)
)
sp = ngc.evaluate_demand(graph, A, C, ngc.EvaluationOptions(shortest_path))
te = ngc.evaluate_demand(graph, A, C, ngc.EvaluationOptions(optimized_te))
print(f"Shortest path: {sp.call_delivered} Gbit/s")  # Shortest path: 10 Gbit/s
print(f"Optimized TE: {te.call_delivered} Gbit/s")  # Optimized TE: 30 Gbit/s
print(te.marginal_cost_curve)  # [(1, 10), (2, 20)]
```

`evaluate_demand` answers one question on a fresh copy of the graph at full capacity: here, with
no offer given, the most A can deliver to C under the routing model. Shortest-path routing
forwards only on the lowest-cost routes and splits traffic equally when several tie. In this
triangle the direct edge is the only lowest-cost route, so nothing is split: A delivers 10, and
the route through B carries nothing because it costs more. The ideal controller (`CentralizedTE`
with optimized source allocation) computes the cheapest flow that delivers the most: 10 over A→C
and 20 over A→B→C. That fills both edges leaving A, 10 + 20, so no routing can deliver more
than 30.

A `marginal_cost_curve` lists `(cost, units)` pairs: the first 10 units cost 1 each, the next 20
cost 2 each.

### Bandwidth units and q

`q` is the bandwidth of one unit: with bandwidth in Gbit/s, `q = 1.0` means 1 Gbit/s per unit and
`q = 0.25` means 0.25 Gbit/s per unit. It is the `resolution` argument of `Graph`. Core has no
default `q`; the caller passes a power of two from `2**-32` to `2**32`.

`Graph` takes capacities already in units, and setting `q` does not convert them. The Quick Start
uses `q = 1.0` only because its capacities are whole Gbit/s, which keeps the numbers easy to
read. For physical amounts, convert with `to_units(value, q)` and back with `as_float(units, q)`.
A value between two whole units goes to the nearest one, and an exact halfway value rounds up;
from then on the engine computes exactly on those integers.

```python
q = 0.25
print(ngc.to_units(10.0, q), ngc.as_float(40, q))  # 40 10.0
print(ngc.to_units(2.625, q), ngc.as_float(11, q))  # 11 2.75
print([ngc.to_units(gbps, q) for gbps in (10.0, 20.0, 20.0)])  # [40, 80, 80]
```

At `q = 0.25`, 10 Gbit/s is 40 units, and 2.625 Gbit/s is 10.5 units, rounded up to 11 units, or
2.75 Gbit/s. Keeping the Quick Start's capacity array `[10, 20, 20]` with `q = 0.25` would
describe a different network with a quarter of the bandwidth. To keep the same network, convert
its capacities again, to `[40, 80, 80]`, and do the same for every traffic amount you pass;
results then come back in units of 0.25 Gbit/s.

[NetGraph](https://github.com/networmix/NetGraph) chooses `q` itself when it prepares a scenario.

## The main ideas

- **Units.** Capacities and traffic are whole numbers of units of the graph's `q`
  ([Bandwidth units and q](#bandwidth-units-and-q)).
- **The graph** (`Graph`) is directed and may have parallel links. Each edge has a capacity and
  an integer cost of at least 1. Edges are numbered (their `EdgeId`) by sorting on source, cost
  and destination, so an edge's number can differ from its position in the input arrays. An edge
  of capacity 0 is an empty link, not a missing one: routing still uses it, and traffic sent there
  is lost.
- **A failure** is an `Eligibility`: the nodes and edges a call may not use.
- **A demand** is traffic from a set of sources to a set of destinations. A `FlowPolicy` routes
  one demand on a `FlowLedger`, the record of what each edge carries; several demands share one
  ledger, and `apply_actions` places them one after another or at the same time.
- **The routing model** decides the routes:
  - `ShortestPaths`: routers forwarding on shortest paths, split equally (ECMP) or by
    weights (WCMP), or a controller restricted to those paths;
  - `DistributedTE` (called T1 in the docs): each source router sets up its own N tunnels of
    fixed bandwidth on the shortest path with enough unreserved bandwidth, as RSVP-TE head ends
    do;
  - `CentralizedTE`: an ideal controller, the cheapest flow that delivers the most; with a
    tunnel count (T2), it lays out the same number of single-path tunnels for every source
    router, all of one planned bandwidth;
  - `Pinned`: forwards a given program, to replay a result or study a failure.
- **Admission** says what happens to traffic that does not fit: `LOSSLESS` admits only what goes
  through without loss, `LOSSY` sends everything and reports where it is lost.
- **A program** (`ForwardingProgram`) is what a model decided: the routes and the integer
  weights that split traffic over them. `place`, `evaluate`, `replace` and `forward` return a
  `PlacementResult`: what was delivered and lost, and where.
- **The evaluators** answer one question each on a fresh ledger at full capacity:
  `evaluate_demand` for one demand (with `CentralizedTE` and the optimized source allocation,
  the exact maximum flow, its cost curve and its minimum cut), `evaluate_pairs` for many
  independent pairs in parallel, and `evaluate_edge_failures` for how much delivery each edge's
  removal takes away.
- **Calls are all or nothing.** A failing call raises one exception and leaves every ledger,
  policy and cache as it was, with two limits: running out of memory inside pybind11's own
  object construction may crash the process, and a supplied shortest-path cache may keep
  completed entries and lookup statistics. Native work releases Python's interpreter lock (the
  GIL); calls on one ledger take effect one at a time, and different ledgers run in parallel.

The exact rules, including tunnel setup and reservations, losses on shared edges and the
shortest-path cache, are in [docs/semantics.md](https://github.com/networmix/NetGraph-Core/blob/main/docs/semantics.md).

## Documentation

- [docs/semantics.md](https://github.com/networmix/NetGraph-Core/blob/main/docs/semantics.md): the exact rules (units and rounding, edge numbering,
  which route a search returns, how each model splits traffic, results, failures, the
  shortest-path cache and the evaluators), with worked examples that the tests run.
- The type stub `netgraph_core/__init__.pyi`: every signature, with its docstring.
- The C++ API: the headers under [include/netgraph/core/](https://github.com/networmix/NetGraph-Core/tree/main/include/netgraph/core); their names
  are listed in the last section of the semantics.
- [proofs/](https://github.com/networmix/NetGraph-Core/blob/main/proofs/README.md): machine-checked Lean proofs of the round-robin splitting rule,
  explained in [docs/formal/theorems/](https://github.com/networmix/NetGraph-Core/blob/main/docs/formal/theorems/round-robin.md). They prove the
  mathematical rule, not the C++ code; nothing in the build or the package needs Lean.

## Development

[CONTRIBUTING.md](https://github.com/networmix/NetGraph-Core/blob/main/CONTRIBUTING.md) covers setup, the test targets and the release process. In
short: `make dev` sets up a venv with the editable build and the pre-commit hooks, and
`make check` rebuilds the extension and runs every check without editing a file.
`make install-native` builds with CPU-specific optimizations (not portable). Performance claims
are timed with the A/B/A tool in [dev/perf/](https://github.com/networmix/NetGraph-Core/blob/main/dev/perf/README.md).

## Environment Variables

| Variable | Effect |
| --- | --- |
| `NGRAPH_CORE_BATCH_THREADS` | Default worker count of `evaluate_pairs`: an integer in [1, 4096]; 1 when unset or invalid, or when a cache is supplied. An explicit `workers` wins. |
| `NGRAPH_CORE_SENSITIVITY_THREADS` | The same for `evaluate_edge_failures`. |

## License

[MIT License](https://github.com/networmix/NetGraph-Core/blob/main/LICENSE)
