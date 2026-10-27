# NetGraph-Core

NetGraph-Core is a C++ engine with Python bindings that computes how traffic crosses a network:
the maximum flow between two sets of nodes, how ECMP or WCMP shortest-path routing, RSVP-TE
tunnels or an ideal controller would carry a demand, and how much a link failure costs. All
amounts are exact integers, so a result never depends on floating-point tolerance. It is a
library for developers who build network analysis tools. For a ready tool, with networks and
failure scenarios in YAML and reports, use [NetGraph](https://github.com/networmix/NetGraph),
which is built on this engine.

## The main ideas

- **Units.** Capacities and traffic are whole numbers of units. Each graph carries its
  resolution `q`, the size of one unit (a power of two, for example `2**-12` Gbit/s);
  `to_units` converts a physical value (rounding half up) and `as_float` converts back.
- **The graph** (`Graph`) is directed and may have parallel links. Each edge has a capacity and
  an integer cost of at least 1. Edges are numbered (their `EdgeId`) by sorting on source, cost
  and destination. An edge of capacity 0 is an empty link, not a missing one: routing still uses
  it, and traffic sent there is lost.
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

## Quick Start

Two paths from node 0 to node 3: 10 Gbit/s at cost 1 per hop, and 5 Gbit/s at cost 2 per hop.
One unit is 1 Gbit/s (`q = 1.0`), so every printed amount is in Gbit/s.

```python
import netgraph_core as ngc

q = 1.0
graph = ngc.Graph(4, [0, 1, 0, 2], [1, 3, 2, 3], [10, 10, 5, 5], [1, 1, 2, 2], q)

# The maximum flow: 10 on the cost-2 path, then 5 on the cost-4 path; the minimum cut is the first
# link of each path (edges 0 and 1).
te = ngc.FlowPolicyConfig(routing=ngc.CentralizedTE(source_allocation=ngc.SourceAllocation.OPTIMIZED))
r = ngc.evaluate_demand(graph, 0, 3, ngc.EvaluationOptions(te, outputs=ngc.OutputOptions(min_cut=True)))
print(r.call_delivered, r.marginal_cost_curve, r.min_cut)  # 15 [(2, 10), (4, 5)] [0, 1]

# ECMP sends 12 on the only shortest path; its first link (edge 0) carries 10 and drops 2.
ecmp = ngc.FlowPolicyConfig(routing=ngc.ShortestPaths(), admission=ngc.Admission.LOSSY)
r = ngc.FlowPolicy(ngc.FlowLedger(graph), 0, 3, ecmp).place(12)
print(r.call_delivered, r.call_edge_drops)  # 10 [(0, 2)]

# T1: two tunnels of 6 for a demand of 12; the first takes the 10 path, the second finds no path
# with 6 unreserved. A reservation does not cap traffic: all 12 go into the tunnel that is up,
# and its path carries 10.
t1 = ngc.FlowPolicyConfig(routing=ngc.DistributedTE(ngc.TunnelCount(2)), admission=ngc.Admission.LOSSY)
r = ngc.FlowPolicy(ngc.FlowLedger(graph), 0, 3, t1).replace(12, design_demand=12)
print([t.status.name for t in r.program.tunnels], r.call_delivered)  # ['UP', 'INSUFFICIENT_BANDWIDTH'] 10
```

A `marginal_cost_curve` lists `(cost, units)` pairs: the first 10 units cost 2 each, the next 5
cost 4 each. With a finer resolution the same calls count fractions exactly:
`ngc.to_units(2.5, 2**-12)` is 10240 units.

## Installation

```bash
pip install netgraph-core
```

From source, `pip install .` builds the shipped configuration; `pip install -e .` is the editable
development build, which also carries private test hooks. The package is typed: it ships
`py.typed` and a type stub generated from the compiled module.

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

## Requirements

- **C++:** C++20 compiler (GCC 10+, Clang 12+, MSVC 2019+)
- **Python:** 3.11+
- **Build:** CMake 3.23+, scikit-build-core
- **Dependencies:** pybind11, NumPy

## License

[MIT License](https://github.com/networmix/NetGraph-Core/blob/main/LICENSE)
