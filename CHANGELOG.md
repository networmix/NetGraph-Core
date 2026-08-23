# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.8.0] - 2026-08-22

### Added

- **Flow Policy**: `set_static_paths(src, dst, paths)` pins a demand to explicit path bundles (MPLS-style routing), finishing a feature that was ported but never bound and whose port bound every flow to the first bundle. One flow per usable bundle, each bound to its own bundle in supply order; bundles are validated against the graph and pruned against the policy's masks, and a bundle with no surviving `src->dst` walk is down (no reroute). `EqualBalanced` spreads over the up bundles only. With static paths the policy neither grows its flow set nor reoptimizes; `max_path_cost`/`-factor`, `min_flow_count` and `reoptimize_flows_on_each_placement` are inert.
- **Shortest Paths**: `PredDAG.from_edges(graph, edges)` builds a single-path `PredDAG` from a contiguous edge-id sequence (`make_path_dag` in C++) -- the missing constructor for operator-defined explicit paths. Extracting it also deduplicates the two identical path-to-DAG conversion blocks in `k_shortest_paths` (outputs unchanged).

### Fixed

- **Max-Flow**: `calc_max_flow` could return less than the true maximum, since the tier loop augments only along forward SPF DAGs and never cancels an earlier placement (the reported `min_cut` then contradicted `total_flow`). Added a residual completion phase with reverse arcs for `Proportional` + `require_capacity` + `shortest_path=false`; results increase where the previous value was suboptimal.
- **Shortest Paths**: Zero-cost edges produced a cyclic predecessor DAG, making `EqualBalanced` placement return 0 flow and `resolve_to_paths` hang. Equal-cost predecessors are now recorded only while a node is unsettled, and `resolve_to_paths` guards against cycles in caller-supplied DAGs.
- **K-Shortest Paths**: Spur enumeration materialized every equal-cost path, exponential in ECMP fan-out (a 22-stage ladder needed ~19 s and ~3.8 GB for `k=3`). Enumeration is now bounded by the number of candidates still acceptable; tie-breaking among equal-cost paths may differ, but path counts and costs are unchanged.
- **Flow Policy**: `place_demand` silently routed a second `(src, dst)` pair over the first pair's paths. It now raises `invalid_argument`; use one policy per demand or call `remove_demand()` first.
- **Graph Construction**: `from_arrays` now rejects a total edge cost at or above 2^62, which overflows the int64 path arithmetic in SPF and silently corrupts results.
- **Flow State**: `compute_min_cut` and max-flow reachability derived placed flow from `capacity - residual`, overstating it for a `FlowState` built with a custom `residual_init`. Both now use `edge_flow`.
- **K-Shortest Paths**: `k_shortest_paths` read `paths.back()` on an empty vector when `max_cost_factor < 1.0` put the cost ceiling below the shortest path, crashing for `k > 1` (the `k == 1` early return masked it). All `k` now return no paths for such a factor.
- **Python Bindings**: `batch_max_flow` rejected a valid int32 `pairs` array on Windows, because the dtype check compared buffer format strings and NumPy spells int32 as `NPY_LONG` on LLP64 but `NPY_INT` on LP64. It now compares dtype equivalence, and rejects non-contiguous `pairs` rather than reading them as if packed.
- **Build**: `NETGRAPH_CORE_SANITIZE` passed its flags as one quoted string, so sanitizer builds never compiled, and the test target was missing them entirely. `make sanitize-test` also no longer sets `detect_leaks=1` on macOS, where it aborts.
- **Build**: Coverage builds now use `-fprofile-update=atomic` on GCC; the default non-atomic counters corrupt under the threads used by `batch_max_flow` and `sensitivity_analysis`.

### Changed

- **Shortest Paths**: Replaced nested per-node predecessor vectors with flat intrusive lists, removing the allocation churn that dominated the hot path (~3-4x faster; output unchanged).
- **Flow Placement**: `place_on_dag` now rebuilds its `(parent, child)` edge groups into a reused arena instead of nested per-node vectors, which was ~70% of its runtime (~2-2.4x faster; output unchanged).
- **Max-Flow**: Parallelized `batch_max_flow` across source/destination pairs using `std::async`; workers claim pairs from a shared counter, so a batch whose costs are unevenly distributed still parallelizes. Thread count from `NGRAPH_CORE_BATCH_THREADS` env or hardware concurrency.
- **Flow Policy**: `place_demand` computes one SPF when seeding initial flows instead of repeating an identical one per flow; seeding cost no longer scales with `min_flow_count`.
- **Flow Policy**: The `EqualBalanced` rebalance recursion (`place_demand` -> `rebalance_demand` -> `place_demand`) is now an iterative loop with numerically identical results (returned leftover may differ by ~1 ulp from floating-point summation order); with many pinned bundles of heterogeneous capacity the recursion depth grew like `U * ln(imbalance/kMinFlow)`, a stack-overflow risk on worker threads.
- **Python Bindings**: A wrong-typed `graph`/`algorithms` argument now raises `TypeError` instead of an opaque `RuntimeError: Unable to cast ... to C++ type '?'`, `Algorithms.spf` raises `TypeError` for a wrong `residual` length (matching every other length check), `Algorithms.ksp` validates `dtype` before running, and `batch_max_flow` rejects out-of-range node ids like the single-pair entry points.
- **Python Typing**: Corrected `_docs.py` (pybind11 enums and classes are not `enum.Enum`/`dataclass`, `MinCut.edges` is an `int32` array, removed a `Path` class that never existed, and added the constructors for `FlowIndex`, `FlowPolicyConfig` and the three enums) and widened the type-checking re-exports from 6 names to all 16, so `Algorithms`, `FlowPolicy`, `StrictMultiDiGraph` and friends are no longer `Unknown` to type checkers despite the shipped `py.typed`. A missing stub constructor is now worse than no stub, so a test asserts every stub constructor matches its binding.
- **Tests**: Several assertions could never fail and now do: a self-loop assertion sat inside `except Exception: pass`, a thread-safety test ended in `assert True` while its workers called `pytest.fail()` from a thread where it cannot fail the test, and two tautologies (`assert size >= 0`) were replaced with the properties actually under test.
- **Docs**: Corrected README/CONTRIBUTING (required CMake is 3.23 not 3.15, `make py-test` does not exist, `make test` does not collect coverage), the `sensitivity_analysis` description in README and `backend.hpp` (it measures flow *lost* on edge removal, not gain from relaxing capacity), the zero-copy and GIL claims, and header contracts for `from_arrays`, `PredDAG`, `FlowSummary.costs`, `calc_max_flow` and `FlowPolicy`.

### Removed

- **Dead code**: `.gcovr.cfg` (never read -- gcovr's default config name has no leading dot, and `make cov` passes every flag explicitly), the unreachable `FlowGraph` include and `expect_flow_conservation()` helper in the C++ test utilities, unused test helpers/fixtures, `_USE_MATH_DEFINES` (no `M_*` macro is used), and unused `<optional>`/`<stack>`/`<unordered_map>`/`<unordered_set>`/`<cmath>` includes across `src/` and `include/`.
- **Python Bindings**: `Algorithms.build_graph_from_arrays`. Nothing used it (no test, and not NetGraph), and the `Graph` handle it returned could not be passed to `FlowGraph` or `FlowPolicy`, so it only ever served the stateless algorithms while looking like a peer of `build_graph`. Build the graph with `StrictMultiDiGraph.from_arrays()` and pass it to `Algorithms.build_graph()`.
- **C++ API**: The `build_graph(std::shared_ptr<const StrictMultiDiGraph>)` overload on `Algorithms` and `Backend`, which existed solely to serve that binding. Removing it also drops a pure virtual that every `Backend` implementation had to provide. Callers wanting the handle to own the graph can construct it directly: `GraphHandle{my_shared_ptr}`.
- **Python Bindings**: The unreachable `Flow` class (bound from `FlowRecord`). No binding constructed or returned one, it was absent from `__all__`, and the name collided with the C++ alias `using Flow = double`.

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
