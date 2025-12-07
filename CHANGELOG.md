# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
