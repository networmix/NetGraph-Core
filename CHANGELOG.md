# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
