/* The private hooks of the development build (internal). CMake option NETGRAPH_CORE_PRIVATE_HOOKS:
   ON for the editable development build and the test suites, OFF for an installed wheel. The
   hooks are test controls and prototype-parity diagnostics that the contract does not name:
   a lowered creation cap, an injected sensitivity model failure and the transaction trace.
   Without the option none of them exists: the engine reads the constant kMaxFlowsPerCall,
   sensitivity has no injection branch, and every trace block is a compile-time-dead branch
   (kTraced), so no trace code is emitted. */
#pragma once

#ifndef NETGRAPH_CORE_PRIVATE_HOOKS
#define NETGRAPH_CORE_PRIVATE_HOOKS 0
#endif

namespace netgraph::core::detail {

// Whether this build traces placement calls (the private hooks).
inline constexpr bool kTraced = NETGRAPH_CORE_PRIVATE_HOOKS != 0;

} // namespace netgraph::core::detail

#if NETGRAPH_CORE_PRIVATE_HOOKS

#include <atomic>
#include <cstdint>
#include <functional>
#include <optional>

#include "netgraph/core/detail/trace.hpp"
#include "netgraph/core/flow_policy.hpp"

namespace netgraph::core::detail {

// The call-wide creation cap every CapacityAware call reads (kMaxFlowsPerCall unless a test
// lowers it to reach the cap on small graphs).
std::atomic<std::int64_t>& max_flows_per_call() noexcept;

// The sensitivity candidate whose deletion raises a model exception (-1: none). No valid
// configuration raises under an edge deletion (every construction check is
// deletion-independent), so this is how the tests reach `unevaluable` (APPROACH 11).
std::atomic<std::int64_t>& sensitivity_fail_edge() noexcept;

// FlowPolicy::place with the transaction's trace (prototype 2's keys). `staged`, when given, runs
// once the result and the trace are final and before anything is published (the commit, the
// epoch state): an exception from it aborts the call with nothing changed (the binding builds
// its Python objects there).
PlacementResult place_traced(FlowPolicy& policy, std::optional<Amount> volume, const OutputOptions& outputs, Trace& trace,
                             const std::function<void()>* staged = nullptr);

} // namespace netgraph::core::detail

#endif
