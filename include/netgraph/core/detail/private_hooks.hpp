/* The private hooks of the development build (internal). CMake option NETGRAPH_CORE_PRIVATE_HOOKS:
   ON for the editable development build and the test suites, OFF for an installed wheel. The
   hooks declared here are test controls and diagnostics the tests compare with their
   independent model of the policy; none of them is public API:
   a lowered creation cap and the transaction trace.
   Without the option neither is built: the engine reads the constant kMaxFlowsPerCall, and every
   trace write is a compile-time-dead branch (kTraced), so only empty trace values are built. */
#pragma once

#include "netgraph/core/detail/hooks_flag.hpp"

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

// FlowPolicy::place with the transaction's trace (a Trace: the call's intermediate values,
// keyed as the tests' model of the policy keys them). `staged`, when given, runs
// once the result and the trace are final and before anything is published (the commit, the
// epoch state): an exception from it aborts the call with nothing changed (the binding builds
// its Python objects there).
PlacementResult place_traced(FlowPolicy& policy, std::optional<Amount> volume, const OutputOptions& outputs, Trace& trace,
                             const std::function<void()>* staged = nullptr);

} // namespace netgraph::core::detail

#endif
