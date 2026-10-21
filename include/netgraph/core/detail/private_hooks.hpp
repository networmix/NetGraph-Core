/* The private hooks of the development build (internal). CMake option NETGRAPH_CORE_PRIVATE_HOOKS:
   ON for the editable development build and the test suites, OFF for an installed wheel. The
   hooks declared here are test controls and diagnostics the tests compare with their
   independent model of the policy; none of them is public API: the transaction trace.
   Without the option it is not built: every trace write is a compile-time-dead branch (kTraced),
   so only empty trace values are built. */
#pragma once

#include "netgraph/core/detail/hooks_flag.hpp"

namespace netgraph::core::detail {

// Whether this build traces placement calls (the private hooks).
inline constexpr bool kTraced = NETGRAPH_CORE_PRIVATE_HOOKS != 0;

} // namespace netgraph::core::detail

#if NETGRAPH_CORE_PRIVATE_HOOKS

#include <optional>

#include "netgraph/core/detail/trace.hpp"
#include "netgraph/core/flow_policy.hpp"

namespace netgraph::core::detail {

// FlowPolicy::place with a native trace. Python conversion is a separate read,
// after the operation has released the ledger gate.
PlacementResult place_traced(FlowPolicy& policy, std::optional<Amount> offer, const OutputOptions& outputs, Trace& trace);

} // namespace netgraph::core::detail

#endif
