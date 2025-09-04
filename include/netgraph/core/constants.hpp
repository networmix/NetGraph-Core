/* Numeric thresholds shared by SPF and flow placement. */
#pragma once

#include "netgraph/core/types.hpp"

namespace netgraph::core {

// Global numeric thresholds used across SPF and flow placement
inline constexpr Cap  kMinCap  = static_cast<Cap>(1.0 / 4096.0);   // minimum effective remaining capacity
inline constexpr Flow kMinFlow = static_cast<Flow>(1.0 / 4096.0);  // minimum meaningful flow when augmenting
inline constexpr double kEpsilon = 1e-12;                           // numeric clamp

} // namespace netgraph::core
