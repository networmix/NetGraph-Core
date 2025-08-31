#pragma once

#include <cstdint>

namespace netgraph::core {

// NetGraph uses signed 32-bit ids everywhere at the boundary and internally.
using NodeId = std::int32_t;
using EdgeId = std::int32_t;
using Cost   = double;
using Cap    = double;

enum class EdgeSelect {
  AllMinCost = 1,
  SingleMinCost = 2,
  AllMinCostWithCapRemaining = 3,
  AllAnyCostWithCapRemaining = 4,
  SingleMinCostWithCapRemaining = 5,
  SingleMinCostWithCapRemainingLoadFactored = 6,
  UserDefined = 99
};

enum class FlowPlacement {
  Proportional = 1,
  EqualBalanced = 2
};

} // namespace netgraph::core
