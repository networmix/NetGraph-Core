/* Overflow-safe helpers for Cost arithmetic. */
#pragma once

#include <cmath>
#include <limits>

#include "netgraph/core/types.hpp"

namespace netgraph::core {

[[nodiscard]] inline constexpr Cost cost_max() noexcept {
  return std::numeric_limits<Cost>::max();
}

// Largest finite cost value (keeps cost_max() reserved as "infinity"/unreachable).
[[nodiscard]] inline constexpr Cost cost_finite_max() noexcept {
  return static_cast<Cost>(std::numeric_limits<Cost>::max() - 1);
}

// Saturating add for path costs.
// - If either argument is "infinity" (cost_max), returns cost_max.
// - If finite addition would overflow, returns cost_finite_max.
[[nodiscard]] inline constexpr Cost saturating_cost_add(Cost a, Cost b) noexcept {
  if (a == cost_max() || b == cost_max()) return cost_max();
  if (a >= cost_finite_max() || b >= cost_finite_max()) return cost_finite_max();
  if (b > 0 && a > cost_finite_max() - b) return cost_finite_max();
  if (b < 0 && a < std::numeric_limits<Cost>::min() - b) return std::numeric_limits<Cost>::min();
  return static_cast<Cost>(a + b);
}

// Saturating conversion of base*factor into Cost to avoid UB from out-of-range
// float->int conversions.
[[nodiscard]] inline Cost saturating_cost_mul_factor(Cost base, double factor) noexcept {
  if (base <= 0 || factor <= 0.0) return static_cast<Cost>(0);
  const long double scaled = static_cast<long double>(base) * static_cast<long double>(factor);
  const long double hi = static_cast<long double>(cost_finite_max());
  if (!std::isfinite(static_cast<double>(scaled)) || scaled >= hi) return cost_finite_max();
  if (scaled <= 0.0L) return static_cast<Cost>(0);
  return static_cast<Cost>(scaled);
}

} // namespace netgraph::core
