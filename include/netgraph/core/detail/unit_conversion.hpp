// Numeric kernels after validation of the resolution and physical argument.
#pragma once

#include <cmath>
#include "netgraph/core/types.hpp"

namespace netgraph::core::detail {

// q = 2^k. The caller has checked finite, nonnegative x and the resolution.
// Keep the rounded value as a double until its caller has checked the unit bound.
inline double rounded_units(double x, int k) noexcept {
  // x / q is exact (q is a power of two) unless it overflows to infinity or falls below the
  // normal range (below 2^-1022 it may lose low bits; it then rounds to 0 units either way).
  const double y = std::ldexp(x, -k);
  constexpr double kTwo52 = 4503599627370496.0;  // 2^52: doubles at or above it are integers
  if (y >= kTwo52) return y;
  const double f = std::floor(y);
  return (y - f >= 0.5) ? f + 1.0 : f;  // y - f is exact below 2^52: round half up
}

inline double scaled_units(Amount x, int k) noexcept {
  return std::ldexp(static_cast<double>(x), k);
}

} // namespace netgraph::core::detail
