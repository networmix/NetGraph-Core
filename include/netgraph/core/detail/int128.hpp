/* Portable unsigned 128-bit products and comparisons for exact rational arithmetic (internal).

   Used by the highest-averages apportionment (quotients w / d compared by cross products, lower
   quotas floor(I * w / W)) and by the Headroom priority (residual / capacity compared by cross
   products). Every operand is below 2^63, so every product is below 2^126 and every quotient
   that is formed fits 64 bits. Compilers with __int128 use it; others use 32-bit limbs.
*/
#pragma once

#include <cstdint>

namespace netgraph::core::detail {

struct U128 {
  std::uint64_t hi {0};
  std::uint64_t lo {0};
};

inline bool operator==(U128 a, U128 b) noexcept { return a.hi == b.hi && a.lo == b.lo; }
inline bool operator<(U128 a, U128 b) noexcept { return a.hi != b.hi ? a.hi < b.hi : a.lo < b.lo; }
inline bool operator>(U128 a, U128 b) noexcept { return b < a; }

#if defined(__SIZEOF_INT128__)

inline U128 mul_u64(std::uint64_t a, std::uint64_t b) noexcept {
  const unsigned __int128 p = static_cast<unsigned __int128>(a) * b;
  return {static_cast<std::uint64_t>(p >> 64), static_cast<std::uint64_t>(p)};
}

// floor(n / d) for d > 0 and a quotient below 2^64 (n.hi < d).
inline std::uint64_t div_u128(U128 n, std::uint64_t d) noexcept {
  const unsigned __int128 x = (static_cast<unsigned __int128>(n.hi) << 64) | n.lo;
  return static_cast<std::uint64_t>(x / d);
}

#else

inline U128 mul_u64(std::uint64_t a, std::uint64_t b) noexcept {
  const std::uint64_t a0 = a & 0xffffffffULL, a1 = a >> 32;
  const std::uint64_t b0 = b & 0xffffffffULL, b1 = b >> 32;
  const std::uint64_t p00 = a0 * b0, p01 = a0 * b1, p10 = a1 * b0, p11 = a1 * b1;
  const std::uint64_t mid = (p00 >> 32) + (p01 & 0xffffffffULL) + (p10 & 0xffffffffULL);
  U128 r;
  r.lo = (mid << 32) | (p00 & 0xffffffffULL);
  r.hi = p11 + (p01 >> 32) + (p10 >> 32) + (mid >> 32);
  return r;
}

// floor(n / d) for d > 0 and a quotient below 2^64 (n.hi < d): restoring division, one
// quotient bit per step.
inline std::uint64_t div_u128(U128 n, std::uint64_t d) noexcept {
  std::uint64_t rem = n.hi;
  std::uint64_t q = 0;
  for (int i = 63; i >= 0; --i) {
    const bool carry = (rem >> 63) != 0;
    rem = (rem << 1) | ((n.lo >> i) & 1U);
    if (carry || rem >= d) {
      rem -= d;
      q |= std::uint64_t{1} << i;
    }
  }
  return q;
}

#endif

} // namespace netgraph::core::detail
