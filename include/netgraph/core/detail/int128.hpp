/* Portable 128-bit values (internal).

   U128, unsigned, for exact rational arithmetic: the highest-averages apportionment (quotients
   w / d compared by cross products, lower quotas floor(I * w / W)), the Headroom priority
   (residual / capacity compared by cross products), the comparisons with a double factor read
   as the exact dyadic rational it is (the Stretch ceiling, the diminishing guard) and the
   cumulative counters (Counter128). Every product operand is below 2^63, so every product is
   below 2^126 and every quotient that is formed fits 64 bits. Compilers with __int128 use it;
   others use 32-bit limbs.

   I128, signed, for the deletion re-solve of the flow core (APPROACH 2.2): working potentials,
   distance keys and sums of tier maxima, which stay within (22 n + 15) C for n vertices and the
   cost sum C < 2^62. Only addition, subtraction and comparison, on two limbs everywhere.
*/
#pragma once

#include <bit>
#include <cmath>
#include <compare>
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
  // A dividend below 2^64 (the common lower quota) takes the 64-bit division: the same quotient
  // without the 128-bit division helper.
  if (n.hi == 0) return n.lo / d;
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

// The number of significant bits of v (0 for 0).
inline int bits(U128 v) noexcept { return v.hi != 0 ? 64 + std::bit_width(v.hi) : std::bit_width(v.lo); }

// v * 2^k for 0 <= k < 128; the caller keeps the result below 2^128 (bits(v) + k <= 128).
inline U128 shl(U128 v, int k) noexcept {
  if (k == 0) return v;
  if (k >= 64) return {v.lo << (k - 64), 0};
  return {(v.hi << k) | (v.lo >> (64 - k)), v.lo << k};
}

// A signed 128-bit integer in two's complement, value = hi * 2^64 + lo. The member order makes
// the defaulted comparison (hi signed, then lo unsigned) the numeric order.
struct I128 {
  std::int64_t hi {0};
  std::uint64_t lo {0};

  constexpr I128() noexcept = default;
  // Implicit, so int64 costs and amounts enter the arithmetic directly.
  constexpr I128(std::int64_t v) noexcept : hi(v < 0 ? -1 : 0), lo(static_cast<std::uint64_t>(v)) {}

  friend constexpr bool operator==(const I128&, const I128&) noexcept = default;
  friend constexpr std::strong_ordering operator<=>(const I128&, const I128&) noexcept = default;

  friend constexpr I128 operator+(I128 a, I128 b) noexcept {
    I128 r;
    r.lo = a.lo + b.lo;
    const std::uint64_t carry = r.lo < a.lo ? 1U : 0U;
    r.hi = static_cast<std::int64_t>(static_cast<std::uint64_t>(a.hi) + static_cast<std::uint64_t>(b.hi) + carry);
    return r;
  }
  friend constexpr I128 operator-(I128 a, I128 b) noexcept {
    I128 r;
    r.lo = a.lo - b.lo;
    const std::uint64_t borrow = a.lo < b.lo ? 1U : 0U;
    r.hi = static_cast<std::int64_t>(static_cast<std::uint64_t>(a.hi) - static_cast<std::uint64_t>(b.hi) - borrow);
    return r;
  }
  constexpr I128& operator+=(I128 b) noexcept { return *this = *this + b; }
};

// A positive finite double as the exact dyadic rational it is: f = m * 2^x with m odd and
// below 2^53.
struct Dyadic {
  std::uint64_t m {0};
  int x {0};
};
inline Dyadic dyadic(double f) noexcept {
  int e = 0;
  const double frac = std::frexp(f, &e);  // f = frac * 2^e, frac in [0.5, 1)
  Dyadic d {static_cast<std::uint64_t>(std::ldexp(frac, 53)), e - 53};
  const int tz = std::countr_zero(d.m);
  d.m >>= tz;
  d.x += tz;
  return d;
}

} // namespace netgraph::core::detail
