/* Portable 128-bit arithmetic (internal).

   U128 holds unsigned values for exact allocation shares, priorities and factor comparisons
   (Counter128 keeps its own two words). Compare fractions a/b and c/d by comparing a*d and
   c*b; compute a lower quota as floor(I*w/W), where I is the total, w one weight and W the sum
   of weights.
   A cost-stretch factor is an exact fraction n/d.
   Product operands stay below 2^63, products below 2^126, and quotients fit 64 bits.
   Use __int128 when available; otherwise split products into 32-bit parts (limbs).

   I128 holds signed costs for deletion repairs: vertex cost offsets (potentials), distance
   keys and sums of each shortest-path search's largest distance (tier maxima). These stay
   within (22 n + 15) C for n vertices and cost sum C < 2^62. Only addition, subtraction and
   comparison are needed, always using two 64-bit parts.
*/
#pragma once

#include <bit>
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

// floor(n / d), with d > 0 and quotient < 2^64 (n.hi < d). Build one quotient bit per step:
// shift in the next dividend bit, then subtract d if the current remainder is large enough.
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

// a + b and a - b (a >= b), exactly; the caller keeps a + b below 2^128.
inline U128 add(U128 a, U128 b) noexcept {
  U128 r {a.hi + b.hi, a.lo + b.lo};
  if (r.lo < a.lo) ++r.hi;
  return r;
}
inline U128 sub(U128 a, U128 b) noexcept {
  U128 r {a.hi - b.hi, a.lo - b.lo};
  if (a.lo < b.lo) --r.hi;
  return r;
}

// The number of significant bits of v (0 for 0).
inline int bits(U128 v) noexcept { return v.hi != 0 ? 64 + std::bit_width(v.hi) : std::bit_width(v.lo); }

// v * 2^k for 0 <= k < 128; the caller keeps the result below 2^128 (bits(v) + k <= 128).
inline U128 shl(U128 v, int k) noexcept {
  if (k == 0) return v;
  if (k >= 64) return {v.lo << (k - 64), 0};
  return {(v.hi << k) | (v.lo >> (64 - k)), v.lo << k};
}

// floor(n / d) for d > 0 and a quotient below 2^64 (a lower quota floor(I*w/W) with w <= W: at
// most I). A divisor below 2^64 takes div_u128 above; a larger one (a weight sum of 2^64 or more)
// the long division, one quotient bit per step from the highest that can be set.
inline std::uint64_t div_u128(U128 n, U128 d) noexcept {
  if (d.hi == 0) return div_u128(n, d.lo);
  std::uint64_t q = 0;
  for (int s = bits(n) - bits(d); s >= 0; --s) {
    const U128 ds = shl(d, s);
    if (!(n < ds)) {
      n = sub(n, ds);
      q |= std::uint64_t {1} << s;
    }
  }
  return q;
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

// cost <= (n / d) * base, exactly (the CostStretch cost ceiling): cost * d <= n * base, for cost and
// base in [0, 2^63) and a reduced factor n / d >= 1 (so n >= d >= 1). Each product is below 2^126.
inline bool stretch_le(std::int64_t cost, std::int64_t base, std::int64_t n, std::int64_t d) noexcept {
  const U128 lhs = mul_u64(static_cast<std::uint64_t>(cost), static_cast<std::uint64_t>(d));
  const U128 rhs = mul_u64(static_cast<std::uint64_t>(base), static_cast<std::uint64_t>(n));
  return !(rhs < lhs);
}

} // namespace netgraph::core::detail
