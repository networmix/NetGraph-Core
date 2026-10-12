// Exact 128-bit integers for the test oracles: Int128 (signed) and UInt128 (unsigned). The
// portable classes wrap modulo 2^128; the built-in signed type must not overflow (that is
// undefined), so a test that wants a wrapped signed result forms it in UInt128 and converts back.
// Compilers with a built-in 128-bit integer (GCC, Clang) use it; others (MSVC) use the two small
// classes below, written for the tests alone and independent of the library's detail/int128.hpp,
// so an arithmetic error there cannot hide in an expected value. The CI job "C++ tests (windows,
// MSVC)" runs every test on the classes, and PortableInt128.PortableClassesMatchTheBuiltInType
// compares them with the built-in type wherever both exist. Defining
// NETGRAPH_CORE_TEST_PORTABLE_INT128 selects the classes on any compiler; no make target or CI job
// sets it. To run every test on them locally, from the repository root (a fresh build directory:
// CMake reads CXXFLAGS only when it first configures one):
//   CXXFLAGS=-DNETGRAPH_CORE_TEST_PORTABLE_INT128 cmake -S . -B build/cpp-tests-portable -DNETGRAPH_CORE_BUILD_TESTS=ON -DNETGRAPH_CORE_PRIVATE_HOOKS=ON
//   cmake --build build/cpp-tests-portable --parallel && ctest --test-dir build/cpp-tests-portable --no-tests=error
#pragma once

#include <compare>
#include <concepts>
#include <cstdint>
#include <cstdlib>
#include <type_traits>

namespace netgraph::core::test {

// Two 64-bit words. Every operation is modulo 2^128; division by zero aborts.
class PortableUInt128 {
 public:
  using U = PortableUInt128;

  constexpr PortableUInt128() noexcept = default;
  // Implicit, like a built-in integral conversion: a negative value wraps modulo 2^128.
  template <std::integral T>
  constexpr PortableUInt128(T v) noexcept : lo_(static_cast<std::uint64_t>(v)) {
    if constexpr (std::is_signed_v<T>) hi_ = v < 0 ? ~std::uint64_t {0} : 0;
  }
  static constexpr U of_words(std::uint64_t hi, std::uint64_t lo) noexcept {
    U r;
    r.hi_ = hi;
    r.lo_ = lo;
    return r;
  }
  [[nodiscard]] constexpr std::uint64_t high() const noexcept { return hi_; }
  [[nodiscard]] constexpr std::uint64_t low() const noexcept { return lo_; }
  // The low bits, like a built-in narrowing conversion.
  template <std::integral T>
  constexpr explicit operator T() const noexcept {
    return static_cast<T>(lo_);
  }

  friend constexpr bool operator==(const U&, const U&) noexcept = default;
  friend constexpr std::strong_ordering operator<=>(const U& a, const U& b) noexcept {
    if (a.hi_ != b.hi_) return a.hi_ <=> b.hi_;
    return a.lo_ <=> b.lo_;
  }

  friend constexpr U operator+(U a, U b) noexcept {
    const std::uint64_t lo = a.lo_ + b.lo_;
    return of_words(a.hi_ + b.hi_ + (lo < a.lo_ ? 1 : 0), lo);
  }
  friend constexpr U operator-(U a, U b) noexcept {
    return of_words(a.hi_ - b.hi_ - (a.lo_ < b.lo_ ? 1 : 0), a.lo_ - b.lo_);
  }
  friend constexpr U operator*(U a, U b) noexcept {
    // The full product of the low words from 32-bit halves; the products with a high word only
    // reach the high word (modulo 2^128).
    const std::uint64_t mask = 0xffffffffU;
    const std::uint64_t a0 = a.lo_ & mask, a1 = a.lo_ >> 32;
    const std::uint64_t b0 = b.lo_ & mask, b1 = b.lo_ >> 32;
    const std::uint64_t p00 = a0 * b0, p01 = a0 * b1, p10 = a1 * b0, p11 = a1 * b1;
    const std::uint64_t middle = (p00 >> 32) + (p01 & mask) + (p10 & mask);
    const std::uint64_t lo = (middle << 32) | (p00 & mask);
    const std::uint64_t hi = p11 + (p01 >> 32) + (p10 >> 32) + (middle >> 32) + a.hi_ * b.lo_ + a.lo_ * b.hi_;
    return of_words(hi, lo);
  }
  // Long division, one quotient bit per step: n = q * d + r with r < d.
  static constexpr void divide(U n, U d, U& q, U& r) noexcept {
    if (d == U {}) std::abort();
    q = U {};
    r = U {};
    for (int i = 127; i >= 0; --i) {
      r = (r << 1) | ((n >> i) & U {1});
      if (r >= d) {
        r = r - d;
        q = q | (U {1} << i);
      }
    }
  }
  friend constexpr U operator/(U a, U b) noexcept {
    U q, r;
    divide(a, b, q, r);
    return q;
  }
  friend constexpr U operator%(U a, U b) noexcept {
    U q, r;
    divide(a, b, q, r);
    return r;
  }
  // Shifts by 0 to 127 bits.
  friend constexpr U operator<<(U a, int s) noexcept {
    if (s == 0) return a;
    if (s >= 64) return of_words(a.lo_ << (s - 64), 0);
    return of_words((a.hi_ << s) | (a.lo_ >> (64 - s)), a.lo_ << s);
  }
  friend constexpr U operator>>(U a, int s) noexcept {
    if (s == 0) return a;
    if (s >= 64) return of_words(0, a.hi_ >> (s - 64));
    return of_words(a.hi_ >> s, (a.lo_ >> s) | (a.hi_ << (64 - s)));
  }
  friend constexpr U operator|(U a, U b) noexcept { return of_words(a.hi_ | b.hi_, a.lo_ | b.lo_); }
  friend constexpr U operator&(U a, U b) noexcept { return of_words(a.hi_ & b.hi_, a.lo_ & b.lo_); }
  constexpr U& operator+=(U b) noexcept { return *this = *this + b; }
  constexpr U& operator-=(U b) noexcept { return *this = *this - b; }
  constexpr U& operator*=(U b) noexcept { return *this = *this * b; }
  constexpr U& operator/=(U b) noexcept { return *this = *this / b; }

 private:
  std::uint64_t hi_ {0};
  std::uint64_t lo_ {0};
};

// The two's-complement reading of a PortableUInt128 bit pattern: signed comparison and division
// (truncating toward zero, like the built-in type); the other operations are those of the bits.
class PortableInt128 {
 public:
  using I = PortableInt128;
  using U = PortableUInt128;

  constexpr PortableInt128() noexcept = default;
  // Implicit, like a built-in integral conversion: the value is kept.
  template <std::integral T>
  constexpr PortableInt128(T v) noexcept : bits_(v) {}
  constexpr explicit PortableInt128(U bits) noexcept : bits_(bits) {}
  constexpr explicit operator U() const noexcept { return bits_; }
  template <std::integral T>
  constexpr explicit operator T() const noexcept {
    return static_cast<T>(bits_);
  }

  friend constexpr bool operator==(const I&, const I&) noexcept = default;
  friend constexpr std::strong_ordering operator<=>(const I& a, const I& b) noexcept {
    const auto ah = static_cast<std::int64_t>(a.bits_.high());
    const auto bh = static_cast<std::int64_t>(b.bits_.high());
    if (ah != bh) return ah <=> bh;
    return a.bits_.low() <=> b.bits_.low();
  }

  friend constexpr I operator+(I a, I b) noexcept { return I(a.bits_ + b.bits_); }
  friend constexpr I operator-(I a, I b) noexcept { return I(a.bits_ - b.bits_); }
  friend constexpr I operator-(I a) noexcept { return I(U {} - a.bits_); }
  friend constexpr I operator*(I a, I b) noexcept { return I(a.bits_ * b.bits_); }
  friend constexpr I operator/(I a, I b) noexcept {
    const I q(magnitude(a) / magnitude(b));
    return (a < 0) != (b < 0) ? -q : q;
  }
  friend constexpr I operator%(I a, I b) noexcept {
    const I r(magnitude(a) % magnitude(b));
    return a < 0 ? -r : r;
  }
  constexpr I& operator+=(I b) noexcept { return *this = *this + b; }
  constexpr I& operator-=(I b) noexcept { return *this = *this - b; }
  constexpr I& operator*=(I b) noexcept { return *this = *this * b; }

 private:
  static constexpr U magnitude(I x) noexcept { return x < 0 ? (-x).bits_ : x.bits_; }
  U bits_ {};
};

#if defined(__SIZEOF_INT128__) && !defined(NETGRAPH_CORE_TEST_PORTABLE_INT128)
using Int128 = __int128;
using UInt128 = unsigned __int128;
#else
using Int128 = PortableInt128;
using UInt128 = PortableUInt128;
#endif

// hi * 2^64 + lo.
constexpr Int128 int128_of(std::int64_t hi, std::uint64_t lo) noexcept {
  return static_cast<Int128>((static_cast<UInt128>(static_cast<std::uint64_t>(hi)) << 64) | UInt128 {lo});
}

// 2^127 - 1.
constexpr Int128 int128_max() noexcept { return static_cast<Int128>((UInt128 {1} << 127) - UInt128 {1}); }

} // namespace netgraph::core::test
