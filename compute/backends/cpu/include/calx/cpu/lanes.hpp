// Vector lanes for the native path of the CPU backend: W lanes of an arithmetic type T, on
// GCC's (and Clang's) vector extensions, so that a kernel written once is compiled to SSE2,
// AVX2, AVX-512 or NEON, whichever the build allows (target.hpp, LANE_BYTES), and to correct
// scalar code where there is none.
//
//   using V = calx::cpu::Lanes<std::int16_t>;     // as many lanes as the build's widest vector
//   V h = V::load(p), e = V::splat(-3);
//   h = max(adds(h, e), V::zero());               // saturating add, then a lane-wise maximum
//   h = h.shift_up(0);                            // lane i takes lane i - 1; lane 0 takes 0
//
// Semantics, for every width (they are checked lane by lane against scalar code):
// - integer +, -, * wrap modulo 2^bits, computed on the unsigned type so that no signed
//   overflow is undefined; floating-point arithmetic is IEEE in each lane;
// - adds and subs saturate at the limits of T (integers of 8 and 16 bits use the
//   instructions PADDS/PADDUS/PSUBS/PSUBUS where the build has them, and the same results
//   from a generic formula elsewhere);
// - comparisons give a Mask, all ones in a lane where they hold;
// - shift_up(x) moves every lane up by one (towards higher addresses), and puts x in lane 0:
//   the shift of Farrar's striped Smith–Waterman (Farrar 2007, Bioinformatics 23, 156–161);
// - hmax and hmin reduce across the lanes; any and all test a mask.
//
// The vectors are values of the extension's types, so an expression of Lanes compiles to the
// instructions themselves; there is no dispatch at run time. A program that wants the host's
// widest lanes is built with -march=native (or the flags of its target).

#pragma once

#include <calx/cpu/target.hpp>

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <type_traits>
#include <utility>

#if defined(__SSE2__)
#include <immintrin.h>
#endif

// Passing a 32- or 64-byte vector by value changes the ABI between builds with and without
// AVX; within one program it is consistent, which is all these inline functions need.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpsabi"

namespace calx::cpu {

namespace detail {
template <std::size_t N> struct signed_of;
template <> struct signed_of<1> { using type = std::int8_t; };
template <> struct signed_of<2> { using type = std::int16_t; };
template <> struct signed_of<4> { using type = std::int32_t; };
template <> struct signed_of<8> { using type = std::int64_t; };
}  // namespace detail

template <class T, std::size_t W>
struct Mask {
  using Int = typename detail::signed_of<sizeof(T)>::type;
  typedef Int Vec __attribute__((vector_size(W * sizeof(T))));
  Vec m;

  bool operator[](std::size_t i) const { return m[i] != 0; }
  friend Mask operator&(Mask a, Mask b) { return {a.m & b.m}; }
  friend Mask operator|(Mask a, Mask b) { return {a.m | b.m}; }
  friend Mask operator^(Mask a, Mask b) { return {a.m ^ b.m}; }
  Mask operator~() const { return {~m}; }

  // Whether any lane is set.
  bool any() const {
#if defined(__AVX512BW__)
    if constexpr (sizeof(Vec) == 64) { const __m512i x = (__m512i)m; return _mm512_test_epi8_mask(x, x) != 0; }
#endif
#if defined(__AVX2__)
    if constexpr (sizeof(Vec) == 32) return _mm256_movemask_epi8((__m256i)m) != 0;
#endif
#if defined(__SSE2__)
    if constexpr (sizeof(Vec) == 16) return _mm_movemask_epi8((__m128i)m) != 0;
#endif
    if constexpr (sizeof(Vec) % 8 == 0) {
      std::uint64_t w[sizeof(Vec) / 8];
      std::memcpy(w, &m, sizeof(Vec));
      std::uint64_t r = 0;
      for (std::uint64_t x : w) r |= x;
      return r != 0;
    } else {
      for (std::size_t i = 0; i < W; i++)
        if (m[i]) return true;
      return false;
    }
  }
  bool all() const { return !(~*this).any(); }
  // Lane i in bit i (W at most 64).
  std::uint64_t bits() const {
    static_assert(W <= 64);
    std::uint64_t b = 0;
    for (std::size_t i = 0; i < W; i++) b |= static_cast<std::uint64_t>(m[i] != 0) << i;
    return b;
  }
};

template <class T, std::size_t W = LANE_BYTES / sizeof(T)>
struct Lanes {
  static_assert(std::is_arithmetic_v<T> && !std::is_same_v<T, bool>, "lanes hold numbers");
  static_assert(W >= 1 && (W & (W - 1)) == 0, "the number of lanes is a power of two");

  static constexpr std::size_t width = W;
  static constexpr std::size_t bytes = W * sizeof(T);
  using value_type = T;
  using MaskType = Mask<T, W>;
  using Int = typename detail::signed_of<sizeof(T)>::type;
  using UInt = std::make_unsigned_t<Int>;
  typedef T Vec __attribute__((vector_size(bytes)));
  typedef UInt Bits __attribute__((vector_size(bytes)));

  Vec v;

  // ---- making and moving -------------------------------------------------------------------
  static Lanes zero() { return {Vec{}}; }
  static Lanes splat(T x) { return {Vec{} + x}; }
  // Lanes start, start + 1, ...
  static Lanes iota(T start = 0) {
    Lanes r{};
    for (std::size_t i = 0; i < W; i++) r.v[i] = static_cast<T>(start + static_cast<T>(i));
    return r;
  }
  // From W values at p, which need not be aligned.
  static Lanes load(const T* p) {
    Lanes r;
    std::memcpy(&r.v, p, bytes);
    return r;
  }
  void store(T* p) const { std::memcpy(p, &v, bytes); }
  T operator[](std::size_t i) const { return v[i]; }
  void set(std::size_t i, T x) { v[i] = x; }

  // ---- arithmetic (integers wrap) --------------------------------------------------------
  friend Lanes operator+(Lanes a, Lanes b) {
    if constexpr (std::is_integral_v<T>) return {(Vec)((Bits)a.v + (Bits)b.v)};
    else return {a.v + b.v};
  }
  friend Lanes operator-(Lanes a, Lanes b) {
    if constexpr (std::is_integral_v<T>) return {(Vec)((Bits)a.v - (Bits)b.v)};
    else return {a.v - b.v};
  }
  friend Lanes operator*(Lanes a, Lanes b) {
    if constexpr (std::is_integral_v<T>) return {(Vec)((Bits)a.v * (Bits)b.v)};
    else return {a.v * b.v};
  }
  friend Lanes operator&(Lanes a, Lanes b) requires std::is_integral_v<T> { return {a.v & b.v}; }
  friend Lanes operator|(Lanes a, Lanes b) requires std::is_integral_v<T> { return {a.v | b.v}; }
  friend Lanes operator^(Lanes a, Lanes b) requires std::is_integral_v<T> { return {a.v ^ b.v}; }
  friend Lanes operator<<(Lanes a, int n) requires std::is_integral_v<T> { return {(Vec)((Bits)a.v << n)}; }
  // Arithmetic for signed lanes, logical for unsigned.
  friend Lanes operator>>(Lanes a, int n) requires std::is_integral_v<T> { return {a.v >> n}; }
  Lanes& operator+=(Lanes b) { return *this = *this + b; }
  Lanes& operator-=(Lanes b) { return *this = *this - b; }

  // ---- comparisons and choices ------------------------------------------------------------
  friend MaskType operator>(Lanes a, Lanes b) { return {a.v > b.v}; }
  friend MaskType operator<(Lanes a, Lanes b) { return {a.v < b.v}; }
  friend MaskType operator>=(Lanes a, Lanes b) { return {a.v >= b.v}; }
  friend MaskType operator<=(Lanes a, Lanes b) { return {a.v <= b.v}; }
  friend MaskType operator==(Lanes a, Lanes b) { return {a.v == b.v}; }
  friend MaskType operator!=(Lanes a, Lanes b) { return {a.v != b.v}; }

  // a where the mask is set, b elsewhere.
  friend Lanes select(MaskType m, Lanes a, Lanes b) { return {m.m ? a.v : b.v}; }
  friend Lanes max(Lanes a, Lanes b) { return {a.v > b.v ? a.v : b.v}; }
  friend Lanes min(Lanes a, Lanes b) { return {a.v < b.v ? a.v : b.v}; }

  // ---- saturating arithmetic ------------------------------------------------------------------
  friend Lanes adds(Lanes a, Lanes b) requires std::is_integral_v<T> {
#if defined(__SSE2__)
    if constexpr (sizeof(T) <= 2) {
      if constexpr (bytes == 16) return {(Vec)native_adds((__m128i)a.v, (__m128i)b.v)};
#if defined(__AVX2__)
      if constexpr (bytes == 32) return {(Vec)native_adds((__m256i)a.v, (__m256i)b.v)};
#endif
#if defined(__AVX512BW__)
      if constexpr (bytes == 64) return {(Vec)native_adds((__m512i)a.v, (__m512i)b.v)};
#endif
    }
#endif
    const Lanes r = a + b;
    if constexpr (std::is_unsigned_v<T>) {
      return select(r < a, splat(std::numeric_limits<T>::max()), r);
    } else {
      // Overflow when a and b have one sign and the sum the other (Hacker's Delight, §2-13).
      const MaskType over = ((a ^ r) & (b ^ r)) < zero();
      return select(over, select(a < zero(), splat(std::numeric_limits<T>::min()), splat(std::numeric_limits<T>::max())), r);
    }
  }
  friend Lanes subs(Lanes a, Lanes b) requires std::is_integral_v<T> {
#if defined(__SSE2__)
    if constexpr (sizeof(T) <= 2) {
      if constexpr (bytes == 16) return {(Vec)native_subs((__m128i)a.v, (__m128i)b.v)};
#if defined(__AVX2__)
      if constexpr (bytes == 32) return {(Vec)native_subs((__m256i)a.v, (__m256i)b.v)};
#endif
#if defined(__AVX512BW__)
      if constexpr (bytes == 64) return {(Vec)native_subs((__m512i)a.v, (__m512i)b.v)};
#endif
    }
#endif
    const Lanes r = a - b;
    if constexpr (std::is_unsigned_v<T>) {
      return select(a > b, r, zero());
    } else {
      const MaskType over = ((a ^ b) & (a ^ r)) < zero();
      return select(over, select(a < zero(), splat(std::numeric_limits<T>::min()), splat(std::numeric_limits<T>::max())), r);
    }
  }

  // ---- across lanes ------------------------------------------------------------------------
  // Lane i takes lane i - 1, and lane 0 takes x.
  Lanes shift_up(T x) const { return {__builtin_shuffle(v, splat(x).v, up_index(std::make_index_sequence<W>{}))}; }
  // Lane i takes lane i + 1, and lane W - 1 takes x.
  Lanes shift_down(T x) const { return {__builtin_shuffle(v, splat(x).v, down_index(std::make_index_sequence<W>{}))}; }
  T hmax() const {
    T m = v[0];
    for (std::size_t i = 1; i < W; i++) m = v[i] > m ? v[i] : m;
    return m;
  }
  T hmin() const {
    T m = v[0];
    for (std::size_t i = 1; i < W; i++) m = v[i] < m ? v[i] : m;
    return m;
  }
  // The sum of the lanes, in T (integers wrap).
  T hsum() const {
    T s = 0;
    for (std::size_t i = 0; i < W; i++) s = static_cast<T>(s + v[i]);
    return s;
  }

 private:
  template <std::size_t... K>
  static constexpr Bits up_index(std::index_sequence<K...>) { return Bits{static_cast<UInt>(K == 0 ? W : K - 1)...}; }
  template <std::size_t... K>
  static constexpr Bits down_index(std::index_sequence<K...>) { return Bits{static_cast<UInt>(K + 1)...}; }

#if defined(__SSE2__)
  static __m128i native_adds(__m128i a, __m128i b) {
    if constexpr (std::is_same_v<T, std::uint8_t>) return _mm_adds_epu8(a, b);
    else if constexpr (std::is_same_v<T, std::int8_t>) return _mm_adds_epi8(a, b);
    else if constexpr (std::is_same_v<T, std::uint16_t>) return _mm_adds_epu16(a, b);
    else return _mm_adds_epi16(a, b);
  }
  static __m128i native_subs(__m128i a, __m128i b) {
    if constexpr (std::is_same_v<T, std::uint8_t>) return _mm_subs_epu8(a, b);
    else if constexpr (std::is_same_v<T, std::int8_t>) return _mm_subs_epi8(a, b);
    else if constexpr (std::is_same_v<T, std::uint16_t>) return _mm_subs_epu16(a, b);
    else return _mm_subs_epi16(a, b);
  }
#endif
#if defined(__AVX2__)
  static __m256i native_adds(__m256i a, __m256i b) {
    if constexpr (std::is_same_v<T, std::uint8_t>) return _mm256_adds_epu8(a, b);
    else if constexpr (std::is_same_v<T, std::int8_t>) return _mm256_adds_epi8(a, b);
    else if constexpr (std::is_same_v<T, std::uint16_t>) return _mm256_adds_epu16(a, b);
    else return _mm256_adds_epi16(a, b);
  }
  static __m256i native_subs(__m256i a, __m256i b) {
    if constexpr (std::is_same_v<T, std::uint8_t>) return _mm256_subs_epu8(a, b);
    else if constexpr (std::is_same_v<T, std::int8_t>) return _mm256_subs_epi8(a, b);
    else if constexpr (std::is_same_v<T, std::uint16_t>) return _mm256_subs_epu16(a, b);
    else return _mm256_subs_epi16(a, b);
  }
#endif
#if defined(__AVX512BW__)
  static __m512i native_adds(__m512i a, __m512i b) {
    if constexpr (std::is_same_v<T, std::uint8_t>) return _mm512_adds_epu8(a, b);
    else if constexpr (std::is_same_v<T, std::int8_t>) return _mm512_adds_epi8(a, b);
    else if constexpr (std::is_same_v<T, std::uint16_t>) return _mm512_adds_epu16(a, b);
    else return _mm512_adds_epi16(a, b);
  }
  static __m512i native_subs(__m512i a, __m512i b) {
    if constexpr (std::is_same_v<T, std::uint8_t>) return _mm512_subs_epu8(a, b);
    else if constexpr (std::is_same_v<T, std::int8_t>) return _mm512_subs_epi8(a, b);
    else if constexpr (std::is_same_v<T, std::uint16_t>) return _mm512_subs_epu16(a, b);
    else return _mm512_subs_epi16(a, b);
  }
#endif
};

}  // namespace calx::cpu

#pragma GCC diagnostic pop
