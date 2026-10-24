// Checks of the CPU backend's target description and lanes: the host runs what the build
// emits; and every operation of Lanes<T, W>, for 8-, 16-, 32- and 64-bit integers and both
// floating-point types, at 16, 32 and 64 bytes, agrees lane by lane with scalar code on
// random values and on the limits of each type (wrapping arithmetic, saturation, shifts
// across lanes, reductions, masks). Built twice by CMake: with the default flags, and with
// -march=native, so that both the generic formulas and the instructions are checked.

#include <calx/cpu/lanes.hpp>
#include <calx/cpu/target.hpp>
#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include <type_traits>

using namespace calx::cpu;

namespace {

template <class T>
const char* type_name() {
  if constexpr (std::is_same_v<T, std::uint8_t>) return "u8";
  else if constexpr (std::is_same_v<T, std::int8_t>) return "i8";
  else if constexpr (std::is_same_v<T, std::uint16_t>) return "u16";
  else if constexpr (std::is_same_v<T, std::int16_t>) return "i16";
  else if constexpr (std::is_same_v<T, std::int32_t>) return "i32";
  else if constexpr (std::is_same_v<T, std::uint32_t>) return "u32";
  else if constexpr (std::is_same_v<T, std::int64_t>) return "i64";
  else if constexpr (std::is_same_v<T, float>) return "f32";
  else return "f64";
}

// A random value of T, often at or near its limits, where saturation and wrapping are tested.
template <class T>
T draw(calx::random::Stream& r) {
  if constexpr (std::is_integral_v<T>) {
    const auto lo = std::numeric_limits<T>::min(), hi = std::numeric_limits<T>::max();
    switch (r.below(6)) {
      case 0: return lo;
      case 1: return hi;
      case 2: return static_cast<T>(lo + static_cast<T>(r.below(3)));
      case 3: return static_cast<T>(hi - static_cast<T>(r.below(3)));
      case 4: return static_cast<T>(static_cast<std::int64_t>(r.below(9)) - 4);
      default: return static_cast<T>(r.next64());
    }
  } else {
    return static_cast<T>(r.normal(0, 100));
  }
}

template <class T>
T wrap_add(T a, T b) { if constexpr (std::is_integral_v<T>) return static_cast<T>(static_cast<std::make_unsigned_t<T>>(a) + static_cast<std::make_unsigned_t<T>>(b)); else return a + b; }
template <class T>
T wrap_sub(T a, T b) { if constexpr (std::is_integral_v<T>) return static_cast<T>(static_cast<std::make_unsigned_t<T>>(a) - static_cast<std::make_unsigned_t<T>>(b)); else return a - b; }
template <class T>
T wrap_mul(T a, T b) {
  if constexpr (std::is_integral_v<T>) {
    using U = std::conditional_t<(sizeof(T) < sizeof(unsigned)), unsigned, std::make_unsigned_t<T>>;
    return static_cast<T>(static_cast<U>(static_cast<std::make_unsigned_t<T>>(a)) * static_cast<U>(static_cast<std::make_unsigned_t<T>>(b)));
  } else {
    return a * b;
  }
}
template <class T>
T saturate(long double x) {
  const long double lo = std::numeric_limits<T>::min(), hi = std::numeric_limits<T>::max();
  return static_cast<T>(std::clamp(x, lo, hi));
}

// Every operation of Lanes<T, W> against scalar code, over many random vectors.
template <class T, std::size_t W>
void lanes(calx::verification::Suite& s, calx::random::Stream& r) {
  using V = Lanes<T, W>;
  const std::string name = std::string(type_name<T>()) + " x " + std::to_string(W);
  bool arith = true, cmp = true, sat = true, shift = true, reduce = true, mask = true, memory = true;
  for (int trial = 0; trial < 400; trial++) {
    T a[W], b[W], out[W];
    for (std::size_t i = 0; i < W; i++) { a[i] = draw<T>(r); b[i] = draw<T>(r); }
    const V x = V::load(a), y = V::load(b);
    x.store(out);
    for (std::size_t i = 0; i < W; i++) memory = memory && out[i] == a[i];
    const V sum = x + y, dif = x - y, pro = x * y, hi = max(x, y), lo = min(x, y);
    const auto gt = x > y, eq = x == y;
    const V sel = select(gt, x, y);
    const T fill = draw<T>(r);
    const V up = x.shift_up(fill), down = x.shift_down(fill);
    bool any_gt = false;
    T m = a[0], n = a[0];
    for (std::size_t i = 0; i < W; i++) {
      arith = arith && sum[i] == wrap_add(a[i], b[i]) && dif[i] == wrap_sub(a[i], b[i]) && pro[i] == wrap_mul(a[i], b[i]);
      arith = arith && hi[i] == std::max(a[i], b[i]) && lo[i] == std::min(a[i], b[i]);
      cmp = cmp && gt[i] == (a[i] > b[i]) && eq[i] == (a[i] == b[i]) && sel[i] == (a[i] > b[i] ? a[i] : b[i]);
      shift = shift && up[i] == (i == 0 ? fill : a[i - 1]) && down[i] == (i + 1 == W ? fill : a[i + 1]);
      any_gt = any_gt || a[i] > b[i];
      m = std::max(m, a[i]);
      n = std::min(n, a[i]);
    }
    reduce = reduce && x.hmax() == m && x.hmin() == n;
    mask = mask && gt.any() == any_gt && (~gt).all() == !any_gt && (gt.bits() != 0) == any_gt;
    if constexpr (std::is_integral_v<T>) {
      const V ps = adds(x, y), ms = subs(x, y);
      for (std::size_t i = 0; i < W; i++)
        sat = sat && ps[i] == saturate<T>(static_cast<long double>(a[i]) + b[i]) && ms[i] == saturate<T>(static_cast<long double>(a[i]) - b[i]);
    }
  }
  s.require(name + ": load and store return the values", memory);
  s.require(name + ": +, -, *, max and min agree with scalar code", arith);
  s.require(name + ": comparisons and select agree with scalar code", cmp);
  s.require(name + ": shift_up and shift_down move every lane by one and fill the end", shift);
  s.require(name + ": hmax and hmin reduce across the lanes", reduce);
  s.require(name + ": masks: any, all and bits", mask);
  if constexpr (std::is_integral_v<T>) s.require(name + ": adds and subs saturate at the limits", sat);
}

template <class T>
void all_widths(calx::verification::Suite& s, calx::random::Stream& r) {
  lanes<T, 16 / sizeof(T)>(s, r);
  lanes<T, 32 / sizeof(T)>(s, r);
  lanes<T, 64 / sizeof(T)>(s, r);
}

}  // namespace

int main() {
  calx::verification::Suite s;

  const Target& t = host();
  s.require("the host is described", !t.architecture.empty() && t.threads >= 1, describe(t));
  s.require("the host runs what the build emits", LANE_BYTES <= std::max<std::size_t>(16, t.widest()), std::string("built for ") + ISA + ", host's best " + t.best());
  s.require("the features are consistent (AVX-512BW implies AVX-512F, AVX2 implies AVX)", (!t.avx512bw || t.avx512f) && (!t.avx2 || t.avx));
  s.require("Lanes<T> defaults to the build's widest vector", Lanes<std::int16_t>::bytes == LANE_BYTES && Lanes<std::uint8_t>::width == LANE_BYTES);

  calx::random::Stream r(20261024, 0x63707521, 0, 0);
  all_widths<std::uint8_t>(s, r);
  all_widths<std::int8_t>(s, r);
  all_widths<std::uint16_t>(s, r);
  all_widths<std::int16_t>(s, r);
  all_widths<std::int32_t>(s, r);
  all_widths<std::uint32_t>(s, r);
  all_widths<std::int64_t>(s, r);
  all_widths<float>(s, r);
  all_widths<double>(s, r);

  // A small kernel: the running maximum of a striped column, as Farrar's lazy-F loop takes it.
  {
    using V = Lanes<std::int16_t, 8>;
    V f = V::iota(-4);  // -4, -3, ..., 3
    V g = max(f, f.shift_up(-100) - V::splat(1));
    s.require("a striped step: max(f, shift_up(f) - 1) keeps the larger of each lane and its neighbour", g[0] == -4 && g[1] == -3 && g[7] == 3);
  }

  return s.done();
}
