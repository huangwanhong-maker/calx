// Philox4x32-10, a counter-based pseudorandom function (Salmon, Moraes, Dror and Shaw 2011,
// "Parallel random numbers: as easy as 1, 2, 3", Proc. SC11, §4).
//
// It is a keyed bijection of 128-bit counters: ten rounds, each two 32-bit
// multiply-high-low products with the constants M0 and M1 and a xor with the key, the key
// advanced between rounds by the Weyl constants W0 and W1 (the golden ratio and sqrt(3) - 1
// in 32-bit fixed point). Distinct counters under one key give independent-looking
// outputs; it passes TestU01's BigCrush. Because the output is a pure function of
// (counter, key), any draw can be made by any thread in any order.
//
// Checked against the known-answer vectors of Random123 (tests/kat_vectors).

#pragma once

#include <array>
#include <cstdint>

namespace calx::random {

struct Philox4x32 {
  using Counter = std::array<std::uint32_t, 4>;
  using Key = std::array<std::uint32_t, 2>;
  static constexpr std::uint32_t M0 = 0xD2511F53u, M1 = 0xCD9E8D57u;
  static constexpr std::uint32_t W0 = 0x9E3779B9u, W1 = 0xBB67AE85u;
  static constexpr int rounds = 10;

  static constexpr Counter round(const Counter& c, const Key& k) {
    const std::uint64_t p0 = static_cast<std::uint64_t>(M0) * c[0];
    const std::uint64_t p1 = static_cast<std::uint64_t>(M1) * c[2];
    const auto hi0 = static_cast<std::uint32_t>(p0 >> 32), lo0 = static_cast<std::uint32_t>(p0);
    const auto hi1 = static_cast<std::uint32_t>(p1 >> 32), lo1 = static_cast<std::uint32_t>(p1);
    return {hi1 ^ c[1] ^ k[0], lo1, hi0 ^ c[3] ^ k[1], lo0};
  }

  static constexpr Counter apply(Counter c, Key k) {
    c = round(c, k);
    for (int r = 1; r < rounds; r++) {
      k[0] += W0;
      k[1] += W1;
      c = round(c, k);
    }
    return c;
  }
};

}  // namespace calx::random
