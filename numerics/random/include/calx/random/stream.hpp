// A stream of random draws, and the sampling of the standard distributions from it.
//
// A Stream is Philox4x32-10 under a 64-bit key (the seed), with three of the four counter
// words fixed by the caller to say what the stream is for (its identity: a purpose, an
// object, an entity, a step; how they are packed is the caller's) and the fourth walking.
// Two streams with different identities never share a counter, so their draws are
// independent; and a draw depends only on (seed, identity, position), never on which
// thread makes it.
//
// Distributions (exact, not approximated, unless stated):
// - uniform on (0, 1): 53 random bits, shifted by half a step so that 0 and 1 never occur;
// - below(n): uniform on {0, ..., n-1}, unbiased, by Lemire's nearly divisionless method
//   (Lemire 2019, "Fast random integer generation in an interval", ACM TOMACS 29, 3);
// - normal: Box and Muller (1958), both values of a pair used;
// - exponential: by inversion;
// - geometric: failures before the first success, by inversion;
// - Poisson: by multiplying uniforms below a mean of 10 (Knuth, TAOCP 2, §3.4.1), and by
//   Hörmann's transformed rejection with squeeze (PTRS; Hörmann 1993, "The transformed
//   rejection method for generating Poisson random variables", Insurance: Mathematics and
//   Economics 12, 39-45) above it, as NumPy's random_poisson_ptrs does.
// Others (binomial, gamma, ...) are added when a caller needs them.

#pragma once

#include <array>
#include <cmath>
#include <cstdint>
#include <numbers>

#include "philox.hpp"

namespace calx::random {

class Stream {
 public:
  Stream(std::uint64_t seed, std::uint32_t a, std::uint32_t b, std::uint32_t c)
      : key_{static_cast<std::uint32_t>(seed), static_cast<std::uint32_t>(seed >> 32)}, id_{a, b, c} {}

  // Moves to position `block` of the stream (each block is four 32-bit draws).
  void seek(std::uint32_t block) { block_ = block; pos_ = 4; has_spare_ = false; }

  std::uint32_t next32() {
    if (pos_ == 4) refill();
    return buf_[pos_++];
  }
  std::uint64_t next64() {
    const std::uint64_t hi = next32();
    return (hi << 32) | next32();
  }

  double uniform() { return (static_cast<double>(next64() >> 11) + 0.5) * 0x1.0p-53; }
  double uniform(double a, double b) { return a + (b - a) * uniform(); }

  std::uint64_t below(std::uint64_t n) {
    std::uint64_t x = next64();
    __uint128_t m = static_cast<__uint128_t>(x) * n;
    auto l = static_cast<std::uint64_t>(m);
    if (l < n) {
      const std::uint64_t t = (0 - n) % n;
      while (l < t) {
        x = next64();
        m = static_cast<__uint128_t>(x) * n;
        l = static_cast<std::uint64_t>(m);
      }
    }
    return static_cast<std::uint64_t>(m >> 64);
  }

  double normal() {
    if (has_spare_) { has_spare_ = false; return spare_; }
    const double u = uniform(), v = uniform();
    const double r = std::sqrt(-2 * std::log(u)), a = 2 * std::numbers::pi * v;
    spare_ = r * std::sin(a);
    has_spare_ = true;
    return r * std::cos(a);
  }
  double normal(double mean, double sd) { return mean + sd * normal(); }
  double lognormal(double mu, double sigma) { return std::exp(normal(mu, sigma)); }
  double exponential() { return -std::log(uniform()); }
  bool bernoulli(double p) { return uniform() < p; }

  std::uint64_t geometric(double p) {
    if (p >= 1) return 0;
    return static_cast<std::uint64_t>(std::floor(std::log(uniform()) / std::log1p(-p)));
  }

  std::int64_t poisson(double lam) {
    if (!(lam > 0)) return 0;
    if (lam < 10) {
      const double enlam = std::exp(-lam);
      std::int64_t k = 0;
      double prod = 1;
      for (;;) {
        prod *= uniform();
        if (prod > enlam) k++;
        else return k;
      }
    }
    const double slam = std::sqrt(lam), loglam = std::log(lam);
    const double b = 0.931 + 2.53 * slam, a = -0.059 + 0.02483 * b;
    const double invalpha = 1.1239 + 1.1328 / (b - 3.4), vr = 0.9277 - 3.6224 / (b - 2);
    for (;;) {
      const double U = uniform() - 0.5, V = uniform();
      const double us = 0.5 - std::fabs(U);
      const auto k = static_cast<std::int64_t>(std::floor((2 * a / us + b) * U + lam + 0.43));
      if (us >= 0.07 && V <= vr) return k;
      if (k < 0 || (us < 0.013 && V > us)) continue;
      if (std::log(V) + std::log(invalpha) - std::log(a / (us * us) + b) <= -lam + k * loglam - std::lgamma(k + 1.0)) return k;
    }
  }

 private:
  void refill() {
    buf_ = Philox4x32::apply({block_, id_[0], id_[1], id_[2]}, key_);
    block_++;
    pos_ = 0;
  }
  Philox4x32::Key key_;
  std::array<std::uint32_t, 3> id_;
  std::uint32_t block_ = 0;
  Philox4x32::Counter buf_{};
  int pos_ = 4;
  double spare_ = 0;
  bool has_spare_ = false;
};

}  // namespace calx::random
