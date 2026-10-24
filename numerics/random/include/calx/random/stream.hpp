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
// - gamma(k): Marsaglia and Tsang's squeeze on a cubed normal for shape k ≥ 1, and Γ(k + 1)
//   times U^{1/k} below it (Marsaglia and Tsang 2000, "A simple method for generating gamma
//   variables", ACM TOMS 26, 363–372, §3 and §6); beta(a, b) as X/(X + Y) of two gammas;
// - binomial(n, p): by inversion of the distribution function when n·min(p, 1 − p) < 10
//   (Kachitvichyanukul and Schmeiser 1988, "Binomial random variate generation", CACM 31,
//   216–222, BINV), and by Hörmann's transformed rejection with squeeze above it (BTRS;
//   Hörmann 1993, "The generation of binomial random variates", J. Statist. Comput. Simul.
//   46, 101–110), with log k! by Stirling's series (thread-safe, unlike the C library's lgamma,
//   which writes the global signgam).
// Others are added when a caller needs them.

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
  // Gamma with shape k > 0 and scale 1.
  double gamma(double k) {
    if (!(k > 0)) return 0;
    if (k < 1) return gamma(k + 1) * std::pow(uniform(), 1 / k);
    const double d = k - 1.0 / 3, c = 1 / std::sqrt(9 * d);
    for (;;) {
      double x, v;
      do {
        x = normal();
        v = 1 + c * x;
      } while (v <= 0);
      v = v * v * v;
      const double u = uniform(), x2 = x * x;
      if (u < 1 - 0.0331 * x2 * x2) return d * v;
      if (std::log(u) < 0.5 * x2 + d * (1 - v + std::log(v))) return d * v;
    }
  }
  double gamma(double k, double scale) { return scale * gamma(k); }
  double beta(double a, double b) {
    const double x = gamma(a), y = gamma(b);
    return x / (x + y);
  }
  // Binomial with n trials of probability p.
  std::int64_t binomial(std::int64_t n, double p) {
    if (n <= 0 || !(p > 0)) return 0;
    if (p >= 1) return n;
    if (p > 0.5) return n - binomial(n, 1 - p);
    const double q = 1 - p;
    if (static_cast<double>(n) * p < 10) {
      // Inversion: walk the distribution function from 0.
      const double s = p / q, a = (static_cast<double>(n) + 1) * s;
      for (;;) {
        double r = std::pow(q, static_cast<double>(n)), u = uniform();
        std::int64_t x = 0;
        while (u > r) {
          u -= r;
          x++;
          if (x > n) break;
          r *= a / static_cast<double>(x) - s;
        }
        if (x <= n) return x;
      }
    }
    const double nd = static_cast<double>(n), spq = std::sqrt(nd * p * q);
    const double b = 1.15 + 2.53 * spq, a = -0.0873 + 0.0248 * b + 0.01 * p, c = nd * p + 0.5;
    const double vr = 0.92 - 4.2 / b, alpha = (2.83 + 5.1 / b) * spq, lpq = std::log(p / q);
    const auto m = static_cast<std::int64_t>(std::floor((nd + 1) * p));
    const double h = log_factorial(m) + log_factorial(n - m);
    for (;;) {
      const double U = uniform() - 0.5, V = uniform();
      const double us = 0.5 - std::fabs(U);
      const auto k = static_cast<std::int64_t>(std::floor((2 * a / us + b) * U + c));
      if (k < 0 || k > n) continue;
      if (us >= 0.07 && V <= vr) return k;
      if (std::log(V * alpha / (a / (us * us) + b)) <= h - log_factorial(k) - log_factorial(n - k) + static_cast<double>(k - m) * lpq) return k;
    }
  }

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

  // log k!, exactly summed below 32 and by Stirling's series above (error below 1e-13).
  static double log_factorial(std::int64_t k) {
    if (k < 32) {
      double s = 0;
      for (std::int64_t i = 2; i <= k; i++) s += std::log(static_cast<double>(i));
      return s;
    }
    const double x = static_cast<double>(k), x2 = x * x;
    return x * std::log(x) - x + 0.5 * std::log(2 * std::numbers::pi * x) + (1.0 / 12 - (1.0 / 360 - 1.0 / (1260 * x2)) / x2) / x;
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
