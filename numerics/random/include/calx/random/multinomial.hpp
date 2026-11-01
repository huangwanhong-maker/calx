// Multinomial variates from a Stream.
//
// multinomial(rng, n, p): counts (k_0, …, k_{K−1}) with Σ k = n of n independent draws among K
// categories with probabilities proportional to p (non-negative, not all zero; normalised here). By
// conditional binomials (Davis 1993, "The computer generation of multinomial random variates",
// Comput. Stat. Data Anal. 16, 205–217, the conditional distribution method): k_0 ~ Binomial(n, p_0),
// then k_j ~ Binomial(n − Σ_{i<j} k_i, p_j / Σ_{i≥j} p_i), each by the Stream's exact binomial; the
// last category takes what is left. Exact: the chain rule of the multinomial law. O(K) binomial
// draws, each O(1) on average.

#pragma once

#include <calx/random/stream.hpp>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <vector>

namespace calx::random {

inline std::vector<std::int64_t> multinomial(Stream& rng, std::int64_t n, std::span<const double> p) {
  if (n < 0) throw std::invalid_argument("calx::random::multinomial: n >= 0");
  std::vector<std::int64_t> k(p.size(), 0);
  // The tail sums Σ_{i≥j} p_i, accumulated from the end so that the last ratio is exactly 1.
  std::vector<double> tail(p.size() + 1, 0.0);
  for (std::size_t j = p.size(); j-- > 0;) {
    if (!(p[j] >= 0) || !std::isfinite(p[j])) throw std::invalid_argument("calx::random::multinomial: probabilities non-negative and finite");
    tail[j] = tail[j + 1] + p[j];
  }
  if (p.empty() || !(tail[0] > 0)) {
    if (n == 0) return k;
    throw std::invalid_argument("calx::random::multinomial: some probability positive");
  }
  std::int64_t left = n;
  for (std::size_t j = 0; j < p.size() && left > 0; j++) {
    if (!(tail[j] > 0)) break;
    const double q = p[j] / tail[j];
    k[j] = q >= 1 ? left : rng.binomial(left, q);
    left -= k[j];
  }
  return k;
}

}  // namespace calx::random
