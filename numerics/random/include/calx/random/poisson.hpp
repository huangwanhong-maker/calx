// Poisson variates of one fixed mean, drawn many times: by inversion of a tabulated distribution
// function.
//
// For a mean lam the table holds F_k = sum_{j <= k} p_j, p_0 = e^{-lam}, p_j = p_{j-1} lam / j,
// from k = 0 until the terms past the mode fall below 2^-60; a variate is the least k with
// U < F_k for one uniform U, found by binary search. So each variate costs one uniform and
// O(log lam) comparisons, where Knuth's product of uniforms (Stream::poisson below a mean of
// 10) costs an exponential and lam + 1 uniforms: the difference for a network's Poisson input,
// one variate per neuron per step. The distribution is exact to the rounding of the table; a U
// beyond its last entry, which happens with a probability below 1e-17, is answered by the
// stream's own sampler, as is every mean above 500 (whose e^{-lam} would lose the table's
// precision): Hörmann's transformed rejection (PTRS).
//
// The table is built once, before the draws, and read by many threads at once.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

#include "stream.hpp"

namespace calx::random {

class PoissonTable {
 public:
  explicit PoissonTable(double lam = 0) : lam_(lam) {
    if (!(lam > 0) || lam > 500) return;
    double p = std::exp(-lam), F = p;
    cdf_.push_back(F);
    for (std::int64_t k = 1;; k++) {
      p *= lam / static_cast<double>(k);
      F += p;
      cdf_.push_back(F);
      if (static_cast<double>(k) > lam && p < 0x1.0p-60) break;
    }
  }

  double mean() const { return lam_; }
  std::size_t size() const { return cdf_.size(); }

  std::int64_t operator()(Stream& z) const {
    if (!(lam_ > 0)) return 0;
    if (cdf_.empty()) return z.poisson(lam_);
    const double u = z.uniform();
    const auto it = std::upper_bound(cdf_.begin(), cdf_.end(), u);
    if (it == cdf_.end()) return z.poisson(lam_);
    return static_cast<std::int64_t>(it - cdf_.begin());
  }

 private:
  double lam_;
  std::vector<double> cdf_;
};

}  // namespace calx::random
