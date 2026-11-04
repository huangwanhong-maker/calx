// Categorical variates from a Stream: one draw from a law given by weights, and many draws from one
// law by Walker's alias method.
//
// categorical(rng, w): an index j with probability w_j / Σ w (weights non-negative and finite, not all
// zero), by inversion: one uniform u, and the least j with w_j > 0 and Σ_{i≤j} w_i > u Σ w (Devroye
// 1986, Non-Uniform Random Variate Generation, ch. III). O(K) a draw; an index of weight zero is never
// drawn.
//
// AliasTable(w): Walker's alias method (Walker 1977, "An efficient method for generating discrete
// random variables with general distributions", ACM TOMS 3, 253–256), built by Vose's stable
// construction (Vose 1991, "A linear algorithm for generating random numbers with a given
// distribution", IEEE TSE 17, 972–975) in O(K). Column i keeps itself with probability prob_i and
// otherwise gives alias_i; a draw takes one bounded integer and one uniform, O(1). The law it draws
// is exactly (prob_j + Σ_{i: alias_i = j} (1 − prob_i)) / K, which equals w_j / Σ w up to the rounding
// of the construction; an index of weight zero is never drawn.

#pragma once

#include <calx/random/stream.hpp>

#include <cmath>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace calx::random {

namespace detail {
inline double categorical_total(std::span<const double> w, const char* who) {
  double total = 0;
  for (const double x : w) {
    if (!(x >= 0) || !std::isfinite(x)) throw std::invalid_argument(std::string(who) + ": weights non-negative and finite");
    total += x;
  }
  if (!(total > 0)) throw std::invalid_argument(std::string(who) + ": at least one positive weight");
  return total;
}
}  // namespace detail

inline std::size_t categorical(Stream& rng, std::span<const double> w) {
  const double total = detail::categorical_total(w, "calx::random::categorical");
  const double u = rng.uniform() * total;
  double sum = 0;
  std::size_t last = 0;
  for (std::size_t j = 0; j < w.size(); j++) {
    if (w[j] == 0) continue;
    sum += w[j], last = j;
    if (sum > u) return j;
  }
  return last;  // u within rounding of the total: the last index of positive weight
}

class AliasTable {
 public:
  AliasTable() = default;
  explicit AliasTable(std::span<const double> w) : prob_(w.size(), 0.0), alias_(w.size(), 0) {
    const double total = detail::categorical_total(w, "calx::random::AliasTable");
    const std::size_t n = w.size();
    std::vector<double> p(n);
    std::vector<std::size_t> small, large;
    std::size_t positive = 0;
    for (std::size_t i = 0; i < n; i++) {
      p[i] = w[i] * static_cast<double>(n) / total;
      (p[i] < 1 ? small : large).push_back(i);
      if (w[i] > 0) positive = i;
    }
    while (!small.empty() && !large.empty()) {
      const std::size_t l = small.back(), g = large.back();
      small.pop_back(), large.pop_back();
      prob_[l] = p[l], alias_[l] = g;
      p[g] = (p[g] + p[l]) - 1;
      (p[g] < 1 ? small : large).push_back(g);
    }
    for (const std::size_t g : large) prob_[g] = 1, alias_[g] = g;
    // Left only by rounding: a column whose mass is within rounding of 1 keeps itself, unless its
    // weight is zero, when it gives an index of positive weight.
    for (const std::size_t l : small) {
      if (w[l] > 0) prob_[l] = 1, alias_[l] = l;
      else prob_[l] = 0, alias_[l] = positive;
    }
  }

  std::size_t size() const { return prob_.size(); }
  std::size_t operator()(Stream& rng) const {
    const std::size_t i = static_cast<std::size_t>(rng.below(prob_.size()));
    return rng.uniform() < prob_[i] ? i : alias_[i];
  }
  // The law the table draws: (prob_j + Σ_{i: alias_i = j} (1 − prob_i)) / K.
  std::vector<double> law() const {
    std::vector<double> q(prob_.size(), 0.0);
    const double k = static_cast<double>(prob_.size());
    for (std::size_t i = 0; i < prob_.size(); i++) q[i] += prob_[i] / k, q[alias_[i]] += (1 - prob_[i]) / k;
    return q;
  }
  const std::vector<double>& probabilities() const { return prob_; }
  const std::vector<std::size_t>& aliases() const { return alias_; }

 private:
  std::vector<double> prob_;
  std::vector<std::size_t> alias_;
};

}  // namespace calx::random
