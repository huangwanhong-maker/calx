// Greedy maximisation of a monotone submodular set function, with what bounds the best (decision 0019).
//
// A set function f on the subsets of a ground set {0, …, n − 1} is monotone if f(A) ≤ f(B) for A ⊆ B, and submodular if the gain of an
// element never grows as the set grows: f(A ∪ {e}) − f(A) ≥ f(B ∪ {e}) − f(B) for A ⊆ B, e ∉ B. Coverage, the people within reach of
// the sites opened, a facility location's best similarities, and the greatest flow into the sinks chosen (a polymatroid's rank) are
// such functions.
//
// UNDER A NUMBER k. The greedy rule adds each time the element of greatest gain. Its value is at least (1 − 1/e) of the best set of k
// elements (Nemhauser, Wolsey and Fisher 1978, Math. Programming 14, 265, Theorem 4.2), and no rule that asks f polynomially often can
// promise more in general (Nemhauser and Wolsey 1978). Lazy evaluation (Minoux 1978, Accelerated greedy algorithms, LNCS 7, 234) keeps
// each element's last gain in a heap: by submodularity a gain found earlier bounds the gain now, so an element whose fresh gain still
// tops the heap is the greedy choice, and the others need not be asked.
//
// WHAT BOUNDS THE BEST. For any set S, the best set of k elements is worth at most f(S) plus the k greatest gains with respect to S
// (Nemhauser, Wolsey and Fisher 1978, the inequality behind Theorem 4.2; Leskovec et al. 2007, KDD, Theorem 3). The bound reported is
// the least of these over the greedy sets S_0 ⊂ S_1 ⊂ …, using gains that are fresh or, being stale, still upper bounds: a certificate
// computed from the instance, usually far tighter than 1 − 1/e.
//
// UNDER A BUDGET. With costs c_e > 0 and a budget B, the better of two greedy runs (by gain, and by gain per cost, each taking what
// still fits) is worth at least (1 − 1/e)/2 of the best within the budget (Leskovec et al. 2007, Theorem 2; Khuller, Moss and Naor
// 1999, Information Processing Letters 70, 39). Its bound is f(S) plus the greatest fractional knapsack of gains within B, which is
// at least the integral one.
//
// f is asked for whole sets, its domain the caller's; f(∅) is asked once. Termination: at most k (or n) rounds of at most n questions.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <queue>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace calx::optimisation {

struct GreedySet {
  std::vector<std::size_t> chosen;  // in the order chosen
  std::vector<double> value;        // f after each choice
  double empty = 0;                 // f(∅)
  double bound = std::numeric_limits<double>::infinity();  // an upper bound on the best set of the same size (or within the budget)
  double guarantee = 0;             // the ratio promised a priori: 1 − 1/e under a number, (1 − 1/e)/2 under a budget
  std::size_t evaluations = 0;      // how often f was asked
  double total() const { return value.empty() ? empty : value.back(); }
};

namespace submodular_detail {

// Lazy greedy; `ratio` divides gains by costs; `fits` says whether an element may still be added. Gains and bounds per element.
inline GreedySet run(std::size_t n, const std::function<double(const std::vector<std::size_t>&)>& f, std::size_t k, const std::vector<double>* cost,
                     double budget, bool ratio) {
  GreedySet out;
  std::vector<std::size_t> S;
  out.empty = f(S);
  out.evaluations = 1;
  double current = out.empty;
  double spent = 0;
  // Upper bounds on gains (initially +inf: not yet asked), with the round at which each was last computed.
  std::vector<double> gain(n, std::numeric_limits<double>::infinity());
  std::vector<std::size_t> fresh(n, static_cast<std::size_t>(-1));
  std::vector<char> taken(n, 0);
  const auto key = [&](std::size_t e) { return ratio ? gain[e] / (*cost)[e] : gain[e]; };
  const auto fits = [&](std::size_t e) { return !cost || spent + (*cost)[e] <= budget * (1 + 1e-12); };
  // The bound with respect to the present set from the gains held (fresh or stale upper bounds), each finite once asked.
  const auto bound_now = [&]() {
    std::vector<std::pair<double, double>> g;  // (gain, cost)
    for (std::size_t e = 0; e < n; e++)
      if (!taken[e]) g.push_back({std::max(0.0, gain[e]), cost ? (*cost)[e] : 1.0});
    if (!cost) {
      std::sort(g.begin(), g.end(), [](const auto& a, const auto& b) { return a.first > b.first; });
      double sum = 0;
      for (std::size_t i = 0; i < std::min(k, g.size()); i++) sum += g[i].first;
      return current + sum;
    }
    std::sort(g.begin(), g.end(), [](const auto& a, const auto& b) { return a.first * b.second > b.first * a.second; });
    double room = budget, sum = 0;
    for (const auto& [v, c] : g) {
      if (room <= 0) break;
      const double take = std::min(1.0, room / c);
      sum += take * v;
      room -= take * c;
    }
    return current + sum;
  };
  // Every gain is asked once with respect to ∅ so that the first bound is finite.
  for (std::size_t e = 0; e < n; e++) {
    S.push_back(e);
    gain[e] = f(S) - current;
    out.evaluations++;
    S.pop_back();
    fresh[e] = 0;
  }
  out.bound = bound_now();
  using Item = std::pair<double, std::size_t>;
  std::priority_queue<Item> heap;
  for (std::size_t e = 0; e < n; e++) heap.push({key(e), e});
  std::size_t round = 0;
  while (S.size() < k && !heap.empty()) {
    auto [kv, e] = heap.top();
    heap.pop();
    if (taken[e] || !fits(e)) continue;
    if (fresh[e] != round) {
      S.push_back(e);
      gain[e] = f(S) - current;
      out.evaluations++;
      S.pop_back();
      fresh[e] = round;
      heap.push({key(e), e});
      continue;
    }
    if (gain[e] <= 0) break;  // nothing gains any more: the set stays as it is
    S.push_back(e);
    taken[e] = 1;
    current += gain[e];
    if (cost) spent += (*cost)[e];
    out.chosen.push_back(e);
    out.value.push_back(current);
    round++;
    out.bound = std::min(out.bound, bound_now());
  }
  return out;
}

inline void check(std::size_t n, const std::function<double(const std::vector<std::size_t>&)>& f, const char* who) {
  if (!f) throw std::invalid_argument(std::string(who) + ": a set function");
  if (n == 0) throw std::invalid_argument(std::string(who) + ": a ground set with elements");
}

}  // namespace submodular_detail

// The greedy set of at most k elements, lazily, with the bound on the best set of k.
inline GreedySet greedy_submodular(std::size_t n, const std::function<double(const std::vector<std::size_t>&)>& f, std::size_t k) {
  const char* who = "calx::optimisation::greedy_submodular";
  submodular_detail::check(n, f, who);
  if (k == 0 || k > n) throw std::invalid_argument(std::string(who) + ": 1 ≤ k ≤ the ground set's size");
  auto out = submodular_detail::run(n, f, k, nullptr, 0, false);
  out.guarantee = 1 - std::exp(-1.0);
  return out;
}

// The better of the greedy sets by gain and by gain per cost within a budget, with the bound on the best within it.
inline GreedySet greedy_submodular_budget(std::size_t n, const std::function<double(const std::vector<std::size_t>&)>& f, const std::vector<double>& cost,
                                          double budget) {
  const char* who = "calx::optimisation::greedy_submodular_budget";
  submodular_detail::check(n, f, who);
  if (cost.size() != n) throw std::invalid_argument(std::string(who) + ": a cost for each element");
  for (double c : cost)
    if (!(c > 0 && std::isfinite(c))) throw std::invalid_argument(std::string(who) + ": costs are finite and > 0");
  if (!(budget > 0 && std::isfinite(budget))) throw std::invalid_argument(std::string(who) + ": a budget > 0");
  auto by_gain = submodular_detail::run(n, f, n, &cost, budget, false);
  auto by_ratio = submodular_detail::run(n, f, n, &cost, budget, true);
  GreedySet out = by_gain.total() >= by_ratio.total() ? by_gain : by_ratio;
  out.bound = std::min(by_gain.bound, by_ratio.bound);
  out.evaluations = by_gain.evaluations + by_ratio.evaluations;
  out.guarantee = (1 - std::exp(-1.0)) / 2;
  return out;
}

}  // namespace calx::optimisation
