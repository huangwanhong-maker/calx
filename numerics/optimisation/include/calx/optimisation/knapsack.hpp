// The 0/1 knapsack problem: the subset of greatest total value within a capacity (decision 0018).
//
// THE PROBLEM. Given items with values v_i and weights w_i ≥ 0 and a capacity C ≥ 0, choose x_i ∈ {0, 1} maximising Σ v_i x_i
// subject to Σ w_i x_i ≤ C. It is NP-hard (Karp 1972), and ordinary instances are solved fast.
//
// THE BOUND. Relaxing x_i to [0, 1] gives a linear programme whose optimum is found by sorting: take the items in the order of
// v_i / w_i, greatest first, until one does not fit, and take the fraction of that one that does (Dantzig 1957, Discrete-variable
// extremum problems, Operations Research 5, 266–277). Its value bounds every 0/1 solution from above.
//
// THE ALGORITHM. Branch and bound, depth first, in the order of the ratios (Horowitz and Sahni 1974, Computing partitions with
// applications to the knapsack problem, J. ACM 21, 277–292; Martello and Toth 1990, Knapsack Problems, Wiley, ch. 2): at each node
// the Dantzig bound of what is still undecided is added to what is taken; a node that cannot beat the best found is left. The
// first leaf is the greedy solution, which is already within the most valuable single item of the optimum. An item of no positive
// value is never taken; one of no weight and positive value always is; one heavier than the capacity never.
//
// THE CERTIFICATE. When the search ends, the solution is optimal (`exact`), and `bound` is its value. When it is stopped at
// `node_limit` nodes, the best solution found is returned with `exact` false and `bound` the Dantzig bound at the root: the
// optimum lies between `value` and `bound`. Values and weights must be finite. Equal values to the last bit may be reached by
// another subset; the value is the optimum either way.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numeric>
#include <span>
#include <stdexcept>
#include <vector>

namespace calx::optimisation {

struct Knapsack {
  std::vector<std::size_t> chosen;  // the items taken, in increasing order
  double value = 0, weight = 0;     // their totals
  double bound = 0;                 // an upper bound on every solution's value: `value` when `exact`
  bool exact = true;                // whether the search ended before the node limit
  std::size_t nodes = 0;            // nodes of the search visited
};

namespace detail {

struct KnapsackSearch {
  std::vector<double> v, w;         // the items worth deciding, in the order of their ratios
  std::size_t limit = 0, nodes = 0;
  bool complete = true;
  double best = 0;
  std::vector<char> taken, best_taken;

  // The Dantzig bound of the items from k on, in a capacity c.
  double relaxed(std::size_t k, double c) const {
    double total = 0;
    for (; k < v.size(); k++) {
      if (w[k] <= c) { total += v[k]; c -= w[k]; }
      else return total + v[k] * (c / w[k]);
    }
    return total;
  }

  void go(std::size_t k, double value, double c) {
    if (nodes >= limit) { complete = false; return; }
    nodes++;
    if (value > best) { best = value; best_taken = taken; }
    if (k == v.size() || value + relaxed(k, c) <= best) return;
    if (w[k] <= c) {
      taken[k] = 1;
      go(k + 1, value + v[k], c - w[k]);
      taken[k] = 0;
    }
    go(k + 1, value, c);
  }
};

}  // namespace detail

// The subset of the items of greatest total value whose total weight is at most `capacity`.
inline Knapsack knapsack(std::span<const double> value, std::span<const double> weight, double capacity, std::size_t node_limit = 2000000) {
  const std::size_t n = value.size();
  if (weight.size() != n) throw std::invalid_argument("calx::optimisation::knapsack: a weight for each value");
  if (!(capacity >= 0) || !std::isfinite(capacity)) throw std::invalid_argument("calx::optimisation::knapsack: the capacity is finite and not negative");
  for (std::size_t i = 0; i < n; i++)
    if (!std::isfinite(value[i]) || !std::isfinite(weight[i]) || weight[i] < 0) throw std::invalid_argument("calx::optimisation::knapsack: values are finite, weights finite and not negative");

  Knapsack out;
  std::vector<std::size_t> open;  // the items worth a decision
  for (std::size_t i = 0; i < n; i++) {
    if (value[i] <= 0 || weight[i] > capacity) continue;
    if (weight[i] == 0) { out.chosen.push_back(i); out.value += value[i]; continue; }
    open.push_back(i);
  }
  std::stable_sort(open.begin(), open.end(), [&](std::size_t a, std::size_t b) { return value[a] * weight[b] > value[b] * weight[a]; });

  detail::KnapsackSearch s;
  for (std::size_t i : open) s.v.push_back(value[i]), s.w.push_back(weight[i]);
  s.limit = std::max<std::size_t>(node_limit, open.size() + 1);  // room at least for the greedy leaf
  s.taken.assign(open.size(), 0);
  s.best_taken = s.taken;
  const double root = s.relaxed(0, capacity);
  s.go(0, 0, capacity);

  for (std::size_t k = 0; k < open.size(); k++)
    if (s.best_taken[k]) out.chosen.push_back(open[k]);
  std::sort(out.chosen.begin(), out.chosen.end());
  const double free_value = out.value;
  out.value = 0;
  for (std::size_t i : out.chosen) out.value += value[i], out.weight += weight[i];
  out.exact = s.complete;
  out.nodes = s.nodes;
  out.bound = s.complete ? out.value : std::max(out.value, free_value + root);
  return out;
}

}  // namespace calx::optimisation
