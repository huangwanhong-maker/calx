// Checks of calx/optimisation/submodular.hpp: greedy maximisation of monotone submodular functions, with a bound on the best.
//
// - Weighted coverage on 200 random instances (n = 14, k = 4): greedy reaches at least 1 − 1/e of the best of every 4-set, listed; the
//   bound it reports is at least the best; lazy evaluation picks what the plain greedy rule picks.
// - Facility location (the sum over clients of their best similarity among the chosen) on 100 instances: the same.
// - Under a budget (n = 12, costs and budget random): the better of the two greedy runs reaches (1 − 1/e)/2 of the best within the
//   budget, listed; its set fits; its bound is at least the best.
// - Lazy evaluation asks f far less than n each round on a coverage of 400 elements.
// - Refusals: k beyond the ground set, a cost that is not positive.

#include <calx/optimisation/submodular.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>
#include <vector>

namespace co = calx::optimisation;

namespace {

bool throws(auto&& f) {
  try {
    f();
  } catch (const std::invalid_argument&) {
    return true;
  }
  return false;
}

// Weighted coverage: each element covers items; f(S) = the weight of the items covered.
struct Coverage {
  std::vector<std::vector<std::size_t>> covers;
  std::vector<double> weight;
  double operator()(const std::vector<std::size_t>& S) const {
    std::vector<char> hit(weight.size(), 0);
    double sum = 0;
    for (std::size_t e : S)
      for (std::size_t i : covers[e])
        if (!hit[i]) {
          hit[i] = 1;
          sum += weight[i];
        }
    return sum;
  }
};

Coverage random_coverage(calx::random::Stream& rng, std::size_t n, std::size_t items, double p) {
  Coverage c;
  c.covers.resize(n);
  for (std::size_t i = 0; i < items; i++) c.weight.push_back(rng.uniform(0.1, 1.0));
  for (std::size_t e = 0; e < n; e++)
    for (std::size_t i = 0; i < items; i++)
      if (rng.uniform() < p) c.covers[e].push_back(i);
  return c;
}

// The plain greedy rule, asking every element each round.
std::vector<std::size_t> plain(std::size_t n, const std::function<double(const std::vector<std::size_t>&)>& f, std::size_t k) {
  std::vector<std::size_t> S;
  std::vector<char> taken(n, 0);
  for (std::size_t r = 0; r < k; r++) {
    double best = -std::numeric_limits<double>::infinity();
    std::size_t arg = 0;
    for (std::size_t e = 0; e < n; e++) {
      if (taken[e]) continue;
      S.push_back(e);
      const double v = f(S);
      S.pop_back();
      if (v > best) best = v, arg = e;
    }
    S.push_back(arg);
    taken[arg] = 1;
  }
  return S;
}

}  // namespace

int main() {
  calx::verification::Suite s;
  const double ratio = 1 - std::exp(-1.0);
  {
    calx::random::Stream rng(19, 0, 1, 2);
    bool guaranteed = true, bounded = true, same = true;
    for (int trial = 0; trial < 200; trial++) {
      const std::size_t n = 14, k = 4;
      const Coverage c = random_coverage(rng, n, 30, 0.12);
      const std::function<double(const std::vector<std::size_t>&)> f = c;
      const auto G = co::greedy_submodular(n, f, k);
      double best = 0;
      for (unsigned mask = 0; mask < (1u << n); mask++) {
        if (__builtin_popcount(mask) != static_cast<int>(k)) continue;
        std::vector<std::size_t> S;
        for (std::size_t e = 0; e < n; e++)
          if (mask >> e & 1u) S.push_back(e);
        best = std::max(best, f(S));
      }
      guaranteed = guaranteed && G.total() >= ratio * best - 1e-12;
      bounded = bounded && G.bound >= best - 1e-9 && G.bound >= G.total();
      std::vector<std::size_t> P = plain(n, f, G.chosen.size());
      same = same && std::fabs(f(P) - G.total()) < 1e-9;
    }
    s.require("weighted coverage, 200 instances: greedy reaches 1 − 1/e of the best 4-set, listed", guaranteed);
    s.require("and the bound it reports is never below the best", bounded);
    s.require("lazy evaluation reaches what the plain greedy rule reaches", same);
  }
  {
    calx::random::Stream rng(23, 3, 4, 5);
    bool guaranteed = true, bounded = true;
    for (int trial = 0; trial < 100; trial++) {
      const std::size_t n = 12, clients = 25, k = 3;
      std::vector<std::vector<double>> sim(clients, std::vector<double>(n));
      for (auto& row : sim)
        for (double& x : row) x = rng.uniform();
      const std::function<double(const std::vector<std::size_t>&)> f = [&](const std::vector<std::size_t>& S) {
        double sum = 0;
        for (const auto& row : sim) {
          double b = 0;
          for (std::size_t e : S) b = std::max(b, row[e]);
          sum += b;
        }
        return sum;
      };
      const auto G = co::greedy_submodular(n, f, k);
      double best = 0;
      for (unsigned mask = 0; mask < (1u << n); mask++) {
        if (__builtin_popcount(mask) != static_cast<int>(k)) continue;
        std::vector<std::size_t> S;
        for (std::size_t e = 0; e < n; e++)
          if (mask >> e & 1u) S.push_back(e);
        best = std::max(best, f(S));
      }
      guaranteed = guaranteed && G.total() >= ratio * best - 1e-12;
      bounded = bounded && G.bound >= best - 1e-9;
    }
    s.require("facility location, 100 instances: greedy reaches 1 − 1/e of the best 3-set, its bound above the best", guaranteed && bounded);
  }
  {
    calx::random::Stream rng(29, 6, 7, 8);
    bool guaranteed = true, bounded = true, fits = true;
    for (int trial = 0; trial < 150; trial++) {
      const std::size_t n = 12;
      const Coverage c = random_coverage(rng, n, 30, 0.15);
      const std::function<double(const std::vector<std::size_t>&)> f = c;
      std::vector<double> cost(n);
      for (double& x : cost) x = rng.uniform(1, 5);
      const double budget = rng.uniform(4, 10);
      const auto G = co::greedy_submodular_budget(n, f, cost, budget);
      double best = 0;
      for (unsigned mask = 0; mask < (1u << n); mask++) {
        double spent = 0;
        std::vector<std::size_t> S;
        for (std::size_t e = 0; e < n; e++)
          if (mask >> e & 1u) S.push_back(e), spent += cost[e];
        if (spent <= budget) best = std::max(best, f(S));
      }
      double spent = 0;
      for (std::size_t e : G.chosen) spent += cost[e];
      fits = fits && spent <= budget * (1 + 1e-12);
      guaranteed = guaranteed && G.total() >= (ratio / 2) * best - 1e-12;
      bounded = bounded && G.bound >= best - 1e-9;
    }
    s.require("under a budget, 150 instances: the set fits, reaches (1 − 1/e)/2 of the best within the budget, its bound above the best",
              fits && guaranteed && bounded);
  }
  {
    calx::random::Stream rng(31, 9, 10, 11);
    const std::size_t n = 400, k = 20;
    const Coverage c = random_coverage(rng, n, 2000, 0.01);
    const auto G = co::greedy_submodular(n, c, k);
    s.require("lazy evaluation asks f fewer than a third of the n k times the plain rule does, on 400 elements and k = 20",
              G.evaluations < n * k / 3 && G.chosen.size() == k);
  }
  {
    const std::function<double(const std::vector<std::size_t>&)> f = [](const std::vector<std::size_t>& S) { return static_cast<double>(S.size()); };
    s.require("refusals: k beyond the ground set, a cost that is not positive",
              throws([&] { co::greedy_submodular(3, f, 4); }) && throws([&] { co::greedy_submodular_budget(2, f, {1, 0}, 3); }));
  }
  return s.done();
}
