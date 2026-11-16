// Checks of calx/optimisation/knapsack.hpp: the 0/1 knapsack problem by branch and bound.
//
// - Brute force: on 500 random problems (n ≤ 14; real values and weights, and small integers with many ties), the value found
//   equals the greatest over every subset that fits, and the subset returned fits and has that value.
// - Dynamic programming: on problems with whole weights (n ≤ 60, capacity ≤ 400) the value equals the table's.
// - The certificate: stopped early on subset-sum problems, where every ratio is equal and the bound prunes nothing, the answer
//   says it is not exact, and the optimum lies between its value and its bound; allowed to finish, it is exact.
// - A case by hand; items of no value, of no weight and too heavy; an empty problem; refusals.

#include <calx/optimisation/knapsack.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

using calx::optimisation::knapsack;

int main() {
  calx::verification::Suite s;

  // ---- A case by hand ------------------------------------------------------------------------------------------------------------
  {
    // Values 60, 100, 120; weights 10, 20, 30; capacity 50. Ratios 6, 5, 4: greedy takes the first two (160) and cannot fit the
    // third; the optimum is the last two, 220. The relaxation takes two and two thirds of the third: 160 + 80 = 240.
    const std::vector<double> v = {60, 100, 120}, w = {10, 20, 30};
    const auto k = knapsack(v, w, 50);
    s.require("three items by hand: the two of least ratio are the best, where greedy would take the two of greatest", k.chosen == std::vector<std::size_t>{1, 2} && k.value == 220 && k.weight == 50 && k.exact && k.bound == 220);
    const auto stopped = knapsack(v, w, 50, 1);
    s.require("stopped at once, it returns greedy's 160, says it is not exact, and bounds the optimum by the relaxation's 240", !stopped.exact && stopped.value == 160 && stopped.bound == 240 && stopped.chosen == std::vector<std::size_t>{0, 1});
  }

  // ---- Brute force ---------------------------------------------------------------------------------------------------------------
  {
    calx::random::Stream rng(20261116, 3, 0, 0);
    double worst = 0;
    bool fits = true, exact = true;
    for (int trial = 0; trial < 500; trial++) {
      const std::size_t n = 1 + static_cast<std::size_t>(rng.uniform(0, 14));
      const bool whole = trial % 2 == 1;
      std::vector<double> v(n), w(n);
      for (std::size_t i = 0; i < n; i++) {
        v[i] = whole ? std::floor(rng.uniform(-1, 6)) : rng.uniform(-0.5, 5);
        w[i] = whole ? std::floor(rng.uniform(0, 5)) : rng.uniform(0, 4);
      }
      double total = 0;
      for (double x : w) total += x;
      const double c = whole ? std::floor(rng.uniform(0, total + 1)) : rng.uniform(0, total);
      double best = 0;
      for (std::size_t mask = 0; mask < (std::size_t{1} << n); mask++) {
        double value = 0, weight = 0;
        for (std::size_t i = 0; i < n; i++)
          if (mask >> i & 1) value += v[i], weight += w[i];
        if (weight <= c) best = std::max(best, value);
      }
      const auto k = knapsack(v, w, c);
      double value = 0, weight = 0;
      for (std::size_t i : k.chosen) value += v[i], weight += w[i];
      worst = std::max(worst, std::abs(k.value - best));
      fits = fits && weight <= c && std::abs(value - k.value) < 1e-12 && std::abs(weight - k.weight) < 1e-12 && std::is_sorted(k.chosen.begin(), k.chosen.end());
      exact = exact && k.exact && k.bound == k.value;
    }
    s.check("the value found equals the greatest over every subset that fits (500 problems of up to 14 items, real and tied whole numbers)", worst, 1e-12);
    s.require("the subset returned fits, has that value and weight, and is listed in order", fits);
    s.require("each search ends, and says its answer is exact", exact);
  }

  // ---- Dynamic programming on whole weights ----------------------------------------------------------------------------------------
  {
    calx::random::Stream rng(20261116, 4, 0, 0);
    double worst = 0;
    for (int trial = 0; trial < 60; trial++) {
      const std::size_t n = 20 + static_cast<std::size_t>(rng.uniform(0, 41));
      const std::size_t c = 50 + static_cast<std::size_t>(rng.uniform(0, 351));
      std::vector<double> v(n), w(n);
      for (std::size_t i = 0; i < n; i++) v[i] = rng.uniform(1, 100), w[i] = 1 + std::floor(rng.uniform(0, 60));
      std::vector<double> table(c + 1, 0);
      for (std::size_t i = 0; i < n; i++)
        for (std::size_t x = c; x >= static_cast<std::size_t>(w[i]); x--) {
          table[x] = std::max(table[x], table[x - static_cast<std::size_t>(w[i])] + v[i]);
          if (x == 0) break;
        }
      const auto k = knapsack(v, w, static_cast<double>(c));
      worst = std::max(worst, std::abs(k.value - table[c]) / table[c]);
    }
    s.check("with whole weights the value equals dynamic programming's (60 problems of 20 to 60 items, relative)", worst, 1e-12);
  }

  // ---- The certificate ----------------------------------------------------------------------------------------------------------------
  {
    // Subset sum: every value equals its weight, so every ratio is 1 and the relaxation always says "the capacity": nothing is
    // pruned until a subset fills it exactly.
    calx::random::Stream rng(20261116, 5, 0, 0);
    bool between = true, said = true, finished = true;
    for (int trial = 0; trial < 20; trial++) {
      const std::size_t n = 20;
      std::vector<double> w(n);
      for (auto& x : w) x = 1000 + std::floor(rng.uniform(0, 9000));
      const double c = 20011;  // a prime no small sum is likely to hit
      std::vector<char> reach(static_cast<std::size_t>(c) + 1, 0);
      reach[0] = 1;
      for (double x : w)
        for (std::size_t y = static_cast<std::size_t>(c); y >= static_cast<std::size_t>(x); y--)
          if (reach[y - static_cast<std::size_t>(x)]) reach[y] = 1;
      double optimum = 0;
      for (std::size_t y = 0; y <= static_cast<std::size_t>(c); y++)
        if (reach[y]) optimum = static_cast<double>(y);
      const auto stopped = knapsack(w, w, c, 200);
      between = between && stopped.value <= optimum && optimum <= stopped.bound && stopped.bound == c;
      said = said && (!stopped.exact || stopped.value == stopped.bound) && stopped.nodes <= 200;
      const auto whole = knapsack(w, w, c, 50000000);
      finished = finished && whole.exact && whole.value == optimum;
    }
    s.require("stopped early on subset-sum problems, the optimum lies between the value returned and its bound, the capacity", between);
    s.require("what it calls exact has its bound at its value, and it visits no more nodes than allowed", said);
    s.require("allowed to finish, it finds the subset sum closest beneath the capacity, as the table of reachable sums has it", finished);
  }

  // ---- Edges and refusals ---------------------------------------------------------------------------------------------------------------
  {
    const std::vector<double> v = {5, -2, 0, 7, 3}, w = {0, 0, 1, 99, 2};
    const auto k = knapsack(v, w, 2);
    s.require("an item of no weight and some value is always taken; one of no value never; one heavier than the capacity never", k.chosen == std::vector<std::size_t>{0, 4} && k.value == 8 && k.weight == 2 && k.exact);
    const auto none = knapsack(std::vector<double>{}, std::vector<double>{}, 10);
    s.require("an empty problem has an empty answer", none.chosen.empty() && none.value == 0 && none.exact && none.bound == 0);
    const auto zero = knapsack(std::vector<double>{4, 4}, std::vector<double>{1, 1}, 0);
    s.require("with no capacity nothing that weighs is taken", zero.chosen.empty() && zero.value == 0);
    const auto all = knapsack(std::vector<double>{1, 2, 3}, std::vector<double>{1, 1, 1}, 10);
    s.require("with room for everything, everything of value is taken", all.chosen.size() == 3 && all.value == 6);
    bool r1 = false, r2 = false, r3 = false, r4 = false;
    try { (void)knapsack(std::vector<double>{1}, std::vector<double>{1, 2}, 1); } catch (const std::invalid_argument&) { r1 = true; }
    try { (void)knapsack(std::vector<double>{1}, std::vector<double>{-1}, 1); } catch (const std::invalid_argument&) { r2 = true; }
    try { (void)knapsack(std::vector<double>{std::numeric_limits<double>::quiet_NaN()}, std::vector<double>{1}, 1); } catch (const std::invalid_argument&) { r3 = true; }
    try { (void)knapsack(std::vector<double>{1}, std::vector<double>{1}, -1); } catch (const std::invalid_argument&) { r4 = true; }
    s.require("values without weights, a negative weight, a value that is no number and a negative capacity are refused", r1 && r2 && r3 && r4);
  }
  return s.done();
}
