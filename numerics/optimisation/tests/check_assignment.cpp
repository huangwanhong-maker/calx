// Checks of calx/optimisation/assignment.hpp: the linear assignment problem by the Hungarian method.
//
// - Brute force: on 500 random problems (n ≤ 6 rows, m ≤ 8 columns; real costs, and integer costs
//   with many ties), the cost found equals the least over every injection.
// - The dual certificate: u_i + v_j ≤ c_ij everywhere, with equality on the assignment, v_j ≤ 0 on
//   the unassigned columns, and Σu + Σv equal to the cost (strong duality), on problems up to 60 × 80.
// - The assignment is injective; a permutation matrix of costs −1 is recovered; maximisation by −s;
//   refusals (more rows than columns, a non-finite cost).

#include <calx/optimisation/assignment.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

using calx::linalg::Rectangular;

namespace {

// The least cost over all injections of rows into columns, by recursion.
double brute(const Rectangular& c, std::size_t i, std::vector<bool>& used) {
  if (i == c.rows) return 0;
  double best = std::numeric_limits<double>::infinity();
  for (std::size_t j = 0; j < c.cols; j++) {
    if (used[j]) continue;
    used[j] = true;
    best = std::min(best, c(i, j) + brute(c, i + 1, used));
    used[j] = false;
  }
  return best;
}

}  // namespace

int main() {
  calx::verification::Suite s;
  calx::random::Stream z(1, 0, 0, 0);

  double worst = 0;
  bool injective = true;
  for (int rep = 0; rep < 500; rep++) {
    const std::size_t n = 1 + z.below(6), m = n + z.below(3);
    Rectangular c(n, m);
    for (auto& v : c.a) v = rep % 2 ? static_cast<double>(z.below(4)) : z.normal(0, 3);
    const auto a = calx::optimisation::assignment(c);
    std::vector<bool> used(m, false);
    worst = std::max(worst, std::fabs(a.cost - brute(c, 0, used)));
    std::vector<std::size_t> cols = a.column;
    std::sort(cols.begin(), cols.end());
    injective = injective && std::adjacent_find(cols.begin(), cols.end()) == cols.end() && cols.back() < m;
  }
  s.check("brute force: the least cost over every injection (500 problems, n ≤ 6, m ≤ 8, real and tied integer costs)", worst, 1e-12);
  s.require("the assignment is injective", injective);

  double dual_gap = 0, feasibility = 0, slackness = 0, unassigned = 0;
  for (int rep = 0; rep < 40; rep++) {
    const std::size_t n = 5 + z.below(56), m = n + z.below(21);
    Rectangular c(n, m);
    for (auto& v : c.a) v = z.uniform(0, 10);
    const auto a = calx::optimisation::assignment(c);
    std::vector<bool> taken(m, false);
    for (std::size_t i = 0; i < n; i++) taken[a.column[i]] = true;
    double dual = 0;
    for (const double u : a.u) dual += u;
    for (const double v : a.v) dual += v;
    dual_gap = std::max(dual_gap, std::fabs(dual - a.cost) / (1 + a.cost));
    for (std::size_t i = 0; i < n; i++)
      for (std::size_t j = 0; j < m; j++) feasibility = std::max(feasibility, a.u[i] + a.v[j] - c(i, j));
    for (std::size_t i = 0; i < n; i++) slackness = std::max(slackness, std::fabs(a.u[i] + a.v[a.column[i]] - c(i, a.column[i])));
    for (std::size_t j = 0; j < m; j++)
      if (!taken[j]) unassigned = std::max(unassigned, a.v[j]);
  }
  s.check("dual certificate: u_i + v_j ≤ c_ij (largest excess; problems up to 60 × 80)", std::max(0.0, feasibility), 1e-10);
  s.check("dual certificate: equality on the assignment (complementary slackness)", slackness, 1e-10);
  s.check("dual certificate: v_j ≤ 0 on the unassigned columns", std::max(0.0, unassigned), 1e-12);
  s.check("dual certificate: Σu + Σv equals the cost (strong duality; relative)", dual_gap, 1e-12);

  // A hidden permutation of cost −1 among costs 0; maximising a score by its negative.
  std::vector<std::size_t> perm(9);
  std::iota(perm.begin(), perm.end(), 0);
  for (std::size_t i = 9; i > 1; i--) std::swap(perm[i - 1], perm[z.below(i)]);
  Rectangular c(9, 9);
  for (std::size_t i = 0; i < 9; i++) c(i, perm[i]) = -1;
  const auto a = calx::optimisation::assignment(c);
  s.require("a hidden permutation of cost −1 recovered", a.column == perm && a.cost == -9);
  bool refused = false;
  try {
    calx::optimisation::assignment(Rectangular(3, 2));
  } catch (const std::invalid_argument&) {
    refused = true;
  }
  bool refused_nan = false;
  try {
    Rectangular bad(2, 2);
    bad(0, 1) = std::numeric_limits<double>::infinity();
    calx::optimisation::assignment(bad);
  } catch (const std::invalid_argument&) {
    refused_nan = true;
  }
  s.require("refusals: more rows than columns, and an infinite cost", refused && refused_nan);
  return s.done();
}
