// Checks of calx/optimisation's simplex method on textbook programmes with known optima:
// Chvátal (1983), ch. 2, maximise 5x₁ + 4x₂ + 3x₃ (optimum 13 at (2, 0, 1), duals (1, 0, 1));
// Chvátal's degenerate example, p. 31, on which the largest-coefficient rule cycles and Bland's
// rule does not (optimum 1 at (1, 0, 1, 0)); the Wyndor Glass problem of Hillier and Lieberman
// (optimum 36 at (2, 6), shadow prices (0, 3/2, 1)); their radiation-therapy problem with an
// equality and a ≥ row (optimum 5.25 at (7.5, 4.5)); an infeasible and an unbounded programme;
// bounds, a free variable and a redundant equality; complementary slackness and strong duality
// at every optimum; and the values of zero-sum games: rock–paper–scissors (0, uniform strategies)
// and a 2 × 2 game with a mixed equilibrium (1/7, at p = (3/7, 4/7)).

#include <calx/linalg/rectangular.hpp>
#include <calx/optimisation/linear_programming.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <limits>
#include <string>
#include <vector>

using namespace calx::optimisation;
using calx::linalg::Rectangular;

namespace {

constexpr double inf = std::numeric_limits<double>::infinity();

// Complementary slackness, dual feasibility and strong duality for x ≥ 0 programmes.
void duality(calx::verification::Suite& s, const std::string& name, const LinearProgram& lp, const LpResult& r) {
  double slack = 0, sign = 0, reduced = 0, dual_value = 0;
  for (std::size_t i = 0; i < lp.b_ub.size(); i++) {
    double ax = 0;
    for (std::size_t j = 0; j < lp.c.size(); j++) ax += lp.a_ub(i, j) * r.x[j];
    slack = std::max(slack, std::fabs(r.duals_ub[i] * (lp.b_ub[i] - ax)));
    sign = std::max(sign, r.duals_ub[i]);  // ≤ 0
    dual_value += lp.b_ub[i] * r.duals_ub[i];
  }
  for (std::size_t i = 0; i < lp.b_eq.size(); i++) dual_value += lp.b_eq[i] * r.duals_eq[i];
  for (std::size_t j = 0; j < lp.c.size(); j++) {
    double d = lp.c[j];
    for (std::size_t i = 0; i < lp.b_ub.size(); i++) d -= lp.a_ub(i, j) * r.duals_ub[i];
    for (std::size_t i = 0; i < lp.b_eq.size(); i++) d -= lp.a_eq(i, j) * r.duals_eq[i];
    reduced = std::max(reduced, std::max(-d, std::fabs(d * r.x[j])));  // d ≥ 0, d x = 0
  }
  s.check(name + ": complementary slackness, max |yᵢ (bᵢ − aᵢᵀx)|", slack, 1e-10);
  s.check(name + ": dual feasibility, y_ub ≤ 0", std::max(sign, 0.0), 1e-12);
  s.check(name + ": reduced costs c − Aᵀy ≥ 0 and complementary to x", reduced, 1e-10);
  s.check(name + ": strong duality, cᵀx = bᵀy", r.value - dual_value, 1e-10);
}

}  // namespace

int main() {
  calx::verification::Suite s;

  // Chvátal (1983), ch. 2: maximise 5x₁ + 4x₂ + 3x₃.
  {
    LinearProgram lp;
    lp.c = {-5, -4, -3};
    lp.a_ub = Rectangular::from_rows({{2, 3, 1}, {4, 1, 2}, {3, 4, 2}});
    lp.b_ub = {5, 11, 8};
    const auto r = simplex(lp);
    s.require("simplex: Chvátal's first example is optimal", r.status == LpStatus::optimal);
    s.check("simplex: its optimum 13", -r.value - 13, 1e-12);
    s.check("simplex: at (2, 0, 1)", std::fabs(r.x[0] - 2) + std::fabs(r.x[1]) + std::fabs(r.x[2] - 1), 1e-12);
    s.check("simplex: its duals (1, 0, 1) for the maximum, −(1, 0, 1) as marginals of the minimum",
            std::fabs(r.duals_ub[0] + 1) + std::fabs(r.duals_ub[1]) + std::fabs(r.duals_ub[2] + 1), 1e-12);
    duality(s, "simplex (Chvátal ch. 2)", lp, r);
  }

  // Chvátal (1983), p. 31: the degenerate example that cycles under the largest-coefficient rule.
  {
    LinearProgram lp;
    lp.c = {-10, 57, 9, 24};
    lp.a_ub = Rectangular::from_rows({{0.5, -5.5, -2.5, 9}, {0.5, -1.5, -0.5, 1}, {1, 0, 0, 0}});
    lp.b_ub = {0, 0, 1};
    const auto r = simplex(lp);
    s.require("simplex: Chvátal's cycling example terminates under Bland's rule", r.status == LpStatus::optimal, std::to_string(r.iterations) + " pivots");
    s.check("simplex: its optimum 1", -r.value - 1, 1e-12);
    s.check("simplex: at (1, 0, 1, 0)", std::fabs(r.x[0] - 1) + std::fabs(r.x[1]) + std::fabs(r.x[2] - 1) + std::fabs(r.x[3]), 1e-12);
    duality(s, "simplex (Chvátal p. 31)", lp, r);
  }

  // Hillier and Lieberman, the Wyndor Glass Co.: maximise 3x₁ + 5x₂.
  {
    LinearProgram lp;
    lp.c = {-3, -5};
    lp.a_ub = Rectangular::from_rows({{1, 0}, {0, 2}, {3, 2}});
    lp.b_ub = {4, 12, 18};
    const auto r = simplex(lp);
    s.check("simplex: Wyndor Glass, optimum 36", -r.value - 36, 1e-12);
    s.check("simplex: at (2, 6)", std::fabs(r.x[0] - 2) + std::fabs(r.x[1] - 6), 1e-12);
    s.check("simplex: shadow prices (0, 3/2, 1)", std::fabs(r.duals_ub[0]) + std::fabs(r.duals_ub[1] + 1.5) + std::fabs(r.duals_ub[2] + 1), 1e-12);
    duality(s, "simplex (Wyndor)", lp, r);
  }

  // Hillier and Lieberman, radiation therapy: minimise 0.4x₁ + 0.5x₂ with ≤, = and ≥ rows.
  {
    LinearProgram lp;
    lp.c = {0.4, 0.5};
    lp.a_ub = Rectangular::from_rows({{0.3, 0.1}, {-0.6, -0.4}});
    lp.b_ub = {2.7, -6};
    lp.a_eq = Rectangular::from_rows({{0.5, 0.5}});
    lp.b_eq = {6};
    const auto r = simplex(lp);
    s.check("simplex: radiation therapy, optimum 5.25", r.value - 5.25, 1e-12);
    s.check("simplex: at (7.5, 4.5)", std::fabs(r.x[0] - 7.5) + std::fabs(r.x[1] - 4.5), 1e-11);
    duality(s, "simplex (radiation therapy)", lp, r);
  }

  {
    LinearProgram lp;
    lp.c = {1, 1};
    lp.a_ub = Rectangular::from_rows({{1, 1}, {-1, -1}});
    lp.b_ub = {1, -3};
    s.require("simplex: x₁ + x₂ ≤ 1 and ≥ 3 is infeasible", simplex(lp).status == LpStatus::infeasible);
  }
  {
    LinearProgram lp;
    lp.c = {-1, -1};
    lp.a_ub = Rectangular::from_rows({{1, -1}});
    lp.b_ub = {1};
    const auto r = simplex(lp);
    s.require("simplex: maximising x₁ + x₂ with x₁ − x₂ ≤ 1 is unbounded", r.status == LpStatus::unbounded && r.value == -inf);
  }

  // Bounds, a free variable and a redundant equality: minimise x₁ − x₂ + x₃ with
  // 1 ≤ x₁ ≤ 3, x₂ ≤ 2 (no lower bound), x₃ free, x₁ + x₃ = 4 twice over (once scaled), x₃ ≥ −1 as a row.
  {
    LinearProgram lp;
    lp.c = {1, -1, 1};
    lp.a_eq = Rectangular::from_rows({{1, 0, 1}, {2, 0, 2}});
    lp.b_eq = {4, 8};
    lp.a_ub = Rectangular::from_rows({{0, 0, -1}});
    lp.b_ub = {1};
    lp.lower = {1, -inf, -inf};
    lp.upper = {3, 2, inf};
    const auto r = simplex(lp);
    // Value = x₁ + x₃ − x₂ = 4 − x₂, least at x₂ = 2: 2, with x₁ + x₃ = 4, x₃ ≥ −1.
    s.require("simplex: bounds, a free variable and a redundant equality: optimal", r.status == LpStatus::optimal);
    s.check("simplex: its optimum 2", r.value - 2, 1e-12);
    s.check("simplex: x₂ at its upper bound 2", r.x[1] - 2, 1e-12);
    s.check("simplex: and the equality held", r.x[0] + r.x[2] - 4, 1e-12);
    s.require("simplex: and the bounds", r.x[0] >= 1 - 1e-12 && r.x[0] <= 3 + 1e-12 && r.x[2] >= -1 - 1e-12);
  }

  // Zero-sum games: the row player maximises v subject to Aᵀp ≥ v 1, Σp = 1, p ≥ 0, v free.
  auto game = [](const std::vector<std::vector<double>>& A) {
    const std::size_t m = A.size(), n = A[0].size();
    LinearProgram lp;
    lp.c.assign(m + 1, 0.0);
    lp.c[m] = -1;  // maximise v
    lp.a_ub = Rectangular(n, m + 1);
    lp.b_ub.assign(n, 0.0);
    for (std::size_t j = 0; j < n; j++) {
      for (std::size_t i = 0; i < m; i++) lp.a_ub(j, i) = -A[i][j];
      lp.a_ub(j, m) = 1;  // v − (Aᵀp)_j ≤ 0
    }
    lp.a_eq = Rectangular(1, m + 1);
    for (std::size_t i = 0; i < m; i++) lp.a_eq(0, i) = 1;
    lp.b_eq = {1};
    lp.lower.assign(m + 1, 0.0);
    lp.upper.assign(m + 1, inf);
    lp.lower[m] = -inf;
    return simplex(lp);
  };
  {
    const auto r = game({{0, -1, 1}, {1, 0, -1}, {-1, 1, 0}});
    s.check("simplex: the value of rock–paper–scissors is 0", r.x[3], 1e-12);
    s.check("simplex: with the uniform strategy", std::fabs(r.x[0] - 1.0 / 3) + std::fabs(r.x[1] - 1.0 / 3) + std::fabs(r.x[2] - 1.0 / 3), 1e-12);
    // The column player's optimal strategy is the dual of the row constraints: −duals_ub.
    s.check("simplex: and the column player's strategy, the duals, uniform too",
            std::fabs(-r.duals_ub[0] - 1.0 / 3) + std::fabs(-r.duals_ub[1] - 1.0 / 3) + std::fabs(-r.duals_ub[2] - 1.0 / 3), 1e-12);
  }
  {
    const auto r = game({{3, -1}, {-2, 1}});
    s.check("simplex: a 2 × 2 game with value 1/7", r.x[2] - 1.0 / 7, 1e-12);
    s.check("simplex: at the mixed strategy (3/7, 4/7)", std::fabs(r.x[0] - 3.0 / 7) + std::fabs(r.x[1] - 4.0 / 7), 1e-12);
  }
  return s.done();
}
