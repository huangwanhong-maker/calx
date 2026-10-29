// Linear programmes by the two-phase simplex method, with Bland's rule against cycling.
//
//   minimise cᵀx  subject to  A_ub x ≤ b_ub,  A_eq x = b_eq,  lower ≤ x ≤ upper,
//
// with the bounds 0 ≤ x < ∞ by default (any bound may be ±∞). The problem is brought to the
// standard form min c̃ᵀz, Ãz = b̃ ≥ 0, z ≥ 0: a variable with a finite lower bound is shifted to
// it, one with only an upper bound is reflected, a free one is split into two non-negative
// parts, a finite upper bound becomes a row; each inequality gains a slack, and a row with a
// negative right-hand side is negated. Phase I minimises the sum of one artificial variable per
// row from the artificial basis; a positive minimum proves the problem infeasible. Artificials
// left in the basis at zero are pivoted out, and a row in which no other column can replace
// one is redundant and dropped. Phase II minimises cᵀx from the feasible basis so found. The
// tableau is dense and updated by Gauss–Jordan pivots (Dantzig 1963, Linear Programming and
// Extensions; Chvátal 1983, Linear Programming, ch. 2–3 and 8).
//
// Pivoting: the entering column is the one of least index with a negative reduced cost, and the
// leaving row the one of least ratio, ties going to the basic variable of least index: Bland's
// rule (Bland 1977, "New finite pivoting rules for the simplex method", Math. Oper. Res. 2,
// 103–107; Chvátal 1983, Theorem 3.3), under which the method cannot cycle and so terminates on
// degenerate problems too. It is slower than Dantzig's largest-coefficient rule on large
// problems, and this solver serves the small and moderate programmes of statistics, games and
// allocation (hundreds of rows and columns); large sparse programmes want a revised or
// interior-point method (planned).
//
// Duals. For the optimal basis B of the standard form, y = B⁻ᵀ c_B (read from the reduced costs
// of the artificial columns, −yᵢ); they are mapped back to the original rows as the marginals
//   duals_ub[i] = ∂(optimal value)/∂b_ub[i] ≤ 0,   duals_eq[i] = ∂(optimal value)/∂b_eq[i],
// the convention of scipy.optimize.linprog's marginals. They satisfy complementary slackness,
// duals_ub[i] (b_ub[i] − a_iᵀx) = 0, and, for the default bounds x ≥ 0, strong duality
// cᵀx = b_ubᵀ duals_ub + b_eqᵀ duals_eq with c − A_ubᵀ duals_ub − A_eqᵀ duals_eq ≥ 0 (Chvátal
// 1983, ch. 5). At a degenerate optimum the duals are one of several valid choices.
//
// Tolerances: entries below `tolerance` (default 10⁻⁹, relative to the scale of the data) count
// as zero in pivoting; the phase-I minimum counts as zero below tolerance (1 + ‖b‖∞).

#pragma once

#include <calx/linalg/rectangular.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace calx::optimisation {

struct LinearProgram {
  std::vector<double> c;                 // the objective, minimised; one entry per variable
  calx::linalg::Rectangular a_ub;        // rows of A_ub x ≤ b_ub (0 rows: none)
  std::vector<double> b_ub;
  calx::linalg::Rectangular a_eq;        // rows of A_eq x = b_eq (0 rows: none)
  std::vector<double> b_eq;
  std::vector<double> lower, upper;      // empty: 0 and +∞; else one per variable, ±∞ allowed
};

enum class LpStatus { optimal, infeasible, unbounded, iteration_limit };

struct LpOptions {
  double tolerance = 1e-9;
  int max_iterations = 0;                // 0: 50 (rows + columns of the standard form) + 1000
};

struct LpResult {
  LpStatus status = LpStatus::infeasible;
  std::vector<double> x;                 // the optimal point (when optimal)
  double value = 0;                      // cᵀx (when optimal; −∞ when unbounded)
  std::vector<double> duals_ub, duals_eq;  // the marginals ∂value/∂b (when optimal)
  int iterations = 0;                    // pivots of both phases
};

inline LpResult simplex(const LinearProgram& lp, const LpOptions& options = {}) {
  constexpr double inf = std::numeric_limits<double>::infinity();
  const std::size_t n = lp.c.size();
  const std::size_t mu = lp.b_ub.size(), me = lp.b_eq.size();
  if (lp.a_ub.rows != mu || (mu > 0 && lp.a_ub.cols != n)) throw std::invalid_argument("calx::optimisation::simplex: A_ub is |b_ub| × |c|");
  if (lp.a_eq.rows != me || (me > 0 && lp.a_eq.cols != n)) throw std::invalid_argument("calx::optimisation::simplex: A_eq is |b_eq| × |c|");
  if (!lp.lower.empty() && lp.lower.size() != n) throw std::invalid_argument("calx::optimisation::simplex: one lower bound per variable");
  if (!lp.upper.empty() && lp.upper.size() != n) throw std::invalid_argument("calx::optimisation::simplex: one upper bound per variable");
  const double tol = options.tolerance;

  // ---- the standard columns: x_j = offset_j + Σ sign · z_col ----
  struct Piece { std::size_t var; double sign; };
  std::vector<Piece> pieces;
  std::vector<double> offset(n, 0.0);
  std::vector<std::pair<std::size_t, double>> bound_rows;  // (standard column, upper limit of z)
  for (std::size_t j = 0; j < n; j++) {
    const double lo = lp.lower.empty() ? 0.0 : lp.lower[j];
    const double hi = lp.upper.empty() ? inf : lp.upper[j];
    if (lo > hi) {
      LpResult r;
      r.status = LpStatus::infeasible;
      return r;
    }
    if (std::isfinite(lo)) {
      offset[j] = lo;
      pieces.push_back({j, 1});
      if (std::isfinite(hi)) bound_rows.push_back({pieces.size() - 1, hi - lo});
    } else if (std::isfinite(hi)) {
      offset[j] = hi;
      pieces.push_back({j, -1});
    } else {
      pieces.push_back({j, 1});
      pieces.push_back({j, -1});
    }
  }
  const std::size_t nz = pieces.size();
  const std::size_t mb = bound_rows.size();
  const std::size_t m = mu + mb + me;     // rows: inequalities, bounds, equalities
  const std::size_t ns = mu + mb;         // slacks
  // Columns: z (nz) | slacks (ns) | artificials (m) | rhs.
  const std::size_t col_slack = nz, col_art = nz + ns, width = nz + ns + m + 1, rhs = width - 1;
  std::vector<double> T(m * width, 0.0);
  auto at = [&](std::size_t i, std::size_t j) -> double& { return T[i * width + j]; };
  std::vector<double> flip(m, 1.0);
  double bscale = 0;
  for (std::size_t i = 0; i < m; i++) {
    double b;
    if (i < mu) {
      b = lp.b_ub[i];
      for (std::size_t k = 0; k < nz; k++) at(i, k) = pieces[k].sign * lp.a_ub(i, pieces[k].var);
      for (std::size_t j = 0; j < n; j++) b -= lp.a_ub(i, j) * offset[j];
      at(i, col_slack + i) = 1;
    } else if (i < mu + mb) {
      const auto& br = bound_rows[i - mu];
      at(i, br.first) = 1;
      b = br.second;
      at(i, col_slack + i) = 1;
    } else {
      const std::size_t e = i - mu - mb;
      b = lp.b_eq[e];
      for (std::size_t k = 0; k < nz; k++) at(i, k) = pieces[k].sign * lp.a_eq(e, pieces[k].var);
      for (std::size_t j = 0; j < n; j++) b -= lp.a_eq(e, j) * offset[j];
    }
    at(i, rhs) = b;
    if (b < 0) {
      flip[i] = -1;
      for (std::size_t j = 0; j < col_art; j++) at(i, j) = -at(i, j);
      at(i, rhs) = -b;
    }
    at(i, col_art + i) = 1;
    bscale = std::max(bscale, std::fabs(b));
  }
  std::vector<std::size_t> basis(m);
  for (std::size_t i = 0; i < m; i++) basis[i] = col_art + i;
  std::vector<bool> alive(m, true);  // rows not dropped as redundant

  LpResult out;
  const int max_iterations = options.max_iterations > 0 ? options.max_iterations : static_cast<int>(50 * (m + width) + 1000);

  std::vector<double> z(width, 0.0);  // the reduced costs, z[rhs] = −(objective value)
  auto price = [&](const std::vector<double>& cost) {
    std::fill(z.begin(), z.end(), 0.0);
    for (std::size_t j = 0; j < rhs; j++) z[j] = cost[j];
    for (std::size_t i = 0; i < m; i++) {
      if (!alive[i]) continue;
      const double cb = cost[basis[i]];
      if (cb == 0) continue;
      for (std::size_t j = 0; j < width; j++) z[j] -= cb * at(i, j);
    }
  };
  auto pivot = [&](std::size_t r, std::size_t c) {
    const double pv = at(r, c);
    for (std::size_t j = 0; j < width; j++) at(r, j) /= pv;
    at(r, c) = 1;
    for (std::size_t i = 0; i < m; i++) {
      if (i == r || !alive[i]) continue;
      const double f = at(i, c);
      if (f == 0) continue;
      for (std::size_t j = 0; j < width; j++) at(i, j) -= f * at(r, j);
      at(i, c) = 0;
    }
    const double f = z[c];
    if (f != 0) {
      for (std::size_t j = 0; j < width; j++) z[j] -= f * at(r, j);
      z[c] = 0;
    }
    basis[r] = c;
    out.iterations++;
  };
  // Runs the simplex on columns [0, last) until optimal (true) or unbounded (false).
  auto run = [&](std::size_t last, bool& unbounded) {
    unbounded = false;
    while (out.iterations < max_iterations) {
      std::size_t enter = width;
      for (std::size_t j = 0; j < last; j++)
        if (z[j] < -tol) {
          enter = j;
          break;
        }
      if (enter == width) return true;
      std::size_t leave = m;
      double best = inf;
      for (std::size_t i = 0; i < m; i++) {
        if (!alive[i] || at(i, enter) <= tol) continue;
        const double ratio = at(i, rhs) / at(i, enter);
        const double tie = tol * std::max(1.0, std::fabs(ratio));
        if (leave == m || ratio < best - tie) {
          best = ratio;
          leave = i;
        } else if (std::fabs(ratio - best) <= tie && basis[i] < basis[leave]) {
          leave = i;  // Bland: among equal ratios, the basic variable of least index leaves
        }
      }
      if (leave == m) {
        unbounded = true;
        return false;
      }
      pivot(leave, enter);
    }
    return false;
  };

  // ---- phase I ----
  std::vector<double> cost1(rhs, 0.0);
  for (std::size_t i = 0; i < m; i++) cost1[col_art + i] = 1;
  price(cost1);
  bool unbounded = false;
  if (!run(rhs, unbounded)) {
    out.status = LpStatus::iteration_limit;
    return out;
  }
  if (-z[rhs] > tol * (1 + bscale)) {
    out.status = LpStatus::infeasible;
    return out;
  }
  // Pivot the artificials out of the basis, or drop their rows as redundant.
  for (std::size_t i = 0; i < m; i++) {
    if (basis[i] < col_art) continue;
    std::size_t c = col_art;
    for (std::size_t j = 0; j < col_art; j++)
      if (std::fabs(at(i, j)) > tol) {
        c = j;
        break;
      }
    if (c < col_art) pivot(i, c);
    else alive[i] = false;
  }

  // ---- phase II ----
  std::vector<double> cost2(rhs, 0.0);
  for (std::size_t k = 0; k < nz; k++) cost2[k] = pieces[k].sign * lp.c[pieces[k].var];
  price(cost2);
  if (!run(col_art, unbounded)) {
    out.status = unbounded ? LpStatus::unbounded : LpStatus::iteration_limit;
    if (unbounded) out.value = -inf;
    return out;
  }

  // ---- the solution and the duals ----
  std::vector<double> zval(col_art, 0.0);
  for (std::size_t i = 0; i < m; i++)
    if (alive[i] && basis[i] < col_art) zval[basis[i]] = at(i, rhs);
  out.x = offset;
  for (std::size_t k = 0; k < nz; k++) out.x[pieces[k].var] += pieces[k].sign * zval[k];
  out.value = 0;
  for (std::size_t j = 0; j < n; j++) out.value += lp.c[j] * out.x[j];
  out.duals_ub.assign(mu, 0.0);
  out.duals_eq.assign(me, 0.0);
  // y = −(reduced costs of the artificials), for every row, a dropped one included: the tableau's
  // rows are combinations M A of the original rows, the artificial columns hold M, and
  // y = Mᵀ c_B solves yᵀ A_B = c_Bᵀ on all rows.
  for (std::size_t i = 0; i < m; i++) {
    const double y = -z[col_art + i] * flip[i];
    if (i < mu) out.duals_ub[i] = y;
    else if (i >= mu + mb) out.duals_eq[i - mu - mb] = y;
  }
  out.status = LpStatus::optimal;
  return out;
}

}  // namespace calx::optimisation
