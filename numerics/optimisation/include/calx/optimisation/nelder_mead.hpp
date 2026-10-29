// Minima of a function of several variables without derivatives: the Nelder–Mead simplex method.
//
// Nelder and Mead (1965), Computer Journal 7, 308–313, in the precise form of Lagarias, Reeds,
// Wright and Wright (1998), "Convergence properties of the Nelder–Mead simplex method in low
// dimensions", SIAM J. Optim. 9, 112–147, §2 (Algorithm NM): the n + 1 vertices are kept
// ordered by value (a new vertex is placed after those of equal value, their non-shrink
// tie-breaking rule); each iteration reflects the worst vertex through the centroid of the
// others (coefficient ρ), expands (χ), contracts outside or inside (γ), or shrinks towards the
// best (σ). The standard coefficients are ρ = 1, χ = 2, γ = 1/2, σ = 1/2; the adaptive ones of
// Gao and Han (2012), "Implementing the Nelder–Mead simplex algorithm with adaptive
// parameters", Comput. Optim. Appl. 51, 259–277: ρ = 1, χ = 1 + 2/n,
// γ = 3/4 − 1/(2n), σ = 1 − 1/n, keep the method effective in higher dimensions.
//
// The initial simplex is x₀ and x₀ + δ_i e_i, with δ_i = initial_step[i] when given (one value
// for all, or one per coordinate), else 5% of x₀_i, or 0.00025 where x₀_i = 0 (the convention of
// MATLAB's fminsearch). The method stops, converged, when the vertices' values lie within
// f_tolerance of the best and the vertices within x_tolerance of the best in the max norm
// (fminsearch's test), or unconverged after max_evaluations evaluations.
//
// Limits. The method has no general convergence theory beyond one dimension (Lagarias et al.
// prove convergence for strictly convex functions in dimension 1 and to a minimiser in some
// cases in dimension 2; McKinnon 1998, SIAM J. Optim. 9, 148–158, gives a smooth convex function
// in two dimensions on which it converges to a non-stationary point). It suits rough or noisy
// objectives of a few variables, and polishing; with derivatives, quasi_newton.hpp converges far
// faster. A minimum is found to about √ε relative in x, at best, as with every method using
// only values. The function is called with a std::span<const double> of length n.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <span>
#include <stdexcept>
#include <vector>

namespace calx::optimisation {

struct NelderMeadOptions {
  double x_tolerance = 1e-8;            // max-norm spread of the vertices around the best
  double f_tolerance = 1e-10;           // spread of the values around the best
  int max_evaluations = 0;              // 0: max(1000, 200 n²)
  std::vector<double> initial_step;     // empty: fminsearch's 5% (0.00025 at zero); one value: for all; n values
  bool adaptive = false;                // Gao and Han's coefficients
};

struct NelderMeadResult {
  std::vector<double> x;
  double value = 0;
  int iterations = 0, evaluations = 0;
  bool converged = false;
};

template <class F>
inline NelderMeadResult nelder_mead(F&& f, std::vector<double> x0, const NelderMeadOptions& options = {}) {
  const std::size_t n = x0.size();
  if (n == 0) throw std::invalid_argument("calx::optimisation::nelder_mead: at least one variable");
  const double dn = static_cast<double>(n);
  const double rho = 1, chi = options.adaptive ? 1 + 2 / dn : 2, gamma = options.adaptive ? 0.75 - 1 / (2 * dn) : 0.5,
               sigma = options.adaptive ? 1 - 1 / dn : 0.5;
  const int max_evaluations = options.max_evaluations > 0 ? options.max_evaluations : std::max(1000, static_cast<int>(200 * n * n));
  NelderMeadResult out;
  auto eval = [&](const std::vector<double>& x) {
    out.evaluations++;
    const double v = f(std::span<const double>(x.data(), x.size()));
    return std::isnan(v) ? std::numeric_limits<double>::infinity() : v;
  };

  // The simplex: vertices v[order[0]] (best) ... v[order[n]] (worst).
  std::vector<std::vector<double>> v(n + 1, x0);
  std::vector<double> fv(n + 1);
  for (std::size_t i = 0; i < n; i++) {
    double step;
    if (options.initial_step.size() == 1) step = options.initial_step[0];
    else if (options.initial_step.size() == n) step = options.initial_step[i];
    else if (!options.initial_step.empty()) throw std::invalid_argument("calx::optimisation::nelder_mead: initial_step has 0, 1 or n values");
    else step = x0[i] != 0 ? 0.05 * x0[i] : 0.00025;
    v[i + 1][i] += step;
  }
  for (std::size_t i = 0; i <= n; i++) fv[i] = eval(v[i]);
  std::vector<std::size_t> order(n + 1);
  std::iota(order.begin(), order.end(), std::size_t{0});
  std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return fv[a] < fv[b]; });

  // Places vertex `slot` (just replaced) after every vertex of value ≤ its own, among the others.
  auto reinsert_worst = [&]() {
    const std::size_t w = order[n];
    std::size_t pos = n;
    while (pos > 0 && fv[order[pos - 1]] > fv[w]) {
      order[pos] = order[pos - 1];
      pos--;
    }
    order[pos] = w;
  };

  std::vector<double> centroid(n), xr(n), xe(n), xc(n);
  auto point = [&](std::vector<double>& out_x, double t) {  // centroid + t (centroid − worst)
    const auto& worst = v[order[n]];
    for (std::size_t j = 0; j < n; j++) out_x[j] = centroid[j] + t * (centroid[j] - worst[j]);
  };

  while (true) {
    // The test of fminsearch.
    const auto& best = v[order[0]];
    double fspread = 0, xspread = 0;
    for (std::size_t i = 1; i <= n; i++) {
      fspread = std::max(fspread, std::fabs(fv[order[i]] - fv[order[0]]));
      for (std::size_t j = 0; j < n; j++) xspread = std::max(xspread, std::fabs(v[order[i]][j] - best[j]));
    }
    if (fspread <= options.f_tolerance && xspread <= options.x_tolerance) {
      out.converged = true;
      break;
    }
    if (out.evaluations >= max_evaluations) break;
    out.iterations++;

    std::fill(centroid.begin(), centroid.end(), 0.0);
    for (std::size_t i = 0; i < n; i++)
      for (std::size_t j = 0; j < n; j++) centroid[j] += v[order[i]][j];
    for (double& c : centroid) c /= dn;

    const double f1 = fv[order[0]], fn = fv[order[n - 1]], fworst = fv[order[n]];
    point(xr, rho);
    const double fr = eval(xr);
    bool shrink = false;
    if (fr < f1) {  // expand
      point(xe, rho * chi);
      const double fe = eval(xe);
      if (fe < fr) {
        v[order[n]] = xe;
        fv[order[n]] = fe;
      } else {
        v[order[n]] = xr;
        fv[order[n]] = fr;
      }
    } else if (fr < fn) {  // accept the reflection
      v[order[n]] = xr;
      fv[order[n]] = fr;
    } else if (fr < fworst) {  // outside contraction
      point(xc, rho * gamma);
      const double fc = eval(xc);
      if (fc <= fr) {
        v[order[n]] = xc;
        fv[order[n]] = fc;
      } else {
        shrink = true;
      }
    } else {  // inside contraction
      point(xc, -gamma);
      const double fc = eval(xc);
      if (fc < fworst) {
        v[order[n]] = xc;
        fv[order[n]] = fc;
      } else {
        shrink = true;
      }
    }
    if (shrink) {
      const auto b = v[order[0]];
      for (std::size_t i = 1; i <= n; i++) {
        auto& x = v[order[i]];
        for (std::size_t j = 0; j < n; j++) x[j] = b[j] + sigma * (x[j] - b[j]);
        fv[order[i]] = eval(x);
      }
      std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t c) { return fv[a] < fv[c]; });
    } else {
      reinsert_worst();
    }
  }
  out.x = v[order[0]];
  out.value = fv[order[0]];
  return out;
}

}  // namespace calx::optimisation
