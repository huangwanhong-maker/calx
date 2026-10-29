// Minima of a smooth function of several variables: the BFGS quasi-Newton method with a line
// search satisfying the strong Wolfe conditions.
//
// BFGS (Broyden, Fletcher, Goldfarb and Shanno, 1970) as Nocedal and Wright (2006), Numerical
// Optimization, 2nd ed., Algorithm 6.1, give it: the search direction p = −H ∇f with H an
// approximation of the inverse Hessian, updated after each step s = x₊ − x, y = ∇f₊ − ∇f by
//   H₊ = (I − ρ s yᵀ) H (I − ρ y sᵀ) + ρ s sᵀ,  ρ = 1 / yᵀs          (N&W eq. 6.17),
// starting from H₀ = (yᵀs / yᵀy) I, set after the first step and before the first update
// (eq. 6.20). The step length comes from the line search of N&W Algorithm 3.5 with zoom
// (Algorithm 3.6), which returns a step satisfying the strong Wolfe conditions
//   f(x + αp) ≤ f(x) + c₁ α ∇fᵀp,   |∇f(x + αp)ᵀp| ≤ c₂ |∇fᵀp|      (eq. 3.7),
// c₁ = 10⁻⁴, c₂ = 0.9 by default, trying α = 1 first, doubling while the function keeps
// decreasing, and interpolating within the bracket by the minimiser of the cubic through the
// two ends' values and slopes (eq. 3.59), safeguarded to the middle 80% of the bracket. Under
// the Wolfe conditions yᵀs > 0, so H stays positive definite.
//
// Convergence. On a twice continuously differentiable function with a positive definite Hessian
// at the minimiser, BFGS converges superlinearly (N&W Theorem 6.6). It stops, converged, when
// ‖∇f‖∞ ≤ gradient_tolerance; unconverged when the line search fails to find a decrease (the
// function is flat to rounding along p, or the gradient is wrong) or after max_iterations.
//
// The inverse Hessian. On a strictly convex quadratic with exact line searches, BFGS ends in
// n steps with H = ∇²f⁻¹ exactly (N&W §6.1); a small c₂ makes the line search exact on quadratics (the cubic interpolant
// of a quadratic is itself). In general H approximates the inverse Hessian only along the
// directions explored, so standard errors of an estimate are better computed from the Hessian
// itself, by hessian_by_differences() below, than from BfgsResult::inverse_hessian.
//
// Derivatives. The gradient is given, as g(x) returning a std::vector<double>, or computed by
// central differences with steps h_i = ε^{1/3} max(|x_i|, 1), whose error is O(h²) truncation plus
// O(ε f / h) rounding, about ε^{2/3} relative (N&W §8.1). The function is called with
// a std::span<const double>.

#pragma once

#include <calx/linalg/dense.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <stdexcept>
#include <vector>

namespace calx::optimisation {

struct BfgsOptions {
  double gradient_tolerance = 1e-6;  // ‖∇f‖∞ at which it stops, converged
  int max_iterations = 0;            // 0: max(200, 50 n)
  double c1 = 1e-4, c2 = 0.9;        // the Wolfe constants, 0 < c1 < c2 < 1
  double max_step = 1e10;            // the largest α tried along a direction
};

struct BfgsResult {
  std::vector<double> x;
  double value = 0;
  double gradient_norm = 0;          // ‖∇f(x)‖∞
  int iterations = 0, evaluations = 0;  // evaluations of f (gradient calls counted apart, not here)
  bool converged = false;
  calx::linalg::Matrix inverse_hessian;  // the final BFGS approximation H
};

// The gradient of f at x by central differences, steps ε^{1/3} max(|x_i|, 1).
template <class F>
inline std::vector<double> gradient_by_differences(F&& f, std::span<const double> x) {
  const double h0 = std::cbrt(std::numeric_limits<double>::epsilon());
  std::vector<double> xx(x.begin(), x.end()), g(x.size());
  for (std::size_t i = 0; i < x.size(); i++) {
    const double h = h0 * std::max(std::fabs(x[i]), 1.0);
    xx[i] = x[i] + h;
    const double fp = f(std::span<const double>(xx));
    xx[i] = x[i] - h;
    const double fm = f(std::span<const double>(xx));
    xx[i] = x[i];
    g[i] = (fp - fm) / (2 * h);
  }
  return g;
}

// The Hessian of f at x by central differences of values (Abramowitz and Stegun 1964, §25.3),
// steps ε^{1/4} max(|x_i|, 1), error about ε^{1/2} relative; symmetric by
// construction. 2n² + 1 evaluations.
template <class F>
inline calx::linalg::Matrix hessian_by_differences(F&& f, std::span<const double> x) {
  const std::size_t n = x.size();
  const double h0 = std::pow(std::numeric_limits<double>::epsilon(), 0.25);
  std::vector<double> xx(x.begin(), x.end()), h(n);
  for (std::size_t i = 0; i < n; i++) h[i] = h0 * std::max(std::fabs(x[i]), 1.0);
  auto at = [&]() { return f(std::span<const double>(xx)); };
  const double f0 = at();
  calx::linalg::Matrix H(n);
  for (std::size_t i = 0; i < n; i++) {
    xx[i] = x[i] + h[i];
    const double fp = at();
    xx[i] = x[i] - h[i];
    const double fm = at();
    xx[i] = x[i];
    H(i, i) = (fp - 2 * f0 + fm) / (h[i] * h[i]);
    for (std::size_t j = 0; j < i; j++) {
      double s = 0;
      for (int a = -1; a <= 1; a += 2)
        for (int b = -1; b <= 1; b += 2) {
          xx[i] = x[i] + a * h[i];
          xx[j] = x[j] + b * h[j];
          s += a * b * at();
        }
      xx[i] = x[i];
      xx[j] = x[j];
      H(i, j) = H(j, i) = s / (4 * h[i] * h[j]);
    }
  }
  return H;
}

namespace detail {

// The minimiser of the cubic through (a, fa, ga) and (b, fb, gb) (Nocedal and Wright, eq. 3.59),
// or the midpoint when it does not exist or falls outside the middle 80% of [a, b].
inline double cubic_step(double a, double fa, double ga, double b, double fb, double gb) {
  const double lo = std::min(a, b), hi = std::max(a, b), width = hi - lo;
  const double d1 = ga + gb - 3 * (fa - fb) / (a - b);
  const double disc = d1 * d1 - ga * gb;
  double t = 0.5 * (a + b);
  if (disc >= 0 && b != a) {
    const double d2 = std::copysign(std::sqrt(disc), b - a);
    const double denom = gb - ga + 2 * d2;
    if (denom != 0) t = b - (b - a) * (gb + d2 - d1) / denom;
  }
  if (!(t >= lo + 0.1 * width && t <= hi - 0.1 * width)) t = 0.5 * (a + b);
  return t;
}

}  // namespace detail

template <class F, class G>
inline BfgsResult bfgs(F&& f, G&& gradient, std::vector<double> x0, const BfgsOptions& options = {}) {
  const std::size_t n = x0.size();
  if (n == 0) throw std::invalid_argument("calx::optimisation::bfgs: at least one variable");
  if (!(0 < options.c1 && options.c1 < options.c2 && options.c2 < 1)) throw std::invalid_argument("calx::optimisation::bfgs: 0 < c1 < c2 < 1");
  const int max_iterations = options.max_iterations > 0 ? options.max_iterations : std::max(200, static_cast<int>(50 * n));
  BfgsResult out;
  auto value = [&](const std::vector<double>& x) {
    out.evaluations++;
    return f(std::span<const double>(x.data(), x.size()));
  };
  auto grad = [&](const std::vector<double>& x) {
    std::vector<double> g = gradient(std::span<const double>(x.data(), x.size()));
    if (g.size() != n) throw std::invalid_argument("calx::optimisation::bfgs: the gradient has one entry per variable");
    return g;
  };
  auto dot = [](const std::vector<double>& a, const std::vector<double>& b) {
    double s = 0;
    for (std::size_t i = 0; i < a.size(); i++) s += a[i] * b[i];
    return s;
  };
  auto inf_norm = [](const std::vector<double>& a) {
    double m = 0;
    for (double v : a) m = std::max(m, std::fabs(v));
    return m;
  };

  std::vector<double> x = std::move(x0), g = grad(x), p(n), xt(n), gt(n);
  double fx = value(x);
  calx::linalg::Matrix H = calx::linalg::Matrix::identity(n);
  bool first = true;
  out.gradient_norm = inf_norm(g);

  for (int it = 0; it < max_iterations && out.gradient_norm > options.gradient_tolerance; it++) {
    for (std::size_t i = 0; i < n; i++) {
      double s = 0;
      for (std::size_t j = 0; j < n; j++) s -= H(i, j) * g[j];
      p[i] = s;
    }
    double slope0 = dot(g, p);
    if (!(slope0 < 0)) {  // not a descent direction (rounding in H): restart from steepest descent
      H = calx::linalg::Matrix::identity(n);
      for (std::size_t i = 0; i < n; i++) p[i] = -g[i];
      slope0 = dot(g, p);
      first = true;
    }
    // ---- the line search (N&W Algorithms 3.5 and 3.6) ----
    auto phi = [&](double alpha, double& slope) {
      for (std::size_t i = 0; i < n; i++) xt[i] = x[i] + alpha * p[i];
      const double v = value(xt);
      if (std::isfinite(v)) {
        gt = grad(xt);
        slope = dot(gt, p);
      } else {
        slope = std::numeric_limits<double>::quiet_NaN();
      }
      return v;
    };
    const double f0 = fx;
    double alpha_found = -1, f_found = 0;
    std::vector<double> g_found;
    auto zoom = [&](double lo, double flo, double glo, double hi, double fhi, double ghi) {
      for (int z = 0; z < 60; z++) {
        const double a = std::isfinite(fhi) && std::isfinite(ghi) ? detail::cubic_step(lo, flo, glo, hi, fhi, ghi) : 0.5 * (lo + hi);
        double ga;
        const double fa = phi(a, ga);
        if (!std::isfinite(fa) || fa > f0 + options.c1 * a * slope0 || fa >= flo) {
          hi = a;
          fhi = fa;
          ghi = ga;
        } else {
          if (std::fabs(ga) <= -options.c2 * slope0) {
            alpha_found = a;
            f_found = fa;
            g_found = gt;
            return;
          }
          if (ga * (hi - lo) >= 0) {
            hi = lo;
            fhi = flo;
            ghi = glo;
          }
          lo = a;
          flo = fa;
          glo = ga;
        }
        if (std::fabs(hi - lo) <= std::numeric_limits<double>::epsilon() * std::max(std::fabs(lo), 1.0)) break;
      }
      // The bracket has shrunk to rounding: accept its low end if it decreased f.
      if (lo > 0 && flo < f0) {
        double gl;
        f_found = phi(lo, gl);
        alpha_found = lo;
        g_found = gt;
      }
    };
    {
      double a_prev = 0, f_prev = f0, g_prev = slope0, a = 1;
      for (int i = 1; i <= 60; i++) {
        double ga;
        const double fa = phi(a, ga);
        if (!std::isfinite(fa) || fa > f0 + options.c1 * a * slope0 || (i > 1 && fa >= f_prev)) {
          zoom(a_prev, f_prev, g_prev, a, fa, ga);
          break;
        }
        if (std::fabs(ga) <= -options.c2 * slope0) {
          alpha_found = a;
          f_found = fa;
          g_found = gt;
          break;
        }
        if (ga >= 0) {
          zoom(a, fa, ga, a_prev, f_prev, g_prev);
          break;
        }
        a_prev = a;
        f_prev = fa;
        g_prev = ga;
        a = std::min(2 * a, options.max_step);
        if (a_prev >= options.max_step) break;
      }
    }
    if (alpha_found <= 0) break;  // no acceptable step: unconverged
    // ---- the step and the update (N&W Algorithm 6.1) ----
    std::vector<double> s(n), y(n);
    for (std::size_t i = 0; i < n; i++) {
      s[i] = alpha_found * p[i];
      y[i] = g_found[i] - g[i];
      x[i] += s[i];
    }
    fx = f_found;
    g = g_found;
    out.iterations = it + 1;
    out.gradient_norm = inf_norm(g);
    const double ys = dot(y, s);
    if (ys > 0) {
      if (first) {
        const double scale = ys / dot(y, y);
        H = calx::linalg::Matrix::identity(n);
        for (std::size_t i = 0; i < n; i++) H(i, i) = scale;
        first = false;
      }
      const double rho = 1 / ys;
      // H₊ = H − ρ (H y sᵀ + s yᵀ H) + (ρ² yᵀH y + ρ) s sᵀ, the expansion of eq. 6.17.
      std::vector<double> Hy(n, 0.0);
      for (std::size_t i = 0; i < n; i++)
        for (std::size_t j = 0; j < n; j++) Hy[i] += H(i, j) * y[j];
      const double yHy = dot(y, Hy);
      for (std::size_t i = 0; i < n; i++)
        for (std::size_t j = 0; j < n; j++) H(i, j) += -rho * (Hy[i] * s[j] + s[i] * Hy[j]) + (rho * rho * yHy + rho) * s[i] * s[j];
    }
  }
  out.converged = out.gradient_norm <= options.gradient_tolerance;
  out.x = std::move(x);
  out.value = fx;
  out.inverse_hessian = std::move(H);
  return out;
}

// BFGS with the gradient by central differences.
template <class F>
inline BfgsResult bfgs(F&& f, std::vector<double> x0, const BfgsOptions& options = {}) {
  return bfgs(f, [&](std::span<const double> x) { return gradient_by_differences(f, x); }, std::move(x0), options);
}

}  // namespace calx::optimisation
