// Roots of a function of one variable, bracketed.
//
// Brent's method (Brent 1973, Algorithms for Minimization without Derivatives, ch. 4, the
// procedure zero; in the form of Press et al. 2007, Numerical Recipes, 3rd ed., §9.3):
// inverse quadratic interpolation where it makes progress, bisection where it does not,
// so that it converges superlinearly on smooth functions and never more slowly than
// bisection. f(a) and f(b) must differ in sign; the root is found to within
// 2 eps |x| + tolerance / 2.

#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace calx::solvers {

struct Root {
  double x = 0;
  int iterations = 0;
  bool converged = false;
};

template <class F>
inline Root brent(F&& f, double a, double b, double tolerance = 1e-14, int max_iterations = 200) {
  constexpr double eps = std::numeric_limits<double>::epsilon();
  double fa = f(a), fb = f(b);
  if ((fa > 0 && fb > 0) || (fa < 0 && fb < 0)) throw std::invalid_argument("calx::solvers::brent: the root is not bracketed");
  double c = b, fc = fb, d = b - a, e = d;
  for (int it = 1; it <= max_iterations; it++) {
    if ((fb > 0 && fc > 0) || (fb < 0 && fc < 0)) {
      c = a; fc = fa; d = b - a; e = d;
    }
    if (std::fabs(fc) < std::fabs(fb)) {
      a = b; b = c; c = a;
      fa = fb; fb = fc; fc = fa;
    }
    const double tol1 = 2 * eps * std::fabs(b) + 0.5 * tolerance;
    const double xm = 0.5 * (c - b);
    if (std::fabs(xm) <= tol1 || fb == 0) return {b, it, true};
    if (std::fabs(e) >= tol1 && std::fabs(fa) > std::fabs(fb)) {
      const double s = fb / fa;
      double p, q;
      if (a == c) {
        p = 2 * xm * s;
        q = 1 - s;
      } else {
        const double qq = fa / fc, r = fb / fc;
        p = s * (2 * xm * qq * (qq - r) - (b - a) * (r - 1));
        q = (qq - 1) * (r - 1) * (s - 1);
      }
      if (p > 0) q = -q;
      p = std::fabs(p);
      if (2 * p < std::min(3 * xm * q - std::fabs(tol1 * q), std::fabs(e * q))) {
        e = d;
        d = p / q;
      } else {
        d = xm;
        e = d;
      }
    } else {
      d = xm;
      e = d;
    }
    a = b;
    fa = fb;
    b += std::fabs(d) > tol1 ? d : std::copysign(tol1, xm);
    fb = f(b);
  }
  return {b, max_iterations, false};
}

}  // namespace calx::solvers
