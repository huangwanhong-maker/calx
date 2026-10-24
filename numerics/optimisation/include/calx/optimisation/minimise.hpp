// Minima of a function of one variable on an interval, without derivatives.
//
// Brent's method (Brent 1973, Algorithms for Minimization without Derivatives, ch. 5, the
// procedure localmin; in the form of Press et al. 2007, Numerical Recipes, 3rd ed., §10.3):
// golden-section steps, which shrink the interval by 0.618 whatever the function, and parabolic
// interpolation through the three best points where it lands inside and moves less than half the
// step before last, so that it converges superlinearly on smooth functions near a minimum and
// never more slowly than golden section. The minimum is a local one: on a unimodal function, the
// minimum on [a, b]; on others, one of them. It is found to within 2 eps |x| + tolerance (an
// interval of half-width tol1 = sqrt(eps) |x| + tolerance / 3 at the end, Brent's eq. 5.2.5 with
// his eps the square root of the machine's; asking for less than sqrt(eps) |x| is asking for
// noise, since f is flat to second order there).
//
// golden_section() is the plain golden-section search, for functions too rough for parabolas and
// as the reference the checks compare the counts of evaluations with.

#pragma once

#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace calx::optimisation {

struct Minimum {
  double x = 0, f = 0;
  int iterations = 0, evaluations = 0;
  bool converged = false;
};

template <class F>
inline Minimum brent_minimise(F&& f, double a, double b, double tolerance = 1e-10, int max_iterations = 500) {
  if (!(a < b)) throw std::invalid_argument("calx::optimisation::brent_minimise: an interval a < b");
  if (!(tolerance > 0)) throw std::invalid_argument("calx::optimisation::brent_minimise: a positive tolerance");
  constexpr double golden = 0.3819660112501051;  // (3 − √5) / 2
  const double root_eps = std::sqrt(std::numeric_limits<double>::epsilon());
  // x: the best point so far; w: the second best; v: the previous w. d: the last step; e: the one before.
  double x = a + golden * (b - a), w = x, v = x;
  double fx = f(x), fw = fx, fv = fx;
  double d = 0, e = 0;
  Minimum m;
  m.evaluations = 1;
  for (int it = 1; it <= max_iterations; it++) {
    const double mid = 0.5 * (a + b);
    const double tol1 = root_eps * std::fabs(x) + tolerance / 3, tol2 = 2 * tol1;
    if (std::fabs(x - mid) <= tol2 - 0.5 * (b - a)) {
      m.x = x;
      m.f = fx;
      m.iterations = it - 1;
      m.converged = true;
      return m;
    }
    bool golden_step = true;
    if (std::fabs(e) > tol1) {
      // The parabola through (v, fv), (w, fw), (x, fx): its vertex at x + p / q.
      const double r = (x - w) * (fx - fv);
      double q = (x - v) * (fx - fw);
      double p = (x - v) * q - (x - w) * r;
      q = 2 * (q - r);
      if (q > 0) p = -p;
      q = std::fabs(q);
      const double before = e;
      e = d;
      // Accepted when it lands inside (a, b) and moves less than half the step before last.
      if (std::fabs(p) < std::fabs(0.5 * q * before) && p > q * (a - x) && p < q * (b - x)) {
        d = p / q;
        const double u = x + d;
        if (u - a < tol2 || b - u < tol2) d = std::copysign(tol1, mid - x);
        golden_step = false;
      }
    }
    if (golden_step) {
      e = x >= mid ? a - x : b - x;
      d = golden * e;
    }
    // Never a step shorter than tol1: f would not tell the points apart.
    const double u = std::fabs(d) >= tol1 ? x + d : x + std::copysign(tol1, d);
    const double fu = f(u);
    m.evaluations++;
    if (fu <= fx) {
      (u >= x ? a : b) = x;
      v = w; fv = fw;
      w = x; fw = fx;
      x = u; fx = fu;
    } else {
      (u < x ? a : b) = u;
      if (fu <= fw || w == x) {
        v = w; fv = fw;
        w = u; fw = fu;
      } else if (fu <= fv || v == x || v == w) {
        v = u; fv = fu;
      }
    }
  }
  m.x = x;
  m.f = fx;
  m.iterations = max_iterations;
  return m;
}

// The maximum, as the minimum of −f; Minimum::f is f's value there.
template <class F>
inline Minimum brent_maximise(F&& f, double a, double b, double tolerance = 1e-10, int max_iterations = 500) {
  Minimum m = brent_minimise([&](double t) { return -f(t); }, a, b, tolerance, max_iterations);
  m.f = -m.f;
  return m;
}

// Golden-section search to an interval of width `tolerance` around a local minimum.
template <class F>
inline Minimum golden_section(F&& f, double a, double b, double tolerance = 1e-10, int max_iterations = 500) {
  if (!(a < b)) throw std::invalid_argument("calx::optimisation::golden_section: an interval a < b");
  constexpr double g = 0.6180339887498949;  // (√5 − 1) / 2
  double c = b - g * (b - a), d = a + g * (b - a);
  double fc = f(c), fd = f(d);
  Minimum m;
  m.evaluations = 2;
  int it = 0;
  while (b - a > tolerance && it < max_iterations) {
    it++;
    if (fc <= fd) {
      b = d; d = c; fd = fc;
      c = b - g * (b - a);
      fc = f(c);
    } else {
      a = c; c = d; fc = fd;
      d = a + g * (b - a);
      fd = f(d);
    }
    m.evaluations++;
  }
  m.converged = b - a <= tolerance;
  m.iterations = it;
  m.x = fc <= fd ? c : d;
  m.f = std::min(fc, fd);
  return m;
}

}  // namespace calx::optimisation
