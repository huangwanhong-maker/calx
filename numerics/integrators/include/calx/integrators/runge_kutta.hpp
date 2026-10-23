// One-step integrators of x' = f(x) over a step h, for states of fixed small dimension.
//
// - euler: forward Euler, order 1.
// - heun: the explicit trapezoidal rule, order 2.
// - rk4: the classical Runge-Kutta method, order 4.
// - dopri5: the embedded Runge-Kutta pair of Dormand and Prince (1980), order 5 with an
//   error estimate of order 4, taking as many substeps of [0, h] as its tolerances ask
//   (Hairer, Nørsett and Wanner 1993, II.4-5); the last accepted substep is returned so
//   that the next call may start from it, as NEST does with the GSL's solver.
//
// f is called as f(x, dx). Exponential integrators (Rush and Larsen 1978) need the
// structure of a model's gates and live with the models.
//
// Each method checks against closed-form solutions and its order of convergence
// (tests/check_integrators.cpp).

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>

namespace calx::integrators {

template <std::size_t N>
using Vec = std::array<double, N>;

template <std::size_t N, class F>
inline void euler(F&& f, Vec<N>& x, double h) {
  Vec<N> k;
  f(x, k);
  for (std::size_t i = 0; i < N; i++) x[i] += h * k[i];
}

template <std::size_t N, class F>
inline void heun(F&& f, Vec<N>& x, double h) {
  Vec<N> k1, k2, y;
  f(x, k1);
  for (std::size_t i = 0; i < N; i++) y[i] = x[i] + h * k1[i];
  f(y, k2);
  for (std::size_t i = 0; i < N; i++) x[i] += 0.5 * h * (k1[i] + k2[i]);
}

template <std::size_t N, class F>
inline void rk4(F&& f, Vec<N>& x, double h) {
  Vec<N> k1, k2, k3, k4, y;
  f(x, k1);
  for (std::size_t i = 0; i < N; i++) y[i] = x[i] + 0.5 * h * k1[i];
  f(y, k2);
  for (std::size_t i = 0; i < N; i++) y[i] = x[i] + 0.5 * h * k2[i];
  f(y, k3);
  for (std::size_t i = 0; i < N; i++) y[i] = x[i] + h * k3[i];
  f(y, k4);
  for (std::size_t i = 0; i < N; i++) x[i] += h / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
}

struct Tolerance {
  double absolute = 1e-6, relative = 1e-6;
};

// Integrates over [0, h]; `sub` is the substep to try first (0: h), and is set to the next
// one to try. Returns the number of accepted substeps.
template <std::size_t N, class F>
inline int dopri5(F&& f, Vec<N>& x, double h, const std::array<double, N>& scale, Tolerance tol, double& sub) {
  constexpr double a21 = 1.0 / 5;
  constexpr double a31 = 3.0 / 40, a32 = 9.0 / 40;
  constexpr double a41 = 44.0 / 45, a42 = -56.0 / 15, a43 = 32.0 / 9;
  constexpr double a51 = 19372.0 / 6561, a52 = -25360.0 / 2187, a53 = 64448.0 / 6561, a54 = -212.0 / 729;
  constexpr double a61 = 9017.0 / 3168, a62 = -355.0 / 33, a63 = 46732.0 / 5247, a64 = 49.0 / 176, a65 = -5103.0 / 18656;
  constexpr double b1 = 35.0 / 384, b3 = 500.0 / 1113, b4 = 125.0 / 192, b5 = -2187.0 / 6784, b6 = 11.0 / 84;
  constexpr double e1 = 71.0 / 57600, e3 = -71.0 / 16695, e4 = 71.0 / 1920, e5 = -17253.0 / 339200, e6 = 22.0 / 525, e7 = -1.0 / 40;
  double t = 0, s = (sub > 0 && sub <= h) ? sub : h;
  int accepted = 0;
  Vec<N> k1, k2, k3, k4, k5, k6, k7, y, z;
  f(x, k1);
  for (int guard = 0; t < h && guard < 100000; guard++) {
    if (t + s > h) s = h - t;
    for (std::size_t i = 0; i < N; i++) y[i] = x[i] + s * a21 * k1[i];
    f(y, k2);
    for (std::size_t i = 0; i < N; i++) y[i] = x[i] + s * (a31 * k1[i] + a32 * k2[i]);
    f(y, k3);
    for (std::size_t i = 0; i < N; i++) y[i] = x[i] + s * (a41 * k1[i] + a42 * k2[i] + a43 * k3[i]);
    f(y, k4);
    for (std::size_t i = 0; i < N; i++) y[i] = x[i] + s * (a51 * k1[i] + a52 * k2[i] + a53 * k3[i] + a54 * k4[i]);
    f(y, k5);
    for (std::size_t i = 0; i < N; i++) y[i] = x[i] + s * (a61 * k1[i] + a62 * k2[i] + a63 * k3[i] + a64 * k4[i] + a65 * k5[i]);
    f(y, k6);
    for (std::size_t i = 0; i < N; i++) z[i] = x[i] + s * (b1 * k1[i] + b3 * k3[i] + b4 * k4[i] + b5 * k5[i] + b6 * k6[i]);
    f(z, k7);
    double err = 0;
    for (std::size_t i = 0; i < N; i++) {
      const double e = s * (e1 * k1[i] + e3 * k3[i] + e4 * k4[i] + e5 * k5[i] + e6 * k6[i] + e7 * k7[i]);
      const double sc = scale[i] * tol.absolute + tol.relative * std::max(std::fabs(x[i]), std::fabs(z[i]));
      err = std::max(err, std::fabs(e) / sc);
    }
    if (err <= 1 || s <= h * 1e-9) {
      t += s;
      x = z;
      k1 = k7;  // first same as last
      accepted++;
      const double grow = err > 0 ? std::min(5.0, 0.9 * std::pow(err, -0.2)) : 5.0;
      if (t < h) s *= grow;
      else sub = std::min(h, s * grow);
    } else {
      s *= std::max(0.2, 0.9 * std::pow(err, -0.2));
    }
  }
  return accepted;
}

}  // namespace calx::integrators
