// Forward differentiation by truncated Taylor arithmetic: numbers that carry their own
// first and second derivatives.
//
// A jet of order K in N variables is the Taylor polynomial of a function at a point,
// truncated after degree K: Taylor<N, 2> holds the value f, the gradient ∂_i f and the
// Hessian ∂_i ∂_j f; Taylor<N, 1> (a dual number) holds the value and the gradient. The
// arithmetic of jets is the arithmetic of such polynomials modulo terms of degree K + 1,
// so a function written once as a template, and evaluated on the jets of the coordinates,
// returns its derivatives exact to rounding (Griewank and Walther 2008, Evaluating
// Derivatives, 2nd ed., SIAM, ch. 3 for the forward mode and ch. 13 for Taylor arithmetic;
// Rall 1981). No step size enters, unlike finite differences.
//
// Rules, for jets a and b and a function φ of one variable (Griewank and Walther, Table 13.1):
//   (ab)_i  = a b_i + b a_i,       (ab)_ij = a b_ij + b a_ij + a_i b_j + a_j b_i;
//   φ(a)_i  = φ'(a) a_i,           φ(a)_ij = φ'(a) a_ij + φ''(a) a_i a_j.
// Every elementary function below is the second rule with its own φ' and φ''; atan2 and the
// power of two jets use the same rule in two variables.
//
// A function to be differentiated is written as a template over its number type and calls
// the elementary functions unqualified, after `using std::sin;` and the like, so that the
// same code serves double and jets (argument-dependent lookup finds the overloads here).
// Branches on a jet compare its value. The cost of a product is O(N^2) for order 2.

#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <type_traits>

namespace calx::differentiation {

template <std::size_t N, unsigned K>
struct Taylor {
  static_assert(K == 1 || K == 2, "calx::differentiation: jets of order 1 or 2");
  struct None {};
  using Hessian = std::conditional_t<K == 2, std::array<std::array<double, N>, N>, None>;

  double v = 0;                 // the value
  std::array<double, N> d{};    // the gradient
  [[no_unique_address]] Hessian h{};  // the Hessian (order 2 only)

  Taylor() = default;
  Taylor(double c) : v(c) {}  // a constant
};

template <std::size_t N>
using Jet = Taylor<N, 2>;
template <std::size_t N>
using Dual = Taylor<N, 1>;

// The value of a number or a jet, for code that branches.
inline double value(double x) { return x; }
template <std::size_t N, unsigned K>
inline double value(const Taylor<N, K>& x) { return x.v; }

// The jets of the coordinates x_i: value x_i, gradient the i-th unit vector.
template <std::size_t N, unsigned K = 2>
inline std::array<Taylor<N, K>, N> variables(const std::array<double, N>& x) {
  std::array<Taylor<N, K>, N> z;
  for (std::size_t i = 0; i < N; i++) {
    z[i].v = x[i];
    z[i].d[i] = 1;
  }
  return z;
}

// φ(a), given φ(a.v), φ'(a.v) and φ''(a.v).
template <std::size_t N, unsigned K>
inline Taylor<N, K> chain(const Taylor<N, K>& a, double f0, double f1, double f2) {
  Taylor<N, K> r;
  r.v = f0;
  for (std::size_t i = 0; i < N; i++) r.d[i] = f1 * a.d[i];
  if constexpr (K == 2)
    for (std::size_t i = 0; i < N; i++)
      for (std::size_t j = 0; j < N; j++) r.h[i][j] = f1 * a.h[i][j] + f2 * a.d[i] * a.d[j];
  return r;
}

// φ(a, b), given φ and its first and second partial derivatives at (a.v, b.v).
template <std::size_t N, unsigned K>
inline Taylor<N, K> chain(const Taylor<N, K>& a, const Taylor<N, K>& b, double f0, double fa, double fb, double faa, double fab, double fbb) {
  Taylor<N, K> r;
  r.v = f0;
  for (std::size_t i = 0; i < N; i++) r.d[i] = fa * a.d[i] + fb * b.d[i];
  if constexpr (K == 2)
    for (std::size_t i = 0; i < N; i++)
      for (std::size_t j = 0; j < N; j++)
        r.h[i][j] = fa * a.h[i][j] + fb * b.h[i][j] + faa * a.d[i] * a.d[j] + fab * (a.d[i] * b.d[j] + b.d[i] * a.d[j]) + fbb * b.d[i] * b.d[j];
  return r;
}

// ---- arithmetic -------------------------------------------------------------------------------------------------------

template <std::size_t N, unsigned K>
inline Taylor<N, K> operator+(const Taylor<N, K>& a) { return a; }
template <std::size_t N, unsigned K>
inline Taylor<N, K> operator-(const Taylor<N, K>& a) {
  Taylor<N, K> r;
  r.v = -a.v;
  for (std::size_t i = 0; i < N; i++) r.d[i] = -a.d[i];
  if constexpr (K == 2)
    for (std::size_t i = 0; i < N; i++)
      for (std::size_t j = 0; j < N; j++) r.h[i][j] = -a.h[i][j];
  return r;
}
template <std::size_t N, unsigned K>
inline Taylor<N, K>& operator+=(Taylor<N, K>& a, const Taylor<N, K>& b) {
  a.v += b.v;
  for (std::size_t i = 0; i < N; i++) a.d[i] += b.d[i];
  if constexpr (K == 2)
    for (std::size_t i = 0; i < N; i++)
      for (std::size_t j = 0; j < N; j++) a.h[i][j] += b.h[i][j];
  return a;
}
template <std::size_t N, unsigned K>
inline Taylor<N, K>& operator-=(Taylor<N, K>& a, const Taylor<N, K>& b) {
  a.v -= b.v;
  for (std::size_t i = 0; i < N; i++) a.d[i] -= b.d[i];
  if constexpr (K == 2)
    for (std::size_t i = 0; i < N; i++)
      for (std::size_t j = 0; j < N; j++) a.h[i][j] -= b.h[i][j];
  return a;
}
template <std::size_t N, unsigned K>
inline Taylor<N, K>& operator*=(Taylor<N, K>& a, double c) {
  a.v *= c;
  for (std::size_t i = 0; i < N; i++) a.d[i] *= c;
  if constexpr (K == 2)
    for (std::size_t i = 0; i < N; i++)
      for (std::size_t j = 0; j < N; j++) a.h[i][j] *= c;
  return a;
}
template <std::size_t N, unsigned K>
inline Taylor<N, K> operator*(const Taylor<N, K>& a, const Taylor<N, K>& b) {
  Taylor<N, K> r;
  r.v = a.v * b.v;
  for (std::size_t i = 0; i < N; i++) r.d[i] = a.v * b.d[i] + b.v * a.d[i];
  if constexpr (K == 2)
    for (std::size_t i = 0; i < N; i++)
      for (std::size_t j = 0; j < N; j++) r.h[i][j] = a.v * b.h[i][j] + b.v * a.h[i][j] + a.d[i] * b.d[j] + a.d[j] * b.d[i];
  return r;
}
template <std::size_t N, unsigned K>
inline Taylor<N, K>& operator*=(Taylor<N, K>& a, const Taylor<N, K>& b) { return a = a * b; }
template <std::size_t N, unsigned K>
inline Taylor<N, K> reciprocal(const Taylor<N, K>& a) {
  const double r = 1 / a.v;
  return chain(a, r, -r * r, 2 * r * r * r);
}
template <std::size_t N, unsigned K>
inline Taylor<N, K> operator/(const Taylor<N, K>& a, const Taylor<N, K>& b) { return a * reciprocal(b); }
template <std::size_t N, unsigned K>
inline Taylor<N, K>& operator/=(Taylor<N, K>& a, const Taylor<N, K>& b) { return a = a / b; }

template <std::size_t N, unsigned K>
inline Taylor<N, K> operator+(Taylor<N, K> a, const Taylor<N, K>& b) { return a += b; }
template <std::size_t N, unsigned K>
inline Taylor<N, K> operator-(Taylor<N, K> a, const Taylor<N, K>& b) { return a -= b; }
template <std::size_t N, unsigned K>
inline Taylor<N, K> operator+(Taylor<N, K> a, double c) { a.v += c; return a; }
template <std::size_t N, unsigned K>
inline Taylor<N, K> operator+(double c, Taylor<N, K> a) { a.v += c; return a; }
template <std::size_t N, unsigned K>
inline Taylor<N, K> operator-(Taylor<N, K> a, double c) { a.v -= c; return a; }
template <std::size_t N, unsigned K>
inline Taylor<N, K> operator-(double c, const Taylor<N, K>& a) { return -a + c; }
template <std::size_t N, unsigned K>
inline Taylor<N, K> operator*(Taylor<N, K> a, double c) { return a *= c; }
template <std::size_t N, unsigned K>
inline Taylor<N, K> operator*(double c, Taylor<N, K> a) { return a *= c; }
template <std::size_t N, unsigned K>
inline Taylor<N, K> operator/(Taylor<N, K> a, double c) { return a *= 1 / c; }
template <std::size_t N, unsigned K>
inline Taylor<N, K> operator/(double c, const Taylor<N, K>& a) { return reciprocal(a) * c; }
template <std::size_t N, unsigned K>
inline Taylor<N, K>& operator+=(Taylor<N, K>& a, double c) { a.v += c; return a; }
template <std::size_t N, unsigned K>
inline Taylor<N, K>& operator-=(Taylor<N, K>& a, double c) { a.v -= c; return a; }
template <std::size_t N, unsigned K>
inline Taylor<N, K>& operator/=(Taylor<N, K>& a, double c) { return a *= 1 / c; }

// Comparisons are of values.
template <std::size_t N, unsigned K> inline bool operator<(const Taylor<N, K>& a, const Taylor<N, K>& b) { return a.v < b.v; }
template <std::size_t N, unsigned K> inline bool operator>(const Taylor<N, K>& a, const Taylor<N, K>& b) { return a.v > b.v; }
template <std::size_t N, unsigned K> inline bool operator<=(const Taylor<N, K>& a, const Taylor<N, K>& b) { return a.v <= b.v; }
template <std::size_t N, unsigned K> inline bool operator>=(const Taylor<N, K>& a, const Taylor<N, K>& b) { return a.v >= b.v; }
template <std::size_t N, unsigned K> inline bool operator<(const Taylor<N, K>& a, double c) { return a.v < c; }
template <std::size_t N, unsigned K> inline bool operator>(const Taylor<N, K>& a, double c) { return a.v > c; }
template <std::size_t N, unsigned K> inline bool operator<=(const Taylor<N, K>& a, double c) { return a.v <= c; }
template <std::size_t N, unsigned K> inline bool operator>=(const Taylor<N, K>& a, double c) { return a.v >= c; }
template <std::size_t N, unsigned K> inline bool operator<(double c, const Taylor<N, K>& a) { return c < a.v; }
template <std::size_t N, unsigned K> inline bool operator>(double c, const Taylor<N, K>& a) { return c > a.v; }

// ---- elementary functions ---------------------------------------------------------------------------------------------

template <std::size_t N, unsigned K>
inline Taylor<N, K> sqrt(const Taylor<N, K>& a) {
  const double s = std::sqrt(a.v);
  return chain(a, s, 0.5 / s, -0.25 / (s * a.v));
}
template <std::size_t N, unsigned K>
inline Taylor<N, K> exp(const Taylor<N, K>& a) {
  const double e = std::exp(a.v);
  return chain(a, e, e, e);
}
template <std::size_t N, unsigned K>
inline Taylor<N, K> log(const Taylor<N, K>& a) { return chain(a, std::log(a.v), 1 / a.v, -1 / (a.v * a.v)); }
template <std::size_t N, unsigned K>
inline Taylor<N, K> sin(const Taylor<N, K>& a) {
  const double s = std::sin(a.v), c = std::cos(a.v);
  return chain(a, s, c, -s);
}
template <std::size_t N, unsigned K>
inline Taylor<N, K> cos(const Taylor<N, K>& a) {
  const double s = std::sin(a.v), c = std::cos(a.v);
  return chain(a, c, -s, -c);
}
template <std::size_t N, unsigned K>
inline Taylor<N, K> tan(const Taylor<N, K>& a) {
  const double t = std::tan(a.v), s = 1 + t * t;
  return chain(a, t, s, 2 * t * s);
}
template <std::size_t N, unsigned K>
inline Taylor<N, K> sinh(const Taylor<N, K>& a) {
  const double s = std::sinh(a.v), c = std::cosh(a.v);
  return chain(a, s, c, s);
}
template <std::size_t N, unsigned K>
inline Taylor<N, K> cosh(const Taylor<N, K>& a) {
  const double s = std::sinh(a.v), c = std::cosh(a.v);
  return chain(a, c, s, c);
}
template <std::size_t N, unsigned K>
inline Taylor<N, K> tanh(const Taylor<N, K>& a) {
  const double t = std::tanh(a.v), s = 1 - t * t;
  return chain(a, t, s, -2 * t * s);
}
template <std::size_t N, unsigned K>
inline Taylor<N, K> atan(const Taylor<N, K>& a) {
  const double q = 1 / (1 + a.v * a.v);
  return chain(a, std::atan(a.v), q, -2 * a.v * q * q);
}
// atan2(y, x), with ∂_y = x / r², ∂_x = -y / r².
template <std::size_t N, unsigned K>
inline Taylor<N, K> atan2(const Taylor<N, K>& y, const Taylor<N, K>& x) {
  const double r2 = x.v * x.v + y.v * y.v, r4 = r2 * r2;
  return chain(y, x, std::atan2(y.v, x.v), x.v / r2, -y.v / r2, -2 * x.v * y.v / r4, (y.v * y.v - x.v * x.v) / r4, 2 * x.v * y.v / r4);
}
template <std::size_t N, unsigned K>
inline Taylor<N, K> abs(const Taylor<N, K>& a) { return a.v < 0 ? -a : a; }
// a^p for a real exponent (a > 0 unless p is an integer).
template <std::size_t N, unsigned K>
inline Taylor<N, K> pow(const Taylor<N, K>& a, double p) {
  if (p == 0) return Taylor<N, K>(1.0);
  const double f0 = std::pow(a.v, p), f1 = p * std::pow(a.v, p - 1), f2 = p * (p - 1) * std::pow(a.v, p - 2);
  return chain(a, f0, f1, f2);
}
// a^n for an integer exponent, by the same rule, defined for every a (a != 0 when n < 0).
template <std::size_t N, unsigned K>
inline Taylor<N, K> pow(const Taylor<N, K>& a, int n) {
  if (n == 0) return Taylor<N, K>(1.0);
  const double f1 = n * std::pow(a.v, n - 1), f2 = n == 1 ? 0.0 : n * (n - 1) * std::pow(a.v, n - 2);
  return chain(a, std::pow(a.v, n), f1, f2);
}
// a^b = exp(b log a), a > 0.
template <std::size_t N, unsigned K>
inline Taylor<N, K> pow(const Taylor<N, K>& a, const Taylor<N, K>& b) { return exp(b * log(a)); }

// ---- derivatives of functions -----------------------------------------------------------------------------------------

// The jet of f at x: f evaluated on the jets of the coordinates.
template <std::size_t N, class F>
inline Jet<N> taylor(F&& f, const std::array<double, N>& x) { return f(variables<N, 2>(x)); }

template <std::size_t N, class F>
inline std::array<double, N> gradient(F&& f, const std::array<double, N>& x) { return f(variables<N, 1>(x)).d; }

template <std::size_t N, class F>
inline std::array<std::array<double, N>, N> hessian(F&& f, const std::array<double, N>& x) { return f(variables<N, 2>(x)).h; }

}  // namespace calx::differentiation
