// The gamma function and its relatives: log Γ, the digamma and trigamma functions, the
// regularised incomplete gamma functions P(a, x) and Q(a, x) with the inverse of P, and the
// regularised incomplete beta function I_x(a, b). With them come the tails and quantiles of
// the χ² distribution, the gamma distribution's quantiles (the discrete-Γ rates of
// phylogenetics), and the binomial, F and Student distributions' tails.
//
// - log_gamma(x), x > 0: Lanczos' approximation with g = 7 and nine coefficients (Lanczos
//   1964, SIAM J. Numer. Anal. B1, 86–96; the coefficients of Godfrey's tabulation as Press
//   et al., Numerical Recipes, 3rd ed. 2007, §6.1, give a set of the kind), relative error
//   below 2e-15 away from the zeros at 1 and 2, where the error is absolute, about 1e-15;
//   reflection Γ(x)Γ(1 − x) = π / sin(πx) for 0 < x < 0.5. It is ours rather than the C
//   library's lgamma, which writes the global signgam and so is not safe on several threads.
// - digamma ψ and trigamma ψ′, x > 0: the recurrences ψ(x) = ψ(x + 1) − 1/x and
//   ψ′(x) = ψ′(x + 1) + 1/x² up to x ≥ 10, then the asymptotic series through x^{-14} and
//   x^{-15} (Abramowitz and Stegun 1964, 6.3.18 and 6.4.12); absolute error about 1e-15.
// - gamma_p, gamma_q: the series for x < a + 1 and Legendre's continued fraction otherwise,
//   evaluated by the modified Lentz method (Numerical Recipes §6.2 and §5.2); relative error
//   near 1e-14 for the moderate a (below 10^4) these functions are used for. Q is computed
//   directly where it is small, never as 1 − P.
// - gamma_p_inverse(a, p): Halley's method from the initial guesses of Numerical Recipes
//   §6.2.1 (after DiDonato and Morris 1986, ACM TOMS 12, 377–393); converged to a relative
//   1e-13 in x.
// - beta_i(a, b, x): the continued fraction of Numerical Recipes §6.4 (Lentz), with the
//   symmetry I_x(a, b) = 1 − I_{1−x}(b, a) where it converges faster.

#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>

namespace calx::special_functions {

inline double log_gamma(double x) {
  if (!(x > 0)) {
    if (std::isnan(x)) return x;
    throw std::domain_error("calx::special_functions::log_gamma: x > 0");
  }
  if (x < 0.5) return std::log(std::numbers::pi / std::sin(std::numbers::pi * x)) - log_gamma(1 - x);
  static constexpr double c[9] = {0.99999999999980993,  676.5203681218851,     -1259.1392167224028,
                                  771.32342877765313,   -176.61502916214059,   12.507343278686905,
                                  -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7};
  const double z = x - 1;
  double a = c[0];
  for (int i = 1; i < 9; i++) a += c[i] / (z + i);
  const double t = z + 7.5;
  return 0.5 * std::log(2 * std::numbers::pi) + (z + 0.5) * std::log(t) - t + std::log(a);
}

// log B(a, b) = log Γ(a) + log Γ(b) − log Γ(a + b).
inline double log_beta(double a, double b) { return log_gamma(a) + log_gamma(b) - log_gamma(a + b); }

// log C(n, k), for real n ≥ k ≥ 0.
inline double log_choose(double n, double k) { return log_gamma(n + 1) - log_gamma(k + 1) - log_gamma(n - k + 1); }

inline double digamma(double x) {
  if (!(x > 0)) throw std::domain_error("calx::special_functions::digamma: x > 0");
  double r = 0;
  while (x < 10) { r -= 1 / x; x += 1; }
  const double f = 1 / (x * x);
  return r + std::log(x) - 0.5 / x -
         f * (1.0 / 12 - f * (1.0 / 120 - f * (1.0 / 252 - f * (1.0 / 240 - f * (1.0 / 132 - f * (691.0 / 32760 - f / 12))))));
}

inline double trigamma(double x) {
  if (!(x > 0)) throw std::domain_error("calx::special_functions::trigamma: x > 0");
  double r = 0;
  while (x < 10) { r += 1 / (x * x); x += 1; }
  const double f = 1 / (x * x);
  return r + 1 / x + f / 2 + f / x * (1.0 / 6 - f * (1.0 / 30 - f * (1.0 / 42 - f * (1.0 / 30 - f * (5.0 / 66 - f * (691.0 / 2730 - f * 7.0 / 6))))));
}

namespace detail {

inline constexpr double EPS = 1e-16;
inline constexpr double TINY = 1e-300;

// P(a, x) by its series, x < a + 1.
inline double gamma_series(double a, double x) {
  double ap = a, del = 1 / a, sum = del;
  for (int n = 0; n < 100000; n++) {
    ap += 1;
    del *= x / ap;
    sum += del;
    if (std::fabs(del) < std::fabs(sum) * EPS) break;
  }
  return sum * std::exp(-x + a * std::log(x) - log_gamma(a));
}

// Q(a, x) by Legendre's continued fraction, x ≥ a + 1.
inline double gamma_fraction(double a, double x) {
  double b = x + 1 - a, c = 1 / TINY, d = 1 / b, h = d;
  for (int i = 1; i < 100000; i++) {
    const double an = -i * (i - a);
    b += 2;
    d = an * d + b;
    if (std::fabs(d) < TINY) d = TINY;
    c = b + an / c;
    if (std::fabs(c) < TINY) c = TINY;
    d = 1 / d;
    const double del = d * c;
    h *= del;
    if (std::fabs(del - 1) < EPS) break;
  }
  return std::exp(-x + a * std::log(x) - log_gamma(a)) * h;
}

}  // namespace detail

// The regularised lower incomplete gamma function P(a, x) = γ(a, x) / Γ(a), a > 0, x ≥ 0.
inline double gamma_p(double a, double x) {
  if (!(a > 0) || x < 0) throw std::domain_error("calx::special_functions::gamma_p: a > 0 and x >= 0");
  if (x == 0) return 0;
  if (std::isinf(x)) return 1;
  return x < a + 1 ? detail::gamma_series(a, x) : 1 - detail::gamma_fraction(a, x);
}

// The regularised upper incomplete gamma function Q(a, x) = 1 − P(a, x).
inline double gamma_q(double a, double x) {
  if (!(a > 0) || x < 0) throw std::domain_error("calx::special_functions::gamma_q: a > 0 and x >= 0");
  if (x == 0) return 1;
  if (std::isinf(x)) return 0;
  return x < a + 1 ? 1 - detail::gamma_series(a, x) : detail::gamma_fraction(a, x);
}

// The x with P(a, x) = p, for 0 ≤ p ≤ 1.
inline double gamma_p_inverse(double a, double p) {
  if (!(a > 0) || p < 0 || p > 1) throw std::domain_error("calx::special_functions::gamma_p_inverse: a > 0 and 0 <= p <= 1");
  if (p == 0) return 0;
  if (p == 1) return std::numeric_limits<double>::infinity();
  const double gln = log_gamma(a), a1 = a - 1;
  double x, lna1 = 0, afac = 0;
  if (a > 1) {
    lna1 = std::log(a1);
    afac = std::exp(a1 * (lna1 - 1) - gln);
    const double pp = p < 0.5 ? p : 1 - p;
    const double t = std::sqrt(-2 * std::log(pp));
    x = (2.30753 + t * 0.27061) / (1 + t * (0.99229 + t * 0.04481)) - t;
    if (p < 0.5) x = -x;
    x = std::max(1e-3, a * std::pow(1 - 1 / (9 * a) - x / (3 * std::sqrt(a)), 3));
  } else {
    const double t = 1 - a * (0.253 + a * 0.12);
    x = p < t ? std::pow(p / t, 1 / a) : 1 - std::log(1 - (p - t) / (1 - t));
  }
  for (int j = 0; j < 100; j++) {
    if (x <= 0) return 0;
    const double err = gamma_p(a, x) - p;
    const double t = a > 1 ? afac * std::exp(-(x - a1) + a1 * (std::log(x) - lna1)) : std::exp(-x + a1 * std::log(x) - gln);
    const double u = err / t;
    const double step = u / (1 - 0.5 * std::min(1.0, u * ((a - 1) / x - 1)));
    x -= step;
    if (x <= 0) x = 0.5 * (x + step);
    if (std::fabs(step) < 1e-13 * x) break;
  }
  return x;
}

// The tail of the χ² distribution with k degrees of freedom above x, and its quantile.
inline double chi_square_q(double x, double k) { return gamma_q(k / 2, x / 2); }
inline double chi_square_quantile(double p, double k) { return 2 * gamma_p_inverse(k / 2, p); }

namespace detail {

// The continued fraction of the incomplete beta function (Numerical Recipes 6.4.5).
inline double beta_fraction(double a, double b, double x) {
  const double qab = a + b, qap = a + 1, qam = a - 1;
  double c = 1, d = 1 - qab * x / qap;
  if (std::fabs(d) < TINY) d = TINY;
  d = 1 / d;
  double h = d;
  for (int m = 1; m < 100000; m++) {
    const int m2 = 2 * m;
    double aa = m * (b - m) * x / ((qam + m2) * (a + m2));
    d = 1 + aa * d;
    if (std::fabs(d) < TINY) d = TINY;
    c = 1 + aa / c;
    if (std::fabs(c) < TINY) c = TINY;
    d = 1 / d;
    h *= d * c;
    aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
    d = 1 + aa * d;
    if (std::fabs(d) < TINY) d = TINY;
    c = 1 + aa / c;
    if (std::fabs(c) < TINY) c = TINY;
    d = 1 / d;
    const double del = d * c;
    h *= del;
    if (std::fabs(del - 1) < EPS) break;
  }
  return h;
}

}  // namespace detail

// The regularised incomplete beta function I_x(a, b), a, b > 0, 0 ≤ x ≤ 1.
inline double beta_i(double a, double b, double x) {
  if (!(a > 0) || !(b > 0) || x < 0 || x > 1) throw std::domain_error("calx::special_functions::beta_i: a, b > 0 and 0 <= x <= 1");
  if (x == 0 || x == 1) return x;
  const double bt = std::exp(-log_beta(a, b) + a * std::log(x) + b * std::log1p(-x));
  return x < (a + 1) / (a + b + 2) ? bt * detail::beta_fraction(a, b, x) / a : 1 - bt * detail::beta_fraction(b, a, 1 - x) / b;
}

// P(X ≤ k) for X binomial with n trials of probability p: I_{1−p}(n − k, k + 1).
inline double binomial_cdf(double k, double n, double p) {
  if (k < 0) return 0;
  if (k >= n) return 1;
  return beta_i(n - k, k + 1, 1 - p);
}

// The two-sided tail of Student's t with ν degrees of freedom beyond |t|: I_{ν/(ν+t²)}(ν/2, 1/2).
inline double student_t_two_sided(double t, double nu) { return beta_i(nu / 2, 0.5, nu / (nu + t * t)); }

// The upper tail of the F distribution with (d1, d2) degrees of freedom above f.
inline double f_q(double f, double d1, double d2) { return beta_i(d2 / 2, d1 / 2, d2 / (d2 + d1 * f)); }

}  // namespace calx::special_functions
