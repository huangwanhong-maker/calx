// Checks of the gamma family: log Γ at the integers (log factorials) and at 1/2, against the C
// library over a range; ψ(1) = −γ, ψ(1/2) = −γ − 2 log 2, ψ′(1) = π²/6, ψ′(1/2) = π²/2, and the
// recurrences; P(1, x) = 1 − e^{−x} and P(1/2, x) = erf √x; P + Q = 1; the χ² quantiles of the
// tables; P inverted; I_x(a, 1) = x^a, I_x(1, b) = 1 − (1 − x)^b, the symmetry, I_{1/2}(a, a) =
// 1/2; and the binomial distribution against its sum.

#include <calx/special_functions/gamma.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>
#include <numbers>

using namespace calx::special_functions;

int main() {
  calx::verification::Suite s;
  constexpr double euler = 0.57721566490153286061;

  {
    double err = 0, factorial = 1;
    for (int n = 1; n <= 20; n++) {
      if (n > 1) factorial *= (n - 1);
      err = std::max(err, std::fabs(log_gamma(n) - std::log(factorial)));
    }
    s.check("log Γ(n) = log (n − 1)! for n = 1 to 20", err, 2e-14);
    s.check("log Γ(1/2) = log √π", log_gamma(0.5) - 0.5 * std::log(std::numbers::pi), 1e-15);
    double rel = 0;
    for (double x = 0.01; x < 300; x *= 1.07) rel = std::max(rel, std::fabs(log_gamma(x) - std::lgamma(x)) / std::max(1.0, std::fabs(std::lgamma(x))));
    s.check("log Γ against the C library on [0.01, 300]", rel, 1e-13);
  }
  {
    s.check("ψ(1) = −γ", digamma(1) + euler, 1e-14);
    s.check("ψ(1/2) = −γ − 2 log 2", digamma(0.5) + euler + 2 * std::log(2.0), 1e-14);
    s.check("ψ′(1) = π²/6", trigamma(1) - std::numbers::pi * std::numbers::pi / 6, 1e-14);
    s.check("ψ′(1/2) = π²/2", trigamma(0.5) - std::numbers::pi * std::numbers::pi / 2, 1e-13);
    double rec = 0;
    for (double x = 0.3; x < 50; x *= 1.3) rec = std::max({rec, std::fabs(digamma(x + 1) - digamma(x) - 1 / x), std::fabs(trigamma(x) - trigamma(x + 1) - 1 / (x * x))});
    s.check("the recurrences ψ(x + 1) = ψ(x) + 1/x and ψ′(x) = ψ′(x + 1) + 1/x²", rec, 1e-13);
  }
  {
    double e1 = 0, eh = 0, sum = 0;
    for (double x = 0.001; x < 40; x *= 1.2) {
      e1 = std::max(e1, std::fabs(gamma_p(1, x) - (-std::expm1(-x))));
      eh = std::max(eh, std::fabs(gamma_p(0.5, x) - std::erf(std::sqrt(x))));
      for (double a : {0.1, 0.5, 2.5, 10.0, 150.0}) sum = std::max(sum, std::fabs(gamma_p(a, x) + gamma_q(a, x) - 1));
    }
    s.check("P(1, x) = 1 − e^{−x}", e1, 1e-14);
    s.check("P(1/2, x) = erf √x", eh, 1e-14);
    s.check("P + Q = 1", sum, 1e-14);
    s.check("Q(1, 30) = e^{−30}, relative (computed directly, not as 1 − P)", gamma_q(1, 30) / std::exp(-30.0) - 1, 1e-13);
  }
  {
    // Upper 5% and 1% points of χ² (any table, e.g. Abramowitz and Stegun 1964, Table 26.8).
    s.check("χ² with 1 degree of freedom: 95% point 3.8414588", chi_square_quantile(0.95, 1) - 3.841458820694124, 1e-10);
    s.check("χ² with 10 degrees of freedom: 95% point 18.307038", chi_square_quantile(0.95, 10) - 18.307038053275146, 1e-9);
    s.check("χ² with 5 degrees of freedom: 99% point 15.086272", chi_square_quantile(0.99, 5) - 15.086272469388987, 1e-9);
    s.check("its tail at that point is 0.01", chi_square_q(15.086272469388987, 5) - 0.01, 1e-13);
    double inv = 0;
    for (double a : {0.05, 0.3, 1.0, 2.0, 7.5, 60.0, 900.0})
      for (double p : {1e-8, 0.001, 0.1, 0.5, 0.9, 0.999})
        inv = std::max(inv, std::fabs(gamma_p(a, gamma_p_inverse(a, p)) - p) / p);
    s.check("P(a, P⁻¹(a, p)) = p, relative, for a from 0.05 to 900", inv, 1e-10);
  }
  {
    double e = 0;
    for (double x = 0.01; x < 1; x += 0.07)
      for (double a : {0.3, 1.0, 4.0, 25.0}) {
        e = std::max(e, std::fabs(beta_i(a, 1, x) - std::pow(x, a)));
        e = std::max(e, std::fabs(beta_i(1, a, x) - (1 - std::pow(1 - x, a))));
        e = std::max(e, std::fabs(beta_i(a, 2.5, x) + beta_i(2.5, a, 1 - x) - 1));
      }
    s.check("I_x(a, 1) = x^a, I_x(1, b) = 1 − (1 − x)^b, and I_x(a, b) = 1 − I_{1−x}(b, a)", e, 1e-13);
    double h = 0;
    for (double a : {0.5, 3.0, 40.0, 400.0}) h = std::max(h, std::fabs(beta_i(a, a, 0.5) - 0.5));
    s.check("I_{1/2}(a, a) = 1/2", h, 1e-12);
    // Binomial(30, 0.3): P(X ≤ k) as a sum of the terms.
    double b = 0, cdf = 0;
    for (int k = 0; k <= 30; k++) {
      cdf += std::exp(log_choose(30, k) + k * std::log(0.3) + (30 - k) * std::log(0.7));
      b = std::max(b, std::fabs(binomial_cdf(k, 30, 0.3) - cdf));
    }
    s.check("the binomial distribution (30, 0.3) against the sum of its terms", b, 1e-13);
    // Student's t with an even number ν of degrees of freedom has a closed form (Abramowitz and
    // Stegun 1964, 26.7.3): P(|T| < t) = sin θ (1 + cos²θ/2 + 1·3 cos⁴θ/(2·4) + ... ), θ = atan(t/√ν).
    double st = 0;
    for (double nu : {2.0, 10.0, 30.0})
      for (double t = 0.1; t < 8; t *= 1.4) {
        const double th = std::atan(t / std::sqrt(nu)), c = std::cos(th);
        double term = 1, total = 1;
        for (int k = 1; k < nu / 2; k++) { term *= (2.0 * k - 1) / (2.0 * k) * c * c; total += term; }
        st = std::max(st, std::fabs(student_t_two_sided(t, nu) - (1 - std::sin(th) * total)));
      }
    s.check("Student's t with 2, 10 and 30 degrees of freedom against the closed form for even ν", st, 1e-13);
  }

  return s.done();
}
