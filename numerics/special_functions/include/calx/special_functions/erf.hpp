// The scaled complementary error function, erfcx(x) = exp(x^2) erfc(x).
//
// It is the bounded form in which the error function enters, for instance, the firing rate
// of an integrate-and-fire neuron driven by noise (Siegert 1951; Brunel 2000), where
// exp(u^2) (1 + erf(u)) = erfcx(-u) must be integrated over arguments at which exp(u^2)
// alone overflows.
//
// - 0 <= x < 10: exp(x^2) erfc(x) with the C library's erfc, whose relative accuracy holds
//   until erfc underflows near x = 26; the rounding of x^2 costs a relative error of about
//   x^2 eps, below 2e-14 here.
// - x >= 10: the asymptotic series erfcx(x) ~ (1 / (x sqrt(pi))) sum_n (-1)^n (2n-1)!! /
//   (2x^2)^n (Abramowitz and Stegun 1964, 7.1.23), summed while its terms decrease and
//   exceed 1e-17 of the sum; at x = 10 it reaches that in fourteen terms.
// - x < 0: erfcx(x) = 2 exp(x^2) - erfcx(-x); it overflows below about -26.6.

#pragma once

#include <cmath>
#include <numbers>

namespace calx::special_functions {

inline double erfcx(double x) {
  if (std::isnan(x)) return x;
  if (x < 0) return 2 * std::exp(x * x) - erfcx(-x);
  if (x < 10) return std::exp(x * x) * std::erfc(x);
  const double y = 1 / (2 * x * x);
  double term = 1, sum = 1;
  for (int n = 1; n < 60; n++) {
    const double next = -term * (2 * n - 1) * y;
    if (std::fabs(next) >= std::fabs(term) || std::fabs(next) < 1e-17 * std::fabs(sum)) {
      if (std::fabs(next) < std::fabs(term)) sum += next;
      break;
    }
    term = next;
    sum += term;
  }
  return sum / (x * std::sqrt(std::numbers::pi));
}

}  // namespace calx::special_functions
