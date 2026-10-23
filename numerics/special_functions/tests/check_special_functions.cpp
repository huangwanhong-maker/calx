// Checks of calx/special_functions: erfcx against values computed from the C library's
// erfc where both are representable, the meeting of its two branches at x = 10, the
// reflection for negative arguments, and the leading terms of its asymptotic series.

#include <calx/special_functions/erf.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <numbers>

using calx::special_functions::erfcx;
using calx::verification::relative;

int main() {
  calx::verification::Suite s;
  double worst = 0;
  for (double x = 0; x <= 5; x += 0.125) worst = std::max(worst, relative(erfcx(x), std::exp(x * x) * std::erfc(x)));
  s.check("special_functions: erfcx(x) = exp(x^2) erfc(x) on [0, 5]", worst, 1e-15);
  s.check("special_functions: erfcx(1) = e erfc(1) = 0.42758357615580700", relative(erfcx(1), 0.42758357615580700), 2e-16);
  // At x = 10 the series and the product meet; the product carries the rounding of x^2.
  const double below = std::exp(100.0) * std::erfc(10.0);
  s.check("special_functions: erfcx's branches agree at x = 10", relative(erfcx(10), below), 3e-14);
  // erfcx(x) ~ 1 / (x sqrt(pi)) (1 - 1/(2x^2) + 3/(4x^4) - 15/(8x^6)), error O(x^-9).
  const double x = 1000, lead = 1 / (x * std::sqrt(std::numbers::pi)) * (1 - 1 / (2 * x * x) + 3 / (4 * std::pow(x, 4)) - 15 / (8 * std::pow(x, 6)));
  // The truncation is below 1e-23; what remains is the rounding of the two sums, a few ulps.
  s.check("special_functions: erfcx(1000) by its asymptotic series, to a few ulps", relative(erfcx(x), lead), 4e-16);
  double refl = 0;
  for (double y = 0.25; y <= 5; y += 0.25) refl = std::max(refl, relative(erfcx(-y), std::exp(y * y) * (1 + std::erf(y))));
  s.check("special_functions: erfcx(-x) = exp(x^2) (1 + erf(x)) on [0.25, 5]", refl, 1e-14);
  s.require("special_functions: erfcx(0) = 1 exactly", erfcx(0) == 1);
  return s.done();
}
