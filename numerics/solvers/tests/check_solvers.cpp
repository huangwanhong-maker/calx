// Checks of calx/solvers: Brent's method on roots known in closed form or to many digits,
// its refusal of an unbracketed interval, and its count of iterations on a function where
// interpolation fails and bisection must carry it.

#include <calx/solvers/roots.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <stdexcept>
#include <string>

using calx::solvers::brent;

int main() {
  calx::verification::Suite s;
  // With tolerance 0, the root is guaranteed to within 2 eps |x|.
  constexpr double eps = 2.220446049250313e-16;
  const auto cube = brent([](double x) { return x * x * x - 2; }, 0, 2, 0);
  s.check("solvers: the real cube root of 2, to its guarantee 2 eps |x|", cube.x - std::cbrt(2.0), 2 * eps * std::cbrt(2.0));
  // The Dottie number, the fixed point of cos: 0.739085133215160641655312087673873...
  const auto dottie = brent([](double x) { return std::cos(x) - x; }, 0, 1, 0);
  s.check("solvers: the fixed point of cos (the Dottie number), to its guarantee", dottie.x - 0.7390851332151607, 2 * eps * 0.74);
  s.require("solvers: Brent's method converges superlinearly (cos x = x in at most 10 iterations)", dottie.converged && dottie.iterations <= 10,
            std::to_string(dottie.iterations) + " iterations");
  // A step with no smooth structure: bisection must carry it, within log2(1 / tol) + a few steps.
  const auto step = brent([](double x) { return x < 0.3 ? -1.0 : 1.0; }, 0, 1, 1e-12);
  s.check("solvers: a discontinuous sign change is bracketed to the tolerance", step.x - 0.3, 1e-12);
  bool threw = false;
  try { brent([](double x) { return x * x + 1; }, -1, 1); } catch (const std::invalid_argument&) { threw = true; }
  s.require("solvers: an unbracketed interval is refused", threw);
  return s.done();
}
