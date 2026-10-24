// Checks of calx/optimisation: Brent's method on minima known in closed form, to its guarantee;
// superlinear convergence (a parabola in a handful of evaluations, far fewer than golden
// section's); a minimum at an end of the interval; a quartic, flat to third order, where the
// tolerance in x is all f can tell; a function with a kink, carried by golden steps; the maximum
// of a log-likelihood; and the refusal of an empty interval.

#include <calx/optimisation/minimise.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <stdexcept>
#include <string>

using namespace calx::optimisation;

int main() {
  calx::verification::Suite s;
  const double pi = 3.141592653589793;

  const auto parabola = brent_minimise([](double x) { return (x - 2) * (x - 2) + 1; }, -3, 7, 1e-10);
  s.check("optimisation: the vertex of (x − 2)² + 1, to the tolerance", parabola.x - 2, 1e-9);
  s.require("optimisation: a parabola is found by interpolation in a few evaluations", parabola.converged && parabola.evaluations <= 12,
            std::to_string(parabola.evaluations) + " evaluations");
  const auto golden = golden_section([](double x) { return (x - 2) * (x - 2) + 1; }, -3, 7, 1e-9);
  s.require("optimisation: golden section needs many more on the same parabola", golden.evaluations > 3 * parabola.evaluations,
            std::to_string(golden.evaluations) + " against " + std::to_string(parabola.evaluations));

  const auto cosine = brent_minimise([](double x) { return std::cos(x); }, 2, 4.5, 1e-12);
  s.check("optimisation: the minimum of cos on [2, 4.5] is π", cosine.x - pi, 2e-8);
  s.check("optimisation: and its value −1", cosine.f + 1, 1e-15);

  // The minimum of x e^x is at x = −1, value −1/e.
  const auto xex = brent_minimise([](double x) { return x * std::exp(x); }, -4, 1, 1e-12);
  s.check("optimisation: x eˣ is least at −1", xex.x + 1, 2e-8);
  s.check("optimisation: where it is −1/e", xex.f + std::exp(-1.0), 1e-15);

  const auto end = brent_minimise([](double x) { return x; }, 0.5, 3, 1e-9);
  // At an end the guarantee is the final half-width doubled: 2 (√eps |x| + tolerance / 3).
  s.check("optimisation: a monotone function is least at the interval's left end, to the guarantee", end.x - 0.5, 2 * (std::sqrt(2.220446049250313e-16) * 0.5 + 1e-9 / 3));

  // (x − 1)⁴ is flat: f changes by eps only beyond eps^(1/4) of the minimum, so x is known only to about 1e-4.
  const auto quartic = brent_minimise([](double x) { return std::pow(x - 1, 4); }, -2, 5, 1e-12);
  s.check("optimisation: a quartic minimum, to what its flatness allows", quartic.x - 1, 2e-3);

  // A kink: no parabola fits; golden steps carry it.
  const auto kink = brent_minimise([](double x) { return std::fabs(x - 0.3); }, -1, 1, 1e-10);
  s.check("optimisation: |x − 0.3|, a kink, to the tolerance", kink.x - 0.3, 1e-9);

  // The maximum of a binomial log-likelihood, 7 successes in 20: p̂ = 0.35.
  const auto mle = brent_maximise([](double p) { return 7 * std::log(p) + 13 * std::log(1 - p); }, 1e-6, 1 - 1e-6, 1e-12);
  s.check("optimisation: the maximum of a binomial log-likelihood is the proportion seen", mle.x - 0.35, 1e-8);

  bool threw = false;
  try { brent_minimise([](double x) { return x; }, 1, 1); } catch (const std::invalid_argument&) { threw = true; }
  s.require("optimisation: an empty interval is refused", threw);
  return s.done();
}
