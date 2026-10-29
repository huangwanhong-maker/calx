// Checks of calx/optimisation's Nelder–Mead method: Rosenbrock's function from the classical
// start (−1.2, 1) to its minimum (1, 1), with the standard and the adaptive coefficients; the
// exact minimum of a positive definite quadratic in three variables; a function of one
// variable against Brent's method; Powell's quartic-like function of four variables, whose
// minimum is 0 at the origin and whose Hessian is singular there; the budget of evaluations,
// honoured; a NaN treated as +∞; an initial simplex of the caller's steps.

#include <calx/optimisation/minimise.hpp>
#include <calx/optimisation/nelder_mead.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <span>
#include <string>
#include <vector>

using namespace calx::optimisation;

int main() {
  calx::verification::Suite s;
  auto rosenbrock = [](std::span<const double> x) { return 100 * (x[1] - x[0] * x[0]) * (x[1] - x[0] * x[0]) + (1 - x[0]) * (1 - x[0]); };

  NelderMeadOptions tight;
  tight.x_tolerance = 1e-10;
  tight.f_tolerance = 1e-14;
  for (bool adaptive : {false, true}) {
    NelderMeadOptions o = tight;
    o.adaptive = adaptive;
    const auto r = nelder_mead(rosenbrock, {-1.2, 1}, o);
    const std::string tag = adaptive ? " (adaptive coefficients)" : " (standard coefficients)";
    s.require("nelder–mead: Rosenbrock from (−1.2, 1) converges" + tag, r.converged, std::to_string(r.evaluations) + " evaluations");
    s.check("nelder–mead: to x₁ = 1" + tag, r.x[0] - 1, 1e-6);
    s.check("nelder–mead: and x₂ = 1" + tag, r.x[1] - 1, 1e-6);
    s.check("nelder–mead: where the value is 0" + tag, r.value, 1e-12);
  }

  // f = (x − a)ᵀ M (x − a) with M positive definite, minimum 0 at a.
  {
    const double a[3] = {1.5, -0.5, 2.0};
    auto quad = [&](std::span<const double> x) {
      const double d0 = x[0] - a[0], d1 = x[1] - a[1], d2 = x[2] - a[2];
      return 4 * d0 * d0 + 3 * d1 * d1 + 2 * d2 * d2 + 2 * d0 * d1 + d1 * d2;
    };
    const auto r = nelder_mead(quad, {0, 0, 0}, tight);
    double err = 0;
    for (int i = 0; i < 3; i++) err = std::max(err, std::fabs(r.x[i] - a[i]));
    s.check("nelder–mead: the minimiser of a positive definite quadratic in three variables", err, 1e-6);
  }

  // One variable: against Brent's minimiser.
  {
    auto f1 = [](std::span<const double> x) { return x[0] * std::exp(x[0]); };
    const auto r = nelder_mead(f1, {0.5}, tight);
    const auto b = brent_minimise([](double t) { return t * std::exp(t); }, -4, 1, 1e-12);
    s.check("nelder–mead: one variable, x eˣ, against Brent's minimiser", r.x[0] - b.x, 1e-6);
  }

  // Powell's singular function (Moré, Garbow and Hillstrom 1981, problem 13): minimum 0 at 0.
  {
    auto powell = [](std::span<const double> x) {
      const double a = x[0] + 10 * x[1], b = x[2] - x[3], c = x[1] - 2 * x[2], d = x[0] - x[3];
      return a * a + 5 * b * b + c * c * c * c + 10 * d * d * d * d;
    };
    NelderMeadOptions o = tight;
    o.adaptive = true;
    o.max_evaluations = 20000;
    const auto r = nelder_mead(powell, {3, -1, 0, 1}, o);
    s.check("nelder–mead: Powell's singular function from (3, −1, 0, 1), its value", r.value, 1e-8);
  }

  {
    NelderMeadOptions o;
    o.max_evaluations = 50;
    const auto r = nelder_mead(rosenbrock, {-1.2, 1}, o);
    s.require("nelder–mead: the budget of evaluations is honoured and reported", !r.converged && r.evaluations <= 52, std::to_string(r.evaluations));
  }

  {
    // NaN outside x > 0: the method stays inside.
    auto f = [](std::span<const double> x) { return x[0] > 0 ? (std::log(x[0]) - 1) * (std::log(x[0]) - 1) : std::nan(""); };
    const auto r = nelder_mead(f, {0.1}, tight);
    s.check("nelder–mead: NaN taken as +∞; the minimum of (log x − 1)² at e", r.x[0] - std::exp(1.0), 1e-6);
  }

  {
    NelderMeadOptions o = tight;
    o.initial_step = {1.0};
    const auto r = nelder_mead(rosenbrock, {0, 0}, o);
    s.check("nelder–mead: Rosenbrock from the origin with unit initial steps", std::hypot(r.x[0] - 1, r.x[1] - 1), 1e-6);
  }
  return s.done();
}
