// Checks of calx/integrators: the observed order of convergence of each fixed-step method
// on a problem with a closed-form solution (the damped oscillator), and the adaptive
// Dormand-Prince pair against its tolerance.

#include <calx/integrators/runge_kutta.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <string>

namespace ci = calx::integrators;

namespace {

// x'' + 2 z w x' + w^2 x = 0 as a first-order system, z = 0.1, w = 2; x(0) = 1, x'(0) = 0.
constexpr double zeta = 0.1, omega = 2;
void field(const ci::Vec<2>& x, ci::Vec<2>& dx) {
  dx[0] = x[1];
  dx[1] = -2 * zeta * omega * x[1] - omega * omega * x[0];
}
double exact(double t) {
  const double wd = omega * std::sqrt(1 - zeta * zeta);
  return std::exp(-zeta * omega * t) * (std::cos(wd * t) + zeta * omega / wd * std::sin(wd * t));
}

template <class Step>
double error_at(Step step, int n, double T = 2) {
  ci::Vec<2> x{1, 0};
  const double h = T / n;
  for (int k = 0; k < n; k++) step(x, h);
  return std::fabs(x[0] - exact(T));
}

}  // namespace

int main() {
  calx::verification::Suite s;
  auto f = [](const ci::Vec<2>& x, ci::Vec<2>& dx) { field(x, dx); };
  struct Method { const char* name; double order; double (*err)(int); };
  const Method methods[] = {
      {"euler", 1, [](int n) { return error_at([](ci::Vec<2>& x, double h) { ci::euler<2>(field, x, h); }, n); }},
      {"heun", 2, [](int n) { return error_at([](ci::Vec<2>& x, double h) { ci::heun<2>(field, x, h); }, n); }},
      {"rk4", 4, [](int n) { return error_at([](ci::Vec<2>& x, double h) { ci::rk4<2>(field, x, h); }, n); }},
  };
  for (const auto& m : methods) {
    const int n = m.order >= 4 ? 100 : 2000;
    const double e1 = m.err(n), e2 = m.err(2 * n);
    s.check(std::string("integrators: ") + m.name + ", observed order minus " + calx::verification::shortest(m.order) + " (damped oscillator)",
            calx::verification::order(e1, e2) - m.order, 0.05);
  }
  {
    ci::Vec<2> x{1, 0};
    double sub = 0;
    int steps = 0;
    for (int k = 0; k < 20; k++) steps += ci::dopri5<2>(f, x, 0.1, {1, 1}, {1e-10, 1e-10}, sub);
    s.check("integrators: dopri5 at tolerance 1e-10, error at t = 2", x[0] - exact(2), 1e-8);
    s.require("integrators: dopri5 takes substeps as its tolerance asks", steps > 20, std::to_string(steps) + " substeps in 20 steps");
  }
  return s.done();
}
