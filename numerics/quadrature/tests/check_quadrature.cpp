// Checks of calx/quadrature: the Gauss-Kronrod pair integrates polynomials to the degrees
// it should (Kronrod 15 to degree 22, Gauss 7 to degree 13) and its weights sum to the
// interval's length; the adaptive scheme reaches its tolerance on smooth and on peaked and
// singular-derivative integrands with known integrals.

#include <calx/quadrature/gauss_kronrod.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <numbers>

namespace cq = calx::quadrature;

int main() {
  calx::verification::Suite s;
  {
    // On [-1, 1], x^22 integrates to 2/23 (K15 exact); x^14 is the first power G7 misses,
    // so |K15 - G7| for x^14 is nonzero while for x^13 (odd) and x^12 it vanishes.
    const auto k22 = cq::gauss_kronrod_15([](double x) { return std::pow(x, 22); }, -1, 1);
    s.check("quadrature: K15 integrates x^22 exactly on [-1, 1]", calx::verification::relative(k22[0], 2.0 / 23), 1e-14);
    const auto k12 = cq::gauss_kronrod_15([](double x) { return std::pow(x, 12); }, -1, 1);
    s.check("quadrature: G7 integrates x^12 exactly (|K15 - G7| vanishes)", k12[1], 1e-15);
    const auto k14 = cq::gauss_kronrod_15([](double x) { return std::pow(x, 14); }, -1, 1);
    s.require("quadrature: G7 does not integrate x^14 exactly", k14[1] > 1e-6);
    const auto one = cq::gauss_kronrod_15([](double) { return 1.0; }, 2, 5);
    s.check("quadrature: the Kronrod weights sum to the length of the interval", one[0] - 3, 1e-14);
  }
  {
    const auto r = cq::integrate([](double x) { return std::sin(x); }, 0, std::numbers::pi);
    s.check("quadrature: the integral of sin over [0, pi] is 2", r.value - 2, 1e-13);
    const auto g = cq::integrate([](double x) { return std::exp(-x * x); }, -6, 6);
    s.check("quadrature: the integral of exp(-x^2) over [-6, 6] is sqrt(pi) erf(6)", g.value - std::sqrt(std::numbers::pi) * std::erf(6.0), 1e-13);
    // A narrow peak: 1 / (1e-4 + x^2) over [-1, 1] = 2 / sqrt(1e-4) atan(1 / sqrt(1e-4)).
    const double e = 1e-4;
    const auto p = cq::integrate([e](double x) { return 1 / (e + x * x); }, -1, 1, 1e-12, 1e-12);
    s.check("quadrature: a narrow Lorentzian peak, relative error", calx::verification::relative(p.value, 2 / std::sqrt(e) * std::atan(1 / std::sqrt(e))), 1e-11);
    // sqrt(x) on [0, 1]: a singular derivative at 0; the integral is 2/3.
    const auto q = cq::integrate([](double x) { return std::sqrt(x); }, 0, 1, 1e-12, 1e-12);
    s.check("quadrature: sqrt(x) on [0, 1], with its singular derivative, is 2/3", q.value - 2.0 / 3, 1e-11);
    s.require("quadrature: each adaptive integral meets its tolerance", r.converged && g.converged && p.converged && q.converged);
  }
  return s.done();
}
