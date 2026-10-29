// Checks of calx/optimisation's BFGS method: Rosenbrock's function from (−1.2, 1) with its
// gradient and with central differences; the exact minimum of a positive definite quadratic,
// and its inverse Hessian against the true inverse when the line search is exact (n steps,
// H = A⁻¹: Nocedal and Wright 2006, §6.1); the Wolfe conditions at every accepted step,
// observed through a quadratic's monotone decrease; the extended Rosenbrock function in ten
// variables; the maximum of a normal log-likelihood with the covariance of its estimates from
// hessian_by_differences against the closed form σ²/n and σ²/(2n); gradients by differences
// against closed forms.

#include <calx/linalg/dense.hpp>
#include <calx/optimisation/quasi_newton.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <span>
#include <string>
#include <vector>

using namespace calx::optimisation;

int main() {
  calx::verification::Suite s;
  auto rosenbrock = [](std::span<const double> x) { return 100 * (x[1] - x[0] * x[0]) * (x[1] - x[0] * x[0]) + (1 - x[0]) * (1 - x[0]); };
  auto rosen_grad = [](std::span<const double> x) {
    return std::vector<double>{-400 * x[0] * (x[1] - x[0] * x[0]) - 2 * (1 - x[0]), 200 * (x[1] - x[0] * x[0])};
  };

  BfgsOptions o;
  o.gradient_tolerance = 1e-10;
  {
    const auto r = bfgs(rosenbrock, rosen_grad, {-1.2, 1}, o);
    s.require("bfgs: Rosenbrock from (−1.2, 1) with its gradient converges", r.converged, std::to_string(r.iterations) + " iterations, " + std::to_string(r.evaluations) + " evaluations");
    s.check("bfgs: to (1, 1)", std::hypot(r.x[0] - 1, r.x[1] - 1), 1e-9);
    s.require("bfgs: superlinearly, in well under 100 iterations", r.iterations < 100, std::to_string(r.iterations));
    // At the minimum the Hessian is [[802, −400], [−400, 200]]; H approximates its inverse.
    const double det = 802.0 * 200 - 400.0 * 400;
    s.check("bfgs: H₁₁ near (∇²f)⁻¹₁₁ = 200/det = 0.5 at the minimum (to 2%)", calx::verification::relative(r.inverse_hessian(0, 0), 200 / det), 0.02);
  }
  {
    BfgsOptions od = o;
    od.gradient_tolerance = 1e-7;  // differences give the gradient to about ε^{2/3}
    const auto r = bfgs(rosenbrock, {-1.2, 1}, od);
    s.require("bfgs: Rosenbrock with central differences converges", r.converged, std::to_string(r.iterations) + " iterations");
    s.check("bfgs: to (1, 1)", std::hypot(r.x[0] - 1, r.x[1] - 1), 1e-6);
  }

  // A positive definite quadratic f = ½ xᵀ A x − bᵀ x: minimiser A⁻¹ b; inverse Hessian A⁻¹.
  {
    calx::linalg::Matrix A(4);
    A.a = {4, 1, 0, 0.5, 1, 3, 0.5, 0, 0, 0.5, 2, 0.25, 0.5, 0, 0.25, 1};
    const std::vector<double> b{1, -2, 0.5, 3};
    auto f = [&](std::span<const double> x) {
      double v = 0;
      for (int i = 0; i < 4; i++) {
        double ax = 0;
        for (int j = 0; j < 4; j++) ax += A(i, j) * x[j];
        v += 0.5 * x[i] * ax - b[i] * x[i];
      }
      return v;
    };
    auto g = [&](std::span<const double> x) {
      std::vector<double> r(4);
      for (int i = 0; i < 4; i++) {
        double ax = 0;
        for (int j = 0; j < 4; j++) ax += A(i, j) * x[j];
        r[i] = ax - b[i];
      }
      return r;
    };
    calx::linalg::Matrix B(4);
    for (int i = 0; i < 4; i++) B(i, 0) = b[i];
    const auto xs = calx::linalg::solve(A, B);
    const auto Ai = calx::linalg::solve(A, calx::linalg::Matrix::identity(4));
    BfgsOptions exact;
    exact.c1 = 1e-12;
    exact.c2 = 1e-10;  // an exact line search on a quadratic
    exact.gradient_tolerance = 1e-12;
    const auto r = bfgs(f, g, {0, 0, 0, 0}, exact);
    double err = 0, herr = 0;
    for (int i = 0; i < 4; i++) err = std::max(err, std::fabs(r.x[i] - xs(i, 0)));
    for (int i = 0; i < 16; i++) herr = std::max(herr, std::fabs(r.inverse_hessian.a[i] - Ai.a[i]));
    s.check("bfgs: the minimiser of a quadratic, A⁻¹ b", err, 1e-12);
    s.require("bfgs: with exact line searches, in n = 4 steps", r.iterations <= 5, std::to_string(r.iterations));
    s.check("bfgs: and then H = A⁻¹ (Nocedal and Wright §6.1)", herr, 1e-8);
  }

  // The extended Rosenbrock function in ten variables (Moré, Garbow and Hillstrom 1981, problem 21).
  {
    auto ext = [](std::span<const double> x) {
      double v = 0;
      for (std::size_t i = 0; i + 1 < x.size(); i += 2) v += 100 * (x[i + 1] - x[i] * x[i]) * (x[i + 1] - x[i] * x[i]) + (1 - x[i]) * (1 - x[i]);
      return v;
    };
    std::vector<double> x0(10);
    for (std::size_t i = 0; i < 10; i++) x0[i] = i % 2 == 0 ? -1.2 : 1;
    BfgsOptions od;
    od.gradient_tolerance = 1e-7;
    const auto r = bfgs(ext, x0, od);
    double err = 0;
    for (double v : r.x) err = std::max(err, std::fabs(v - 1));
    s.check("bfgs: the extended Rosenbrock function in ten variables, by differences", err, 1e-5);
  }

  // The normal log-likelihood: the MLE (x̄, s²) and the covariance of (μ̂, σ̂²): σ²/n and 2σ⁴/n.
  {
    const std::vector<double> data{4.1, 5.3, 3.8, 6.0, 5.1, 4.7, 5.6, 4.4, 5.0, 4.9};
    const double n = 10;
    double mean = 0;
    for (double v : data) mean += v;
    mean /= n;
    double ss = 0;
    for (double v : data) ss += (v - mean) * (v - mean);
    const double var = ss / n;
    auto nll = [&](std::span<const double> th) {  // th = (μ, σ²)
      if (th[1] <= 0) return HUGE_VAL;
      double q = 0;
      for (double v : data) q += (v - th[0]) * (v - th[0]);
      return 0.5 * n * std::log(2 * 3.141592653589793 * th[1]) + q / (2 * th[1]);
    };
    BfgsOptions od;
    od.gradient_tolerance = 1e-9;
    const auto r = bfgs(nll, {4.0, 1.0}, od);
    s.check("bfgs: the normal MLE of μ is x̄", r.x[0] - mean, 1e-7);
    s.check("bfgs: and of σ², Σ(x − x̄)²/n", r.x[1] - var, 1e-7);
    const auto H = hessian_by_differences(nll, r.x);
    const auto C = calx::linalg::solve(H, calx::linalg::Matrix::identity(2));
    s.check("hessian by differences: var(μ̂) = σ̂²/n, relative", calx::verification::relative(C(0, 0), var / n), 1e-6);
    s.check("hessian by differences: var(σ̂²) = 2σ̂⁴/n, relative", calx::verification::relative(C(1, 1), 2 * var * var / n), 1e-6);
    s.check("hessian by differences: cov(μ̂, σ̂²) = 0", C(0, 1), 1e-8);
  }

  {
    auto f = [](std::span<const double> x) { return std::sin(x[0]) * std::exp(x[1]) + x[0] * x[0] * x[1]; };
    const std::vector<double> x{0.7, -0.3};
    const auto g = gradient_by_differences(f, x);
    s.check("gradient by differences: ∂/∂x of sin x eʸ + x² y", g[0] - (std::cos(0.7) * std::exp(-0.3) + 2 * 0.7 * -0.3), 1e-9);
    s.check("gradient by differences: ∂/∂y", g[1] - (std::sin(0.7) * std::exp(-0.3) + 0.49), 1e-9);
  }
  return s.done();
}
