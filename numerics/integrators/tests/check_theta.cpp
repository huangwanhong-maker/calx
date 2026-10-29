// Checks of the θ-method (theta.hpp) against closed forms: one step on a diagonal system multiplies each
// component by R(hλ) exactly; on a coupled system the error against the matrix exponential falls as h
// (θ = 0, 1) or h² (θ = ½); on a stiff system θ ≥ ½ is bounded where explicit Euler grows; a constant
// forcing is reached in one implicit step of great length; a row without its diagonal is refused.

#include <calx/integrators/theta.hpp>
#include <calx/linalg/dense.hpp>
#include <calx/linalg/expm.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

namespace ci = calx::integrators;
namespace cl = calx::linalg;

namespace {

struct System {
  cl::Csr<std::uint32_t> A;
  std::vector<double> a;
};

// A symmetric matrix from its dense entries, every diagonal entry stored.
System sparse(const std::vector<std::vector<double>>& dense) {
  std::vector<std::uint32_t> r, c;
  std::vector<double> v;
  for (std::size_t i = 0; i < dense.size(); i++)
    for (std::size_t j = 0; j < dense.size(); j++)
      if (dense[i][j] != 0 || i == j) {
        r.push_back(static_cast<std::uint32_t>(i));
        c.push_back(static_cast<std::uint32_t>(j));
        v.push_back(dense[i][j]);
      }
  System s;
  std::vector<std::uint64_t> order;
  s.A = cl::compress<std::uint32_t>(dense.size(), dense.size(), r, c, &order);
  s.a.resize(v.size());
  for (std::size_t k = 0; k < v.size(); k++) s.a[k] = v[order[k]];
  return s;
}

}  // namespace

int main() {
  calx::verification::Suite s;
  const std::vector<double> thetas{0.0, 0.5, 1.0};

  // One step on a diagonal system: R(hλ) per component.
  {
    const auto S = sparse({{1, 0, 0}, {0, 4, 0}, {0, 0, 9}});
    const std::vector<double> m{1, 2, 3}, y{1, -2, 0.5}, zero(3, 0.0);
    const double h = 0.1;
    double worst = 0;
    for (double th : thetas) {
      const auto step = ci::theta_step(S.A, std::span<const double>(S.a), m, y, zero, h, th, 1e-15);
      const std::vector<double> lambda{1.0 / 1, 4.0 / 2, 9.0 / 3};
      for (int i = 0; i < 3; i++) worst = std::max(worst, std::fabs(step.y[i] / y[i] - ci::theta_factor(h * lambda[i], th)));
    }
    s.check("theta: one step on a diagonal system multiplies each component by R(hλ), θ = 0, ½, 1", worst, 1e-14);
  }

  // A coupled system M ẏ = −A y against y(T) = exp(−T M⁻¹A) y(0): the order of the error.
  {
    const std::vector<std::vector<double>> A{{2, -1, 0}, {-1, 2, -1}, {0, -1, 2}};
    const auto S = sparse(A);
    const std::vector<double> m{1, 0.5, 2}, y0{1, 0, -1}, zero(3, 0.0);
    const double T = 1;
    cl::Matrix G(3);
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) G(i, j) = -T * A[i][j] / m[i];
    const cl::Matrix E = cl::expm(G);
    std::vector<double> exact(3, 0.0);
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) exact[i] += E(i, j) * y0[j];
    auto error = [&](double th, int steps) {
      std::vector<double> y = y0;
      for (int k = 0; k < steps; k++) y = ci::theta_step(S.A, std::span<const double>(S.a), m, y, zero, T / steps, th, 1e-15).y;
      double e = 0;
      for (int i = 0; i < 3; i++) e = std::max(e, std::fabs(y[i] - exact[i]));
      return e;
    };
    for (double th : thetas) {
      const double order = std::log2(error(th, 200) / error(th, 400));
      const double expected = th == 0.5 ? 2 : 1;
      s.check("theta: the observed order of θ = " + std::to_string(th).substr(0, 3) + " against the matrix exponential, h = 1/200 → 1/400",
              order - expected, 0.02);
    }
  }

  // A stiff system, hλ = 10⁶: bounded for θ ≥ ½, growing for θ = 0.
  {
    const auto S = sparse({{1e6}});
    const std::vector<double> m{1}, y{1}, zero(1, 0.0);
    const auto cn = ci::theta_step(S.A, std::span<const double>(S.a), m, y, zero, 1.0, 0.5);
    const auto ie = ci::theta_step(S.A, std::span<const double>(S.a), m, y, zero, 1.0, 1.0);
    const auto ee = ci::theta_step(S.A, std::span<const double>(S.a), m, y, zero, 1.0, 0.0);
    s.require("theta: on a stiff mode Crank–Nicolson and implicit Euler are bounded, explicit Euler grows",
              std::fabs(cn.y[0]) <= 1 && std::fabs(ie.y[0]) < 1e-5 && std::fabs(ee.y[0]) > 1e5);
  }

  // A constant forcing: one implicit step of great length reaches A⁻¹ c.
  {
    const auto S = sparse({{2, -1}, {-1, 2}});
    const std::vector<double> m{1, 1}, y{0, 0}, c{1, 0};
    const auto step = ci::theta_step(S.A, std::span<const double>(S.a), m, y, c, 1e12, 1.0, 1e-15);
    s.check("theta: one implicit step of length 10¹² reaches the steady state A⁻¹c = (⅔, ⅓)", std::max(std::fabs(step.y[0] - 2.0 / 3), std::fabs(step.y[1] - 1.0 / 3)), 1e-11);
    s.require("theta: the solve converged, its iterations counted", step.converged && step.iterations >= 1 && step.relative_residual <= 1e-14);
  }

  // A row without its diagonal entry is refused.
  {
    const std::vector<std::uint32_t> r{0, 1}, c{1, 0};
    const auto A = cl::compress<std::uint32_t>(2, 2, r, c);
    const std::vector<double> a{-1, -1}, m{1, 1}, y{1, 1}, f{0, 0};
    bool threw = false;
    try { (void)ci::theta_step(A, std::span<const double>(a), m, y, f, 0.1, 0.5); } catch (const std::invalid_argument&) { threw = true; }
    s.require("theta: a matrix without a stored diagonal entry is refused", threw);
  }
  return s.done();
}
