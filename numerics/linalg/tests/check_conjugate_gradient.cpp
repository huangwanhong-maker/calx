// Checks of conjugate gradients (conjugate_gradient.hpp).
//
// A random symmetric positive definite system of order 60 solved as calx's dense solve solves it;
// the discrete Laplacian of a 200-point path with fixed ends against the exact solution of
// −u″ = 1 (the parabola x(1 − x)/2, which the second difference reproduces exactly); the singular
// Laplacian of a ring with a right-hand side of zero sum, whose least-norm solution has zero sum
// and agrees with the dense solve of the system made definite by the kernel; a zero right-hand
// side gives zero at once.

#include <calx/linalg/conjugate_gradient.hpp>
#include <calx/linalg/dense.hpp>
#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <cstdint>
#include <span>
#include <vector>

namespace cl = calx::linalg;

namespace {

struct Triplets {
  std::vector<std::uint32_t> r, c;
  std::vector<double> v;
  void add(std::uint32_t i, std::uint32_t j, double x) {
    r.push_back(i);
    c.push_back(j);
    v.push_back(x);
  }
  std::pair<cl::Csr<>, std::vector<double>> csr(std::size_t n) const {
    std::vector<std::uint64_t> order;
    auto A = cl::compress<std::uint32_t>(n, n, r, c, &order);
    std::vector<double> a(order.size());
    for (std::size_t k = 0; k < order.size(); k++) a[k] = v[order[k]];
    return {A, a};
  }
};

}  // namespace

int main() {
  calx::verification::Suite s;
  calx::random::Stream z(20261024, 0x63677366u, 0, 0);  // "cgsf"

  // A random SPD system: B Bᵀ + n I.
  {
    const std::size_t n = 60;
    cl::Matrix B(n), M(n), rhs(n);
    for (double& x : B.a) x = z.uniform() - 0.5;
    Triplets T;
    for (std::size_t i = 0; i < n; i++)
      for (std::size_t j = 0; j < n; j++) {
        double s = i == j ? static_cast<double>(n) * 0.1 : 0.0;
        for (std::size_t k = 0; k < n; k++) s += B(i, k) * B(j, k);
        M(i, j) = s;
        T.add(static_cast<std::uint32_t>(i), static_cast<std::uint32_t>(j), s);
      }
    std::vector<double> b(n);
    for (std::size_t i = 0; i < n; i++) rhs(i, 0) = b[i] = z.uniform() - 0.5;
    const auto X = cl::solve(M, rhs);
    const auto [A, a] = T.csr(n);
    const auto R = cl::conjugate_gradient(A, std::span<const double>(a), std::span<const double>(b), 1e-13);
    double worst = 0;
    for (std::size_t i = 0; i < n; i++) worst = std::max(worst, std::fabs(R.x[i] - X(i, 0)));
    s.require("a random SPD system of order 60: converged", R.converged, std::to_string(R.iterations) + " iterations");
    s.check("... to the dense solve", worst, 1e-10);
  }
  // −u″ = 1 on (0, 1) with u(0) = u(1) = 0: the parabola x(1 − x)/2, exact for second differences.
  {
    const std::size_t n = 200;
    const double h = 1.0 / static_cast<double>(n + 1);
    Triplets T;
    std::vector<double> b(n, h * h);
    for (std::size_t i = 0; i < n; i++) {
      T.add(static_cast<std::uint32_t>(i), static_cast<std::uint32_t>(i), 2);
      if (i > 0) T.add(static_cast<std::uint32_t>(i), static_cast<std::uint32_t>(i - 1), -1);
      if (i + 1 < n) T.add(static_cast<std::uint32_t>(i), static_cast<std::uint32_t>(i + 1), -1);
    }
    const auto [A, a] = T.csr(n);
    const auto R = cl::conjugate_gradient(A, std::span<const double>(a), std::span<const double>(b), 1e-12);
    double worst = 0;
    for (std::size_t i = 0; i < n; i++) {
      const double x = h * static_cast<double>(i + 1);
      worst = std::max(worst, std::fabs(R.x[i] - x * (1 - x) / 2));
    }
    s.require("the Laplacian of a path, −u″ = 1: converged within n steps", R.converged && R.iterations <= n + 1, std::to_string(R.iterations) + " iterations");
    s.check("... to the parabola x(1 − x)/2", worst, 1e-11);
  }
  // A singular Laplacian: a ring of 50, with a right-hand side of zero sum.
  {
    const std::size_t n = 50;
    Triplets T;
    std::vector<double> b(n);
    double mean = 0;
    for (double& x : b) mean += x = z.uniform() - 0.5;
    mean /= static_cast<double>(n);
    for (double& x : b) x -= mean;
    for (std::size_t i = 0; i < n; i++) {
      const auto u = static_cast<std::uint32_t>(i), v = static_cast<std::uint32_t>((i + 1) % n);
      const double w = 0.5 + z.uniform();
      T.add(u, u, w), T.add(v, v, w), T.add(u, v, -w), T.add(v, u, -w);
    }
    const auto [A, a] = T.csr(n);
    const std::vector<std::vector<double>> kernel{std::vector<double>(n, 1 / std::sqrt(static_cast<double>(n)))};
    const auto R = cl::conjugate_gradient(A, std::span<const double>(a), std::span<const double>(b), 1e-13, 0, kernel);
    double sum = 0;
    for (double x : R.x) sum += x;
    s.require("a singular ring Laplacian, b of zero sum: converged", R.converged, std::to_string(R.relative_residual));
    s.check("... to the solution of least norm (zero sum)", sum, 1e-10);
    // Against the dense solve of A + 11ᵀ/n, definite, with the same solution in the complement.
    cl::Matrix D(n), rhs(n);
    for (std::size_t i = 0; i < n; i++)
      for (std::uint64_t k = A.row[i]; k < A.row[i + 1]; k++) D(i, A.col[k]) += a[k];
    for (std::size_t i = 0; i < n; i++) {
      for (std::size_t j = 0; j < n; j++) D(i, j) += 1.0 / static_cast<double>(n);
      rhs(i, 0) = b[i];
    }
    const auto X = cl::solve(D, rhs);
    double worst = 0;
    for (std::size_t i = 0; i < n; i++) worst = std::max(worst, std::fabs(R.x[i] - X(i, 0)));
    s.check("... and agrees with the dense solve made definite by the kernel", worst, 1e-9);
    const std::vector<double> zero(n, 0.0);
    const auto Z = cl::conjugate_gradient(A, std::span<const double>(a), std::span<const double>(zero));
    s.require("a zero right-hand side gives zero at once", Z.converged && Z.iterations == 0);
  }
  return s.done();
}
