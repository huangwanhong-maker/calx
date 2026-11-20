// Checks of calx/linalg/band: the band LU's solutions against the dense solver on random band matrices of several shapes;
// a matrix whose diagonal is zero (pivoting needed, and its fill reaching ku + kl); a singular one refused; the transpose
// and the product against their definitions; a tridiagonal system against the closed form of the discrete Poisson problem.

#include <calx/linalg/band.hpp>
#include <calx/linalg/dense.hpp>
#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <string>

using namespace calx::linalg;

int main() {
  calx::verification::Suite s;
  calx::random::Stream r(20, 1, 0, 0);
  const std::size_t shapes[][3] = {{40, 1, 1}, {60, 2, 3}, {50, 4, 1}, {30, 0, 2}, {35, 3, 0}, {80, 5, 5}};
  for (const auto& sh : shapes) {
    const std::size_t n = sh[0], kl = sh[1], ku = sh[2];
    Band B(n, kl, ku);
    Matrix D(n);
    for (std::size_t i = 0; i < n; i++)
      for (std::size_t j = (i >= kl ? i - kl : 0); j <= std::min(n - 1, i + ku); j++) {
        const double v = r.uniform() * 2 - 1 + (i == j ? 0.2 : 0.0);
        B(i, j) = v, D(i, j) = v;
      }
    Matrix bm(n);
    std::vector<double> b(n);
    for (std::size_t i = 0; i < n; i++) b[i] = r.uniform() - 0.5, bm(i, 0) = b[i];
    const auto x = band_solve(B, b);
    const Matrix X = solve(D, bm);
    double worst = 0;
    for (std::size_t i = 0; i < n; i++) worst = std::max(worst, std::fabs(x[i] - X(i, 0)) / (1 + std::fabs(X(i, 0))));
    s.check("a random band matrix (n " + std::to_string(n) + ", kl " + std::to_string(kl) + ", ku " + std::to_string(ku) + "): the band LU against the dense solver", worst, 1e-9);
    // The residual relative to ‖A‖∞‖x‖∞ (a random triangular matrix is badly conditioned and its solution large).
    const auto y = multiply(B, x);
    double res = 0, xn = 0, an = 0;
    for (std::size_t i = 0; i < n; i++) {
      res = std::max(res, std::fabs(y[i] - b[i])), xn = std::max(xn, std::fabs(x[i]));
      double row = 0;
      for (std::size_t j = 0; j < n; j++) row += std::fabs(D(i, j));
      an = std::max(an, row);
    }
    s.check("... and its residual relative to the size of A and x", res / (an * xn), 1e-14);
  }
  {
    // Zero diagonal: every step must pivot.
    const std::size_t n = 12;
    Band B(n, 1, 1);
    Matrix D(n);
    for (std::size_t i = 0; i < n; i++) {
      if (i > 0) B(i, i - 1) = 1 + 0.1 * i, D(i, i - 1) = 1 + 0.1 * i;
      if (i + 1 < n) B(i, i + 1) = 2 - 0.05 * i, D(i, i + 1) = 2 - 0.05 * i;
    }
    std::vector<double> b(n);
    Matrix bm(n);
    for (std::size_t i = 0; i < n; i++) b[i] = std::sin(1.0 + i), bm(i, 0) = b[i];
    const auto F = band_lu(B);
    std::size_t swaps = 0;
    for (std::size_t k = 0; k < n; k++) swaps += F.pivot[k] != k;
    const auto x = band_solve(F, b);
    const Matrix X = solve(D, bm);
    double worst = 0;
    for (std::size_t i = 0; i < n; i++) worst = std::max(worst, std::fabs(x[i] - X(i, 0)));
    s.require("a zero diagonal: the factorisation pivots", swaps > 0);
    s.check("... and solves as the dense solver does", worst, 1e-12);
  }
  {
    Band B(4, 1, 1);
    B(0, 0) = 1, B(0, 1) = 2, B(1, 0) = 2, B(1, 1) = 4, B(2, 2) = 1, B(3, 3) = 1;
    bool refused = false;
    try {
      band_lu(B);
    } catch (const std::runtime_error&) {
      refused = true;
    }
    s.require("a singular band matrix is refused", refused);
  }
  {
    Band B(7, 2, 1);
    for (std::size_t i = 0; i < 7; i++)
      for (std::size_t j = (i >= 2 ? i - 2 : 0); j <= std::min<std::size_t>(6, i + 1); j++) B(i, j) = 10 * i + j + 1;
    const Band T = transpose(B);
    const Band& A = B;  // read through a const view: zero outside the band
    bool same = T.kl == 1 && T.ku == 2;
    for (std::size_t i = 0; i < 7; i++)
      for (std::size_t j = 0; j < 7; j++) same = same && T(i, j) == A(j, i);
    s.require("the transpose: Aᵀ(i, j) = A(j, i), bandwidths exchanged", same);
  }
  {
    // −u'' = 1 on (0, 1), u(0) = u(1) = 0, at n interior points: the discrete solution equals x(1 − x)/2 exactly.
    const std::size_t n = 99;
    const double h = 1.0 / (n + 1);
    Band B(n, 1, 1);
    for (std::size_t i = 0; i < n; i++) {
      B(i, i) = 2 / (h * h);
      if (i > 0) B(i, i - 1) = -1 / (h * h);
      if (i + 1 < n) B(i, i + 1) = -1 / (h * h);
    }
    const auto u = band_solve(B, std::vector<double>(n, 1.0));
    double worst = 0;
    for (std::size_t i = 0; i < n; i++) {
      const double x = (i + 1) * h;
      worst = std::max(worst, std::fabs(u[i] - x * (1 - x) / 2));
    }
    s.check("the discrete Poisson problem: x(1 − x)/2 at the nodes, exactly", worst, 1e-11);
  }
  return s.done();
}
