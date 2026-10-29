// Checks of calx/linalg's Cholesky factorisations: L Lᵀ = A on random positive definite
// matrices and the solve against Gaussian elimination; the Hilbert matrix of order 5, whose
// determinant ∏ L_ii² is 1/266716800000 exactly; a matrix that is not positive definite,
// refused; the band
// factorisation of a pentadiagonal matrix of the Hodrick–Prescott kind against the dense solve,
// and of the second-difference matrix against its exact inverse.

#include <calx/linalg/cholesky.hpp>
#include <calx/linalg/dense.hpp>
#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <vector>

using namespace calx::linalg;

int main() {
  calx::verification::Suite s;

  // ---- dense ----
  {
    calx::random::Stream rng(20261024, 9, 0, 0);
    const std::size_t n = 12;
    Matrix B(n);
    for (double& v : B.a) v = rng.normal();
    Matrix A(n);
    for (std::size_t i = 0; i < n; i++)
      for (std::size_t j = 0; j < n; j++) {
        double v = i == j ? 1.0 : 0.0;
        for (std::size_t k = 0; k < n; k++) v += B(i, k) * B(j, k);
        A(i, j) = v;
      }
    const Cholesky f = cholesky(A);
    s.require("cholesky: B Bᵀ + I is positive definite", f.positive_definite);
    double recon = 0;
    for (std::size_t i = 0; i < n; i++)
      for (std::size_t j = 0; j < n; j++) {
        double v = 0;
        for (std::size_t k = 0; k <= std::min(i, j); k++) v += f.l(i, k) * f.l(j, k);
        recon = std::max(recon, std::fabs(v - A(i, j)));
      }
    s.check("cholesky: L Lᵀ = A on a random positive definite matrix of order 12", recon, 1e-12);
    std::vector<double> b(n);
    Matrix Bm(n);
    for (std::size_t i = 0; i < n; i++) b[i] = Bm(i, 0) = std::sin(static_cast<double>(i));
    const auto x = cholesky_solve(f, b);
    const Matrix xs = solve(A, Bm);
    double diff = 0;
    for (std::size_t i = 0; i < n; i++) diff = std::max(diff, std::fabs(x[i] - xs(i, 0)));
    s.check("cholesky: the solve against Gaussian elimination", diff, 1e-12);
  }

  // ---- the Hilbert matrix of order 5: det = 1 / 266716800000 ----
  {
    Matrix H(5);
    for (std::size_t i = 0; i < 5; i++)
      for (std::size_t j = 0; j < 5; j++) H(i, j) = 1.0 / static_cast<double>(i + j + 1);
    const Cholesky f = cholesky(H);
    double det = 1;
    for (std::size_t i = 0; i < 5; i++) det *= f.l(i, i) * f.l(i, i);
    s.require("cholesky: the Hilbert matrix of order 5 is positive definite", f.positive_definite);
    s.check("cholesky: its determinant, the square of ∏ L_ii, against 1/266716800000", calx::verification::relative(det, 1.0 / 266716800000.0), 1e-10);
  }

  // ---- not positive definite ----
  {
    Matrix A(3);
    A.a = {1, 2, 0, 2, 1, 0, 0, 0, 1};  // eigenvalues 3, −1, 1
    s.require("cholesky: an indefinite matrix is refused", !cholesky(A).positive_definite);
  }

  // ---- band: the Hodrick–Prescott matrix I + λ DᵀD, half-bandwidth 2 ----
  {
    const std::size_t n = 40;
    const double lambda = 1600;
    SymmetricBand A(n, 2);
    Matrix D(n);
    // DᵀD for the second differences: rows (1, −2, 1).
    for (std::size_t t = 0; t + 2 < n; t++) {
      const double d[3] = {1, -2, 1};
      for (int a = 0; a < 3; a++)
        for (int b = 0; b < 3; b++) D(t + a, t + b) += lambda * d[a] * d[b];
    }
    for (std::size_t i = 0; i < n; i++) {
      D(i, i) += 1;
      for (std::size_t j = (i >= 2 ? i - 2 : 0); j <= i; j++) A(i, j) = D(i, j);
    }
    std::vector<double> y(n);
    Matrix Y(n);
    for (std::size_t i = 0; i < n; i++) y[i] = Y(i, 0) = std::cos(0.3 * static_cast<double>(i)) + 0.01 * static_cast<double>(i * i);
    const auto tau = band_solve(A, y);
    const Matrix ts = solve(D, Y);
    double diff = 0, scale = 0;
    for (std::size_t i = 0; i < n; i++) {
      diff = std::max(diff, std::fabs(tau[i] - ts(i, 0)));
      scale = std::max(scale, std::fabs(ts(i, 0)));
    }
    s.check("band: the pentadiagonal system I + 1600 DᵀD against the dense solve (relative)", diff / scale, 1e-11);
    s.require("band: the entry outside the band reads as zero", A.at(10, 2) == 0 && A.at(2, 4) == A(4, 2));
  }

  // ---- band: the second-difference matrix, inverse (T⁻¹)_ij = min(i, j)(n + 1 − max(i, j)) / (n + 1) ----
  {
    const std::size_t n = 50;
    SymmetricBand T(n, 1);
    for (std::size_t i = 0; i < n; i++) {
      T(i, i) = 2;
      if (i > 0) T(i, i - 1) = -1;
    }
    std::vector<double> e(n, 0.0);
    e[19] = 1;  // column 20 of the inverse (1-based)
    const auto col = band_solve(T, e);
    double err = 0;
    for (std::size_t i = 0; i < n; i++) {
      const double ii = static_cast<double>(i + 1), jj = 20, nn = static_cast<double>(n);
      err = std::max(err, std::fabs(col[i] - std::min(ii, jj) * (nn + 1 - std::max(ii, jj)) / (nn + 1)));
    }
    s.check("band: a column of the inverse of the second-difference matrix against its closed form", err, 1e-12);
  }
  return s.done();
}
