// Checks of calx/linalg's QR factorisation and least squares: the factors of random matrices
// (Q orthonormal, R upper triangular with non-increasing diagonal, A P = Q R); least squares
// against a line fitted in closed form, against the normal equations on a well-conditioned
// problem, and the residual's orthogonality to the columns; Läuchli's matrix, on which the
// normal equations are singular in floating point and QR is not; a consistent square system
// against Gaussian elimination; the numerical rank of a matrix with a dependent column;
// (AᵀA)⁻¹ against the inverse of the Gram matrix; products and transposes of Rectangular.

#include <calx/linalg/qr.hpp>
#include <calx/linalg/rectangular.hpp>
#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <string>
#include <vector>

using namespace calx::linalg;

namespace {

Rectangular random_matrix(std::size_t m, std::size_t n, std::uint32_t id) {
  calx::random::Stream rng(20261024, 7, id, 0);
  Rectangular A(m, n);
  for (double& v : A.a) v = rng.normal();
  return A;
}

double max_abs(const std::vector<double>& v) {
  double m = 0;
  for (double x : v) m = std::max(m, std::fabs(x));
  return m;
}

}  // namespace

int main() {
  calx::verification::Suite s;

  // ---- the factors ----
  for (auto [m, n, id] : {std::tuple<std::size_t, std::size_t, std::uint32_t>{30, 5, 1}, {12, 12, 2}, {6, 9, 3}}) {
    const Rectangular A = random_matrix(m, n, id);
    const Qr f = qr(A);
    const std::size_t t = std::min(m, n);
    const Rectangular Q = thin_q(f, t);
    double orth = 0;
    for (std::size_t a = 0; a < t; a++)
      for (std::size_t b = 0; b < t; b++) {
        double d = 0;
        for (std::size_t i = 0; i < m; i++) d += Q(i, a) * Q(i, b);
        orth = std::max(orth, std::fabs(d - (a == b ? 1.0 : 0.0)));
      }
    const std::string shape = std::to_string(m) + " × " + std::to_string(n);
    s.check("qr: Qᵀ Q = I for a random " + shape + " matrix", orth, 1e-14);
    double recon = 0;
    for (std::size_t i = 0; i < m; i++)
      for (std::size_t j = 0; j < n; j++) {
        double v = 0;
        for (std::size_t k = 0; k <= std::min(j, t - 1); k++) v += Q(i, k) * f.factors(k, j);
        recon = std::max(recon, std::fabs(v - A(i, f.permutation[j])));
      }
    s.check("qr: A P = Q R for the " + shape + " matrix", recon, 1e-13);
    bool decreasing = true;
    for (std::size_t k = 1; k < t; k++) decreasing = decreasing && std::fabs(f.factors(k, k)) <= std::fabs(f.factors(k - 1, k - 1)) * (1 + 1e-12);
    s.require("qr: column pivoting makes |R_kk| non-increasing (" + shape + ")", decreasing);
    s.require("qr: the rank of a random " + shape + " matrix is min(m, n)", f.rank == t, std::to_string(f.rank));
  }

  // ---- least squares: a line in closed form ----
  {
    const std::vector<double> x{0, 1, 2, 3, 4, 5}, y{1.0, 3.1, 4.9, 7.2, 8.8, 11.1};
    Rectangular X(6, 2);
    double sx = 0, sy = 0;
    for (std::size_t i = 0; i < 6; i++) {
      X(i, 0) = 1;
      X(i, 1) = x[i];
      sx += x[i];
      sy += y[i];
    }
    const double mx = sx / 6, my = sy / 6;
    double sxx = 0, sxy = 0;
    for (std::size_t i = 0; i < 6; i++) {
      sxx += (x[i] - mx) * (x[i] - mx);
      sxy += (x[i] - mx) * (y[i] - my);
    }
    const double b = sxy / sxx, a = my - b * mx;
    const auto beta = least_squares(X, y);
    s.check("least squares: the slope of a fitted line, against Sxy / Sxx", beta[1] - b, 1e-14);
    s.check("least squares: its intercept, against ȳ − b x̄", beta[0] - a, 1e-14);
  }

  // ---- against the normal equations, and the residual's orthogonality ----
  {
    const Rectangular A = random_matrix(40, 6, 4);
    calx::random::Stream rng(20261024, 8, 0, 0);
    std::vector<double> b(40);
    for (double& v : b) v = rng.normal();
    const auto x = least_squares(A, b);
    const Matrix G = gram(A);
    const auto atb = multiply_transpose(A, b);
    Matrix B(6);
    for (std::size_t i = 0; i < 6; i++) B(i, 0) = atb[i];
    const Matrix xn = solve(G, B);
    double diff = 0;
    for (std::size_t i = 0; i < 6; i++) diff = std::max(diff, std::fabs(x[i] - xn(i, 0)));
    s.check("least squares: QR against the normal equations on a well-conditioned 40 × 6 problem", diff, 1e-13);
    auto r = multiply(A, x);
    for (std::size_t i = 0; i < 40; i++) r[i] = b[i] - r[i];
    s.check("least squares: the residual is orthogonal to the columns, ‖Aᵀ r‖∞", max_abs(multiply_transpose(A, r)), 1e-13);
  }

  // ---- Läuchli's matrix: AᵀA is singular in floating point, A is not ----
  {
    const double d = 1e-8;  // 1 + d² rounds to 1
    const Rectangular A = Rectangular::from_rows({{1, 1}, {d, 0}, {0, d}});
    const std::vector<double> b{2, d, d};  // A (1, 1) = b exactly
    const auto x = least_squares(A, b);
    s.check("least squares: Läuchli's matrix (δ = 1e-8), x₁ = 1, by QR", x[0] - 1, 1e-7);
    s.check("least squares: Läuchli's matrix, x₂ = 1", x[1] - 1, 1e-7);
    const Matrix G = gram(A);
    s.require("least squares: its Gram matrix is singular in floating point (1 + δ² = 1)", G(0, 0) * G(1, 1) - G(0, 1) * G(1, 0) == 0);
  }

  // ---- a square system, against Gaussian elimination ----
  {
    const Rectangular A = random_matrix(8, 8, 5);
    std::vector<double> b(8);
    for (std::size_t i = 0; i < 8; i++) b[i] = static_cast<double>(i) - 3.5;
    const auto x = least_squares(A, b);
    Matrix B(8);
    for (std::size_t i = 0; i < 8; i++) B(i, 0) = b[i];
    const Matrix xs = solve(to_square(A), B);
    double diff = 0;
    for (std::size_t i = 0; i < 8; i++) diff = std::max(diff, std::fabs(x[i] - xs(i, 0)));
    s.check("least squares: a square system, against Gaussian elimination", diff, 1e-12);
  }

  // ---- the numerical rank ----
  {
    Rectangular A = random_matrix(20, 4, 6);
    for (std::size_t i = 0; i < 20; i++) A(i, 3) = 2 * A(i, 0) - A(i, 2);  // dependent
    const Qr f = qr(A);
    s.require("qr: a column that is a combination of two others leaves rank 3", f.rank == 3, std::to_string(f.rank));
    const Qr g = qr(A, true, 0.0);  // tolerance 0: rounding counts as rank
    s.require("qr: with tolerance 0 the rounding error is taken for rank", g.rank == 4, std::to_string(g.rank));
    // The basic solution fits a consistent right-hand side exactly.
    std::vector<double> coef{1, -2, 0.5, 0}, b = multiply(A, coef);
    const auto x = least_squares(f, b);
    auto r = multiply(A, x);
    for (std::size_t i = 0; i < 20; i++) r[i] -= b[i];
    s.check("least squares: a rank-deficient but consistent system is fitted exactly", max_abs(r), 1e-12);
  }

  // ---- (AᵀA)⁻¹ ----
  {
    const Rectangular A = random_matrix(25, 5, 7);
    const Matrix C = normal_inverse(qr(A)), Gi = inverse(gram(A));
    double diff = 0;
    for (std::size_t i = 0; i < 25; i++) diff = std::max(diff, std::fabs(C.a[i] - Gi.a[i]));
    s.check("qr: (AᵀA)⁻¹ = R⁻¹R⁻ᵀ in the original order, against the inverse of the Gram matrix", diff, 1e-13);
  }

  // ---- Rectangular ----
  {
    const Rectangular A = Rectangular::from_rows({{1, 2, 3}, {4, 5, 6}});
    const Rectangular B = Rectangular::from_columns({{1, 0, -1}, {2, 1, 0}});  // 3 × 2
    const Rectangular C = A * B;  // [[−2, 4], [−2, 13]]
    s.require("rectangular: a 2 × 3 by 3 × 2 product", C.rows == 2 && C.cols == 2 && C(0, 0) == -2 && C(0, 1) == 4 && C(1, 0) == -2 && C(1, 1) == 13);
    const Rectangular At = A.transpose();
    s.require("rectangular: the transpose", At.rows == 3 && At.cols == 2 && At(2, 1) == 6 && At(0, 1) == 4);
    s.require("rectangular: rows and columns", A.row(1) == std::vector<double>{4, 5, 6} && A.column(2) == std::vector<double>{3, 6});
  }
  return s.done();
}
