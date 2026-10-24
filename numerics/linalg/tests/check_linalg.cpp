// Checks of calx/linalg: the solve of a dense system by its residual; the exponential
// against closed forms (diagonal, nilpotent, a rotation with damping at small and large
// norm, so that scaling and squaring is exercised) and against exp(A) exp(-A) = I; and
// compressed sparse rows against the dense product; determinants against triangular matrices,
// permutations, multiplicativity, Vandermonde's product and the Hilbert matrix.

#include <calx/linalg/dense.hpp>
#include <calx/linalg/expm.hpp>
#include <calx/linalg/sparse.hpp>
#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <cstdint>
#include <vector>

using calx::linalg::Matrix;

namespace {

double max_abs(const Matrix& m) {
  double b = 0;
  for (double v : m.a) b = std::max(b, std::fabs(v));
  return b;
}

// [[b, -w], [w, b]] t has exponential e^(b t) [[cos w t, -sin w t], [sin w t, cos w t]].
double rotation_error(double b, double w, double t) {
  Matrix A(2);
  A(0, 0) = b * t; A(0, 1) = -w * t; A(1, 0) = w * t; A(1, 1) = b * t;
  const Matrix E = calx::linalg::expm(A);
  const double g = std::exp(b * t), c = std::cos(w * t), s = std::sin(w * t);
  const double ref[4] = {g * c, -g * s, g * s, g * c};
  double e = 0;
  for (int i = 0; i < 4; i++) e = std::max(e, std::fabs(E.a[i] - ref[i]) / g);
  return e;
}

}  // namespace

int main() {
  calx::verification::Suite s;
  calx::random::Stream r(1, 0, 0, 0);

  // ---- solve ------------------------------------------------------------------------------------
  {
    const std::size_t n = 7;
    Matrix A(n), B(n);
    for (double& v : A.a) v = r.normal();
    for (double& v : B.a) v = r.normal();
    const Matrix X = calx::linalg::solve(A, B);
    s.check("linalg: a dense 7 x 7 system solved, the largest residual of A X - B", max_abs(A * X - B), 1e-12);
  }
  // ---- the exponential ---------------------------------------------------------------------------
  {
    Matrix D(3);
    D(0, 0) = -2; D(1, 1) = 0.5; D(2, 2) = 3;
    const Matrix E = calx::linalg::expm(D);
    const double e = std::max({std::fabs(E(0, 0) - std::exp(-2)) / std::exp(-2), std::fabs(E(1, 1) - std::exp(0.5)) / std::exp(0.5),
                               std::fabs(E(2, 2) - std::exp(3)) / std::exp(3), std::fabs(E(0, 1)), std::fabs(E(2, 0))});
    s.check("linalg: exp of a diagonal matrix, relative error", e, 1e-14);
    Matrix N(3);
    N(0, 1) = 1; N(1, 2) = 1;  // exp(N) = I + N + N^2 / 2
    const Matrix F = calx::linalg::expm(N);
    s.check("linalg: exp of a nilpotent matrix is its finite series", std::max({std::fabs(F(0, 0) - 1), std::fabs(F(0, 1) - 1), std::fabs(F(0, 2) - 0.5), std::fabs(F(1, 0)), std::fabs(F(2, 2) - 1)}), 1e-15);
    s.check("linalg: exp of a damped rotation, small norm (no squaring)", rotation_error(-0.3, 1.1, 1), 1e-14);
    s.check("linalg: exp of a damped rotation, norm 70 (seven squarings)", rotation_error(-0.5, 50, 1), 1e-11);
    Matrix A(5);
    for (double& v : A.a) v = r.normal();
    const Matrix P = calx::linalg::expm(A) * calx::linalg::expm(-1.0 * A);
    s.check("linalg: exp(A) exp(-A) = I for a random 5 x 5 matrix", max_abs(P - Matrix::identity(5)), 1e-12);
  }
  // ---- compressed sparse rows -------------------------------------------------------------------------
  {
    const std::size_t rows = 50, cols = 40;
    std::vector<std::uint32_t> ri, ci;
    std::vector<double> vi;
    Matrix dense(std::max(rows, cols));
    for (int k = 0; k < 600; k++) {
      const auto i = static_cast<std::uint32_t>(r.below(rows)), j = static_cast<std::uint32_t>(r.below(cols));
      const double v = r.normal();
      ri.push_back(i); ci.push_back(j); vi.push_back(v);
      dense(i, j) += v;  // repeated coordinates add, as in the product
    }
    std::vector<std::uint64_t> order;
    const auto A = calx::linalg::compress<std::uint32_t>(rows, cols, ri, ci, &order);
    std::vector<double> values(order.size());
    for (std::size_t k = 0; k < order.size(); k++) values[k] = vi[order[k]];
    bool stable = true;
    for (std::size_t i = 0; i < rows; i++)
      for (auto k = A.begin(i); k + 1 < A.end(i); k++) stable &= order[k] < order[k + 1];
    s.require("linalg: compressed rows keep the given order within each row", stable && A.entries() == 600 && A.row.back() == 600);
    std::vector<double> x(cols), y(rows);
    for (double& v : x) v = r.normal();
    calx::linalg::multiply<std::uint32_t, double>(A, values, x, y);
    double err = 0;
    for (std::size_t i = 0; i < rows; i++) {
      double d = 0;
      for (std::size_t j = 0; j < cols; j++) d += dense(i, j) * x[j];
      err = std::max(err, std::fabs(d - y[i]));
    }
    s.check("linalg: the product by compressed rows equals the dense product", err, 1e-12);
  }
  {
    // Determinants: a triangular matrix's is the product of its diagonal; a permutation's is its
    // sign; det(AB) = det A det B; Vandermonde's product; the Hilbert matrix of order 6,
    // 1 / 186313420339200000 (det Hₙ = cₙ⁴ / c₂ₙ, cₙ = ∏_{i<n} i!); a singular matrix.
    calx::linalg::Matrix T(4);
    double diag = 1;
    for (std::size_t i = 0; i < 4; i++)
      for (std::size_t j = i; j < 4; j++) T(i, j) = (i == j) ? (diag *= -1.5 - static_cast<double>(i), -1.5 - static_cast<double>(i)) : 0.3 * static_cast<double>(i + j);
    s.check("linalg: a triangular matrix's determinant is its diagonal's product", calx::linalg::determinant(T) - diag, 1e-12);
    calx::linalg::Matrix P(3);
    P(0, 1) = P(1, 2) = P(2, 0) = 1;  // a 3-cycle: even
    calx::linalg::Matrix Q(3);
    Q(0, 1) = Q(1, 0) = Q(2, 2) = 1;  // a transposition: odd
    s.require("linalg: a permutation matrix's determinant is its sign", calx::linalg::determinant(P) == 1 && calx::linalg::determinant(Q) == -1);
    calx::random::Stream r(20261024, 0x64657420, 0, 0);
    double worst = 0;
    for (int rep = 0; rep < 50; rep++) {
      calx::linalg::Matrix A(5), B(5);
      for (auto& v : A.a) v = r.uniform() * 2 - 1;
      for (auto& v : B.a) v = r.uniform() * 2 - 1;
      const double da = calx::linalg::determinant(A), db = calx::linalg::determinant(B), dab = calx::linalg::determinant(A * B);
      worst = std::max(worst, std::fabs(dab - da * db) / std::max(1e-300, std::fabs(da * db)));
    }
    s.check("linalg: det(AB) = det A det B on 50 random pairs of order 5 (relative)", worst, 1e-11);
    const double xs[5] = {0.5, -1.25, 2, 3.5, -0.75};
    calx::linalg::Matrix V(5);
    double vandermonde = 1;
    for (std::size_t i = 0; i < 5; i++) {
      for (std::size_t j = 0; j < 5; j++) V(i, j) = std::pow(xs[i], static_cast<double>(j));
      for (std::size_t j = i + 1; j < 5; j++) vandermonde *= xs[j] - xs[i];
    }
    s.check("linalg: Vandermonde's determinant is ∏_{i<j} (x_j − x_i)", calx::linalg::determinant(V) / vandermonde - 1, 1e-12);
    calx::linalg::Matrix H(6);
    for (std::size_t i = 0; i < 6; i++)
      for (std::size_t j = 0; j < 6; j++) H(i, j) = 1.0 / static_cast<double>(i + j + 1);
    const auto lh = calx::linalg::log_determinant(H);
    s.check("linalg: the Hilbert matrix of order 6, through its logarithm (ill-conditioned, κ ≈ 1.5e7)", lh.log_abs + std::log(186313420339200000.0), 1e-8);
    s.require("linalg: and positive", lh.sign == 1);
    calx::linalg::Matrix S(3);
    S(0, 0) = 1; S(0, 1) = 2; S(0, 2) = 3; S(1, 0) = 2; S(1, 1) = 4; S(1, 2) = 6; S(2, 0) = 1; S(2, 1) = 0; S(2, 2) = 1;
    s.require("linalg: a singular matrix's determinant is 0, its log −∞ with sign 0", calx::linalg::determinant(S) == 0 && calx::linalg::log_determinant(S).sign == 0);
  }
  return s.done();
}
