// Checks of the symmetric eigenproblem: the eigenvalues of the second-difference matrix,
// 2 − 2 cos(kπ/(n + 1)), and its eigenvectors, the sines; A V = V Λ and Vᵀ V = I on random
// symmetric matrices of orders 1 to 120; repeated eigenvalues (the identity, a rank-one
// matrix); a diagonal matrix returns its diagonal sorted; and the exponential of a reversible
// rate matrix by its spectral decomposition against calx's expm.

#include <calx/linalg/eigen.hpp>
#include <calx/linalg/expm.hpp>
#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>
#include <numbers>
#include <string>
#include <vector>

using namespace calx::linalg;

namespace {

// max |A V − V Λ| / ‖A‖ and max |Vᵀ V − I|.
std::pair<double, double> residuals(const Matrix& A, const SymmetricEigen& E) {
  const std::size_t n = A.n;
  double scale = 0;
  for (double v : A.a) scale = std::max(scale, std::fabs(v));
  scale = std::max(scale, 1e-300);
  double r = 0, o = 0;
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t k = 0; k < n; k++) {
      double av = 0, vtv = 0;
      for (std::size_t j = 0; j < n; j++) {
        av += A(i, j) * E.vectors(j, k);
        vtv += E.vectors(j, i) * E.vectors(j, k);
      }
      r = std::max(r, std::fabs(av - E.vectors(i, k) * E.values[k]) / scale);
      o = std::max(o, std::fabs(vtv - (i == k ? 1 : 0)));
    }
  return {r, o};
}

}  // namespace

int main() {
  calx::verification::Suite s;

  // The second-difference matrix tridiag(−1, 2, −1) of order 50.
  {
    const std::size_t n = 50;
    Matrix A(n);
    for (std::size_t i = 0; i < n; i++) {
      A(i, i) = 2;
      if (i + 1 < n) A(i, i + 1) = A(i + 1, i) = -1;
    }
    const auto E = symmetric_eigen(A);
    double err = 0, vec = 0;
    for (std::size_t k = 0; k < n; k++) {
      err = std::max(err, std::fabs(E.values[k] - (2 - 2 * std::cos((k + 1) * std::numbers::pi / (n + 1)))));
      // The eigenvector of the k-th value is sin(j (k + 1) π / (n + 1)), normalised, up to sign.
      double dot = 0;
      for (std::size_t j = 0; j < n; j++) dot += E.vectors(j, k) * std::sqrt(2.0 / (n + 1)) * std::sin((j + 1) * (k + 1) * std::numbers::pi / (n + 1));
      vec = std::max(vec, 1 - std::fabs(dot));
    }
    s.check("second differences: the eigenvalues are 2 − 2 cos(kπ/51)", err, 1e-13);
    s.check("second differences: the eigenvectors are the sines", vec, 1e-12);
  }

  // Random symmetric matrices.
  calx::random::Stream r(20261024, 0x65696721, 0, 0);
  double worst_r = 0, worst_o = 0;
  bool sorted = true;
  for (std::size_t n : {1, 2, 3, 5, 8, 17, 40, 120}) {
    Matrix A(n);
    for (std::size_t i = 0; i < n; i++)
      for (std::size_t j = 0; j <= i; j++) A(i, j) = A(j, i) = r.normal();
    const auto E = symmetric_eigen(A);
    const auto [res, orth] = residuals(A, E);
    worst_r = std::max(worst_r, res);
    worst_o = std::max(worst_o, orth);
    sorted = sorted && std::is_sorted(E.values.begin(), E.values.end());
  }
  s.check("random symmetric matrices of orders 1 to 120: max |AV − VΛ| / max |A|", worst_r, 1e-12);
  s.check("random symmetric matrices: max |VᵀV − I|", worst_o, 1e-12);
  s.require("the eigenvalues are ascending", sorted);

  // Repeated eigenvalues: the identity, and the rank-one matrix of ones (n, 0, ..., 0).
  {
    const auto I = symmetric_eigen(Matrix::identity(6));
    const auto [ri, oi] = residuals(Matrix::identity(6), I);
    s.check("the identity: its eigenpairs", std::max(ri, oi), 1e-15);
    Matrix J(7);
    for (double& v : J.a) v = 1;
    const auto E = symmetric_eigen(J);
    double err = std::fabs(E.values.back() - 7);
    for (std::size_t k = 0; k + 1 < 7; k++) err = std::max(err, std::fabs(E.values[k]));
    const auto [rj, oj] = residuals(J, E);
    s.check("the matrix of ones: eigenvalues 0 (six times) and 7", err, 1e-13);
    s.check("the matrix of ones: an orthonormal basis of eigenvectors", std::max(rj, oj), 1e-13);
  }

  // A diagonal matrix.
  {
    Matrix D(4);
    D(0, 0) = 3; D(1, 1) = -1; D(2, 2) = 2; D(3, 3) = 0.5;
    const auto E = symmetric_eigen(D);
    s.require("a diagonal matrix: its diagonal, sorted", E.values == std::vector<double>{-1, 0.5, 2, 3});
  }

  // exp(Qt) for the reversible HKY85 rate matrix (κ = 4, π = 0.1, 0.2, 0.3, 0.4), by
  // Q = Π^{-1/2} S Π^{1/2} with S symmetric, against the scaling-and-squaring expm.
  {
    const double pi[4] = {0.1, 0.2, 0.3, 0.4}, kappa = 4;
    Matrix Q(4);
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++)
        if (i != j) Q(i, j) = pi[j] * (((i ^ j) == 2) ? kappa : 1);  // A=0 C=1 G=2 T=3: A↔G and C↔T are transitions
    for (int i = 0; i < 4; i++) {
      double row = 0;
      for (int j = 0; j < 4; j++)
        if (j != i) row += Q(i, j);
      Q(i, i) = -row;
    }
    Matrix S(4);
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++) S(i, j) = std::sqrt(pi[i]) * Q(i, j) / std::sqrt(pi[j]);
    const auto E = symmetric_eigen(S);
    const double t = 0.37;
    Matrix P(4);
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++) {
        double v = 0;
        for (int k = 0; k < 4; k++) v += E.vectors(i, k) * std::exp(E.values[k] * t) * E.vectors(j, k);
        P(i, j) = v * std::sqrt(pi[j]) / std::sqrt(pi[i]);
      }
    const Matrix R = expm(t * Q);
    double err = 0;
    for (std::size_t k = 0; k < 16; k++) err = std::max(err, std::fabs(P.a[k] - R.a[k]));
    s.check("exp(Qt) of HKY85 by its spectral decomposition, against expm", err, 1e-13);
    s.check("the largest eigenvalue of a rate matrix is 0", E.values.back(), 1e-14);
  }

  return s.done();
}
