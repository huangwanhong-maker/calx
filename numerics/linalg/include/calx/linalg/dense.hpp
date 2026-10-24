// Dense matrices of small order, in double precision.
//
// A Matrix is square, of run-time order, row-major. It serves the small systems of models
// (a propagator of order below twenty, computed once and applied many times); large dense
// problems will go to blocked kernels or a borrowed BLAS behind this module's interface
// (docs/mathematics.md of the infrastructure, §5.3).
//
// solve() is Gaussian elimination with partial pivoting (Golub and Van Loan 2013,
// Matrix Computations, 4th ed., Algorithm 3.4.1), backward stable in practice; O(n^3).
// determinant() is the product of the pivots of the same elimination, with the sign of its row
// exchanges (Golub and Van Loan 2013, §3.4); log_determinant() its logarithm, which neither
// overflows nor underflows where the product would.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

namespace calx::linalg {

struct Matrix {
  std::size_t n = 0;
  std::vector<double> a;

  Matrix() = default;
  explicit Matrix(std::size_t order) : n(order), a(order * order, 0.0) {}
  static Matrix identity(std::size_t order) {
    Matrix m(order);
    for (std::size_t i = 0; i < order; i++) m(i, i) = 1;
    return m;
  }
  double& operator()(std::size_t i, std::size_t j) { return a[i * n + j]; }
  double operator()(std::size_t i, std::size_t j) const { return a[i * n + j]; }
};

inline Matrix operator*(const Matrix& x, const Matrix& y) {
  Matrix z(x.n);
  for (std::size_t i = 0; i < x.n; i++)
    for (std::size_t k = 0; k < x.n; k++) {
      const double v = x(i, k);
      if (v == 0) continue;
      for (std::size_t j = 0; j < x.n; j++) z(i, j) += v * y(k, j);
    }
  return z;
}
inline Matrix operator+(Matrix x, const Matrix& y) {
  for (std::size_t i = 0; i < x.a.size(); i++) x.a[i] += y.a[i];
  return x;
}
inline Matrix operator-(Matrix x, const Matrix& y) {
  for (std::size_t i = 0; i < x.a.size(); i++) x.a[i] -= y.a[i];
  return x;
}
inline Matrix operator*(double s, Matrix x) {
  for (double& v : x.a) v *= s;
  return x;
}

// The 1-norm: the largest column sum of absolute values.
inline double norm1(const Matrix& m) {
  double best = 0;
  for (std::size_t j = 0; j < m.n; j++) {
    double s = 0;
    for (std::size_t i = 0; i < m.n; i++) s += std::fabs(m(i, j));
    best = std::max(best, s);
  }
  return best;
}

// X = A^-1 B.
inline Matrix solve(Matrix A, Matrix B) {
  const std::size_t n = A.n;
  for (std::size_t c = 0; c < n; c++) {
    std::size_t p = c;
    for (std::size_t r = c + 1; r < n; r++)
      if (std::fabs(A(r, c)) > std::fabs(A(p, c))) p = r;
    if (A(p, c) == 0) throw std::runtime_error("calx::linalg::solve: singular matrix");
    if (p != c)
      for (std::size_t j = 0; j < n; j++) {
        std::swap(A(c, j), A(p, j));
        std::swap(B(c, j), B(p, j));
      }
    for (std::size_t r = c + 1; r < n; r++) {
      const double f = A(r, c) / A(c, c);
      if (f == 0) continue;
      for (std::size_t j = c; j < n; j++) A(r, j) -= f * A(c, j);
      for (std::size_t j = 0; j < n; j++) B(r, j) -= f * B(c, j);
    }
  }
  Matrix X(n);
  for (std::size_t j = 0; j < n; j++)
    for (std::size_t r = n; r-- > 0;) {
      double s = B(r, j);
      for (std::size_t k = r + 1; k < n; k++) s -= A(r, k) * X(k, j);
      X(r, j) = s / A(r, r);
    }
  return X;
}

// The LU factors of A with partial pivoting, in place: U on and above the diagonal, the
// multipliers of L below it; the sign of the row permutation. A zero pivot leaves the rest unfactored.
struct Lu {
  Matrix lu;
  int sign = 1;
  bool singular = false;
};
inline Lu lu(Matrix A) {
  const std::size_t n = A.n;
  Lu f;
  for (std::size_t c = 0; c < n; c++) {
    std::size_t p = c;
    for (std::size_t r = c + 1; r < n; r++)
      if (std::fabs(A(r, c)) > std::fabs(A(p, c))) p = r;
    if (A(p, c) == 0) {
      f.singular = true;
      break;
    }
    if (p != c) {
      for (std::size_t j = 0; j < n; j++) std::swap(A(c, j), A(p, j));
      f.sign = -f.sign;
    }
    for (std::size_t r = c + 1; r < n; r++) {
      const double m = A(r, c) / A(c, c);
      A(r, c) = m;
      if (m == 0) continue;
      for (std::size_t j = c + 1; j < n; j++) A(r, j) -= m * A(c, j);
    }
  }
  f.lu = std::move(A);
  return f;
}

inline double determinant(const Matrix& A) {
  const Lu f = lu(A);
  if (f.singular) return 0;
  double d = f.sign;
  for (std::size_t i = 0; i < A.n; i++) d *= f.lu(i, i);
  return d;
}

// log |det A| and the sign of det A (0 when singular, with log −∞).
struct LogDeterminant {
  double log_abs = 0;
  int sign = 1;
};
inline LogDeterminant log_determinant(const Matrix& A) {
  const Lu f = lu(A);
  if (f.singular) return {-HUGE_VAL, 0};
  LogDeterminant r{0, f.sign};
  for (std::size_t i = 0; i < A.n; i++) {
    const double u = f.lu(i, i);
    r.log_abs += std::log(std::fabs(u));
    if (u < 0) r.sign = -r.sign;
  }
  return r;
}

}  // namespace calx::linalg
