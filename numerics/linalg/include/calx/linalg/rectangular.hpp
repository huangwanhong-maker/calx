// Dense rectangular matrices in double precision.
//
// A Rectangular is an m × n matrix of run-time shape, row-major: the design matrices of
// regressions (one row per observation, one column per regressor), the constraint matrices of
// linear programmes, the data of several series side by side. `Matrix` (dense.hpp) stays the
// square matrix of small order on which solve, lu, expm and the eigenproblems work; the two meet
// in gram() (AᵀA, square) and in from_square() / to_square().
//
// Products are the plain triple loops, O(mnk), adequate for the moderate sizes of statistics
// (thousands of rows, tens of columns); large dense problems will go to blocked kernels or a
// borrowed BLAS behind this interface (docs/mathematics.md of the infrastructure, §10.3).

#pragma once

#include <calx/linalg/dense.hpp>

#include <cstddef>
#include <span>
#include <stdexcept>
#include <vector>

namespace calx::linalg {

struct Rectangular {
  std::size_t rows = 0, cols = 0;
  std::vector<double> a;  // row-major: entry (i, j) at a[i * cols + j]

  Rectangular() = default;
  Rectangular(std::size_t r, std::size_t c) : rows(r), cols(c), a(r * c, 0.0) {}

  double& operator()(std::size_t i, std::size_t j) { return a[i * cols + j]; }
  double operator()(std::size_t i, std::size_t j) const { return a[i * cols + j]; }

  // The matrix whose columns are the given vectors, all of one length.
  static Rectangular from_columns(const std::vector<std::vector<double>>& columns) {
    Rectangular m(columns.empty() ? 0 : columns[0].size(), columns.size());
    for (std::size_t j = 0; j < m.cols; j++) {
      if (columns[j].size() != m.rows) throw std::invalid_argument("calx::linalg::Rectangular::from_columns: columns of one length");
      for (std::size_t i = 0; i < m.rows; i++) m(i, j) = columns[j][i];
    }
    return m;
  }
  // The matrix whose rows are the given vectors, all of one length.
  static Rectangular from_rows(const std::vector<std::vector<double>>& rows) {
    Rectangular m(rows.size(), rows.empty() ? 0 : rows[0].size());
    for (std::size_t i = 0; i < m.rows; i++) {
      if (rows[i].size() != m.cols) throw std::invalid_argument("calx::linalg::Rectangular::from_rows: rows of one length");
      for (std::size_t j = 0; j < m.cols; j++) m(i, j) = rows[i][j];
    }
    return m;
  }
  static Rectangular from_square(const Matrix& s) {
    Rectangular m(s.n, s.n);
    m.a = s.a;
    return m;
  }

  std::vector<double> column(std::size_t j) const {
    std::vector<double> c(rows);
    for (std::size_t i = 0; i < rows; i++) c[i] = (*this)(i, j);
    return c;
  }
  std::vector<double> row(std::size_t i) const { return std::vector<double>(a.begin() + static_cast<std::ptrdiff_t>(i * cols), a.begin() + static_cast<std::ptrdiff_t>((i + 1) * cols)); }

  Rectangular transpose() const {
    Rectangular t(cols, rows);
    for (std::size_t i = 0; i < rows; i++)
      for (std::size_t j = 0; j < cols; j++) t(j, i) = (*this)(i, j);
    return t;
  }
};

inline Matrix to_square(const Rectangular& r) {
  if (r.rows != r.cols) throw std::invalid_argument("calx::linalg::to_square: a square matrix");
  Matrix m(r.rows);
  m.a = r.a;
  return m;
}

inline Rectangular operator*(const Rectangular& x, const Rectangular& y) {
  if (x.cols != y.rows) throw std::invalid_argument("calx::linalg: the shapes of a product do not compose");
  Rectangular z(x.rows, y.cols);
  for (std::size_t i = 0; i < x.rows; i++)
    for (std::size_t k = 0; k < x.cols; k++) {
      const double v = x(i, k);
      if (v == 0) continue;
      for (std::size_t j = 0; j < y.cols; j++) z(i, j) += v * y(k, j);
    }
  return z;
}

// A x.
inline std::vector<double> multiply(const Rectangular& A, std::span<const double> x) {
  if (x.size() != A.cols) throw std::invalid_argument("calx::linalg::multiply: x has one entry per column");
  std::vector<double> y(A.rows, 0.0);
  for (std::size_t i = 0; i < A.rows; i++) {
    double s = 0;
    for (std::size_t j = 0; j < A.cols; j++) s += A(i, j) * x[j];
    y[i] = s;
  }
  return y;
}

// Aᵀ x.
inline std::vector<double> multiply_transpose(const Rectangular& A, std::span<const double> x) {
  if (x.size() != A.rows) throw std::invalid_argument("calx::linalg::multiply_transpose: x has one entry per row");
  std::vector<double> y(A.cols, 0.0);
  for (std::size_t i = 0; i < A.rows; i++) {
    const double v = x[i];
    if (v == 0) continue;
    for (std::size_t j = 0; j < A.cols; j++) y[j] += A(i, j) * v;
  }
  return y;
}

// A x for a square Matrix.
inline std::vector<double> multiply(const Matrix& A, std::span<const double> x) {
  if (x.size() != A.n) throw std::invalid_argument("calx::linalg::multiply: x has one entry per column");
  std::vector<double> y(A.n, 0.0);
  for (std::size_t i = 0; i < A.n; i++) {
    double s = 0;
    for (std::size_t j = 0; j < A.n; j++) s += A(i, j) * x[j];
    y[i] = s;
  }
  return y;
}

// The Gram matrix AᵀA, square of order A.cols. Forming it squares the condition number; least
// squares go through qr.hpp instead, and gram() serves where the cross-products themselves are
// wanted (covariances of data, the meat of a sandwich).
inline Matrix gram(const Rectangular& A) {
  Matrix g(A.cols);
  for (std::size_t i = 0; i < A.rows; i++)
    for (std::size_t j = 0; j < A.cols; j++) {
      const double v = A(i, j);
      if (v == 0) continue;
      for (std::size_t k = j; k < A.cols; k++) g(j, k) += v * A(i, k);
    }
  for (std::size_t j = 0; j < A.cols; j++)
    for (std::size_t k = 0; k < j; k++) g(j, k) = g(k, j);
  return g;
}

// The transpose of a square Matrix.
inline Matrix transpose(const Matrix& m) {
  Matrix t(m.n);
  for (std::size_t i = 0; i < m.n; i++)
    for (std::size_t j = 0; j < m.n; j++) t(j, i) = m(i, j);
  return t;
}

// A⁻¹ by Gaussian elimination with partial pivoting (solve against the identity).
inline Matrix inverse(const Matrix& A) { return solve(A, Matrix::identity(A.n)); }

}  // namespace calx::linalg
