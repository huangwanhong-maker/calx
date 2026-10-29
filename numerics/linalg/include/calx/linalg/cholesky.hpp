// The Cholesky factorisation A = L Lᵀ of a symmetric positive definite matrix, and the
// factorisation of symmetric positive definite band matrices.
//
// cholesky(): the gaxpy Cholesky of Golub and Van Loan (2013), Matrix Computations, 4th ed.,
// Algorithm 4.2.1, reading the lower triangle of A. It is backward stable without pivoting
// (Higham 2002, Theorem 10.3): the computed L is exact for A + E with |E| ≤ c n ε |L||Lᵀ|. A
// non-positive pivot ends it with positive_definite = false (and L filled up to there); that
// is the cheapest test of definiteness there is. O(n³/3).
//
// Uses: the orthogonalised impulse responses of a vector autoregression (Σ = P Pᵀ), sampling
// from a multivariate normal, solving normal equations that are known to be well conditioned.
//
// Band matrices. A symmetric positive definite matrix with half-bandwidth p (a_ij = 0 for
// |i − j| > p) keeps a Cholesky factor of the same band (Golub and Van Loan, §4.3,
// band Cholesky and the band triangular solves): O(n p²) flops and O(n p) memory, so that the
// pentadiagonal system of the Hodrick–Prescott filter, (I + λ DᵀD) τ = y with p = 2, is solved in
// linear time for any length. Accuracy as for the dense factorisation.

#pragma once

#include <calx/linalg/dense.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <vector>

namespace calx::linalg {

struct Cholesky {
  Matrix l;  // lower triangular, A = L Lᵀ
  bool positive_definite = false;
};

inline Cholesky cholesky(const Matrix& A) {
  const std::size_t n = A.n;
  Cholesky f;
  f.l = Matrix(n);
  Matrix& L = f.l;
  for (std::size_t j = 0; j < n; j++) {
    double d = A(j, j);
    for (std::size_t k = 0; k < j; k++) d -= L(j, k) * L(j, k);
    if (!(d > 0)) return f;
    const double ljj = std::sqrt(d);
    L(j, j) = ljj;
    for (std::size_t i = j + 1; i < n; i++) {
      double s = A(i, j);
      for (std::size_t k = 0; k < j; k++) s -= L(i, k) * L(j, k);
      L(i, j) = s / ljj;
    }
  }
  f.positive_definite = true;
  return f;
}

// x with A x = b from A's Cholesky factor: L y = b, then Lᵀ x = y.
inline std::vector<double> cholesky_solve(const Cholesky& f, std::span<const double> b) {
  if (!f.positive_definite) throw std::runtime_error("calx::linalg::cholesky_solve: the matrix is not positive definite");
  const std::size_t n = f.l.n;
  if (b.size() != n) throw std::invalid_argument("calx::linalg::cholesky_solve: b has one entry per row");
  std::vector<double> x(b.begin(), b.end());
  for (std::size_t i = 0; i < n; i++) {
    for (std::size_t k = 0; k < i; k++) x[i] -= f.l(i, k) * x[k];
    x[i] /= f.l(i, i);
  }
  for (std::size_t i = n; i-- > 0;) {
    for (std::size_t k = i + 1; k < n; k++) x[i] -= f.l(k, i) * x[k];
    x[i] /= f.l(i, i);
  }
  return x;
}

// A symmetric band matrix of order n and half-bandwidth p, its lower band stored by rows:
// entry (i, j), i − p ≤ j ≤ i, at a[i (p + 1) + (j − i + p)].
struct SymmetricBand {
  std::size_t n = 0, p = 0;
  std::vector<double> a;

  SymmetricBand() = default;
  SymmetricBand(std::size_t order, std::size_t half_bandwidth) : n(order), p(half_bandwidth), a(order * (half_bandwidth + 1), 0.0) {}
  // Entry (i, j) of the lower band, i − p ≤ j ≤ i.
  double& operator()(std::size_t i, std::size_t j) { return a[i * (p + 1) + (j + p - i)]; }
  double operator()(std::size_t i, std::size_t j) const { return a[i * (p + 1) + (j + p - i)]; }
  // Any entry, zero outside the band.
  double at(std::size_t i, std::size_t j) const {
    if (j > i) std::swap(i, j);
    return i - j > p ? 0.0 : (*this)(i, j);
  }
};

// The band Cholesky factor (Golub and Van Loan, §4.3), in the same storage: L's band.
struct BandCholesky {
  SymmetricBand l;
  bool positive_definite = false;
};

inline BandCholesky band_cholesky(const SymmetricBand& A) {
  const std::size_t n = A.n, p = A.p;
  BandCholesky f;
  f.l = A;
  SymmetricBand& L = f.l;
  for (std::size_t j = 0; j < n; j++) {
    const std::size_t k0 = j > p ? j - p : 0;
    double d = L(j, j);
    for (std::size_t k = k0; k < j; k++) d -= L(j, k) * L(j, k);
    if (!(d > 0)) return f;
    const double ljj = std::sqrt(d);
    L(j, j) = ljj;
    for (std::size_t i = j + 1; i <= std::min(n - 1, j + p); i++) {
      double s = L(i, j);
      const std::size_t kk = i > p ? i - p : 0;
      for (std::size_t k = std::max(k0, kk); k < j; k++) s -= L(i, k) * L(j, k);
      L(i, j) = s / ljj;
    }
  }
  f.positive_definite = true;
  return f;
}

// x with A x = b from A's band Cholesky factor; O(n p).
inline std::vector<double> band_cholesky_solve(const BandCholesky& f, std::span<const double> b) {
  if (!f.positive_definite) throw std::runtime_error("calx::linalg::band_cholesky_solve: the matrix is not positive definite");
  const std::size_t n = f.l.n, p = f.l.p;
  if (b.size() != n) throw std::invalid_argument("calx::linalg::band_cholesky_solve: b has one entry per row");
  std::vector<double> x(b.begin(), b.end());
  for (std::size_t i = 0; i < n; i++) {
    for (std::size_t k = i > p ? i - p : 0; k < i; k++) x[i] -= f.l(i, k) * x[k];
    x[i] /= f.l(i, i);
  }
  for (std::size_t i = n; i-- > 0;) {
    for (std::size_t k = i + 1; k <= std::min(n - 1, i + p); k++) x[i] -= f.l(k, i) * x[k];
    x[i] /= f.l(i, i);
  }
  return x;
}

// One call: x with A x = b for A symmetric positive definite and banded.
inline std::vector<double> band_solve(const SymmetricBand& A, std::span<const double> b) {
  const BandCholesky f = band_cholesky(A);
  if (!f.positive_definite) throw std::runtime_error("calx::linalg::band_solve: the matrix is not positive definite");
  return band_cholesky_solve(f, b);
}

}  // namespace calx::linalg
