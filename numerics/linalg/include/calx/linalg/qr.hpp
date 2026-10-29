// The QR factorisation by Householder reflections with column pivoting, and linear least squares.
//
// A P = Q R for an m × n matrix A: P a permutation of the columns, Q = H_0 H_1 ⋯ H_{t−1}
// orthogonal (t = min(m, n) reflections H_k = I − β_k v_k v_kᵀ), R upper triangular with
// |R_00| ≥ |R_11| ≥ ⋯ (approximately). The reflections are those of Golub and Van Loan (2013),
// Matrix Computations, 4th ed., Algorithm 5.1.1 (house, with v_k(k) = 1 and R_kk = ‖x‖ ≥ 0,
// computed on the column scaled by its largest entry, so that nothing overflows); the pivoting
// is Businger and Golub's (1965; Golub and Van Loan, Algorithm 5.4.1): at step k the column of
// largest remaining norm is brought forward. The remaining norms are downdated, and recomputed
// when the downdate has lost more than half the digits (the safeguard of LAPACK's xLAQP2,
// after Drmač and Bujanović 2008, ACM TOMS 35, 12).
//
// Accuracy. Householder QR is backward stable: the computed R is the exact factor of A + E with
// ‖E(:, j)‖ ≤ c m n ε ‖A(:, j)‖ column by column (Higham 2002, Accuracy and Stability of
// Numerical Algorithms, 2nd ed., Theorem 19.4), so a least-squares solution is as accurate as
// the problem's conditioning allows: relative error about κ(A) ε when the residual is small,
// κ(A)² ε when it is not (Golub and Van Loan, §5.3). Forming AᵀA instead squares κ always.
// Column pivoting reveals the numerical rank: rank = the number of |R_kk| above
// tolerance · |R_00| (default max(m, n) ε). Scale the columns to comparable norms first if their
// units differ greatly, since the rank decision is not invariant under column scaling.
// O(mn²) flops (2mn² − 2n³/3 without accumulating Q).
//
// least_squares() returns the basic solution of min ‖Ax − b‖₂: the coordinates of the columns
// beyond the rank are set to zero (Golub and Van Loan, §5.4); for a matrix of full column
// rank it is the unique minimiser.

#pragma once

#include <calx/linalg/dense.hpp>
#include <calx/linalg/rectangular.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace calx::linalg {

struct Qr {
  // R on and above the diagonal (its first min(m, n) rows); below the diagonal of column k,
  // the essential part of v_k (its entry k is 1 and not stored).
  Rectangular factors;
  std::vector<double> beta;              // H_k = I − β_k v_k v_kᵀ
  std::vector<std::size_t> permutation;  // column k of A P is column permutation[k] of A
  std::size_t rank = 0;                  // numerical rank, by the tolerance below
  double tolerance = 0;                  // relative tolerance used for the rank
};

// tolerance < 0: max(m, n) ε. pivoting = false keeps the columns in order (P = I); the rank is
// then the count of leading |R_kk| above the tolerance before the first that is not.
inline Qr qr(Rectangular A, bool pivoting = true, double tolerance = -1) {
  const std::size_t m = A.rows, n = A.cols, t = std::min(m, n);
  constexpr double eps = std::numeric_limits<double>::epsilon();
  Qr f;
  f.tolerance = tolerance < 0 ? static_cast<double>(std::max(m, n)) * eps : tolerance;
  f.beta.assign(t, 0.0);
  f.permutation.resize(n);
  std::iota(f.permutation.begin(), f.permutation.end(), std::size_t{0});

  auto column_norm = [&](std::size_t j, std::size_t from) {
    double scale = 0, ssq = 1;  // the scaled sum of squares of LAPACK's xNRM2, without overflow
    for (std::size_t i = from; i < m; i++) {
      const double v = std::fabs(A(i, j));
      if (v == 0) continue;
      if (scale < v) {
        ssq = 1 + ssq * (scale / v) * (scale / v);
        scale = v;
      } else {
        ssq += (v / scale) * (v / scale);
      }
    }
    return scale * std::sqrt(ssq);
  };
  std::vector<double> vn1(n), vn2(n);  // current and reference norms of the remaining columns
  for (std::size_t j = 0; j < n; j++) vn1[j] = vn2[j] = column_norm(j, 0);
  const double tol3z = std::sqrt(eps);

  for (std::size_t k = 0; k < t; k++) {
    if (pivoting) {
      std::size_t p = k;
      for (std::size_t j = k + 1; j < n; j++)
        if (vn1[j] > vn1[p]) p = j;
      if (p != k) {
        for (std::size_t i = 0; i < m; i++) std::swap(A(i, k), A(i, p));
        std::swap(f.permutation[k], f.permutation[p]);
        std::swap(vn1[k], vn1[p]);
        std::swap(vn2[k], vn2[p]);
      }
    }
    // house(x), x = A(k:m, k), on x scaled by its largest entry (Golub and Van Loan, Alg. 5.1.1).
    double s = 0;
    for (std::size_t i = k; i < m; i++) s = std::max(s, std::fabs(A(i, k)));
    double beta = 0, mu = 0;
    if (s > 0) {
      const double x1 = A(k, k) / s;
      double sigma = 0;
      for (std::size_t i = k + 1; i < m; i++) {
        const double v = A(i, k) / s;
        sigma += v * v;
      }
      if (sigma == 0) {
        beta = x1 >= 0 ? 0 : 2;
        mu = std::fabs(x1);
        for (std::size_t i = k + 1; i < m; i++) A(i, k) = 0;
      } else {
        mu = std::sqrt(x1 * x1 + sigma);
        const double v1 = x1 <= 0 ? x1 - mu : -sigma / (x1 + mu);
        beta = 2 * v1 * v1 / (sigma + v1 * v1);
        for (std::size_t i = k + 1; i < m; i++) A(i, k) = (A(i, k) / s) / v1;
      }
      mu *= s;
    }
    f.beta[k] = beta;
    // Apply H_k to the remaining columns: A ← A − β v (vᵀ A).
    if (beta != 0)
      for (std::size_t j = k + 1; j < n; j++) {
        double w = A(k, j);
        for (std::size_t i = k + 1; i < m; i++) w += A(i, k) * A(i, j);
        w *= beta;
        A(k, j) -= w;
        for (std::size_t i = k + 1; i < m; i++) A(i, j) -= w * A(i, k);
      }
    A(k, k) = mu;
    // Downdate the remaining norms (LAPACK xLAQP2).
    if (pivoting)
      for (std::size_t j = k + 1; j < n; j++) {
        if (vn1[j] == 0) continue;
        double temp = std::fabs(A(k, j)) / vn1[j];
        temp = std::max(0.0, 1 - temp * temp);
        const double ratio = vn1[j] / vn2[j];
        if (temp * ratio * ratio <= tol3z) {
          vn1[j] = vn2[j] = column_norm(j, k + 1);
        } else {
          vn1[j] *= std::sqrt(temp);
        }
      }
  }
  // The numerical rank.
  const double r00 = t > 0 ? std::fabs(A(0, 0)) : 0;
  std::size_t r = 0;
  while (r < t && r00 > 0 && std::fabs(A(r, r)) > f.tolerance * r00) r++;
  f.rank = r;
  f.factors = std::move(A);
  return f;
}

// Qᵀ b, b of length m.
inline std::vector<double> apply_qt(const Qr& f, std::span<const double> b) {
  const std::size_t m = f.factors.rows;
  if (b.size() != m) throw std::invalid_argument("calx::linalg::apply_qt: b has one entry per row");
  std::vector<double> y(b.begin(), b.end());
  for (std::size_t k = 0; k < f.beta.size(); k++) {
    if (f.beta[k] == 0) continue;
    double w = y[k];
    for (std::size_t i = k + 1; i < m; i++) w += f.factors(i, k) * y[i];
    w *= f.beta[k];
    y[k] -= w;
    for (std::size_t i = k + 1; i < m; i++) y[i] -= w * f.factors(i, k);
  }
  return y;
}

// Q y, y of length m.
inline std::vector<double> apply_q(const Qr& f, std::span<const double> yin) {
  const std::size_t m = f.factors.rows;
  if (yin.size() != m) throw std::invalid_argument("calx::linalg::apply_q: y has one entry per row");
  std::vector<double> y(yin.begin(), yin.end());
  for (std::size_t k = f.beta.size(); k-- > 0;) {
    if (f.beta[k] == 0) continue;
    double w = y[k];
    for (std::size_t i = k + 1; i < m; i++) w += f.factors(i, k) * y[i];
    w *= f.beta[k];
    y[k] -= w;
    for (std::size_t i = k + 1; i < m; i++) y[i] -= w * f.factors(i, k);
  }
  return y;
}

// The first `columns` columns of Q (default: the rank), an m × columns matrix with orthonormal
// columns spanning the range of A's first `columns` pivoted columns. Its rows give the leverages
// h_ii = ‖Q(i, 0:r)‖² of a regression.
inline Rectangular thin_q(const Qr& f, std::size_t columns = static_cast<std::size_t>(-1)) {
  const std::size_t m = f.factors.rows;
  const std::size_t c = std::min(columns == static_cast<std::size_t>(-1) ? f.rank : columns, std::min(m, f.factors.cols));
  Rectangular Q(m, c);
  std::vector<double> e(m);
  for (std::size_t j = 0; j < c; j++) {
    std::fill(e.begin(), e.end(), 0.0);
    e[j] = 1;
    const auto q = apply_q(f, e);
    for (std::size_t i = 0; i < m; i++) Q(i, j) = q[i];
  }
  return Q;
}

// The leading r × r block of R (default r = rank), as a square Matrix in pivoted order.
inline Matrix r_factor(const Qr& f, std::size_t r = static_cast<std::size_t>(-1)) {
  const std::size_t k = std::min(r == static_cast<std::size_t>(-1) ? f.rank : r, std::min(f.factors.rows, f.factors.cols));
  Matrix R(k);
  for (std::size_t i = 0; i < k; i++)
    for (std::size_t j = i; j < k; j++) R(i, j) = f.factors(i, j);
  return R;
}

// x with U x = b for U upper triangular (back substitution), using U's leading order b.size().
inline std::vector<double> solve_upper(const Matrix& U, std::span<const double> b) {
  const std::size_t k = b.size();
  std::vector<double> x(k);
  for (std::size_t i = k; i-- > 0;) {
    double s = b[i];
    for (std::size_t j = i + 1; j < k; j++) s -= U(i, j) * x[j];
    if (U(i, i) == 0) throw std::runtime_error("calx::linalg::solve_upper: a zero on the diagonal");
    x[i] = s / U(i, i);
  }
  return x;
}

// x with Uᵀ x = b for U upper triangular (forward substitution).
inline std::vector<double> solve_upper_transpose(const Matrix& U, std::span<const double> b) {
  const std::size_t k = b.size();
  std::vector<double> x(k);
  for (std::size_t i = 0; i < k; i++) {
    double s = b[i];
    for (std::size_t j = 0; j < i; j++) s -= U(j, i) * x[j];
    if (U(i, i) == 0) throw std::runtime_error("calx::linalg::solve_upper_transpose: a zero on the diagonal");
    x[i] = s / U(i, i);
  }
  return x;
}

// U⁻¹ for U upper triangular, itself upper triangular.
inline Matrix upper_inverse(const Matrix& U) {
  const std::size_t k = U.n;
  Matrix X(k);
  for (std::size_t j = 0; j < k; j++) {
    if (U(j, j) == 0) throw std::runtime_error("calx::linalg::upper_inverse: a zero on the diagonal");
    X(j, j) = 1 / U(j, j);
    for (std::size_t i = j; i-- > 0;) {
      double s = 0;
      for (std::size_t l = i + 1; l <= j; l++) s += U(i, l) * X(l, j);
      X(i, j) = -s / U(i, i);
    }
  }
  return X;
}

// The basic least-squares solution of min ‖A x − b‖₂, in A's original column order: the
// coordinates of the columns beyond the rank are zero.
inline std::vector<double> least_squares(const Qr& f, std::span<const double> b) {
  const std::size_t r = f.rank;
  const auto c = apply_qt(f, b);
  const auto z = solve_upper(r_factor(f, r), std::span<const double>(c.data(), r));
  std::vector<double> x(f.factors.cols, 0.0);
  for (std::size_t j = 0; j < r; j++) x[f.permutation[j]] = z[j];
  return x;
}

// One call: x minimising ‖A x − b‖₂ (the basic solution when A is rank-deficient).
inline std::vector<double> least_squares(const Rectangular& A, std::span<const double> b) { return least_squares(qr(A), b); }

// (AᵀA)⁻¹ = P R⁻¹ R⁻ᵀ Pᵀ in A's original column order, for A of full column rank; the
// unscaled covariance of least-squares coefficients. With rank r < n, the block of the first r
// pivoted columns is filled and the rest is zero.
inline Matrix normal_inverse(const Qr& f) {
  const std::size_t n = f.factors.cols, r = f.rank;
  const Matrix Ri = upper_inverse(r_factor(f, r));
  Matrix C(n);
  for (std::size_t i = 0; i < r; i++)
    for (std::size_t j = 0; j < r; j++) {
      double s = 0;
      for (std::size_t l = std::max(i, j); l < r; l++) s += Ri(i, l) * Ri(j, l);
      C(f.permutation[i], f.permutation[j]) = s;
    }
  return C;
}

}  // namespace calx::linalg
