// The symmetric eigenproblem: A = V Λ Vᵀ for a real symmetric matrix A, with the eigenvalues
// Λ in ascending order and the orthonormal eigenvectors as the columns of V.
//
// Householder reduction to tridiagonal form, accumulating the transformations, then the
// implicit QL algorithm with Wilkinson's shifts on the tridiagonal matrix: the procedures
// tred2 and tql2 of Bowdler, Martin, Reinsch and Wilkinson (Wilkinson and Reinsch 1971,
// Handbook for Automatic Computation II, Linear Algebra, contributions II/2 and II/3), as
// EISPACK has them and as the public-domain JAMA library (Hicks, Moler et al., NIST and
// MathWorks) transcribes them; Golub and Van Loan (2013), Matrix Computations, 4th ed.,
// §8.3. Backward stable: the computed pairs are exact for A + E with ‖E‖ about ε‖A‖, so an
// eigenvalue is found to within about ε‖A‖ and an eigenvector to within ε‖A‖ / gap. O(n³);
// it serves orders up to a few thousand (the principal components of a genotype matrix, the
// rate matrix of a substitution model); only the lower triangle of A is read.
//
// Uses: the spectral decomposition of a reversible rate matrix Q = Π^{-1/2} S Π^{1/2}, with S
// symmetric, gives exp(Qt) for every t from one decomposition; principal components are the
// eigenvectors of a Gram or covariance matrix.

#pragma once

#include <calx/linalg/dense.hpp>

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace calx::linalg {

struct SymmetricEigen {
  std::vector<double> values;  // ascending
  Matrix vectors;              // column k is the eigenvector of values[k]
};

inline SymmetricEigen symmetric_eigen(const Matrix& A) {
  const std::size_t n = A.n;
  SymmetricEigen out;
  out.vectors = Matrix(n);
  if (n == 0) return out;
  Matrix& V = out.vectors;
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t j = 0; j <= i; j++) V(i, j) = V(j, i) = A(i, j);
  std::vector<double> d(n), e(n);

  // ---- tred2: Householder reduction to tridiagonal form --------------------------------------
  for (std::size_t j = 0; j < n; j++) d[j] = V(n - 1, j);
  for (std::size_t i = n - 1; i > 0; i--) {
    double scale = 0, h = 0;
    for (std::size_t k = 0; k < i; k++) scale += std::fabs(d[k]);
    if (scale == 0) {
      e[i] = d[i - 1];
      for (std::size_t j = 0; j < i; j++) {
        d[j] = V(i - 1, j);
        V(i, j) = 0;
        V(j, i) = 0;
      }
    } else {
      for (std::size_t k = 0; k < i; k++) {
        d[k] /= scale;
        h += d[k] * d[k];
      }
      double f = d[i - 1];
      double g = std::sqrt(h);
      if (f > 0) g = -g;
      e[i] = scale * g;
      h -= f * g;
      d[i - 1] = f - g;
      for (std::size_t j = 0; j < i; j++) e[j] = 0;
      for (std::size_t j = 0; j < i; j++) {
        f = d[j];
        V(j, i) = f;
        g = e[j] + V(j, j) * f;
        for (std::size_t k = j + 1; k + 1 <= i; k++) {
          g += V(k, j) * d[k];
          e[k] += V(k, j) * f;
        }
        e[j] = g;
      }
      f = 0;
      for (std::size_t j = 0; j < i; j++) {
        e[j] /= h;
        f += e[j] * d[j];
      }
      const double hh = f / (h + h);
      for (std::size_t j = 0; j < i; j++) e[j] -= hh * d[j];
      for (std::size_t j = 0; j < i; j++) {
        f = d[j];
        g = e[j];
        for (std::size_t k = j; k + 1 <= i; k++) V(k, j) -= (f * e[k] + g * d[k]);
        d[j] = V(i - 1, j);
        V(i, j) = 0;
      }
    }
    d[i] = h;
  }
  // Accumulate the transformations.
  for (std::size_t i = 0; i + 1 < n; i++) {
    V(n - 1, i) = V(i, i);
    V(i, i) = 1;
    const double h = d[i + 1];
    if (h != 0) {
      for (std::size_t k = 0; k <= i; k++) d[k] = V(k, i + 1) / h;
      for (std::size_t j = 0; j <= i; j++) {
        double g = 0;
        for (std::size_t k = 0; k <= i; k++) g += V(k, i + 1) * V(k, j);
        for (std::size_t k = 0; k <= i; k++) V(k, j) -= g * d[k];
      }
    }
    for (std::size_t k = 0; k <= i; k++) V(k, i + 1) = 0;
  }
  for (std::size_t j = 0; j < n; j++) {
    d[j] = V(n - 1, j);
    V(n - 1, j) = 0;
  }
  V(n - 1, n - 1) = 1;
  e[0] = 0;

  // ---- tql2: the implicit QL algorithm on the tridiagonal matrix (d, e) ---------------------------
  for (std::size_t i = 1; i < n; i++) e[i - 1] = e[i];
  e[n - 1] = 0;
  double f = 0, tst1 = 0;
  const double eps = std::ldexp(1.0, -52);
  for (std::size_t l = 0; l < n; l++) {
    tst1 = std::max(tst1, std::fabs(d[l]) + std::fabs(e[l]));
    std::size_t m = l;
    while (m < n && std::fabs(e[m]) > eps * tst1) m++;
    if (m == n) m = n - 1;  // e[n - 1] is zero, so the search ends there
    if (m > l) {
      int iterations = 0;
      do {
        if (++iterations > 64) throw std::runtime_error("calx::linalg::symmetric_eigen: the QL iteration did not converge");
        double g = d[l];
        double p = (d[l + 1] - g) / (2 * e[l]);
        double r = std::hypot(p, 1.0);
        if (p < 0) r = -r;
        d[l] = e[l] / (p + r);
        d[l + 1] = e[l] * (p + r);
        const double dl1 = d[l + 1];
        double h = g - d[l];
        for (std::size_t i = l + 2; i < n; i++) d[i] -= h;
        f += h;
        p = d[m];
        double c = 1, c2 = 1, c3 = 1, s = 0, s2 = 0;
        const double el1 = e[l + 1];
        for (std::size_t i = m; i-- > l;) {
          c3 = c2;
          c2 = c;
          s2 = s;
          g = c * e[i];
          h = c * p;
          r = std::hypot(p, e[i]);
          e[i + 1] = s * r;
          s = e[i] / r;
          c = p / r;
          p = c * d[i] - s * g;
          d[i + 1] = h + s * (c * g + s * d[i]);
          for (std::size_t k = 0; k < n; k++) {
            h = V(k, i + 1);
            V(k, i + 1) = s * V(k, i) + c * h;
            V(k, i) = c * V(k, i) - s * h;
          }
        }
        p = -s * s2 * c3 * el1 * e[l] / dl1;
        e[l] = s * p;
        d[l] = c * p;
      } while (std::fabs(e[l]) > eps * tst1);
    }
    d[l] += f;
    e[l] = 0;
  }

  // Ascending order, the vectors with their values.
  for (std::size_t i = 0; i + 1 < n; i++) {
    std::size_t k = i;
    double p = d[i];
    for (std::size_t j = i + 1; j < n; j++)
      if (d[j] < p) { k = j; p = d[j]; }
    if (k != i) {
      d[k] = d[i];
      d[i] = p;
      for (std::size_t j = 0; j < n; j++) std::swap(V(j, i), V(j, k));
    }
  }
  out.values = std::move(d);
  return out;
}

}  // namespace calx::linalg
