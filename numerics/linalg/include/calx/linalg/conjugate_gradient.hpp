// Conjugate gradients for symmetric positive semidefinite systems in CSR.
//
// For A symmetric positive definite, the conjugate-gradient method (Hestenes and Stiefel 1952)
// minimises the A-norm of the error over growing Krylov spaces and reaches A x = b in at most n
// steps in exact arithmetic; with the diagonal of A as preconditioner (Jacobi) it is Saad's
// Algorithm 9.1 (Iterative Methods for Sparse Linear Systems, 2nd ed., 2003). For A positive
// semidefinite and b in its range (the graph and Hodge Laplacians, whose kernels are the locally
// constant functions and the harmonic forms), the method started from 0 converges to a solution
// (Kaasschieter 1988, J. Comput. Appl. Math. 24, 265); when an orthonormal basis of the kernel is
// given, b and the iterate are projected off it, so that the solution returned is the one of
// least norm. The iteration stops when ‖b − A x‖ ≤ tolerance · ‖b‖.
//
// The matrix is a calx Csr with its values stored beside the columns, as multiply() reads them;
// duplicate entries of a row add up.

#pragma once

#include <calx/linalg/sparse.hpp>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <vector>

namespace calx::linalg {

struct ConjugateGradientResult {
  std::vector<double> x;
  std::size_t iterations = 0;
  double relative_residual = 0;  // ‖b − A x‖ / ‖b‖
  bool converged = false;
};

namespace cg_detail {

inline double dot(std::span<const double> a, std::span<const double> b) {
  double s = 0;
  for (std::size_t i = 0; i < a.size(); i++) s += a[i] * b[i];
  return s;
}

inline void project_off(std::vector<double>& v, const std::vector<std::vector<double>>& kernel) {
  for (const auto& k : kernel) {
    const double c = dot(v, k);
    for (std::size_t i = 0; i < v.size(); i++) v[i] -= c * k[i];
  }
}

}  // namespace cg_detail

template <class Index, class Value>
ConjugateGradientResult conjugate_gradient(const Csr<Index>& A, std::span<const Value> a, std::span<const double> b, double tolerance = 1e-12,
                                           std::size_t max_iterations = 0, const std::vector<std::vector<double>>& kernel = {}) {
  using cg_detail::dot;
  const std::size_t n = A.rows;
  if (A.cols != n || b.size() != n || a.size() != A.entries()) throw std::invalid_argument("calx::linalg::conjugate_gradient: a square matrix, its values and a right-hand side of its order");
  if (max_iterations == 0) max_iterations = 10 * n + 100;
  // The Jacobi preconditioner: 1 / diag(A), 1 where the diagonal vanishes.
  std::vector<double> inverse_diagonal(n, 0.0);
  for (std::size_t i = 0; i < n; i++)
    for (std::uint64_t k = A.row[i]; k < A.row[i + 1]; k++)
      if (A.col[k] == i) inverse_diagonal[i] += static_cast<double>(a[k]);
  for (double& d : inverse_diagonal) d = d > 0 ? 1 / d : 1;
  ConjugateGradientResult R;
  R.x.assign(n, 0.0);
  std::vector<double> r(b.begin(), b.end());
  cg_detail::project_off(r, kernel);
  const double norm_b = std::sqrt(dot(r, r));
  if (norm_b == 0) {
    R.converged = true;
    return R;
  }
  std::vector<double> z(n), p(n), q(n);
  for (std::size_t i = 0; i < n; i++) z[i] = inverse_diagonal[i] * r[i];
  p = z;
  double rz = dot(r, z);
  for (R.iterations = 0; R.iterations < max_iterations; R.iterations++) {
    multiply(A, a, std::span<const double>(p), std::span<double>(q));
    const double pq = dot(p, q);
    if (!(pq > 0)) break;  // p in the kernel, or A not positive semidefinite
    const double alpha = rz / pq;
    for (std::size_t i = 0; i < n; i++) {
      R.x[i] += alpha * p[i];
      r[i] -= alpha * q[i];
    }
    R.relative_residual = std::sqrt(dot(r, r)) / norm_b;
    if (R.relative_residual <= tolerance) {
      R.converged = true;
      R.iterations++;
      break;
    }
    for (std::size_t i = 0; i < n; i++) z[i] = inverse_diagonal[i] * r[i];
    const double rz_next = dot(r, z);
    const double beta = rz_next / rz;
    rz = rz_next;
    for (std::size_t i = 0; i < n; i++) p[i] = z[i] + beta * p[i];
  }
  cg_detail::project_off(R.x, kernel);
  // The residual of the solution returned.
  std::vector<double> Ax(n);
  multiply(A, a, std::span<const double>(R.x), std::span<double>(Ax));
  std::vector<double> res(b.begin(), b.end());
  cg_detail::project_off(res, kernel);
  for (std::size_t i = 0; i < n; i++) res[i] -= Ax[i];
  R.relative_residual = std::sqrt(dot(res, res)) / norm_b;
  R.converged = R.relative_residual <= 10 * tolerance;
  return R;
}

}  // namespace calx::linalg
