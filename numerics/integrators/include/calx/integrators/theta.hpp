// The θ-method for large linear systems M ẏ = −A y + c (decision 0011).
//
// M is diagonal and positive (a mass, a heat capacity per cell), A symmetric positive semidefinite and
// sparse (a stiffness, a conductance Laplacian), c a forcing held over the step. One step of size h is
//
//   (M + θ h A) y′ = (M − (1 − θ) h A) y + h c,
//
// explicit Euler at θ = 0 (first order; stable for h λ_max ≤ 2, λ_max the largest eigenvalue of M⁻¹A),
// Crank–Nicolson at θ = ½ (second order, A-stable) and implicit Euler at θ = 1 (first order, L-stable)
// (Hairer and Wanner 1996, §IV.3; Quarteroni and Valli 1994, §11.3). On an eigenvector of M⁻¹A of
// eigenvalue λ it multiplies by R(z) = (1 − (1 − θ) z)/(1 + θ z), z = h λ.
//
// The step is solved for the increment, (M + θ h A)(y′ − y) = h (c − A y), by conjugate gradients from
// zero (calx linalg, Jacobi-preconditioned): the increment is small when the step is, which makes zero a
// good start. At θ = 0 the matrix is M and the step is explicit. A must store every diagonal entry.

#pragma once

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <vector>

#include <calx/linalg/conjugate_gradient.hpp>
#include <calx/linalg/sparse.hpp>

namespace calx::integrators {

struct ThetaStep {
  std::vector<double> y;         // y′
  std::size_t iterations = 0;    // of conjugate gradients (0 at θ = 0)
  double relative_residual = 0;  // of the linear solve
  bool converged = true;
};

// The values of M + θ h A on A's pattern; throws if a row lacks its diagonal entry.
template <class Index>
std::vector<double> theta_matrix(const linalg::Csr<Index>& A, std::span<const double> a, std::span<const double> m, double h, double theta) {
  std::vector<double> v(a.size());
  for (std::size_t i = 0; i < A.rows; i++) {
    bool diagonal = false;
    for (std::uint64_t k = A.row[i]; k < A.row[i + 1]; k++) {
      v[k] = theta * h * a[k];
      if (A.col[k] == i && !diagonal) {
        v[k] += m[i];
        diagonal = true;
      }
    }
    if (!diagonal) throw std::invalid_argument("calx::integrators::theta_step: every row of A stores its diagonal entry");
  }
  return v;
}

template <class Index>
ThetaStep theta_step(const linalg::Csr<Index>& A, std::span<const double> a, std::span<const double> m, std::span<const double> y,
                     std::span<const double> c, double h, double theta, double tolerance = 1e-12) {
  const std::size_t n = A.rows;
  if (A.cols != n || m.size() != n || y.size() != n || c.size() != n || a.size() != A.entries())
    throw std::invalid_argument("calx::integrators::theta_step: a square A, its values, and M, y and c of its order");
  if (!(theta >= 0 && theta <= 1)) throw std::invalid_argument("calx::integrators::theta_step: θ in [0, 1]");
  if (!(h > 0)) throw std::invalid_argument("calx::integrators::theta_step: a positive step");
  for (double mi : m)
    if (!(mi > 0)) throw std::invalid_argument("calx::integrators::theta_step: M positive");
  // r = h (c − A y).
  std::vector<double> r(n);
  linalg::multiply(A, a, y, std::span<double>(r));
  for (std::size_t i = 0; i < n; i++) r[i] = h * (c[i] - r[i]);
  ThetaStep S;
  S.y.assign(y.begin(), y.end());
  if (theta == 0) {
    for (std::size_t i = 0; i < n; i++) S.y[i] += r[i] / m[i];
    return S;
  }
  const std::vector<double> v = theta_matrix(A, a, m, h, theta);
  const auto R = linalg::conjugate_gradient(A, std::span<const double>(v), std::span<const double>(r), tolerance);
  for (std::size_t i = 0; i < n; i++) S.y[i] += R.x[i];
  S.iterations = R.iterations;
  S.relative_residual = R.relative_residual;
  S.converged = R.converged;
  return S;
}

// R(z) = (1 − (1 − θ) z)/(1 + θ z): the factor of one step on an eigenvector of M⁻¹A, z = h λ.
inline double theta_factor(double z, double theta) { return (1 - (1 - theta) * z) / (1 + theta * z); }

}  // namespace calx::integrators
