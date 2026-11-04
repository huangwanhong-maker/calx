// Eigenvalues and eigenvectors of a real symmetric 3 × 3 matrix in closed form, for the many small problems
// of image analysis: the Hessian or structure tensor at every voxel of a volume (decision 0015).
//
// Eigenvalues (Smith 1961, Eigenvalues of a symmetric 3 × 3 matrix, Commun. ACM 4(4), 168). With q = tr(A)/3,
// B = (A − qI)/p and p² = tr((A − qI)²)/6, the characteristic polynomial of B is λ³ − 3λ − det(B) = 0, whose
// three real roots are 2 cos(φ + 2πk/3) with φ = arccos(det(B)/2)/3 ∈ [0, π/3]. Hence
//   λ₁ = q + 2p cos φ ≥ λ₂ ≥ λ₃ = q + 2p cos(φ + 2π/3),  λ₂ = 3q − λ₁ − λ₃.
// det(B)/2 is clamped to [−1, 1] against rounding. A matrix with p = 0 is qI.
//
// Eigenvectors. For a simple eigenvalue λ, the rows r₀, r₁, r₂ of A − λI span a plane, and the largest of
// their cross products r_i × r_j is normal to it: the eigenvector (Kopp 2008, Efficient numerical
// diagonalization of hermitian 3 × 3 matrices, Int. J. Mod. Phys. C 19, 523, §2.2). Where λ is double the
// rows are parallel; the third vector is then completed orthonormally from the other two, so that the
// result is always an orthonormal, right-handed frame.
//
// Accuracy. Near a double eigenvalue the trigonometric form loses half the digits (arccos has an infinite
// slope at ±1: errors of about √ε ‖A‖), although a symmetric matrix's eigenvalues are well conditioned
// (Weyl). So where two eigenvalues come within 10⁻⁴ ‖A‖, the cyclic Jacobi method finishes the problem
// (Golub and Van Loan 2013, Matrix Computations, 4th ed., §8.5.2), backward stable: eigenvalues to a few ulp
// of ‖A‖ everywhere. Eigenvectors of nearly equal eigenvalues are as ill-determined as the problem makes them.
// Checks: tests/check_symmetric3.cpp.

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <numbers>

namespace calx::linalg {

// The symmetric matrix [[a, d, e], [d, b, f], [e, f, c]] by its six entries.
struct Symmetric3 {
  double a = 0, b = 0, c = 0, d = 0, e = 0, f = 0;  // xx, yy, zz, xy, xz, yz
};

// Eigenvalues by the trigonometric form alone, descending (accurate to √ε ‖A‖ near a double eigenvalue).
inline std::array<double, 3> eigenvalues3_trigonometric(const Symmetric3& A) {
  const double q = (A.a + A.b + A.c) / 3;
  const double off = A.d * A.d + A.e * A.e + A.f * A.f;
  const double p2 = ((A.a - q) * (A.a - q) + (A.b - q) * (A.b - q) + (A.c - q) * (A.c - q) + 2 * off) / 6;
  if (!(p2 > 0)) return {q, q, q};
  const double p = std::sqrt(p2);
  const double ba = (A.a - q) / p, bb = (A.b - q) / p, bc = (A.c - q) / p, bd = A.d / p, be = A.e / p, bf = A.f / p;
  double r = (ba * (bb * bc - bf * bf) - bd * (bd * bc - bf * be) + be * (bd * bf - bb * be)) / 2;
  r = std::clamp(r, -1.0, 1.0);
  const double phi = std::acos(r) / 3;
  const double l1 = q + 2 * p * std::cos(phi);
  const double l3 = q + 2 * p * std::cos(phi + 2 * std::numbers::pi / 3);
  return {l1, 3 * q - l1 - l3, l3};
}

struct Eigen3 {
  std::array<double, 3> values{};                 // descending
  std::array<std::array<double, 3>, 3> vectors{};  // vectors[k] belongs to values[k]; orthonormal, right-handed
};

namespace detail {
// The cyclic Jacobi method on the 3 × 3 matrix: eigenvalues descending, eigenvectors orthonormal and right-handed.
inline void jacobi3(const Symmetric3& S, std::array<double, 3>& values, std::array<std::array<double, 3>, 3>& vectors) {
  double a[3][3] = {{S.a, S.d, S.e}, {S.d, S.b, S.f}, {S.e, S.f, S.c}};
  double v[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  const int pairs[3][2] = {{0, 1}, {0, 2}, {1, 2}};
  for (int sweep = 0; sweep < 32; sweep++) {
    const double off = a[0][1] * a[0][1] + a[0][2] * a[0][2] + a[1][2] * a[1][2];
    const double diag = a[0][0] * a[0][0] + a[1][1] * a[1][1] + a[2][2] * a[2][2];
    if (off == 0 || off <= 1e-36 * diag) break;
    for (const auto& pq : pairs) {
      const int p = pq[0], q = pq[1];
      if (a[p][q] == 0) continue;
      const double tau = (a[q][q] - a[p][p]) / (2 * a[p][q]);
      const double t = (tau >= 0 ? 1.0 : -1.0) / (std::fabs(tau) + std::sqrt(1 + tau * tau));
      const double c = 1 / std::sqrt(1 + t * t), s = t * c;
      for (int k = 0; k < 3; k++) {
        const double akp = a[k][p], akq = a[k][q];
        a[k][p] = c * akp - s * akq, a[k][q] = s * akp + c * akq;
      }
      for (int k = 0; k < 3; k++) {
        const double apk = a[p][k], aqk = a[q][k];
        a[p][k] = c * apk - s * aqk, a[q][k] = s * apk + c * aqk;
      }
      for (int k = 0; k < 3; k++) {
        const double vkp = v[k][p], vkq = v[k][q];
        v[k][p] = c * vkp - s * vkq, v[k][q] = s * vkp + c * vkq;
      }
    }
  }
  int order[3] = {0, 1, 2};
  std::sort(order, order + 3, [&](int i, int j) { return a[i][i] > a[j][j]; });
  for (int k = 0; k < 3; k++) {
    values[k] = a[order[k]][order[k]];
    for (int i = 0; i < 3; i++) vectors[k][i] = v[i][order[k]];
  }
  const auto& u = vectors[0];
  const auto& w = vectors[1];
  const auto& x = vectors[2];
  const double det = u[0] * (w[1] * x[2] - w[2] * x[1]) - u[1] * (w[0] * x[2] - w[2] * x[0]) + u[2] * (w[0] * x[1] - w[1] * x[0]);
  if (det < 0)
    for (auto& z : vectors[2]) z = -z;
}
// Whether two eigenvalues are within 10⁻⁴ of the matrix's scale: then the trigonometric form is not accurate enough.
inline bool near_double(const std::array<double, 3>& l) {
  const double scale = std::max(std::fabs(l[0]), std::fabs(l[2]));
  return scale > 0 && std::min(l[0] - l[1], l[1] - l[2]) < 1e-4 * scale;
}
}  // namespace detail

// Eigenvalues in descending order, to a few ulp of ‖A‖.
inline std::array<double, 3> eigenvalues3(const Symmetric3& A) {
  auto l = eigenvalues3_trigonometric(A);
  if (detail::near_double(l) && (A.d != 0 || A.e != 0 || A.f != 0)) {
    std::array<std::array<double, 3>, 3> V;
    detail::jacobi3(A, l, V);
  }
  return l;
}

namespace detail {
inline std::array<double, 3> cross3(const std::array<double, 3>& u, const std::array<double, 3>& v) {
  return {u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]};
}
inline double norm2(const std::array<double, 3>& u) { return u[0] * u[0] + u[1] * u[1] + u[2] * u[2]; }
// The eigenvector of λ from the largest cross product of the rows of A − λI; false when the rows are parallel.
inline bool vector_of(const Symmetric3& A, double l, std::array<double, 3>& out, double scale) {
  const std::array<double, 3> r0{A.a - l, A.d, A.e}, r1{A.d, A.b - l, A.f}, r2{A.e, A.f, A.c - l};
  const auto c01 = cross3(r0, r1), c02 = cross3(r0, r2), c12 = cross3(r1, r2);
  const double n01 = norm2(c01), n02 = norm2(c02), n12 = norm2(c12);
  const auto& best = n01 >= n02 && n01 >= n12 ? c01 : n02 >= n12 ? c02 : c12;
  const double n = std::max({n01, n02, n12});
  if (!(n > 1e-24 * scale * scale * scale * scale)) return false;
  const double s = 1 / std::sqrt(n);
  out = {best[0] * s, best[1] * s, best[2] * s};
  return true;
}
// A unit vector orthogonal to u (unit).
inline std::array<double, 3> orthogonal(const std::array<double, 3>& u) {
  const std::array<double, 3> t = std::fabs(u[0]) < 0.6 ? std::array<double, 3>{1, 0, 0} : std::array<double, 3>{0, 1, 0};
  auto v = cross3(u, t);
  const double s = 1 / std::sqrt(norm2(v));
  return {v[0] * s, v[1] * s, v[2] * s};
}
}  // namespace detail

inline Eigen3 eigen3(const Symmetric3& A) {
  Eigen3 E;
  E.values = eigenvalues3_trigonometric(A);
  if (detail::near_double(E.values) && (A.d != 0 || A.e != 0 || A.f != 0)) {
    detail::jacobi3(A, E.values, E.vectors);
    return E;
  }
  const double scale = std::max({std::fabs(E.values[0]), std::fabs(E.values[2]), 1e-300});
  std::array<std::array<double, 3>, 3>& V = E.vectors;
  const bool ok0 = detail::vector_of(A, E.values[0], V[0], scale);
  const bool ok2 = detail::vector_of(A, E.values[2], V[2], scale);
  if (ok0 && ok2) {
    V[1] = detail::cross3(V[2], V[0]);
  } else if (ok0) {
    V[1] = detail::orthogonal(V[0]);
    V[2] = detail::cross3(V[0], V[1]);
  } else if (ok2) {
    V[0] = detail::orthogonal(V[2]);
    V[1] = detail::cross3(V[2], V[0]);
  } else {
    V = {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};  // a multiple of the identity: every frame is a frame of eigenvectors
  }
  return E;
}

}  // namespace calx::linalg
