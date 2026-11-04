// Checks of the closed-form symmetric 3 × 3 eigenproblem (symmetric3.hpp; decision 0015).
//
// - 10,000 random symmetric matrices, entries uniform in [−1, 1]: the eigenvalues agree with the Householder–QL
//   solver (eigen.hpp) to 10⁻¹³ of the matrix's norm; every pair satisfies ‖A v − λ v‖ ≤ 10⁻¹² ‖A‖; the
//   vectors are orthonormal and right-handed.
// - Matrices with a double eigenvalue, Q diag(λ, λ, μ) Qᵀ, and a multiple of the identity: the values exact
//   to rounding, the frame still orthonormal, A V = V Λ.
// - Scale: a matrix multiplied by 10⁶ or 10⁻⁶ has its eigenvalues multiplied alike, relatively to rounding.
// - The Hessian of a tube along x, diag(0, −1, −1) turned by a rotation: the eigenvector of the eigenvalue
//   near 0 is the tube's axis.

#include <algorithm>
#include <cmath>
#include <vector>

#include <calx/linalg/eigen.hpp>
#include <calx/linalg/symmetric3.hpp>
#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

using calx::verification::Suite;
namespace cl = calx::linalg;

namespace {

double norm(const cl::Symmetric3& A) {
  return std::sqrt(A.a * A.a + A.b * A.b + A.c * A.c + 2 * (A.d * A.d + A.e * A.e + A.f * A.f));
}
std::array<double, 3> apply(const cl::Symmetric3& A, const std::array<double, 3>& v) {
  return {A.a * v[0] + A.d * v[1] + A.e * v[2], A.d * v[0] + A.b * v[1] + A.f * v[2], A.e * v[0] + A.f * v[1] + A.c * v[2]};
}
// Residual of the pairs, orthonormality and handedness of the frame.
void frame_errors(const cl::Symmetric3& A, const cl::Eigen3& E, double& residual, double& ortho, double& hand) {
  residual = ortho = hand = 0;
  for (int k = 0; k < 3; k++) {
    const auto Av = apply(A, E.vectors[k]);
    double r = 0;
    for (int i = 0; i < 3; i++) r += (Av[i] - E.values[k] * E.vectors[k][i]) * (Av[i] - E.values[k] * E.vectors[k][i]);
    residual = std::max(residual, std::sqrt(r));
    for (int j = 0; j < 3; j++) {
      double d = 0;
      for (int i = 0; i < 3; i++) d += E.vectors[k][i] * E.vectors[j][i];
      ortho = std::max(ortho, std::fabs(d - (k == j ? 1.0 : 0.0)));
    }
  }
  const auto& u = E.vectors[0];
  const auto& v = E.vectors[1];
  const auto& w = E.vectors[2];
  const double det = u[0] * (v[1] * w[2] - v[2] * w[1]) - u[1] * (v[0] * w[2] - v[2] * w[0]) + u[2] * (v[0] * w[1] - v[1] * w[0]);
  hand = std::fabs(det - 1);
}
// Q diag(l) Qᵀ for the rotation about a unit axis by an angle.
cl::Symmetric3 rotated(std::array<double, 3> l, std::array<double, 3> axis, double angle) {
  const double n = std::sqrt(axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]);
  for (auto& x : axis) x /= n;
  const double c = std::cos(angle), s = std::sin(angle), t = 1 - c;
  const double x = axis[0], y = axis[1], z = axis[2];
  const double Q[3][3] = {{t * x * x + c, t * x * y - s * z, t * x * z + s * y},
                          {t * x * y + s * z, t * y * y + c, t * y * z - s * x},
                          {t * x * z - s * y, t * y * z + s * x, t * z * z + c}};
  double M[3][3] = {};
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      for (int k = 0; k < 3; k++) M[i][j] += Q[i][k] * l[k] * Q[j][k];
  return {M[0][0], M[1][1], M[2][2], M[0][1], M[0][2], M[1][2]};
}

}  // namespace

int main() {
  Suite s;
  {
    calx::random::Stream z(15, 1, 0, 0);
    double worst_value = 0, worst_residual = 0, worst_ortho = 0, worst_hand = 0;
    for (int trial = 0; trial < 10000; trial++) {
      const cl::Symmetric3 A{z.uniform(-1, 1), z.uniform(-1, 1), z.uniform(-1, 1), z.uniform(-1, 1), z.uniform(-1, 1), z.uniform(-1, 1)};
      cl::Matrix M(3);
      M(0, 0) = A.a, M(1, 1) = A.b, M(2, 2) = A.c;
      M(0, 1) = M(1, 0) = A.d, M(0, 2) = M(2, 0) = A.e, M(1, 2) = M(2, 1) = A.f;
      const auto ref = cl::symmetric_eigen(M);
      const auto E = cl::eigen3(A);
      const double n = norm(A);
      for (int k = 0; k < 3; k++) worst_value = std::max(worst_value, std::fabs(E.values[k] - ref.values[2 - k]) / n);
      double r, o, h;
      frame_errors(A, E, r, o, h);
      worst_residual = std::max(worst_residual, r / n), worst_ortho = std::max(worst_ortho, o), worst_hand = std::max(worst_hand, h);
    }
    s.check("10,000 random matrices: the eigenvalues as Householder–QL finds them, relative to ‖A‖", worst_value, 1e-13);
    s.check("every pair: ‖A v − λ v‖ / ‖A‖", worst_residual, 1e-12);
    s.check("the vectors orthonormal", worst_ortho, 1e-12);
    s.check("and right-handed (det − 1)", worst_hand, 1e-12);
  }
  {
    double worst_value = 0, worst_residual = 0, worst_ortho = 0;
    for (const auto& l : {std::array<double, 3>{2, 2, -1}, std::array<double, 3>{-3, 0.5, 0.5}, std::array<double, 3>{1e-3, 1e-3, 7}}) {
      const auto A = rotated(l, {1, 2, 3}, 0.9);
      const auto E = cl::eigen3(A);
      auto sorted = l;
      std::sort(sorted.begin(), sorted.end(), std::greater<>());
      for (int k = 0; k < 3; k++) worst_value = std::max(worst_value, std::fabs(E.values[k] - sorted[k]));
      double r, o, h;
      frame_errors(A, E, r, o, h);
      worst_residual = std::max(worst_residual, r), worst_ortho = std::max(worst_ortho, o);
    }
    s.check("double eigenvalues: the values exact to rounding", worst_value, 1e-13);
    s.check("and A V = V Λ with an orthonormal frame", std::max(worst_residual, worst_ortho), 1e-12);
    const auto I = cl::eigen3({4, 4, 4, 0, 0, 0});
    s.require("a multiple of the identity: three equal values and the standard frame", I.values == std::array<double, 3>{4, 4, 4} && I.vectors[0][0] == 1);
  }
  {
    const cl::Symmetric3 A{0.3, -0.7, 0.2, 0.5, -0.1, 0.4};
    const auto E = cl::eigenvalues3(A);
    double worst = 0;
    for (double k : {1e6, 1e-6}) {
      const auto F = cl::eigenvalues3({A.a * k, A.b * k, A.c * k, A.d * k, A.e * k, A.f * k});
      for (int i = 0; i < 3; i++) worst = std::max(worst, std::fabs(F[i] / (k * E[i]) - 1));
    }
    s.check("scaled by 10⁶ and 10⁻⁶, the eigenvalues scale alike, relatively", worst, 1e-13);
  }
  {
    const auto H = rotated({0, -1, -1.2}, {0.3, -1, 0.5}, 1.1);  // a bright tube along the rotated x axis
    const auto E = cl::eigen3(H);
    // The axis: the rotated x axis, Q e_x.
    const auto A = rotated({1, 0, 0}, {0.3, -1, 0.5}, 1.1);  // Q diag(1,0,0) Qᵀ = (Q e_x)(Q e_x)ᵀ
    const std::array<double, 3> axis{std::sqrt(A.a), A.d / std::sqrt(A.a), A.e / std::sqrt(A.a)};
    const auto& v = E.vectors[0];  // the eigenvalue 0 is the largest
    s.check("a tube's Hessian: the eigenvector of the eigenvalue near 0 is its axis", 1 - std::fabs(v[0] * axis[0] + v[1] * axis[1] + v[2] * axis[2]), 1e-12);
  }
  return s.done();
}
