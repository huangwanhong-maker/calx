// The exponential of a matrix, by scaling and squaring with the diagonal Padé approximant
// of degree 13 (Higham 2005, "The scaling and squaring method for the matrix exponential
// revisited", SIAM J. Matrix Anal. Appl. 26, 1179-1193, Algorithm 2.3 with m = 13).
//
// A is scaled by 2^-s so that its 1-norm is at most theta_13 = 5.371920351148152, the
// approximant r_13(X) = q_13(X)^-1 p_13(X) is formed with p_13(X) = sum_j c_j X^j,
// q_13(X) = p_13(-X), c_j = (2m - j)! m! / ((2m)! j! (m - j)!), and the result is squared s
// times. Higham's bound gives a backward error below the unit roundoff in exact
// arithmetic; the checks compare with closed forms to about 1e-14. O(n^3 log ||A||).
//
// Simplification: degree 13 is used for every norm, where Higham takes lower degrees for
// small norms to save products; the accuracy is the same and the matrices here are small.

#pragma once

#include <array>
#include <cmath>

#include "dense.hpp"

namespace calx::linalg {

inline Matrix expm(const Matrix& A) {
  constexpr int m = 13;
  constexpr double theta = 5.371920351148152;
  std::array<double, m + 1> c{};
  c[0] = 1;
  for (int j = 1; j <= m; j++) c[j] = c[j - 1] * (m - j + 1) / (static_cast<double>(j) * (2 * m - j + 1));
  const double nrm = norm1(A);
  int s = 0;
  if (nrm > theta) s = static_cast<int>(std::ceil(std::log2(nrm / theta)));
  const Matrix X = std::ldexp(1.0, -s) * A;
  const std::size_t n = A.n;
  Matrix P = Matrix::identity(n), Q = Matrix::identity(n), Xk = Matrix::identity(n);
  for (int j = 1; j <= m; j++) {
    Xk = Xk * X;
    P = P + c[j] * Xk;
    Q = Q + ((j % 2) ? -c[j] : c[j]) * Xk;
  }
  Matrix R = solve(Q, P);
  for (int k = 0; k < s; k++) R = R * R;
  return R;
}

}  // namespace calx::linalg
