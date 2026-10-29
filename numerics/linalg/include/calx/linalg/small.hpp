// Small dense systems of a size fixed at compile time, on std::array, without allocation: the linear
// systems of a Newton step in a few unknowns (an implicit integrator's stage, a constraint's multipliers).
//
// solve(A, b) is Gaussian elimination with partial pivoting (Golub and Van Loan 2013, Matrix
// Computations, 4th ed., §3.4); a pivot that is exactly zero is reported, as dense.hpp's solve does.

#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>

namespace calx::linalg {

template <std::size_t N>
using SmallVector = std::array<double, N>;
template <std::size_t N>
using SmallMatrix = std::array<std::array<double, N>, N>;

// x = A⁻¹ b.
template <std::size_t N>
inline SmallVector<N> solve(SmallMatrix<N> A, SmallVector<N> b) {
  for (std::size_t c = 0; c < N; c++) {
    std::size_t piv = c;
    for (std::size_t r = c + 1; r < N; r++)
      if (std::fabs(A[r][c]) > std::fabs(A[piv][c])) piv = r;
    if (A[piv][c] == 0) throw std::runtime_error("calx::linalg::solve: singular matrix");
    std::swap(A[c], A[piv]);
    std::swap(b[c], b[piv]);
    for (std::size_t r = c + 1; r < N; r++) {
      const double f = A[r][c] / A[c][c];
      if (f == 0) continue;
      for (std::size_t k = c; k < N; k++) A[r][k] -= f * A[c][k];
      b[r] -= f * b[c];
    }
  }
  SmallVector<N> x{};
  for (std::size_t r = N; r-- > 0;) {
    double s = b[r];
    for (std::size_t k = r + 1; k < N; k++) s -= A[r][k] * x[k];
    x[r] = s / A[r][r];
  }
  return x;
}

}  // namespace calx::linalg
