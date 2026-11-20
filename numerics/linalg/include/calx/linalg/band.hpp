// General band matrices and their LU factorisation with partial pivoting (decision 0020).
//
// A matrix with kl diagonals below the main one and ku above is stored by rows, each row the 2kl + ku + 1 entries from
// column i − kl to column i + ku + kl: the last kl diagonals hold nothing at first and receive the fill that row
// interchanges bring (Golub and Van Loan 2013, Matrix Computations, 4th ed., §4.3.5: with partial pivoting U has upper
// bandwidth ku + kl, and L is kept as Gauss transforms of lower bandwidth kl). Elimination at column k pivots among rows
// k … k + kl, swaps the two rows on columns k … k + ku + kl, and stores each multiplier where it eliminated; the solve
// applies the interchanges and the transforms in the same order, then back-substitutes with U. O(n·kl·(ku + kl)) to
// factorise, O(n·(2kl + ku)) to solve.
//
// The matrices of Markov generators on a grid of cells (a density moved by drift and diffusion, a value equation's
// implicit step) are of this kind, and nonsymmetric, so that Cholesky and conjugate gradients do not apply.
#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <vector>

namespace calx::linalg {

struct Band {
  std::size_t n = 0, kl = 0, ku = 0;
  std::vector<double> a;  // row i, column j at a[i·width() + (j − i + kl)], for −kl ≤ j − i ≤ ku + kl

  Band() = default;
  Band(std::size_t order, std::size_t lower, std::size_t upper) : n(order), kl(lower), ku(upper), a(order * (2 * lower + upper + 1), 0.0) {}
  std::size_t width() const { return 2 * kl + ku + 1; }
  // Whether (i, j) lies within the declared band (not the fill).
  bool inside(std::size_t i, std::size_t j) const { return i < n && j < n && j + kl >= i && j <= i + ku; }
  double& operator()(std::size_t i, std::size_t j) {
    if (!inside(i, j)) throw std::out_of_range("calx::linalg::Band: an entry outside the band");
    return a[i * width() + (j + kl - i)];
  }
  double operator()(std::size_t i, std::size_t j) const { return inside(i, j) ? a[i * width() + (j + kl - i)] : 0.0; }
};

// y = A x.
inline std::vector<double> multiply(const Band& A, std::span<const double> x) {
  if (x.size() != A.n) throw std::invalid_argument("calx::linalg::multiply: a vector of the band matrix's order");
  std::vector<double> y(A.n, 0.0);
  for (std::size_t i = 0; i < A.n; i++) {
    const std::size_t j0 = i >= A.kl ? i - A.kl : 0, j1 = std::min(A.n - 1, i + A.ku);
    double s = 0;
    for (std::size_t j = j0; j <= j1; j++) s += A(i, j) * x[j];
    y[i] = s;
  }
  return y;
}

// Aᵀ, whose lower and upper bandwidths are A's upper and lower.
inline Band transpose(const Band& A) {
  Band T(A.n, A.ku, A.kl);
  for (std::size_t i = 0; i < A.n; i++) {
    const std::size_t j0 = i >= A.kl ? i - A.kl : 0, j1 = std::min(A.n - 1, i + A.ku);
    for (std::size_t j = j0; j <= j1; j++) T(j, i) = A(i, j);
  }
  return T;
}

struct BandLu {
  Band f;                          // U on and above the diagonal (to ku + kl), the multipliers below
  std::vector<std::size_t> pivot;  // at step k, row k was exchanged with row pivot[k]
};

inline BandLu band_lu(Band A) {
  const std::size_t n = A.n, kl = A.kl, ku = A.ku, w = A.width();
  BandLu out;
  out.pivot.resize(n);
  auto at = [&](std::size_t i, std::size_t j) -> double& { return A.a[i * w + (j + kl - i)]; };
  for (std::size_t k = 0; k < n; k++) {
    const std::size_t last = std::min(n - 1, k + kl), reach = std::min(n - 1, k + ku + kl);
    std::size_t p = k;
    for (std::size_t i = k + 1; i <= last; i++)
      if (std::fabs(at(i, k)) > std::fabs(at(p, k))) p = i;
    if (at(p, k) == 0) throw std::runtime_error("calx::linalg::band_lu: singular matrix");
    out.pivot[k] = p;
    if (p != k)
      for (std::size_t j = k; j <= reach; j++) std::swap(at(k, j), at(p, j));
    const double d = at(k, k);
    for (std::size_t i = k + 1; i <= last; i++) {
      const double l = at(i, k) / d;
      at(i, k) = l;
      if (l == 0) continue;
      for (std::size_t j = k + 1; j <= reach; j++) at(i, j) -= l * at(k, j);
    }
  }
  out.f = std::move(A);
  return out;
}

// x with A x = b, from A's factorisation.
inline std::vector<double> band_solve(const BandLu& F, std::span<const double> b) {
  const Band& A = F.f;
  const std::size_t n = A.n, kl = A.kl, ku = A.ku, w = A.width();
  if (b.size() != n) throw std::invalid_argument("calx::linalg::band_solve: a right-hand side of the matrix's order");
  std::vector<double> x(b.begin(), b.end());
  auto at = [&](std::size_t i, std::size_t j) { return A.a[i * w + (j + kl - i)]; };
  for (std::size_t k = 0; k < n; k++) {
    if (F.pivot[k] != k) std::swap(x[k], x[F.pivot[k]]);
    const std::size_t last = std::min(n - 1, k + kl);
    for (std::size_t i = k + 1; i <= last; i++) x[i] -= at(i, k) * x[k];
  }
  for (std::size_t i = n; i-- > 0;) {
    const std::size_t reach = std::min(n - 1, i + ku + kl);
    double s = x[i];
    for (std::size_t j = i + 1; j <= reach; j++) s -= at(i, j) * x[j];
    x[i] = s / at(i, i);
  }
  return x;
}

inline std::vector<double> band_solve(const Band& A, std::span<const double> b) { return band_solve(band_lu(A), b); }

}  // namespace calx::linalg
