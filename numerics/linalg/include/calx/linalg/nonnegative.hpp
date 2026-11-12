// Non-negative matrices: the Perron root and its vectors, their sensitivities, and the matrices of a
// process that leaves its states (decision 0016).
//
// A square matrix A ≥ 0 has a real eigenvalue λ = ρ(A), its spectral radius, with a non-negative right
// vector w and left vector v (Perron 1907; Frobenius 1912). When A is irreducible (its graph of positive
// entries is strongly connected) λ is simple and w, v are positive and unique up to scale; when A is
// moreover primitive (irreducible with period 1, the gcd of its cycles' lengths) every other eigenvalue
// is smaller in modulus, so that A^t x / λ^t → w (vᵀx) for any x ≥ 0 with vᵀx > 0 (Seneta 2006,
// Non-negative Matrices and Markov Chains, Theorems 1.1, 1.5 and 1.2; Horn and Johnson 2013, Matrix
// Analysis, 2nd ed., §8.4–8.5).
//
//   irreducible(A), period(A), primitive(A)   the graph of A's positive entries: strong connection by
//                                             a search forwards and backwards from one vertex; the period
//                                             as the gcd over edges u→v of level(u) + 1 − level(v), levels
//                                             from a breadth-first search (Denardo 1977, Management Sci.
//                                             23: 1029–1035; Jarvis and Shier 1999)
//   perron(A)                                 λ = ρ(A), w (summing to one), v (with vᵀw = 1), the damping
//                                             ratio λ / |λ₂|, and whether λ is simple
//   sensitivities(P), elasticities(A, P)      ∂λ/∂a_ij = v_i w_j / vᵀw, valid for a simple λ (Caswell 2001,
//                                             Matrix Population Models, 2nd ed., ch. 9);
//                                             e_ij = (a_ij / λ) ∂λ/∂a_ij, which sum to one (de Kroon et al.
//                                             1986, Ecology 67: 1427–1431; Caswell ch. 9)
//   fundamental(U)                            N = (I − U)⁻¹ for a substochastic U with ρ(U) < 1: N_ij the
//                                             expected number of steps spent in state i starting from j, so
//                                             that the column sums of N are the expected times to leaving
//                                             (Kemeny and Snell 1960, Finite Markov Chains, ch. III;
//                                             Caswell ch. 5). Columns index the state one starts from, as a
//                                             projection n(t+1) = A n(t) acts on columns.
//   next_generation(U, F)                     R₀ = ρ(F N), the dominant eigenvalue of the next-generation
//                                             matrix: R₀ > 1 exactly when ρ(U + F) > 1 (Cushing and Zhou
//                                             1994, Nat. Resour. Model. 8: 297–333; Caswell ch. 5).
//
// Method. λ is found as calx's spectral radius (balancing, Hessenberg reduction and the Francis double-
// shift QR: backward stable, linalg/eigenvalues.hpp). w and v are then found by inverse iteration with
// the shift σ = λ(1 + 10⁻¹⁰) + 10⁻³⁰⁰, nearer to λ than to any other eigenvalue (they lie in the closed
// disc of radius λ, so their distance from σ exceeds σ − λ unless they equal λ): each iteration
// multiplies the error by at most (σ − λ)/min|σ − λ_k|, so two or three iterations reach rounding
// (Wilkinson 1965, The Algebraic Eigenvalue Problem, §9.31–9.40; Golub and Van Loan 2013, §7.6.1).
// The nearly singular A − σI is what makes inverse iteration converge; its LU factors are formed once.
// This serves every non-negative matrix, primitive, imprimitive or reducible, where power iteration
// would fail to converge on an imprimitive one. A simple λ is decided from the gap between λ and the
// next real part within 10⁻⁸ λ of it; when λ is not simple the vectors returned are one non-negative
// choice among many and the sensitivities are not defined (refused).
//
// Accuracy. λ to within about ε‖A‖ (λ is simple for an irreducible A, and its condition 1/s(λ) =
// ‖v‖‖w‖/vᵀw is modest for the stage matrices of populations); w, v to within about ε‖A‖/gap. O(n³);
// the matrices of structured populations, Markov chains and next-generation operators are small.

#pragma once

#include <calx/linalg/dense.hpp>
#include <calx/linalg/eigenvalues.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace calx::linalg {

namespace nonnegative_detail {

inline void require_nonnegative(const Matrix& A, const char* where) {
  for (double x : A.a)
    if (!(x >= 0)) throw std::invalid_argument(std::string("calx::linalg::") + where + ": a non-negative matrix, with finite entries");
}

// Breadth-first levels from vertex 0 along the edges i → j where A(j, i) > 0 (column j to row i: the
// direction in which n(t+1) = A n(t) carries mass) or, reversed, A(i, j) > 0.
inline std::vector<long> levels(const Matrix& A, bool forward) {
  const std::size_t n = A.n;
  std::vector<long> level(n, -1);
  if (n == 0) return level;
  std::vector<std::size_t> queue{0};
  level[0] = 0;
  for (std::size_t q = 0; q < queue.size(); q++) {
    const std::size_t u = queue[q];
    for (std::size_t x = 0; x < n; x++) {
      const bool edge = forward ? A(x, u) > 0 : A(u, x) > 0;
      if (edge && level[x] < 0) {
        level[x] = level[u] + 1;
        queue.push_back(x);
      }
    }
  }
  return level;
}

// Solves (LU) x = b in place, with the row exchanges as lu() made them: lu() swaps whole rows as it
// goes, so the exchanges are replayed on b from the recorded pivots.
struct Factored {
  Matrix lu;
  std::vector<std::size_t> pivot;
  bool singular = false;
};

inline Factored factor(Matrix A) {
  const std::size_t n = A.n;
  Factored f;
  f.pivot.assign(n, 0);
  for (std::size_t c = 0; c < n; c++) {
    std::size_t p = c;
    for (std::size_t r = c + 1; r < n; r++)
      if (std::fabs(A(r, c)) > std::fabs(A(p, c))) p = r;
    f.pivot[c] = p;
    if (A(p, c) == 0) {
      f.singular = true;
      A(p, c) = 1e-300;  // an exactly singular shift: perturb the pivot, as inverse iteration allows
    }
    if (p != c)
      for (std::size_t j = 0; j < n; j++) std::swap(A(c, j), A(p, j));
    for (std::size_t r = c + 1; r < n; r++) {
      const double m = A(r, c) / A(c, c);
      A(r, c) = m;
      if (m == 0) continue;
      for (std::size_t j = c + 1; j < n; j++) A(r, j) -= m * A(c, j);
    }
  }
  f.lu = std::move(A);
  return f;
}

inline void solve_in_place(const Factored& f, std::vector<double>& b) {
  const std::size_t n = f.lu.n;
  for (std::size_t c = 0; c < n; c++)
    if (f.pivot[c] != c) std::swap(b[c], b[f.pivot[c]]);
  for (std::size_t r = 0; r < n; r++)
    for (std::size_t k = 0; k < r; k++) b[r] -= f.lu(r, k) * b[k];
  for (std::size_t r = n; r-- > 0;) {
    for (std::size_t k = r + 1; k < n; k++) b[r] -= f.lu(r, k) * b[k];
    b[r] /= f.lu(r, r);
  }
}

inline Matrix transpose(const Matrix& A) {
  Matrix T(A.n);
  for (std::size_t i = 0; i < A.n; i++)
    for (std::size_t j = 0; j < A.n; j++) T(j, i) = A(i, j);
  return T;
}

// A non-negative null vector of A − σI by inverse iteration from the all-ones vector, scaled to sum one.
inline std::vector<double> inverse_iteration(const Matrix& A, double sigma, const char* where) {
  const std::size_t n = A.n;
  Matrix S = A;
  for (std::size_t i = 0; i < n; i++) S(i, i) -= sigma;
  const Factored f = factor(std::move(S));
  std::vector<double> x(n, 1.0 / static_cast<double>(n));
  for (int it = 0; it < 6; it++) {
    solve_in_place(f, x);
    double scale = 0;
    for (double y : x) scale += y;
    if (!std::isfinite(scale) || scale == 0) {
      double big = 0;
      for (double y : x) big = std::max(big, std::fabs(y));
      if (!(big > 0) || !std::isfinite(big)) throw std::runtime_error(std::string("calx::linalg::") + where + ": inverse iteration lost the vector");
      for (double& y : x) y /= big;
      continue;
    }
    for (double& y : x) y /= scale;
  }
  // The Perron vector is non-negative: entries of the wrong sign are rounding about zero.
  double total = 0, negative = 0;
  for (double& y : x) {
    if (y < 0) negative = std::max(negative, -y), y = 0;
    total += y;
  }
  if (negative > 1e-8) throw std::runtime_error(std::string("calx::linalg::") + where + ": the vector found has entries of both signs");
  for (double& y : x) y /= total;
  return x;
}

}  // namespace nonnegative_detail

// Whether the graph of A's positive entries is strongly connected (a 1 × 1 matrix is irreducible).
inline bool irreducible(const Matrix& A) {
  nonnegative_detail::require_nonnegative(A, "irreducible");
  if (A.n <= 1) return true;
  const auto f = nonnegative_detail::levels(A, true), b = nonnegative_detail::levels(A, false);
  for (std::size_t i = 0; i < A.n; i++)
    if (f[i] < 0 || b[i] < 0) return false;
  return true;
}

// The period of an irreducible A: the gcd of the lengths of the cycles of its graph. 0 for a reducible
// matrix (its period is not defined), and for the 1 × 1 zero matrix.
inline long period(const Matrix& A) {
  if (!irreducible(A)) return 0;
  const auto level = nonnegative_detail::levels(A, true);
  long g = 0;
  for (std::size_t u = 0; u < A.n; u++)
    for (std::size_t x = 0; x < A.n; x++)
      if (A(x, u) > 0) g = std::gcd(g, std::labs(level[u] + 1 - level[x]));
  return g;
}

inline bool primitive(const Matrix& A) { return period(A) == 1; }

struct Perron {
  double lambda = 0;           // ρ(A)
  std::vector<double> right;   // w ≥ 0, Σ w = 1
  std::vector<double> left;    // v ≥ 0, vᵀw = 1
  double damping = HUGE_VAL;   // λ / |λ₂|, the largest modulus among the other eigenvalues (∞ for none)
  bool simple = true;          // λ a simple eigenvalue: the vectors are unique and the sensitivities defined
  bool irreducible = true;
  long period = 1;             // 0 when reducible
};

inline Perron perron(const Matrix& A) {
  nonnegative_detail::require_nonnegative(A, "perron");
  const std::size_t n = A.n;
  if (n == 0) throw std::invalid_argument("calx::linalg::perron: an empty matrix");
  Perron P;
  P.irreducible = irreducible(A);
  P.period = period(A);
  const auto values = eigenvalues(A);
  double lambda = 0;
  for (const auto& z : values) lambda = std::max(lambda, std::abs(z));
  P.lambda = lambda;
  // λ is real and an eigenvalue; the others: those within 1e−8 λ of λ itself count as λ repeated.
  int at_lambda = 0;
  double second = 0;
  for (const auto& z : values) {
    if (std::abs(z - std::complex<double>(lambda, 0)) <= 1e-8 * std::max(lambda, 1e-300)) at_lambda++;
    else second = std::max(second, std::abs(z));
  }
  P.simple = at_lambda <= 1 || P.irreducible;  // an irreducible matrix's λ is simple (Perron–Frobenius)
  if (P.irreducible) {
    // For an irreducible matrix the count above may include λ found twice by rounding; recompute the
    // second modulus over all values but one copy of λ.
    second = 0;
    bool skipped = false;
    for (const auto& z : values) {
      if (!skipped && std::abs(z - std::complex<double>(lambda, 0)) <= 1e-8 * std::max(lambda, 1e-300)) { skipped = true; continue; }
      second = std::max(second, std::abs(z));
    }
  }
  P.damping = second > 0 ? lambda / second : HUGE_VAL;
  const double sigma = lambda * (1 + 1e-10) + 1e-300;
  P.right = nonnegative_detail::inverse_iteration(A, sigma, "perron");
  P.left = nonnegative_detail::inverse_iteration(nonnegative_detail::transpose(A), sigma, "perron");
  double vw = 0;
  for (std::size_t i = 0; i < n; i++) vw += P.left[i] * P.right[i];
  if (vw > 0)
    for (double& x : P.left) x /= vw;
  else
    P.simple = false;  // vᵀw = 0 only for a λ that is not simple
  return P;
}

// ∂λ/∂a_ij = v_i w_j / vᵀw (with vᵀw = 1 as perron() scales them). Refused when λ is not simple.
inline Matrix sensitivities(const Perron& P) {
  if (!P.simple) throw std::invalid_argument("calx::linalg::sensitivities: λ is not a simple eigenvalue, so it is not differentiable in the entries");
  const std::size_t n = P.right.size();
  Matrix S(n);
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t j = 0; j < n; j++) S(i, j) = P.left[i] * P.right[j];
  return S;
}

// e_ij = (a_ij / λ) ∂λ/∂a_ij: the proportional change of λ for a proportional change of a_ij. Σ e_ij = 1.
inline Matrix elasticities(const Matrix& A, const Perron& P) {
  if (!(P.lambda > 0)) throw std::invalid_argument("calx::linalg::elasticities: λ = 0");
  Matrix E = sensitivities(P);
  for (std::size_t i = 0; i < A.n; i++)
    for (std::size_t j = 0; j < A.n; j++) E(i, j) *= A(i, j) / P.lambda;
  return E;
}

// N = (I − U)⁻¹ for a non-negative U with ρ(U) < 1 (a transient part: survival and passage among states).
inline Matrix fundamental(const Matrix& U) {
  nonnegative_detail::require_nonnegative(U, "fundamental");
  if (!(spectral_radius(U) < 1)) throw std::invalid_argument("calx::linalg::fundamental: ρ(U) < 1 is needed, so that every state is left in finite expected time");
  return solve(Matrix::identity(U.n) - U, Matrix::identity(U.n));
}

// The expected number of steps before leaving, from each state: the column sums of N.
inline std::vector<double> expected_time_to_leave(const Matrix& U) {
  const Matrix N = fundamental(U);
  std::vector<double> t(U.n, 0.0);
  for (std::size_t j = 0; j < U.n; j++)
    for (std::size_t i = 0; i < U.n; i++) t[j] += N(i, j);
  return t;
}

// R₀ = ρ(F (I − U)⁻¹), for a transient part U and a part F that makes new members (births, infections).
inline double next_generation(const Matrix& U, const Matrix& F) {
  nonnegative_detail::require_nonnegative(F, "next_generation");
  if (F.n != U.n) throw std::invalid_argument("calx::linalg::next_generation: U and F of one order");
  return spectral_radius(F * fundamental(U));
}

}  // namespace calx::linalg
