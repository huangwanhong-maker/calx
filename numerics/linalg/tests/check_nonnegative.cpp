// Checks of calx/linalg's non-negative matrices (nonnegative.hpp): the 2 × 2 Perron root and vectors
// in closed form; random positive matrices (λ the spectral radius, A w = λ w and vᵀA = λ vᵀ, positive
// vectors); the sensitivities against differences of λ, and the elasticities summing to one; an
// imprimitive matrix (period 2: power iteration oscillates, inverse iteration does not); a reducible
// matrix; the period of graphs with known cycles; a Leslie matrix, whose λ solves the Euler–Lotka
// equation Σ f_x l_x λ^{−x−1} = 1 and whose R₀ is Σ f_x l_x; the fundamental matrix of one transient
// state, 1/(1 − s), and of a chain against N (I − U) = I; R₀ > 1 exactly when λ > 1.

#include <calx/linalg/dense.hpp>
#include <calx/linalg/eigenvalues.hpp>
#include <calx/linalg/nonnegative.hpp>
#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>
#include <vector>

using namespace calx::linalg;

namespace {

double residual_right(const Matrix& A, const Perron& P) {
  double r = 0;
  for (std::size_t i = 0; i < A.n; i++) {
    double s = -P.lambda * P.right[i];
    for (std::size_t j = 0; j < A.n; j++) s += A(i, j) * P.right[j];
    r = std::max(r, std::fabs(s));
  }
  return r;
}

double residual_left(const Matrix& A, const Perron& P) {
  double r = 0;
  for (std::size_t j = 0; j < A.n; j++) {
    double s = -P.lambda * P.left[j];
    for (std::size_t i = 0; i < A.n; i++) s += P.left[i] * A(i, j);
    r = std::max(r, std::fabs(s));
  }
  return r;
}

}  // namespace

int main() {
  calx::verification::Suite s;
  calx::random::Stream rng(16, 0, 0, 0);

  // ---- 2 × 2 in closed form ----------------------------------------------------------------------------
  {
    const double a = 0.3, b = 2.0, c = 0.4, d = 0.6;
    Matrix A(2);
    A(0, 0) = a, A(0, 1) = b, A(1, 0) = c, A(1, 1) = d;
    const double disc = std::sqrt((a - d) * (a - d) / 4 + b * c);
    const double l1 = (a + d) / 2 + disc, l2 = (a + d) / 2 - disc;
    const auto P = perron(A);
    s.check("2 × 2: λ = (a + d)/2 + √(((a − d)/2)² + bc)", P.lambda - l1, 1e-14);
    // w ∝ (b, λ − a), summing to one.
    const double w0 = b / (b + l1 - a);
    s.check("2 × 2: w ∝ (b, λ − a)", P.right[0] - w0, 1e-13);
    // v ∝ (c, λ − a) as a row: vᵀA = λvᵀ gives v₀ a + v₁ c = λ v₀.
    s.check("2 × 2: v₁/v₀ = (λ − a)/c", P.left[1] / P.left[0] - (l1 - a) / c, 1e-12);
    s.check("2 × 2: damping ratio λ₁/|λ₂|", P.damping - l1 / std::fabs(l2), 1e-12);
    s.require("2 × 2: positive entries make A primitive", primitive(A) && P.irreducible && P.period == 1 && P.simple);
  }

  // ---- random positive matrices ------------------------------------------------------------------------
  {
    double worst_lambda = 0, worst_res = 0, worst_sum = 0, worst_fd = 0, worst_vw = 0;
    bool positive = true, all_primitive = true;
    for (int t = 0; t < 200; t++) {
      const std::size_t n = 2 + rng.below(9);
      Matrix A(n);
      for (double& x : A.a) x = rng.uniform() < 0.3 ? 0.0 : rng.uniform() * 3;
      for (std::size_t i = 0; i < n; i++) A(i, (i + 1) % n) += 0.1;  // a cycle through every state: irreducible
      A(0, 0) += 0.1;                                                   // and a loop: primitive
      const auto P = perron(A);
      worst_lambda = std::max(worst_lambda, std::fabs(P.lambda - spectral_radius(A)) / P.lambda);
      worst_res = std::max({worst_res, residual_right(A, P) / P.lambda, residual_left(A, P) / P.lambda});
      double vw = 0;
      for (std::size_t i = 0; i < n; i++) {
        vw += P.left[i] * P.right[i];
        positive = positive && P.right[i] > 0 && P.left[i] > 0;
      }
      worst_vw = std::max(worst_vw, std::fabs(vw - 1));
      const Matrix E = elasticities(A, P);
      double sum = 0;
      for (double e : E.a) sum += e;
      worst_sum = std::max(worst_sum, std::fabs(sum - 1));
      // ∂λ/∂a_ij by a central difference at a positive entry.
      const std::size_t i = rng.below(n), j = (i + 1) % n;
      const double h = 1e-6;
      Matrix Ap = A, Am = A;
      Ap(i, j) += h, Am(i, j) -= h;
      const double fd = (spectral_radius(Ap) - spectral_radius(Am)) / (2 * h);
      worst_fd = std::max(worst_fd, std::fabs(fd - sensitivities(P)(i, j)));
      all_primitive = all_primitive && primitive(A);
    }
    s.require("random: the cycle and the loop make A primitive", all_primitive);
    s.check("random: λ is the spectral radius (relative)", worst_lambda, 1e-13);
    s.check("random: A w = λ w and vᵀA = λ vᵀ (relative)", worst_res, 1e-12);
    s.check("random: vᵀw = 1", worst_vw, 1e-12);
    s.require("random: w and v positive for irreducible A", positive);
    s.check("random: elasticities sum to one", worst_sum, 1e-12);
    s.check("random: sensitivities against central differences", worst_fd, 1e-6);
  }

  // ---- imprimitive: period 2 ----------------------------------------------------------------------------
  {
    Matrix A(2);
    A(0, 1) = 2, A(1, 0) = 3;  // eigenvalues ±√6
    const auto P = perron(A);
    s.check("imprimitive: λ = √6", P.lambda - std::sqrt(6.0), 1e-14);
    s.require("imprimitive: irreducible, period 2, not primitive", irreducible(A) && period(A) == 2 && !primitive(A));
    s.check("imprimitive: A w = λ w although power iteration would oscillate", residual_right(A, P), 1e-13);
    s.check("imprimitive: damping ratio 1", P.damping - 1, 1e-12);
  }

  // ---- reducible ----------------------------------------------------------------------------------------
  {
    Matrix A(2);
    A(0, 0) = 2, A(0, 1) = 1, A(1, 1) = 1;  // upper triangular: state 1 feeds state 0, never back
    const auto P = perron(A);
    s.require("reducible: not irreducible, period 0", !irreducible(A) && period(A) == 0);
    s.check("reducible: λ = 2", P.lambda - 2, 1e-14);
    s.check("reducible: w = (1, 0)", P.right[0] - 1, 1e-12);
    s.check("reducible: A w = λ w", residual_right(A, P), 1e-12);
    s.check("reducible: vᵀA = λ vᵀ", residual_left(A, P), 1e-12);
  }

  // ---- periods of known graphs -------------------------------------------------------------------------
  {
    // Cycles of lengths 2 and 3 through vertex 0: gcd 1.
    Matrix A(4);
    A(1, 0) = A(0, 1) = 1;                    // 0 → 1 → 0
    A(2, 0) = A(3, 2) = A(0, 3) = 1;          // 0 → 2 → 3 → 0
    s.require("period: cycles of 2 and 3 give period 1", period(A) == 1);
    Matrix B(4);
    B(1, 0) = B(0, 1) = 1;                    // 0 → 1 → 0
    B(2, 0) = B(3, 2) = B(1, 3) = 1;          // 0 → 2 → 3 → 1 → 0: length 4
    s.require("period: cycles of 2 and 4 give period 2", period(B) == 2);
  }

  // ---- a Leslie matrix --------------------------------------------------------------------------------
  {
    // Ages 0..3; survivals s_x into the next age, fertilities f_x (the newborns per member of age x per step).
    const std::vector<double> f{0.0, 1.2, 1.5, 0.8}, sv{0.6, 0.7, 0.5};
    Matrix L(4);
    for (std::size_t j = 0; j < 4; j++) L(0, j) = f[j];
    for (std::size_t j = 0; j < 3; j++) L(j + 1, j) = sv[j];
    const auto P = perron(L);
    // Euler–Lotka: Σ_x f_x l_x λ^{−x−1} = 1 with l_0 = 1, l_{x+1} = l_x s_x.
    double l = 1, e = 0;
    for (std::size_t x = 0; x < 4; x++) {
      e += f[x] * l * std::pow(P.lambda, -static_cast<double>(x) - 1);
      if (x < 3) l *= sv[x];
    }
    s.check("Leslie: λ solves Euler–Lotka, Σ f_x l_x λ^(−x−1) = 1", e - 1, 1e-12);
    // R₀ by the next-generation matrix: U the survivals, F the fertilities.
    Matrix U(4), F(4);
    for (std::size_t j = 0; j < 3; j++) U(j + 1, j) = sv[j];
    for (std::size_t j = 0; j < 4; j++) F(0, j) = f[j];
    double r0 = 0;
    l = 1;
    for (std::size_t x = 0; x < 4; x++) {
      r0 += f[x] * l;
      if (x < 3) l *= sv[x];
    }
    s.check("Leslie: R₀ = ρ(F N) = Σ f_x l_x", next_generation(U, F) - r0, 1e-12);
    s.require("Leslie: R₀ > 1 exactly when λ > 1", (r0 > 1) == (P.lambda > 1));
    // The stable age structure: w_x ∝ l_x λ^{−x}.
    l = 1;
    std::vector<double> stable(4);
    double total = 0;
    for (std::size_t x = 0; x < 4; x++) {
      stable[x] = l * std::pow(P.lambda, -static_cast<double>(x));
      total += stable[x];
      if (x < 3) l *= sv[x];
    }
    double worst = 0;
    for (std::size_t x = 0; x < 4; x++) worst = std::max(worst, std::fabs(P.right[x] - stable[x] / total));
    s.check("Leslie: the stable structure w_x ∝ l_x λ^(−x)", worst, 1e-12);
  }

  // ---- the fundamental matrix -------------------------------------------------------------------------
  {
    Matrix U(1);
    U(0, 0) = 0.9;
    s.check("fundamental: one state kept with probability s is left after 1/(1 − s) steps", expected_time_to_leave(U)[0] - 10, 1e-12);
    // Stages that grow and stay: N (I − U) = I.
    Matrix V(3);
    V(0, 0) = 0.5, V(1, 0) = 0.3, V(1, 1) = 0.7, V(2, 1) = 0.2, V(2, 2) = 0.9;
    const Matrix N = fundamental(V);
    const Matrix I = N * (Matrix::identity(3) - V);
    double worst = 0;
    for (std::size_t i = 0; i < 3; i++)
      for (std::size_t j = 0; j < 3; j++) worst = std::max(worst, std::fabs(I(i, j) - (i == j ? 1.0 : 0.0)));
    s.check("fundamental: N (I − U) = I", worst, 1e-13);
    // By first-step analysis: t₂ = 1/(1 − 0.9) = 10, t₁ = (1 + 0.2 t₂)/(1 − 0.7) = 10,
    // t₀ = (1 + 0.3 t₁)/(1 − 0.5) = 8.
    const auto t = expected_time_to_leave(V);
    s.check("fundamental: expected times by first-step analysis (8, 10, 10)", std::fabs(t[0] - 8) + std::fabs(t[1] - 10) + std::fabs(t[2] - 10), 1e-12);
    bool refused = false;
    try {
      Matrix W(1);
      W(0, 0) = 1.0;
      (void)fundamental(W);
    } catch (const std::invalid_argument&) { refused = true; }
    s.require("fundamental: a state never left is refused", refused);
  }

  // ---- the zero matrix: a part that makes no new members ----------------------------------------------
  {
    Matrix U(3), Z(3);
    U(1, 0) = 0.5, U(2, 1) = 0.5;
    s.check("zero F: R₀ = 0 (the eigenvalues of the zero matrix, which the QR iteration once did not finish)", next_generation(U, Z), 0);
    s.check("zero matrix: λ = 0", perron(Z).lambda, 0);
  }

  // ---- refusals ---------------------------------------------------------------------------------------
  {
    bool refused = false;
    try {
      Matrix A(2);
      A(0, 1) = -1;
      (void)perron(A);
    } catch (const std::invalid_argument&) { refused = true; }
    s.require("a negative entry is refused", refused);
    // Two equal blocks: λ = 1 twice, sensitivities undefined.
    Matrix B = Matrix::identity(2);
    const auto P = perron(B);
    bool undefined = false;
    try { (void)sensitivities(P); } catch (const std::invalid_argument&) { undefined = true; }
    s.require("a repeated λ: sensitivities are refused", !P.simple && undefined);
  }

  return s.done();
}
