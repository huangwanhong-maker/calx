// Checks of calx/linalg's eigenvalues of general real matrices: the companion matrix of a
// polynomial with known roots (real, and a complex pair); a rotation, whose eigenvalues are
// e^{±iθ}; a random symmetric matrix against the symmetric eigenproblem; the trace and the
// determinant of random nonsymmetric matrices as the sum and product of the eigenvalues; the
// conjugate symmetry of complex pairs; a badly scaled matrix that balancing makes harmless; a
// Jordan block, to the accuracy a defective eigenvalue allows; the spectral radius of a stable
// and an unstable autoregression's companion matrix.

#include <calx/linalg/dense.hpp>
#include <calx/linalg/eigen.hpp>
#include <calx/linalg/eigenvalues.hpp>
#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <string>
#include <vector>

using namespace calx::linalg;
using cd = std::complex<double>;

namespace {

// The companion matrix of z^n + c[0] z^{n−1} + … + c[n−1].
Matrix companion(const std::vector<double>& c) {
  const std::size_t n = c.size();
  Matrix C(n);
  for (std::size_t j = 0; j < n; j++) C(0, j) = -c[j];
  for (std::size_t i = 1; i < n; i++) C(i, i - 1) = 1;
  return C;
}

// The largest distance from each expected value to the nearest computed one.
double match(std::vector<cd> got, const std::vector<cd>& want) {
  double worst = 0;
  for (const cd& w : want) {
    auto it = std::min_element(got.begin(), got.end(), [&](const cd& a, const cd& b) { return std::abs(a - w) < std::abs(b - w); });
    worst = std::max(worst, std::abs(*it - w));
    got.erase(it);
  }
  return worst;
}

}  // namespace

int main() {
  calx::verification::Suite s;

  // (z − 1)(z − 2)(z − 3) = z³ − 6z² + 11z − 6.
  s.check("eigenvalues: the companion matrix of (z − 1)(z − 2)(z − 3)", match(eigenvalues(companion({-6, 11, -6})), {1, 2, 3}), 1e-12);
  // (z² + 1)(z − 0.5) = z³ − 0.5 z² + z − 0.5: roots ±i and 0.5.
  s.check("eigenvalues: a complex pair ±i and 0.5", match(eigenvalues(companion({-0.5, 1, -0.5})), {cd(0, 1), cd(0, -1), cd(0.5, 0)}), 1e-13);

  {
    const double th = 0.7;
    Matrix R(2);
    R.a = {std::cos(th), -std::sin(th), std::sin(th), std::cos(th)};
    s.check("eigenvalues: a rotation by θ has e^{±iθ}", match(eigenvalues(R), {std::polar(1.0, th), std::polar(1.0, -th)}), 1e-15);
  }

  {
    calx::random::Stream rng(20261024, 10, 0, 0);
    const std::size_t n = 15;
    Matrix A(n);
    for (std::size_t i = 0; i < n; i++)
      for (std::size_t j = 0; j <= i; j++) A(i, j) = A(j, i) = rng.normal();
    const auto sym = symmetric_eigen(A);
    std::vector<cd> want(sym.values.begin(), sym.values.end());
    const auto got = eigenvalues(A);
    double imag = 0;
    for (const cd& l : got) imag = std::max(imag, std::fabs(l.imag()));
    s.check("eigenvalues: a random symmetric matrix of order 15 against the symmetric eigenproblem", match(got, want), 1e-12);
    s.check("eigenvalues: and they are real", imag, 1e-12);
  }

  for (std::size_t n : {5u, 10u, 30u}) {
    calx::random::Stream rng(20261024, 11, static_cast<std::uint32_t>(n), 0);
    Matrix A(n);
    for (double& v : A.a) v = rng.normal();
    const auto l = eigenvalues(A);
    cd sum = 0, prod = 1;
    double trace = 0;
    for (std::size_t i = 0; i < n; i++) trace += A(i, i);
    for (const cd& v : l) {
      sum += v;
      prod *= v;
    }
    const std::string tag = " (random nonsymmetric, order " + std::to_string(n) + ")";
    s.check("eigenvalues: Σλ = trace" + tag, std::abs(sum - trace), 1e-12 * static_cast<double>(n));
    s.check("eigenvalues: ∏λ = det, relative" + tag, std::abs(prod - determinant(A)) / std::fabs(determinant(A)), 1e-11);
    // Complex eigenvalues come in exactly conjugate pairs.
    bool paired = true;
    for (const cd& v : l)
      if (v.imag() != 0) paired = paired && std::any_of(l.begin(), l.end(), [&](const cd& w) { return w == std::conj(v); });
    s.require("eigenvalues: complex ones in exact conjugate pairs" + tag, paired);
  }

  {
    // Badly scaled: D A D⁻¹ with D = diag(1, 1e6, 1e-6) has A's eigenvalues 1, 2, 3.
    Matrix A(3);
    A.a = {1, 1, 0, 0, 2, 1, 0, 0, 3};
    const double d[3] = {1, 1e6, 1e-6};
    Matrix B(3);
    for (std::size_t i = 0; i < 3; i++)
      for (std::size_t j = 0; j < 3; j++) B(i, j) = d[i] * A(i, j) / d[j];
    s.check("eigenvalues: a matrix scaled by 1e±6, balanced back", match(eigenvalues(B), {1, 2, 3}), 1e-10);
  }

  {
    Matrix J(2);
    J.a = {1, 1, 0, 1};
    s.check("eigenvalues: a Jordan block of order 2, to ε^{1/2}", match(eigenvalues(J), {1, 1}), 1e-7);
  }

  // x_t = 0.5 x_{t−1} + 0.3 x_{t−2}: stable; 1.2 and −0.1: not (a root 1.1099…).
  s.check("eigenvalues: a stable AR(2)'s companion, spectral radius (√0.25+1.2 + 0.5)/2", spectral_radius(companion({-0.5, -0.3})) - (0.5 + std::sqrt(0.25 + 1.2)) / 2, 1e-14);
  s.require("eigenvalues: an explosive AR(2)'s companion has spectral radius > 1", spectral_radius(companion({-1.2, 0.1})) > 1);
  return s.done();
}
