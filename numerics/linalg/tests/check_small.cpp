// Checks of small dense systems (linalg/small.hpp; decision 0012).
//
// - Random systems of orders 1 to 8 against dense.hpp's solve, and their residuals.
// - A system that needs a row exchange (a zero first pivot) solved exactly.
// - A singular matrix refused.

#include <calx/linalg/dense.hpp>
#include <calx/linalg/small.hpp>
#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>

using calx::verification::Suite;
namespace cl = calx::linalg;

template <std::size_t N>
double worst_against_dense(calx::random::Stream& z, double& residual) {
  double worst = 0;
  for (int trial = 0; trial < 50; trial++) {
    cl::SmallMatrix<N> A{};
    cl::SmallVector<N> b{};
    cl::Matrix D(N), B(N);
    for (std::size_t i = 0; i < N; i++) {
      for (std::size_t j = 0; j < N; j++) A[i][j] = D(i, j) = z.normal();
      b[i] = B(i, 0) = z.normal();
    }
    const auto x = cl::solve<N>(A, b);
    const auto y = cl::solve(D, B);
    for (std::size_t i = 0; i < N; i++) {
      worst = std::max(worst, std::fabs(x[i] - y(i, 0)) / std::max(1.0, std::fabs(y(i, 0))));
      double r = -b[i];
      for (std::size_t j = 0; j < N; j++) r += A[i][j] * x[j];
      residual = std::max(residual, std::fabs(r));
    }
  }
  return worst;
}

int main() {
  Suite s;
  calx::random::Stream z(20261029, 0x736d616cu, 0, 0);  // "smal"
  double residual = 0;
  const double worst = std::max({worst_against_dense<1>(z, residual), worst_against_dense<2>(z, residual), worst_against_dense<3>(z, residual),
                                 worst_against_dense<5>(z, residual), worst_against_dense<8>(z, residual)});
  s.check("small: random systems of orders 1 to 8 agree with dense.hpp's solve (relative)", worst, 1e-10);
  s.check("small: and their residuals |A x − b|", residual, 1e-12);
  const auto x = cl::solve<3>({{{0, 1, 0}, {0, 0, 1}, {2, 0, 0}}}, {3, 5, 4});
  s.check("small: a zero first pivot, by a row exchange, solved exactly", std::max({std::fabs(x[0] - 2), std::fabs(x[1] - 3), std::fabs(x[2] - 5)}), 0.0);
  bool refused = false;
  try {
    (void)cl::solve<2>({{{1, 2}, {2, 4}}}, {1, 1});
  } catch (const std::runtime_error&) {
    refused = true;
  }
  s.require("small: a singular matrix is refused", refused);
  return s.done();
}
