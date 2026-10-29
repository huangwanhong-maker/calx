# 0008. Optimisation in Several Variables, Linear Programming, and the Matrices of Statistics

- **Status:** Accepted
- **Date:** 2026-10-24
- **Decided by:** the user approved a new foundation library for inference, *tyche* (The Science Infrastructures' [ADR 0048](https://github.com/huangwanhong-maker/The-Science-Infrastructures/tree/main/development/decisions)), and a second, *nomos*, for econometrics, measurement and games, under the infrastructure's rule that a general method is written once, in the foundation where it belongs. Their regressions, likelihoods, time series and games need methods that name no statistical object: minima of functions of several variables, linear programmes, least squares, Cholesky and band factorisations, and the eigenvalues of non-symmetric matrices. Claude placed them in calx and wrote this record.

## 1. Context

`optimisation` held Brent's minimiser in one variable (decision 0006); `linalg` held square dense matrices with LU, solves and determinants, the symmetric eigenproblem, the exponential, compressed rows and conjugate gradients. A regression needs a rectangular design matrix and a stable least-squares solver; a likelihood needs a minimiser in several variables and the Hessian at its minimum; the value of a zero-sum game is a linear programme; the Hodrick–Prescott filter is a pentadiagonal positive definite system; the stability of a vector autoregression is the spectral radius of a non-symmetric companion matrix; orthogonalised impulse responses need a Cholesky factor.

## 2. Decision

### 2.1 linalg

- **`rectangular.hpp`:** `Rectangular`, an m × n row-major matrix (`from_columns`, `from_rows`, `column`, `row`, `transpose`), its products with matrices and vectors (`multiply`, `multiply_transpose`), the Gram matrix `gram(A)` = AᵀA as a square `Matrix`, and for square matrices `transpose`, `inverse` and `multiply`. `Matrix` stays square: changing its shape would silently break every caller that reads `n` as its order.
- **`qr.hpp`:** Householder QR with column pivoting (Golub and Van Loan 2013, Algorithms 5.1.1, 5.2.1, 5.4.1; LAPACK's safeguarded norm downdating), the numerical rank, Qᵀb and Qy, the thin Q, R, triangular solves and inverses, the basic least-squares solution, and (AᵀA)⁻¹ = R⁻¹R⁻ᵀ.
- **`cholesky.hpp`:** the Cholesky factorisation and its solve, and the band Cholesky factorisation of symmetric positive definite band matrices (`SymmetricBand`, `band_cholesky`, `band_solve`), O(np²).
- **`eigenvalues.hpp`:** the eigenvalues of a general real matrix: balancing (Parlett and Reinsch), Householder reduction to Hessenberg form (orthes) and the Francis double-shift QR iteration (hqr, as EISPACK and JAMA have them); `spectral_radius`.

### 2.2 optimisation

- **`nelder_mead.hpp`:** the Nelder–Mead method in the form of Lagarias, Reeds, Wright and Wright (1998), with the adaptive coefficients of Gao and Han (2012) optional; `NelderMeadResult{x, value, iterations, evaluations, converged}`.
- **`quasi_newton.hpp`:** BFGS (Nocedal and Wright 2006, Algorithm 6.1, update 6.17, initial scaling 6.20) with the strong-Wolfe line search of Algorithms 3.5 and 3.6 (cubic interpolation, 3.59), the gradient given or by central differences; `BfgsResult{x, value, gradient_norm, iterations, evaluations, converged, inverse_hessian}`; `gradient_by_differences` and `hessian_by_differences`, the latter for the standard errors of estimates.
- **`linear_programming.hpp`:** the two-phase dense simplex method with Bland's rule, for min cᵀx under ≤ and = rows and any bounds (free variables split, finite upper bounds as rows), returning the status (optimal, infeasible, unbounded, iteration limit), the solution, the value and the duals as marginals ∂value/∂b.

All additive: no existing header changes, except that `calx::optimisation` now links `calx::linalg`.

## 3. Alternatives Considered

- **Borrowing.** Eigen (MPL 2.0) has QR, Cholesky and general eigenvalues, and could sit behind these interfaces later; GLPK is GPL, HiGHS (MIT) is large. For the small and moderate problems of statistics, a few hundred lines each that are checked against published references are worth more than a dependency, and writing them is how they are understood (docs/mathematics.md of the infrastructure, §10.3).
- **Least squares in tyche.** Least squares names no statistical object: the same QR serves curve fitting, interpolation and the domains, so it belongs here.
- **The normal equations.** Cheaper, but they square the condition number; NIST's Longley and Filip data are solved to 11 and 7 digits by the QR and lose most of them by the normal equations.
- **Stability by power iteration.** Unreliable when the dominant eigenvalues are a complex pair of equal modulus, which is common in autoregressions; the QR algorithm is not.
- **An interior-point LP solver.** Better for large sparse programmes; planned when one is needed.

## 4. Consequences

- New checks: `check_calx_qr` (27), `check_calx_cholesky` (9), `check_calx_eigenvalues` (18), `check_calx_nelder_mead` (14), `check_calx_quasi_newton` (17), `check_calx_linear_programming` (40): random factorisations, Läuchli's matrix, the Hilbert matrix, polynomial companions, a balanced badly scaled matrix, Rosenbrock's and Powell's functions, BFGS's inverse Hessian equal to A⁻¹ on a quadratic with exact line searches, Chvátal's examples (including the one that cycles without Bland's rule), Hillier and Lieberman's Wyndor and radiation problems with their shadow prices, infeasible and unbounded programmes, complementary slackness and strong duality, and the values of zero-sum games.
- tyche's regression, GLMs, ARIMA, VAR and HP filter, and nomos's estimators and games, use them.
