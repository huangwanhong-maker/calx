# numerics: Numerical Methods

*Status: draft (2026-10-23): the first modules are written (below); the rest are planned.* This is the upper layer of calx. It may use [`compute/`](../compute/) for bulk kernels, but `compute/` never depends on it.

## 1. Purpose

This layer provides general numerical algorithms on floating-point numbers and arrays. The algorithms know nothing of physical meaning. Of algebraic structure, they know only what arrays and functions require.

## 2. Planned Modules

| Module | Scope | Status |
|---|---|---|
| `number_systems` | IEEE formats, double-double, interval arithmetic, compensated summation | planned |
| [`linalg`](linalg/) | Numerical linear algebra: dense (BLAS kernels, factorisations, eigenproblems, least squares), sparse (compressed rows and columns, iterative solvers, preconditioners), functions of matrices, structured systems | draft: small dense and rectangular matrices, LU, solves and determinants, QR with column pivoting and least squares, Cholesky and band Cholesky, the exponential, compressed sparse rows, conjugate gradients, the symmetric eigenproblem and the eigenvalues of general matrices |
| [`random`](random/) | Counter-based generators, streams keyed by purpose and entity, sampling of the standard distributions; quasi-random sequences | draft: Philox4x32-10 and the common distributions, gamma, beta and binomial among them |
| `approximation` | Polynomial, Chebyshev, rational and spline approximation; interpolation in one and many dimensions | planned |
| [`special_functions`](special_functions/) | Error and gamma functions, orthogonal polynomials, Legendre, Bessel, hypergeometric and elliptic functions | draft: erfcx; log Γ, ψ, ψ′, the incomplete gamma and beta functions and their tails; Schmidt's associated Legendre functions |
| [`differentiation`](differentiation/) | Derivatives of functions written as code: forward mode by truncated Taylor arithmetic; reverse mode and transformation of the IR (planned) | draft: jets of order 1 and 2 |
| [`transforms`](transforms/) | Fast Fourier, cosine, wavelet and spherical-harmonic transforms | draft: the FFT of any length |
| [`signal`](signal/) | Digital filters, windows, convolution, resampling, periodograms | draft: the Hann window, Welch's spectral density |
| [`quadrature`](quadrature/) | Integration on intervals, on the sphere, in many dimensions | draft: adaptive Gauss–Kronrod |
| [`integrators`](integrators/) | Ordinary differential equations: explicit, adaptive, implicit and IMEX; differential-algebraic, delay and stochastic equations | draft: Euler, Heun, RK4, Dormand–Prince 5(4) |
| [`solvers`](solvers/) | Roots, nonlinear systems, fixed points, continuation | draft: Brent |
| [`optimisation`](optimisation/) | Unconstrained, constrained, linear, quadratic and conic programming; stochastic and derivative-free methods; least squares | draft: Brent's minimiser in one variable; Nelder–Mead; BFGS with a strong-Wolfe line search; gradients and Hessians by differences; the simplex method for linear programmes |
| [`verification`](verification/) | Checks against references within declared limits; convergence orders; the report of every C++ check | draft |

Some related methods live in sibling libraries:

- Exact integers, rationals and polynomials belong to morphe.
- PDE discretisations belong to physis. So do Lie-group and other geometric integrators that need group structure.
- Dynamical systems as objects (flows, events, hybrid systems) belong to physis (`dynamics`); they integrate by `integrators`.
- Probability and statistics belong to tyche, the inference library; `random` samples, tyche describes.

The place of each module in the whole of the project's mathematics is set in the infrastructure's [map of the mathematics](https://github.com/huangwanhong-maker/The-Science-Infrastructures/blob/main/docs/mathematics.md).

## 3. Candidate External Implementations

Check every licence before adopting a library.

- **Permissive licences**, compatible with MIT: Eigen (MPL 2.0), SUNDIALS (BSD), pocketfft (BSD).
- **Copyleft licences**, which need an explicit decision: FFTW (GPL), SHTns (CeCILL).

## 4. Open Question

Automatic differentiation as a transformation of the IR is an open decision ([CLAUDE.md §8](../CLAUDE.md)); forward differentiation of functions written in C++ is in `differentiation` (decision 0005).
