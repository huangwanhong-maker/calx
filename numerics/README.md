# numerics: Numerical Methods

*Status: planned.* This is the upper layer of calx. It may use [`compute/`](../compute/) for bulk kernels, but `compute/` never depends on it.

## 1. Purpose

This layer provides general numerical algorithms on floating-point numbers and arrays. The algorithms know nothing of physical meaning. Of algebraic structure, they know only what arrays and functions require.

## 2. Planned Modules

| Module | Scope |
|---|---|
| `number_systems` | IEEE formats, double-double, interval arithmetic, compensated summation |
| `linalg` | Dense and sparse linear algebra, eigenproblems |
| `approximation` | Chebyshev series, polynomial and spline interpolation |
| `special_functions` | Associated Legendre, Bessel and elliptic functions |
| `transforms` | Fast Fourier transforms, spherical-harmonic transforms |
| `quadrature` | Numerical integration, including on the sphere |
| `integrators` | ODE solvers: adaptive, implicit and symplectic |
| `solvers` | Root finding, optimisation, nonlinear systems |

Some related methods live in sibling libraries:

- Exact integers, rationals and polynomials belong to morphe.
- PDE discretisations belong to physis. So do Lie-group and other geometric integrators that need group structure.

## 3. Candidate External Implementations

Check every licence before adopting a library.

- **Permissive licences**, compatible with MIT: Eigen (MPL 2.0), SUNDIALS (BSD), pocketfft (BSD).
- **Copyleft licences**, which need an explicit decision: FFTW (GPL), SHTns (CeCILL).

## 4. Open Question

Automatic differentiation is an open decision ([CLAUDE.md §8](../CLAUDE.md)).
