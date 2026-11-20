# 0020. General Band Matrices and Their LU Factorisation

- **Status:** Accepted, provisionally (2026-11-21).
- **Date:** 2026-11-21
- **Decided by:** Claude, for kinesis's densities moved by a generator (its decision 0009), at the infrastructure's Stage 7 of evolutionary allocation (its ADR 0139).

## 1. Context

A density on a grid of cells moved by drift and diffusion, in types between which it switches, is moved by a Markov generator whose matrix has bandwidth equal to the number of types; so is the implicit step of a value equation solved by the upwind scheme. These matrices are nonsymmetric. calx factorised symmetric band matrices (band Cholesky) and solved symmetric sparse systems (conjugate gradients), and solved general systems densely; a general band matrix of a few thousand rows had no solver but the dense one, at O(n³).

## 2. Decision

`linalg/band.hpp`:

- `Band{n, kl, ku}`: rows stored with room for the fill that row interchanges bring (2kl + ku + 1 entries a row), entries addressed within the declared band, read as zero outside it;
- `multiply`, `transpose`;
- `band_lu`: Gaussian elimination with partial pivoting within the band, the Gauss transforms stored where they eliminated and U of upper bandwidth ku + kl (Golub and Van Loan 2013, §4.3.5), refusing a singular matrix;
- `band_solve`, from a factorisation or a matrix.

## 3. Alternatives Considered

- **A sparse LU or GMRES.** More general, and more than a banded generator needs; a sparse direct solver is a project of its own, and an iterative one needs preconditioning and tolerances that a banded LU does not.
- **Writing the solver in kinesis.** It names no density and no law, and would be written again by the next library that met a banded system.

## 4. Consequences

`check_calx_band` (17 checks): random band matrices of six shapes against the dense solver, with their residuals relative to the size of A and x; a zero diagonal, where every step pivots; a singular matrix refused; the transpose; the discrete Poisson problem against x(1 − x)/2 exactly.
