# 0015. The Symmetric 3 × 3 Eigenproblem in Closed Form

- **Status:** Accepted, provisionally (trigonometric form with a Jacobi fallback near double eigenvalues)
- **Date:** 2026-11-04
- **Decided by:** the new image foundation eikon (its decision 0002, the infrastructure's ADR 0070) takes the eigenvalues and eigenvectors of the Hessian at every voxel of a volume, to tell tubes from blobs and sheets and to find a vessel's axis. The problem names no image: it is a small symmetric eigenproblem, solved millions of times. Claude wrote it here, with this record, as the user asked for the reconstruction of vessels from images.

## 1. Context

`linalg` held the symmetric eigenproblem of any order by Householder tridiagonalisation and QL (decision 0006): backward stable, but allocating and O(n³) with a large constant, too slow for one problem per voxel of a volume of 10⁷ voxels at several scales. The 3 × 3 case has a closed form: the roots of the characteristic cubic by trigonometry (Smith 1961), and eigenvectors from cross products of the rows of A − λI (Kopp 2008). The closed form loses half the digits where two eigenvalues nearly coincide, which in a Hessian is common (a tube's two cross-sectional eigenvalues are equal), while a symmetric matrix's eigenvalues are well conditioned there.

## 2. Decision

- **`linalg/symmetric3.hpp`.** `Symmetric3{a, b, c, d, e, f}` (the six entries); `eigenvalues3(A)`, descending; `eigen3(A)` → `Eigen3{values, vectors}`, an orthonormal right-handed frame. The trigonometric form serves where the eigenvalues are separated by more than 10⁻⁴ ‖A‖; nearer, the cyclic Jacobi method finishes the problem (Golub and Van Loan 2013, §8.5.2), so that eigenvalues are accurate to a few ulp of ‖A‖ everywhere. Eigenvectors of a double eigenvalue are completed orthonormally from the third.
- No allocation, no exceptions: one call per voxel.

## 3. Alternatives Considered

- **`symmetric_eigen` at order 3.** Rejected for its cost per call (allocation of a `Matrix` and the general iteration).
- **The trigonometric form alone.** Rejected: on matrices with a double eigenvalue it erred by 4 × 10⁻⁸ relative, which the checks caught; tubes' Hessians are such matrices.
- **The solver in eikon.** Rejected: it names no image, and the structure tensor, the inertia of a body and stress tensors need it as much.

## 4. Consequences

- **`check_calx_symmetric3`** (9 checks): 10,000 random matrices against Householder–QL to 10⁻¹³ ‖A‖; every pair's residual to 10⁻¹² ‖A‖; orthonormal right-handed frames; double eigenvalues exact to rounding with A V = V Λ; a multiple of the identity; scaling by 10^±6; a tube's Hessian giving its axis.
- **Interfaces.** Additive; nothing existing changes.
