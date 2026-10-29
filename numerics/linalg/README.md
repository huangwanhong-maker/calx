# linalg: Numerical Linear Algebra

*Status: draft (2026-10-24): small dense matrices, the exponential, compressed sparse rows, conjugate gradients, the symmetric eigenproblem; rectangular matrices, QR with column pivoting and least squares, Cholesky and band Cholesky factorisations, the eigenvalues of general matrices ([decision 0008](../../development/decisions/0008-optimisation-and-least-squares.md)).* The rest of the module's plan (factorisations, eigenproblems, iterative solvers, structured systems) is in the infrastructure's [map of the mathematics](https://github.com/huangwanhong-maker/The-Science-Infrastructures/blob/main/docs/mathematics.md), §4.

## 1. Purpose

Linear algebra in floating point, dense and sparse. Exact linear algebra over the integers and finite fields belongs to morphe.

## 2. Content

| Header | Content |
|---|---|
| `small.hpp` | `SmallMatrix<N>`, `SmallVector<N>` on `std::array`, and `solve<N>` by Gaussian elimination with partial pivoting, without allocation: the systems of a Newton step in a few unknowns ([decision 0012](../../development/decisions/0012-symplectic-integrators-and-small-systems.md)) |
| `dense.hpp` | `Matrix` of run-time order, row-major; products, sums, the 1-norm; `solve` by Gaussian elimination with partial pivoting; `lu` (the factors in place, the permutation's sign), `determinant` and `log_determinant` (log \|det\| and the sign) |
| `expm.hpp` | The exponential, by scaling and squaring with the Padé approximant of degree 13 (Higham 2005) |
| `eigen.hpp` | `symmetric_eigen`: the eigenvalues (ascending) and orthonormal eigenvectors of a symmetric matrix, by Householder tridiagonalisation and the implicit QL algorithm with Wilkinson's shifts (tred2 and tql2: Wilkinson and Reinsch 1971; EISPACK; JAMA); backward stable, O(n³) ([decision 0006](../../development/decisions/0006-lanes-task-graphs-and-the-genome.md)) |
| `sparse.hpp` | `Csr`: compressed sparse rows whose values live in arrays the caller keeps beside the columns; `compress` by a stable counting sort, with its permutation; `multiply` |
| `rectangular.hpp` | `Rectangular`, m × n row-major (`from_columns`, `from_rows`, `column`, `row`, `transpose`); products, `multiply`, `multiply_transpose`, `gram` (AᵀA); for square matrices `transpose`, `inverse`, `multiply` ([decision 0008](../../development/decisions/0008-optimisation-and-least-squares.md)) |
| `qr.hpp` | `qr(A, pivoting, tolerance)` → `Qr{factors, beta, permutation, rank}`: Householder QR with column pivoting (Golub and Van Loan 2013, Algorithms 5.1.1, 5.2.1, 5.4.1; LAPACK's norm downdating); `apply_qt`, `apply_q`, `thin_q`, `r_factor`, `solve_upper`, `solve_upper_transpose`, `upper_inverse`, `least_squares` (the basic solution), `normal_inverse` ((AᵀA)⁻¹); backward stable, O(mn²) |
| `cholesky.hpp` | `cholesky(A)` → `Cholesky{l, positive_definite}`, `cholesky_solve` (Golub and Van Loan, Algorithm 4.2.1); `SymmetricBand`, `band_cholesky`, `band_cholesky_solve`, `band_solve` for symmetric positive definite band matrices, O(np²) |
| `eigenvalues.hpp` | `eigenvalues(A)`: the complex eigenvalues of a real matrix by balancing, Hessenberg reduction and the Francis double-shift QR iteration (Parlett and Reinsch 1969; orthes and hqr of Wilkinson and Reinsch 1971, as EISPACK and JAMA have them); `balance`, `hessenberg`, `spectral_radius` |
| `conjugate_gradient.hpp` | `conjugate_gradient(A, a, b, tolerance, max_iterations, kernel)`: Jacobi-preconditioned conjugate gradients for symmetric positive semidefinite systems on a `Csr` (Hestenes and Stiefel 1952; Saad 2003, Alg. 9.1), the least-norm solution of a singular consistent system when its kernel is given; `ConjugateGradientResult` ([decision 0007](../../development/decisions/0007-conjugate-gradients.md)) |

Large dense problems will go to blocked kernels, or to a borrowed BLAS and LAPACK, behind this interface.

## 3. Checks

A dense solve by its residual; determinants against triangular matrices, permutations, det(AB) = det A det B, Vandermonde's product and the Hilbert matrix of order 6 (through its logarithm), and a singular matrix; the eigenvalues and sine eigenvectors of the second-difference matrix, AV = VΛ and VᵀV = I on random symmetric matrices of orders 1 to 120, repeated eigenvalues, and exp(Qt) of a reversible rate matrix by its spectral decomposition against `expm`; the exponential against closed forms (diagonal, nilpotent, a damped rotation at small and large norm) and against exp(A) exp(−A) = I; compressed rows against the dense product, and the stability of their order. Conjugate gradients (`check_conjugate_gradient.cpp`, 7 checks): a random positive definite system of order 60 against the dense solve; the path Laplacian against the parabola of −u″ = 1, within n steps; a singular ring Laplacian's least-norm solution against the dense solve made definite by its kernel; a zero right-hand side. QR (`check_qr.cpp`, 27): QᵀQ = I, AP = QR and a non-increasing |R_kk| on random matrices of three shapes; least squares against a line in closed form, the normal equations and the residual's orthogonality; Läuchli's matrix (singular normal equations, exact QR solution); a square system against elimination; the rank of a matrix with a dependent column; (AᵀA)⁻¹. Cholesky (`check_cholesky.cpp`, 9): L Lᵀ = A and the solve; the Hilbert matrix's determinant 1/266716800000; an indefinite matrix refused; the Hodrick–Prescott band against the dense solve; the inverse of the second-difference matrix in closed form. Eigenvalues (`check_eigenvalues.cpp`, 18): companion matrices with real and complex roots; a rotation; a symmetric matrix against `symmetric_eigen`; trace and determinant of random matrices of orders 5–30, conjugate pairs exact; a matrix scaled by 10^±6; a Jordan block; autoregressions' spectral radii. Small systems (`check_small.cpp`, 4): random systems of orders 1 to 8 against the dense solve and their residuals; a zero first pivot exchanged; a singular matrix refused.

## 4. Dependencies

None (the checks use `random`).

## 5. References

- Golub, G. H. and Van Loan, C. F. (2013), *Matrix Computations*, 4th ed., Johns Hopkins.
- Higham, N. J. (2005), The scaling and squaring method for the matrix exponential revisited, *SIAM J. Matrix Anal. Appl.* 26, 1179–1193.
- Wilkinson, J. H. and Reinsch, C. (1971), *Handbook for Automatic Computation II: Linear Algebra*, Springer, contributions II/2 (tred2) and II/3 (tql2).
- Saad, Y. (2003), *Iterative Methods for Sparse Linear Systems*, 2nd ed., SIAM.
- Businger, P. and Golub, G. H. (1965), Linear least squares solutions by Householder transformations, *Numer. Math.* 7, 269–276.
- Drmač, Z. and Bujanović, Z. (2008), On the failure of rank-revealing QR factorization software, *ACM TOMS* 35, 12.
- Higham, N. J. (2002), *Accuracy and Stability of Numerical Algorithms*, 2nd ed., SIAM.
- Parlett, B. N. and Reinsch, C. (1969), Balancing a matrix for calculation of eigenvalues and eigenvectors, *Numer. Math.* 13, 293–304.

## 6. Open Questions

1. Balancing before the exponential (Parlett and Reinsch) for badly scaled matrices; physis's propagators scale by the components' physical sizes instead.
2. When to borrow Eigen or a BLAS for large dense problems.
