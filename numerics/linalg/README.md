# linalg: Numerical Linear Algebra

*Status: draft (2026-10-24): small dense matrices, the exponential, compressed sparse rows, conjugate gradients, the symmetric eigenproblem.* The rest of the module's plan (factorisations, eigenproblems, iterative solvers, structured systems) is in the infrastructure's [map of the mathematics](https://github.com/huangwanhong-maker/The-Science-Infrastructures/blob/main/docs/mathematics.md), §4.

## 1. Purpose

Linear algebra in floating point, dense and sparse. Exact linear algebra over the integers and finite fields belongs to morphe.

## 2. Content

| Header | Content |
|---|---|
| `dense.hpp` | `Matrix` of run-time order, row-major; products, sums, the 1-norm; `solve` by Gaussian elimination with partial pivoting; `lu` (the factors in place, the permutation's sign), `determinant` and `log_determinant` (log \|det\| and the sign) |
| `expm.hpp` | The exponential, by scaling and squaring with the Padé approximant of degree 13 (Higham 2005) |
| `eigen.hpp` | `symmetric_eigen`: the eigenvalues (ascending) and orthonormal eigenvectors of a symmetric matrix, by Householder tridiagonalisation and the implicit QL algorithm with Wilkinson's shifts (tred2 and tql2: Wilkinson and Reinsch 1971; EISPACK; JAMA); backward stable, O(n³) ([decision 0006](../../development/decisions/0006-lanes-task-graphs-and-the-genome.md)) |
| `sparse.hpp` | `Csr`: compressed sparse rows whose values live in arrays the caller keeps beside the columns; `compress` by a stable counting sort, with its permutation; `multiply` |
| `conjugate_gradient.hpp` | `conjugate_gradient(A, a, b, tolerance, max_iterations, kernel)`: Jacobi-preconditioned conjugate gradients for symmetric positive semidefinite systems on a `Csr` (Hestenes and Stiefel 1952; Saad 2003, Alg. 9.1), the least-norm solution of a singular consistent system when its kernel is given; `ConjugateGradientResult` ([decision 0007](../../development/decisions/0007-conjugate-gradients.md)) |

Large dense problems will go to blocked kernels, or to a borrowed BLAS and LAPACK, behind this interface.

## 3. Checks

A dense solve by its residual; determinants against triangular matrices, permutations, det(AB) = det A det B, Vandermonde's product and the Hilbert matrix of order 6 (through its logarithm), and a singular matrix; the eigenvalues and sine eigenvectors of the second-difference matrix, AV = VΛ and VᵀV = I on random symmetric matrices of orders 1 to 120, repeated eigenvalues, and exp(Qt) of a reversible rate matrix by its spectral decomposition against `expm`; the exponential against closed forms (diagonal, nilpotent, a damped rotation at small and large norm) and against exp(A) exp(−A) = I; compressed rows against the dense product, and the stability of their order. Conjugate gradients (`check_conjugate_gradient.cpp`, 7 checks): a random positive definite system of order 60 against the dense solve; the path Laplacian against the parabola of −u″ = 1, within n steps; a singular ring Laplacian's least-norm solution against the dense solve made definite by its kernel; a zero right-hand side.

## 4. Dependencies

None (the checks use `random`).

## 5. References

- Golub, G. H. and Van Loan, C. F. (2013), *Matrix Computations*, 4th ed., Johns Hopkins.
- Higham, N. J. (2005), The scaling and squaring method for the matrix exponential revisited, *SIAM J. Matrix Anal. Appl.* 26, 1179–1193.
- Wilkinson, J. H. and Reinsch, C. (1971), *Handbook for Automatic Computation II: Linear Algebra*, Springer, contributions II/2 (tred2) and II/3 (tql2).
- Saad, Y. (2003), *Iterative Methods for Sparse Linear Systems*, 2nd ed., SIAM.

## 6. Open Questions

1. Balancing before the exponential (Parlett and Reinsch) for badly scaled matrices; physis's propagators scale by the components' physical sizes instead.
2. When to borrow Eigen or a BLAS for large dense problems.
