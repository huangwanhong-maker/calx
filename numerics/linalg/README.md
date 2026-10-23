# linalg: Numerical Linear Algebra

*Status: draft (2026-10-23): small dense matrices, the exponential, compressed sparse rows.* The rest of the module's plan (factorisations, eigenproblems, iterative solvers, structured systems) is in the infrastructure's [map of the mathematics](https://github.com/huangwanhong-maker/The-Science-Infrastructures/blob/main/docs/mathematics.md), §4.

## 1. Purpose

Linear algebra in floating point, dense and sparse. Exact linear algebra over the integers and finite fields belongs to morphe.

## 2. Content

| Header | Content |
|---|---|
| `dense.hpp` | `Matrix` of run-time order, row-major; products, sums, the 1-norm; `solve` by Gaussian elimination with partial pivoting |
| `expm.hpp` | The exponential, by scaling and squaring with the Padé approximant of degree 13 (Higham 2005) |
| `sparse.hpp` | `Csr`: compressed sparse rows whose values live in arrays the caller keeps beside the columns; `compress` by a stable counting sort, with its permutation; `multiply` |

Large dense problems will go to blocked kernels, or to a borrowed BLAS and LAPACK, behind this interface.

## 3. Checks

A dense solve by its residual; the exponential against closed forms (diagonal, nilpotent, a damped rotation at small and large norm) and against exp(A) exp(−A) = I; compressed rows against the dense product, and the stability of their order.

## 4. Dependencies

None (the checks use `random`).

## 5. References

- Golub, G. H. and Van Loan, C. F. (2013), *Matrix Computations*, 4th ed., Johns Hopkins.
- Higham, N. J. (2005), The scaling and squaring method for the matrix exponential revisited, *SIAM J. Matrix Anal. Appl.* 26, 1179–1193.
- Saad, Y. (2003), *Iterative Methods for Sparse Linear Systems*, 2nd ed., SIAM.

## 6. Open Questions

1. Balancing before the exponential (Parlett and Reinsch) for badly scaled matrices; physis's propagators scale by the components' physical sizes instead.
2. When to borrow Eigen or a BLAS for large dense problems.
