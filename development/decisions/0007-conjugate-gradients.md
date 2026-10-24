# 0007. Conjugate Gradients on Compressed Sparse Rows

- **Status:** Accepted
- **Date:** 2026-10-24
- **Decided by:** the user asked for an economics domain that investigates the circulation of value on maps and in spacetime, under the infrastructure's rule that domains reuse the foundations. Its Hodge decompositions (physis `discretisation`) reach complexes of thousands of cells, beyond dense solves. Claude placed the solver here, on the existing `Csr`, after agreeing the files with the session working on calx's linear algebra, and wrote this record; the composition is The Science Infrastructures' [ADR 0044](https://github.com/huangwanhong-maker/The-Science-Infrastructures/blob/main/development/decisions/0044-value-in-space-and-time.md).

## 1. Context

The Laplacians of a complex (on vertices, d₀ᵀW d₀; on 2-cells, d₁W⁻¹d₁ᵀ) are sparse, symmetric and positive semidefinite, with known kernels: the locally constant functions, and the 2-cycles. A dense solve costs O(n³) and n² memory, too much at n of a few thousand. `linalg` had compressed sparse rows (`Csr`, decision 0004, which the neuroscience engine uses for synapses) and its product, but no solver.

## 2. Decision

`conjugate_gradient.hpp` adds `conjugate_gradient(A, a, b, tolerance, max_iterations, kernel)` on a `Csr` with its values beside the columns.

- **Method.** The preconditioned conjugate-gradient method with the Jacobi preconditioner (Hestenes and Stiefel 1952; Saad 2003, Algorithm 9.1). It stops at a relative residual ‖b − Ax‖ ≤ tolerance · ‖b‖.
- **Semidefinite systems.** For a consistent b, the method started from 0 converges to a solution (Kaasschieter 1988). Given an orthonormal basis of the kernel, b and the solution are projected off it, so that the solution returned is the one of least norm.
- **The result.** `ConjugateGradientResult` carries the solution, the iterations, the true relative residual of the solution returned, and whether it converged.
- **Purely additive.** `sparse.hpp` is unchanged.

## 3. Alternatives Considered

- **A borrowed solver (Eigen's ConjugateGradient).** Eigen is licence-compatible and could sit behind this interface later; for one method on calx's own CSR, the dependency is not yet worth it.
- **A sparse Cholesky factorisation.** Better for many right-hand sides on one matrix; it needs an ordering (AMD, nested dissection) and is planned.
- **The solver in physis.** It names no structure of physics and other domains will need it, so it belongs here (calx's scope: numerical methods on arrays).

## 4. Consequences

- `linalg` gains the header and `tests/check_conjugate_gradient.cpp` (7 checks): a random symmetric positive definite system of order 60 against the dense solve; the path Laplacian against the exact parabola of −u″ = 1, converged within n steps; a singular ring Laplacian, its least-norm solution of zero sum against the dense solve made definite by the kernel; a zero right-hand side.
- physis's Hodge decomposition uses it beyond 1500 vertices or 2-cells (physis decision 0008).
