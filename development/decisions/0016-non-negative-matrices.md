# 0016. Non-Negative Matrices: the Perron Root, Its Sensitivities, and Processes That Leave Their States

- **Status:** Accepted
- **Date:** 2026-11-12
- **Decided by:** the infrastructure's ecology domain, founded for urban forests (the infrastructure's ADR 0101), projects a population of trees by size class with a stage matrix, as the demography domain projects a population by age with a Leslie matrix. Both need the dominant eigenvalue of a non-negative matrix and its vectors; the demography domain had found them by power iteration inside its own `projection`. The pattern having appeared twice, the method went down to calx, which it names nothing of: Claude wrote it here, with this record, as the user asked for the ecology domain.

## 1. Context

A non-negative square matrix A has a real eigenvalue λ = ρ(A) with non-negative right and left vectors (Perron–Frobenius). Structured populations (stages, ages, sizes), Markov chains, input–output models and next-generation operators of epidemics are such matrices, and what is asked of them is the same: λ, the stable structure w, the reproductive values v, how λ answers each entry (∂λ/∂a_ij = v_i w_j / vᵀw), and, for the part of a process that keeps its members, how long they stay ((I − U)⁻¹) and how many new members each makes (R₀).

Power iteration, which the demography domain used, converges only when λ strictly dominates the other eigenvalues in modulus: it oscillates for an imprimitive matrix (a population whose members reproduce only at even ages, a cyclic chain) and may stall on a reducible one (post-reproductive stages).

## 2. Decision

- **`linalg/nonnegative.hpp`**, in `calx::linalg`:
  - `irreducible`, `period`, `primitive`, from the graph of positive entries (strong connection; the period as the gcd of level differences along edges);
  - `perron(A)` → `Perron{lambda, right, left, damping, simple, irreducible, period}`: λ as calx's spectral radius (balancing, Hessenberg, Francis QR), then the vectors by inverse iteration with a shift just above λ, from LU factors formed once; w sums to one and vᵀw = 1;
  - `sensitivities(P)`, `elasticities(A, P)`, refused when λ is not simple;
  - `fundamental(U)`, `expected_time_to_leave(U)` and `next_generation(U, F)` = ρ(F (I − U)⁻¹), refusing ρ(U) ≥ 1.
- Columns index the state one comes from, as n(t+1) = A n(t) acts.

## 3. Alternatives Considered

- **Power iteration, as the demography domain had it.** Rejected as the general method: it fails on imprimitive matrices and gives no left vector without a second run; inverse iteration after the QR's λ serves every case in two or three solves.
- **The method in the ecology domain, beside demography's.** Rejected: it would be written twice, and it names no population.
- **In physis's `dynamics`, as the linear map of a discrete dynamical system.** Rejected: nothing here is physical; it is linear algebra of a sign-structured matrix.

## 4. Consequences

- **`check_calx_nonnegative`** (35 checks): closed forms of the 2 × 2 case; 200 random irreducible matrices (residuals 10⁻¹⁴ relative, elasticities summing to one, sensitivities against differences); an imprimitive and a reducible matrix; periods of known graphs; a Leslie matrix against Euler–Lotka, its stable structure and R₀; the fundamental matrix against first-step analysis; the zero matrix; refusals.
- **A repair found while checking.** `eigenvalues` never finished on the zero matrix: its test for a negligible subdiagonal compares with ε‖H‖, which is zero there (the test as JAMA transcribes it, which calx followed). It now returns the zeros at once. The next-generation matrix of a population that makes no new members (street trees, which are planted) is such a matrix.
- **Interfaces.** Additive; nothing existing changes. The demography domain's `projection` may later take its λ from here (not changed in this record).
