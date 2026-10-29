# 0011. The θ-Method for Large Linear Systems

- **Status:** Accepted, provisionally (linear systems with a diagonal mass; a nonlinear θ-step with Newton's method is planned)
- **Date:** 2026-10-29
- **Decided by:** the user asked for fields that evolve "in structured way, like PDE", the evolution "separated" from the field. physis's first balance law, heat conduction, needed a stiff integrator for a large sparse system. Claude placed it here, where the method names no physics, and wrote this record. Its use is physis's [decision 0012](https://github.com/huangwanhong-maker/physis/blob/main/development/decisions/0012-fields-as-sections.md).

## 1. Context

`integrators` held explicit one-step methods for small systems (Euler, Heun, RK4, Dormand–Prince). A diffusion law on a cellulation is a large linear system $M\dot y=-Ay+c$:

- $M$ is diagonal (heat capacities);
- $A$ is symmetric positive semidefinite and sparse (a conductance Laplacian);
- its stiffness grows as the square of the cells' number per axis.

An explicit step is stable only below $2/\lambda_\text{max}$, and nothing in calx took an implicit step on compressed sparse rows.

## 2. Decision

`numerics/integrators/include/calx/integrators/theta.hpp`:

- **`theta_step(A, a, m, y, c, h, θ, tol)`** solves
$$(M+\theta hA)\,y'=(M-(1-\theta)hA)\,y+hc$$
for the increment, $(M+\theta hA)(y'-y)=h(c-Ay)$. The solve is by linalg's Jacobi-preconditioned conjugate gradients from zero, the increment being small when the step is. At $\theta=0$ the matrix is $M$, and the step is explicit. It returns the new state, the iterations and the relative residual.
- **`theta_matrix`** assembles $M+\theta hA$ on $A$'s pattern. A matrix without a stored diagonal entry is refused.
- **`theta_factor(θ, z)`** is the stability function $R(z)=(1-(1-\theta)z)/(1+\theta z)$.

Its properties are standard (Hairer and Wanner 1996, §IV.3; Quarteroni and Valli 1994, §11.3):

| θ | Method | Order | Stability |
|---|---|---|---|
| 0 | explicit Euler | 1 | stable for $h\lambda_\text{max}\le2$ |
| ½ | Crank–Nicolson | 2 | A-stable |
| 1 | implicit Euler | 1 | L-stable |

*Rationale.* Solving for the increment rather than the state makes a loose tolerance cost little where the step changes little. A law that must conserve exactly (physis's heat law) does not rely on the solve for it: it updates its conserved quantity by the fluxes of the solved state, so that conservation holds to rounding whatever the tolerance.

## 3. Alternatives Considered

- **A sparse factorisation** (Cholesky on the pattern). Deferred: conjugate gradients suffice for the first laws, and many right-hand sides on one matrix (where a factorisation pays) are not yet needed.
- **SUNDIALS behind the interface.** Deferred: the interface is small, and a borrowed library would be adopted when adaptive implicit methods or DAEs are needed (BSD-licensed, compatible).
- **The method in physis.** Rejected: it names no physics.

## 4. Consequences

- calx has 27 check programs. `check_calx_theta` (8 checks) holds:
  - one step on a diagonal system against $R(h\lambda)$ at θ = 0, ½, 1;
  - the observed orders of the three;
  - their behaviour on a stiff mode: bounded at θ ≥ ½, growing explicitly;
  - one implicit step of length $10^{12}$ reaching the steady state $A^{-1}c$;
  - the solve's convergence;
  - a matrix without a stored diagonal entry refused.
- `integrators` depends on `linalg`.
