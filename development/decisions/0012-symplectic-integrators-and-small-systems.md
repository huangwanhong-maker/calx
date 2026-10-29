# 0012. Symplectic Integrators, Free and Constrained, and Small Dense Systems

- **Status:** Accepted, provisionally (separable systems with a diagonal mass; holonomic constraints of full rank)
- **Date:** 2026-10-29
- **Decided by:** the user asked for physical agents and their dynamics in one world, "different dynamics system" for different needs; physis's discrete mechanics (its decision 0013) and the robotics domain's vehicles needed a symplectic step and a constrained one. Claude wrote these and this record, for the infrastructure's [ADR 0054](https://github.com/huangwanhong-maker/The-Science-Infrastructures/blob/main/development/decisions/0054-bodies-dynamics-and-knowledge.md).

## 1. Context

calx's integrators were explicit Runge–Kutta methods and the θ-method (decision 0011). They know nothing of a system's structure, so over long times they let a conserved energy drift. A body held to a surface (a rover on the Moon) needs its constraint kept to rounding, and its constraint force known, since that force is the normal force its wheels' friction depends on. Newton's method in a few unknowns (an implicit step, a constraint's multipliers) needs a solve of a small system without allocation.

## 2. Decision

### 2.1 `integrators/symplectic.hpp`

- `stormer_verlet`: velocity Verlet for $H = \tfrac12 p^{\mathsf T}M^{-1}p + V(q)$, with a diagonal mass. It is symplectic, symmetric and of order 2, and its energy error is bounded (Hairer, Lubich and Wanner 2006).
- `rattle`: RATTLE for $m$ holonomic constraints of full rank. Its first half is solved by Newton's method on the multipliers. Its second projects the momentum onto the constraint's tangent space. It returns both multipliers, the constraint's forces over the two halves.

A momentum-dependent force is the caller's to split off: an exact half step of it on either side, as robotics's Coulomb kicks are.

### 2.2 `linalg/small.hpp`

`SmallMatrix<N>`, `SmallVector<N>` and `solve<N>`, Gaussian elimination with partial pivoting on `std::array`. It refuses a zero pivot, as dense.hpp's `solve` does.

## 3. Alternatives Considered

- **Leave structure-preserving integration to physis.** physis's discrete mechanics derives its steps from an action (decision 0013 there). The trapezoidal rule's is Störmer–Verlet, which physis's checks hold to calx's here. A separable system needs no action or jets, and RATTLE is a method that names no physics, so both belong here.
- **SHAKE without the momentum projection.** Rejected: the momentum would leave the tangent space, and the second multiplier, half the normal force, would be lost.
- **Eigen for small systems.** Unnecessary for systems of order 1 to 12. calx's own keeps no dependency and no allocation.

## 4. Consequences

- `check_calx_symplectic` (10):
  - Störmer–Verlet's modified energy on the oscillator conserved to $10^{-13}$, and its phase $10^5\theta$ with $\cos\theta = 1 - h^2\omega^2/2$;
  - order 2 on the pendulum against Dormand–Prince, and area preserved;
  - an orbit of $e = 0.6$ over a thousand periods without drift;
  - RATTLE on the spherical pendulum: on the sphere and tangent to it to $10^{-15}$, its multipliers the normal force at order 2, reversible.
- `check_calx_small` (4): random systems against dense.hpp, a zero pivot exchanged, a singular matrix refused.
- The robotics domain's dynamic vehicle and physis's variational and rigid-body integrators use them.
