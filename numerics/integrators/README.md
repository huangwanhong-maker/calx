# integrators: Ordinary Differential Equations

*Status: draft (2026-10-29): explicit one-step methods; the θ-method for large linear systems ([decision 0011](../../development/decisions/0011-the-theta-method.md)); symplectic integrators, free and constrained ([decision 0012](../../development/decisions/0012-symplectic-integrators-and-small-systems.md)).* Nonlinear implicit, IMEX, differential-algebraic, delay and stochastic methods are planned (the infrastructure's map of the mathematics, §4).

## 1. Purpose

Methods for x' = f(x) on vectors of fixed small dimension, knowing nothing of f; and for large linear systems M ẏ = −A y + c with a diagonal mass and a sparse symmetric stiffness, the semi-discrete form of diffusion on a cellulation. A system as an object, with its events and invariants, is physis's (`physis/dynamics`), and a balance law is physis's (`physis/continuum`); both integrate by these methods.

## 2. Content

`runge_kutta.hpp`: `euler` (order 1), `heun` (order 2), `rk4` (order 4), and `dopri5`, the embedded pair of Dormand and Prince, of order 5(4), with error control per component, taking as many substeps of a step as its tolerance asks and carrying the last substep forward.

`symplectic.hpp`: `stormer_verlet(M⁻¹, ∇V, q, p, h)`, velocity Verlet for a separable Hamiltonian with a diagonal mass: symplectic, symmetric, of order 2, its energy bounded; `rattle(M⁻¹, ∇V, g, q, p, h, tol)`, RATTLE for holonomic constraints g(q) = 0 of full rank, the first half solved by Newton's method on the multipliers and the momentum projected onto the constraint's tangent space; `RattleStep`, the multipliers of both halves (the constraint's forces), the iterations and the residual.

`theta.hpp`: `theta_step(A, a, m, y, c, h, θ, tol)`, one step of $(M+\theta hA)\,y'=(M-(1-\theta)hA)\,y+hc$ solved for the increment by linalg's preconditioned conjugate gradients (explicit Euler at θ = 0, Crank–Nicolson at ½, implicit Euler at 1), with its iterations and residual; `theta_matrix`, $M+\theta hA$ on $A$'s pattern; `theta_factor`, the stability function $R(z)=(1-(1-\theta)z)/(1+\theta z)$.

## 3. Checks

`check_integrators.cpp`: the observed order of each fixed-step method on the damped oscillator; Dormand–Prince against its tolerance.

`check_symplectic.cpp` (10 checks): Störmer–Verlet's modified energy on the oscillator, conserved to 10⁻¹³ over 10⁵ steps, and its phase; order 2 on the pendulum against Dormand–Prince; area kept; an orbit of e = 0.6 over a thousand periods without drift; RATTLE on the spherical pendulum, on the sphere and tangent to it to 10⁻¹⁵, its multipliers the normal force at order 2, reversible.

`check_theta.cpp` (8 checks): a step against $R(h\lambda)$ at θ = 0, ½, 1; their observed orders; a stiff mode bounded at θ ≥ ½ and growing explicitly; one implicit step of $10^{12}$ reaching the steady state; the solve's convergence; a matrix without its diagonal refused.

## 4. Dependencies

`linalg` (compressed sparse rows, conjugate gradients), for the θ-method; its small dense systems, for RATTLE's multipliers.

## 5. References

- Dormand, J. R. and Prince, P. J. (1980), A family of embedded Runge–Kutta formulae, *J. Comput. Appl. Math.* 6, 19–26.
- Andersen, H. C. (1983), Rattle: a "velocity" version of the Shake algorithm for molecular dynamics calculations, *J. Comput. Phys.* 52, 24–34.
- Hairer, E., Lubich, C. and Wanner, G. (2006), *Geometric Numerical Integration*, 2nd ed., Springer.
- Hairer, E., Nørsett, S. P. and Wanner, G. (1993), *Solving Ordinary Differential Equations I*, 2nd ed., Springer.
- Hairer, E. and Wanner, G. (1996), *Solving Ordinary Differential Equations II*, 2nd ed., Springer.
- Quarteroni, A. and Valli, A. (1994), *Numerical Approximation of Partial Differential Equations*, Springer.
