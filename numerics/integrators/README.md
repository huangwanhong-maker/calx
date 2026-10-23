# integrators: Ordinary Differential Equations

*Status: draft (2026-10-23): explicit one-step methods.* Implicit, IMEX, differential-algebraic, delay and stochastic methods are planned (the infrastructure's map of the mathematics, §4).

## 1. Purpose

Methods for x' = f(x) on vectors of fixed small dimension, knowing nothing of f. A system as an object, with its events and invariants, is physis's (`physis/dynamics`), which integrates by these methods.

## 2. Content

`runge_kutta.hpp`: `euler` (order 1), `heun` (order 2), `rk4` (order 4), and `dopri5`, the embedded pair of Dormand and Prince, of order 5(4), with error control per component, taking as many substeps of a step as its tolerance asks and carrying the last substep forward.

## 3. Checks

The observed order of each fixed-step method on the damped oscillator; Dormand–Prince against its tolerance.

## 4. Dependencies

None.

## 5. References

- Dormand, J. R. and Prince, P. J. (1980), A family of embedded Runge–Kutta formulae, *J. Comput. Appl. Math.* 6, 19–26.
- Hairer, E., Nørsett, S. P. and Wanner, G. (1993), *Solving Ordinary Differential Equations I*, 2nd ed., Springer.
