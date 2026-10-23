# solvers: Roots and Nonlinear Systems

*Status: draft (2026-10-23): bracketed roots.* Newton and Broyden for systems, fixed points and continuation are planned; optimisation is a module of its own (`optimisation`, planned).

## 1. Content

`roots.hpp`: `brent`, Brent's method for a root bracketed by a sign change, superlinear on smooth functions and never slower than bisection, to within 2 eps |x| + tolerance / 2.

## 2. Checks

Roots known in closed form or to many digits, to the method's guarantee; a discontinuous sign change; the refusal of an unbracketed interval.

## 3. Dependencies

None.

## 4. References

- Brent, R. P. (1973), *Algorithms for Minimization without Derivatives*, Prentice-Hall, ch. 4.
- Press, W. H., Teukolsky, S. A., Vetterling, W. T. and Flannery, B. P. (2007), *Numerical Recipes*, 3rd ed., Cambridge, §9.3.
