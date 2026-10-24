# optimisation: Minima

*Status: draft (2026-10-24): minima of a function of one variable on an interval ([decision 0006](../../development/decisions/0006-lanes-task-graphs-and-the-genome.md)).* Minima of several variables (Nelder–Mead, quasi-Newton with line searches, trust regions), constrained problems and the expectation–maximisation of mixtures are planned.

## 1. Concepts

A minimum without derivatives is found by narrowing an interval that holds it. Golden section narrows it by 0.618 at each evaluation, whatever the function; a parabola through the three best points jumps to the minimum of a smooth function. Brent's method takes the parabola's step where it is safe (inside the interval, and shorter than half the step before last) and golden section's otherwise, so it converges superlinearly near a smooth minimum and never more slowly than golden section. Near a minimum f is flat to second order, so x is known only to about √eps relative to its size.

## 2. Content

| Header | Content |
|---|---|
| `minimise.hpp` | `Minimum{x, f, iterations, evaluations, converged}`; `brent_minimise(f, a, b, tolerance)` (Brent 1973, ch. 5, localmin), `brent_maximise`; `golden_section(f, a, b, tolerance)` |

## 3. Checks

`check_calx_optimisation`: minima known in closed form (a parabola, cos on [2, 4.5], x eˣ) to the method's guarantee; a parabola in a few evaluations, against golden section's many; a minimum at an end; a quartic, to what its flatness allows; a kink; a binomial log-likelihood's maximum; the refusal of an empty interval.

## 4. Dependencies

None.

## 5. References

- Brent, R. P. (1973), *Algorithms for Minimization without Derivatives*, Prentice-Hall, ch. 5.
- Kiefer, J. (1953), Sequential minimax search for a maximum, *Proc. AMS* 4, 502–506.
- Press, W. H., Teukolsky, S. A., Vetterling, W. T. and Flannery, B. P. (2007), *Numerical Recipes*, 3rd ed., Cambridge, §10.2–10.3.
