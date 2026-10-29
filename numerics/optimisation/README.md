# optimisation: Minima and Linear Programmes

*Status: draft (2026-10-24): minima of a function of one variable on an interval ([decision 0006](../../development/decisions/0006-lanes-task-graphs-and-the-genome.md)); minima of several variables without and with derivatives, and linear programmes ([decision 0008](../../development/decisions/0008-optimisation-and-least-squares.md)).* Trust regions, constrained nonlinear problems (SQP, interior point), quadratic and conic programming, nonlinear least squares (Levenberg–Marquardt) and stochastic methods are planned.

## 1. Concepts

A minimum without derivatives is found by narrowing an interval that holds it. Golden section narrows it by 0.618 at each evaluation, whatever the function; a parabola through the three best points jumps to the minimum of a smooth function. Brent's method takes the parabola's step where it is safe (inside the interval, and shorter than half the step before last) and golden section's otherwise, so it converges superlinearly near a smooth minimum and never more slowly than golden section. Near a minimum f is flat to second order, so x is known only to about √ε relative to its size.

In several variables, the Nelder–Mead method moves a simplex of n + 1 points by reflection, expansion, contraction and shrinking, using values only: robust on rough functions of a few variables, with no convergence theory beyond one dimension. With a gradient, a quasi-Newton method builds an approximation H of the inverse Hessian from the changes of the gradient along its steps (the BFGS update, which keeps H positive definite when the step satisfies the Wolfe conditions), and converges superlinearly. At the minimum of a negative log-likelihood, the inverse of the Hessian is the asymptotic covariance of the estimates; it is better computed afresh by differences than taken from BFGS's H, which knows the curvature only along the directions explored.

A linear programme minimises a linear function over a polyhedron; an optimum, if there is one, is at a vertex, and the simplex method walks from vertex to vertex (basis to basis) while the objective decreases. Phase I finds a first vertex or proves there is none; Bland's rule (least index entering and leaving) prevents the cycling that degenerate vertices can cause. The duals y, read from the final basis, price the constraints: ∂(optimal value)/∂b, zero on a slack constraint (complementary slackness) and with bᵀy equal to the optimal value (strong duality). The value of a zero-sum game is such a programme, and its two optimal mixed strategies are a primal and dual pair.

## 2. Content

| Header | Content |
|---|---|
| `minimise.hpp` | `Minimum{x, f, iterations, evaluations, converged}`; `brent_minimise(f, a, b, tolerance)` (Brent 1973, ch. 5, localmin), `brent_maximise`; `golden_section(f, a, b, tolerance)` |
| `nelder_mead.hpp` | `nelder_mead(f, x0, NelderMeadOptions{x_tolerance, f_tolerance, max_evaluations, initial_step, adaptive})` → `NelderMeadResult{x, value, iterations, evaluations, converged}`: Lagarias, Reeds, Wright and Wright (1998); Gao and Han's (2012) adaptive coefficients; fminsearch's initial simplex and stopping test |
| `quasi_newton.hpp` | `bfgs(f, gradient, x0, BfgsOptions{gradient_tolerance, max_iterations, c1, c2, max_step})` and `bfgs(f, x0, options)` (central differences) → `BfgsResult{x, value, gradient_norm, iterations, evaluations, converged, inverse_hessian}`: Nocedal and Wright (2006), Algorithm 6.1 with the strong-Wolfe line search of Algorithms 3.5–3.6; `gradient_by_differences(f, x)`, `hessian_by_differences(f, x)` |
| `linear_programming.hpp` | `simplex(LinearProgram{c, a_ub, b_ub, a_eq, b_eq, lower, upper}, LpOptions)` → `LpResult{status, x, value, duals_ub, duals_eq, iterations}`: two-phase dense simplex with Bland's rule (Dantzig 1963; Bland 1977; Chvátal 1983); `LpStatus{optimal, infeasible, unbounded, iteration_limit}` |

Functions of several variables are called with a `std::span<const double>`; a gradient returns a `std::vector<double>`.

## 3. Checks

- `check_calx_optimisation`: minima known in closed form (a parabola, cos on [2, 4.5], x eˣ) to the method's guarantee; a parabola in a few evaluations, against golden section's many; a minimum at an end; a quartic, to what its flatness allows; a kink; a binomial log-likelihood's maximum; the refusal of an empty interval.
- `check_calx_nelder_mead` (14): Rosenbrock's function from (−1.2, 1) with the standard and the adaptive coefficients; a quadratic in three variables; x eˣ against Brent; Powell's singular function in four variables; the budget of evaluations; NaN as +∞; the caller's initial steps.
- `check_calx_quasi_newton` (17): Rosenbrock with its gradient (superlinear, to 10⁻⁹) and by differences; the extended Rosenbrock function in ten variables; a quadratic's minimiser, with H = A⁻¹ after n exact line searches (Nocedal and Wright §6.1); the normal log-likelihood's maximum and the covariance of (μ̂, σ̂²) from `hessian_by_differences` against σ²/n and 2σ⁴/n; gradients by differences against closed forms.
- `check_calx_linear_programming` (40): Chvátal's (1983) first example (13 at (2, 0, 1), duals (1, 0, 1)) and his cycling example (1 at (1, 0, 1, 0)); Hillier and Lieberman's Wyndor Glass problem (36 at (2, 6), shadow prices (0, 3/2, 1)) and radiation therapy (5.25 at (7.5, 4.5), with = and ≥ rows); an infeasible and an unbounded programme; bounds, a free variable and a redundant equality; complementary slackness, dual feasibility, reduced costs and strong duality at every optimum; rock–paper–scissors (value 0, uniform strategies for both players) and a 2 × 2 game (value 1/7 at (3/7, 4/7)).

## 4. Dependencies

`linalg` (the inverse Hessian is a `Matrix`; the constraint matrices are `Rectangular`).

## 5. References

- Bland, R. G. (1977), New finite pivoting rules for the simplex method, *Math. Oper. Res.* 2, 103–107.
- Brent, R. P. (1973), *Algorithms for Minimization without Derivatives*, Prentice-Hall, ch. 5.
- Chvátal, V. (1983), *Linear Programming*, Freeman.
- Dantzig, G. B. (1963), *Linear Programming and Extensions*, Princeton.
- Gao, F. and Han, L. (2012), Implementing the Nelder–Mead simplex algorithm with adaptive parameters, *Comput. Optim. Appl.* 51, 259–277.
- Hillier, F. S. and Lieberman, G. J., *Introduction to Operations Research*, McGraw-Hill (the Wyndor Glass and radiation-therapy examples, ch. 3–4 and 6).
- Kiefer, J. (1953), Sequential minimax search for a maximum, *Proc. AMS* 4, 502–506.
- Lagarias, J. C., Reeds, J. A., Wright, M. H. and Wright, P. E. (1998), Convergence properties of the Nelder–Mead simplex method in low dimensions, *SIAM J. Optim.* 9, 112–147.
- Moré, J. J., Garbow, B. S. and Hillstrom, K. E. (1981), Testing unconstrained optimization software, *ACM TOMS* 7, 17–41.
- Nelder, J. A. and Mead, R. (1965), A simplex method for function minimization, *Computer J.* 7, 308–313.
- Nocedal, J. and Wright, S. J. (2006), *Numerical Optimization*, 2nd ed., Springer.
- Press, W. H., Teukolsky, S. A., Vetterling, W. T. and Flannery, B. P. (2007), *Numerical Recipes*, 3rd ed., Cambridge, §10.2–10.3.

## 6. Open Questions

1. Bounds and constraints for the minimisers (L-BFGS-B, or a transformation of variables), which likelihoods with positive parameters would use.
2. A revised or interior-point method for large sparse linear programmes.
