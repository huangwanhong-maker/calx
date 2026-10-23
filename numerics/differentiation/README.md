# differentiation: Derivatives of Functions Written as Code

*Status: draft (2026-10-23): forward differentiation of order 1 and 2 by truncated Taylor arithmetic* ([decision 0005](../../development/decisions/0005-forward-differentiation.md)). Reverse mode, higher orders, and differentiation as a transformation of the IR (open decision 1) are planned.

## 1. Purpose

The derivatives of a function that is written once, as a C++ template over its number type, exact to rounding and without a step size. physis uses them for the curvature of metrics (`geometry`); optimisation, sensitivities and Jacobians of models will use them too.

## 2. Content

`jet.hpp`:

| Name | Content |
|---|---|
| `Taylor<N, K>`, `Jet<N>`, `Dual<N>` | A truncated Taylor polynomial in N variables of order K = 2 (value, gradient, Hessian) or K = 1 (value, gradient) |
| arithmetic | `+ - * /` with jets and doubles, comparisons of values, `value(x)` |
| elementary functions | `sqrt`, `exp`, `log`, `sin`, `cos`, `tan`, `sinh`, `cosh`, `tanh`, `atan`, `atan2`, `abs`, `pow` (real, integer and jet exponents), and `chain` for any other function given its first two derivatives |
| `variables<N, K>(x)` | The jets of the coordinates at x |
| `taylor(f, x)`, `gradient(f, x)`, `hessian(f, x)` | The jet, gradient (by dual numbers) and Hessian of a generic function at x |

A function to be differentiated calls the elementary functions unqualified after `using std::sin;` and the like, so the same code serves doubles and jets.

## 3. Checks

`tests/check_differentiation.cpp`: Rosenbrock's function against its closed-form gradient and Hessian; a composition of every elementary function against its closed-form gradient, and its Hessian against central differences of that gradient; atan2's gradient and its harmonicity; the mixed derivative of x^y; integer powers at negative arguments; Euler's theorem for a homogeneous function, at first and second order.

## 4. Dependencies

None.

## 5. References

- Griewank, A. and Walther, A. (2008), *Evaluating Derivatives: Principles and Techniques of Algorithmic Differentiation*, 2nd ed., SIAM, ch. 3 and 13.
- Rall, L. B. (1981), *Automatic Differentiation: Techniques and Applications*, Lecture Notes in Computer Science 120, Springer.
