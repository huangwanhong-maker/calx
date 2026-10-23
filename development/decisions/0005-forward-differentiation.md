# 0005. Forward Differentiation by Truncated Taylor Arithmetic

- **Status:** Accepted, provisionally (open decision 1, automatic differentiation as a transformation of the IR, stays open)
- **Date:** 2026-10-23
- **Decided by:** the user asked for spacetime models and the infrastructure of field theory in physis (its decision 0005). Their curvature needs the first and second derivatives of a metric written as code. Claude placed the method here and wrote this record; the composition is The Science Infrastructures' [ADR 0039](https://github.com/huangwanhong-maker/The-Science-Infrastructures/blob/main/development/decisions/0039-complexes-symmetry-and-fields.md).

## 1. Context

The Riemann tensor of a metric g(x) needs ∂g and ∂∂g. Finite differences lose half the digits for the first derivatives and more for the second; symbolic differentiation needs a computer algebra system that the infrastructure does not have. A metric written once as a C++ template, and evaluated on numbers that carry their own derivatives, gives every derivative to rounding (Griewank and Walther 2008, ch. 3 and 13). Differentiation is a numerical method on functions, which is calx's (the map of the mathematics, MSC 65D25).

## 2. Decision

1. **A new module, `numerics/differentiation`,** header-only, status `draft`: `Jet<N>`, a truncated multivariate Taylor polynomial of order 2 in N variables (value, gradient and Hessian), with its arithmetic and the elementary functions (sqrt, exp, log, sin, cos, tan, sinh, cosh, tanh, atan, pow); `Dual<N>`, its first-order part; `variables<N>(x)`, the jets of the coordinates; and `gradient`, `hessian` of a function written as a template.
2. **Operator overloading, not a transformation of the IR.** It serves functions written directly in C++, such as metrics and potentials; open decision 1 (automatic differentiation as a transformation of the IR, for kernels on devices and reverse mode for optimisation and learning) stays open, and the two will share their checks.
3. **Cost.** A product of two jets costs O(N²) operations; for the metrics of spacetime (N = 4) this is a few hundred operations per product, against a symbolic expression that would have to be derived for each metric.

## 3. Alternatives Considered

- **Finite differences with Richardson extrapolation:** about 10⁻⁸ relative in second derivatives at best, too coarse to check that a vacuum solution is Ricci-flat.
- **Derivatives written by hand for each metric:** error-prone, and against the aim that a metric be stated once.
- **A borrowed library** (autodiff, CppAD): licence-compatible ones exist, but a second-order jet is a page of code, and writing it is how the method is understood (the map, §10.3).

## 4. Consequences

- physis `geometry` computes Christoffel symbols and curvature from any metric written as a template.
- calx's README and CLAUDE.md list the module; its checks compare jets with closed-form derivatives.
