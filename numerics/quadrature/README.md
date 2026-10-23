# quadrature: Numerical Integration

*Status: draft (2026-10-23): adaptive Gauss–Kronrod on an interval.* Integration on the sphere and in many dimensions is planned.

## 1. Content

`gauss_kronrod.hpp`: the pair (7, 15) on an interval, and `integrate`, which bisects the subinterval of largest estimated error until the summed estimate meets its absolute or relative tolerance (QUADPACK's QAG without extrapolation; nodes and weights from qk15).

## 2. Checks

The degrees of exactness of the pair; integrals with closed forms, smooth, peaked, and with a singular derivative.

## 3. Dependencies

None.

## 4. References

- Kronrod, A. S. (1965), *Nodes and Weights of Quadrature Formulas*, Consultants Bureau.
- Piessens, R., de Doncker-Kapenga, E., Überhuber, C. W. and Kahaner, D. K. (1983), *QUADPACK*, Springer.
