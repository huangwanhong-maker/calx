# special_functions: Special Functions

*Status: draft (2026-10-23): the scaled complementary error function.* Gamma functions, orthogonal polynomials, and the Legendre, Bessel, hypergeometric and elliptic functions are planned.

## 1. Content

`erf.hpp`: `erfcx(x) = exp(x²) erfc(x)`, from the C library's erfc below 10, the asymptotic series above, and reflection for negative arguments; relative accuracy about 2e-14 below 10 and 1e-16 above.

## 2. Checks

Agreement with exp(x²) erfc(x) on [0, 5]; the meeting of the two branches at 10; the asymptotic series at 1000; the reflection; erfcx(0) = 1.

## 3. Dependencies

None.

## 4. References

- Abramowitz, M. and Stegun, I. A. (1964), *Handbook of Mathematical Functions*, National Bureau of Standards, 7.1.23.
