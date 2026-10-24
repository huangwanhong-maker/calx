# special_functions: Special Functions

*Status: draft (2026-10-24): the scaled complementary error function; the gamma family ([decision 0006](../../development/decisions/0006-lanes-task-graphs-and-the-genome.md)).* Orthogonal polynomials, and the Legendre, Bessel, hypergeometric and elliptic functions are planned.

## 1. Content

`erf.hpp`: `erfcx(x) = exp(x²) erfc(x)`, from the C library's erfc below 10, the asymptotic series above, and reflection for negative arguments; relative accuracy about 2e-14 below 10 and 1e-16 above.

`gamma.hpp`: `log_gamma` (Lanczos, g = 7, nine coefficients; thread-safe, unlike the C library's lgamma, which writes the global signgam), `log_beta`, `log_choose`, `digamma` and `trigamma` (recurrence to x ≥ 10, then the asymptotic series); the regularised incomplete gamma functions `gamma_p` and `gamma_q` (series and Legendre's continued fraction by Lentz's method, Q computed directly where small) and `gamma_p_inverse` (Halley's method from the starting points of Numerical Recipes §6.2.1); `chi_square_q` and `chi_square_quantile`; the regularised incomplete beta function `beta_i` (its continued fraction, with the symmetry); `binomial_cdf`, `student_t_two_sided`, `f_q`.

## 2. Checks

Agreement with exp(x²) erfc(x) on [0, 5]; the meeting of the two branches at 10; the asymptotic series at 1000; the reflection; erfcx(0) = 1. The gamma family (`check_calx_gamma`): log Γ at the integers and at 1/2 and against the C library on [0.01, 300]; ψ(1), ψ(1/2), ψ′(1), ψ′(1/2) and the recurrences; P(1, x) = 1 − e^{−x}, P(1/2, x) = erf √x, P + Q = 1, Q small and exact; the χ² points of the tables; P inverted from a = 0.05 to 900; I_x(a, 1), I_x(1, b), the symmetry, I_{1/2}(a, a) = 1/2; the binomial distribution against the sum of its terms; Student's t against the closed form for even degrees of freedom (Abramowitz and Stegun 26.7.3).

## 3. Dependencies

None.

## 4. References

- Abramowitz, M. and Stegun, I. A. (1964), *Handbook of Mathematical Functions*, National Bureau of Standards, 6.3.18, 6.4.12, 7.1.23, 26.7.3, Table 26.8.
- Lanczos, C. (1964), A precision approximation of the gamma function, *SIAM J. Numer. Anal. B* 1, 86–96.
- Press, W. H., Teukolsky, S. A., Vetterling, W. T. and Flannery, B. P. (2007), *Numerical Recipes*, 3rd ed., §6.1, §6.2, §6.4.
- DiDonato, A. R. and Morris, A. H. (1986), Computation of the incomplete gamma function ratios and their inverse, *ACM TOMS* 12, 377–393.
