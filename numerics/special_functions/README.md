# special_functions: Special Functions

*Status: draft (2026-10-24): the scaled complementary error function; the gamma family ([decision 0006](../../development/decisions/0006-lanes-task-graphs-and-the-genome.md)); 2026-10-28, the associated Legendre functions in Schmidt's semi-normalisation, in extended range ([decision 0010](../../development/decisions/0010-kernel-awareness-and-extended-range-legendre.md)).* Orthogonal polynomials, the fully normalised Legendre functions with extended range, and the Bessel, hypergeometric and elliptic functions are planned.

## 1. Content

`erf.hpp`: `erfcx(x) = exp(x²) erfc(x)`, from the C library's erfc below 10, the asymptotic series above, and reflection for negative arguments; relative accuracy about 2e-14 below 10 and 1e-16 above.

`gamma.hpp`: `log_gamma` (Lanczos, g = 7, nine coefficients; thread-safe, unlike the C library's lgamma, which writes the global signgam), `log_beta`, `log_choose`, `digamma` and `trigamma` (recurrence to x ≥ 10, then the asymptotic series); the regularised incomplete gamma functions `gamma_p` and `gamma_q` (series and Legendre's continued fraction by Lentz's method, Q computed directly where small) and `gamma_p_inverse` (Halley's method from the starting points of Numerical Recipes §6.2.1); `chi_square_q` and `chi_square_quantile`; the regularised incomplete beta function `beta_i` (its continued fraction, with the symmetry); `binomial_cdf`, `student_t_two_sided`, `f_q`.

`legendre.hpp`: `schmidt_legendre(degree, φ)`, the associated Legendre functions without the Condon–Shortley phase in Schmidt's semi-normalisation, $\breve P_n^m(\sin\varphi)$ for every $n \le$ degree and $m \le n$, with their derivatives in the latitude and $\breve P_n^m/\cos\varphi$, finite at the poles; by recurrences on the normalised functions themselves, without factorials. They are the functions of geomagnetism's expansions (the geophysics domain's WMM2025 uses them); geodesy's fully normalised functions are $\sqrt{2n+1}$ times them. They were written in the geophysics domain and moved here, since they name no object in nature. Each order's column is kept in extended range (Fukushima 2012), a double and a binary exponent of its own, renormalised past $2^{480}$, so that the sectoral values $\propto\cos^m\varphi$, which underflow below the smallest double near the poles at high degree, are not lost. The functions hold at every latitude, the poles included, to degrees far beyond the crustal models' (449 for the Moon); the addition theorem $\sum_m(\breve P_n^m)^2=1$ holds to $n^2\varepsilon$.

## 2. Checks

Agreement with exp(x²) erfc(x) on [0, 5]; the meeting of the two branches at 10; the asymptotic series at 1000; the reflection; erfcx(0) = 1. The gamma family (`check_calx_gamma`): log Γ at the integers and at 1/2 and against the C library on [0.01, 300]; ψ(1), ψ(1/2), ψ′(1), ψ′(1/2) and the recurrences; P(1, x) = 1 − e^{−x}, P(1/2, x) = erf √x, P + Q = 1, Q small and exact; the χ² points of the tables; P inverted from a = 0.05 to 900; I_x(a, 1), I_x(1, b), the symmetry, I_{1/2}(a, a) = 1/2; the binomial distribution against the sum of its terms; Student's t against the closed form for even degrees of freedom (Abramowitz and Stegun 26.7.3). The Legendre functions (`check_calx_legendre`): their normalisation over $[-1, 1]$ by Gauss–Kronrod to degree 12; against the standard library's `std::assoc_legendre` (C++17), an independent implementation; $\breve P_3^m$ in closed form; the derivative against the identity in the degree and against central differences; $\breve P_n^m/\cos\varphi$ and its limits at the poles; the addition theorem to $n^2\varepsilon$ at degrees 449, 1000 and 2700, from the equator to $2\times10^{-6}$ rad of each pole; the derivatives by the identity and by central differences at degree 449 near the poles.

## 3. Dependencies

None; the checks of the Legendre functions use `quadrature`.

## 4. References

- Abramowitz, M. and Stegun, I. A. (1964), *Handbook of Mathematical Functions*, National Bureau of Standards, 6.3.18, 6.4.12, 7.1.23, 26.7.3, Table 26.8.
- Lanczos, C. (1964), A precision approximation of the gamma function, *SIAM J. Numer. Anal. B* 1, 86–96.
- Press, W. H., Teukolsky, S. A., Vetterling, W. T. and Flannery, B. P. (2007), *Numerical Recipes*, 3rd ed., §6.1, §6.2, §6.4.
- Chulliat, A. et al. (2025), *The US/UK World Magnetic Model for 2025–2030: Technical Report*, NOAA National Centers for Environmental Information, §1.2 (eqs. 5, 6, 16).
- Fukushima, T. (2012), Numerical computation of spherical harmonics of arbitrary degree and order by extending exponent of floating point numbers, *Journal of Geodesy* 86, 271–285.
- Holmes, S. A. and Featherstone, W. E. (2002), A unified approach to the Clenshaw summation and the recursive computation of very high degree and order normalised associated Legendre functions, *Journal of Geodesy* 76, 279–299.
- DiDonato, A. R. and Morris, A. H. (1986), Computation of the incomplete gamma function ratios and their inverse, *ACM TOMS* 12, 377–393.
