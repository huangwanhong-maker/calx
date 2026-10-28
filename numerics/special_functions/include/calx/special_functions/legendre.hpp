// Associated Legendre functions in Schmidt's semi-normalisation, of sin φ for a latitude φ, with
// their derivatives in φ: the functions of the spherical-harmonic expansions of potential fields
// (geomagnetism's WMM and IGRF; with a factor √(2n + 1), geodesy's fully normalised functions).
//
// For a latitude φ (radians; the colatitude θ = π/2 − φ, so that sin φ = cos θ and d/dθ = −d/dφ)
// write μ = sin φ and ν = cos φ. The functions are
//   P̆ₙᵐ(μ) = √(2 (n − m)!/(n + m)!) Pₙ,ₘ(μ) for m > 0,   P̆ₙ⁰(μ) = Pₙ,₀(μ),
// with Pₙ,ₘ = ν^m dᵐPₙ/dμᵐ the associated Legendre functions without the Condon–Shortley phase
// (Pₙ,ₘ = (−1)ᵐ Pₙᵐ of Abramowitz and Stegun 1972, §8.6; the convention of geomagnetism, as the
// WMM2025 report has it, Chulliat et al. 2025, §1.2, eq. 5). Their mean square over the sphere
// with cos mλ is 1/(2n + 1) for every m: the coefficients of one degree are on one scale.
//
// They are computed by recurrences on P̆ themselves, which need no factorials and stay within a
// few units of rounding to the degrees of the geomagnetic and gravity models:
//   P̆₀⁰ = 1,  P̆₁¹ = ν,  P̆ₘᵐ = √((2m − 1)/(2m)) ν P̆ₘ₋₁ᵐ⁻¹ (m ≥ 2),
//   P̆ₙᵐ = [(2n − 1) μ P̆ₙ₋₁ᵐ − √((n − 1)² − m²) P̆ₙ₋₂ᵐ] / √(n² − m²)   (n > m),
// the second the three-term recurrence in the degree (Abramowitz and Stegun 1972, 8.5.3) rescaled by
// the normalisation. Differentiating each step in φ (dμ/dφ = ν, dν/dφ = −μ) gives the derivatives
// with no division by ν, so they hold at the poles. For m ≥ 1 every P̆ₙᵐ carries the factor νᵐ, and
// the same recurrences started from P̆₁¹/ν = 1 give P̆ₙᵐ/ν, which the east component of a potential
// field needs and which stays finite at the poles: its limit there is √(n(n + 1)/2) for m = 1
// (times (−1)ⁿ⁺¹ at the south pole) and 0 for m > 1.
//
// Caveat. The sectoral recurrence multiplies by ν at each order, so P̆ₘᵐ underflows near the poles
// for orders of several hundred; the extended-range arithmetic of Holmes and Featherstone (2002)
// is needed beyond degree some 1800 everywhere and is not written. The checks hold the functions
// to degree 13.

#pragma once

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace calx::special_functions {

// The functions of every degree n ≤ degree and order m ≤ n at one latitude, stored by (n, m).
struct SchmidtLegendre {
  int degree = 0;
  std::vector<double> p;           // P̆ₙᵐ(sin φ)
  std::vector<double> dp;          // dP̆ₙᵐ(sin φ)/dφ
  std::vector<double> p_over_cos;  // P̆ₙᵐ(sin φ)/cos φ for m ≥ 1 (finite at the poles); 0 for m = 0

  // The position of (n, m), 0 ≤ m ≤ n, in the arrays: n(n + 1)/2 + m.
  static constexpr std::size_t index(int n, int m) { return std::size_t(n) * std::size_t(n + 1) / 2 + std::size_t(m); }
  double operator()(int n, int m) const { return p[index(n, m)]; }
};

// P̆ₙᵐ(sin φ), its derivative in φ and P̆ₙᵐ/cos φ for every n ≤ degree and m ≤ n, at the latitude φ
// (radians).
inline SchmidtLegendre schmidt_legendre(int degree, double latitude_rad) {
  if (degree < 0) throw std::invalid_argument("calx::special_functions::schmidt_legendre: a negative degree");
  SchmidtLegendre s;
  s.degree = degree;
  const std::size_t size = SchmidtLegendre::index(degree, degree) + 1;
  s.p.assign(size, 0.0);
  s.dp.assign(size, 0.0);
  s.p_over_cos.assign(size, 0.0);
  const double mu = std::sin(latitude_rad), nu = std::cos(latitude_rad);
  auto at = [](int n, int m) { return SchmidtLegendre::index(n, m); };
  s.p[0] = 1;
  for (int m = 0; m <= degree; m++) {
    if (m >= 1) {
      // The sectoral function from the one below it.
      const double k = m == 1 ? 1.0 : std::sqrt((2.0 * m - 1) / (2.0 * m));
      const std::size_t below = at(m - 1, m - 1), here = at(m, m);
      s.p[here] = k * nu * s.p[below];
      s.dp[here] = k * (nu * s.dp[below] - mu * s.p[below]);
      s.p_over_cos[here] = m == 1 ? 1.0 : k * nu * s.p_over_cos[below];
    }
    // Up the degrees at this order.
    for (int n = m + 1; n <= degree; n++) {
      const double root = std::sqrt(double(n) * n - double(m) * m);
      const double c1 = (2.0 * n - 1) / root;
      const double c2 = std::sqrt(double(n - 1) * (n - 1) - double(m) * m) / root;  // 0 when n = m + 1
      const std::size_t here = at(n, m), one = at(n - 1, m);
      const double p2 = n - 2 >= m ? s.p[at(n - 2, m)] : 0.0;
      const double dp2 = n - 2 >= m ? s.dp[at(n - 2, m)] : 0.0;
      const double q2 = n - 2 >= m ? s.p_over_cos[at(n - 2, m)] : 0.0;
      s.p[here] = c1 * mu * s.p[one] - c2 * p2;
      s.dp[here] = c1 * (nu * s.p[one] + mu * s.dp[one]) - c2 * dp2;
      if (m >= 1) s.p_over_cos[here] = c1 * mu * s.p_over_cos[one] - c2 * q2;
    }
  }
  return s;
}

}  // namespace calx::special_functions
