// Associated Legendre functions in Schmidt's semi-normalisation, of sin φ for a latitude φ, with
// their derivatives in φ: the functions of the spherical-harmonic expansions of potential fields
// (geomagnetism's WMM and IGRF, the crustal fields of the Moon and Mars; with a factor √(2n + 1),
// geodesy's fully normalised functions).
//
// For a latitude φ (radians; the colatitude θ = π/2 − φ, so that sin φ = cos θ and d/dθ = −d/dφ)
// write μ = sin φ and ν = cos φ. The functions are
//   P̆ₙᵐ(μ) = √(2 (n − m)!/(n + m)!) Pₙ,ₘ(μ) for m > 0,   P̆ₙ⁰(μ) = Pₙ,₀(μ),
// with Pₙ,ₘ = ν^m dᵐPₙ/dμᵐ the associated Legendre functions without the Condon–Shortley phase
// (Pₙ,ₘ = (−1)ᵐ Pₙᵐ of Abramowitz and Stegun 1972, §8.6; the convention of geomagnetism, as the
// WMM2025 report has it, Chulliat et al. 2025, §1.2, eq. 5). Their mean square over the sphere
// with cos mλ is 1/(2n + 1) for every m, and by the addition theorem Σₘ P̆ₙᵐ(μ)² = 1 at every μ:
// the coefficients of one degree are on one scale, and no function exceeds 1 in magnitude.
//
// They are computed by recurrences on P̆ themselves, which need no factorials:
//   P̆₀⁰ = 1,  P̆₁¹ = ν,  P̆ₘᵐ = √((2m − 1)/(2m)) ν P̆ₘ₋₁ᵐ⁻¹ (m ≥ 2),
//   P̆ₙᵐ = [(2n − 1) μ P̆ₙ₋₁ᵐ − √((n − 1)² − m²) P̆ₙ₋₂ᵐ] / √(n² − m²)   (n > m),
// the second the three-term recurrence in the degree (Abramowitz and Stegun 1972, 8.5.3) rescaled by
// the normalisation. Differentiating each step in φ (dμ/dφ = ν, dν/dφ = −μ) gives the derivatives
// with no division by ν, so they hold at the poles. For m ≥ 1 every P̆ₙᵐ carries the factor νᵐ, and
// the same recurrences started from P̆₁¹/ν = 1 give P̆ₙᵐ/ν, which the east component of a potential
// field needs and which stays finite at the poles: its limit there is √(n(n + 1)/2) for m = 1
// (times (−1)ⁿ⁺¹ at the south pole) and 0 for m > 1.
//
// Extended range. The sectoral seed P̆ₘᵐ carries νᵐ, which underflows the doubles near the poles for
// orders of a few hundred (ν = 0.2 and m = 449 give 10⁻³¹⁴), though the functions of high degree it
// seeds are not negligible there. Each order's column is therefore computed in scaled arithmetic, the
// X-numbers of Fukushima (2012): its values are held as doubles times 2^e with one integer exponent
// e for the column, renormalised when they grow, and multiplied out by 2^e only at the end, where a
// value below the smallest double is truly negligible beside the others of its degree (whose squares
// sum to 1). The functions then hold at every degree and latitude to the rounding of the recurrence,
// which grows as n²ε at most (some 10⁻¹⁰ at degree 2700); the checks hold them to degree 2700, at the
// equator and within 2 × 10⁻⁶ rad of the poles.

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

  // The sectoral functions, scaled: (sp, sd, sq) · 2^se are P̆ₘᵐ, dP̆ₘᵐ/dφ and P̆ₘᵐ/ν, with sp kept in
  // [1/2, 1) by moving powers of two into se.
  double sp = 1, sd = 0, sq = 0;
  long se = 0;
  constexpr int big = 480;  // renormalise a column once its values pass 2^480
  const double up = std::ldexp(1.0, big), down = std::ldexp(1.0, -big);
  for (int m = 0; m <= degree; m++) {
    if (m >= 1) {
      const double k = m == 1 ? 1.0 : std::sqrt((2.0 * m - 1) / (2.0 * m));
      const double p0 = sp, d0 = sd;
      sp = k * nu * p0;
      sd = k * (nu * d0 - mu * p0);
      sq = m == 1 ? 1.0 : k * nu * sq;  // P̆₁¹/ν = 1; the scale is still 2⁰ at m = 1
      if (sp != 0) {
        int e = 0;
        std::frexp(sp, &e);
        sp = std::ldexp(sp, -e);
        sd = std::ldexp(sd, -e);
        sq = std::ldexp(sq, -e);
        se += e;
      }
    }
    // The column of order m, on its own scale e: y = P̆/2^e.
    long e = se;
    double y1 = sp, d1 = sd, q1 = sq;      // degree m
    double y2 = 0, d2 = 0, q2 = 0;         // degree m − 1 (absent)
    auto store = [&](int n, double y, double d, double q) {
      const std::size_t k = at(n, m);
      const int ei = e < -100000 ? -100000 : e > 100000 ? 100000 : int(e);
      s.p[k] = std::ldexp(y, ei);
      s.dp[k] = std::ldexp(d, ei);
      if (m >= 1) s.p_over_cos[k] = std::ldexp(q, ei);
    };
    store(m, y1, d1, q1);
    for (int n = m + 1; n <= degree; n++) {
      const double root = std::sqrt(double(n) * n - double(m) * m);
      const double c1 = (2.0 * n - 1) / root;
      const double c2 = std::sqrt(double(n - 1) * (n - 1) - double(m) * m) / root;  // 0 when n = m + 1
      const double y = c1 * mu * y1 - c2 * y2;
      const double d = c1 * (nu * y1 + mu * d1) - c2 * d2;
      const double q = m >= 1 ? c1 * mu * q1 - c2 * q2 : 0.0;
      y2 = y1, d2 = d1, q2 = q1;
      y1 = y, d1 = d, q1 = q;
      if (std::fabs(y1) > up || std::fabs(d1) > up || std::fabs(q1) > up) {
        y1 *= down, d1 *= down, q1 *= down, y2 *= down, d2 *= down, q2 *= down;
        e += big;
      }
      store(n, y1, d1, q1);
    }
  }
  return s;
}

}  // namespace calx::special_functions
