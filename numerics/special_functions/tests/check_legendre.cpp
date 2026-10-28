// Checks of the Schmidt semi-normalised associated Legendre functions: their normalisation over
// [−1, 1] by Gauss–Kronrod, against the standard library's associated Legendre functions (C++17,
// [sf.cmath.assoc_legendre], which carry no Condon–Shortley phase) to degree 12, P̆₃ᵐ against its
// closed forms, the derivative in latitude against the identity of the degree and against central
// differences, and P̆ₙᵐ/cos φ with its limits at the poles.

#include <calx/quadrature/gauss_kronrod.hpp>
#include <calx/special_functions/legendre.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <numbers>

using namespace calx::special_functions;

namespace {
double factorial(int n) { return std::tgamma(n + 1.0); }
}  // namespace

int main() {
  calx::verification::Suite s;
  constexpr double pi = std::numbers::pi;

  // Normalisation: ∫₋₁¹ P̆ₙᵐ(μ)² dμ = 2 (2 − δₘ₀)/(2n + 1), from the orthogonality of the associated
  // Legendre functions (Abramowitz and Stegun 1972, 8.14.13), to degree 12.
  {
    double worst = 0;
    for (int n = 0; n <= 12; n++)
      for (int m = 0; m <= n; m++) {
        auto f = [&](double mu) { const double v = schmidt_legendre(12, std::asin(mu))(n, m); return v * v; };
        const double exact = 2.0 * (m == 0 ? 1 : 2) / (2 * n + 1);
        worst = std::max(worst, std::fabs(calx::quadrature::integrate(f, -1, 1, 1e-15, 1e-14).value / exact - 1));
      }
    s.check("Schmidt functions: their normalisation to degree 12, by Gauss-Kronrod", worst, 1e-12);
  }

  // An independent implementation: the standard library's Pₙ,ₘ times the normalisation.
  {
    double worst = 0;
    for (const double lat : {-1.4, -0.7, -0.05, 0.3, 0.9, 1.5}) {
      const SchmidtLegendre sf = schmidt_legendre(12, lat);
      for (int n = 0; n <= 12; n++)
        for (int m = 0; m <= n; m++) {
          const double norm = m ? std::sqrt(2 * factorial(n - m) / factorial(n + m)) : 1.0;
          worst = std::max(worst, std::fabs(sf(n, m) - norm * std::assoc_legendre(unsigned(n), unsigned(m), std::sin(lat))));
        }
    }
    s.check("Schmidt functions: against std::assoc_legendre to degree 12", worst, 1e-13);
  }

  // P₃,ₘ in closed form (as the WMM2025 report prints them in its eq. 6, from Pₙ,ₘ = νᵐ dᵐPₙ/dμᵐ),
  // times √(2 (3 − m)!/(3 + m)!).
  {
    double worst = 0;
    for (const double lat : {-1.2, -0.4, 0.1, 0.7, 1.4}) {
      const double sp = std::sin(lat), cp = std::cos(lat);
      const std::array<double, 4> p3{0.5 * sp * (5 * sp * sp - 3), -1.5 * cp * (1 - 5 * sp * sp), 15 * sp * (1 - sp * sp), 15 * cp * cp * cp};
      const SchmidtLegendre sf = schmidt_legendre(3, lat);
      for (int m = 0; m <= 3; m++) worst = std::max(worst, std::fabs(sf(3, m) - (m ? std::sqrt(2 * factorial(3 - m) / factorial(3 + m)) : 1.0) * p3[m]));
    }
    s.check("Schmidt functions: P̆₃ᵐ against their closed forms", worst, 1e-14);
  }

  // The derivative in latitude, against the identity in the degree
  // dP̆ₙᵐ/dφ = (n + 1) tan φ P̆ₙᵐ − √((n + 1)² − m²) sec φ P̆ₙ₊₁ᵐ (the WMM2025 report's eq. 16, from
  // Abramowitz and Stegun 1972, 8.5.4, rescaled), and against central differences.
  {
    double identity = 0, fd = 0;
    for (const double lat : {-1.3, -0.6, 0.004, 0.5, 1.2}) {
      const SchmidtLegendre sf = schmidt_legendre(13, lat), up = schmidt_legendre(12, lat + 1e-6), down = schmidt_legendre(12, lat - 1e-6);
      for (int n = 1; n <= 12; n++)
        for (int m = 0; m <= n; m++) {
          const std::size_t k = SchmidtLegendre::index(n, m);
          const double r = (n + 1) * std::tan(lat) * sf.p[k] - std::sqrt((n + 1.0) * (n + 1) - double(m) * m) / std::cos(lat) * sf(n + 1, m);
          identity = std::max(identity, std::fabs(sf.dp[k] - r));
          fd = std::max(fd, std::fabs(sf.dp[k] - (up.p[k] - down.p[k]) / 2e-6));
        }
    }
    s.check("Schmidt functions: dP̆ₙᵐ/dφ against the identity in the degree", identity, 1e-11);
    s.check("Schmidt functions: dP̆ₙᵐ/dφ against central differences", fd, 1e-8);
  }

  // P̆ₙᵐ/cos φ, and its limits at the poles: √(n(n + 1)/2) for m = 1 at the north pole, (−1)ⁿ⁺¹ times
  // that at the south pole, and 0 for m > 1.
  {
    double ratio = 0, limits = 0;
    for (const double lat : {-1.5, -0.2, 0.9, 1.55}) {
      const SchmidtLegendre sf = schmidt_legendre(12, lat);
      for (int n = 1; n <= 12; n++)
        for (int m = 1; m <= n; m++) ratio = std::max(ratio, std::fabs(sf.p_over_cos[SchmidtLegendre::index(n, m)] * std::cos(lat) - sf(n, m)));
    }
    for (const double pole : {1.0, -1.0}) {
      const SchmidtLegendre sf = schmidt_legendre(12, pole * pi / 2);
      for (int n = 1; n <= 12; n++)
        for (int m = 1; m <= n; m++) {
          const double expected = m == 1 ? (pole > 0 || n % 2 == 1 ? 1 : -1) * std::sqrt(n * (n + 1) / 2.0) : 0.0;
          limits = std::max(limits, std::fabs(sf.p_over_cos[SchmidtLegendre::index(n, m)] - expected));
        }
    }
    s.check("Schmidt functions: P̆ₙᵐ/cos φ times cos φ is P̆ₙᵐ", ratio, 1e-14);
    s.check("Schmidt functions: P̆ₙᵐ/cos φ at the poles, √(n(n + 1)/2) for m = 1 and 0 beyond", limits, 1e-12);
  }
  s.require("Schmidt functions: a negative degree is refused", [] {
    try { schmidt_legendre(-1, 0.0); } catch (const std::invalid_argument&) { return true; }
    return false;
  }());

  return s.done();
}
