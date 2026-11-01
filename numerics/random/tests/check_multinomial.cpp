// Checks of calx/random/multinomial.hpp: multinomial variates by conditional binomials.
//
// - Every draw sums to n, with no count in a category of probability 0.
// - Means n p_j and variances n p_j (1 − p_j) of 20 000 draws, within four standard errors; the
//   covariance of two categories −n p_i p_j.
// - Pearson's χ² of one draw of a large n against its expectation within the χ² distribution's bulk.
// - Probabilities that do not sum to 1 are normalised; refusals (negative n, no positive probability).

#include <calx/random/multinomial.hpp>

#include <calx/verification/check.hpp>

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

int main() {
  calx::verification::Suite s;
  calx::random::Stream rng(7, 1, 2, 3);
  const std::vector<double> p{0.1, 0.0, 0.25, 0.05, 0.6};  // sums to 1
  const std::int64_t n = 50;
  const int draws = 20000;
  std::vector<double> mean(p.size(), 0.0), sq(p.size(), 0.0);
  double cross = 0;
  bool sums = true, zero = true;
  for (int d = 0; d < draws; d++) {
    const auto k = calx::random::multinomial(rng, n, p);
    std::int64_t t = 0;
    for (std::size_t j = 0; j < p.size(); j++) {
      t += k[j];
      mean[j] += static_cast<double>(k[j]) / draws;
      sq[j] += static_cast<double>(k[j]) * static_cast<double>(k[j]) / draws;
    }
    cross += static_cast<double>(k[0]) * static_cast<double>(k[4]) / draws;
    sums = sums && t == n;
    zero = zero && k[1] == 0;
  }
  s.require("every draw sums to n", sums);
  s.require("a category of probability 0 is never drawn", zero);
  double worst_mean = 0, worst_var = 0;
  for (std::size_t j = 0; j < p.size(); j++) {
    if (p[j] == 0) continue;
    const double v = n * p[j] * (1 - p[j]);
    worst_mean = std::max(worst_mean, std::fabs(mean[j] - n * p[j]) / std::sqrt(v / draws));
    // The sample variance's standard error, about v √(2/draws) for counts this large.
    worst_var = std::max(worst_var, std::fabs(sq[j] - mean[j] * mean[j] - v) / (v * std::sqrt(2.0 / draws) * 1.5));
  }
  s.check("means n p_j (in standard errors)", worst_mean, 4);
  s.check("variances n p_j (1 − p_j) (in standard errors)", worst_var, 4);
  const double cov = cross - mean[0] * mean[4], expect = -n * p[0] * p[4];
  s.check("covariance of two categories −n p_i p_j (relative)", cov / expect - 1, 0.1);

  // Pearson's χ² of one large draw: 9 degrees of freedom; its 0.1% and 99.9% quantiles are 1.15 and 27.9.
  std::vector<double> q(10);
  for (std::size_t j = 0; j < 10; j++) q[j] = 1.0 + static_cast<double>(j);  // unnormalised
  double total = 0;
  for (const double v : q) total += v;
  const std::int64_t N = 1000000;
  const auto k = calx::random::multinomial(rng, N, q);
  double chi2 = 0;
  for (std::size_t j = 0; j < 10; j++) {
    const double e = N * q[j] / total;
    chi2 += (static_cast<double>(k[j]) - e) * (static_cast<double>(k[j]) - e) / e;
  }
  s.require("Pearson's χ² of a draw of 10⁶ with unnormalised probabilities lies in (1.15, 27.9)", chi2 > 1.15 && chi2 < 27.9, std::to_string(chi2));
  bool refused = false, refused_zero = false;
  try {
    calx::random::multinomial(rng, -1, p);
  } catch (const std::invalid_argument&) {
    refused = true;
  }
  try {
    const std::vector<double> none{0, 0};
    calx::random::multinomial(rng, 3, none);
  } catch (const std::invalid_argument&) {
    refused_zero = true;
  }
  s.require("refusals: n < 0, and no positive probability", refused && refused_zero);
  return s.done();
}
