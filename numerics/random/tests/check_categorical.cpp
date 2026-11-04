// Checks of calx/random/categorical.hpp: categorical variates by inversion and by Walker's alias table.
//
// - The alias table's law, read back from its columns, equals the normalised weights to rounding,
//   for even, skewed and sparse weights.
// - Pearson's χ² of 10⁶ draws by each method against the probabilities lies within the χ²
//   distribution's bulk; an index of weight zero is never drawn.
// - One weight alone is always drawn; a draw depends only on the stream (the same seed, the same
//   draws).
// - Refusals: negative, infinite or not-a-number weights, and no positive weight.

#include <calx/random/categorical.hpp>

#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

double chi2(const std::vector<double>& count, const std::vector<double>& w, double n) {
  double total = 0, x = 0;
  for (const double v : w) total += v;
  for (std::size_t j = 0; j < w.size(); j++) {
    if (w[j] == 0) continue;
    const double e = n * w[j] / total;
    x += (count[j] - e) * (count[j] - e) / e;
  }
  return x;
}

}  // namespace

int main() {
  calx::verification::Suite s;
  using calx::random::AliasTable;
  using calx::random::categorical;
  using calx::random::Stream;

  // ---- the alias table's law --------------------------------------------------------------------------------------
  const std::vector<std::vector<double>> laws = {
      {1, 1, 1, 1}, {0.7, 0.2, 0.05, 0.05}, {0, 3, 0, 1, 0, 0, 6, 0}, {1e-9, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}};
  double worst = 0;
  for (const auto& w : laws) {
    double total = 0;
    for (const double v : w) total += v;
    const auto q = AliasTable(w).law();
    for (std::size_t j = 0; j < w.size(); j++) worst = std::max(worst, std::fabs(q[j] - w[j] / total));
  }
  s.check("the alias table's law equals the normalised weights", worst, 1e-14);

  // ---- frequencies --------------------------------------------------------------------------------------------------
  // Eight categories, two of weight zero: χ² on 5 degrees of freedom, whose 0.1% and 99.9% quantiles are 0.21 and 20.5.
  const std::vector<double> w{0.5, 0, 2, 1, 0, 3, 0.25, 1.25};
  const double n = 1e6;
  Stream rng(20261104, 0x63617467u, 0, 0);  // "catg"
  std::vector<double> by_inversion(w.size(), 0.0), by_alias(w.size(), 0.0);
  const AliasTable table(w);
  for (int k = 0; k < static_cast<int>(n); k++) by_inversion[categorical(rng, w)] += 1, by_alias[table(rng)] += 1;
  const double x1 = chi2(by_inversion, w, n), x2 = chi2(by_alias, w, n);
  s.require("by inversion: Pearson's χ² of 10⁶ draws lies in (0.21, 20.5)", x1 > 0.21 && x1 < 20.5, std::to_string(x1));
  s.require("by the alias table: Pearson's χ² of 10⁶ draws lies in (0.21, 20.5)", x2 > 0.21 && x2 < 20.5, std::to_string(x2));
  s.require("an index of weight zero is never drawn, by either method",
            by_inversion[1] == 0 && by_inversion[4] == 0 && by_alias[1] == 0 && by_alias[4] == 0);

  // ---- one weight; the same stream, the same draws ------------------------------------------------------------------
  {
    const std::vector<double> one{0, 0, 5, 0};
    const AliasTable t(one);
    bool always = true;
    for (int k = 0; k < 1000; k++) always = always && categorical(rng, one) == 2 && t(rng) == 2;
    s.require("one positive weight is always drawn", always);
    Stream a(9, 1, 2, 3), b(9, 1, 2, 3);
    bool same = true;
    for (int k = 0; k < 1000; k++) same = same && categorical(a, w) == categorical(b, w) && table(a) == table(b);
    s.require("the same seed and identity give the same draws", same);
  }

  // ---- refusals -----------------------------------------------------------------------------------------------------
  int refused = 0;
  for (const auto& bad : std::vector<std::vector<double>>{{1, -1}, {0, 0}, {1, INFINITY}, {NAN, 1}, {}}) {
    try {
      categorical(rng, bad);
    } catch (const std::invalid_argument&) {
      refused++;
    }
    try {
      AliasTable t(bad);
    } catch (const std::invalid_argument&) {
      refused++;
    }
  }
  s.require("refusals: a negative, infinite or undefined weight, and no positive weight", refused == 10, std::to_string(refused));
  return s.done();
}
