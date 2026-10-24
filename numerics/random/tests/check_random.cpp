// Checks of calx/random: Philox4x32-10 against Random123's known-answer vectors, and the
// moments of each distribution against their exact values, as z-scores of the sample mean
// (or variance) over a million draws; a correct sampler exceeds |z| = 5 with probability
// below 1e-6, and the seeds are fixed, so the checks are deterministic.

#include <calx/random/philox.hpp>
#include <calx/random/poisson.hpp>
#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

using calx::random::Philox4x32;
using calx::random::Stream;
using calx::random::PoissonTable;

namespace {

// Welford's running mean and variance.
struct Moments {
  double n = 0, mean = 0, m2 = 0;
  void add(double x) {
    n++;
    const double d = x - mean;
    mean += d / n;
    m2 += d * (x - mean);
  }
  double variance() const { return m2 / (n - 1); }
};

}  // namespace

int main() {
  calx::verification::Suite s;

  // ---- Philox4x32-10: Random123's kat_vectors (DEShawResearch/random123, tests/) --------------
  struct Kat { Philox4x32::Counter ctr; Philox4x32::Key key; Philox4x32::Counter out; };
  const Kat kats[] = {
      {{0, 0, 0, 0}, {0, 0}, {0x6627e8d5, 0xe169c58d, 0xbc57ac4c, 0x9b00dbd8}},
      {{0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff}, {0xffffffff, 0xffffffff}, {0x408f276d, 0x41c83b0e, 0xa20bc7c6, 0x6d5451fd}},
      {{0x243f6a88, 0x85a308d3, 0x13198a2e, 0x03707344}, {0xa4093822, 0x299f31d0}, {0xd16cfe09, 0x94fdcceb, 0x5001e420, 0x24126ea1}},
  };
  int matched = 0;
  for (const auto& k : kats) matched += Philox4x32::apply(k.ctr, k.key) == k.out;
  s.require("random: Philox4x32-10 gives Random123's three known answers", matched == 3, std::to_string(matched) + " of 3");
  static_assert(Philox4x32::apply({0, 0, 0, 0}, {0, 0})[0] == 0x6627e8d5u, "Philox is evaluated at compile time");

  // ---- streams ----------------------------------------------------------------------------------
  {
    Stream a(42, 1, 2, 3), b(42, 1, 2, 3), c(42, 1, 2, 4), d(43, 1, 2, 3);
    bool same = true, differ_c = false, differ_d = false;
    for (int i = 0; i < 1000; i++) {
      const auto x = a.next32(), y = b.next32(), z = c.next32(), w = d.next32();
      same &= x == y;
      differ_c |= x != z;
      differ_d |= x != w;
    }
    s.require("random: a stream is a function of its seed and identity", same && differ_c && differ_d);
    Stream e(7, 0, 0, 0), f(7, 0, 0, 0);
    for (int i = 0; i < 40; i++) e.next32();
    f.seek(10);
    s.require("random: seek(10) goes to the eleventh block of four draws", e.next32() == f.next32());
  }

  const int N = 1000000;
  // ---- uniform ----------------------------------------------------------------------------------
  {
    Stream r(2026, 10, 0, 0);
    Moments m;
    double lo = 1, hi = 0;
    for (int i = 0; i < N; i++) {
      const double u = r.uniform();
      m.add(u);
      lo = std::min(lo, u);
      hi = std::max(hi, u);
    }
    s.check("random: uniform, z of the mean (1/2)", (m.mean - 0.5) / std::sqrt(1.0 / 12 / N), 5);
    s.check("random: uniform, z of the variance (1/12)", (m.variance() - 1.0 / 12) / std::sqrt(1.0 / 180 / N), 5);
    s.require("random: uniform draws lie in the open interval (0, 1)", lo > 0 && hi < 1);
  }
  // ---- bounded integers: a chi-square over 10 cells ---------------------------------------------
  {
    Stream r(2026, 11, 0, 0);
    std::vector<double> cells(10, 0);
    for (int i = 0; i < N; i++) cells[r.below(10)]++;
    double chi2 = 0;
    for (double c : cells) chi2 += (c - N / 10.0) * (c - N / 10.0) / (N / 10.0);
    // 9 degrees of freedom: P(chi2 > 33.7) = 1e-4.
    s.check("random: below(10), chi-square over its ten values (9 degrees of freedom)", chi2, 33.7);
  }
  // ---- normal -----------------------------------------------------------------------------------
  {
    Stream r(2026, 12, 0, 0);
    Moments m;
    double m4 = 0;
    for (int i = 0; i < N; i++) {
      const double x = r.normal();
      m.add(x);
      m4 += x * x * x * x;
    }
    m4 /= N;
    s.check("random: normal, z of the mean (0)", m.mean / std::sqrt(1.0 / N), 5);
    s.check("random: normal, z of the variance (1)", (m.variance() - 1) / std::sqrt(2.0 / N), 5);
    s.check("random: normal, z of the fourth moment (3)", (m4 - 3) / std::sqrt(96.0 / N), 5);
  }
  // ---- exponential and geometric ------------------------------------------------------------------
  {
    Stream r(2026, 13, 0, 0);
    Moments e, g;
    const double p = 0.2;
    for (int i = 0; i < N; i++) {
      e.add(r.exponential());
      g.add(static_cast<double>(r.geometric(p)));
    }
    s.check("random: exponential, z of the mean (1)", (e.mean - 1) / std::sqrt(1.0 / N), 5);
    const double gm = (1 - p) / p, gv = (1 - p) / (p * p);
    s.check("random: geometric (p = 0.2), z of the mean (4)", (g.mean - gm) / std::sqrt(gv / N), 5);
  }
  // ---- Poisson, by multiplication and by PTRS ---------------------------------------------------
  for (double lam : {0.3, 3.0, 9.9, 10.0, 30.0, 300.0}) {
    Stream r(2026, 14, static_cast<std::uint32_t>(lam * 10), 0);
    Moments m;
    for (int i = 0; i < N; i++) m.add(static_cast<double>(r.poisson(lam)));
    const std::string name = "random: Poisson (mean " + calx::verification::shortest(lam) + ")";
    s.check(name + ", z of the mean", (m.mean - lam) / std::sqrt(lam / N), 5);
    s.check(name + ", z of the variance", (m.variance() - lam) / std::sqrt((lam + 2 * lam * lam) / N), 5);
  }
  // ---- Poisson by a tabulated distribution function: its frequencies against the probabilities ----
  for (double lam : {0.2, 2.0, 7.5, 40.0}) {
    const PoissonTable table(lam);
    Stream r(2026, 15, static_cast<std::uint32_t>(lam * 10), 0);
    const int n = 1000000;
    std::vector<double> count(static_cast<std::size_t>(lam + 12 * std::sqrt(lam) + 15), 0.0);
    Moments m;
    for (int i = 0; i < n; i++) {
      const auto k = table(r);
      m.add(static_cast<double>(k));
      count[std::min<std::size_t>(static_cast<std::size_t>(k), count.size() - 1)] += 1;
    }
    // Pearson's chi-square over the cells expecting at least 5 draws, the rest pooled.
    double chi2 = 0, p = std::exp(-lam), pooled_p = 0, pooled_c = 0;
    int cells = 0;
    for (std::size_t k = 0; k < count.size(); k++) {
      if (k > 0) p *= lam / static_cast<double>(k);
      const double e = n * p;
      if (e >= 5) { chi2 += (count[k] - e) * (count[k] - e) / e; cells++; }
      else { pooled_p += p; pooled_c += count[k]; }
    }
    if (pooled_p * n >= 5) { chi2 += (pooled_c - n * pooled_p) * (pooled_c - n * pooled_p) / (n * pooled_p); cells++; }
    // Against the chi-square law of cells - 1 degrees of freedom: z of its mean and spread.
    const double df = cells - 1, z = (chi2 - df) / std::sqrt(2 * df);
    const std::string name = "random: tabulated Poisson (mean " + calx::verification::shortest(lam) + ")";
    s.check(name + ", chi-square of its frequencies against the probabilities (z)", z, 5);
    s.check(name + ", z of the mean", (m.mean - lam) / std::sqrt(lam / n), 5);
  }
  {
    const PoissonTable none(0), big(800);
    Stream r(2026, 16, 0, 0);
    double sum = 0;
    for (int i = 0; i < 10000; i++) sum += static_cast<double>(big(r));
    s.require("random: a tabulated Poisson of mean 0 draws 0, and one above 500 goes to PTRS", none(r) == 0 && big.size() == 0 && std::fabs(sum / 10000 - 800) < 5 * std::sqrt(800.0 / 10000));
  }
  // ---- gamma, beta and binomial: moments as z-scores, and binomial frequencies by chi-square ----------
  {
    Stream r(77, 0x67616d6d, 0, 0);
    const double N = 1e6;
    for (double k : {0.3, 1.0, 2.5, 40.0}) {
      Moments m;
      for (int i = 0; i < N; i++) m.add(r.gamma(k));
      s.check("random: gamma (shape " + calx::verification::shortest(k) + "), z of the mean (k)", (m.mean - k) / std::sqrt(k / N), 5);
      // Var of the sample variance: (μ4 − σ⁴)/N with μ4 = 3k² + 6k for the gamma law.
      s.check("random: gamma (shape " + calx::verification::shortest(k) + "), z of the variance (k)", (m.variance() - k) / std::sqrt((3 * k * k + 6 * k - k * k) / N), 5);
    }
    Moments bm;
    for (int i = 0; i < N; i++) bm.add(r.beta(2, 5));
    s.check("random: beta (2, 5), z of the mean (2/7)", (bm.mean - 2.0 / 7) / std::sqrt((10.0 / (49 * 8)) / N), 5);
    for (auto [n, p] : {std::pair<std::int64_t, double>{20, 0.3}, {1000, 0.004}, {1000, 0.4}, {100000, 0.73}, {7, 0.9}}) {
      Moments m;
      std::vector<double> count(static_cast<std::size_t>(n) + 1, 0);
      const int draws = 400000;
      for (int i = 0; i < draws; i++) {
        const auto x = r.binomial(n, p);
        m.add(static_cast<double>(x));
        if (x >= 0 && x <= n) count[static_cast<std::size_t>(x)]++;
      }
      const double mean = n * p, var = n * p * (1 - p);
      const std::string name = "random: binomial (" + std::to_string(n) + ", " + calx::verification::shortest(p) + ")";
      s.check(name + ", z of the mean", (m.mean - mean) / std::sqrt(var / draws), 5);
      // Pearson's chi-square over cells expecting at least 5, against its law (z of mean and spread).
      double chi2 = 0, pooled_obs = 0, pooled_exp = 0;
      int cells = 0;
      for (std::int64_t k = 0; k <= n; k++) {
        const double logp = Stream::log_factorial(n) - Stream::log_factorial(k) - Stream::log_factorial(n - k) + k * std::log(p) + (n - k) * std::log1p(-p);
        const double e = draws * std::exp(logp);
        if (e >= 5) { chi2 += (count[static_cast<std::size_t>(k)] - e) * (count[static_cast<std::size_t>(k)] - e) / e; cells++; }
        else { pooled_obs += count[static_cast<std::size_t>(k)]; pooled_exp += e; }
      }
      if (pooled_exp >= 5) { chi2 += (pooled_obs - pooled_exp) * (pooled_obs - pooled_exp) / pooled_exp; cells++; }
      const double dof = cells - 1;
      s.check(name + ", chi-square of its frequencies against the probabilities (z)", (chi2 - dof) / std::sqrt(2 * dof), 5);
    }
  }

  return s.done();
}
