// Checks of calx/random: Philox4x32-10 against Random123's known-answer vectors, and the
// moments of each distribution against their exact values, as z-scores of the sample mean
// (or variance) over a million draws; a correct sampler exceeds |z| = 5 with probability
// below 1e-6, and the seeds are fixed, so the checks are deterministic.

#include <calx/random/philox.hpp>
#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

using calx::random::Philox4x32;
using calx::random::Stream;

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
  return s.done();
}
