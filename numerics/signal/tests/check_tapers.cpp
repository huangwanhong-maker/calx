// Checks of calx/signal/tapers.hpp: the discrete prolate spheroidal sequences.
//
// - Against scipy.signal.windows.dpss (make_fixtures.py; LAPACK's bisection and inverse iteration):
//   the tapers, with scipy's signs, and their concentrations, for n = 64, 1000, 4096 with NW = 2 and 4,
//   K = 2NW − 1, and n = 10⁵ with NW = 4 (every 1000th sample), timed.
// - Their defining properties: unit energy and mutual orthogonality; even tapers symmetric, odd ones
//   antisymmetric; concentrations decreasing in (0, 1); the leading taper's concentration equal to its
//   energy fraction in |f| ≤ W measured directly on a fine grid of frequencies (brute force).

#include <calx/signal/tapers.hpp>
#include <calx/verification/check.hpp>

#include <chrono>
#include <cmath>
#include <fstream>
#include <map>
#include <numbers>
#include <sstream>
#include <string>
#include <vector>

namespace cs = calx::signal;

namespace {

std::map<std::string, std::vector<double>> read_fixtures(const std::string& path) {
  std::map<std::string, std::vector<double>> out;
  std::ifstream in(path);
  std::string line, key;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream f(line);
    f >> key;
    std::vector<double> v;
    double x;
    while (f >> x) v.push_back(x);
    out[key] = std::move(v);
  }
  return out;
}

}  // namespace

int main() {
  calx::verification::Suite s;
  const auto ref = read_fixtures(CALX_SIGNAL_FIXTURES "/spectra-reference.txt");
  s.require("the reference values from scipy are there", ref.count("dpss/64/2/ratios") == 1);

  // ---- against scipy ----
  struct Case { std::size_t n; int nw; double taper_limit, ratio_limit; };
  for (const Case c : {Case{64, 2, 1e-12, 1e-13}, Case{64, 4, 1e-12, 1e-13}, Case{1000, 2, 1e-11, 1e-12}, Case{1000, 4, 1e-11, 1e-12},
                       Case{4096, 2, 1e-11, 1e-13}, Case{4096, 4, 1e-11, 1e-13}, Case{100000, 4, 1e-9, 1e-12}}) {
    const std::size_t K = static_cast<std::size_t>(2 * c.nw - 1);
    const auto t0 = std::chrono::steady_clock::now();
    const cs::Tapers tp = cs::dpss(c.n, c.nw, K);
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    const std::string base = "dpss/" + std::to_string(c.n) + "/" + std::to_string(c.nw) + "/";
    double worst = 0, worst_ratio = 0;
    bool signs = true;
    for (std::size_t k = 0; k < K; k++) {
      const auto& r = ref.at(base + "taper" + std::to_string(k));
      const auto stride = static_cast<std::size_t>(r[0]);
      for (std::size_t i = 1; i < r.size(); i++) {
        const double mine = tp.taper[k][(i - 1) * stride];
        worst = std::max(worst, std::fabs(mine - r[i]));
        if (std::fabs(r[i]) > 1e-3 && (mine > 0) != (r[i] > 0)) signs = false;
      }
      worst_ratio = std::max(worst_ratio, std::fabs(tp.concentration[k] - ref.at(base + "ratios")[k]));
    }
    const std::string what = "n = " + std::to_string(c.n) + ", NW = " + std::to_string(c.nw) + ", K = " + std::to_string(K);
    s.check("dpss: the tapers equal scipy's (" + what + "; largest difference)", worst, c.taper_limit);
    s.require("dpss: with scipy's signs (" + what + ")", signs);
    s.check("dpss: the concentrations λ_k equal scipy's (" + what + ")", worst_ratio, c.ratio_limit);
    if (c.n == 100000) s.check("dpss: n = 10⁵, NW = 4, K = 7 computed in O(nK) (seconds)", seconds, 10.0, "s");

    // ---- defining properties ----
    double gram = 0, par = 0;
    for (std::size_t i = 0; i < K; i++)
      for (std::size_t j = 0; j <= i; j++) {
        double p = 0;
        for (std::size_t t = 0; t < c.n; t++) p += tp.taper[i][t] * tp.taper[j][t];
        gram = std::max(gram, std::fabs(p - (i == j ? 1.0 : 0.0)));
      }
    for (std::size_t k = 0; k < K; k++)
      for (std::size_t t = 0; t < c.n; t++) par = std::max(par, std::fabs(tp.taper[k][t] - (k % 2 ? -1.0 : 1.0) * tp.taper[k][c.n - 1 - t]));
    s.check("dpss: unit energy and orthogonal, Σ h_i h_j = δ_ij (" + what + ")", gram, 1e-12);
    s.check("dpss: even tapers symmetric, odd antisymmetric (" + what + ")", par, 1e-15);
    bool ordered = tp.concentration[0] < 1 + 1e-15 && tp.concentration[K - 1] > 0;
    for (std::size_t k = 1; k < K; k++) ordered = ordered && tp.concentration[k] < tp.concentration[k - 1];
    s.require("dpss: 1 > λ_0 > λ_1 > … > 0 (" + what + ")", ordered);
  }

  // ---- the concentration as an energy fraction (brute force) ----
  {
    // λ = ∫_{−W}^{W} |H(f)|² df for the taper of order 3, n = 64, NW = 2 (λ well below 1), by the midpoint
    // rule on 20 000 frequencies in [−W, W], H(f) = Σ_t h_t e^{−2πift} summed directly.
    const std::size_t n = 64;
    const double nw = 2, W = nw / static_cast<double>(n);
    const cs::Tapers tp = cs::dpss(n, nw, 4);
    const int M = 20000;
    double lam = 0;
    for (int i = 0; i < M; i++) {
      const double f = -W + (i + 0.5) * 2 * W / M;
      std::complex<double> H = 0;
      for (std::size_t t = 0; t < n; t++) H += tp.taper[3][t] * std::polar(1.0, -2 * std::numbers::pi * f * static_cast<double>(t));
      lam += std::norm(H) * 2 * W / M;
    }
    s.check("dpss: λ_3 = ∫_{−W}^{W} |H_3(f)|² df (n = 64, NW = 2; midpoint rule on 20 000 frequencies)", lam - tp.concentration[3], 1e-8);
  }
  return s.done();
}
