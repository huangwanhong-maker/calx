// Checks of calx/signal/spectrum.hpp: the periodogram, Thomson's multitaper estimate, cross spectra and
// coherence.
//
// - The periodogram against scipy.signal.periodogram (make_fixtures.py): density, one-sided, n = 1000
//   and 999 (an odd length has no Nyquist bin), rectangular with and without demeaning, and under the
//   periodic Hann window (not of unit energy: the normalisation by Σh² is checked); Parseval in closed
//   form: Σ_j S_j fs/n = Σ h² x² / Σ h² to rounding.
// - The multitaper estimate: the plain estimate is the mean of the eigenspectra and each eigenspectrum is
//   the periodogram under its taper (exactly); white noise lies flat at 2σ²/fs (simulation, 20 series),
//   and the fraction of bins outside the 95 % interval of χ²_2K/2K is 5 %; a sinusoid's line lies at its
//   bin, its power A²/2 is all there (Parseval), and the band ±W holds at least the mean concentration; adaptive
//   weights remove the broadband leakage that tapers of poor concentration (K = 12 > 2NW at NW = 4) put
//   into the plain average of an AR(1) of ρ = 0.95 (simulation, 300 series, against the exact spectrum);
//   the degrees of freedom, 2K plain and K at the ends, and between 2 and 2K adaptive.
// - Cross spectra and coherence against scipy.signal.welch, csd and coherence (L = 256 and 255); the
//   phase convention in closed form: y lags x by 3 samples, so arg S_xy = −2πf·3/fs.

#include <calx/random/stream.hpp>
#include <calx/signal/spectrum.hpp>
#include <calx/special_functions/gamma.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <complex>
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

// The largest difference relative to the largest reference value.
double scaled_difference(const std::vector<double>& a, const std::vector<double>& b) {
  if (a.size() != b.size()) return INFINITY;
  double m = 0, top = 0;
  for (std::size_t i = 0; i < a.size(); i++) {
    m = std::max(m, std::fabs(a[i] - b[i]));
    top = std::max(top, std::fabs(b[i]));
  }
  return m / top;
}

}  // namespace

int main() {
  calx::verification::Suite s;
  const auto ref = read_fixtures(CALX_SIGNAL_FIXTURES "/spectra-reference.txt");
  s.require("the reference values from scipy are there", ref.count("series/x1000") == 1);

  // ---- the periodogram against scipy ----
  for (const std::string name : {"x1000", "x999"}) {
    const std::vector<double>& x = ref.at("series/" + name);
    const double fs = 10;
    const auto box = cs::periodogram(x, fs);
    const auto raw = cs::periodogram(x, fs, {}, false);
    const auto w = cs::hann(x.size());
    const auto han = cs::periodogram(x, fs, w);
    s.check("periodogram: rectangular, demeaned, = scipy (" + name + "; relative to the peak)", scaled_difference(box.density, ref.at("periodogram/" + name + "/boxcar")), 1e-13);
    s.check("periodogram: rectangular, not demeaned, = scipy (" + name + ")", scaled_difference(raw.density, ref.at("periodogram/" + name + "/nodemean")), 1e-13);
    s.check("periodogram: under the periodic Hann window, = scipy (" + name + ")", scaled_difference(han.density, ref.at("periodogram/" + name + "/hann")), 1e-13);
    s.require("periodogram: ⌊n/2⌋ + 1 frequencies j fs/n (" + name + ")", box.frequency.size() == x.size() / 2 + 1 && std::fabs(box.frequency[1] - fs / static_cast<double>(x.size())) < 1e-15);
    // Parseval: Σ S_j fs/n = Σ h²(x − x̄)² / Σ h².
    double mean = 0;
    for (double v : x) mean += v;
    mean /= static_cast<double>(x.size());
    double power = 0, hh = 0, total = 0, total_box = 0, var = 0;
    for (std::size_t t = 0; t < x.size(); t++) {
      power += w[t] * w[t] * (x[t] - mean) * (x[t] - mean);
      hh += w[t] * w[t];
      var += (x[t] - mean) * (x[t] - mean);
    }
    var /= static_cast<double>(x.size());
    for (double p : han.density) total += p * fs / static_cast<double>(x.size());
    for (double p : box.density) total_box += p * fs / static_cast<double>(x.size());
    s.check("periodogram: Parseval, Σ S_j fs/n = Σ h²(x − x̄)²/Σ h² (Hann; " + name + "; relative)", total / (power / hh) - 1, 1e-13);
    s.check("periodogram: Parseval, the rectangular periodogram integrates to the variance (" + name + ")", total_box / var - 1, 1e-13);
  }

  // ---- the multitaper estimate: its algebra ----
  {
    const std::vector<double>& x = ref.at("series/x1000");
    const double fs = 10;
    const auto plain = cs::multitaper(x, fs, 4, 7, false);
    const auto tp = cs::dpss(x.size(), 4, 7);
    double mean_err = 0, eigen_err = 0;
    std::vector<double> xd = x;
    double m = 0;
    for (double v : x) m += v;
    for (double& v : xd) v -= m / static_cast<double>(x.size());
    for (std::size_t j = 0; j < plain.density.size(); j++) {
      double a = 0;
      for (std::size_t k = 0; k < 7; k++) a += plain.eigen[k][j];
      mean_err = std::max(mean_err, std::fabs(a / 7 - plain.density[j]) / plain.density[j]);
    }
    for (std::size_t k = 0; k < 7; k++) eigen_err = std::max(eigen_err, scaled_difference(plain.eigen[k], cs::periodogram(xd, fs, tp.taper[k], false).density));
    s.check("multitaper: the plain estimate is the mean of the eigenspectra", mean_err, 1e-15);
    s.check("multitaper: each eigenspectrum is the periodogram under its taper", eigen_err, 1e-15);
    s.require("multitaper: plain degrees of freedom 2K, K at 0 and n/2", plain.dof[1] == 14 && plain.dof[0] == 7 && plain.dof.back() == 7);
    s.check("multitaper: bandwidth W = NW fs / n", plain.bandwidth - 4 * fs / 1000, 1e-15);
  }

  // ---- white noise: flat, within its χ² spread ----
  {
    const std::size_t n = 4096, K = 7;
    const double fs = 2, sd = 1.5, level = 2 * sd * sd / fs;
    const double lo = calx::special_functions::chi_square_quantile(0.025, 2.0 * K) / (2.0 * K);
    const double hi = calx::special_functions::chi_square_quantile(0.975, 2.0 * K) / (2.0 * K);
    double mean = 0, mean_adaptive = 0, outside = 0, count = 0, dof_lo = 1e9, dof_hi = 0;
    for (std::uint32_t r = 0; r < 20; r++) {
      calx::random::Stream z(7, r, 0, 0);
      std::vector<double> x(n);
      for (auto& v : x) v = z.normal(0, sd);
      const auto mt = cs::multitaper(x, fs, 4, K, false);
      const auto ad = cs::multitaper(x, fs, 4, K, true);
      for (std::size_t j = 8; j + 8 < mt.density.size(); j++) {
        const double q = mt.density[j] / level;
        mean += q;
        mean_adaptive += ad.density[j] / level;
        outside += (q < lo || q > hi) ? 1 : 0;
        count++;
        dof_lo = std::min(dof_lo, ad.dof[j]);
        dof_hi = std::max(dof_hi, ad.dof[j]);
      }
    }
    s.check("multitaper: white noise flat at 2σ²/fs (plain; relative, 20 series of 4096)", mean / count - 1, 0.01);
    s.check("multitaper: and adaptive (relative)", mean_adaptive / count - 1, 0.01);
    s.check("multitaper: 5 % of the bins outside the 95 % interval of χ²_2K/2K (fraction − 0.05)", outside / count - 0.05, 0.01);
    s.require("multitaper: adaptive degrees of freedom between 2 and 2K on white noise", dof_lo >= 2 && dof_hi <= 14 + 1e-9, "from " + std::to_string(dof_lo) + " to " + std::to_string(dof_hi));
  }

  // ---- a sinusoid: its line and its power ----
  {
    const std::size_t n = 2048;
    const double fs = 1, A = 3, nw = 4;
    const std::size_t bin = 300;
    std::vector<double> x(n);
    for (std::size_t t = 0; t < n; t++) x[t] = A * std::cos(2 * std::numbers::pi * static_cast<double>(bin * t) / static_cast<double>(n) + 0.8);
    const auto mt = cs::multitaper(x, fs, nw, 7, false);
    std::size_t peak = 0;
    for (std::size_t j = 0; j < mt.density.size(); j++) if (mt.density[j] > mt.density[peak]) peak = j;
    double total = 0, band = 0;
    for (std::size_t j = 0; j < mt.density.size(); j++) {
      total += mt.density[j] * fs / static_cast<double>(n);
      if (std::fabs(mt.frequency[j] - mt.frequency[bin]) <= mt.bandwidth) band += mt.density[j] * fs / static_cast<double>(n);
    }
    const auto tp = cs::dpss(n, nw, 7);
    double lam = 0, beat = 0;
    for (double l : tp.concentration) lam += l / 7;
    // Parseval under each taper: Σ_t h_t² x_t² = (A²/2)(1 + Σ_t h_t² cos(4πf₀t + 2φ)), the second term the
    // transform of h² at 2f₀ (some 10⁻⁷ here).
    for (const auto& h : tp.taper)
      for (std::size_t t = 0; t < n; t++) beat += h[t] * h[t] * std::cos(4 * std::numbers::pi * static_cast<double>(bin * t) / static_cast<double>(n) + 1.6) / 7;
    s.require("multitaper: a sinusoid's line at its bin", peak == bin, "at bin " + std::to_string(peak));
    s.check("multitaper: its power (A²/2)(1 + mean_k Σ h_k² cos(4πf₀t + 2φ)) in the estimate (relative)", total / (A * A / 2 * (1 + beat)) - 1, 1e-12);
    // The bins within ±W hold at least the mean concentration (the band's continuous fraction).
    s.require("multitaper: and within ±W of the line at least the tapers' mean concentration of it", band / (A * A / 2) >= lam - 1e-9 && band <= total * (1 + 1e-12),
              "fraction " + std::to_string(band / (A * A / 2)) + ", mean λ " + std::to_string(lam));
  }

  // ---- adaptive weights against leakage: an AR(1) of ρ = 0.95 ----
  {
    const std::size_t n = 1024, K = 12;
    const double fs = 1, rho = 0.95, nw = 4;
    double plain_ratio = 0, adaptive_ratio = 0, count = 0;
    for (std::uint32_t r = 0; r < 300; r++) {
      calx::random::Stream z(11, r, 0, 0);
      std::vector<double> x(n);
      double v = z.normal() / std::sqrt(1 - rho * rho);
      for (std::size_t t = 0; t < n; t++) {
        v = rho * v + z.normal();
        x[t] = v;
      }
      const auto pl = cs::multitaper(x, fs, nw, K, false);
      const auto ad = cs::multitaper(x, fs, nw, K, true);
      for (std::size_t j = 0; j < pl.density.size(); j++) {
        const double f = pl.frequency[j];
        if (f < 0.375 * fs || j == n / 2) continue;
        const double truth = 2 / fs / std::norm(1.0 - rho * std::polar(1.0, -2 * std::numbers::pi * f / fs));
        plain_ratio += pl.density[j] / truth;
        adaptive_ratio += ad.density[j] / truth;
        count++;
      }
    }
    plain_ratio /= count;
    adaptive_ratio /= count;
    s.require("multitaper: K = 12 tapers at NW = 4 leak into the plain average above 3fs/8 (bias above 5 %)", plain_ratio > 1.05, "mean ratio to the truth " + std::to_string(plain_ratio));
    s.check("multitaper: adaptive weights remove that bias (mean ratio to the truth − 1, 300 series)", adaptive_ratio - 1, 0.02);
  }

  // ---- cross spectra and coherence against scipy ----
  {
    const std::vector<double>& x = ref.at("series/x4000");
    const std::vector<double>& y = ref.at("series/y4000");
    const double fs = 2;
    for (const std::size_t L : {std::size_t{256}, std::size_t{255}}) {
      const std::string key = "csd/" + std::to_string(L) + "/";
      const auto c = cs::cross_spectrum(x, y, fs, L);
      std::vector<double> re, im;
      for (const auto& v : c.sxy) {
        re.push_back(v.real());
        im.push_back(v.imag());
      }
      const std::string what = " (L = " + std::to_string(L) + ")";
      s.check("cross_spectrum: S_xx = scipy.signal.welch" + what, scaled_difference(c.sxx, ref.at(key + "pxx")), 1e-13);
      s.check("cross_spectrum: S_yy = scipy.signal.welch" + what, scaled_difference(c.syy, ref.at(key + "pyy")), 1e-13);
      s.check("cross_spectrum: Re S_xy = scipy.signal.csd" + what, scaled_difference(re, ref.at(key + "re")), 1e-13);
      s.check("cross_spectrum: Im S_xy = scipy.signal.csd (the phase convention conj(X) Y)" + what, scaled_difference(im, ref.at(key + "im")), 1e-13);
      s.check("coherence: = scipy.signal.coherence" + what, scaled_difference(cs::coherence(c), ref.at(key + "coh")), 1e-12);
      const auto w = cs::welch(x, fs, L);
      s.check("cross_spectrum: S_xx is welch()'s" + what, scaled_difference(c.sxx, w.density), 1e-15);
    }
    // y_t = 0.8 x_{t−3} + noise: the phase −2πfτ with τ = 3/fs, its slope fitted by least squares over the
    // bins of coherence ≥ 0.8 (the phase unwrapped about the expected line).
    const auto c = cs::cross_spectrum(x, y, fs, 256);
    const auto coh = cs::coherence(c);
    double sfp = 0, sff = 0;
    for (std::size_t k = 1; k < c.frequency.size(); k++) {
      if (coh[k] < 0.8) continue;
      const double w2 = 2 * std::numbers::pi * c.frequency[k];
      const double phase = -w2 * 3 / fs + std::remainder(std::arg(c.sxy[k]) + w2 * 3 / fs, 2 * std::numbers::pi);
      sfp += w2 * phase;
      sff += w2 * w2;
    }
    s.check("cross_spectrum: y lagging x by 3 samples has arg S_xy = −2πfτ: τ fitted (samples − 3; a wrong convention gives −6)", -sfp / sff * fs - 3, 0.05);
  }
  return s.done();
}
