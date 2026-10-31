// Checks of calx/signal/bispectrum.hpp: the bispectrum and the bicoherence.
//
// - A deterministic triad on bins, rectangular window, one segment, in closed form: x = cos(2πj₁t/L + φ₁)
//   + cos(2πj₂t/L + φ₂) + cos(2π(j₁+j₂)t/L + φ₁ + φ₂) has X(j) = (L/2) e^{iφ} at its bins, so
//   B(j₁, j₂) = (L/2)³ exactly (real and positive: the conjugate falls on the sum) and b² = 1; and B is
//   zero, to rounding, at pairs of bins that are not a triad.
// - Against a brute-force transform (the DFT summed directly) of Hann-windowed, demeaned, overlapping
//   segments of noise: B and b² at every pair of the principal domain.
// - Quadratic phase coupling (simulation, after Kim and Powers 1979): random phases φ₁, φ₂ in each of 64
//   segments, the third component at φ₁ + φ₂, give b² ≈ 1 at (j₁, j₂); an independent third phase gives
//   b² of mean 1/S over 800 records.
// - Gaussian noise: b² of mean 1/S over the interior of the principal domain (S = 64, no overlap).
// - b² ∈ [0, 1]; at(j, k) = at(k, j); bins and the frequencies with max_bin.

#include <calx/random/stream.hpp>
#include <calx/signal/bispectrum.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <complex>
#include <numbers>
#include <vector>

namespace cs = calx::signal;

int main() {
  calx::verification::Suite s;
  const double two_pi = 2 * std::numbers::pi;

  // ---- a deterministic triad, in closed form ----
  {
    const std::size_t L = 64, j1 = 9, j2 = 5;
    const double p1 = 0.7, p2 = -1.9, fs = 4;
    std::vector<double> x(L);
    for (std::size_t t = 0; t < L; t++) {
      const double u = two_pi * static_cast<double>(t) / static_cast<double>(L);
      x[t] = std::cos(u * j1 + p1) + std::cos(u * j2 + p2) + std::cos(u * (j1 + j2) + p1 + p2);
    }
    const auto b = cs::bispectrum(x, fs, L, 0.0, false);
    const double expected = std::pow(L / 2.0, 3);
    s.require("bispectrum: one segment of L, L/2 + 1 bins at j fs/L", b.segments == 1 && b.bins == L / 2 + 1 && std::fabs(b.frequency[3] - 3 * fs / L) < 1e-15);
    s.check("bispectrum: a phase-coupled triad on bins has B(j₁, j₂) = (L/2)³ (relative)", std::abs(b.at(j1, j2) - expected) / expected, 1e-13);
    s.check("bispectrum: its bicoherence is 1", b.b2(j1, j2) - 1, 1e-13);
    s.check("bispectrum: B(j₂, j₁) = B(j₁, j₂)", std::abs(b.at(j2, j1) - b.at(j1, j2)), 0.0);
    s.check("bispectrum: B vanishes at pairs that are no triad, e.g. (j₁, j₁) (relative)", std::abs(b.at(j1, j1)) / expected, 1e-13);
  }

  // ---- against the transform summed directly ----
  {
    const std::size_t n = 300, L = 40;
    const double fs = 1, overlap = 0.5;
    calx::random::Stream z(5, 0, 0, 0);
    std::vector<double> x(n);
    for (auto& v : x) v = z.normal() + 0.3 * z.normal() * z.normal();
    const auto b = cs::bispectrum(x, fs, L, overlap, true);
    const auto w = cs::hann(L);
    const std::size_t step = L / 2, half = L / 2;
    std::vector<std::complex<double>> num((half + 1) * (half + 1));
    std::vector<double> d1((half + 1) * (half + 1)), d2(half + 1);
    std::size_t S = 0;
    for (std::size_t s0 = 0; s0 + L <= n; s0 += step, S++) {
      double mean = 0;
      for (std::size_t t = 0; t < L; t++) mean += x[s0 + t];
      mean /= static_cast<double>(L);
      std::vector<std::complex<double>> X(L);
      for (std::size_t j = 0; j < L; j++)
        for (std::size_t t = 0; t < L; t++) X[j] += (x[s0 + t] - mean) * w[t] * std::polar(1.0, -two_pi * static_cast<double>(j * t) / static_cast<double>(L));
      for (std::size_t q = 0; q <= half; q++) d2[q] += std::norm(X[q]);
      for (std::size_t j = 0; j <= half; j++)
        for (std::size_t k = 0; k <= j && j + k <= half; k++) {
          num[j * (half + 1) + k] += X[j] * X[k] * std::conj(X[j + k]);
          d1[j * (half + 1) + k] += std::norm(X[j] * X[k]);
        }
    }
    double worst_b = 0, worst_b2 = 0, scale = 0;
    for (std::size_t j = 0; j <= half; j++)
      for (std::size_t k = 0; k <= j && j + k <= half; k++) {
        const std::size_t i = j * (half + 1) + k;
        const auto B = num[i] / static_cast<double>(S);
        scale = std::max(scale, std::abs(B));
        worst_b = std::max(worst_b, std::abs(b.at(j, k) - B));
        const double den = d1[i] * d2[j + k];
        worst_b2 = std::max(worst_b2, std::fabs(b.b2(j, k) - (den > 0 ? std::norm(num[i]) / den : 0.0)));
      }
    s.require("bispectrum: 14 half-overlapping segments of 40 in 300 samples", b.segments == S && S == 14);
    s.check("bispectrum: B = (1/S) Σ X(j) X(k) conj X(j+k) with the DFT summed directly (relative to max |B|)", worst_b / scale, 1e-12);
    s.check("bispectrum: b² as Kim and Powers's eq. 13, with the DFT summed directly", worst_b2, 1e-12);
    double outside = 0;
    for (std::size_t j = 0; j <= half; j++)
      for (std::size_t k = 0; k <= half; k++)
        if (k > j || j + k > half) outside = std::max(outside, std::abs(b.value[j * b.bins + k]) + b.bicoherence[j * b.bins + k]);
    s.check("bispectrum: zero outside the principal domain 0 ≤ k ≤ j, j + k ≤ L/2", outside, 0.0);
    const auto c = cs::bispectrum(x, fs, L, overlap, true, 6);
    double same = 0;
    for (std::size_t j = 0; j < 7; j++)
      for (std::size_t k = 0; k <= j; k++) same = std::max(same, std::abs(c.at(j, k) - b.at(j, k)) + std::fabs(c.b2(j, k) - b.b2(j, k)));
    s.require("bispectrum: max_bin = 6 keeps 7 bins", c.bins == 7 && c.value.size() == 49);
    s.check("bispectrum: and the same values there (relative to max |B|)", same / scale, 1e-15);
  }

  // ---- quadratic phase coupling ----
  {
    const std::size_t L = 128, S = 64, j1 = 23, j2 = 10;
    auto record = [&](calx::random::Stream& z, bool coupled) {
      std::vector<double> x;
      x.reserve(L * S);
      for (std::size_t seg = 0; seg < S; seg++) {
        const double p1 = z.uniform(0, two_pi), p2 = z.uniform(0, two_pi), p3 = z.uniform(0, two_pi);
        const double p = coupled ? p1 + p2 : p3;
        for (std::size_t t = 0; t < L; t++) {
          const double u = two_pi * static_cast<double>(t) / static_cast<double>(L);
          x.push_back(std::cos(u * j1 + p1) + std::cos(u * j2 + p2) + std::cos(u * (j1 + j2) + p) + 0.3 * z.normal());
        }
      }
      return x;
    };
    calx::random::Stream z(9, 0, 0, 0);
    const auto coupled = cs::bispectrum(record(z, true), 1.0, L, 0.0);
    s.require("bispectrum: 64 segments of 128 without overlap", coupled.segments == S);
    s.check("bispectrum: quadratic phase coupling gives b²(j₁, j₂) ≈ 1 (1 − b²)", 1 - coupled.b2(j1, j2), 0.03);
    double mean = 0;
    const int R = 800;
    for (int r = 0; r < R; r++) {
      calx::random::Stream zr(10, static_cast<std::uint32_t>(r), 0, 0);
      mean += cs::bispectrum(record(zr, false), 1.0, L, 0.0).b2(j1, j2) / R;
    }
    s.check("bispectrum: an independent third phase gives b² of mean 1/S (S b² − 1, 800 records)", mean * S - 1, 0.12);
  }

  // ---- Gaussian noise ----
  {
    const std::size_t L = 128, S = 64;
    double mean = 0, count = 0, top = 0, low = 1;
    for (std::uint32_t r = 0; r < 5; r++) {
      calx::random::Stream z(12, r, 0, 0);
      std::vector<double> x(L * S);
      for (auto& v : x) v = z.normal();
      const auto b = cs::bispectrum(x, 1.0, L, 0.0);
      for (std::size_t j = 1; j < L / 2; j++)
        for (std::size_t k = 1; k <= j && j + k < L / 2; k++) {
          mean += b.b2(j, k);
          count++;
        }
      for (double v : b.bicoherence) {
        top = std::max(top, v);
        low = std::min(low, v);
      }
    }
    s.check("bispectrum: Gaussian noise has b² of mean 1/S (S b̄² − 1, interior of the principal domain)", mean / count * S - 1, 0.08);
    s.require("bispectrum: b² ∈ [0, 1]", low >= 0 && top <= 1);
  }
  return s.done();
}
