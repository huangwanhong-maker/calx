// Windows and spectra of sampled signals.
//
// The Hann window is taken in its DFT-even (periodic) form, w_j = (1 - cos(2 pi j / L)) / 2,
// j = 0 .. L-1, the form for spectral analysis (Harris 1978, §4).
//
// Welch's (1967) estimate of a power spectral density: the signal, sampled at fs, is cut into
// segments of length L that overlap by a fraction of it; each is detrended by its mean,
// windowed, and transformed; and the squared moduli are averaged over the segments. Scaled as a
// one-sided density, P_k = c |X_k|^2 / (fs sum_j w_j^2), c = 2 except at zero frequency and at
// the Nyquist frequency (c = 1), at frequencies f_k = k fs / L, k = 0 .. L/2. So scaled, the
// sum of P_k fs / L over k is the mean, over segments, of the windowed segments' power
// (sum_j w_j^2 x_j^2 / sum_j w_j^2), which for a stationary signal is its variance; and white
// noise of variance s^2 has the flat density 2 s^2 / fs.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <stdexcept>
#include <vector>

#include <calx/transforms/fft.hpp>

namespace calx::signal {

inline std::vector<double> hann(std::size_t L) {
  std::vector<double> w(L);
  for (std::size_t j = 0; j < L; j++) w[j] = 0.5 * (1 - std::cos(2 * std::numbers::pi * static_cast<double>(j) / static_cast<double>(L)));
  return w;
}

struct Spectrum {
  std::vector<double> frequency;  // Hz (if fs is in Hz)
  std::vector<double> density;    // the signal's unit squared per Hz
  std::size_t segments = 0;
};

// Welch's estimate with segments of L samples (L is shortened to the signal's length if longer)
// overlapping by `overlap` of a segment, Hann-windowed.
inline Spectrum welch(const std::vector<double>& x, double fs, std::size_t L, double overlap = 0.5) {
  if (!(fs > 0)) throw std::invalid_argument("calx::signal::welch: the sampling rate must be positive");
  if (!(overlap >= 0 && overlap < 1)) throw std::invalid_argument("calx::signal::welch: the overlap lies in [0, 1)");
  L = std::min(L, x.size());
  if (L < 2) throw std::invalid_argument("calx::signal::welch: too few samples");
  const auto step = std::max<std::size_t>(1, static_cast<std::size_t>(std::llround(static_cast<double>(L) * (1 - overlap))));
  const std::vector<double> w = hann(L);
  double ww = 0;
  for (double v : w) ww += v * v;
  Spectrum out;
  const std::size_t K = L / 2 + 1;
  out.density.assign(K, 0.0);
  std::vector<transforms::Complex> seg(L);
  for (std::size_t s = 0; s + L <= x.size(); s += step) {
    double mean = 0;
    for (std::size_t j = 0; j < L; j++) mean += x[s + j];
    mean /= static_cast<double>(L);
    for (std::size_t j = 0; j < L; j++) seg[j] = (x[s + j] - mean) * w[j];
    transforms::fft(seg);
    for (std::size_t k = 0; k < K; k++) out.density[k] += std::norm(seg[k]);
    out.segments++;
  }
  for (std::size_t k = 0; k < K; k++) {
    const double c = (k == 0 || (L % 2 == 0 && k == L / 2)) ? 1.0 : 2.0;
    out.density[k] *= c / (fs * ww * static_cast<double>(out.segments));
  }
  out.frequency.resize(K);
  for (std::size_t k = 0; k < K; k++) out.frequency[k] = static_cast<double>(k) * fs / static_cast<double>(L);
  return out;
}

}  // namespace calx::signal
