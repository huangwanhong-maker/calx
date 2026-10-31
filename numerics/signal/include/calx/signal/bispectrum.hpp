// The bispectrum and the bicoherence of a sampled signal, by averaging over segments.
//
// The law. The bispectrum is the Fourier transform of the third-order cumulant; for a stationary signal
// it is estimated from the transforms X_s of S segments (each of L samples, demeaned and, by default,
// Hann-windowed as welch() windows them; spectrum.hpp) as
//   B(j, k) = (1/S) Σ_s X_s(j) X_s(k) conj X_s(j + k),
// at the frequencies f_j = j fs/L, f_k = k fs/L and f_j + f_k (Kim and Powers 1979, IEEE Trans. Plasma
// Sci. 7, 120–131; Hinich 1982, J. Time Ser. Anal. 3, 169–176). It vanishes for a Gaussian
// signal and for components whose phases are independent; it does not when a component at f_j + f_k is
// phase-coupled to components at f_j and f_k (quadratic phase coupling: θ_{j+k} = θ_j + θ_k, as a
// quadratic nonlinearity makes it). Its normalised modulus, the bicoherence
//   b²(j, k) = |Σ_s X_s(j) X_s(k) conj X_s(j+k)|² / (Σ_s |X_s(j) X_s(k)|² · Σ_s |X_s(j+k)|²)   (K&P eq. 13),
// lies in [0, 1] by the Cauchy–Schwarz inequality: 1 when the phase θ_j + θ_k − θ_{j+k} is the same in
// every segment, about 1/S when it is uniform and independent from segment to segment, and with mean
// near 1/S for Gaussian noise (2S b² is then asymptotically χ²₂: Elgar and Guza 1988, IEEE Trans.
// Acoust. Speech Signal Process. 36, 1667–1668; the level is tyche's).
//
// The domain. B(j, k) = B(k, j), and for a real signal B(j, k) = conj B(−j, −k); the principal domain
// 0 ≤ k ≤ j, j + k ≤ L/2 determines the rest. Only it is computed, with j, k ≤ max_bin when max_bin > 0
// (bins = min(L/2, max_bin) + 1), so that memory is O(bins²) rather than O(L²). B is in the units of x³
// (the transforms are not scaled); b² is dimensionless.
//
// Cost O(S (L log L + bins²/2)). Checks: tests/check_bispectrum.cpp (a deterministic triad in closed form,
// a brute-force transform, coupled and uncoupled phases, Gaussian noise).

#pragma once

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include <calx/signal/spectrum.hpp>
#include <calx/transforms/fft.hpp>

namespace calx::signal {

struct Bispectrum {
  std::vector<double> frequency;              // j fs/L for j = 0…bins−1
  std::size_t bins = 0;                       // min(L/2, max_bin) + 1
  std::size_t segments = 0;
  std::vector<std::complex<double>> value;    // row-major bins × bins; zero outside the principal domain
  std::vector<double> bicoherence;            // b², same layout

  // B(j, k), by B(j, k) = B(k, j); zero outside the principal domain.
  std::complex<double> at(std::size_t j, std::size_t k) const {
    if (k > j) std::swap(j, k);
    if (j >= bins) throw std::out_of_range("calx::signal::Bispectrum::at: a bin below bins");
    return value[j * bins + k];
  }
  double b2(std::size_t j, std::size_t k) const {
    if (k > j) std::swap(j, k);
    if (j >= bins) throw std::out_of_range("calx::signal::Bispectrum::b2: a bin below bins");
    return bicoherence[j * bins + k];
  }
};

// Segments of L samples (shortened to the signal's length if longer) overlapping by `overlap` of a
// segment, demeaned and Hann-windowed (rectangular if hann_window is false).
inline Bispectrum bispectrum(const std::vector<double>& x, double fs, std::size_t L, double overlap = 0.5, bool hann_window = true, std::size_t max_bin = 0) {
  if (!(fs > 0)) throw std::invalid_argument("calx::signal::bispectrum: the sampling rate must be positive");
  if (!(overlap >= 0 && overlap < 1)) throw std::invalid_argument("calx::signal::bispectrum: the overlap lies in [0, 1)");
  L = std::min(L, x.size());
  if (L < 4) throw std::invalid_argument("calx::signal::bispectrum: segments of at least 4 samples");
  const auto step = std::max<std::size_t>(1, static_cast<std::size_t>(std::llround(static_cast<double>(L) * (1 - overlap))));
  const std::size_t half = L / 2;
  const std::size_t top = max_bin > 0 ? std::min(half, max_bin) : half;
  Bispectrum out;
  out.bins = top + 1;
  const std::size_t B = out.bins;
  out.frequency.resize(B);
  for (std::size_t j = 0; j < B; j++) out.frequency[j] = static_cast<double>(j) * fs / static_cast<double>(L);
  out.value.assign(B * B, {0.0, 0.0});
  out.bicoherence.assign(B * B, 0.0);  // Σ_s |X(j) X(k)|² while summing
  std::vector<double> power(half + 1, 0.0);  // Σ_s |X(q)|², q = j + k ≤ L/2
  const std::vector<double> w = hann_window ? hann(L) : std::vector<double>(L, 1.0);
  std::vector<transforms::Complex> X(L);
  for (std::size_t s = 0; s + L <= x.size(); s += step) {
    double mean = 0;
    for (std::size_t t = 0; t < L; t++) mean += x[s + t];
    mean /= static_cast<double>(L);
    for (std::size_t t = 0; t < L; t++) X[t] = (x[s + t] - mean) * w[t];
    transforms::fft(X);
    for (std::size_t q = 0; q <= half; q++) power[q] += std::norm(X[q]);
    for (std::size_t j = 0; j < B; j++)
      for (std::size_t k = 0; k <= j && j + k <= half; k++) {
        const std::complex<double> p = X[j] * X[k];
        out.value[j * B + k] += p * std::conj(X[j + k]);
        out.bicoherence[j * B + k] += std::norm(p);
      }
    out.segments++;
  }
  const double S = static_cast<double>(out.segments);
  for (std::size_t j = 0; j < B; j++)
    for (std::size_t k = 0; k <= j && j + k <= half; k++) {
      const std::size_t i = j * B + k;
      const double den = out.bicoherence[i] * power[j + k];
      out.bicoherence[i] = den > 0 ? std::min(1.0, std::norm(out.value[i]) / den) : 0.0;
      out.value[i] /= S;
    }
  return out;
}

}  // namespace calx::signal
