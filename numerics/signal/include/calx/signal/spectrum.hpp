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
//
// Every estimate here is a one-sided density in this one convention (scipy.signal's
// scaling='density', return_onesided=True): a bin j of 0 < j < n/2 holds the power of the
// frequencies ±f_j, so c_j = 2; the bins at 0 and, for even n, at n/2 are their own mirror images,
// so c_j = 1, and the expectation there is half the continuous one-sided density.
//
// - periodogram: the (tapered) periodogram S_j = c_j |Σ_t h_t x_t e^{−2πi j t/n}|² / (fs Σ_t h_t²),
//   f_j = j fs/n, j = 0 … ⌊n/2⌋; the rectangular taper when none is given (h_t = 1); the series
//   demeaned first unless asked not to (scipy.signal.periodogram with detrend='constant' or False).
//   For a taper of unit energy this is (c_j/fs) |Σ h_t x_t e^{−2πijt/n}|²; Σ_j S_j fs/n is the
//   tapered power Σ h_t² x_t² / Σ h_t² (Parseval), the variance for the rectangular taper.
// - multitaper: Thomson's (1982, Proc. IEEE 70, 1055–1096) estimate from the K leading discrete
//   prolate spheroidal sequences (tapers.hpp): the eigenspectra Ŝ_k are the periodograms under each
//   taper of the demeaned series; the plain estimate is their mean, with 2K equivalent degrees of
//   freedom (K at 0 and at n/2, where the eigencoefficients are real: each Ŝ_k is then χ²₁ and not
//   χ²₂ about its mean); the adaptive estimate weights them by
//     d_k(f) = √λ_k S(f) / (λ_k S(f) + (1 − λ_k) σ² c_j/fs)      (Percival and Walden 1993, eq. 368a),
//     S(f) = Σ_k d_k² Ŝ_k / Σ_k d_k²                                (eq. 369a),
//   iterated from the mean of the first two eigenspectra until the largest relative change is below
//   10⁻¹⁰ (at most 100 iterations), σ² the series' variance (divisor n) and (1 − λ_k) σ² c_j/fs the
//   broadband bias of the k-th eigenspectrum of white noise of that variance; its degrees of freedom
//   ν(f) = 2 (Σ_k d_k²)² / Σ_k d_k⁴ (eq. 370a; halved at 0 and n/2). The adaptive weights let the
//   tapers that leak (λ_k well below 1) count only where the spectrum stands high above its
//   broadband leakage.
// - cross_spectrum: Welch's averages of the segments' conj(X) X, conj(Y) Y and conj(X) Y, Hann-
//   windowed and demeaned as welch() (scipy.signal.welch and csd). The cross density's phase is
//   arg Y − arg X: if y lags x by τ (y_t = x_{t−τ}), arg S_xy(f) = −2πfτ.
// - coherence: |S_xy|² / (S_xx S_yy) ∈ [0, 1] (scipy.signal.coherence); with K independent segments of
//   independent Gaussian series it is Beta(1, K − 1) at each frequency 0 < f < fs/2.

#pragma once

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <numbers>
#include <span>
#include <stdexcept>
#include <vector>

#include <calx/signal/tapers.hpp>
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

namespace detail {

// c_j: 1 at 0 and, for even n, at n/2; 2 between.
inline double one_sided(std::size_t j, std::size_t n) { return (j == 0 || (n % 2 == 0 && j == n / 2)) ? 1.0 : 2.0; }

inline std::vector<double> frequencies(std::size_t n, double fs) {
  std::vector<double> f(n / 2 + 1);
  for (std::size_t j = 0; j < f.size(); j++) f[j] = static_cast<double>(j) * fs / static_cast<double>(n);
  return f;
}

}  // namespace detail

// ---- the periodogram ----

inline Spectrum periodogram(const std::vector<double>& x, double fs, std::span<const double> taper = {}, bool demean = true) {
  const std::size_t n = x.size();
  if (!(fs > 0)) throw std::invalid_argument("calx::signal::periodogram: the sampling rate must be positive");
  if (n < 2) throw std::invalid_argument("calx::signal::periodogram: too few samples");
  if (!taper.empty() && taper.size() != n) throw std::invalid_argument("calx::signal::periodogram: a taper as long as the series");
  double mean = 0;
  if (demean) {
    for (double v : x) mean += v;
    mean /= static_cast<double>(n);
  }
  double hh = 0;
  std::vector<transforms::Complex> z(n);
  for (std::size_t t = 0; t < n; t++) {
    const double h = taper.empty() ? 1.0 : taper[t];
    hh += h * h;
    z[t] = h * (x[t] - mean);
  }
  if (!(hh > 0)) throw std::invalid_argument("calx::signal::periodogram: a taper that is not zero");
  transforms::fft(z);
  Spectrum out;
  out.frequency = detail::frequencies(n, fs);
  out.density.resize(out.frequency.size());
  for (std::size_t j = 0; j < out.density.size(); j++) out.density[j] = detail::one_sided(j, n) * std::norm(z[j]) / (fs * hh);
  out.segments = 1;
  return out;
}

// ---- Thomson's multitaper estimate ----

struct Multitaper {
  std::vector<double> frequency;              // j fs/n, j = 0…⌊n/2⌋
  std::vector<double> density;                // one-sided estimate
  std::vector<std::vector<double>> eigen;     // the K eigenspectra (one-sided), as periodogram() under each taper
  std::vector<double> dof;                    // equivalent degrees of freedom per frequency: 2K for the plain
                                              // average; Percival and Walden eq. 370a for the adaptive weights
                                              // (K and half of eq. 370a at 0 and n/2)
  double nw = 0;
  std::size_t tapers = 0;
  double bandwidth = 0;                       // W in frequency units: nw · fs / n
};

inline Multitaper multitaper(const std::vector<double>& x, double fs, double nw, std::size_t k, bool adaptive = true) {
  const std::size_t n = x.size();
  if (!(fs > 0)) throw std::invalid_argument("calx::signal::multitaper: the sampling rate must be positive");
  if (n < 2) throw std::invalid_argument("calx::signal::multitaper: too few samples");
  const Tapers tp = dpss(n, nw, k);
  double mean = 0;
  for (double v : x) mean += v;
  mean /= static_cast<double>(n);
  std::vector<double> xd(n);
  double variance = 0;
  for (std::size_t t = 0; t < n; t++) {
    xd[t] = x[t] - mean;
    variance += xd[t] * xd[t];
  }
  variance /= static_cast<double>(n);

  Multitaper out;
  out.frequency = detail::frequencies(n, fs);
  out.nw = nw;
  out.tapers = k;
  out.bandwidth = nw * fs / static_cast<double>(n);
  const std::size_t J = out.frequency.size();
  out.eigen.resize(k);
  for (std::size_t i = 0; i < k; i++) out.eigen[i] = periodogram(xd, fs, tp.taper[i], false).density;
  out.density.assign(J, 0.0);
  out.dof.resize(J);
  for (std::size_t j = 0; j < J; j++) {
    for (std::size_t i = 0; i < k; i++) out.density[j] += out.eigen[i][j];
    out.density[j] /= static_cast<double>(k);
    out.dof[j] = static_cast<double>(k) * detail::one_sided(j, n);
  }
  if (!adaptive || k == 1 || !(variance > 0)) return out;

  // ---- adaptive weights (Percival and Walden eqs. 368a, 369a, 370a) ----
  const std::vector<double>& lambda = tp.concentration;
  std::vector<double> S(J);
  for (std::size_t j = 0; j < J; j++) S[j] = 0.5 * (out.eigen[0][j] + out.eigen[1][j]);
  std::vector<double> next(J), d2(k);
  for (int it = 0; it < 100; it++) {
    double change = 0;
    for (std::size_t j = 0; j < J; j++) {
      const double broadband = variance * detail::one_sided(j, n) / fs;
      double num = 0, den = 0;
      for (std::size_t i = 0; i < k; i++) {
        const double dk = std::sqrt(lambda[i]) * S[j] / (lambda[i] * S[j] + (1 - lambda[i]) * broadband);
        num += dk * dk * out.eigen[i][j];
        den += dk * dk;
      }
      next[j] = den > 0 ? num / den : out.density[j];
      if (S[j] > 0) change = std::max(change, std::fabs(next[j] - S[j]) / S[j]);
    }
    S.swap(next);
    if (change < 1e-10) break;
  }
  for (std::size_t j = 0; j < J; j++) {
    const double broadband = variance * detail::one_sided(j, n) / fs;
    double s2 = 0, s4 = 0;
    for (std::size_t i = 0; i < k; i++) {
      const double dk = std::sqrt(lambda[i]) * S[j] / (lambda[i] * S[j] + (1 - lambda[i]) * broadband);
      s2 += dk * dk;
      s4 += dk * dk * dk * dk;
    }
    out.density[j] = S[j];
    if (s4 > 0) out.dof[j] = 2 * s2 * s2 / s4 * (detail::one_sided(j, n) / 2);
  }
  return out;
}

// ---- cross spectra and coherence ----

struct CrossSpectrum {
  std::vector<double> frequency;              // j fs/L, j = 0…L/2
  std::vector<double> sxx, syy;               // one-sided densities (Welch, Hann, as welch())
  std::vector<std::complex<double>> sxy;      // one-sided cross density, scipy.signal.csd's convention:
                                              // averages of conj(X_k) Y_k — its phase is arg Y − arg X
  std::size_t segments = 0;
};

// Welch's estimates of the two densities and of the cross density, segments of L samples (shortened to
// the series' length if longer) overlapping by `overlap` of a segment, each demeaned and Hann-windowed.
inline CrossSpectrum cross_spectrum(const std::vector<double>& x, const std::vector<double>& y, double fs, std::size_t L, double overlap = 0.5) {
  if (x.size() != y.size()) throw std::invalid_argument("calx::signal::cross_spectrum: two series of one length");
  if (!(fs > 0)) throw std::invalid_argument("calx::signal::cross_spectrum: the sampling rate must be positive");
  if (!(overlap >= 0 && overlap < 1)) throw std::invalid_argument("calx::signal::cross_spectrum: the overlap lies in [0, 1)");
  L = std::min(L, x.size());
  if (L < 2) throw std::invalid_argument("calx::signal::cross_spectrum: too few samples");
  const auto step = std::max<std::size_t>(1, static_cast<std::size_t>(std::llround(static_cast<double>(L) * (1 - overlap))));
  const std::vector<double> w = hann(L);
  double ww = 0;
  for (double v : w) ww += v * v;
  CrossSpectrum out;
  const std::size_t K = L / 2 + 1;
  out.sxx.assign(K, 0.0);
  out.syy.assign(K, 0.0);
  out.sxy.assign(K, {0.0, 0.0});
  std::vector<transforms::Complex> a(L), b(L);
  for (std::size_t s = 0; s + L <= x.size(); s += step) {
    double mx = 0, my = 0;
    for (std::size_t j = 0; j < L; j++) {
      mx += x[s + j];
      my += y[s + j];
    }
    mx /= static_cast<double>(L);
    my /= static_cast<double>(L);
    for (std::size_t j = 0; j < L; j++) {
      a[j] = (x[s + j] - mx) * w[j];
      b[j] = (y[s + j] - my) * w[j];
    }
    transforms::fft(a);
    transforms::fft(b);
    for (std::size_t k = 0; k < K; k++) {
      out.sxx[k] += std::norm(a[k]);
      out.syy[k] += std::norm(b[k]);
      out.sxy[k] += std::conj(a[k]) * b[k];
    }
    out.segments++;
  }
  for (std::size_t k = 0; k < K; k++) {
    const double c = detail::one_sided(k, L) / (fs * ww * static_cast<double>(out.segments));
    out.sxx[k] *= c;
    out.syy[k] *= c;
    out.sxy[k] *= c;
  }
  out.frequency = detail::frequencies(L, fs);
  return out;
}

// The magnitude-squared coherence; zero where either density vanishes.
inline std::vector<double> coherence(const CrossSpectrum& c) {
  std::vector<double> out(c.sxy.size(), 0.0);
  for (std::size_t k = 0; k < out.size(); k++) {
    const double den = c.sxx[k] * c.syy[k];
    if (den > 0) out[k] = std::min(1.0, std::norm(c.sxy[k]) / den);
  }
  return out;
}

}  // namespace calx::signal
