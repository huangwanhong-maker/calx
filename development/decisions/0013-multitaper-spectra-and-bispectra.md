# 0013. Slepian Tapers, Multitaper Spectra, Cross Spectra and Bispectra

- **Status:** Accepted, provisionally (one-dimensional, real, evenly sampled series; spectra in one one-sided convention)
- **Date:** 2026-10-31
- **Decided by:** the infrastructure's street analysis ("street oceanography": activity on streets read as a forced system with daily, weekly, lunar and annual rhythms, their nonlinear intermodulation and red-noise nulls) needs spectral estimators of known variance and bias, the coupling of rhythms, and the coherence of two series. Claude wrote the methods and this record at the user's request; the statistics on them are tyche's (its decision 0005).

## 1. Context

`signal` held the Hann window and Welch's estimate (decision 0004). Welch's estimate trades resolution for variance by cutting the series into segments, which a record of hourly counts over years, whose lines lie a year, a week and a day apart, can ill afford: Thomson's multitaper method keeps the whole record and averages over orthogonal tapers instead, with a bias controlled by the tapers' concentration. Lines at the sums and differences of the forcing frequencies (a daily rhythm modulated weekly) are the mark of a quadratic nonlinearity, which the bispectrum and the bicoherence detect. The relation of two series (counts at two places, counts and weather) needs their cross spectrum and coherence.

## 2. Decision

### 2.1 `signal/tapers.hpp`

`dpss(n, nw, k)` → `Tapers{n, nw, taper, concentration}`: the discrete prolate spheroidal sequences, as the leading eigenvectors of the tridiagonal matrix that commutes with the concentration problem (Slepian 1978; Percival and Walden 1993, eq. 378). Eigenvalues by bisection on Sturm sequences, vectors by inverse iteration with a pivoted tridiagonal LU, each projected onto its parity and orthogonalised against the earlier ones. Concentrations from the tapers' autocorrelations by FFT. Signs as scipy's.

- **Cost.** O(nK) memory, about 50 Sturm sweeps per taper: 0.38 s for n = 10⁵, K = 7 on the development machine. The dense n × n matrix is never formed.
- **Accuracy.** Agrees with scipy (LAPACK's bisection and inverse iteration) to 10⁻¹⁵ at n = 64, 4 × 10⁻¹³ at n = 4096 and 2 × 10⁻¹¹ at n = 10⁵.

### 2.2 `signal/spectrum.hpp`

Beside `hann` and `welch`:

- `periodogram(x, fs, taper, demean)`: the (tapered) periodogram, normalised by Σh², as scipy.signal.periodogram.
- `multitaper(x, fs, nw, k, adaptive)` → `Multitaper{frequency, density, eigen, dof, nw, tapers, bandwidth}`: the plain average of the eigenspectra (2K degrees of freedom), or Thomson's adaptive weights iterated as Percival and Walden eqs. 368a and 369a, with the degrees of freedom of eq. 370a.
- `cross_spectrum(x, y, fs, L, overlap)` → `CrossSpectrum{frequency, sxx, syy, sxy, segments}` and `coherence(c)`: Welch's averages, as scipy.signal.csd and coherence. The cross density is the average of conj(X) Y, so its phase is arg Y − arg X (y lagging x by τ gives −2πfτ).

One convention throughout: one-sided densities, c_j = 2 between 0 and the Nyquist frequency and 1 at either end, so that every estimate integrates to the variance.

### 2.3 `signal/bispectrum.hpp`

`bispectrum(x, fs, L, overlap, hann_window, max_bin)` → `Bispectrum{frequency, bins, segments, value, bicoherence, at(), b2()}`: B(j, k) = (1/S) Σ_s X_s(j) X_s(k) conj X_s(j + k) and Kim and Powers's bicoherence on the principal domain 0 ≤ k ≤ j, j + k ≤ L/2. `max_bin` bounds memory at O(bins²).

## 3. Alternatives Considered

- **The dense eigenproblem of the sinc matrix** (calx's symmetric eigensolver). Rejected: its eigenvalues cluster at 1, so its vectors are ill-conditioned, and it costs O(n²) memory and O(n³) time, impossible at n = 10⁵.
- **LAPACK's dstemr (MRRR) or a borrowed DPSS routine.** Unnecessary: the tridiagonal's leading eigenvalues are separated by O(1) for every n, which bisection and inverse iteration resolve.
- **Interpolating tapers from a shorter length**, as some packages do for speed. Rejected: the exact tapers are cheap enough.
- **Adaptive weights by a fixed number of iterations.** Rejected for a convergence test, so that the estimate does not depend on an arbitrary count.
- **The inference on spectra here.** Rejected: confidence intervals, degrees of freedom of Welch's estimate, the levels of coherence and bicoherence and red-noise nulls are statistics, and belong to tyche (its decision 0005), which uses these estimators.

## 4. Consequences

- `check_calx_tapers` (45):
  - against scipy at n = 64, 1000 and 4096 with NW = 2 and 4, and at n = 10⁵: tapers, signs and concentrations;
  - orthonormality and parity;
  - a concentration recomputed as an energy fraction by brute force.
- `check_calx_spectra` (39):
  - the periodogram against scipy for even and odd lengths, rectangular (demeaned and not) and Hann, and Parseval;
  - the multitaper estimate's algebra; white noise flat within its χ² spread; a sinusoid's line and its power;
  - adaptive weights removing the leakage of poorly concentrated tapers on an AR(1);
  - cross spectra and coherence against scipy, and the phase of a lag.
- `check_calx_bispectrum` (16): a triad in closed form, a brute-force transform, quadratic phase coupling against independent phases, Gaussian noise at 1/S.
- The fixtures are made by `numerics/signal/tests/make_fixtures.py` from scipy 1.18.1.
- tyche's `time_series` builds its spectral inference on these estimators (tyche decision 0005).
- Open: filters, convolution and resampling (still planned); the multitaper F-test for lines (Thomson 1982); spectra of irregularly sampled series (tyche's Lomb–Scargle for now); spectra in several dimensions.
