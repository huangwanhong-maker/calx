# signal: Signal Processing

*Status: draft (2026-10-23): the Hann window and Welch's spectral density; (2026-10-31, [decision 0013](../../development/decisions/0013-multitaper-spectra-and-bispectra.md)) Slepian tapers, periodograms, Thomson's multitaper estimate, cross spectra and coherence, the bispectrum and bicoherence.* Digital filters, convolution and resampling are planned.

## 1. Concepts

A spectral estimate spreads a sampled signal's variance over frequency. Every estimate here is a one-sided density in one convention, scipy.signal's: a bin between zero and the Nyquist frequency holds the power of ±f (c = 2), the bins at zero and at the Nyquist frequency their own (c = 1), so that each integrates to the variance and white noise of variance s² lies flat at 2s²/fs.

The periodogram of a finite record sees the true spectrum through a window, whose sidelobes let strong frequencies leak into weak ones; a taper lowers the sidelobes, and averaging lowers the variance. Welch averages tapered segments, giving up resolution. Thomson keeps the whole record and averages over several orthogonal tapers, the discrete prolate spheroidal sequences, which concentrate as much of their energy as a sequence of their length can in a band ±W: K ≈ 2NW tapers give 2K degrees of freedom at resolution 2W, with their leakage measured by their concentrations λ_k. Adaptive weights let a taper that leaks count only where the spectrum stands well above its leakage.

Two series are related, frequency by frequency, by their cross spectrum, whose phase is the lag of one behind the other and whose normalised modulus, the coherence, is the fraction of either's variance that a linear filter of the other explains. A quadratic nonlinearity couples three frequencies f₁, f₂ and f₁ + f₂ in phase; the bispectrum, a third-order spectrum, is the average of that phase coupling, and the bicoherence its normalised modulus.

## 2. Content

| Header | Content |
|---|---|
| `spectrum.hpp` | `hann(L)` (DFT-even); `Spectrum{frequency, density, segments}`; `welch(x, fs, L, overlap)`; `periodogram(x, fs, taper, demean)` (normalised by Σh², as scipy); `multitaper(x, fs, nw, k, adaptive)` → `Multitaper{frequency, density, eigen, dof, nw, tapers, bandwidth}` (plain: 2K degrees of freedom, K at 0 and n/2; adaptive: Percival and Walden eqs. 368a, 369a, 370a); `cross_spectrum(x, y, fs, L, overlap)` → `CrossSpectrum{frequency, sxx, syy, sxy, segments}` (S_xy = mean conj(X) Y, phase arg Y − arg X); `coherence(c)` |
| `tapers.hpp` | `dpss(n, nw, k)` → `Tapers{n, nw, taper, concentration}`: bisection and inverse iteration on Slepian's tridiagonal matrix, O(nK) memory (n = 10⁵ in under half a second); scipy's signs; concentrations from the autocorrelations |
| `bispectrum.hpp` | `bispectrum(x, fs, L, overlap, hann_window, max_bin)` → `Bispectrum{frequency, bins, segments, value, bicoherence, at(j, k), b2(j, k)}` on the principal domain 0 ≤ k ≤ j, j + k ≤ L/2 |

## 3. Checks

- `tests/check_signal.cpp` (`check_calx_signal`, 6): the density of white noise flat at 2s²/fs and integrating to s² (within 0.4 %); a sinusoid's peak at its frequency, its power A²/2 in the window's main lobe to rounding; the window's symmetry.
- `tests/check_tapers.cpp` (`check_calx_tapers`, 45): tapers, signs and concentrations against scipy.signal.windows.dpss at n = 64, 1000, 4096 (NW = 2, 4) and 10⁵ (NW = 4), to 10⁻¹⁵ – 2 × 10⁻¹¹; orthonormality, parity, the concentrations' order; a concentration recomputed as an energy fraction by brute force.
- `tests/check_spectra.cpp` (`check_calx_spectra`, 39): the periodogram against scipy.signal.periodogram (n = 1000 and 999; rectangular with and without demeaning, Hann) and Parseval; the multitaper estimate's algebra, white noise flat within its χ² spread, a sinusoid's line and its power; adaptive weights removing the leakage of poorly concentrated tapers (K = 12 at NW = 4) on an AR(1) of ρ = 0.95, against its exact spectrum; cross spectra and coherence against scipy.signal.welch, csd and coherence (L = 256, 255); the phase of a lag of three samples.
- `tests/check_bispectrum.cpp` (`check_calx_bispectrum`, 16): a phase-coupled triad in closed form, B = (L/2)³; B and b² against a brute-force transform; quadratic phase coupling (b² ≈ 1) against independent phases (b² ≈ 1/S); Gaussian noise at 1/S; the principal domain and max_bin.

The reference values are in `tests/fixtures/spectra-reference.txt`, written by `tests/make_fixtures.py` from scipy 1.18.1 (numpy 2.3.5) with the input series it drew.

## 4. Dependencies

`transforms` (`random` and `special_functions` in the checks).

## 5. References

- Harris, F. J. (1978), On the use of windows for harmonic analysis with the discrete Fourier transform, *Proc. IEEE* 66, 51–83.
- Hinich, M. J. (1982), Testing for Gaussianity and linearity of a stationary time series, *J. Time Ser. Anal.* 3, 169–176.
- Kim, Y. C. and Powers, E. J. (1979), Digital bispectral analysis and its applications to nonlinear wave interactions, *IEEE Trans. Plasma Sci.* 7, 120–131.
- Percival, D. B. and Walden, A. T. (1993), *Spectral Analysis for Physical Applications: Multitaper and Conventional Univariate Techniques*, Cambridge University Press.
- Slepian, D. (1978), Prolate spheroidal wave functions, Fourier analysis, and uncertainty — V: the discrete case, *Bell Syst. Tech. J.* 57, 1371–1430.
- Thomson, D. J. (1982), Spectrum estimation and harmonic analysis, *Proc. IEEE* 70, 1055–1096.
- Welch, P. D. (1967), The use of fast Fourier transform for the estimation of power spectra: a method based on time averaging over short, modified periodograms, *IEEE Trans. Audio Electroacoust.* 15, 70–73.

## 6. Open Questions

1. Thomson's F-test for lines in coloured noise, and the reshaping of the spectrum about them.
2. Filters, convolution and resampling.
3. Spectra of irregularly sampled series (tyche's Lomb–Scargle meanwhile) and in several dimensions.
