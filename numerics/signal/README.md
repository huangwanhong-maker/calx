# signal: Signal Processing

*Status: draft (2026-10-23): the Hann window and Welch's spectral density.* Digital filters, convolution and resampling are planned.

## 1. Content

`spectrum.hpp`: the Hann window in its DFT-even form; Welch's estimate of a power spectral density, from segments that overlap, each detrended by its mean, windowed and transformed (`transforms`), their squared moduli averaged and scaled as a one-sided density, so that it integrates to the signal's variance and white noise of variance s² lies flat at 2s²/fs.

## 2. Checks

`tests/check_signal.cpp`: the density of white noise flat at 2s²/fs and integrating to s² (within 0.4 %); a sinusoid's peak at its frequency, its power A²/2 in the window's main lobe to rounding; the window's symmetry.

## 3. Dependencies

`transforms` (`random` in the checks).

## 4. References

- Welch, P. D. (1967), The use of fast Fourier transform for the estimation of power spectra: a method based on time averaging over short, modified periodograms, *IEEE Trans. Audio Electroacoust.* 15, 70–73.
- Harris, F. J. (1978), On the use of windows for harmonic analysis with the discrete Fourier transform, *Proc. IEEE* 66, 51–83.
