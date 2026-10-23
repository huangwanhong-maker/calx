# transforms: Fast Transforms

*Status: draft (2026-10-23): the fast Fourier transform of any length.* Cosine, wavelet and spherical-harmonic transforms are planned.

## 1. Content

`fft.hpp`: the discrete Fourier transform of a complex sequence, forward and inverse (divided by n), in O(n log n) for every length: the iterative radix-2 algorithm of Cooley and Tukey for powers of two, with twiddle factors computed directly rather than by recurrence; Bluestein's chirp-z identity for any other length, as a convolution of a power-of-two length, the chirp's phase reduced modulo 2π in integers. A real sequence's transform is given in full.

## 2. Checks

`tests/check_transforms.cpp`: against the defining sum, computed in long double, for every length from 1 to 64 and for nine more up to the prime 1009 (relative error 1.4 × 10⁻¹⁵); the inverse undoing the forward; Parseval's identity for a power of two and a prime; a pure tone in its two bins and nowhere else.

## 3. Dependencies

None (`random` in the checks).

## 4. References

- Cooley, J. W. and Tukey, J. W. (1965), An algorithm for the machine calculation of complex Fourier series, *Math. Comp.* 19, 297–301.
- Bluestein, L. (1970), A linear filtering approach to the computation of discrete Fourier transform, *IEEE Trans. Audio Electroacoust.* 18, 451–455.
