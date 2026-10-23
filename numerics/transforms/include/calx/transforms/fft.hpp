// Fast Fourier transforms of complex sequences of any length.
//
// The forward transform is X_k = sum_{j=0}^{n-1} x_j exp(-2 pi i j k / n); the inverse has the
// opposite sign and is divided by n, so that inverse(forward(x)) = x.
//
// A length that is a power of two is transformed by the iterative radix-2 algorithm of Cooley
// and Tukey (1965), in place after a bit-reversal permutation; the twiddle factors are computed
// directly, cos and sin of 2 pi k / n, and not by recurrence, so that their error does not grow
// with n. Any other length goes by Bluestein's (1970) chirp-z identity
// j k = (j^2 + k^2 - (k - j)^2) / 2, which makes the transform a convolution of length
// m >= 2n - 1, a power of two, done by three transforms of length m; the chirp's phase
// pi j^2 / n is reduced modulo 2 pi exactly, in integers, before its cosine is taken.
// Either way the cost is O(n log n), and the error O(eps log n) relative to the norm of x.

#pragma once

#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <utility>
#include <vector>

namespace calx::transforms {

using Complex = std::complex<double>;

inline bool power_of_two(std::size_t n) { return n != 0 && (n & (n - 1)) == 0; }

namespace detail {

// In place, without the division by n of the inverse.
inline void radix2(std::vector<Complex>& a, bool inverse) {
  const std::size_t n = a.size();
  for (std::size_t i = 1, j = 0; i < n; i++) {
    std::size_t bit = n >> 1;
    for (; j & bit; bit >>= 1) j ^= bit;
    j ^= bit;
    if (i < j) std::swap(a[i], a[j]);
  }
  std::vector<Complex> w(n / 2);
  const double sign = inverse ? 1.0 : -1.0;
  for (std::size_t k = 0; k < n / 2; k++) {
    const double t = 2 * std::numbers::pi * static_cast<double>(k) / static_cast<double>(n);
    w[k] = {std::cos(t), sign * std::sin(t)};
  }
  for (std::size_t len = 2; len <= n; len <<= 1) {
    const std::size_t half = len / 2, step = n / len;
    for (std::size_t i = 0; i < n; i += len)
      for (std::size_t k = 0; k < half; k++) {
        const Complex u = a[i + k], v = a[i + k + half] * w[k * step];
        a[i + k] = u + v;
        a[i + k + half] = u - v;
      }
  }
}

// The forward transform of any length by Bluestein's chirp: with b_j = exp(i pi j^2 / n),
// X_k = conj(b_k) sum_j (x_j conj(b_j)) b_{k-j}, a convolution with the even sequence b.
inline void bluestein(std::vector<Complex>& x) {
  const std::size_t n = x.size();
  std::size_t m = 1;
  while (m < 2 * n - 1) m <<= 1;
  std::vector<Complex> b(n);
  for (std::size_t j = 0; j < n; j++) {
    const auto r = static_cast<std::uint64_t>(j) * j % (2 * static_cast<std::uint64_t>(n));
    const double t = std::numbers::pi * static_cast<double>(r) / static_cast<double>(n);
    b[j] = {std::cos(t), std::sin(t)};
  }
  std::vector<Complex> a(m), c(m);
  for (std::size_t j = 0; j < n; j++) a[j] = x[j] * std::conj(b[j]);
  c[0] = b[0];
  for (std::size_t j = 1; j < n; j++) c[j] = c[m - j] = b[j];
  radix2(a, false);
  radix2(c, false);
  for (std::size_t k = 0; k < m; k++) a[k] *= c[k];
  radix2(a, true);
  const double scale = 1.0 / static_cast<double>(m);
  for (std::size_t k = 0; k < n; k++) x[k] = std::conj(b[k]) * a[k] * scale;
}

}  // namespace detail

// The transform of x, in place; inverse = true for the inverse (divided by n).
inline void fft(std::vector<Complex>& x, bool inverse = false) {
  const std::size_t n = x.size();
  if (n <= 1) return;
  if (power_of_two(n)) {
    detail::radix2(x, inverse);
  } else if (!inverse) {
    detail::bluestein(x);
  } else {
    // The inverse by conjugation: conj(F(conj x)).
    for (auto& v : x) v = std::conj(v);
    detail::bluestein(x);
    for (auto& v : x) v = std::conj(v);
  }
  if (inverse) {
    const double scale = 1.0 / static_cast<double>(n);
    for (auto& v : x) v *= scale;
  }
}

// The transform of a real sequence (all n coefficients; X_{n-k} = conj(X_k)).
inline std::vector<Complex> fft(const std::vector<double>& x) {
  std::vector<Complex> out(x.begin(), x.end());
  fft(out);
  return out;
}

}  // namespace calx::transforms
