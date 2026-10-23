// Checks of calx/transforms: the fast transform against the defining sum, computed in long
// double, for every length from 1 to 64 and for lengths prime, composite and large; the
// inverse undoing the forward; Parseval's identity; and a pure tone landing in its one bin.

#include <calx/transforms/fft.hpp>
#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <numbers>
#include <string>
#include <vector>

namespace ct = calx::transforms;
using ct::Complex;

namespace {

std::vector<Complex> direct(const std::vector<Complex>& x) {
  const std::size_t n = x.size();
  std::vector<Complex> out(n);
  for (std::size_t k = 0; k < n; k++) {
    std::complex<long double> s = 0;
    for (std::size_t j = 0; j < n; j++) {
      const long double t = -2 * std::numbers::pi_v<long double> * static_cast<long double>((j * k) % n) / static_cast<long double>(n);
      s += std::complex<long double>(x[j].real(), x[j].imag()) * std::complex<long double>(std::cos(t), std::sin(t));
    }
    out[k] = {static_cast<double>(s.real()), static_cast<double>(s.imag())};
  }
  return out;
}

std::vector<Complex> noise(std::size_t n, std::uint32_t id) {
  calx::random::Stream z(17, id, 0, 0);
  std::vector<Complex> x(n);
  for (auto& v : x) v = {z.normal(), z.normal()};
  return x;
}

double worst(const std::vector<Complex>& a, const std::vector<Complex>& b) {
  double d = 0, m = 0;
  for (std::size_t i = 0; i < a.size(); i++) { d = std::max(d, std::abs(a[i] - b[i])); m = std::max(m, std::abs(b[i])); }
  return m > 0 ? d / m : d;
}

}  // namespace

int main() {
  calx::verification::Suite s;
  {
    double err = 0, back = 0;
    std::vector<std::size_t> lengths;
    for (std::size_t n = 1; n <= 64; n++) lengths.push_back(n);
    for (std::size_t n : {97u, 100u, 127u, 128u, 243u, 255u, 256u, 1000u, 1009u}) lengths.push_back(n);
    for (std::size_t n : lengths) {
      const auto x = noise(n, static_cast<std::uint32_t>(n));
      auto X = x;
      ct::fft(X);
      err = std::max(err, worst(X, direct(x)));
      ct::fft(X, true);
      back = std::max(back, worst(X, x));
    }
    s.check("transforms: the FFT against the defining sum, every length 1-64 and nine more up to 1009 (relative to max |X|)", err, 5e-14);
    s.check("transforms: the inverse undoes the forward (relative)", back, 5e-14);
  }
  {
    // Parseval: sum |x|^2 = sum |X|^2 / n, for a power of two and for a prime.
    double worst_p = 0;
    for (std::size_t n : {4096u, 4093u}) {
      const auto x = noise(n, 9000 + static_cast<std::uint32_t>(n));
      auto X = x;
      ct::fft(X);
      double a = 0, b = 0;
      for (std::size_t i = 0; i < n; i++) { a += std::norm(x[i]); b += std::norm(X[i]); }
      worst_p = std::max(worst_p, std::fabs(b / static_cast<double>(n) / a - 1));
    }
    s.check("transforms: Parseval's identity, n = 4096 and the prime 4093 (relative)", worst_p, 1e-13);
  }
  {
    // cos(2 pi 5 j / 64) has X_5 = X_59 = 32 and nothing else.
    std::vector<double> x(64);
    for (std::size_t j = 0; j < 64; j++) x[j] = std::cos(2 * std::numbers::pi * 5 * static_cast<double>(j) / 64);
    const auto X = ct::fft(x);
    double off = 0;
    for (std::size_t k = 0; k < 64; k++) if (k != 5 && k != 59) off = std::max(off, std::abs(X[k]));
    s.check("transforms: a pure tone lands in its two bins (|X_5| - 32)", std::abs(X[5]) - 32, 1e-13);
    s.check("transforms: and nowhere else", off, 1e-13);
  }
  return s.done();
}
