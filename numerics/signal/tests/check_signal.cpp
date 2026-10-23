// Checks of calx/signal: Welch's estimate of white noise is flat at 2 s^2 / fs and integrates to
// its variance; a sinusoid's power, A^2 / 2, lies in the main lobe at its frequency; and the
// Hann window is the DFT-even one.

#include <calx/signal/spectrum.hpp>
#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <numbers>
#include <vector>

namespace cs = calx::signal;

int main() {
  calx::verification::Suite s;
  const double fs = 1000;
  {
    const double sd = 3;
    calx::random::Stream z(23, 0, 0, 0);
    std::vector<double> x(1 << 18);
    for (auto& v : x) v = z.normal(0, sd);
    const auto P = cs::welch(x, fs, 1024);
    double mean = 0, total = 0;
    for (std::size_t k = 1; k + 1 < P.density.size(); k++) mean += P.density[k];
    mean /= static_cast<double>(P.density.size() - 2);
    for (double p : P.density) total += p * fs / 1024;
    s.check("signal: Welch's density of white noise is flat at 2 s^2 / fs (relative, mean over the bins)", mean / (2 * sd * sd / fs) - 1, 0.01);
    s.check("signal: and integrates to its variance (relative)", total / (sd * sd) - 1, 0.01);
    s.require("signal: 511 half-overlapping segments of 1024 in 2^18 samples", P.segments == 511);
  }
  {
    // A sin(2 pi f t), f on a bin (f = 125 Hz, L = 1024 at 1 kHz: bin 128): the Hann main lobe
    // spans the bin and its two neighbours, which hold A^2 / 2 between them.
    const double A = 2, f = 125;
    std::vector<double> x(1 << 16);
    for (std::size_t j = 0; j < x.size(); j++) x[j] = A * std::sin(2 * std::numbers::pi * f * static_cast<double>(j) / fs);
    const auto P = cs::welch(x, fs, 1024);
    std::size_t peak = 0;
    for (std::size_t k = 0; k < P.density.size(); k++) if (P.density[k] > P.density[peak]) peak = k;
    double lobe = 0;
    for (std::size_t k = peak - 1; k <= peak + 1; k++) lobe += P.density[k] * fs / 1024;
    s.check("signal: a sinusoid's spectrum peaks at its frequency (Hz)", P.frequency[peak] - f, 1e-12);
    s.check("signal: its main lobe holds A^2 / 2 (relative)", lobe / (A * A / 2) - 1, 1e-10);
  }
  {
    const auto w = cs::hann(8);
    s.check("signal: the DFT-even Hann window: w_0 = 0, w_{L/2} = 1, w_j = w_{L-j}", std::fabs(w[0]) + std::fabs(w[4] - 1) + std::fabs(w[1] - w[7]), 1e-15);
  }
  return s.done();
}
