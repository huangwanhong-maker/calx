// Discrete prolate spheroidal sequences (Slepian sequences), the tapers of Thomson's multitaper method.
//
// The law. Among the sequences h_0 … h_{n−1} of unit energy, the one whose transform
// H(f) = Σ_t h_t e^{−2πift} keeps the largest fraction λ = ∫_{−W}^{W} |H(f)|² df of its energy in the
// band |f| ≤ W (cycles per sample) is the leading eigenvector of the n × n matrix
// A_ts = sin 2πW(t − s) / π(t − s), A_tt = 2W, and λ is its eigenvalue; the k-th eigenvector, orthogonal
// to the first k, is the k-th discrete prolate spheroidal sequence, with 1 > λ_0 > λ_1 > … > 0, and about
// the first 2NW (NW = nW, the time–half-bandwidth product) have λ_k near 1 (Slepian 1978, Bell Syst. Tech.
// J. 57, 1371–1430; Percival and Walden 1993, Spectral Analysis for Physical Applications, §8.3,
// eqs. 377 and 378). Even k are symmetric in t ↦ n − 1 − t, odd k antisymmetric.
//
// The computation. A's leading eigenvalues cluster at 1, so its eigenvectors are ill-conditioned there.
// The symmetric tridiagonal matrix B that commutes with A has the same eigenvectors in the same order
// and gaps of order one between its leading eigenvalues (Slepian 1978; Percival and Walden eq. 378):
//   B_tt = ((n − 1)/2 − t)² cos 2πW,  B_{t−1,t} = B_{t,t−1} = t(n − t)/2,  t = 0 … n − 1.
// Its K largest eigenvalues are found by bisection on Sturm sequences, the count of the negative pivots
// of B − xI = LDLᵀ being the number of eigenvalues below x by Sylvester's law of inertia (Wilkinson 1965,
// The Algebraic Eigenvalue Problem, §5.37–5.39; Golub and Van Loan 2013, Matrix Computations, 4th ed.,
// §8.4.1), with a zero pivot moved to −pivmin as LAPACK's dstebz does; each eigenvector by inverse
// iteration (Golub and Van Loan §8.2.2) with B − θI factorised by Gaussian elimination with partial
// pivoting in O(n) (LAPACK's dgttrf and dgtts2), four solves from a generic start, each followed by the
// projection onto the vector's parity (exact for B, which is persymmetric) and Gram–Schmidt against the
// vectors already found. Memory O(nK); time O(nK) per bisection step and per solve: about 50 steps for
// each eigenvalue. The dense n × n matrix is never formed.
//
// The concentrations, λ_k = Σ_{|τ|<n} r_k(τ) sin(2πWτ)/(πτ) (the τ = 0 term 2W r_k(0)), with r_k the
// taper's autocorrelation Σ_t h_t h_{t+τ} computed by FFT (Percival and Walden §8.3, p. 390), as
// scipy.signal.windows.dpss(…, return_ratios=True) computes them.
//
// Signs, as scipy.signal.windows.dpss (after Percival and Walden p. 379): an even taper has a positive
// sum; an odd taper is positive at the first t where h_t² exceeds max(10⁻⁷, 1/n) (its first lobe).
//
// Accuracy. B is stored to a relative ε in each entry (its off-diagonals exactly), so its eigenvectors are
// determined to about ε·(n²/8)/gap, the gaps between B's leading eigenvalues being O(1) and independent of
// n (about 5 for NW = 2, 7 to 12 for NW = 4): some 10⁻¹¹ at n = 4096 and 10⁻⁸ at n = 10⁵, of no
// consequence for a spectral estimate. Checked against scipy (LAPACK's dstebz and dstein) at n = 64,
// 1000, 4096 and 10⁵ (tests/check_tapers.cpp).

#pragma once

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <vector>

#include <calx/transforms/fft.hpp>

namespace calx::signal {

struct Tapers {
  std::size_t n = 0;                        // length
  double nw = 0;                            // time–half-bandwidth product NW (W = NW/n cycles per sample)
  std::vector<std::vector<double>> taper;   // K tapers of length n, unit energy Σ h² = 1
  std::vector<double> concentration;        // λ_k ∈ (0, 1): the fraction of the taper's energy in |f| ≤ W
};

namespace detail {

// ---- the tridiagonal matrix B and its Sturm sequences ----

// The number of eigenvalues of the symmetric tridiagonal matrix (d, e) below x; e2[t] = e_t² couples
// t − 1 and t (e2[0] unused).
inline std::size_t sturm_count(const std::vector<double>& d, const std::vector<double>& e2, double x, double pivmin) {
  std::size_t count = 0;
  double q = d[0] - x;
  if (std::fabs(q) < pivmin) q = -pivmin;
  if (q < 0) count++;
  for (std::size_t t = 1; t < d.size(); t++) {
    q = (d[t] - x) - e2[t] / q;
    if (std::fabs(q) < pivmin) q = -pivmin;
    if (q < 0) count++;
  }
  return count;
}

// The eigenvalue of index m (ascending) by bisection in [lo, hi], where fewer than m + 1 eigenvalues lie
// below lo and at least m + 1 below hi; to a relative 4ε of its magnitude.
inline double bisect(const std::vector<double>& d, const std::vector<double>& e2, std::size_t m, double lo, double hi, double pivmin) {
  constexpr double eps = std::numeric_limits<double>::epsilon();
  for (int it = 0; it < 200; it++) {
    const double mid = 0.5 * (lo + hi);
    if (hi - lo <= 4 * eps * std::max(std::fabs(lo), std::fabs(hi)) + pivmin || mid == lo || mid == hi) break;
    if (sturm_count(d, e2, mid, pivmin) >= m + 1) hi = mid;
    else lo = mid;
  }
  return 0.5 * (lo + hi);
}

// ---- inverse iteration ----

// B − θI = P L U for B symmetric tridiagonal (d, e), by Gaussian elimination with partial pivoting
// (LAPACK's dgttrf): U has the diagonal u0, the superdiagonals u1 and u2; L the multipliers l; a pivot
// smaller than `tiny` is replaced by ±tiny, as inverse iteration requires near an eigenvalue.
struct ShiftedLu {
  std::vector<double> u0, u1, u2, l;
  std::vector<unsigned char> swapped;
};

inline ShiftedLu shifted_lu(const std::vector<double>& d, const std::vector<double>& e, double theta, double tiny) {
  const std::size_t n = d.size();
  ShiftedLu f;
  f.u0.resize(n);
  f.u1.assign(n > 0 ? n - 1 : 0, 0.0);
  f.u2.assign(n > 1 ? n - 2 : 0, 0.0);
  f.l.assign(n > 0 ? n - 1 : 0, 0.0);
  f.swapped.assign(n > 0 ? n - 1 : 0, 0);
  for (std::size_t t = 0; t < n; t++) f.u0[t] = d[t] - theta;
  for (std::size_t t = 0; t + 1 < n; t++) f.u1[t] = e[t + 1];
  auto guard = [tiny](double& p) {
    if (std::fabs(p) < tiny) p = p < 0 ? -tiny : tiny;
  };
  for (std::size_t t = 0; t + 1 < n; t++) {
    const double below = e[t + 1];  // the subdiagonal entry under the pivot
    if (std::fabs(f.u0[t]) >= std::fabs(below)) {
      guard(f.u0[t]);
      f.l[t] = below / f.u0[t];
      f.u0[t + 1] -= f.l[t] * f.u1[t];
    } else {
      // Rows t and t + 1 exchanged: row t becomes (e, u0[t+1], u1[t+1]).
      f.l[t] = f.u0[t] / below;
      f.u0[t] = below;
      const double next = f.u0[t + 1];
      f.u0[t + 1] = f.u1[t] - f.l[t] * next;
      f.u1[t] = next;
      if (t + 2 < n) {
        f.u2[t] = f.u1[t + 1];
        f.u1[t + 1] = -f.l[t] * f.u1[t + 1];
      }
      f.swapped[t] = 1;
    }
  }
  guard(f.u0[n - 1]);
  return f;
}

// x with (B − θI) x = r, in place.
inline void shifted_solve(const ShiftedLu& f, std::vector<double>& r) {
  const std::size_t n = r.size();
  for (std::size_t t = 0; t + 1 < n; t++) {
    if (f.swapped[t]) std::swap(r[t], r[t + 1]);
    r[t + 1] -= f.l[t] * r[t];
  }
  r[n - 1] /= f.u0[n - 1];
  if (n >= 2) r[n - 2] = (r[n - 2] - f.u1[n - 2] * r[n - 1]) / f.u0[n - 2];
  for (std::size_t t = n - 2; t-- > 0;) r[t] = (r[t] - f.u1[t] * r[t + 1] - f.u2[t] * r[t + 2]) / f.u0[t];
}

inline void normalise(std::vector<double>& v) {
  double s = 0;
  for (double x : v) s += x * x;
  s = 1 / std::sqrt(s);
  for (double& x : v) x *= s;
}

// The projection onto the symmetric (even) or antisymmetric (odd) sequences in t ↦ n − 1 − t.
inline void parity(std::vector<double>& v, bool odd) {
  const std::size_t n = v.size();
  for (std::size_t t = 0; t < n / 2; t++) {
    const double a = v[t], b = v[n - 1 - t];
    v[t] = 0.5 * (odd ? a - b : a + b);
    v[n - 1 - t] = odd ? -v[t] : v[t];
  }
  if (odd && n % 2 == 1) v[n / 2] = 0;
}

}  // namespace detail

// The k leading discrete prolate spheroidal sequences of length n with time–half-bandwidth product nw,
// 0 < nw < n/2, 1 ≤ k ≤ n, and their concentrations.
inline Tapers dpss(std::size_t n, double nw, std::size_t k) {
  if (n < 2) throw std::invalid_argument("calx::signal::dpss: a length of at least 2");
  if (!(nw > 0 && nw < 0.5 * static_cast<double>(n))) throw std::invalid_argument("calx::signal::dpss: 0 < NW < n/2");
  if (k < 1 || k > n) throw std::invalid_argument("calx::signal::dpss: between 1 and n tapers");
  const double dn = static_cast<double>(n), W = nw / dn, cw = std::cos(2 * std::numbers::pi * W);

  // ---- B ----
  std::vector<double> d(n), e(n, 0.0), e2(n, 0.0);
  double emax = 0;
  for (std::size_t t = 0; t < n; t++) {
    const double c = 0.5 * (dn - 1) - static_cast<double>(t);
    d[t] = c * c * cw;
    if (t > 0) {
      e[t] = 0.5 * static_cast<double>(t) * (dn - static_cast<double>(t));
      e2[t] = e[t] * e[t];
      emax = std::max(emax, e[t]);
    }
  }
  double lo = d[0], hi = d[0], norm = 0;
  for (std::size_t t = 0; t < n; t++) {
    const double r = (t > 0 ? e[t] : 0.0) + (t + 1 < n ? e[t + 1] : 0.0);
    lo = std::min(lo, d[t] - r);
    hi = std::max(hi, d[t] + r);
    norm = std::max(norm, std::fabs(d[t]) + r);
  }
  const double pivmin = std::numeric_limits<double>::min() * std::max(1.0, emax * emax);
  const double span = std::max(1.0, hi - lo) * 1e-14;
  lo -= span;
  hi += span;

  // ---- eigenvalues: the k-th largest first, then the others above it ----
  std::vector<double> theta(k);
  theta[k - 1] = detail::bisect(d, e2, n - k, lo, hi, pivmin);
  const double bottom = theta[k - 1] - 8 * std::numeric_limits<double>::epsilon() * std::fabs(theta[k - 1]) - pivmin;
  for (std::size_t j = 0; j + 1 < k; j++) theta[j] = detail::bisect(d, e2, n - 1 - j, bottom, hi, pivmin);

  // ---- eigenvectors ----
  Tapers out;
  out.n = n;
  out.nw = nw;
  out.taper.resize(k);
  const double tiny = std::numeric_limits<double>::epsilon() * norm;
  std::uint64_t state = 0x9E3779B97F4A7C15ull;
  for (std::size_t j = 0; j < k; j++) {
    const bool odd = j % 2 == 1;
    std::vector<double>& v = out.taper[j];
    v.resize(n);
    for (std::size_t t = 0; t < n; t++) {
      state ^= state << 13, state ^= state >> 7, state ^= state << 17;  // xorshift64: a generic start
      v[t] = 0.5 + static_cast<double>(state >> 11) * 0x1.0p-53;
    }
    const detail::ShiftedLu f = detail::shifted_lu(d, e, theta[j], tiny);
    for (int it = 0; it < 4; it++) {
      detail::parity(v, odd);
      detail::normalise(v);
      detail::shifted_solve(f, v);
      detail::parity(v, odd);
      for (std::size_t i = odd ? 1 : 0; i < j; i += 2) {
        double p = 0;
        for (std::size_t t = 0; t < n; t++) p += v[t] * out.taper[i][t];
        for (std::size_t t = 0; t < n; t++) v[t] -= p * out.taper[i][t];
      }
      detail::normalise(v);
    }
    // Signs as scipy's.
    if (!odd) {
      double s = 0;
      for (double x : v) s += x;
      if (s < 0) for (double& x : v) x = -x;
    } else {
      const double thresh = std::max(1e-7, 1 / dn);
      for (double x : v)
        if (x * x > thresh) {
          if (x < 0) for (double& y : v) y = -y;
          break;
        }
    }
  }

  // ---- concentrations from the autocorrelations ----
  std::size_t m = 1;
  while (m < 2 * n - 1) m <<= 1;
  std::vector<double> kernel(n);
  kernel[0] = 2 * W;
  for (std::size_t tau = 1; tau < n; tau++) {
    const double a = static_cast<double>(tau);
    kernel[tau] = 2 * std::sin(2 * std::numbers::pi * W * a) / (std::numbers::pi * a);
  }
  out.concentration.resize(k);
  std::vector<transforms::Complex> z(m);
  for (std::size_t j = 0; j < k; j++) {
    std::fill(z.begin(), z.end(), transforms::Complex{});
    for (std::size_t t = 0; t < n; t++) z[t] = out.taper[j][t];
    transforms::fft(z);
    for (auto& c : z) c = std::norm(c);
    transforms::fft(z, true);
    double s = 0;
    for (std::size_t tau = n; tau-- > 0;) s += z[tau].real() * kernel[tau];
    out.concentration[j] = s;
  }
  return out;
}

}  // namespace calx::signal
