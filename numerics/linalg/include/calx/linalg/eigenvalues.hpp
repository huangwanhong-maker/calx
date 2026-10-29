// The eigenvalues of a real square matrix, in general complex.
//
// Three steps, each a similarity: balancing by powers of two (Parlett and Reinsch 1969, Numer.
// Math. 13, 293–304; the procedure balanc of Press et al. 2007, Numerical Recipes, 3rd ed.,
// §11.6.1), which is exact in binary arithmetic and makes the rows and columns of comparable norm;
// Householder reduction to upper Hessenberg form (the procedure orthes of Martin and Wilkinson,
// in Wilkinson and Reinsch 1971, Handbook for Automatic Computation II, contribution II/13); and
// the Francis double-shift QR iteration on the Hessenberg matrix, with Wilkinson's exceptional
// shifts after 10 and 30 iterations without deflation (the procedure hqr of Martin, Peters and
// Wilkinson, contribution II/14, as EISPACK has it and as the public-domain JAMA library of
// NIST and MathWorks transcribes it in hqr2, here without accumulating the Schur vectors).
// Golub and Van Loan (2013), Matrix Computations, 4th ed., Algorithms 7.4.2 and 7.5.1–7.5.2.
//
// Accuracy. The QR algorithm is backward stable: the computed eigenvalues are exact for A + E
// with ‖E‖ about ε‖A‖. A simple eigenvalue λ is then found to within about ε‖A‖ / s(λ), s the
// cosine between its left and right eigenvectors (Golub and Van Loan, §7.2.2); for a normal
// matrix s = 1, for a defective one the error grows like ε^{1/m} for a block of order m. A
// complex-conjugate pair comes out exactly conjugate. O(10 n³) flops; it serves the small
// matrices of statistics and dynamics: the companion matrix of a vector autoregression (order
// K p) and of an AR or MA polynomial, whose spectral radius decides stability and invertibility.
//
// The eigenvalues are returned in the order in which they deflate (from the bottom of the
// Hessenberg matrix upwards), not sorted.

#pragma once

#include <calx/linalg/dense.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace calx::linalg {

// Balancing in place (Parlett and Reinsch 1969): D⁻¹ A D with D a diagonal of powers of two.
inline void balance(Matrix& A) {
  constexpr double radix = 2, sqrdx = radix * radix;
  const std::size_t n = A.n;
  bool done = false;
  while (!done) {
    done = true;
    for (std::size_t i = 0; i < n; i++) {
      double r = 0, c = 0;
      for (std::size_t j = 0; j < n; j++)
        if (j != i) {
          c += std::fabs(A(j, i));
          r += std::fabs(A(i, j));
        }
      if (c == 0 || r == 0) continue;
      double g = r / radix, f = 1;
      const double s = c + r;
      while (c < g) {
        f *= radix;
        c *= sqrdx;
      }
      g = r * radix;
      while (c > g) {
        f /= radix;
        c /= sqrdx;
      }
      if ((c + r) / f < 0.95 * s) {
        done = false;
        g = 1 / f;
        for (std::size_t j = 0; j < n; j++) A(i, j) *= g;
        for (std::size_t j = 0; j < n; j++) A(j, i) *= f;
      }
    }
  }
}

// Reduction to upper Hessenberg form by Householder similarities (orthes), in place; the
// entries below the subdiagonal are set to zero.
inline void hessenberg(Matrix& H) {
  const std::size_t n = H.n;
  if (n < 3) return;
  std::vector<double> ort(n, 0.0);
  const std::size_t low = 0, high = n - 1;
  for (std::size_t m = low + 1; m + 1 <= high; m++) {
    double scale = 0;
    for (std::size_t i = m; i <= high; i++) scale += std::fabs(H(i, m - 1));
    if (scale == 0) continue;
    double h = 0;
    for (std::size_t i = high + 1; i-- > m;) {
      ort[i] = H(i, m - 1) / scale;
      h += ort[i] * ort[i];
    }
    double g = std::sqrt(h);
    if (ort[m] > 0) g = -g;
    h -= ort[m] * g;
    ort[m] -= g;
    // H = (I − u uᵀ / h) H (I − u uᵀ / h).
    for (std::size_t j = m; j < n; j++) {
      double f = 0;
      for (std::size_t i = high + 1; i-- > m;) f += ort[i] * H(i, j);
      f /= h;
      for (std::size_t i = m; i <= high; i++) H(i, j) -= f * ort[i];
    }
    for (std::size_t i = 0; i <= high; i++) {
      double f = 0;
      for (std::size_t j = high + 1; j-- > m;) f += ort[j] * H(i, j);
      f /= h;
      for (std::size_t j = m; j <= high; j++) H(i, j) -= f * ort[j];
    }
    ort[m] *= scale;
    H(m, m - 1) = scale * g;
  }
  for (std::size_t i = 2; i < n; i++)
    for (std::size_t j = 0; j + 1 < i; j++) H(i, j) = 0;
}

// The eigenvalues of A (balanced, reduced to Hessenberg form, then the double-shift QR
// iteration). Throws std::runtime_error if an eigenvalue fails to deflate in 30 n iterations,
// which does not happen in practice.
inline std::vector<std::complex<double>> eigenvalues(Matrix A) {
  const int nn = static_cast<int>(A.n);
  std::vector<std::complex<double>> out(A.n);
  if (nn == 0) return out;
  balance(A);
  hessenberg(A);
  Matrix& H = A;
  std::vector<double> d(A.n, 0.0), e(A.n, 0.0);
  const int low = 0;
  int n = nn - 1;
  const double eps = std::ldexp(1.0, -52);
  double exshift = 0, p = 0, q = 0, r = 0, s = 0, z = 0, w, x, y;
  double norm = 0;
  for (int i = 0; i < nn; i++)
    for (int j = std::max(i - 1, 0); j < nn; j++) norm += std::fabs(H(i, j));
  int iter = 0, total = 0;
  while (n >= low) {
    // A single small subdiagonal element.
    int l = n;
    while (l > low) {
      s = std::fabs(H(l - 1, l - 1)) + std::fabs(H(l, l));
      if (s == 0) s = norm;
      if (std::fabs(H(l, l - 1)) < eps * s) break;
      l--;
    }
    if (l == n) {  // one root
      H(n, n) += exshift;
      d[n] = H(n, n);
      e[n] = 0;
      n--;
      iter = 0;
    } else if (l == n - 1) {  // two roots
      w = H(n, n - 1) * H(n - 1, n);
      p = (H(n - 1, n - 1) - H(n, n)) / 2;
      q = p * p + w;
      z = std::sqrt(std::fabs(q));
      H(n, n) += exshift;
      H(n - 1, n - 1) += exshift;
      x = H(n, n);
      if (q >= 0) {  // a real pair
        z = p >= 0 ? p + z : p - z;
        d[n - 1] = x + z;
        d[n] = z != 0 ? x - w / z : d[n - 1];
        e[n - 1] = e[n] = 0;
      } else {  // a complex pair
        d[n - 1] = d[n] = x + p;
        e[n - 1] = z;
        e[n] = -z;
      }
      n -= 2;
      iter = 0;
    } else {  // no convergence yet: form the shift
      x = H(n, n);
      y = 0;
      w = 0;
      if (l < n) {
        y = H(n - 1, n - 1);
        w = H(n, n - 1) * H(n - 1, n);
      }
      if (iter == 10) {  // Wilkinson's exceptional shift
        exshift += x;
        for (int i = low; i <= n; i++) H(i, i) -= x;
        s = std::fabs(H(n, n - 1)) + std::fabs(H(n - 1, n - 2));
        x = y = 0.75 * s;
        w = -0.4375 * s * s;
      }
      if (iter == 30) {  // MATLAB's exceptional shift (as JAMA has it)
        s = (y - x) / 2;
        s = s * s + w;
        if (s > 0) {
          s = std::sqrt(s);
          if (y < x) s = -s;
          s = x - w / ((y - x) / 2 + s);
          for (int i = low; i <= n; i++) H(i, i) -= s;
          exshift += s;
          x = y = w = 0.964;
        }
      }
      iter++;
      if (++total > 30 * nn) throw std::runtime_error("calx::linalg::eigenvalues: the QR iteration did not converge");
      // Two consecutive small subdiagonal elements.
      int m = n - 2;
      while (m >= l) {
        z = H(m, m);
        r = x - z;
        s = y - z;
        p = (r * s - w) / H(m + 1, m) + H(m, m + 1);
        q = H(m + 1, m + 1) - z - r - s;
        r = H(m + 2, m + 1);
        s = std::fabs(p) + std::fabs(q) + std::fabs(r);
        p /= s;
        q /= s;
        r /= s;
        if (m == l) break;
        if (std::fabs(H(m, m - 1)) * (std::fabs(q) + std::fabs(r)) <
            eps * (std::fabs(p) * (std::fabs(H(m - 1, m - 1)) + std::fabs(z) + std::fabs(H(m + 1, m + 1)))))
          break;
        m--;
      }
      for (int i = m + 2; i <= n; i++) {
        H(i, i - 2) = 0;
        if (i > m + 2) H(i, i - 3) = 0;
      }
      // The double QR step on rows l..n and columns m..n.
      for (int k = m; k <= n - 1; k++) {
        const bool notlast = k != n - 1;
        if (k != m) {
          p = H(k, k - 1);
          q = H(k + 1, k - 1);
          r = notlast ? H(k + 2, k - 1) : 0;
          x = std::fabs(p) + std::fabs(q) + std::fabs(r);
          if (x == 0) continue;
          p /= x;
          q /= x;
          r /= x;
        }
        s = std::sqrt(p * p + q * q + r * r);
        if (p < 0) s = -s;
        if (s != 0) {
          if (k != m) H(k, k - 1) = -s * x;
          else if (l != m) H(k, k - 1) = -H(k, k - 1);
          p += s;
          x = p / s;
          y = q / s;
          z = r / s;
          q /= p;
          r /= p;
          for (int j = k; j < nn; j++) {  // row modification
            p = H(k, j) + q * H(k + 1, j);
            if (notlast) {
              p += r * H(k + 2, j);
              H(k + 2, j) -= p * z;
            }
            H(k, j) -= p * x;
            H(k + 1, j) -= p * y;
          }
          for (int i = 0; i <= std::min(n, k + 3); i++) {  // column modification
            p = x * H(i, k) + y * H(i, k + 1);
            if (notlast) {
              p += z * H(i, k + 2);
              H(i, k + 2) -= p * r;
            }
            H(i, k) -= p;
            H(i, k + 1) -= p * q;
          }
        }
      }
    }
  }
  for (std::size_t i = 0; i < A.n; i++) out[i] = {d[i], e[i]};
  return out;
}

// The largest modulus of the eigenvalues of A.
inline double spectral_radius(const Matrix& A) {
  double rho = 0;
  for (const auto& l : eigenvalues(A)) rho = std::max(rho, std::abs(l));
  return rho;
}

}  // namespace calx::linalg
