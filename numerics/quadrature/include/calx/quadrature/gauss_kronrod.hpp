// Integration on an interval by the Gauss-Kronrod pair (7, 15), globally adaptive.
//
// On [a, b] the 15-point Kronrod rule K15 extends the 7-point Gauss rule G7 by eight nodes
// (Kronrod 1965); K15 integrates polynomials of degree 22 exactly and |K15 - G7| estimates
// its error. The adaptive scheme keeps the subintervals with their estimates and bisects
// the worst until the summed estimate meets max(absolute, relative |I|) or the number of
// subintervals reaches its limit (Piessens et al. 1983, QUADPACK, QAG without
// extrapolation). The estimate |K15 - G7| is conservative: QUADPACK's scaled estimate
// (200 |K - G|)^1.5 is sharper for smooth integrands, and may be adopted later.
//
// Nodes and weights from QUADPACK's qk15, to 33 digits.

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <queue>
#include <vector>

namespace calx::quadrature {

struct Estimate {
  double value = 0, error = 0;
  int evaluations = 0, intervals = 0;
  bool converged = false;
};

namespace detail {

inline constexpr std::array<double, 8> XGK{
    0.991455371120812639206854697526329, 0.949107912342758524526189684047851, 0.864864423359769072789712788640926,
    0.741531185599394439863864773280788, 0.586087235467691130294144845693013, 0.405845151377397166906606412076961,
    0.207784955007898467600689403773245, 0.000000000000000000000000000000000};
inline constexpr std::array<double, 8> WGK{
    0.022935322010529224963732008058970, 0.063092092629978553290700663189204, 0.104790010322250183839876322541518,
    0.140653259715525918745189590510238, 0.169004726639267902826583426598550, 0.190350578064785409913256402421014,
    0.204432940075298892414161999234649, 0.209482141084727828012999174891714};
// Gauss weights of the nodes XGK[1], XGK[3], XGK[5], XGK[7].
inline constexpr std::array<double, 4> WG{
    0.129484966168869693270611432679082, 0.279705391489276667901467771423780, 0.381830050505118944950369775488975,
    0.417959183673469387755102040816327};

}  // namespace detail

// K15 and G7 on [a, b]: {kronrod, |kronrod - gauss|}.
template <class F>
inline std::array<double, 2> gauss_kronrod_15(F&& f, double a, double b) {
  const double c = 0.5 * (a + b), h = 0.5 * (b - a);
  const double fc = f(c);
  double k = fc * detail::WGK[7], g = fc * detail::WG[3];
  for (int j = 0; j < 7; j++) {
    const double dx = h * detail::XGK[j];
    const double s = f(c - dx) + f(c + dx);
    k += detail::WGK[j] * s;
    if (j % 2 == 1) g += detail::WG[j / 2] * s;
  }
  return {k * h, std::fabs((k - g) * h)};
}

template <class F>
inline Estimate integrate(F&& f, double a, double b, double absolute = 1e-12, double relative = 1e-10, int max_intervals = 1000) {
  struct Piece { double a, b, value, error; };
  auto cmp = [](const Piece& x, const Piece& y) { return x.error < y.error; };
  std::priority_queue<Piece, std::vector<Piece>, decltype(cmp)> heap(cmp);
  Estimate out;
  auto [v0, e0] = gauss_kronrod_15(f, a, b);
  heap.push({a, b, v0, e0});
  out.evaluations = 15;
  double value = v0, error = e0;
  while (error > std::max(absolute, relative * std::fabs(value)) && static_cast<int>(heap.size()) < max_intervals) {
    const Piece p = heap.top();
    heap.pop();
    const double m = 0.5 * (p.a + p.b);
    auto [v1, e1] = gauss_kronrod_15(f, p.a, m);
    auto [v2, e2] = gauss_kronrod_15(f, m, p.b);
    out.evaluations += 30;
    value += v1 + v2 - p.value;
    error += e1 + e2 - p.error;
    heap.push({p.a, m, v1, e1});
    heap.push({m, p.b, v2, e2});
  }
  // Sum afresh, so that the running updates leave no rounding behind.
  value = 0;
  error = 0;
  out.intervals = static_cast<int>(heap.size());
  while (!heap.empty()) { value += heap.top().value; error += heap.top().error; heap.pop(); }
  out.value = value;
  out.error = error;
  out.converged = error <= std::max(absolute, relative * std::fabs(value));
  return out;
}

}  // namespace calx::quadrature
