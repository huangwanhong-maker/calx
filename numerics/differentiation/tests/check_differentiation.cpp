// Checks of calx/differentiation: jets of order 2 against closed-form gradients and Hessians
// of the Rosenbrock function, of a composition of every elementary function, and of the
// power and arctangent functions of two jets; the symmetry of Hessians; dual numbers against
// the first-order part of jets; and Euler's theorem for a homogeneous function.

#include <calx/differentiation/jet.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <array>
#include <cmath>

namespace cd = calx::differentiation;

namespace {

template <std::size_t N>
double worst(const std::array<double, N>& a, const std::array<double, N>& b) {
  double w = 0;
  for (std::size_t i = 0; i < N; i++) w = std::max(w, std::fabs(a[i] - b[i]) / std::max(1.0, std::fabs(b[i])));
  return w;
}
template <std::size_t N>
double worst(const std::array<std::array<double, N>, N>& a, const std::array<std::array<double, N>, N>& b) {
  double w = 0;
  for (std::size_t i = 0; i < N; i++) w = std::max(w, worst(a[i], b[i]));
  return w;
}

}  // namespace

int main() {
  calx::verification::Suite s;

  // ---- Rosenbrock's function f = (1 - x)^2 + 100 (y - x^2)^2 ------------------------------------------------------------
  {
    auto f = [](const auto& z) { return (1.0 - z[0]) * (1.0 - z[0]) + 100.0 * (z[1] - z[0] * z[0]) * (z[1] - z[0] * z[0]); };
    const std::array<double, 2> x{-1.2, 1.0};
    const auto j = cd::taylor(f, x);
    const double a = x[0], b = x[1];
    const std::array<double, 2> g{-2 * (1 - a) - 400 * a * (b - a * a), 200 * (b - a * a)};
    const std::array<std::array<double, 2>, 2> H{{{2 - 400 * b + 1200 * a * a, -400 * a}, {-400 * a, 200}}};
    s.check("differentiation: Rosenbrock's function, value", j.v - f(x), 0);
    s.check("differentiation: Rosenbrock's gradient against its closed form (relative)", worst(j.d, g), 1e-15);
    s.check("differentiation: Rosenbrock's Hessian against its closed form (relative)", worst(j.h, H), 1e-15);
    s.check("differentiation: the gradient of a dual number equals the jet's", worst(cd::gradient(f, x), g), 1e-15);
  }
  // ---- a composition of every elementary function -----------------------------------------------------------------------
  {
    // f(x, y, z) = exp(sin x) * log(y) + sqrt(z) / cosh(x y) + atan(tanh(z) - tan(x)) + sinh(y) * pow(z, 1.5)
    auto f = [](const auto& w) {
      using std::atan; using std::cosh; using std::exp; using std::log; using std::pow; using std::sin; using std::sinh;
      using std::sqrt; using std::tan; using std::tanh;
      const auto& x = w[0];
      const auto& y = w[1];
      const auto& z = w[2];
      return exp(sin(x)) * log(y) + sqrt(z) / cosh(x * y) + atan(tanh(z) - tan(x)) + sinh(y) * pow(z, 1.5);
    };
    const std::array<double, 3> p{0.3, 1.7, 0.9};
    const auto j = cd::taylor(f, p);
    // Central differences of the closed-form gradient, computed independently.
    auto grad = [&](std::array<double, 3> q) {
      const double x = q[0], y = q[1], z = q[2];
      const double u = std::tanh(z) - std::tan(x), du = 1 / (1 + u * u);
      const double ch = std::cosh(x * y), sh = std::sinh(x * y);
      return std::array<double, 3>{
          std::cos(x) * std::exp(std::sin(x)) * std::log(y) - std::sqrt(z) * y * sh / (ch * ch) - du / (std::cos(x) * std::cos(x)),
          std::exp(std::sin(x)) / y - std::sqrt(z) * x * sh / (ch * ch) + std::cosh(y) * std::pow(z, 1.5),
          0.5 / (std::sqrt(z) * ch) + du * (1 - std::tanh(z) * std::tanh(z)) + 1.5 * std::sinh(y) * std::sqrt(z)};
    };
    s.check("differentiation: a composition of every elementary function, gradient against its closed form (relative)", worst(j.d, grad(p)), 1e-14);
    std::array<std::array<double, 3>, 3> H;
    const double h = 1e-5;
    for (int i = 0; i < 3; i++) {
      auto a = p, b = p;
      a[i] += h;
      b[i] -= h;
      const auto ga = grad(a), gb = grad(b);
      for (int k = 0; k < 3; k++) H[k][i] = (ga[k] - gb[k]) / (2 * h);
    }
    s.check("differentiation: its Hessian against central differences of the closed-form gradient (relative)", worst(j.h, H), 1e-8);
    double asym = 0;
    for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++) asym = std::max(asym, std::fabs(j.h[i][k] - j.h[k][i]));
    s.check("differentiation: the Hessian is symmetric", asym, 1e-14);
  }
  // ---- atan2, the power of two jets, and integer powers ----------------------------------------------------------------
  {
    // θ = atan2(y, x) is harmonic: θ_xx + θ_yy = 0, with ∇θ = (-y, x) / r².
    auto theta = [](const auto& w) { using std::atan2; return atan2(w[1], w[0]); };
    const std::array<double, 2> p{-0.8, 0.6};
    const auto j = cd::taylor(theta, p);
    s.check("differentiation: atan2's gradient (-y, x) / r^2", worst(j.d, std::array<double, 2>{-0.6, -0.8}), 1e-15);
    s.check("differentiation: atan2 is harmonic", j.h[0][0] + j.h[1][1], 1e-14);
    // f = x^y: f_x = y x^(y-1), f_y = x^y log x, f_xy = x^(y-1) (1 + y log x).
    auto f = [](const auto& w) { using std::pow; return pow(w[0], w[1]); };
    const std::array<double, 2> q{1.3, 2.4};
    const auto k = cd::taylor(f, q);
    const double x = q[0], y = q[1];
    s.check("differentiation: x^y, mixed derivative x^(y-1) (1 + y log x) (relative)", k.h[0][1] / (std::pow(x, y - 1) * (1 + y * std::log(x))) - 1, 1e-14);
    auto g = [](const auto& w) { using std::pow; return pow(w[0], 3) * pow(w[1], -2); };
    const auto m = cd::taylor(g, std::array<double, 2>{-1.5, 0.7});
    s.check("differentiation: x^3 y^-2 at negative x, d2/dx2 = 6 x y^-2 (relative)", m.h[0][0] / (6 * -1.5 / (0.7 * 0.7)) - 1, 1e-15);
  }
  // ---- Euler's theorem: a function homogeneous of degree k has x . ∇f = k f, and x^T H x = k (k - 1) f ----------------------
  {
    auto f = [](const auto& w) { using std::sqrt; return w[0] * w[1] * w[2] / sqrt(w[0] * w[0] + w[1] * w[1] + w[3] * w[3]) + w[3] * w[3]; };  // degree 2
    const std::array<double, 4> x{0.7, -1.1, 2.3, 0.4};
    const auto j = cd::taylor(f, x);
    double xg = 0, xHx = 0;
    for (int i = 0; i < 4; i++) {
      xg += x[i] * j.d[i];
      for (int k = 0; k < 4; k++) xHx += x[i] * j.h[i][k] * x[k];
    }
    s.check("differentiation: Euler's theorem, x . grad f = 2 f (relative)", xg / (2 * j.v) - 1, 1e-15);
    s.check("differentiation: and x^T H x = 2 f (relative)", xHx / (2 * j.v) - 1, 1e-14);
  }
  return s.done();
}
