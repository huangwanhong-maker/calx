// Checks of the symplectic integrators (integrators/symplectic.hpp; decision 0012).
//
// - Störmer–Verlet on the harmonic oscillator: its map is linear with determinant 1, turns by θ with
//   cos θ = 1 − h²ω²/2, and conserves the modified energy p²/2m + (k/2)(1 − h²ω²/4) q² exactly
//   (from the map's matrix [[c, h/m], [−hk(1 − h²ω²/4), c]], c = 1 − h²ω²/2: for a matrix [[c, β], [−γ, c]]
//   with c² + βγ = 1 the form γq² + βp² is invariant).
// - Order 2 on the pendulum against Dormand–Prince at a tight tolerance; area preservation of the
//   nonlinear map; the energy of an eccentric Kepler orbit bounded over a thousand orbits, without drift.
// - RATTLE on the spherical pendulum: positions on the sphere and momenta tangent to it at every step;
//   the multipliers the normal force (m v²/r − m g ẑ·r̂)/(2r) for g(q) = |q|² − r², at order 2;
//   reversibility; the energy bounded.

#include <array>
#include <cmath>
#include <numbers>

#include <calx/integrators/runge_kutta.hpp>
#include <calx/integrators/symplectic.hpp>
#include <calx/verification/check.hpp>

using calx::verification::Suite;
namespace ci = calx::integrators;

int main() {
  Suite s;
  const double pi = std::numbers::pi;

  // ---- the harmonic oscillator --------------------------------------------------------------------------------
  {
    const double m = 2.0, k = 8.0, w = std::sqrt(k / m), h = 0.05;  // ω = 2, hω = 0.1
    const ci::Vec<1> inv{1 / m};
    const auto grad = [&](const ci::Vec<1>& q, ci::Vec<1>& g) { g[0] = k * q[0]; };
    ci::Vec<1> q{1.0}, p{0.0};
    const auto Eh = [&](double qq, double pp) { return pp * pp / (2 * m) + 0.5 * k * (1 - h * h * w * w / 4) * qq * qq; };
    const double E0 = Eh(q[0], p[0]);
    double worst = 0;
    const int steps = 100000;
    for (int i = 0; i < steps; i++) {
      ci::stormer_verlet<1>(inv, grad, q, p, h);
      worst = std::max(worst, std::fabs(Eh(q[0], p[0]) - E0) / E0);
    }
    s.check("verlet: the oscillator's modified energy p²/2m + (k/2)(1 − h²ω²/4) q² is conserved (10⁵ steps, relative)", worst, 1e-12);
    // The phase: (q, p) turns by θ each step, cos θ = 1 − h²ω²/2, along the invariant ellipse.
    const double theta = std::acos(1 - h * h * w * w / 2);
    const double scale = std::sqrt(0.5 * k * (1 - h * h * w * w / 4) / (1 / (2 * m)));  // p-axis over q-axis of the ellipse
    const double phase = std::atan2(-p[0] / scale, q[0]);
    const double expected = std::remainder(steps * theta, 2 * pi);
    s.check("verlet: after 10⁵ steps its phase is 10⁵ θ, cos θ = 1 − h²ω²/2 (rad)", std::fabs(std::remainder(phase - expected, 2 * pi)), 1e-8);
  }

  // ---- the pendulum: order, area -----------------------------------------------------------------------------
  {
    const ci::Vec<1> inv{1.0};
    const auto grad = [](const ci::Vec<1>& q, ci::Vec<1>& g) { g[0] = std::sin(q[0]); };
    const auto verlet_to = [&](double T, int n) {
      ci::Vec<1> q{1.2}, p{0.0};
      for (int i = 0; i < n; i++) ci::stormer_verlet<1>(inv, grad, q, p, T / n);
      return std::array<double, 2>{q[0], p[0]};
    };
    // The reference: Dormand–Prince on (q, p) at a tight tolerance.
    std::array<double, 2> x{1.2, 0.0};
    const auto f = [](const std::array<double, 2>& y, std::array<double, 2>& dy) { dy[0] = y[1]; dy[1] = -std::sin(y[0]); };
    double hs = 1e-3;
    for (int i = 0; i < 100; i++) ci::dopri5<2>(f, x, 0.1, {1.0, 1.0}, ci::Tolerance{1e-13, 1e-13}, hs);
    const auto e = [&](int n) { const auto y = verlet_to(10.0, n); return std::hypot(y[0] - x[0], y[1] - x[1]); };
    const double e1 = e(400), e2 = e(800), e3 = e(1600);
    s.check("verlet: second order on the pendulum, log₂ of the error's fall − 2 (two halvings)", std::max(std::fabs(std::log2(e1 / e2) - 2), std::fabs(std::log2(e2 / e3) - 2)), 0.05);
    // Area: the Jacobian of one step, by central differences of step 10⁻⁶, has determinant 1.
    const auto step = [&](double q0, double p0) { ci::Vec<1> q{q0}, p{p0}; ci::stormer_verlet<1>(inv, grad, q, p, 0.3); return std::array<double, 2>{q[0], p[0]}; };
    const double d = 1e-6;
    const auto a = step(1.0 + d, 0.4), b = step(1.0 - d, 0.4), c = step(1.0, 0.4 + d), dd = step(1.0, 0.4 - d);
    const double J = ((a[0] - b[0]) * (c[1] - dd[1]) - (a[1] - b[1]) * (c[0] - dd[0])) / (4 * d * d);
    s.check("verlet: one step of the pendulum keeps area (the Jacobian's determinant − 1, by differences)", std::fabs(J - 1), 1e-8);
  }

  // ---- Kepler: the energy over a thousand orbits ---------------------------------------------------------------
  {
    const ci::Vec<2> inv{1.0, 1.0};
    const auto grad = [](const ci::Vec<2>& q, ci::Vec<2>& g) { const double r3 = std::pow(q[0] * q[0] + q[1] * q[1], 1.5); g[0] = q[0] / r3; g[1] = q[1] / r3; };
    const double e = 0.6;
    ci::Vec<2> q{1 - e, 0.0}, p{0.0, std::sqrt((1 + e) / (1 - e))};
    const auto H = [&]() { return 0.5 * (p[0] * p[0] + p[1] * p[1]) - 1 / std::hypot(q[0], q[1]); };
    const double H0 = H(), period = 2 * pi;  // a = 1
    const int per = 2000;
    double first = 0, last = 0;
    for (int orbit = 0; orbit < 1000; orbit++)
      for (int i = 0; i < per; i++) {
        ci::stormer_verlet<2>(inv, grad, q, p, period / per);
        const double err = std::fabs(H() - H0);
        if (orbit < 100) first = std::max(first, err);
        if (orbit >= 900) last = std::max(last, err);
      }
    // Bounded, not drifting: a linear drift would make the last hundred orbits' error some nine times the first's.
    s.check("verlet: an orbit of e = 0.6 over 1000 periods: the energy's error stays near the first hundred orbits' (last over first − 1)", std::fabs(last / first - 1), 0.1);
  }

  // ---- RATTLE: the spherical pendulum ----------------------------------------------------------------------------
  {
    const double m = 1.5, r = 2.0, gz = 9.81;
    const ci::Vec<3> inv{1 / m, 1 / m, 1 / m};
    const auto grad = [&](const ci::Vec<3>&, ci::Vec<3>& g) { g = {0, 0, m * gz}; };
    const auto sphere = [&](const ci::Vec<3>& q, ci::Vec<1>& c, std::array<ci::Vec<3>, 1>& G) {
      c[0] = q[0] * q[0] + q[1] * q[1] + q[2] * q[2] - r * r;
      G[0] = {2 * q[0], 2 * q[1], 2 * q[2]};
    };
    const auto start = [&](ci::Vec<3>& q, ci::Vec<3>& p) {
      const double th = 0.9;  // from the downward vertical
      q = {r * std::sin(th), 0, -r * std::cos(th)};
      p = {0, m * 2.5, 0};
    };
    const auto energy = [&](const ci::Vec<3>& q, const ci::Vec<3>& p) { return (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]) / (2 * m) + m * gz * q[2]; };
    ci::Vec<3> q, p;
    start(q, p);
    const double E0 = energy(q, p);
    double on = 0, tangent = 0, dE = 0;
    for (int i = 0; i < 20000; i++) {
      const auto st = ci::rattle<3, 1>(inv, grad, sphere, q, p, 1e-3);
      on = std::max(on, std::fabs(std::hypot(q[0], q[1], q[2]) - r));
      tangent = std::max(tangent, std::fabs(q[0] * p[0] + q[1] * p[1] + q[2] * p[2]) / (r * m));
      dE = std::max(dE, std::fabs(energy(q, p) - E0) / std::fabs(E0));
      (void)st;
    }
    s.check("rattle: on the sphere at every step (20 000 steps, m)", on, 1e-12, "m");
    s.check("rattle: the momentum tangent to it, |q · p|/(r m)", tangent, 1e-12, "m/s");
    s.check("rattle: the energy bounded (relative, h = 1 ms)", dE, 1e-5);
    // The multipliers: the normal force (m v²/r − m g ẑ·r̂)/(2r), at order 2 in h.
    const auto multiplier_error = [&](double h) {
      ci::Vec<3> qq, pp;
      start(qq, pp);
      double worst = 0;
      for (int i = 0; i < int(std::lround(0.5 / h)); i++) {
        const ci::Vec<3> q0 = qq, p0 = pp;
        const auto st = ci::rattle<3, 1>(inv, grad, sphere, qq, pp, h);
        // At the step's midpoint, the speed and the height from the average of its ends.
        const double v2 = (p0[0] * p0[0] + p0[1] * p0[1] + p0[2] * p0[2] + pp[0] * pp[0] + pp[1] * pp[1] + pp[2] * pp[2]) / (2 * m * m);
        const double zr = (q0[2] + qq[2]) / (2 * r);
        const double expected = (m * v2 / r - m * gz * zr) / (2 * r);
        worst = std::max(worst, std::fabs(0.5 * (st.lambda[0] + st.mu[0]) - expected) / std::fabs(expected));
      }
      return worst;
    };
    const double e1 = multiplier_error(4e-3), e2 = multiplier_error(2e-3);
    s.check("rattle: its multipliers are the sphere's normal force, at order 2 (log₂ of the fall − 2)", std::fabs(std::log2(e1 / e2) - 2), 0.15);
    // Reversibility: a step forward and one back, the momentum reversed, return to the start.
    ci::Vec<3> qa, pa;
    start(qa, pa);
    const ci::Vec<3> q0 = qa, p0 = pa;
    ci::rattle<3, 1>(inv, grad, sphere, qa, pa, 0.01);
    for (auto& x : pa) x = -x;
    ci::rattle<3, 1>(inv, grad, sphere, qa, pa, 0.01);
    s.check("rattle: symmetric: a step and its reverse return to the start (m and kg m/s)",
            std::max({std::fabs(qa[0] - q0[0]), std::fabs(qa[1] - q0[1]), std::fabs(qa[2] - q0[2]), std::fabs(-pa[0] - p0[0]), std::fabs(-pa[1] - p0[1]), std::fabs(-pa[2] - p0[2])}), 1e-12);
  }
  return s.done();
}
