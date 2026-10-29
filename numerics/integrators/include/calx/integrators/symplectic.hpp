// Symplectic integrators of separable Hamiltonian systems, free and with holonomic constraints
// (decision 0012).
//
// A separable system has H(q, p) = ½ pᵀ M⁻¹ p + V(q), with a diagonal mass M (one mass for each
// coordinate) and a potential whose gradient the caller computes. Störmer–Verlet (the velocity
// Verlet form) is
//
//   p_{k+½} = p_k − (h/2) ∇V(q_k),   q_{k+1} = q_k + h M⁻¹ p_{k+½},   p_{k+1} = p_{k+½} − (h/2) ∇V(q_{k+1}):
//
// symplectic, symmetric and of order 2, with an energy error that stays bounded for exponentially
// long times on analytic systems (Hairer, Lubich and Wanner 2006, Geometric Numerical Integration,
// 2nd ed., §I.1.4, §VI.3, §IX.8). It is the variational integrator of the trapezoidal discrete
// Lagrangian (Marsden and West 2001, Acta Numerica 10, §2.6), which physis's discrete mechanics
// checks.
//
// With m holonomic constraints g(q) = 0 whose Jacobian G = ∂g/∂q has full rank, RATTLE (Andersen 1983;
// Leimkuhler and Skeel 1994; Hairer, Lubich and Wanner 2006, §VII.1.4) is
//
//   p_{k+½} = p_k − (h/2)(∇V(q_k) + G(q_k)ᵀ λ),   q_{k+1} = q_k + h M⁻¹ p_{k+½},   g(q_{k+1}) = 0,
//   p_{k+1} = p_{k+½} − (h/2)(∇V(q_{k+1}) + G(q_{k+1})ᵀ μ),   G(q_{k+1}) M⁻¹ p_{k+1} = 0:
//
// symplectic on the constraint manifold, symmetric and of order 2; positions stay on g = 0 and momenta
// on its tangent space to the tolerance of the solve. λ is found by Newton's method on g(q_{k+1}(λ)) = 0,
// whose Jacobian is −(h²/2) G(q_{k+1}) M⁻¹ G(q_k)ᵀ; μ by the linear system (h/2) G M⁻¹ Gᵀ μ = G M⁻¹ p*,
// p* = p_{k+½} − (h/2) ∇V(q_{k+1}). The multipliers are the constraint's forces: −Gᵀλ over the first half
// step and −Gᵀμ over the second, a body's normal force on the surface it is held to.
//
// A force that is not the gradient of a potential (a motor, friction) enters as part of "∇V" with the
// opposite sign; one that depends on the momentum is the caller's to split off (an exact damping step on
// either side, for instance), since folding it into the half steps would lose the symmetry.

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>

#include <calx/integrators/runge_kutta.hpp>  // Vec
#include <calx/linalg/small.hpp>

namespace calx::integrators {

// One Störmer–Verlet step. gradient(q, g) writes ∇V(q) into g.
template <std::size_t N, class Gradient>
inline void stormer_verlet(const Vec<N>& inverse_mass, Gradient&& gradient, Vec<N>& q, Vec<N>& p, double h) {
  Vec<N> g;
  gradient(q, g);
  for (std::size_t i = 0; i < N; i++) p[i] -= 0.5 * h * g[i];
  for (std::size_t i = 0; i < N; i++) q[i] += h * inverse_mass[i] * p[i];
  gradient(q, g);
  for (std::size_t i = 0; i < N; i++) p[i] -= 0.5 * h * g[i];
}

// The outcome of a RATTLE step: the multipliers of its two halves, the Newton iterations of the first,
// and the constraint's residual at the end, |g(q_{k+1})|∞.
template <std::size_t M>
struct RattleStep {
  Vec<M> lambda{}, mu{};
  std::size_t iterations = 0;
  double residual = 0;
};


// One RATTLE step. gradient(q, g) writes ∇V(q); constraint(q, c, G) writes g(q) into c (M values) and its
// Jacobian into G (M rows of N). The step's first half is solved to |g| ≤ tolerance.
template <std::size_t N, std::size_t M, class Gradient, class Constraint>
inline RattleStep<M> rattle(const Vec<N>& inverse_mass, Gradient&& gradient, Constraint&& constraint, Vec<N>& q, Vec<N>& p, double h,
                            double tolerance = 1e-13, std::size_t max_iterations = 50) {
  using Jacobian = std::array<Vec<N>, M>;
  RattleStep<M> out;
  Vec<N> g;
  Vec<M> c;
  Jacobian G0, G1;
  gradient(q, g);
  constraint(q, c, G0);
  // The unconstrained half kick and drift, from which the constraint force is taken away.
  Vec<N> half, free;
  for (std::size_t i = 0; i < N; i++) half[i] = p[i] - 0.5 * h * g[i];
  for (std::size_t i = 0; i < N; i++) free[i] = q[i] + h * inverse_mass[i] * half[i];
  Vec<M> lambda{};
  Vec<N> q1 = free;
  for (std::size_t it = 0;; it++) {
    for (std::size_t i = 0; i < N; i++) {
      double s = 0;
      for (std::size_t a = 0; a < M; a++) s += G0[a][i] * lambda[a];
      q1[i] = free[i] - 0.5 * h * h * inverse_mass[i] * s;
    }
    constraint(q1, c, G1);
    double worst = 0;
    for (std::size_t a = 0; a < M; a++) worst = std::max(worst, std::fabs(c[a]));
    out.iterations = it;
    if (worst <= tolerance) break;
    if (it >= max_iterations) throw std::runtime_error("calx::integrators::rattle: the constraint was not met within the iterations allowed");
    // Newton on λ: ∂g(q₁)/∂λ = −(h²/2) G(q₁) M⁻¹ G(q_k)ᵀ.
    linalg::SmallMatrix<M> J{};
    for (std::size_t a = 0; a < M; a++)
      for (std::size_t b = 0; b < M; b++) {
        double s = 0;
        for (std::size_t i = 0; i < N; i++) s += G1[a][i] * inverse_mass[i] * G0[b][i];
        J[a][b] = -0.5 * h * h * s;
      }
    Vec<M> minus;
    for (std::size_t a = 0; a < M; a++) minus[a] = -c[a];
    const Vec<M> d = linalg::solve<M>(J, minus);
    for (std::size_t a = 0; a < M; a++) lambda[a] += d[a];
  }
  // The first half kick with its constraint force, and the drift.
  for (std::size_t i = 0; i < N; i++) {
    double s = 0;
    for (std::size_t a = 0; a < M; a++) s += G0[a][i] * lambda[a];
    half[i] -= 0.5 * h * s;
  }
  q = q1;
  // The second half kick, its momentum projected onto the constraint's tangent space.
  gradient(q, g);
  Vec<N> star;
  for (std::size_t i = 0; i < N; i++) star[i] = half[i] - 0.5 * h * g[i];
  linalg::SmallMatrix<M> A{};
  Vec<M> rhs{};
  for (std::size_t a = 0; a < M; a++) {
    for (std::size_t i = 0; i < N; i++) rhs[a] += G1[a][i] * inverse_mass[i] * star[i];
    for (std::size_t b = 0; b < M; b++) {
      double s = 0;
      for (std::size_t i = 0; i < N; i++) s += G1[a][i] * inverse_mass[i] * G1[b][i];
      A[a][b] = 0.5 * h * s;
    }
  }
  const Vec<M> mu = linalg::solve<M>(A, rhs);
  for (std::size_t i = 0; i < N; i++) {
    double s = 0;
    for (std::size_t a = 0; a < M; a++) s += G1[a][i] * mu[a];
    p[i] = star[i] - 0.5 * h * s;
  }
  out.lambda = lambda;
  out.mu = mu;
  for (std::size_t a = 0; a < M; a++) out.residual = std::max(out.residual, std::fabs(c[a]));
  return out;
}

}  // namespace calx::integrators
