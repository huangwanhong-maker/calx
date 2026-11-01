// The linear assignment problem: an injection of rows into columns of least total cost.
//
// THE PROBLEM. Given costs c_ij for n rows and m ≥ n columns, find σ: rows → columns, injective,
// minimising Σ_i c_{iσ(i)}. Its linear programme (rows assigned once, columns at most once, x ≥ 0)
// has integral vertices, and its dual is: maximise Σ_i u_i + Σ_j v_j subject to u_i + v_j ≤ c_ij,
// v_j ≤ 0 for the columns left unassigned when m > n. An assignment is optimal exactly when some
// feasible (u, v) is tight on it (complementary slackness) — the certificate returned with it.
//
// THE ALGORITHM. The Hungarian method (Kuhn 1955, "The Hungarian method for the assignment problem",
// Naval Res. Logist. Q. 2, 83–97; Munkres 1957, J. SIAM 5, 32–38) in its O(n² m) form by shortest
// augmenting paths over reduced costs c_ij − u_i − v_j ≥ 0 (Burkard, Dell'Amico and Martello 2012,
// Assignment Problems, rev. ed., SIAM, ch. 4; the dual-potential form of Jonker and Volgenant 1987,
// Computing 38, 325–340, without their initialisation heuristics): rows are added one at a time, and
// each is placed by a Dijkstra search for the cheapest alternating path from it to a free column,
// the potentials updated so that every reduced cost stays non-negative and the matched ones zero.
// Exact in arithmetic; in floating point the cost is optimal to within about n ε max |c_ij|. Costs
// must be finite. Maximising a score s_ij is minimising −s_ij.

#pragma once

#include <calx/linalg/rectangular.hpp>

#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace calx::optimisation {

struct Assignment {
  std::vector<std::size_t> column;  // σ(i): the column assigned to row i
  double cost = 0;                  // Σ_i c_{iσ(i)}
  std::vector<double> u, v;         // dual potentials: u_i + v_j ≤ c_ij, with equality on the assignment
};

// The least-cost assignment of the rows of `cost` (n × m, n ≤ m) to distinct columns.
inline Assignment assignment(const calx::linalg::Rectangular& cost) {
  const std::size_t n = cost.rows, m = cost.cols;
  if (n > m) throw std::invalid_argument("calx::optimisation::assignment: no more rows than columns (transpose the costs)");
  for (const double c : cost.a)
    if (!std::isfinite(c)) throw std::invalid_argument("calx::optimisation::assignment: finite costs");
  const double inf = std::numeric_limits<double>::infinity();
  // 1-based potentials; p[j] is the row (1..n) matched to column j (0: none); column 0 is the root.
  std::vector<double> u(n + 1, 0.0), v(m + 1, 0.0), minv(m + 1);
  std::vector<std::size_t> p(m + 1, 0), way(m + 1, 0);
  std::vector<char> used(m + 1);
  for (std::size_t i = 1; i <= n; i++) {
    p[0] = i;
    std::size_t j0 = 0;
    std::fill(minv.begin(), minv.end(), inf);
    std::fill(used.begin(), used.end(), 0);
    do {
      used[j0] = 1;
      const std::size_t i0 = p[j0];
      double delta = inf;
      std::size_t j1 = 0;
      for (std::size_t j = 1; j <= m; j++) {
        if (used[j]) continue;
        const double reduced = cost(i0 - 1, j - 1) - u[i0] - v[j];
        if (reduced < minv[j]) {
          minv[j] = reduced;
          way[j] = j0;
        }
        if (minv[j] < delta) {
          delta = minv[j];
          j1 = j;
        }
      }
      for (std::size_t j = 0; j <= m; j++) {
        if (used[j]) {
          u[p[j]] += delta;
          v[j] -= delta;
        } else {
          minv[j] -= delta;
        }
      }
      j0 = j1;
    } while (p[j0] != 0);
    // Augment along the path back to the root.
    do {
      const std::size_t j1 = way[j0];
      p[j0] = p[j1];
      j0 = j1;
    } while (j0 != 0);
  }
  Assignment a;
  a.column.assign(n, 0);
  for (std::size_t j = 1; j <= m; j++)
    if (p[j] != 0) a.column[p[j] - 1] = j - 1;
  for (std::size_t i = 0; i < n; i++) a.cost += cost(i, a.column[i]);
  a.u.assign(u.begin() + 1, u.end());
  a.v.assign(v.begin() + 1, v.end());
  return a;
}

}  // namespace calx::optimisation
