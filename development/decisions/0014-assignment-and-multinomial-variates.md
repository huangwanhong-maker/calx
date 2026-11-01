# 0014. The Assignment Problem and Multinomial Variates

- **Status:** Accepted, provisionally (dense costs; multinomial variates by conditional binomials)
- **Date:** 2026-11-01
- **Decided by:** tyche's non-negative matrix factorisation (its decision 0006), written for the infrastructure's genome analyses (its ADR 0064: mutational signatures extracted from cancer genomes), clusters the components of many replicate factorisations so that each replicate gives one component to each cluster, and resamples a catalogue multinomially within its columns. Both steps name no statistical object: they are an optimisation over permutations and a standard distribution. Claude wrote them here, with this record, at the user's request.

## 1. Context

`optimisation` held minima in one and several variables and linear programmes (decisions 0006, 0008); `random` the common distributions, the binomial among them (decision 0006). The constrained clustering of signatures matches each replicate's k components to k centroids by the assignment of greatest total similarity (SigProfilerExtractor's use of the Hungarian algorithm). The linear programme of an assignment could be given to the simplex method, but at O(k⁶) and with degenerate vertices; the Hungarian method is O(k³) and returns a dual certificate of optimality. The bootstrap of Alexandrov et al. (2013) redraws each sample's mutations multinomially with its total kept.

## 2. Decision

- **`optimisation/assignment.hpp`.** `assignment(cost)` → `Assignment{column, cost, u, v}`: the least-cost injection of the rows of an n × m matrix (n ≤ m) into its columns, by the Hungarian method (Kuhn 1955; Munkres 1957) in its O(n² m) form by shortest augmenting paths over reduced costs (Burkard, Dell'Amico and Martello 2012, ch. 4). It returns the dual potentials, so that a caller can verify optimality: u_i + v_j ≤ c_ij everywhere, with equality on the assignment, and v_j = 0 on the columns left unassigned.
- **`random/multinomial.hpp`.** `multinomial(rng, n, p)`: counts of n draws among categories with probabilities proportional to p, by conditional binomials (the chain rule of the multinomial law), each by the Stream's exact binomial.

## 3. Alternatives Considered

- **The assignment as a linear programme** on `simplex`. Rejected: slower by orders of magnitude and without the structure's certificate.
- **The assignment in morphe**, as combinatorics. Rejected: morphe holds exact structures; an optimisation with floating-point costs is calx's.
- **Multinomial variates in tyche**, beside their use. Rejected: sampling the standard distributions is `random`'s (its README, §1), and the domains' simulators may need it too.
- **A member function of `Stream`.** Deferred: a free function over a Stream keeps `stream.hpp` unchanged; it may move into the class if more multivariate distributions join it.

## 4. Consequences

- **`check_calx_assignment`.** The least cost over every injection (brute force, 500 problems up to 6 × 8, real costs and tied integers); the dual certificate (feasibility, complementary slackness, v ≤ 0 off the assignment, strong duality) on problems up to 60 × 80; a hidden permutation recovered; refusals.
- **`check_calx_multinomial`.** The draws' totals; a category of probability 0 never drawn; means, variances and a covariance over 20 000 draws; Pearson's χ² of a draw of 10⁶ with unnormalised probabilities; refusals.
- **Interfaces.** Both are additive; nothing existing changes.
