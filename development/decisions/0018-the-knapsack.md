# 0018. The 0/1 Knapsack Problem by Branch and Bound, with Its Certificate

- **Status:** Accepted (2026-11-16)
- **Date:** 2026-11-16
- **Decided by:** Claude, for tyche's information bought under a budget (its decision 0026), at the infrastructure's request for decisions on evidence in space and time (its ADR 0123).

## 1. Context

tyche needed to choose which tests to run before several separate decisions, within a budget: a subset of greatest total worth whose costs fit, the 0/1 knapsack problem. It names no statistical object, so it is calx's (decision 0008: domain-free numerics are calx's).

## 2. Decision

`optimisation/knapsack.hpp`: `knapsack(value, weight, capacity, node_limit)` → `Knapsack{chosen, value, weight, bound, exact, nodes}`.

- Branch and bound, depth first in the order of the ratios v/w, each node bounded by Dantzig's relaxation (Horowitz and Sahni 1974; Martello and Toth 1990, ch. 2). The first leaf is the greedy solution.
- Items of no positive value are never taken, of no weight always taken, heavier than the capacity never.
- **A certificate.** A search that ends is exact, its bound its value. One stopped at the node limit returns its best with `exact` false and the root's Dantzig bound: the optimum lies between them.

## 3. Alternatives Considered

- **Dynamic programming over the capacity.** Exact and simple for whole weights, useless for real ones (costs in money, worths in harm) without scaling, which changes the problem. It is kept as a check.
- **calx's simplex with the integrality dropped.** That is Dantzig's bound, not a solution.
- **A general mixed-integer solver.** Not needed for one constraint; a dependency of that size waits until a problem needs several.

## 4. Consequences

`check_calx_knapsack` (14 checks): every subset of 500 problems up to 14 items, real and tied whole numbers; dynamic programming's value on 60 problems of up to 60 items; the certificate on subset-sum problems, where every ratio is equal and the bound prunes nothing; edges and refusals. calx 976 checks.
