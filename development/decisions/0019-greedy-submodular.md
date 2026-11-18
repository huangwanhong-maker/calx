# 0019. Greedy Maximisation of Submodular Functions

- **Status:** Accepted, provisionally (2026-11-18).
- **Date:** 2026-11-18
- **Decided by:** Claude, for the infrastructure's urban and disaster risk (its ADR 0129): where to open shelters so that the most people reach one in time, and which roads to reinforce within a budget.

## 1. Context

Choosing sites from a list so that the people within reach of them are the most is maximising a monotone submodular set function: coverage, a facility location's best similarities, and the greatest flow into the sinks chosen (the rank of a polymatroid) are all of this kind. The best set is hard to find; the greedy rule is within 1 − 1/e of it, and an instance gives a much tighter bound of its own. The markets domain's `space` chooses sites greedily for one function of its own; the method names no market, no shelter and no road, and so belongs here.

## 2. Decision

`optimisation/submodular.hpp`:

- `greedy_submodular(n, f, k)`: the greedy set of at most k elements with lazy evaluation (Minoux 1978), and the bound on the best set of k, the least over the greedy sets S of f(S) plus the k greatest gains with respect to S (Nemhauser, Wolsey and Fisher 1978), stale gains being upper bounds by submodularity;
- `greedy_submodular_budget(n, f, cost, budget)`: the better of greedy by gain and by gain per cost, each taking what still fits, at least (1 − 1/e)/2 of the best within the budget (Leskovec et al. 2007), its bound f(S) plus a fractional knapsack of gains within the budget.

f is any callable on sets of indices, asked for whole sets. That f is monotone and submodular is the caller's to know; the guarantees and the bound hold only then, and the header says so.

## 3. Consequences

`check_calx_submodular` (7 checks) holds the guarantees and the bound against every set listed, on random coverages, facility locations and budgets. The markets domain may move its own greedy choice onto this one; that is its decision.
