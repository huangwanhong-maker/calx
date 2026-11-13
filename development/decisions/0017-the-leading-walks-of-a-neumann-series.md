# 0017. The Leading Walks of a Neumann Series

- **Status:** Accepted
- **Date:** 2026-11-13
- **Decided by:** the infrastructure's industrial-ecology domain (its ADR 0106), whose life cycle assessment and environmentally extended input–output analysis both ask which supply chains carry a product's impact. Both reduce to the same series: an impact eᵀ(I − Z)⁻¹y with Z ≥ 0 the requirements of one process or sector on another. The series names no object in nature, so Claude wrote it here, as the user asked for the domain.

## 1. Context

For a non-negative Z with ρ(Z) < 1, the Neumann series (I − Z)⁻¹ = Σ_k Z^k expands a total eᵀ(I − Z)⁻¹y into a sum over walks in the graph of Z's positive entries. Structural path analysis lists the largest of these walks (Defourny and Thorbecke 1984; Lenzen 2002; Peters and Hertwich 2006). Each walk is a supply chain, read from the final demand that starts it to the process or sector whose pressure it reaches.

The walks are infinitely many when the graph has a cycle, and grow as nᵏ with depth k, so they must be pruned. Pruning by depth alone misses large walks deep in the series. Pruning by a walk's own value misses the large walks that extend a small one.

## 2. Decision

- **In `linalg/nonnegative.hpp`**, `leading_walks(Z, e, y, cutoff, depth, keep, budget)` returns `Walks`:
  - `walks`, largest first;
  - `total = eᵀ(I − Z)⁻¹y`;
  - `covered`, the sum of every walk visited;
  - `by_depth`, eᵀZᵏy computed exactly by powers;
  - `visited`, the number of walks visited;
  - `cut`, the number of branches cut;
  - `complete`, whether the bound alone ended the search.
- **The bound.** A depth-first search carries each walk's weight v. Everything that extends a walk ending at j sums to v m_j, with m = eᵀ(I − Z)⁻¹. A branch is therefore cut only when v m_j is below cutoff × total, which is exact for non-negative data:
  - no walk above the cut is missed;
  - when the search is complete, total − covered < cut × cutoff × total.
- **Refusals.** Negative entries of Z, e or y are refused, since the bound would not hold. So is ρ(Z) ≥ 1, for which the series diverges.

## 3. Alternatives Considered

- **Pruning by depth, or by a walk's own value.** Rejected: either misses walks above the cut.
- **In the domain.** Rejected: LCA's technology matrix, a national input–output table and an absorbing Markov chain all ask the same question, and the method names none of them.
- **In morphe, as paths in a weighted graph.** The sums are of floating-point products with a numerical bound; morphe's walks are exact. The method stands with calx's other non-negative matrices.

## 4. Consequences

Nine checks:

- an acyclic chain is covered exactly, its total by hand;
- the series by depth sums to the total;
- the bounded search finds every walk above the cut that brute force finds to depth six;
- the walks kept are above the cut, largest first;
- what is covered never exceeds the total;
- negative weights are refused.

The industrial-ecology domain's checks also hold the bound total − covered < cut × cutoff × total.
