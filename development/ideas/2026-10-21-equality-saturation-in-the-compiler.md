---
title: Equality Saturation in the Compiler
status: seed
date: 2026-10-21
related: [compute/compiler, compute/schedule]
---

# Equality Saturation in the Compiler

## 1. Idea

Compiler optimisation is a search for an equivalent program that runs better.

- **Equality saturation** keeps many equivalent programs at once in an *e-graph*, a congruence-closed set of terms. It then extracts the best program according to a cost model (egg; Willsey et al., 2021).
- **Rewrite rules become declarative.** They are applied until the e-graph is saturated.
- **Phase ordering disappears.** The question of which optimisation pass should run first no longer arises.

## 2. Relation to morphe

E-graphs, term rewriting, congruence closure and program equivalence (including bisimulation) are subjects of the structures library morphe. calx has two options:

- **Its own engine.** Implement a small rewriting engine inside calx, which keeps calx independent.
- **morphe's engine.** Depend on morphe for rewriting. This needs a new decision record, because calx and morphe are independent by design ([decision 0002](../decisions/0002-scope-and-boundaries.md)).

## 3. Question

Which rewrites does the compiler need first? Are they simple enough to implement locally?
