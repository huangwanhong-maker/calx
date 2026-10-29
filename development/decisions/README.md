# Decision Log

Each record captures one decision that shapes the structure, interfaces or tooling of calx.

- Records are numbered and never rewritten.
- When a decision changes, write a new record that supersedes the old one, and point the old record's status to the new one.

## 1. Index

| No. | Decision | Status |
|---|---|---|
| [0001](0001-record-architecture-decisions.md) | Record architecture decisions | Accepted |
| [0002](0002-scope-and-boundaries.md) | Scope and boundaries | Accepted |
| [0003](0003-retargetable-computation-stack.md) | Retargetable computation stack | Accepted |
| [0004](0004-first-modules-and-the-build.md) | The first modules and the build | Accepted, provisionally |
| [0005](0005-forward-differentiation.md) | Forward differentiation by truncated Taylor arithmetic | Accepted, provisionally |
| [0006](0006-lanes-task-graphs-and-the-genome.md) | Lanes, task graphs, and what the genome asks of calx | Accepted |
| [0007](0007-conjugate-gradients.md) | Conjugate gradients on compressed sparse rows | Accepted |
| [0008](0008-optimisation-and-least-squares.md) | Optimisation in several variables, linear programming, and the matrices of statistics | Accepted |
| [0009](0009-real-time-scheduling.md) | Real-time scheduling: tasks with deadlines, their analyses, servers, simulation and dispatch | Accepted, provisionally |
| [0010](0010-kernel-awareness-and-extended-range-legendre.md) | The kernel's real-time classes, asked; Legendre functions in extended range | Accepted, provisionally |

## 2. Template

```markdown
# NNNN. Title

- **Status:** Proposed | Accepted | Superseded by NNNN
- **Date:** YYYY-MM-DD
- **Decided by:** …

## 1. Context
## 2. Decision
## 3. Alternatives Considered
## 4. Consequences
```
