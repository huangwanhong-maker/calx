# schedule: Mapping a Computation onto Hardware

*Status: planned.*

## 1. Purpose

A schedule is a set of transformations that change how a computation runs without changing what it computes. Following Halide and TVM, a schedule is a separate object that is applied to a computation. One computation can therefore have many schedules, for different devices or problem sizes.

## 2. Planned Content

- **Primitives:**
  - split, tile, reorder and fuse;
  - vectorise, parallelise and unroll;
  - cache and memory placement;
  - precision changes, within the tolerance that the computation declares.
- **Semantic checks.** Every primitive is checked to preserve the computation's meaning.
- **Later: automatic scheduling** guided by cost models.

*Remark.* `schedule/` is not [`scheduler/`](../scheduler/), which decides when and where many computations run.
