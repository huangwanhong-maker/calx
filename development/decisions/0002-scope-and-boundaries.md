# 0002. Scope and Boundaries

- **Status:** Accepted
- **Date:** 2026-10-21
- **Decided by:** the user proposed separate libraries. Claude drew the boundaries, recorded in [laboratory ADR 0008](https://github.com/huangwanhong-maker/Integrative_Astronomy_Geography_Oceanography_Laboratory/blob/main/development/decisions/0008-independent-libraries.md).

## 1. Context

The laboratory first held the computation stack and the numerical methods inside its own `core/`. These parts know nothing about astronomy, oceans or the atmosphere. They are also useful to other work, and they change on a different rhythm from the scientific code. The laboratory therefore extracted them into calx, alongside two sibling libraries: morphe (structures) and physis (geometry and physics).

## 2. Decision

1. **Contents.** calx contains:
   - the retargetable computation stack (`compute/`);
   - general numerical methods on floating-point numbers and arrays (`numerics/`).
2. **Exclusions.**
   - Exact arithmetic and algebraic, topological or categorical structures belong to morphe.
   - Geometry, units and physical laws belong to physis.
   - Discretisations that rely on such structure belong to physis. Examples are discrete exterior calculus and Lie-group integrators.
3. **Dependencies.** calx depends on no other project library. In particular, calx and morphe are independent of each other.
4. **Language.** calx is written in C++, following [laboratory ADR 0003](https://github.com/huangwanhong-maker/Integrative_Astronomy_Geography_Oceanography_Laboratory/blob/main/development/decisions/0003-cpp-computational-core.md).

## 3. Alternatives Considered

- **Including exact arithmetic.** Rejected. Rings and fields are algebraic structures, and keeping them in morphe lets each of calx and morphe be used without the other.
- **Including PDE discretisation in general.** Rejected. Discrete exterior calculus needs topology from morphe, which would make calx depend on morphe. Such discretisation belongs where structures meet physics, which is physis.

## 4. Consequences

- physis depends on calx.
- Any feature that would need both calx and morphe goes through an optional adapter, decided by its own record. An example is homology computed at scale on a GPU.
- A dependency of calx on morphe would need a new record (see the [equality-saturation idea](../ideas/2026-10-21-equality-saturation-in-the-compiler.md)).
