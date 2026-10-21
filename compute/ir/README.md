# ir: Computation Representation

*Status: planned.*

## 1. Purpose

The IR is the device-neutral description of *what* is computed. It is the contract between the scientific code and every backend.

## 2. Planned Content

- **High-level IR.** A dataflow graph of operators on typed arrays and fields. Each array or field records its shape, element type, precision and layout.
- **Low-level IR.** Loop nests, index expressions and memory scopes, produced by lowering the high-level IR.
- **Text format.** A versioned, serialisable form. A computation can then be stored, inspected, and compiled later by a backend that has not yet been written.
- **Verifier.** Checks types and shapes.

## 3. Open Questions

- **Typed annotations.** Should the IR offer a generic annotation mechanism, through which physis attaches physical dimensions and reference frames to types? Any such mechanism must keep calx domain-agnostic ([laboratory idea](https://github.com/huangwanhong-maker/Integrative_Astronomy_Geography_Oceanography_Laboratory/blob/main/development/ideas/2026-10-21-units-and-frames-as-types.md)).
- **Relation to MLIR.** Is the low-level IR an MLIR dialect, or does it lower to one?
