# compiler: Passes and Lowering

*Status: planned.*

## 1. Purpose

The compiler turns a computation and its schedule into target-specific IR. It then hands that IR to a backend's code generator.

## 2. Planned Content

- **Pass infrastructure:** a pass manager and a pass interface.
- **Verification passes.**
- **Lowering** from the high-level to the low-level IR.
- **Target-independent optimisation.**
- **Legalisation** against the capabilities a backend declares. The compiler either transforms the computation to fit those capabilities, or rejects it with a clear error. It never silently changes results.
