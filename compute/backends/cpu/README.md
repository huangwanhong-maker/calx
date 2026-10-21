# cpu: CPU Backend

*Status: planned. This is the first backend.*

## 1. Purpose

The CPU is the primary execution target. This backend also defines correct results for every other backend.

## 2. Planned Content

- **Reference interpreter.** It executes the IR directly and simply, with no optimisation. Clarity takes priority over speed.
- **Native path.** It compiles the IR for performance, using multithreading and SIMD. Whether it emits C++, LLVM IR or something else is decided when it is built.
