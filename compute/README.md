# compute: Retargetable Computation Stack

*Status: planned.* This is the lower layer of calx, and it depends on nothing. For the design, see [architecture §3–6](../docs/architecture.md) and [decision 0003](../development/decisions/0003-retargetable-computation-stack.md).

## 1. Purpose

The stack lets a computation, written once, run on the CPU today and on other devices later, including devices that do not yet exist. Moving to a new device never requires changing the code that uses the computation.

## 2. Components

| Directory | Role | Status |
|---|---|---|
| [`ir/`](ir/) | What is computed: a device-neutral representation | planned |
| [`schedule/`](schedule/) | How one computation is mapped onto hardware | planned |
| [`compiler/`](compiler/) | Passes and lowering | planned |
| [`runtime/`](runtime/) | Device API, memory, queues, modules | planned |
| [`scheduler/`](scheduler/) | When and where task graphs run | planned |
| [`backends/`](backends/) | One plug-in per target | planned |

## 3. Invariants

The library-wide invariants are in [CLAUDE.md §2](../CLAUDE.md). For this layer, they mean:

1. **No scientific concepts.** Nothing here refers to a scientific concept.
2. **One directory per device.** A new device means a new directory under `backends/`, and no other change.
3. **Reference results.** The CPU reference interpreter defines correct results.
4. **Declared precision.** Each computation declares its precision, and precision is never silently reduced.
5. **Versioned IR.** The IR text format is versioned.

## 4. External Systems

LLVM, MLIR, TVM, Kokkos or SYCL implementations may be used inside backends. None of them becomes part of the interfaces defined here.
