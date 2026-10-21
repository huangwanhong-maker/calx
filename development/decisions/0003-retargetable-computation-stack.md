# 0003. Retargetable Computation Stack

- **Status:** Accepted
- **Date:** 2026-10-21
- **Decided by:** the user set the direction: a conservative, TVM-like separation, CPU first, open to future devices. Claude designed the components. This record carries forward [laboratory ADR 0004](https://github.com/huangwanhong-maker/Integrative_Astronomy_Geography_Oceanography_Laboratory/blob/main/development/decisions/0004-retargetable-computation-stack.md).

## 1. Context

The user mainly computes on CPUs, but wants computations to be deployable to many devices, including devices that do not yet exist. Hardware changes faster than scientific code. A design that names today's devices in its core abstractions will need rewriting with each new device.

## 2. Decision

calx's computation stack has six separate components ([architecture §3](../../docs/architecture.md)):

1. a device-neutral, versioned representation (`ir/`);
2. schedules, kept separate from computations (`schedule/`);
3. a compiler (`compiler/`);
4. a minimal runtime with a C ABI and DLPack arrays (`runtime/`);
5. a task scheduler (`scheduler/`);
6. backends as plug-ins (`backends/`).

The following rules also apply:

- **The CPU is the first backend.** Its reference interpreter is the correctness oracle.
- **The interfaces belong to calx.** External systems may be used inside backends, but never become the contract.
- **Scalar algorithms stay ordinary C++.** Only kernels that need several devices go through the IR.
- **Precision is declared per computation.** It is never silently reduced.

## 3. Alternatives Considered

- **Kokkos or SYCL as the abstraction.** Rejected as the contract, because it ties the representation to C++ templates and to the devices of today. Either may serve inside a backend.
- **Adopting TVM or MLIR wholesale.** Rejected as the contract. Both are large dependencies that keep changing, and TVM is oriented towards machine-learning tensors. Either may serve inside a backend.
- **Plain C++ only.** Rejected: there would be no path to new devices except rewriting.

## 4. Consequences

- The engineering effort is substantial, so the stack starts small: an IR, a CPU interpreter and one kernel.
- Every backend can be tested against the reference interpreter.
- A second backend tests whether the abstraction really holds.
