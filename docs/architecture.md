# Architecture

*Status: draft, 2026-10-21.* This document describes how calx is structured. The decisions behind it are in [development/decisions/](../development/decisions/).

## 1. Principles

### 1.1 Domain-Agnostic

calx knows only numbers, arrays and devices. It knows nothing about physics, units, the Earth or any other scientific concept. Other libraries hold everything else:

- exact arithmetic and algebraic structures belong to morphe;
- physical laws, and discretisations that preserve structure, belong to physis.

### 1.2 Separation of What from How

Following Halide, TVM and MLIR, calx keeps three things separate:

- a computation, which says *what* is computed;
- its schedule, which says *how* it is mapped onto hardware;
- the backend, which realises it on a device.

### 1.3 Conservative Interfaces

calx's own interfaces are small, versioned and stable: the IR, schedule primitives, the pass interface, the device API, the task-graph interface and the backend interface. Implementations behind them may use LLVM, MLIR, TVM, Kokkos or SYCL, and may be replaced without disturbing callers.

### 1.4 Correctness Before Speed

A reference interpreter defines correct results before any optimised path exists.

### 1.5 Growth by Need

The stack begins with an IR, a CPU interpreter and one kernel. It grows only as real computations demand.

## 2. Structure

```
calx
├── numerics/   numerical methods             (uses compute/ for bulk kernels)
└── compute/    retargetable computation stack (depends on nothing)
```

## 3. Computation Stack

```
computation (what) ──┐
                     ├──► compiler ──► backend code generator ──► executable module
schedule (how) ──────┘                                                  │
                                                                        ▼
task graph ──► scheduler ──► runtime (device API · memory · queues) ──► device
```

| Component | Responsibility | Stable contract |
|---|---|---|
| `ir/` | A device-neutral representation of computations. It has two levels: a dataflow of operators on typed arrays and fields, and a loop-nest IR beneath it. It also has a versioned, serialisable text form. | IR specification and text format |
| `schedule/` | Transformations that change how a computation runs but not what it computes: tiling, fusion, vectorisation, parallelisation, memory placement. Kept separate from the computation. | Schedule primitives |
| `compiler/` | Pass infrastructure: verification, lowering from the high-level to the low-level IR, target-independent optimisation, and hand-off to a backend | Pass interface |
| `runtime/` | A minimal, device-neutral API covering devices, buffers, queues, events, modules and kernel launch. It is exposed as a C ABI, and arrays are exchanged as DLPack. | Device API |
| `scheduler/` | Runs task graphs over one or more devices: tracks dependencies, moves data, decides placement | Task-graph interface |
| `backends/` | One directory per target. Each provides a target description, a code generator and a device-API implementation, and registers as a plug-in. | Backend interface |

*Remark.* `schedule/` maps one computation onto hardware, in the Halide/TVM sense. `scheduler/` decides when and where many computations run. The names are close, but the roles are distinct.

## 4. Adding a Device

To support a new device, add one directory under `compute/backends/` that implements three interfaces:

- a target description;
- a code generator;
- the device API.

Nothing else changes. The backend declares its capabilities. When a device cannot express a construct, the compiler either legalises the computation or rejects it with a clear error. It never silently changes results.

## 5. CPU Reference Backend

The CPU is the first and primary backend. It has two parts:

- **a reference interpreter.** It executes the IR directly and simply. It is the correctness oracle against which every other backend is tested.
- **a native path for performance.** Whether it emits C++, LLVM IR or something else is decided when it is built.

## 6. Precision Policy

- Each computation declares the precision it requires.
- Each backend declares which precisions it supports natively, and at what cost.
- The compiler may emulate higher precision, for example with double-double arithmetic or compensated summation.
- The compiler may never silently reduce precision.

## 7. Numerical Methods and the Stack

Not every computation passes through the stack.

- **Ordinary C++.** Scalar and small-scale algorithms in `numerics/` are ordinary C++ functions.
- **IR kernels.** Bulk data-parallel kernels that must run on more than one device are expressed in the IR. Examples are large linear-algebra operations, transforms and grid updates.
- **Migration.** A kernel may begin as plain C++. It moves into the IR when it first needs a second device.

## 8. Interfaces Offered to Dependent Libraries

- **C++ API.** Public headers live under `include/calx/`.
- **Runtime boundary.** The runtime exposes a C ABI, and arrays cross it as DLPack. NumPy, PyTorch and JAX also understand DLPack.
- **IR text format.** A stable artefact that other tools can produce or consume.
- **Python bindings.** Planned.

## 9. Module Layout

```
<module>/
├── README.md               purpose, status, dependencies, references, open questions
├── include/calx/<module>/  public headers
├── src/                    implementation
├── tests/                  unit tests
├── bindings/python/        Python bindings (optional)
└── benchmarks/             performance measurements (optional)
```
