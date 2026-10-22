# calx

*A retargetable numerical computing library.* The name comes from *calx*, the Latin word for the small stone once used for counting, from which "calculus" derives.

## 1. Purpose

calx lets a numerical computation be written once and run on the CPU today, and on other devices later, including devices that do not yet exist. It provides two things:

- **a computation stack.** It keeps six things separate: representation, schedule, compiler, runtime, task scheduler and device backends. The design follows TVM, Halide and MLIR.
- **general numerical methods**, built on that stack. These cover floating-point number systems, linear algebra, transforms, quadrature, integrators and solvers.

calx knows nothing about any scientific domain. It is one of the foundation libraries of [The Science Infrastructures](https://github.com/huangwanhong-maker/The-Science-Infrastructures), but it stands on its own.

## 2. Design Principles

- **Domain-agnostic.** calx works with numbers, arrays and devices, and knows nothing of physics, units or the Earth.
- **What is separate from how.**
  - A computation says *what* to compute.
  - A schedule says *how* the computation is mapped onto hardware.
  - A backend *realises* the result on a device.
- **New devices are additive.** Supporting a new device means adding one directory under `compute/backends/`. Nothing else changes.
- **Correctness before speed.** A simple CPU interpreter of the IR defines correct results, and every other backend is tested against it.
- **Declared precision.** Each computation declares the precision it needs, and calx never silently reduces it.
- **Own the interfaces, borrow the implementations.** LLVM, MLIR, TVM, Kokkos or SYCL implementations may be used inside backends, but never become part of calx's interfaces.

## 3. Organisation

| Path | Contents |
|---|---|
| [`compute/`](compute/) | The retargetable stack: `ir`, `schedule`, `compiler`, `runtime`, `scheduler`, `backends` |
| [`numerics/`](numerics/) | Number systems, linear algebra, approximation, special functions, transforms, quadrature, integrators, solvers |
| [`docs/`](docs/) | Architecture |
| [`development/`](development/) | Decisions and ideas |

## 4. Relation to Sibling Libraries

| Library | Scope | Relation to calx |
|---|---|---|
| [morphe](https://github.com/huangwanhong-maker/morphe) | Mathematical structures and their equivalences | Independent of calx |
| [physis](https://github.com/huangwanhong-maker/physis) | Geometry and theoretical physics | Depends on calx |

calx depends on no other library of the project.

## 5. Status

**Planned (2026-10-21).** Only the design documents exist so far; there is no code yet.

## 6. License

calx is released under the MIT License. See [`LICENSE`](LICENSE).
