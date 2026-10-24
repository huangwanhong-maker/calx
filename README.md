# calx

*A retargetable numerical computing library.* The name comes from *calx*, the Latin word for the small stone once used for counting, from which "calculus" derives.

## 1. Purpose

calx lets a numerical computation be written once and run on the CPU today, and on other devices later, including devices that do not yet exist. It provides two things:

- **a computation stack.** It keeps six things separate: representation, schedule, compiler, runtime, task scheduler and device backends. The design follows TVM, Halide and MLIR.
- **general numerical methods**, built on that stack. These cover floating-point number systems, linear algebra, random numbers, approximation, special functions, transforms, signal processing, quadrature, integrators, solvers, optimisation, and the verification of results.

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
| [`numerics/`](numerics/) | Number systems, linear algebra, random numbers, approximation, special functions, transforms, signal processing, quadrature, integrators, solvers, optimisation, verification |
| [`docs/`](docs/) | Architecture |
| [`development/`](development/) | Decisions and ideas |

## 4. Relation to Sibling Libraries

| Library | Scope | Relation to calx |
|---|---|---|
| [morphe](https://github.com/huangwanhong-maker/morphe) | Mathematical structures and their equivalences | Independent of calx |
| [physis](https://github.com/huangwanhong-maker/physis) | Geometry and theoretical physics | Depends on calx |

calx depends on no other library of the project.

## 5. Status

**Draft (2026-10-23).** The first modules are written in C++20 and checked, as the neuroscience domain's spiking networks need them ([decision 0004](development/decisions/0004-first-modules-and-the-build.md)): `random`, `linalg` (dense, the exponential, compressed rows), `integrators`, `special_functions` (erfcx), `differentiation` (jets of order 1 and 2, [decision 0005](development/decisions/0005-forward-differentiation.md)), `transforms` (the FFT of any length), `signal` (Welch's spectral density), `quadrature`, `solvers` (Brent), `verification`, and `compute/scheduler` (a pool of host threads). For the genome ([decision 0006](development/decisions/0006-lanes-task-graphs-and-the-genome.md)): the CPU backend's target description and vector lanes (`calx::cpu`), task graphs and dynamic loops in the scheduler, the symmetric eigenproblem and determinants, the gamma family of special functions, gamma, beta and binomial variates, and Brent's minimiser, the first content of `optimisation`. The rest of the stack and of the numerics is designed and planned.

Build and check with CMake 3.20 or later:

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

## 6. License

calx is released under the MIT License. See [`LICENSE`](LICENSE).
