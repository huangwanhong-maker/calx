# 0004. The First Modules and the Build

- **Status:** Accepted, provisionally (the build settles §7's provisional choices until The Science Infrastructures' Phase 0 confirms them)
- **Date:** 2026-10-23
- **Decided by:** the user asked that mathematics common to the sciences live in the foundations, designed as a whole and written as the first simulation, the neuroscience domain's spiking networks, needs it; that mature libraries be borrowed where they already serve and C++ be written where it allows better optimisation and heterogeneous computing. Claude designed calx's part and wrote this record. The whole map is The Science Infrastructures' [ADR 0035](https://github.com/huangwanhong-maker/The-Science-Infrastructures/blob/main/development/decisions/0035-the-mathematical-foundations.md) and [map of the mathematics](https://github.com/huangwanhong-maker/The-Science-Infrastructures/blob/main/docs/mathematics.md).

## 1. Context

calx was designed (decisions 0002 and 0003) with no code. The spiking networks need random numbers reproducible across threads, the exponential of small matrices, one-step integrators, compressed sparse rows, parallel phases on the CPU, and, for their mean-field theory, the scaled complementary error function, adaptive quadrature and bracketed roots. Their checks must report as the infrastructure's other checks do.

## 2. Decision

1. **Three new modules in the plan of `numerics/`:** `random` (counter-based generators and sampling), `signal` (digital filters, windows, resampling, periodograms) and `optimisation` (split from `solvers`, which keeps roots, nonlinear systems, fixed points and continuation); and `verification` (checks against references and their report).
2. **Written now, header-only, status `draft`:** `numerics/random` (Philox4x32-10 and the common distributions), `numerics/linalg` (small dense matrices, LU solve, the exponential by Padé 13 with scaling and squaring, compressed sparse rows), `numerics/integrators` (Euler, Heun, RK4, Dormand–Prince 5(4)), `numerics/special_functions` (erfcx), `numerics/quadrature` (adaptive Gauss–Kronrod 7–15), `numerics/solvers` (Brent), `numerics/verification`, and `compute/scheduler` (a pool of host threads for bulk-synchronous phases, the scheduler's simplest placement). Everything else stays planned.
3. **Layout and names:** a module is `<area>/<module>/` with `include/calx/<module>/`, the namespace `calx::<module>`, and the CMake target `calx::<module>`.
4. **The build:** CMake 3.20 or later and C++20 (no compiler extensions), one project with a target per module; built on its own, calx builds its checks (`CALX_CHECKS`); added to a larger project, it builds them when asked. Checks are programs on `calx::verification`, made by `calx_add_check` and registered with CTest.
5. **Our own implementations,** because each need here is of a kind a general library does not meet: matrices of order below twenty exponentiated once and applied millions of times; random draws identified by purpose and entity, independent of threads; a report shared with the infrastructure's checks. Larger dense linear algebra, transforms and stiff systems may borrow Eigen, a BLAS, pocketfft or SUNDIALS behind these interfaces when they are needed.

## 3. Alternatives Considered

- **Writing the numerics inside the neuroscience domain** and moving them later: rejected by the user's direction, since every simulating domain would need them.
- **A borrowed generator** (std::mt19937, PCG): sequential generators make results depend on the order of draws, and therefore on the number of threads; Philox's draws are functions of their counters.
- **A borrowed test framework** (doctest, Catch2): the checks must print values against limits as the infrastructure's harness does, so that one runner adds them up.
- **Compiled libraries now:** the modules are templates and small inline functions; `src/` comes when code is worth compiling once.

## 4. Consequences

- calx has code, 7 check programs and 61 checks, run by the infrastructure's runner through CMake (locally, or in WSL on Windows).
- The computation stack (IR, schedules, compiler, runtime, backends) is untouched and still planned; the scheduler's thread pool is the CPU placement it will generalise, and kernels written in plain C++ now move into the IR when they first need a second device (architecture §7).
