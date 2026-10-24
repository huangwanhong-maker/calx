# CLAUDE.md

This file guides Claude Code in the calx repository.

calx is one of the foundation libraries of [The Science Infrastructures](https://github.com/huangwanhong-maker/The-Science-Infrastructures), where it is included as the git submodule `foundations/calx`. It is nevertheless a project in its own right, and must make sense without the infrastructure.

## 1. Purpose and Scope

calx is a retargetable numerical computing library: a computation written once runs on the CPU and on other devices, including devices that do not yet exist. It contains:

- `compute/`, the computation stack:
  - representation (`ir/`);
  - schedule (`schedule/`);
  - compiler (`compiler/`);
  - runtime (`runtime/`);
  - task scheduler (`scheduler/`);
  - backends (`backends/`).
- `numerics/`, general numerical methods on floating-point numbers and arrays.

The following are out of scope, with where each belongs instead:

| Out of scope | Belongs to |
|---|---|
| Exact arithmetic and algebraic, topological or categorical structures | morphe |
| Geometry, units and physical laws | physis |
| Discretisations that rely on such structure, e.g. discrete exterior calculus, Lie-group integrators | physis |
| Anything about particular bodies, datasets or observations | the domains of The Science Infrastructures |

## 2. Invariants

1. calx depends on no other project library and knows nothing of any scientific domain.
2. Supporting a new device means adding one directory under `compute/backends/`. Nothing else changes.
3. The CPU reference interpreter defines correct results. Every other backend is tested against it.
4. Each computation declares its precision. Precision is never silently reduced.
5. The IR text format is versioned. Old versions remain readable, or are upgraded by an explicit migration.
6. External systems may be used inside backends, but never become part of calx's interfaces.

## 3. Design Philosophy

calx shares its governing principle with its sibling repositories: **modular, extensible, revisable**.

- **Status.** Every module carries one of these statuses:
  - `planned`: named, with no content yet;
  - `seed`: an idea, possibly wrong;
  - `draft`: being worked out;
  - `stable`: tested and relied upon;
  - `deprecated`: superseded, with a pointer to its replacement.
- **Small units.** Prefer small, composable units. Introduce an abstraction only after its pattern has appeared twice.
- **Conservative change.** calx sits at the base of the dependency graph, so changing its interfaces is expensive for every library that depends on it.
  - Prefer additive change.
  - Deprecate before removing.
  - Record every interface change in the decision log.

## 4. Repository Layout

```
calx/
├── compute/       ir · schedule · compiler · runtime · scheduler · backends/cpu
├── numerics/      number systems · linear algebra · random · approximation · special functions ·
│                  transforms · signal · quadrature · integrators · solvers · optimisation · verification
├── docs/          architecture
└── development/   decisions · ideas · release notes
```

- A subdirectory is created when its first content is written. Until then it appears as a `planned` row in its parent's README.
- The component READMEs in `compute/` are an exception. They specify the stack's interfaces.

## 5. Working Procedures

### 5.1 Module Contract

Each module has:

- a `README.md` giving its purpose, status, dependencies, references and open questions;
- public headers under `include/calx/<module>/`;
- its implementation in `src/`;
- unit tests in `tests/`;
- optionally, Python bindings in `bindings/python/` and benchmarks in `benchmarks/`.

### 5.2 Decisions

- Each decision is a numbered record in `development/decisions/`.
- Records are never rewritten. A changed decision gets a new record that supersedes the old one.

### 5.3 Ideas

- Each idea is a file `development/ideas/YYYY-MM-DD-slug.md`.
- Its status is `seed`, `promoted` or `abandoned`.
- Abandoned ideas are kept, with one sentence saying why.

### 5.4 Release Notes

At each commit, add to `development/release-notes/YYYY-MM-DD.md` (one file per day, one short section per commit) a concise summary in plain prose of what is new, what one can now do and what works better, not a list of changes (the git log keeps those); not between commits. This convention is shared with the sibling repositories.

### 5.5 Change Discipline

- Before restructuring, or changing an invariant (§2), propose the change and wait for the user's agreement.
- Keep each change to one concept.

### 5.6 Shared Conventions

calx, morphe, physis and the infrastructure keep these conventions identical:

- the writing style;
- the status vocabulary;
- the decision-log, idea and release-note conventions.

When one of them changes, change it in every repository of [The Science Infrastructures](https://github.com/huangwanhong-maker/The-Science-Infrastructures) together: the foundation libraries and the infrastructure's repository.

### 5.7 Work Inside the Infrastructure

When calx is edited inside the infrastructure's `foundations/calx`:

1. commit in calx first;
2. then commit the updated submodule pointer in the infrastructure.

## 6. Conventions

### 6.1 Numerical Conventions

- Cite the source of every algorithm: a paper or book, with the equation or algorithm number.
- Docstrings state each algorithm's accuracy, stability, complexity and range of validity.
- For a fixed backend, schedule and input, results are reproducible. Results from different backends agree within their declared tolerance.

### 6.2 Writing Style

- Headings are declarative, numbered and academic. Never use question forms.
- Present concepts before structure and usage.
- Label asides with *Rationale.*, *Remark.* or *Caveat.*
- Write in English.

## 7. Technology

- **Language: C++.** This was decided in ADR 0003 of The Science Infrastructures and adopted in decision 0002 here.
- **Provisional choices** (the build is in use since decision 0004):
  - C++20, header-only targets `calx::<module>` while the modules are templates;
  - CMake 3.20 or later; checks are programs on `calx::verification`, registered with CTest;
  - GCC 11 or later (in WSL on the development machine), Clang 14 or later;
  - a C ABI at the runtime boundary;
  - DLPack for exchanging arrays;
  - nanobind for Python bindings.
- **Licences.** calx is MIT-licensed. Check every dependency for compatibility before adopting it. FFTW (GPL) and SHTns (CeCILL), for example, need an explicit decision.
- **Development machine:**
  - Windows 11, Intel i7-11370H, 32 GB;
  - NVIDIA RTX 3050 Ti Laptop GPU (4 GB, CUDA 12.9) and Intel Iris Xe;
  - no C++ compiler or CMake on the Windows side; GCC 11 and CMake 3.28 in WSL.

## 8. Open Decisions

1. **Automatic differentiation.** Its natural home is a transformation on the IR, as Enzyme does on LLVM IR. Decide when to build it. Forward differentiation by operator overloading (jets of order 1 and 2) is written for functions in C++ (decision 0005).
2. **Relation to MLIR.** Is the low-level IR an MLIR dialect, or does it lower to one?
3. **Native CPU path.** Should it emit C++, LLVM IR or something else?
4. **Typed annotations in the IR.** A generic mechanism through which physis can attach units and frames to IR types ([idea of the infrastructure](https://github.com/huangwanhong-maker/The-Science-Infrastructures/blob/main/development/ideas/2026-10-21-units-and-frames-as-types.md)).
5. **Rewriting in the compiler.** Should calx have its own engine, or use morphe's? See the [idea](development/ideas/2026-10-21-equality-saturation-in-the-compiler.md).
6. **Precision policy details.**
7. **Toolchain and build.** This decision is shared with the infrastructure.

## 9. Current State

As of 2026-10-24, the first modules are written and checked (decision 0004): `numerics/random`, `linalg`, `integrators`, `special_functions`, `differentiation` (decision 0005), `transforms`, `signal`, `quadrature`, `solvers`, `verification`, and `compute/scheduler`. For the genome (decision 0006): the CPU backend's target description and vector lanes (`compute/backends/cpu`), task graphs and dynamic loops in the scheduler, the symmetric eigenproblem, LU and determinants, the gamma family of special functions and its variates, and `numerics/optimisation` with Brent's minimiser. For the circulation of value (decision 0007): conjugate gradients on compressed sparse rows. 17 check programs. The rest of the stack is designed and planned. The repository is published at `https://github.com/huangwanhong-maker/calx`.
