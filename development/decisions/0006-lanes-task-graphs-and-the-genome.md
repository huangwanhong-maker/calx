# 0006. Lanes, Task Graphs, and What the Genome Asks of calx

- **Status:** Accepted
- **Date:** 2026-10-24
- **Decided by:** the user asked for a genome analysis infrastructure in The Science Infrastructures, with algorithms "at low level, classical and advanced, high speed and high quality", and with "distributed parallel structures on heterogeneous devices", noting that calx already follows a TVM-like pipeline. Claude placed the general parts here and wrote this record; the composition is The Science Infrastructures' [ADR 0043](https://github.com/huangwanhong-maker/The-Science-Infrastructures/blob/main/development/decisions/0043-the-genome-infrastructure.md).

## 1. Context

The fast algorithms of sequence analysis run on vector lanes (Farrar's striped Smith–Waterman, bit-parallel and wavefront alignment, k-mer hashing), on many threads (counting k-mers, mapping reads, bootstrap replicates, all-against-all distances), and across processes and machines. Their inferences need a symmetric eigensolver (principal components of genotypes, the spectral decomposition of reversible rate matrices), the incomplete gamma and beta functions (χ² tails, discrete-Γ rates, binomial and Student tails), and gamma, beta and binomial variates (read lengths, ancestry, Wright–Fisher sampling). None of this names an object in nature, and the domains may not reimplement it (the infrastructure's rule that domains reuse the foundations), so it belongs here. calx's stack was planned but for the scheduler's thread pool.

## 2. Decision

1. **The CPU backend's first content, `compute/backends/cpu` → `calx::cpu`** (status `draft`):
   - `target.hpp`: the target description the backend contract asks for (architecture §4): the host's capabilities read from CPUID (SSE4.1, AVX2, AVX-512F/BW/VL, POPCNT, BMI2, with the operating system's saving of the vector state checked by XGETBV), its threads, and the build's own: `LANE_BYTES` and `ISA`, the widest vector the compiler was allowed.
   - `lanes.hpp`: `Lanes<T, W>`, W lanes of an arithmetic type on GCC's and Clang's vector extensions, with wrapping and saturating arithmetic, comparisons as masks, select, max, min, shifts by a lane, reductions. One kernel compiles to SSE2, AVX2, AVX-512 or NEON; saturation uses PADDS/PSUBS where the build has them and a generic formula elsewhere, with the same results. This is the native path's vocabulary (cpu README §2); the code generator from the IR remains planned.
2. **Task graphs in `compute/scheduler`** (`task_graph.hpp`): tasks with costs and predecessors, acyclic by construction, run on the `ThreadPool` with the ready task of highest upward rank first (HEFT's priority, Topcuoglu, Hariri and Wu 2002), and a trace of where and when each ran. And `parallel_for.hpp`: dynamic chunks from an atomic counter for irregular loops. Determinism is the caller's: tasks write their own slots, merges read them in order.
3. **`numerics/linalg/eigen.hpp`:** the symmetric eigenproblem by Householder tridiagonalisation and the implicit QL algorithm (tred2 and tql2 of the Handbook, as EISPACK and JAMA have them).
4. **`numerics/special_functions/gamma.hpp`:** log Γ (Lanczos; thread-safe, unlike the C library's lgamma), ψ and ψ′, the regularised incomplete gamma functions and the inverse of P, χ², the regularised incomplete beta function and the binomial, Student and F tails.
5. **`numerics/random`:** gamma variates (Marsaglia and Tsang 2000), beta as a ratio of gammas, and binomial variates by inversion below n·p = 10 and Hörmann's BTRS above.
6. **The first content of `numerics/optimisation`** (`minimise.hpp`): Brent's method for the minimum of a function of one variable on an interval (Brent 1973, ch. 5, localmin), with golden section beside it; the substitution models of phylogenetics fit their parameters with it. And in `linalg/dense.hpp`, the LU factors with partial pivoting, the determinant and its logarithm (LogDet distances).
7. **Distribution across processes and machines stays outside calx.** calx is pure: a task graph runs on this process's threads. The application that needs more (the Genome Atlas's cluster) cuts a job into shards that are pure functions of a few numbers and runs them on engine processes here and elsewhere; calx gives it the deterministic shards and the traces.

## 3. Alternatives Considered

- **Intrinsics in the domains.** Rejected: the domains never address a device directly (the infrastructure's architecture §4), and a kernel per instruction set triples the code.
- **A borrowed SIMD library** (xsimd, Highway, std::experimental::simd). Licence-compatible ones exist; but the vector extensions are already in the compiler, the kernels need a dozen operations, and writing them is how the method is understood (the map of the mathematics, §10.3). The interface is ours, so a borrowed implementation can go behind it later.
- **Run-time dispatch** by function multiversioning. Deferred: the infrastructure builds its engines on the machine that runs them (-march=native), and the target description records what each was built for.
- **A work-stealing scheduler.** Deferred: the shards of the genome's jobs are coarse (thousands of reads, a replicate, a block of a matrix), and a shared ready queue by rank suffices; per-thread deques come when tasks become fine.

## 4. Consequences

- A new target, `calx::cpu`, built with its checks twice (generic and -march=native), and new checks for task graphs, eigenproblems, the gamma family and the new variates.
- The CPU backend's README, the scheduler's, linalg's, special_functions' and random's record what is written; calx's README lists them.
- The genomics and evolution domains of The Science Infrastructures use all of it.
