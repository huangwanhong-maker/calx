# cpu: CPU Backend

*Status: draft (2026-10-24): the target description and vector lanes ([decision 0006](../../../development/decisions/0006-lanes-task-graphs-and-the-genome.md)). The reference interpreter and the code generator are planned.* This is the first backend.

## 1. Purpose

The CPU is the primary execution target. This backend also defines correct results for every other backend.

## 2. Content

| Header | Content |
|---|---|
| `include/calx/cpu/target.hpp` | `Target`, `host()`, `describe()`: the host's architecture, vendor and name, threads, and capabilities read from CPUID (SSE2, SSE4.1/4.2, POPCNT, AVX, AVX2, BMI2, AVX-512F/BW/VL), AVX and AVX-512 counted only when the operating system saves their registers (XGETBV); the build's `LANE_BYTES` and `ISA`, the widest vector the compiler was allowed |
| `include/calx/cpu/lanes.hpp` | `Lanes<T, W>`: W lanes of an arithmetic type on GCC's and Clang's vector extensions: splat, load and store (unaligned), wrapping +, −, ×, bitwise operations and shifts, comparisons giving a `Mask`, select, max, min, saturating `adds` and `subs` (the PADDS/PSUBS family where the build has it, a generic formula with the same results elsewhere), `shift_up` and `shift_down` by one lane, `hmax`, `hmin`, `hsum`; `Mask` with `any`, `all`, `bits` |

A kernel written on `Lanes<T>` takes the build's widest vector: 64 bytes under -march=native on an AVX-512 machine, 32 with AVX2, 16 with SSE2 or NEON, and correct (emulated) code where there is none. There is no dispatch at run time; a program that wants the host's widest lanes is built for it, and `host()` against `ISA` says whether it was.

## 3. Planned Content

- **Reference interpreter.** It executes the IR directly and simply, with no optimisation. Clarity takes priority over speed.
- **Native path.** It compiles the IR for performance, using multithreading and the lanes. Whether it emits C++, LLVM IR or something else is decided when it is built.

## 4. Checks

`check_calx_cpu` and `check_calx_cpu_native` (the same program, built with the default flags and with -march=native): the host is described and runs what the build emits; every operation of `Lanes<T, W>` for 8-, 16-, 32- and 64-bit integers, signed and unsigned, and both floating-point types, at 16, 32 and 64 bytes, agrees lane by lane with scalar code on random values and on the limits of each type (wrapping, saturation, shifts across lanes, reductions, masks); a step of Farrar's striped kernel.

## 5. References

- Intel 64 and IA-32 Architectures Software Developer's Manual, vol. 2A, CPUID and XGETBV; AMD64 Architecture Programmer's Manual, vol. 3, appendix E.
- GCC manual, *Using Vector Instructions through Built-in Functions* (vector extensions, `__builtin_shuffle`).
- Farrar, M. (2007), Striped Smith–Waterman speeds database searches six times over other SIMD implementations, *Bioinformatics* 23, 156–161 (the shift across lanes).
- Warren, H. S. (2013), *Hacker's Delight*, 2nd ed., §2-13 (overflow of signed addition).
