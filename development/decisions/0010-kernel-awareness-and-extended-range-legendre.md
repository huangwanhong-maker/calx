# 0010. The Kernel's Real-Time Classes, Asked; Legendre Functions in Extended Range

- **Status:** Accepted, provisionally (Linux's classes only; SCHED_DEADLINE is reported, not yet used for dispatch)
- **Date:** 2026-10-28
- **Decided by:** the user asked of Cosmopoiesis's acquisition "whether currently the scheduler algorithm is RT kernel aware", and asked for the Moon's magnetic field, published as an expansion of degree 449, to be read by its virtual instruments. Claude placed the general methods here and wrote this record. The composition is the infrastructure's [ADR 0052](https://github.com/huangwanhong-maker/The-Science-Infrastructures/blob/main/development/decisions/0052-fields-known-in-part.md); the potential fields are physis's [decision 0010](https://github.com/huangwanhong-maker/physis/blob/main/development/decisions/0010-potential-fields.md).

## 1. Context

**The scheduler did not know the kernel.** Decision 0009's dispatcher orders jobs by EDF, fixed priority with aging or FIFO on a pool of threads, and measures their lateness. Whether another process takes a worker's processor is the operating system's decision, and nothing asked it: the workers ran in the time-sharing class, and the caveat of 0009 ("a guarantee needs a real-time kernel beneath") stayed a remark.

**The Legendre functions failed at high degree near the poles.** The Schmidt functions of 0009, by recurrences on the normalised functions, start each order from the sectoral value $\breve P_m^m\propto\cos^m\varphi$. Near a pole at degree 449 this underflows below the smallest double: at $\varphi=89.99°$, $\cos^{449}\varphi\approx10^{-1500}$. The order's column is then lost, though its values at higher degree return to normal size.

## 2. Decision

### 2.1 The Kernel Asked

`compute/scheduler/include/calx/scheduler/kernel.hpp`:

- **`probe_kernel()`** returns a `KernelReport`: the system, its release and version; whether it is WSL; whether it is PREEMPT_RT (from `/sys/kernel/realtime` and the version string); whether SCHED_FIFO and SCHED_DEADLINE would be granted; SCHED_FIFO's range of priorities; RLIMIT_RTPRIO; and the reason for any refusal, from `errno` and the process's limits.
  - Each class is tried on a short-lived thread of its own, raised and given back, so the probe leaves the caller as it was.
  - SCHED_DEADLINE is tried through the `sched_setattr` system call, policy 6, with a runtime of 1 ms in 10 ms.
- **`raise_to_fifo(priority)`**, **`lower_to_time_sharing()`** and **`current_class()`** act on the calling thread.

A dispatcher's pool raises its workers on request, through the pool's own `run`. Off Linux the probe says that calx offers no real-time class there, and that lateness is measured, not bounded.

On the development machine, WSL 2 refuses both classes: the process lacks CAP_SYS_NICE and its RLIMIT_RTPRIO is 0. Its Linux kernel is itself a guest scheduled by Windows.

### 2.2 Legendre Functions in Extended Range

`schmidt_legendre` keeps each order's column in extended range (Fukushima 2012): a value as a double and a binary exponent of its own.

- The sectoral seeds are carried normalised by `frexp`.
- Each column is renormalised when a value exceeds $2^{480}$.
- Values return to doubles with `ldexp`, clamped to zero where they are below the smallest double.

The functions now hold at every latitude, the poles included, at degrees far beyond any crustal model's. The addition theorem $\sum_m(\breve P_n^m)^2=1$ holds to $n^2\varepsilon$, which error analysis of the three-term recurrence predicts: each of $n$ steps loses $O(n\varepsilon)$ in relative terms at worst. It is checked at degrees 449, 1000 and 2700, at latitudes from the equator to $2\times10^{-6}$ rad of each pole. The derivatives are checked by the identity in the degree and by central differences at degree 449 near the poles.

## 3. Alternatives Considered

- **SCHED_FIFO for every pool by default.** Rejected. It needs privilege, and a real-time thread that spins starves the system. It is asked for explicitly, and reported.
- **SCHED_DEADLINE as the dispatcher.** Deferred. It is EDF with constant bandwidth servers in the kernel, the natural backend for decision 0009's model, but it needs each task's runtime, deadline and period declared to the kernel and a kernel that grants it. The probe says whether one does.
- **Legendre functions by Clenshaw's summation with scaling** (Holmes and Featherstone 2002). It is equivalent for sums at a point, but does not give the individual functions, which a grid's synthesis by orders needs. Extended range keeps the functions and costs little.
- **Long double.** Its range ends near $10^{\pm4932}$, beyond the underflow at degree 449 but not at degree 2700 near the poles, and it is not portable. It is not used.

## 4. Consequences

- A check program is added, `check_calx_kernel`. It holds that:
  - the probe leaves the caller unchanged;
  - a raise succeeds exactly when the probe says it would;
  - lowering works;
  - every refusal has a reason.
- `check_calx_legendre` gains the sum rule at high degree and the derivatives near the poles.
- Cosmopoiesis's sessions report the kernel and may raise their workers (its decision 0005).
- physis `analysis` evaluates the Moon's field of degree 449 at every latitude on these functions.
