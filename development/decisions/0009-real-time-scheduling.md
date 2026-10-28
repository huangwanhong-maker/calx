# 0009. Real-Time Scheduling: Tasks with Deadlines, Their Analyses, Servers, Simulation and Dispatch

- **Status:** Accepted, provisionally (dispatch is non-preemptive and measures lateness; a backend on a real-time kernel is planned)
- **Date:** 2026-10-28
- **Decided by:** the user asked Cosmopoiesis, a work of The Science Infrastructures, for instruments whose readings are acquired at the ticks of a clock and on events, from physical and virtual sensors, "prioritised and parallel, batched, with events that do not wait for the tick, without conflict and safely, and without the starvation a priority queue that lets jobs jump ahead can cause", noting that these are well-studied questions of computer science. Claude placed the general methods here, since they name no instrument and no science, and wrote this record; the composition is the infrastructure's [ADR 0050](https://github.com/huangwanhong-maker/The-Science-Infrastructures/blob/main/development/decisions/0050-instruments-and-their-acquisition.md).

## 1. Context

calx's scheduler had a pool of threads for bulk-synchronous phases, task graphs ordered by HEFT's rank, and dynamic loops (decision 0006). It had no notion of time: of a job released at an instant and due by another, of a stream of work that must be given a share and no more, or of a test that says beforehand whether a set of periodic duties can be met. An acquisition runtime needs all three, and they are textbook: Liu and Layland (1973) for periodic tasks and their bounds, the response-time analysis of Joseph and Pandya (1986), the processor-demand criterion of Baruah, Rosier and Howell (1990), Abeni and Buttazzo's constant bandwidth server (1998), Buttazzo's elastic model (2002).

## 2. Decision

`compute/scheduler/include/calx/scheduler/realtime.hpp` adds, without changing the existing headers:

- **The task model.** `PeriodicTask` (C, T, D, φ, priority, criticality, server); `Job` (release, cost, deadline); `periodic_jobs` releasing a set's jobs within a horizon.
- **Analyses.** `utilisation`, `density`, `liu_layland_bound`, `rate_monotonic` (and deadline-monotonic) orders, `response_times` with blocking, `busy_period`, `demand_bound`, `edf_feasible`.
- **Servers.** `ConstantBandwidthServer` with the four rules of Abeni and Buttazzo.
- **Degradation.** `compress`, elastic compression of periods to a desired utilisation.
- **Simulation.** `simulate`, deterministic, on m processors, preemptive or not, by EDF, fixed priority (with aging) or FIFO, with servers; ties go to the earlier release, the lower task index, the earlier job.
- **Dispatch.** `dispatch`, released jobs on a `ThreadPool` in the policy's order, run to their end, jobs of one group never at once, with their starts, ends and lateness.

Times are seconds of a declared clock, logical or of the wall; the module does not care which. The policy is a choice of the caller, not of the module.

*Rationale.* The analyses and the simulation are two computations of one schedule, so the checks hold each against the other on random task sets (the response-time analysis against the simulated critical instant; the demand criterion against the simulated synchronous schedule), as calx's convention prefers an independent computation to a remembered table. Starvation is answered by two mechanisms the user may choose between: servers, which give each stream a guaranteed share and isolate the rest from it (analysable), and aging, which bounds a ready job's wait under fixed priorities (simple, but it breaks the exact analysis).

## 3. Alternatives Considered

- **Leave scheduling to the kernel.** Linux's `SCHED_DEADLINE` (EDF with constant bandwidth servers) and `PREEMPT_RT` are the right foundation for guarantees, and a backend on them is planned. But they schedule threads, not the jobs inside a thread, they need privileges and Linux, and the development machine runs Windows; the analyses and the simulation are needed whatever the kernel.
- **A utility function deciding everything.** A weighted sum of a request's importance, its data's age, its expected information and its cost mixes quantities of different kinds with arbitrary weights and gives no guarantee to the duties that must be met. Criticality and deadlines are kept for what must be met (with servers and admission), and utility is left to the caller's choice of what to release (a policy above this module).
- **A scheduler of our own in each work.** It would be the third (after the neuroscience engine's phases and the genome workbench's jobs) and would drift; this module is general.

## 4. Consequences

- `tests/check_realtime.cpp` (34 checks) holds the analyses to the simulation and the rules to cases worked by hand.
- Cosmopoiesis's acquisition runtime schedules its instruments' jobs by EDF with constant bandwidth servers, checks a set of instruments for admission by `edf_feasible`, degrades under overload by `compress`, and runs its virtual instruments' jobs by `dispatch`.
- Planned: preemption points within long jobs; a backend that maps servers to `SCHED_DEADLINE` reservations where the kernel offers them; mixed-criticality analysis (Vestal 2007) if a work needs more than degradation by compression.
