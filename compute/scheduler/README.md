# scheduler: Task Graphs Across Devices

*Status: draft (2026-10-24): a pool of host threads for bulk-synchronous phases; task graphs on the host's threads by rank, with their traces; dynamic loops ([decision 0006](../../development/decisions/0006-lanes-task-graphs-and-the-genome.md)); real-time scheduling, its analyses, servers, simulation and dispatch (2026-10-28, [decision 0009](../../development/decisions/0009-real-time-scheduling.md)).* Placement on other devices and data movement are planned.

## 1. Purpose

The scheduler decides when and where computations run. A task graph is a set of kernels connected by data dependencies. The scheduler places the graph on one or more devices and moves data between them as needed.

*Remark.* The scheduler is distinct from [`schedule/`](../schedule/), which decides how one computation is mapped onto one device.

## 2. Content

`include/calx/scheduler/thread_pool.hpp`: `ThreadPool`, persistent threads that run one task on every thread at once and wait for all, the caller being thread 0; workers spin briefly before sleeping, so that the many short phases of a simulation do not each pay for a wake-up; an exception in any thread reaches the caller. `part(n, parts, t)` divides a range into balanced contiguous parts. This is the simplest placement, all on the CPU, of one bulk-synchronous phase at a time; its first user is the neuroscience domain's spiking networks, which update their neurons and deliver their spikes in alternating phases.

Checked by `tests/check_scheduler.cpp`: every thread runs every phase once; the parts tile the range; work by parts gives the serial answer; an exception reaches the caller and the pool works on.

`include/calx/scheduler/task_graph.hpp`: `TaskGraph`, tasks with a declared cost and the tasks they wait for (only tasks added before them, so the graph is acyclic by construction and the order of addition is topological); `upward_rank()`, the declared cost of the longest path from each task to the end, and `critical_path()`; `run(pool)` starts every task as soon as its predecessors have ended, a free thread taking the ready task of highest rank (HEFT's priority, Topcuoglu, Hariri and Wu 2002, with one kind of processor), and returns a `Trace`: each task's thread, start and end, the makespan, the time in tasks and the utilisation. An exception in a task stops new tasks, lets those running end, and reaches the caller. Determinism is the caller's: tasks write their own outputs and merges read them in order.

`include/calx/scheduler/parallel_for.hpp`: `parallel_for(pool, n, grain, body)`, [0, n) in chunks of `grain` taken from an atomic counter (OpenMP's dynamic schedule), for iterations of unpredictable cost; the body learns its thread for per-thread scratch space.

Checked by `tests/check_task_graph.cpp`: the ranks and critical path of a diamond; a task cannot wait for one not yet added; on 1, 2, 4 and 8 threads every task of a random graph runs once and after its predecessors, with an ordered trace and utilisation in (0, 1]; the task on the critical path goes first; 64 shards and a merge give the serial sum on every thread count; an exception reaches the caller and its successors do not run; `parallel_for` covers its range once in chunks of the grain.

`include/calx/scheduler/realtime.hpp`: real-time scheduling of periodic and sporadic tasks (Liu and Layland 1973; Mok 1983), for work that must be done by deadlines: the readings of instruments, the steps of a controller, the frames of a display. A `PeriodicTask` is (C, T, D, φ) with a priority, a criticality and an optional server.

- **Analyses.** `utilisation` and `density`; `liu_layland_bound(n)`, the sufficient bound of rate-monotonic priorities; `response_times`, the exact response-time analysis of preemptive fixed priorities (Joseph and Pandya 1986; Audsley et al. 1993) with an optional blocking term; `busy_period`, `demand_bound` and `edf_feasible`, EDF's utilisation test for implicit deadlines and the processor-demand criterion otherwise (Baruah, Rosier and Howell 1990).
- **Servers.** `ConstantBandwidthServer`, the CBS of Abeni and Buttazzo (1998): a stream of jobs of unknown cost served with a budget Q every T takes no more than Q/T from the others, whatever its jobs do, which is how a source of sporadic events is kept from starving the periodic work and how a guaranteed share is given.
- **Degradation.** `compress`, the elastic compression of periods to a desired utilisation (Buttazzo, Lipari, Caccamo and Abeni 2002, §3).
- **Simulation.** `simulate`, a deterministic discrete-event simulation of jobs on m identical processors, preemptive or not, by EDF, fixed priority (with or without aging, which bounds a ready job's wait) or FIFO, with servers; it returns every slice of execution and every job's start, finish and lateness.
- **Dispatch.** `dispatch`, released jobs run on a `ThreadPool`, a free thread taking the ready job first in the policy's order and running it to its end, jobs of one group never at once, with each job's start, end and lateness on the host's clock.

*Caveat.* The analyses are of the model and hold for the costs declared. User code on a general-purpose operating system cannot be preempted by the library and can be delayed by the kernel; `dispatch` measures lateness, it does not guarantee it. A guarantee needs a real-time kernel beneath (Linux's `SCHED_DEADLINE`, which is EDF with constant bandwidth servers in the kernel; or `PREEMPT_RT` with fixed priorities).

`include/calx/scheduler/kernel.hpp` asks the kernel what it grants ([decision 0010](../../development/decisions/0010-kernel-awareness-and-extended-range-legendre.md)):

- **`probe_kernel()`** gives a `KernelReport`: the system and its release, WSL, PREEMPT_RT, whether SCHED_FIFO and SCHED_DEADLINE would be granted, SCHED_FIFO's priorities, RLIMIT_RTPRIO, and why a class is refused. Each class is tried on a short-lived thread and given back, so the caller is unchanged.
- **`raise_to_fifo(priority)`**, **`lower_to_time_sharing()`** and **`current_class()`** act on the calling thread; a pool raises its workers through its own `run`.

Off Linux it says that calx offers no real-time class there.

Checked by `tests/check_kernel.cpp`: the probe leaves the caller unchanged; a raise succeeds exactly when the probe says it would, and lowering works; every refusal is explained. On WSL 2 both classes are refused (no CAP_SYS_NICE, RLIMIT_RTPRIO 0).

Checked by `tests/check_realtime.cpp`: the Liu–Layland bound and its limit; response-time analysis against the simulated schedule of 400 random sets (the first job's response equals Rᵢ, no job's exceeds it, the first rejected task misses); EDF's tests against the simulation (U ≤ 1 meets every implicit deadline and U > 1 misses; the demand criterion agrees on 400 constrained-deadline sets; EDF schedules every set rate-monotonic priorities do); the CBS's rules worked by hand and its isolation beside a stream demanding eight times its share; aging ending a starvation at the step it should; non-preemptive runs in one slice; global EDF on two processors work-conserving; elastic compression; dispatch's order, its once-only runs, its groups and its lateness.

## 3. Planned Content

- Placement policies, starting with all-CPU placement.
- Planning of data movement.
- Later:
  - heterogeneous placement driven by cost models;
  - distributed execution across nodes.
