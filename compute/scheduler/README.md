# scheduler: Task Graphs Across Devices

*Status: draft (2026-10-24): a pool of host threads for bulk-synchronous phases; task graphs on the host's threads by rank, with their traces; dynamic loops ([decision 0006](../../development/decisions/0006-lanes-task-graphs-and-the-genome.md)).* Placement on other devices and data movement are planned.

## 1. Purpose

The scheduler decides when and where computations run. A task graph is a set of kernels connected by data dependencies. The scheduler places the graph on one or more devices and moves data between them as needed.

*Remark.* The scheduler is distinct from [`schedule/`](../schedule/), which decides how one computation is mapped onto one device.

## 2. Content

`include/calx/scheduler/thread_pool.hpp`: `ThreadPool`, persistent threads that run one task on every thread at once and wait for all, the caller being thread 0; workers spin briefly before sleeping, so that the many short phases of a simulation do not each pay for a wake-up; an exception in any thread reaches the caller. `part(n, parts, t)` divides a range into balanced contiguous parts. This is the simplest placement, all on the CPU, of one bulk-synchronous phase at a time; its first user is the neuroscience domain's spiking networks, which update their neurons and deliver their spikes in alternating phases.

Checked by `tests/check_scheduler.cpp`: every thread runs every phase once; the parts tile the range; work by parts gives the serial answer; an exception reaches the caller and the pool works on.

`include/calx/scheduler/task_graph.hpp`: `TaskGraph`, tasks with a declared cost and the tasks they wait for (only tasks added before them, so the graph is acyclic by construction and the order of addition is topological); `upward_rank()`, the declared cost of the longest path from each task to the end, and `critical_path()`; `run(pool)` starts every task as soon as its predecessors have ended, a free thread taking the ready task of highest rank (HEFT's priority, Topcuoglu, Hariri and Wu 2002, with one kind of processor), and returns a `Trace`: each task's thread, start and end, the makespan, the time in tasks and the utilisation. An exception in a task stops new tasks, lets those running end, and reaches the caller. Determinism is the caller's: tasks write their own outputs and merges read them in order.

`include/calx/scheduler/parallel_for.hpp`: `parallel_for(pool, n, grain, body)`, [0, n) in chunks of `grain` taken from an atomic counter (OpenMP's dynamic schedule), for iterations of unpredictable cost; the body learns its thread for per-thread scratch space.

Checked by `tests/check_task_graph.cpp`: the ranks and critical path of a diamond; a task cannot wait for one not yet added; on 1, 2, 4 and 8 threads every task of a random graph runs once and after its predecessors, with an ordered trace and utilisation in (0, 1]; the task on the critical path goes first; 64 shards and a merge give the serial sum on every thread count; an exception reaches the caller and its successors do not run; `parallel_for` covers its range once in chunks of the grain.

## 3. Planned Content

- Placement policies, starting with all-CPU placement.
- Planning of data movement.
- Later:
  - heterogeneous placement driven by cost models;
  - distributed execution across nodes.
