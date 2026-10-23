# scheduler: Task Graphs Across Devices

*Status: draft (2026-10-23): a pool of host threads for bulk-synchronous phases.* Task graphs, placement and data movement are planned.

## 1. Purpose

The scheduler decides when and where computations run. A task graph is a set of kernels connected by data dependencies. The scheduler places the graph on one or more devices and moves data between them as needed.

*Remark.* The scheduler is distinct from [`schedule/`](../schedule/), which decides how one computation is mapped onto one device.

## 2. Content

`include/calx/scheduler/thread_pool.hpp`: `ThreadPool`, persistent threads that run one task on every thread at once and wait for all, the caller being thread 0; workers spin briefly before sleeping, so that the many short phases of a simulation do not each pay for a wake-up; an exception in any thread reaches the caller. `part(n, parts, t)` divides a range into balanced contiguous parts. This is the simplest placement, all on the CPU, of one bulk-synchronous phase at a time; its first user is the neuroscience domain's spiking networks, which update their neurons and deliver their spikes in alternating phases.

Checked by `tests/check_scheduler.cpp`: every thread runs every phase once; the parts tile the range; work by parts gives the serial answer; an exception reaches the caller and the pool works on.

## 3. Planned Content

- A representation of task graphs.
- Placement policies, starting with all-CPU placement.
- Planning of data movement.
- Later:
  - heterogeneous placement driven by cost models;
  - distributed execution across nodes.
