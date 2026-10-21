# scheduler: Task Graphs Across Devices

*Status: planned.*

## 1. Purpose

The scheduler decides when and where computations run. A task graph is a set of kernels connected by data dependencies. The scheduler places the graph on one or more devices and moves data between them as needed.

*Remark.* The scheduler is distinct from [`schedule/`](../schedule/), which decides how one computation is mapped onto one device.

## 2. Planned Content

- A representation of task graphs.
- Placement policies, starting with all-CPU placement.
- Planning of data movement.
- Later:
  - heterogeneous placement driven by cost models;
  - distributed execution across nodes.
