# backends: Target Plug-ins

*Status: planned.*

## 1. Contract

Each backend is one directory that provides three things:

1. **A target description:** capabilities, precision support and cost hints.
2. **A code generator:** from low-level IR to executable form.
3. **A device-API implementation:** an implementation of [`runtime/`](../runtime/).

A backend registers as a plug-in and is tested against the CPU reference interpreter.

## 2. Planned Backends

| Backend | Target | Status |
|---|---|---|
| [`cpu/`](cpu/) | Reference interpreter and native CPU code | draft: the target description and vector lanes (decision 0006); the interpreter and code generator planned |
| `cuda` | NVIDIA GPUs (development machine: RTX 3050 Ti) | planned |
| `sycl` | Intel GPUs via Level Zero (development machine: Iris Xe), and other SYCL devices | planned |
| `opencl` | Portable accelerators | planned |
| `vulkan` | SPIR-V compute | planned |
| `metal` | Apple GPUs | planned |
| `webgpu` | Browsers (WGSL) | planned |
| `fpga` | FPGAs, through high-level synthesis | planned |

Devices not listed here are added the same way, including devices that do not exist yet.
