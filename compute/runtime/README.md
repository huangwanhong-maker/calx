# runtime: Device-Neutral Execution API

*Status: planned.*

## 1. Purpose

The runtime is the minimal interface through which compiled modules run on a device. Every backend implements it, so it is kept small and stable.

## 2. Planned Content

- **Devices:** enumeration and capabilities.
- **Memory:** buffers, allocation, and transfers between host and device.
- **Queues and events:** ordering and synchronisation.
- **Modules:** loading compiled kernels and launching them.
- **A C ABI at the boundary**, for stability across compilers and languages.
- **Arrays exchanged as DLPack.** NumPy, PyTorch and JAX also understand this format.
