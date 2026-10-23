---
title: Spiking networks and quantum trajectories as targets of the compiler
status: seed
date: 2026-10-23
---

# Spiking Networks and Quantum Trajectories as Targets of the Compiler

*Recorded on 2026-10-23 at the user's direction. The first programs that calx's compute stack must compile are the spiking networks of the infrastructure's neuroscience domain, which are simulated today by C++ on the CPU's threads. Compiling them to devices is planned, not built. The domain's side is in its idea [Spiking networks compiled to devices](https://github.com/huangwanhong-maker/neuroscience/blob/main/development/ideas/2026-10-23-spiking-networks-on-devices.md).*

## 1. What the IR Must Express

- **Maps over many small systems:** the same small update for each element of a population, which is a propagator applied, an integrator's step, or a closed-form exponential of a 2 × 2 complex matrix.
- **Events:**
  - a test after each update, which emits records (time, index) into a buffer;
  - jumps conditioned on the test;
  - silent jumps applied unconditionally.
- **Sparse scatter at delays:** each event read once, its row of targets walked, and weights added to ring buffers indexed by arrival time.
- **Event-driven state on edges:** lazy updates of synapses over a history, which are plasticity's traces and closed-form integrals.
- **Counter-based randomness:** Philox streams as a first-class operation, keyed by purpose, entity and step, so that every device draws the same numbers.
- **Declared precision,** fixed point among it, with quantisation explicit.

## 2. Backends Beyond the Usual

- **Neuromorphic chips** (Loihi 2, SpiNNaker 2, BrainScaleS-2). Their native operations are neurons, synapses and spikes, so the backend's task is to map a network's declaration onto them, not to emit instructions. A device-faithful model of the chip's neuron, in the domain's library, is the reference the CPU interpreter checks it against.
- **Quantum processors.** A quantum trajectory becomes a circuit with mid-circuit measurement and classical feedback:
  - pulses become rotations;
  - detection becomes measurement;
  - a reset follows on the outcome.

  The reference is the master equation, which physis's `quantum` module integrates exactly.

## 3. Consequence for the Design of the IR

Events, sparse scatter and stateful edges are not tensor operations. They belong in the IR from the start, as dialects beside the dense ones, rather than being added to a tensor IR later. This bears on the [equality-saturation idea](2026-10-21-equality-saturation-in-the-compiler.md): its rewrites must respect the order of events.
