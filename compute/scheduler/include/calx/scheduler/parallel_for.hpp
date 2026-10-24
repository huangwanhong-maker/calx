// A loop of irregular iterations on the pool: [0, n) cut into chunks of `grain` iterations,
// which the threads take in turn from an atomic counter until none are left (dynamic
// scheduling, as OpenMP's schedule(dynamic, grain)). It suits work whose cost varies from one
// iteration to the next and cannot be predicted, such as mapping reads of different lengths
// or aligning pairs of different similarity; for uniform work the static parts of part() in
// thread_pool.hpp cost less.
//
//   parallel_for(pool, reads.size(), 64, [&](std::size_t begin, std::size_t end, std::size_t thread) {
//     for (std::size_t i = begin; i < end; i++) out[i] = map(reads[i], scratch[thread]);
//   });
//
// Which thread runs which chunk depends on timing, so a deterministic result writes each
// iteration's answer to its own slot (out[i]) and never accumulates across chunks in the
// order they finish. The third argument is the thread's index, for per-thread scratch space.

#pragma once

#include <calx/scheduler/thread_pool.hpp>

#include <algorithm>
#include <atomic>
#include <cstddef>

namespace calx::scheduler {

template <class Body>
void parallel_for(ThreadPool& pool, std::size_t n, std::size_t grain, Body&& body) {
  if (n == 0) return;
  grain = std::max<std::size_t>(1, grain);
  const std::size_t chunks = (n + grain - 1) / grain;
  if (pool.size() == 1 || chunks == 1) {
    for (std::size_t c = 0; c < chunks; c++) body(c * grain, std::min(n, (c + 1) * grain), std::size_t{0});
    return;
  }
  std::atomic<std::size_t> next{0};
  pool.run([&](std::size_t thread) {
    for (;;) {
      const std::size_t c = next.fetch_add(1, std::memory_order_relaxed);
      if (c >= chunks) return;
      body(c * grain, std::min(n, (c + 1) * grain), thread);
    }
  });
}

}  // namespace calx::scheduler
