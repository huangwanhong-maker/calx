// A pool of persistent threads that run one task on every thread at once and wait for all
// of them: the bulk-synchronous step of a simulation (update the neurons of each thread,
// then deliver the spikes to the targets of each thread). The caller is thread 0.
//
// Workers spin briefly on the pool's generation counter before sleeping, so that the many
// short phases of a simulation do not each pay for a wake-up through the kernel. They spin
// on the processor's pause instruction, not by yielding: a yield is a system call, and many
// threads yielding cost more than the phases they wait for. An exception in any thread is
// carried to the caller.
//
// The first content of the task scheduler: all-CPU placement of one bulk-synchronous phase
// at a time. Task graphs and devices come with the rest of the stack.

#pragma once

#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <cstddef>
#include <exception>
#include <functional>
#include <mutex>
#include <thread>
#include <vector>

namespace calx::scheduler {

// A hint to the processor that this is a spin-wait loop.
inline void relax() {
#if defined(__x86_64__) || defined(__i386__)
  __builtin_ia32_pause();
#elif defined(__aarch64__)
  asm volatile("yield");
#endif
}

// The t-th of T contiguous parts of [0, n), the first n % T one longer.
struct Range {
  std::size_t begin = 0, end = 0;
  std::size_t size() const { return end - begin; }
};
inline Range part(std::size_t n, std::size_t parts, std::size_t t) {
  const std::size_t q = n / parts, r = n % parts;
  const std::size_t b = t * q + std::min(t, r);
  return {b, b + q + (t < r ? 1 : 0)};
}

class ThreadPool {
 public:
  explicit ThreadPool(std::size_t threads) : threads_(std::max<std::size_t>(1, threads)) {
    for (std::size_t t = 1; t < threads_; t++) workers_.emplace_back([this, t] { work(t); });
  }
  ~ThreadPool() {
    {
      std::lock_guard<std::mutex> lock(mutex_);
      stop_ = true;
      generation_.fetch_add(1, std::memory_order_release);
    }
    wake_.notify_all();
    for (auto& w : workers_) w.join();
  }
  ThreadPool(const ThreadPool&) = delete;
  ThreadPool& operator=(const ThreadPool&) = delete;

  std::size_t size() const { return threads_; }

  // Run task(t) on every thread t and return when all have finished.
  void run(const std::function<void(std::size_t)>& task) {
    if (threads_ == 1) { task(0); return; }
    {
      std::lock_guard<std::mutex> lock(mutex_);
      task_ = &task;
      error_ = nullptr;
      pending_.store(threads_ - 1, std::memory_order_relaxed);
      generation_.fetch_add(1, std::memory_order_release);
    }
    wake_.notify_all();
    try { task(0); } catch (...) { record(std::current_exception()); }
    // Wait for the workers, spinning first.
    for (int i = 0; pending_.load(std::memory_order_acquire) != 0; i++) {
      relax();
      if (i > (1 << 14)) {
        std::unique_lock<std::mutex> lock(mutex_);
        done_.wait(lock, [&] { return pending_.load(std::memory_order_acquire) == 0; });
        break;
      }
    }
    task_ = nullptr;
    if (error_) std::rethrow_exception(error_);
  }

 private:
  void record(std::exception_ptr e) {
    std::lock_guard<std::mutex> lock(error_mutex_);
    if (!error_) error_ = e;
  }
  void work(std::size_t t) {
    std::size_t seen = 0;
    for (;;) {
      std::size_t g = generation_.load(std::memory_order_acquire);
      for (int i = 0; g == seen && i < (1 << 14); i++) {
        relax();
        g = generation_.load(std::memory_order_acquire);
      }
      if (g == seen) {
        std::unique_lock<std::mutex> lock(mutex_);
        wake_.wait(lock, [&] { return generation_.load(std::memory_order_acquire) != seen; });
        g = generation_.load(std::memory_order_acquire);
      }
      seen = g;
      if (stop_) return;
      try { (*task_)(t); } catch (...) { record(std::current_exception()); }
      if (pending_.fetch_sub(1, std::memory_order_acq_rel) == 1) {
        std::lock_guard<std::mutex> lock(mutex_);
        done_.notify_all();
      }
    }
  }

  std::size_t threads_;
  std::vector<std::thread> workers_;
  std::mutex mutex_, error_mutex_;
  std::condition_variable wake_, done_;
  std::atomic<std::size_t> generation_{0}, pending_{0};
  const std::function<void(std::size_t)>* task_ = nullptr;
  std::exception_ptr error_;
  bool stop_ = false;
};

}  // namespace calx::scheduler
