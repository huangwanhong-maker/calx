// Checks of calx/compute/scheduler: every thread runs every phase once; the parts of a
// range tile it; work done by parts gives the serial answer; an exception in a worker
// reaches the caller and the pool goes on working.

#include <calx/scheduler/thread_pool.hpp>
#include <calx/verification/check.hpp>

#include <atomic>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

using calx::scheduler::part;
using calx::scheduler::ThreadPool;

int main() {
  calx::verification::Suite s;
  bool tiles = true;
  for (std::size_t n : {0, 1, 7, 100, 1001})
    for (std::size_t p : {1, 3, 8}) {
      std::size_t next = 0;
      for (std::size_t t = 0; t < p; t++) {
        const auto r = part(n, p, t);
        tiles &= r.begin == next && r.end >= r.begin && r.size() + 1 >= n / p && r.size() <= n / p + 1;
        next = r.end;
      }
      tiles &= next == n;
    }
  s.require("scheduler: the parts of [0, n) tile it in order, balanced to within one", tiles);

  ThreadPool pool(6);
  std::vector<int> seen(6, 0);
  for (int phase = 0; phase < 1000; phase++) pool.run([&](std::size_t t) { seen[t]++; });
  bool all = true;
  for (int v : seen) all &= v == 1000;
  s.require("scheduler: six threads each run every one of 1000 phases once", all);

  const std::size_t n = 1000003;
  std::vector<double> partial(pool.size(), 0);
  pool.run([&](std::size_t t) {
    const auto r = part(n, pool.size(), t);
    for (std::size_t i = r.begin; i < r.end; i++) partial[t] += 1.0 / (1.0 + static_cast<double>(i));
  });
  double sum = 0, serial = 0;
  for (double v : partial) sum += v;
  for (std::size_t i = 0; i < n; i++) serial += 1.0 / (1.0 + static_cast<double>(i));
  s.check("scheduler: a sum done by parts equals the serial sum (to rounding)", sum - serial, 1e-12);

  bool caught = false;
  try {
    pool.run([](std::size_t t) { if (t == 3) throw std::runtime_error("worker 3"); });
  } catch (const std::runtime_error& e) { caught = std::string(e.what()) == "worker 3"; }
  std::atomic<int> after{0};
  pool.run([&](std::size_t) { after++; });
  s.require("scheduler: a worker's exception reaches the caller, and the pool works on", caught && after == 6);
  return s.done();
}
