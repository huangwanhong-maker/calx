// Checks of task graphs and dynamic loops: every task runs once, after all its predecessors;
// the upward ranks and the critical path of a known graph; the ready task of highest rank goes
// first; a reduction by shards and a merge gives the serial answer on 1, 2, 4 and 8 threads;
// the trace's times are ordered and its utilisation lies in (0, 1]; an exception in a task
// reaches the caller and the pool works on; parallel_for covers its range exactly once in
// chunks of the grain, on every thread count.

#include <calx/scheduler/parallel_for.hpp>
#include <calx/scheduler/task_graph.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <mutex>
#include <numeric>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

using namespace calx::scheduler;

int main() {
  calx::verification::Suite s;

  // A diamond with a tail: a -> {b, c} -> d -> e, costs 1, 5, 2, 1, 3.
  {
    TaskGraph g;
    const auto a = g.add("a", nullptr, {}, 1);
    const auto b = g.add("b", nullptr, {a}, 5);
    const auto c = g.add("c", nullptr, {a}, 2);
    const auto d = g.add("d", nullptr, {b, c}, 1);
    g.add("e", nullptr, {d}, 3);
    const auto rank = g.upward_rank();
    s.require("upward ranks of a diamond: 10, 9, 6, 4, 3", rank == std::vector<double>{10, 9, 6, 4, 3});
    s.check("its critical path is a, b, d, e: 10", g.critical_path() - 10, 0);
    bool threw = false;
    try { g.add("loop", nullptr, {7}); } catch (const std::invalid_argument&) { threw = true; }
    s.require("a task cannot wait for one not yet added, so the graph stays acyclic", threw);
  }

  // Order: every task after its predecessors, each exactly once, on several thread counts.
  for (std::size_t threads : {1, 2, 4, 8}) {
    ThreadPool pool(threads);
    TaskGraph g;
    const std::size_t n = 400;
    std::vector<std::atomic<int>> runs(n);
    std::vector<std::atomic<int>> done(n);
    std::atomic<bool> order{true};
    for (std::size_t i = 0; i < n; i++) {
      std::vector<std::size_t> after;
      // Each task waits for up to three earlier ones, chosen by a fixed rule.
      for (std::size_t k = 1; k <= 3; k++)
        if (i >= 7 * k) after.push_back((i * 31 + k * 17) % (i - 6 * k));
      g.add("t" + std::to_string(i), [&, i, after](std::size_t) {
        for (std::size_t p : after)
          if (!done[p].load()) order = false;
        runs[i]++;
        // A little work of uneven size.
        volatile double x = 0;
        for (std::size_t k = 0; k < 50 * (i % 7); k++) x = x + std::sqrt(static_cast<double>(k));
        done[i] = 1;
      }, after, 1 + static_cast<double>(i % 7));
    }
    const Trace t = g.run(pool);
    bool once = true;
    for (auto& r : runs) once = once && r.load() == 1;
    s.require("on " + std::to_string(threads) + " threads: every task ran once", once && t.records.size() == n);
    s.require("on " + std::to_string(threads) + " threads: every task started after its predecessors ended", order.load());
    bool ordered = true;
    for (const auto& r : t.records) ordered = ordered && r.start <= r.end && r.end <= t.seconds && r.thread < threads;
    s.require("on " + std::to_string(threads) + " threads: the trace's times are ordered, within the run, on the pool's threads", ordered);
    s.require("on " + std::to_string(threads) + " threads: utilisation lies in (0, 1]", t.utilisation() > 0 && t.utilisation() <= 1 + 1e-9,
              std::to_string(t.utilisation()));
  }

  // Priority: on one thread, of two ready tasks the one of higher rank starts first.
  {
    ThreadPool pool(1);
    TaskGraph g;
    std::vector<std::string> seen;
    g.add("cheap", [&](std::size_t) { seen.push_back("cheap"); }, {}, 1);
    const auto heavy = g.add("heavy", [&](std::size_t) { seen.push_back("heavy"); }, {}, 2);
    g.add("tail", [&](std::size_t) { seen.push_back("tail"); }, {heavy}, 5);
    g.run(pool);
    s.require("the ready task on the critical path goes first", seen == std::vector<std::string>{"heavy", "tail", "cheap"});
  }

  // Shards and a merge: the sum of squares below 10^6 by 64 shards, merged in order.
  {
    const std::uint64_t N = 1'000'000, shards = 64;
    std::uint64_t serial = 0;
    for (std::uint64_t i = 0; i < N; i++) serial += i * i;
    bool same = true;
    for (std::size_t threads : {1, 2, 4, 8}) {
      ThreadPool pool(threads);
      TaskGraph g;
      std::vector<std::uint64_t> part(shards, 0);
      std::vector<std::size_t> ids;
      for (std::uint64_t k = 0; k < shards; k++)
        ids.push_back(g.add("shard", [&, k](std::size_t) {
          const auto r = calx::scheduler::part(N, shards, k);
          for (std::uint64_t i = r.begin; i < r.end; i++) part[k] += i * i;
        }));
      std::uint64_t total = 0;
      g.add("merge", [&](std::size_t) { for (auto p : part) total += p; }, ids);
      g.run(pool);
      same = same && total == serial;
    }
    s.require("64 shards and a merge give the serial sum on 1, 2, 4 and 8 threads", same);
  }

  // An exception reaches the caller; the pool works on.
  {
    ThreadPool pool(4);
    TaskGraph g;
    std::atomic<int> after_failure{0};
    const auto bad = g.add("bad", [](std::size_t) { throw std::runtime_error("in a task"); });
    g.add("after", [&](std::size_t) { after_failure++; }, {bad});
    bool caught = false;
    try { g.run(pool); } catch (const std::runtime_error& e) { caught = std::string(e.what()) == "in a task"; }
    s.require("an exception in a task reaches the caller, and its successors do not run", caught && after_failure.load() == 0);
    TaskGraph h;
    int ran = 0;
    h.add("fine", [&](std::size_t) { ran++; });
    h.run(pool);
    s.require("the pool works on after a failed graph", ran == 1);
  }

  // parallel_for: every iteration once, in chunks of the grain.
  for (std::size_t threads : {1, 3, 8}) {
    ThreadPool pool(threads);
    for (std::size_t n : {0, 1, 7, 1000, 4097}) {
      std::vector<int> hits(n, 0);
      std::atomic<bool> chunks_ok{true};
      std::mutex m;
      parallel_for(pool, n, 64, [&](std::size_t b, std::size_t e, std::size_t t) {
        if (b % 64 != 0 || e - b > 64 || e > n || t >= threads) chunks_ok = false;
        std::lock_guard<std::mutex> lock(m);
        for (std::size_t i = b; i < e; i++) hits[i]++;
      });
      bool once = std::all_of(hits.begin(), hits.end(), [](int h) { return h == 1; });
      s.require("parallel_for on " + std::to_string(threads) + " threads covers [0, " + std::to_string(n) + ") once, in chunks of the grain", once && chunks_ok.load());
    }
  }

  return s.done();
}
