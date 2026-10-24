// Task graphs on the host: a set of tasks, each a piece of work with a declared cost, and the
// tasks it must wait for; run on the ThreadPool with every task started as soon as all its
// predecessors have finished, and a trace of where and when each one ran.
//
// The graph is acyclic by construction: a task may wait only for tasks added before it, so
// the order of addition is a topological order.
//
// Among the tasks that are ready, a free thread takes the one of highest upward rank, the
// declared cost of the longest path from it to the end of the graph: the priority of HEFT
// (Topcuoglu, Hariri and Wu 2002, "Performance-effective and low-complexity task scheduling
// for heterogeneous computing", IEEE TPDS 13, 260-274, §4.1), here with one kind of
// processor. Ties go to the task added first. With costs that are only estimates the
// schedule is a good one, not an optimal one (the problem is NP-hard); it is never wrong,
// since dependencies are always respected.
//
//   TaskGraph g;
//   std::vector<Partial> part(shards);
//   std::vector<std::size_t> ids;
//   for (std::size_t k = 0; k < shards; k++)
//     ids.push_back(g.add("shard " + std::to_string(k), [&, k](std::size_t) { part[k] = work(k); }, {}, cost(k)));
//   g.add("merge", [&](std::size_t) { result = merge(part); }, ids);
//   Trace trace = g.run(pool);
//
// Determinism is the caller's: a task that writes only its own outputs, and a merge that reads
// them in a fixed order, give the same result on any number of threads, whatever the trace.
// An exception in a task stops new tasks from starting, lets those running finish, and
// reaches the caller of run().
//
// This is the first representation of task graphs in the scheduler (README §3): all tasks on
// the host's threads. Placement on other devices adds a cost per task and device, as HEFT
// does, and data movement between them.

#pragma once

#include <calx/scheduler/thread_pool.hpp>

#include <algorithm>
#include <chrono>
#include <condition_variable>
#include <cstddef>
#include <exception>
#include <functional>
#include <mutex>
#include <queue>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace calx::scheduler {

// Where and when one task ran: seconds from the start of TaskGraph::run.
struct TaskRecord {
  std::size_t task = 0, thread = 0;
  double start = 0, end = 0;
};

struct Trace {
  std::size_t threads = 1;
  double seconds = 0;               // from the start of the run to the end of the last task
  std::vector<TaskRecord> records;  // in the order the tasks started

  // The time the threads spent in tasks.
  double busy() const {
    double b = 0;
    for (const auto& r : records) b += r.end - r.start;
    return b;
  }
  // The fraction of the threads' time spent in tasks.
  double utilisation() const { return seconds > 0 ? busy() / (seconds * static_cast<double>(threads)) : 0; }
  // The time in tasks of one thread.
  double busy(std::size_t thread) const {
    double b = 0;
    for (const auto& r : records)
      if (r.thread == thread) b += r.end - r.start;
    return b;
  }
};

class TaskGraph {
 public:
  using Work = std::function<void(std::size_t thread)>;

  // Adds a task that runs after every task in `after`, all of which must already be in the
  // graph; returns its index. The cost is in any unit, the same for every task.
  std::size_t add(std::string name, Work work, std::vector<std::size_t> after = {}, double cost = 1) {
    const std::size_t id = tasks_.size();
    for (std::size_t p : after)
      if (p >= id) throw std::invalid_argument("calx::scheduler::TaskGraph::add: a task waits only for tasks added before it");
    if (!(cost >= 0)) throw std::invalid_argument("calx::scheduler::TaskGraph::add: a cost is not negative");
    std::sort(after.begin(), after.end());
    after.erase(std::unique(after.begin(), after.end()), after.end());
    tasks_.push_back({std::move(name), std::move(work), std::move(after), cost});
    return id;
  }

  std::size_t size() const { return tasks_.size(); }
  const std::string& name(std::size_t i) const { return tasks_.at(i).name; }
  const std::vector<std::size_t>& after(std::size_t i) const { return tasks_.at(i).after; }
  double cost(std::size_t i) const { return tasks_.at(i).cost; }

  // The tasks that wait for each task.
  std::vector<std::vector<std::size_t>> successors() const {
    std::vector<std::vector<std::size_t>> s(tasks_.size());
    for (std::size_t i = 0; i < tasks_.size(); i++)
      for (std::size_t p : tasks_[i].after) s[p].push_back(i);
    return s;
  }

  // rank(i) = cost(i) + the largest rank of a successor: the declared cost of the longest path
  // from i to the end. One pass in reverse order of addition, a reverse topological order.
  std::vector<double> upward_rank() const {
    const auto succ = successors();
    std::vector<double> rank(tasks_.size(), 0);
    for (std::size_t k = tasks_.size(); k-- > 0;) {
      double longest = 0;
      for (std::size_t j : succ[k]) longest = std::max(longest, rank[j]);
      rank[k] = tasks_[k].cost + longest;
    }
    return rank;
  }

  // The declared cost of the longest path through the graph, a lower bound on the time of any
  // schedule when costs are times.
  double critical_path() const {
    const auto rank = upward_rank();
    return rank.empty() ? 0 : *std::max_element(rank.begin(), rank.end());
  }

  // Runs every task once, each after its predecessors, on the pool's threads.
  Trace run(ThreadPool& pool) const {
    Trace trace;
    trace.threads = pool.size();
    const std::size_t n = tasks_.size();
    if (n == 0) return trace;
    const auto succ = successors();
    const auto rank = upward_rank();
    std::vector<std::size_t> waiting(n);
    for (std::size_t i = 0; i < n; i++) waiting[i] = tasks_[i].after.size();

    // Ready tasks: the highest rank first, then the first added.
    const auto lower = [&](std::size_t a, std::size_t b) { return rank[a] != rank[b] ? rank[a] < rank[b] : a > b; };
    std::priority_queue<std::size_t, std::vector<std::size_t>, decltype(lower)> ready(lower);
    for (std::size_t i = 0; i < n; i++)
      if (waiting[i] == 0) ready.push(i);

    std::mutex mutex;
    std::condition_variable changed;
    std::size_t unfinished = n;
    std::exception_ptr error;
    trace.records.reserve(n);
    const auto t0 = std::chrono::steady_clock::now();
    const auto now = [&] { return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(); };

    pool.run([&](std::size_t thread) {
      std::unique_lock<std::mutex> lock(mutex);
      for (;;) {
        changed.wait(lock, [&] { return !ready.empty() || unfinished == 0 || error; });
        if (unfinished == 0 || error) return;
        const std::size_t task = ready.top();
        ready.pop();
        const std::size_t slot = trace.records.size();
        trace.records.push_back({task, thread, now(), 0});
        lock.unlock();
        std::exception_ptr failure;
        try {
          if (tasks_[task].work) tasks_[task].work(thread);
        } catch (...) {
          failure = std::current_exception();
        }
        const double end = now();
        lock.lock();
        trace.records[slot].end = end;
        if (failure && !error) error = failure;
        unfinished--;
        for (std::size_t j : succ[task])
          if (--waiting[j] == 0) ready.push(j);
        changed.notify_all();
      }
    });
    for (const auto& r : trace.records) trace.seconds = std::max(trace.seconds, r.end);
    if (error) std::rethrow_exception(error);
    return trace;
  }

 private:
  struct Task {
    std::string name;
    Work work;
    std::vector<std::size_t> after;
    double cost = 1;
  };
  std::vector<Task> tasks_;
};

}  // namespace calx::scheduler
