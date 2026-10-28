// Real-time scheduling: periodic and sporadic tasks with deadlines, the analyses that say whether
// a set of them can meet every deadline, servers that give a stream of work a guaranteed share of
// the processor and no more, a discrete-event simulator of a schedule in logical time, and a
// dispatcher that runs released jobs on the host's threads in the order of a policy.
//
// The model (Liu and Layland 1973, "Scheduling algorithms for multiprogramming in a hard-real-time
// environment", JACM 20, 46-61; Mok 1983 for sporadic tasks; Buttazzo 2011, Hard Real-Time
// Computing Systems, 3rd ed., Springer, ch. 4-6). A task τᵢ = (Cᵢ, Tᵢ, Dᵢ, φᵢ) releases a job every
// Tᵢ from its phase φᵢ (a sporadic task: at least Tᵢ apart), each needing Cᵢ of processor time
// before its absolute deadline, the release plus Dᵢ. Times are seconds of one declared clock,
// logical or of the wall; the module does not care which.
//
// What is here, with its source:
//   utilisation U = Σ Cᵢ/Tᵢ and density Σ Cᵢ/min(Dᵢ, Tᵢ);
//   the Liu–Layland bound n(2^{1/n} − 1) for rate-monotonic priorities with Dᵢ = Tᵢ (Liu and
//     Layland 1973, Theorem 5): sufficient, not necessary;
//   response-time analysis for fixed priorities, Rᵢ = Cᵢ + Bᵢ + Σ_{j ∈ hp(i)} ⌈Rᵢ/Tⱼ⌉ Cⱼ iterated
//     to its least fixed point (Joseph and Pandya 1986, The Computer Journal 29, 390-395; Audsley,
//     Burns, Richardson, Tindell and Wellings 1993, Software Engineering Journal 8, 284-292), with
//     Bᵢ the blocking by a lower-priority job that cannot be preempted: exact for preemptive
//     fixed priorities with Dᵢ ≤ Tᵢ and synchronous release, the critical instant (Liu and Layland
//     1973, Theorem 1);
//   EDF's feasibility on one processor: U ≤ 1 when every Dᵢ = Tᵢ (Liu and Layland 1973,
//     Theorem 7), and otherwise the processor-demand criterion dbf(t) ≤ t at every absolute
//     deadline t up to the synchronous busy period (Baruah, Rosier and Howell 1990, Real-Time
//     Systems 2, 301-324; Spuri 1996, INRIA RR-2772, for the busy period);
//   the constant bandwidth server (Abeni and Buttazzo 1998, "Integrating multimedia applications in
//     hard real-time systems", RTSS, 4-13): a stream of jobs of unknown cost served with a budget
//     Q every period T demands at most Q/T of the processor whatever the jobs do, so that EDF with
//     servers meets every hard deadline when Σ Cᵢ/Tᵢ + Σ Qₛ/Tₛ ≤ 1 (their Theorem 1, isolation);
//   aging for fixed priorities, the effective priority rising by a level for each `aging_step`
//     waited, which bounds how long a ready job can wait (no starvation) at the price of the exact
//     analysis above;
//   elastic compression of periods to bring a set's utilisation to a desired value (Buttazzo,
//     Lipari, Caccamo and Abeni 2002, "Elastic scheduling for flexible workload management", IEEE
//     Transactions on Computers 51, 289-302, §3), the way a set is degraded under overload;
//   simulate(): a discrete-event simulation of jobs on m identical processors, preemptive or not,
//     by EDF, fixed priority (with or without aging) or FIFO, with servers; deterministic, ties
//     going to the earlier release, then the lower task index, then the earlier job;
//   dispatch(): released jobs run on a ThreadPool, a free thread taking the ready job first in the
//     policy's order and running it to its end (user code cannot be preempted), with a record of
//     each job's start, end and lateness.
//
// Caveat. The analyses are of the model: they hold for the costs declared. On a general-purpose
// operating system a thread can be delayed by the kernel, so dispatch() measures lateness rather
// than guaranteeing it; a guarantee needs a real-time kernel (Linux's SCHED_DEADLINE, which is EDF
// with constant bandwidth servers, or PREEMPT_RT with fixed priorities) beneath it.

#pragma once

#include <calx/scheduler/thread_pool.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <cstddef>
#include <exception>
#include <functional>
#include <limits>
#include <mutex>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace calx::scheduler {

inline constexpr std::size_t no_server = std::numeric_limits<std::size_t>::max();

// A periodic or sporadic task. deadline 0 means the period (an implicit deadline); priority is for
// fixed-priority policies, a smaller number first; criticality is the caller's (a larger number
// more critical), carried to jobs and to the records for degradation and reports.
struct PeriodicTask {
  std::string name;
  double cost = 0;      // C, s
  double period = 0;    // T, s (a sporadic task's least separation)
  double deadline = 0;  // D, s; 0: D = T
  double phase = 0;     // φ, s: the first release
  int priority = 0;
  int criticality = 0;
  std::size_t server = no_server;  // served by this server (EDF with servers), or not
  double relative_deadline() const { return deadline > 0 ? deadline : period; }
};

// ---- analysis ---------------------------------------------------------------------------------------

inline double utilisation(const std::vector<PeriodicTask>& tasks) {
  double u = 0;
  for (const auto& t : tasks) u += t.cost / t.period;
  return u;
}

inline double density(const std::vector<PeriodicTask>& tasks) {
  double u = 0;
  for (const auto& t : tasks) u += t.cost / std::min(t.relative_deadline(), t.period);
  return u;
}

// n (2^{1/n} − 1): the utilisation below which rate-monotonic priorities meet every implicit
// deadline of n tasks (Liu and Layland 1973, Theorem 5); 1 for one task, ln 2 as n grows.
inline double liu_layland_bound(std::size_t n) {
  if (n == 0) return 1;
  return static_cast<double>(n) * (std::pow(2.0, 1.0 / static_cast<double>(n)) - 1.0);
}

// The tasks' indices in rate-monotonic order (the shorter period first; Liu and Layland 1973) or
// deadline-monotonic order (the shorter relative deadline first; Leung and Whitehead 1982), ties to
// the lower index.
inline std::vector<std::size_t> rate_monotonic(const std::vector<PeriodicTask>& tasks, bool by_deadline = false) {
  std::vector<std::size_t> order(tasks.size());
  for (std::size_t i = 0; i < order.size(); i++) order[i] = i;
  std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
    const double ka = by_deadline ? tasks[a].relative_deadline() : tasks[a].period;
    const double kb = by_deadline ? tasks[b].relative_deadline() : tasks[b].period;
    return ka < kb;
  });
  return order;
}

// The worst-case response time of each task under preemptive fixed priorities, the tasks given in
// priority order (the first highest), with an optional blocking term per task (the longest job of
// lower priority that cannot be preempted, for non-preemptive or limited-preemptive schedules).
// Rᵢ is the least fixed point of R = Cᵢ + Bᵢ + Σ_{j<i} ⌈R/Tⱼ⌉ Cⱼ (Joseph and Pandya 1986), iterated
// from Cᵢ + Bᵢ; when R exceeds `limit` × Dᵢ the iteration stops and the task is reported
// unschedulable (nullopt for that task). Exact for D ≤ T with synchronous release.
inline std::vector<std::optional<double>> response_times(const std::vector<PeriodicTask>& by_priority,
                                                         const std::vector<double>& blocking = {}, double limit = 1.0) {
  std::vector<std::optional<double>> out(by_priority.size());
  for (std::size_t i = 0; i < by_priority.size(); i++) {
    const auto& ti = by_priority[i];
    const double b = i < blocking.size() ? blocking[i] : 0.0;
    const double bound = limit * ti.relative_deadline();
    double r = ti.cost + b;
    for (int iteration = 0; iteration < 100000; iteration++) {
      double next = ti.cost + b;
      for (std::size_t j = 0; j < i; j++) next += std::ceil(r / by_priority[j].period - 1e-12) * by_priority[j].cost;
      if (next > bound + 1e-12) { r = -1; break; }
      if (std::fabs(next - r) <= 1e-12 * std::max(1.0, r)) { r = next; break; }
      r = next;
    }
    if (r >= 0) out[i] = r;
  }
  return out;
}

// The length of the synchronous busy period, the least L > 0 with L = Σ ⌈L/Tᵢ⌉ Cᵢ (Spuri 1996),
// when U < 1 (or U = 1 with the hyperperiod as its bound); nullopt when U > 1.
inline std::optional<double> busy_period(const std::vector<PeriodicTask>& tasks) {
  const double u = utilisation(tasks);
  if (tasks.empty()) return 0.0;
  if (u > 1 + 1e-12) return std::nullopt;
  double l = 0;
  for (const auto& t : tasks) l += t.cost;
  for (int iteration = 0; iteration < 1000000; iteration++) {
    double next = 0;
    for (const auto& t : tasks) next += std::ceil(l / t.period - 1e-12) * t.cost;
    if (std::fabs(next - l) <= 1e-12 * std::max(1.0, l)) return next;
    l = next;
  }
  return std::nullopt;
}

// The processor demand of the tasks' jobs released and due within [0, t] under synchronous release:
// dbf(t) = Σᵢ max(0, ⌊(t − Dᵢ)/Tᵢ⌋ + 1) Cᵢ (Baruah, Rosier and Howell 1990).
inline double demand_bound(const std::vector<PeriodicTask>& tasks, double t) {
  double d = 0;
  for (const auto& task : tasks) {
    const double k = std::floor((t - task.relative_deadline()) / task.period + 1e-12) + 1;
    if (k > 0) d += k * task.cost;
  }
  return d;
}

// Whether preemptive EDF on one processor meets every deadline of the tasks, phases ignored (the
// synchronous case is the worst). Implicit deadlines: U ≤ 1 (Liu and Layland 1973, Theorem 7).
// Otherwise U ≤ 1 and dbf(t) ≤ t at every absolute deadline t within the busy period (Baruah,
// Rosier and Howell 1990, Theorem 3.1 with Spuri's bound).
struct EdfVerdict {
  bool feasible = false;
  double utilisation = 0;
  double horizon = 0;           // the deadlines tested lie in (0, horizon]
  double first_failure = -1;    // the first t with dbf(t) > t, or −1
  std::size_t points = 0;       // the number of deadlines tested
};

inline EdfVerdict edf_feasible(const std::vector<PeriodicTask>& tasks, std::size_t max_points = 2000000) {
  EdfVerdict v;
  v.utilisation = utilisation(tasks);
  if (v.utilisation > 1 + 1e-12) return v;
  bool implicit = true;
  for (const auto& t : tasks) implicit = implicit && std::fabs(t.relative_deadline() - t.period) <= 1e-15 * t.period;
  if (implicit) { v.feasible = true; return v; }
  const auto l = busy_period(tasks);
  if (!l) return v;
  v.horizon = *l;
  // Every absolute deadline kTᵢ + Dᵢ ≤ L, in increasing order, merged across the tasks.
  std::vector<double> next(tasks.size());
  for (std::size_t i = 0; i < tasks.size(); i++) next[i] = tasks[i].relative_deadline();
  for (;;) {
    double t = std::numeric_limits<double>::infinity();
    for (double x : next) t = std::min(t, x);
    if (!(t <= v.horizon + 1e-12)) break;
    if (++v.points > max_points) return v;  // too many to test: not shown feasible
    if (demand_bound(tasks, t) > t + 1e-9 * std::max(1.0, t)) { v.first_failure = t; return v; }
    for (std::size_t i = 0; i < tasks.size(); i++)
      if (next[i] <= t + 1e-12) next[i] += tasks[i].period;
  }
  v.feasible = true;
  return v;
}

// ---- elastic compression ----------------------------------------------------------------------------

// A task's elasticity: its nominal period, the longest it may be stretched to, and its elastic
// coefficient E ≥ 0 (0: rigid). compress() returns periods with Σ C/T = desired, each task's
// utilisation reduced in proportion to its coefficient among those still free to move, those that
// reach their longest period fixed there and the rest compressed again (Buttazzo et al. 2002, §3,
// the iterative algorithm). nullopt when the desired utilisation cannot be reached even with every
// elastic task at its longest period.
struct Elastic {
  double cost = 0, period = 0, longest = 0, coefficient = 1;
};

inline std::optional<std::vector<double>> compress(const std::vector<Elastic>& tasks, double desired) {
  const std::size_t n = tasks.size();
  std::vector<double> u(n), u_min(n);
  double total = 0;
  for (std::size_t i = 0; i < n; i++) {
    u[i] = tasks[i].cost / tasks[i].period;
    u_min[i] = tasks[i].coefficient > 0 ? tasks[i].cost / std::max(tasks[i].longest, tasks[i].period) : u[i];
    total += u[i];
  }
  if (total <= desired + 1e-15) {
    std::vector<double> periods(n);
    for (std::size_t i = 0; i < n; i++) periods[i] = tasks[i].period;
    return periods;
  }
  std::vector<bool> fixed(n, false);
  std::vector<double> current = u;
  for (std::size_t round = 0; round <= n; round++) {
    double u_fixed = 0, u_free = 0, e_free = 0;
    for (std::size_t i = 0; i < n; i++) {
      if (fixed[i] || tasks[i].coefficient <= 0) u_fixed += fixed[i] ? current[i] : u[i];
      else { u_free += u[i]; e_free += tasks[i].coefficient; }
    }
    if (e_free <= 0) return std::nullopt;
    const double excess = u_free + u_fixed - desired;
    bool again = false;
    for (std::size_t i = 0; i < n; i++) {
      if (fixed[i] || tasks[i].coefficient <= 0) continue;
      current[i] = u[i] - excess * tasks[i].coefficient / e_free;
      if (current[i] < u_min[i] - 1e-15) { current[i] = u_min[i]; fixed[i] = true; again = true; }
    }
    if (!again) {
      std::vector<double> periods(n);
      for (std::size_t i = 0; i < n; i++) periods[i] = tasks[i].cost / current[i];
      return periods;
    }
  }
  return std::nullopt;
}

// ---- the constant bandwidth server ---------------------------------------------------------------

// A CBS with budget Q every period T (Abeni and Buttazzo 1998, §3, rules (1)-(4)). The server has
// a budget c and a deadline d. A job arriving at r when the server has no pending work keeps (c, d)
// if c < (d − r) Q/T, and otherwise gets d = r + T, c = Q (rule 2); the server's jobs are scheduled
// by EDF with the server's deadline; executing consumes c (rule 3), and when c reaches zero it is
// recharged to Q and d postponed by T (rule 4), so that the work it serves is pushed later rather
// than taking more than Q/T. The jobs it serves are served first come, first served.
class ConstantBandwidthServer {
 public:
  ConstantBandwidthServer(double budget = 0, double period = 1) : q_(budget), t_(period) {
    if (!(budget > 0) || !(period > 0) || budget > period) throw std::invalid_argument("calx::scheduler::ConstantBandwidthServer: 0 < Q ≤ T");
  }
  double budget() const { return q_; }
  double period() const { return t_; }
  double bandwidth() const { return q_ / t_; }
  double remaining() const { return c_; }
  double deadline() const { return d_; }

  // A job arrives at r; `idle` says whether the server had no pending work. Returns the deadline.
  double arrive(double r, bool idle) {
    if (idle && c_ >= (d_ - r) * bandwidth()) {
      d_ = r + t_;
      c_ = q_;
    }
    return d_;
  }
  // The server's jobs executed for dt; returns the (possibly postponed) deadline.
  double consume(double dt) {
    c_ -= dt;
    while (c_ <= 1e-15 * q_) {
      c_ += q_;
      d_ += t_;
    }
    return d_;
  }

 private:
  double q_, t_;
  double c_ = 0, d_ = 0;
};

// ---- simulation -------------------------------------------------------------------------------------

enum class Policy { edf, fixed_priority, fifo };

inline const char* policy_name(Policy p) {
  switch (p) {
    case Policy::edf: return "edf";
    case Policy::fixed_priority: return "fixed-priority";
    case Policy::fifo: return "fifo";
  }
  return "?";
}

// One job: released at `release`, needing `cost`, due at `deadline` (absolute). Jobs of a task
// carry its index; others (sporadic arrivals given directly) a task of their own choosing.
struct Job {
  std::size_t task = 0, index = 0;
  double release = 0, cost = 0, deadline = 0;
  int priority = 0, criticality = 0;
  std::size_t server = no_server;
};

// An interval in which a job ran on a processor.
struct Slice {
  std::size_t job = 0, processor = 0;
  double start = 0, end = 0;
};

// What became of a job.
struct Outcome {
  Job job;
  double start = -1, finish = -1;   // −1: not started, not finished by the horizon
  double deadline = 0;              // the deadline it was scheduled by (a server's may move)
  bool finished() const { return finish >= 0; }
  double response() const { return finished() ? finish - job.release : -1; }
  double lateness() const { return finished() ? finish - job.deadline : -1; }  // against its own deadline
  bool missed(double horizon) const {
    return finished() ? finish > job.deadline + 1e-9 * std::max(1.0, std::fabs(job.deadline)) : job.deadline < horizon;
  }
};

struct SimulationOptions {
  Policy policy = Policy::edf;
  std::size_t processors = 1;
  bool preemptive = true;
  double aging_step = 0;                         // fixed priority: a level for each step waited; 0: none
  std::vector<ConstantBandwidthServer> servers;  // EDF: the jobs with a server take its deadlines
};

struct Simulation {
  std::vector<Slice> slices;       // in order of start
  std::vector<Outcome> outcomes;   // one per job, in the order of the jobs given
  double horizon = 0;
  std::size_t misses = 0;          // jobs past their own deadline (or due before the horizon and unfinished)
  double max_lateness = -std::numeric_limits<double>::infinity();
  double busy = 0;                 // processor time spent in jobs
  // Processor time given to the jobs of a server within [0, horizon].
  double served(std::size_t server) const {
    double s = 0;
    for (const auto& sl : slices)
      if (outcomes[sl.job].job.server == server) s += sl.end - sl.start;
    return s;
  }
};

// The jobs of periodic tasks released in [0, horizon), each task at phase + k T.
inline std::vector<Job> periodic_jobs(const std::vector<PeriodicTask>& tasks, double horizon) {
  std::vector<Job> jobs;
  for (std::size_t i = 0; i < tasks.size(); i++) {
    const auto& t = tasks[i];
    if (!(t.period > 0)) throw std::invalid_argument("calx::scheduler::periodic_jobs: a period is positive");
    for (std::size_t k = 0;; k++) {
      const double r = t.phase + static_cast<double>(k) * t.period;
      if (r >= horizon - 1e-12) break;
      jobs.push_back({i, k, r, t.cost, r + t.relative_deadline(), t.priority, t.criticality, t.server});
    }
  }
  return jobs;
}

// Simulates the jobs on the processors until every job has finished or the horizon is reached.
inline Simulation simulate(const std::vector<Job>& jobs, double horizon, SimulationOptions options = {}) {
  if (options.processors == 0) throw std::invalid_argument("calx::scheduler::simulate: at least one processor");
  for (const auto& j : jobs)
    if (j.server != no_server && (options.policy != Policy::edf || j.server >= options.servers.size()))
      throw std::invalid_argument("calx::scheduler::simulate: a job's server is one of the options' servers, under EDF");
  Simulation sim;
  sim.horizon = horizon;
  const std::size_t n = jobs.size();
  sim.outcomes.resize(n);
  for (std::size_t i = 0; i < n; i++) { sim.outcomes[i].job = jobs[i]; sim.outcomes[i].deadline = jobs[i].deadline; }

  std::vector<std::size_t> order(n);  // by release, then task, then index
  for (std::size_t i = 0; i < n; i++) order[i] = i;
  std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
    const Job &x = jobs[a], &y = jobs[b];
    if (x.release != y.release) return x.release < y.release;
    if (x.task != y.task) return x.task < y.task;
    return x.index < y.index;
  });
  std::vector<double> remaining(n);
  for (std::size_t i = 0; i < n; i++) remaining[i] = jobs[i].cost;
  std::vector<std::size_t> ready;                       // released, unfinished, not running
  std::vector<std::size_t> running(options.processors, n);  // n: idle
  std::vector<double> since(options.processors, 0);         // when the current slice began
  // Servers: their pending jobs in arrival order (the head is the one the server offers).
  std::vector<std::vector<std::size_t>> queue(options.servers.size());
  auto& servers = options.servers;

  const auto key_less = [&](std::size_t a, std::size_t b, double t) {
    const Job &x = jobs[a], &y = jobs[b];
    if (options.policy == Policy::edf) {
      const double da = sim.outcomes[a].deadline, db = sim.outcomes[b].deadline;
      if (da != db) return da < db;
    } else if (options.policy == Policy::fixed_priority) {
      double pa = x.priority, pb = y.priority;
      if (options.aging_step > 0) {
        pa -= std::floor((t - x.release) / options.aging_step + 1e-12);
        pb -= std::floor((t - y.release) / options.aging_step + 1e-12);
      }
      if (pa != pb) return pa < pb;
    }
    if (x.release != y.release) return x.release < y.release;
    if (x.task != y.task) return x.task < y.task;
    return x.index < y.index;
  };
  // A server's job is eligible only at the head of its server's queue.
  const auto eligible = [&](std::size_t j) {
    const std::size_t s = jobs[j].server;
    return s == no_server || (!queue[s].empty() && queue[s].front() == j);
  };

  std::size_t next_release = 0, finished = 0;
  double t = 0;
  const auto close_slice = [&](std::size_t p, double now) {
    const std::size_t j = running[p];
    if (j == n || now <= since[p]) return;
    sim.slices.push_back({j, p, since[p], now});
    sim.busy += now - since[p];
    since[p] = now;
  };
  const auto charge = [&](std::size_t j, double dt) {
    remaining[j] -= dt;
    const std::size_t s = jobs[j].server;
    if (s != no_server && dt > 0) {
      const double d = servers[s].consume(dt);
      for (std::size_t k : queue[s]) sim.outcomes[k].deadline = d;
    }
  };

  for (int guard = 0; guard < 50000000; guard++) {
    // Releases at t.
    while (next_release < n && jobs[order[next_release]].release <= t + 1e-12) {
      const std::size_t j = order[next_release++];
      const std::size_t s = jobs[j].server;
      if (s != no_server) {
        const bool idle = queue[s].empty();
        queue[s].push_back(j);
        sim.outcomes[j].deadline = idle ? servers[s].arrive(jobs[j].release, true) : servers[s].deadline();
      }
      ready.push_back(j);
    }
    // Completions at t.
    for (std::size_t p = 0; p < running.size(); p++) {
      const std::size_t j = running[p];
      if (j != n && remaining[j] <= 1e-12 * std::max(1.0, jobs[j].cost)) {
        close_slice(p, t);
        finished++;
        sim.outcomes[j].finish = t;
        running[p] = n;
        const std::size_t s = jobs[j].server;
        if (s != no_server) {
          queue[s].erase(queue[s].begin());
          if (!queue[s].empty()) {
            // The next pending job takes the server's current deadline (it is not idle).
            for (std::size_t k : queue[s]) sim.outcomes[k].deadline = servers[s].deadline();
          }
        }
      }
    }
    if (finished == n || t >= horizon - 1e-12) break;

    // Choose what runs: the eligible jobs first in the policy's order.
    std::vector<std::size_t> candidates;
    for (std::size_t j : ready)
      if (eligible(j)) candidates.push_back(j);
    if (options.preemptive)
      for (std::size_t p = 0; p < running.size(); p++)
        if (running[p] != n) candidates.push_back(running[p]);
    std::sort(candidates.begin(), candidates.end(), [&](std::size_t a, std::size_t b) { return key_less(a, b, t); });
    if (options.preemptive) {
      std::vector<std::size_t> chosen(candidates.begin(), candidates.begin() + std::min(candidates.size(), running.size()));
      // Keep a chosen job on its processor; preempt the rest.
      std::vector<bool> keep(running.size(), false);
      for (std::size_t p = 0; p < running.size(); p++)
        if (running[p] != n && std::find(chosen.begin(), chosen.end(), running[p]) != chosen.end()) keep[p] = true;
      for (std::size_t p = 0; p < running.size(); p++)
        if (running[p] != n && !keep[p]) { close_slice(p, t); ready.push_back(running[p]); running[p] = n; }
      for (std::size_t j : chosen) {
        if (std::find(running.begin(), running.end(), j) != running.end()) continue;
        for (std::size_t p = 0; p < running.size(); p++)
          if (running[p] == n) {
            running[p] = j; since[p] = t;
            if (sim.outcomes[j].start < 0) sim.outcomes[j].start = t;
            ready.erase(std::find(ready.begin(), ready.end(), j));
            break;
          }
      }
    } else {
      std::size_t c = 0;
      for (std::size_t p = 0; p < running.size() && c < candidates.size(); p++)
        if (running[p] == n) {
          const std::size_t j = candidates[c++];
          running[p] = j; since[p] = t;
          if (sim.outcomes[j].start < 0) sim.outcomes[j].start = t;
          ready.erase(std::find(ready.begin(), ready.end(), j));
        }
    }

    // The next event: a release, a completion, a server's budget running out, an aging step.
    double next = horizon;
    if (next_release < n) next = std::min(next, jobs[order[next_release]].release);
    for (std::size_t p = 0; p < running.size(); p++) {
      const std::size_t j = running[p];
      if (j == n) continue;
      next = std::min(next, t + remaining[j]);
      const std::size_t s = jobs[j].server;
      if (s != no_server && options.preemptive) next = std::min(next, t + servers[s].remaining());
    }
    if (options.policy == Policy::fixed_priority && options.aging_step > 0 && options.preemptive)
      for (std::size_t j : ready) {
        const double waited = t - jobs[j].release;
        next = std::min(next, jobs[j].release + (std::floor(waited / options.aging_step + 1e-12) + 1) * options.aging_step);
      }
    if (next <= t) next = t + 1e-12 * std::max(1.0, t);  // guard against no progress
    // Run until then.
    const double dt = next - t;
    for (std::size_t p = 0; p < running.size(); p++)
      if (running[p] != n) charge(running[p], dt);
    t = next;
    // A server whose budget ran out postponed its deadline in charge(); its running job is reconsidered
    // at the next decision, as a preemptive EDF would.
  }
  for (std::size_t p = 0; p < running.size(); p++) close_slice(p, std::min(t, horizon));
  std::sort(sim.slices.begin(), sim.slices.end(), [](const Slice& a, const Slice& b) { return a.start != b.start ? a.start < b.start : a.processor < b.processor; });
  for (const auto& o : sim.outcomes) {
    if (o.missed(horizon)) sim.misses++;
    if (o.finished()) sim.max_lateness = std::max(sim.max_lateness, o.lateness());
  }
  return sim;
}

// ---- dispatch on the host's threads ----------------------------------------------------------------

// A job to run now: its deadline and priority in the policy's terms (seconds of the dispatch's clock,
// which starts at the call), its work, and a group: the jobs of one group run in their order on one
// thread at a time (for work that shares a resource and must not overlap).
struct LiveJob {
  std::string name;
  double deadline = std::numeric_limits<double>::infinity();
  int priority = 0, criticality = 0;
  std::size_t group = no_server;
  std::function<void(std::size_t thread)> work;
};

struct LiveRecord {
  std::size_t job = 0, thread = 0;
  double start = 0, end = 0;
  double lateness(const LiveJob& j) const { return end - j.deadline; }
};

struct Dispatch {
  std::size_t threads = 1;
  double seconds = 0;
  std::vector<LiveRecord> records;  // in the order the jobs started
  std::size_t late = 0;             // jobs that ended after their deadline
};

// Runs every job once on the pool, a free thread taking the ready job first in the policy's order
// (EDF: the earliest deadline; fixed priority: the smallest priority; FIFO: the order given), ties to
// the order given, and running it to its end. `epoch` is the instant the deadlines count from (by
// default the call). An exception in a job stops new jobs, lets those running finish, and reaches the
// caller.
inline Dispatch dispatch(ThreadPool& pool, const std::vector<LiveJob>& jobs, Policy policy = Policy::edf,
                         std::optional<std::chrono::steady_clock::time_point> epoch = std::nullopt) {
  Dispatch out;
  out.threads = pool.size();
  const std::size_t n = jobs.size();
  if (n == 0) return out;
  std::vector<std::size_t> order(n);
  for (std::size_t i = 0; i < n; i++) order[i] = i;
  std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
    if (policy == Policy::edf && jobs[a].deadline != jobs[b].deadline) return jobs[a].deadline < jobs[b].deadline;
    if (policy == Policy::fixed_priority && jobs[a].priority != jobs[b].priority) return jobs[a].priority < jobs[b].priority;
    return false;
  });
  std::mutex mutex;
  std::vector<bool> taken(n, false);
  std::vector<bool> busy_group;  // groups running now
  std::size_t left = n;
  std::exception_ptr error;
  std::condition_variable changed;
  const auto t0 = epoch.value_or(std::chrono::steady_clock::now());
  const auto now = [&] { return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(); };
  std::size_t groups = 0;
  for (const auto& j : jobs)
    if (j.group != no_server) groups = std::max(groups, j.group + 1);
  busy_group.assign(groups, false);
  out.records.reserve(n);
  pool.run([&](std::size_t thread) {
    std::unique_lock<std::mutex> lock(mutex);
    for (;;) {
      std::size_t pick = n;
      changed.wait(lock, [&] {
        if (left == 0 || error) return true;
        for (std::size_t k : order)
          if (!taken[k] && (jobs[k].group == no_server || !busy_group[jobs[k].group])) { pick = k; return true; }
        return false;
      });
      if (left == 0 || error || pick == n) return;
      taken[pick] = true;
      if (jobs[pick].group != no_server) busy_group[jobs[pick].group] = true;
      const std::size_t slot = out.records.size();
      out.records.push_back({pick, thread, now(), 0});
      lock.unlock();
      std::exception_ptr failure;
      try {
        if (jobs[pick].work) jobs[pick].work(thread);
      } catch (...) {
        failure = std::current_exception();
      }
      const double end = now();
      lock.lock();
      out.records[slot].end = end;
      if (jobs[pick].group != no_server) busy_group[jobs[pick].group] = false;
      if (failure && !error) error = failure;
      left--;
      changed.notify_all();
    }
  });
  for (const auto& r : out.records) {
    out.seconds = std::max(out.seconds, r.end);
    if (r.end > jobs[r.job].deadline) out.late++;
  }
  if (error) std::rethrow_exception(error);
  return out;
}

}  // namespace calx::scheduler
