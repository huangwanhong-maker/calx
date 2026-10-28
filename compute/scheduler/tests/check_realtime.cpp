// Checks of real-time scheduling: each analysis against the simulation, which is an independent
// computation of the same schedule.
//
//   the Liu–Layland bound at n = 1, 2 and its limit ln 2;
//   response-time analysis against preemptive fixed priorities simulated over the hyperperiod from a
//     synchronous release, on random task sets: the first job's response equals Rᵢ (the critical
//     instant, Liu and Layland 1973, Theorem 1), no job's exceeds it, and the first task the
//     analysis rejects misses its first deadline;
//   EDF: implicit-deadline sets with U ≤ 1 meet every deadline, sets with U > 1 miss one within the
//     hyperperiod; the processor-demand criterion agrees with the simulation on random
//     constrained-deadline sets; EDF meets every deadline of every set rate-monotonic priorities do;
//   the constant bandwidth server's rules on a sequence worked by hand, and its isolation: hard tasks
//     keep their deadlines beside a server whose jobs demand eight times its budget, and the server
//     receives no more than its bandwidth; without the server the same jobs make the hard tasks miss;
//   aging: a low-priority job beside a stream that fills the processor never runs without aging, and
//     runs when its effective priority reaches the stream's (step 1, level 5: at t = 5);
//   non-preemptive runs give each job one slice; global EDF on two processors never runs more than
//     two jobs and never leaves a processor idle while a job waits (work-conserving);
//   elastic compression reaches the desired utilisation within the periods' bounds, leaves rigid tasks
//     alone, clamps a task at its longest period and refuses what cannot be reached;
//   dispatch on threads: one thread runs the jobs in EDF's order; four run each job once; jobs of one
//     group never overlap; jobs due before they start are counted late.

#include <calx/scheduler/realtime.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <random>
#include <string>
#include <vector>

using namespace calx::scheduler;

namespace {

// Integers from the engine's own output, so that the sets are the same on every platform (the
// standard fixes mt19937_64's sequence but not its distributions').
struct Draw {
  std::mt19937_64 engine;
  explicit Draw(std::uint64_t seed) : engine(seed) {}
  std::uint64_t below(std::uint64_t n) { return engine() % n; }
};

std::int64_t lcm(std::int64_t a, std::int64_t b) { return a / std::gcd(a, b) * b; }

// A random set of n tasks with integer periods from a set whose hyperperiod is at most 120, costs of
// at least 1, implicit deadlines; the utilisation is whatever it comes to.
std::vector<PeriodicTask> random_set(Draw& d, std::size_t n, double max_share = 0.5) {
  static const int periods[] = {2, 3, 4, 5, 6, 8, 10, 12, 15, 20, 24, 30, 40, 60, 120};
  std::vector<PeriodicTask> set;
  for (std::size_t i = 0; i < n; i++) {
    const int t = periods[d.below(std::size(periods))];
    const int c_max = std::max(1, static_cast<int>(std::floor(t * max_share)));
    const int c = 1 + static_cast<int>(d.below(static_cast<std::uint64_t>(c_max)));
    PeriodicTask task;
    task.name = "t" + std::to_string(i);
    task.cost = c;
    task.period = t;
    set.push_back(task);
  }
  return set;
}

double hyperperiod(const std::vector<PeriodicTask>& set) {
  std::int64_t h = 1;
  for (const auto& t : set) h = lcm(h, static_cast<std::int64_t>(t.period));
  return static_cast<double>(h);
}

}  // namespace

int main() {
  calx::verification::Suite s;

  // ---- the Liu–Layland bound -------------------------------------------------------------------------
  s.check("Liu–Layland bound for one task is 1", liu_layland_bound(1) - 1, 1e-15);
  s.check("for two, 2(√2 − 1)", liu_layland_bound(2) - 2 * (std::sqrt(2.0) - 1), 1e-15);
  s.check("as n grows, ln 2 (n = 10⁶)", std::fabs(liu_layland_bound(1000000) - std::log(2.0)), 1e-6);

  // ---- response-time analysis against the simulation ----------------------------------------------
  {
    Draw d(20261028);
    std::size_t sets = 0, accepted = 0, rejected = 0, disagree = 0, worse = 0;
    double worst_gap = 0;
    for (int trial = 0; trial < 400; trial++) {
      auto set = random_set(d, 2 + d.below(5), 0.45);
      // Rate-monotonic priorities, written into the tasks.
      const auto order = rate_monotonic(set);
      std::vector<PeriodicTask> by_priority;
      for (std::size_t k = 0; k < order.size(); k++) { by_priority.push_back(set[order[k]]); by_priority.back().priority = static_cast<int>(k); }
      const auto r = response_times(by_priority);
      const double h = hyperperiod(by_priority);
      SimulationOptions o;
      o.policy = Policy::fixed_priority;
      const auto sim = simulate(periodic_jobs(by_priority, h), h + 4 * h, o);
      sets++;
      // The first task the analysis rejects, if any.
      std::size_t first_bad = by_priority.size();
      for (std::size_t i = 0; i < r.size(); i++)
        if (!r[i]) { first_bad = i; break; }
      for (std::size_t i = 0; i < first_bad; i++) {
        double first = -1, most = 0;
        for (const auto& out : sim.outcomes) {
          if (out.job.task != i) continue;
          if (out.job.index == 0) first = out.response();
          most = std::max(most, out.response());
        }
        worst_gap = std::max(worst_gap, std::fabs(first - *r[i]));
        if (most > *r[i] + 1e-9) worse++;
      }
      if (first_bad < by_priority.size()) {
        rejected++;
        // Its first job, released with every higher-priority task's, finishes after its deadline.
        for (const auto& out : sim.outcomes)
          if (out.job.task == first_bad && out.job.index == 0 && !(out.finish > out.job.deadline + 1e-9)) disagree++;
      } else {
        accepted++;
        if (sim.misses) disagree++;
      }
    }
    s.require("random sets of 2 to 6 tasks, some accepted and some rejected by the analysis", accepted > 50 && rejected > 50,
              std::to_string(accepted) + " accepted, " + std::to_string(rejected) + " rejected of " + std::to_string(sets));
    s.check("the first job's response equals Rᵢ (the critical instant)", worst_gap, 1e-9);
    s.require("no job's response exceeds Rᵢ", worse == 0, std::to_string(worse) + " exceed it");
    s.require("the analysis's verdict is the simulation's", disagree == 0, std::to_string(disagree) + " disagree");
  }
  // A set worked by hand: (C, T) = (1, 4), (2, 6), (3, 10): R = 1, 3, 10.
  {
    std::vector<PeriodicTask> set{{"a", 1, 4}, {"b", 2, 6}, {"c", 3, 10}};
    const auto r = response_times(set);
    s.require("(1,4), (2,6), (3,10): R = 1, 3, 10", r[0] && r[1] && r[2] && *r[0] == 1 && *r[1] == 3 && *r[2] == 10);
    s.require("its utilisation 0.883 exceeds the Liu–Layland bound 0.780 and the set is still schedulable",
              utilisation(set) > liu_layland_bound(3) && r[2]);
  }

  // ---- EDF ---------------------------------------------------------------------------------------------------
  {
    Draw d(7);
    std::size_t under = 0, over = 0, bad_under = 0, bad_over = 0;
    for (int trial = 0; trial < 300; trial++) {
      auto set = random_set(d, 2 + d.below(5), 0.6);
      const double h = hyperperiod(set), u = utilisation(set);
      // Jobs released within the hyperperiod, followed until they all finish: every deadline is ≤ h.
      const auto sim = simulate(periodic_jobs(set, h), 3 * h, {});
      if (u <= 1 + 1e-12) { under++; if (sim.misses) bad_under++; }
      else { over++; if (!sim.misses) bad_over++; }
    }
    s.require("EDF: implicit deadlines and U ≤ 1 meet every deadline", under > 50 && bad_under == 0, std::to_string(bad_under) + " of " + std::to_string(under));
    s.require("EDF: U > 1 misses within the hyperperiod", over > 50 && bad_over == 0, std::to_string(bad_over) + " of " + std::to_string(over));
  }
  {
    Draw d(11);
    std::size_t feasible = 0, infeasible = 0, disagree = 0;
    for (int trial = 0; trial < 400; trial++) {
      auto set = random_set(d, 2 + d.below(4), 0.5);
      for (auto& t : set) {
        const int c = static_cast<int>(t.cost), p = static_cast<int>(t.period);
        t.deadline = c + static_cast<double>(d.below(static_cast<std::uint64_t>(p - c + 1)));  // C ≤ D ≤ T
      }
      const auto v = edf_feasible(set);
      const double h = hyperperiod(set);
      double longest = 0;
      for (const auto& t : set) longest = std::max(longest, t.relative_deadline());
      const auto sim = simulate(periodic_jobs(set, h), h + longest, {});
      // Jobs released within the hyperperiod, all of them due by h + longest.
      std::size_t missed = 0;
      for (const auto& out : sim.outcomes)
        if (!out.finished() || out.finish > out.job.deadline + 1e-9) missed++;
      if (v.feasible) { feasible++; if (missed) disagree++; }
      else { infeasible++; if (!missed) disagree++; }
    }
    s.require("the processor-demand criterion and the simulation agree on constrained deadlines",
              feasible > 40 && infeasible > 40 && disagree == 0,
              std::to_string(feasible) + " feasible, " + std::to_string(infeasible) + " not, " + std::to_string(disagree) + " disagree");
  }
  {
    Draw d(13);
    std::size_t rm_ok = 0, edf_fails = 0;
    for (int trial = 0; trial < 200; trial++) {
      auto set = random_set(d, 2 + d.below(5), 0.45);
      const auto order = rate_monotonic(set);
      std::vector<PeriodicTask> by_priority;
      for (std::size_t k : order) by_priority.push_back(set[k]);
      const auto r = response_times(by_priority);
      if (std::all_of(r.begin(), r.end(), [](const auto& x) { return x.has_value(); })) {
        rm_ok++;
        const double h = hyperperiod(set);
        if (simulate(periodic_jobs(set, h), 3 * h, {}).misses) edf_fails++;
      }
    }
    s.require("EDF meets every deadline of every set rate-monotonic priorities meet", rm_ok > 40 && edf_fails == 0);
  }

  // ---- the constant bandwidth server ----------------------------------------------------------------
  {
    ConstantBandwidthServer cbs(2, 5);  // Q = 2, T = 5, U = 0.4
    const double d1 = cbs.arrive(0, true);
    const double d2 = cbs.consume(2);
    const double d3 = cbs.arrive(3, true);   // c = 2 < (10 − 3) 0.4 = 2.8: keep (c, d)
    const double d4 = cbs.arrive(9, true);   // c = 2 ≥ (10 − 9) 0.4 = 0.4: d = 14, c = 2
    s.require("CBS rules (Q 2, T 5): arrival at 0 → d 5; budget spent → d 10; arrival at 3 keeps 10; at 9 → 14",
              d1 == 5 && d2 == 10 && d3 == 10 && d4 == 14 && cbs.remaining() == 2);
  }
  {
    // Hard tasks (1, 4), (2, 8) and (1, 4), U = 0.75; a server (1, 4), U = 0.25, so that the
    // processor is full; jobs every unit that each need 2: the stream asks for eight times the
    // server's share. (With idle time left, the soft CBS would also take the idle time: its bound
    // is on what it takes from the others, not on what it takes.)
    std::vector<PeriodicTask> hard{{"a", 1, 4}, {"b", 2, 8}, {"c", 1, 4}};
    const double h = 400;
    auto jobs = periodic_jobs(hard, h);
    std::vector<Job> stream;
    for (int k = 0; k < 400; k++) stream.push_back({9, static_cast<std::size_t>(k), static_cast<double>(k), 2.0, k + 1.0, 0, 0, 0});
    jobs.insert(jobs.end(), stream.begin(), stream.end());
    SimulationOptions o;
    o.servers.emplace_back(1, 4);
    const auto sim = simulate(jobs, h, o);
    std::size_t hard_missed = 0;
    for (const auto& out : sim.outcomes)
      if (out.job.server == no_server && out.missed(h)) hard_missed++;
    s.require("CBS isolation: the hard tasks meet every deadline beside a greedy stream", hard_missed == 0, std::to_string(hard_missed) + " missed");
    s.check("the server receives no more than U h + Q of the processor (100 + 1 in 400)", std::max(0.0, sim.served(0) - (0.25 * h + 1)), 1e-9);
    s.check("and does receive its share (at least U h − Q)", std::max(0.0, (0.25 * h - 1) - sim.served(0)), 1e-9);
    // The same stream without a server: its jobs take EDF deadlines of their own and crowd out the tasks.
    auto bare = jobs;
    for (auto& j : bare) j.server = no_server;
    const auto unserved = simulate(bare, h, {});
    std::size_t crowded = 0;
    for (const auto& out : unserved.outcomes)
      if (out.job.task != 9 && out.missed(h)) crowded++;
    s.require("without the server the same stream makes the hard tasks miss", crowded > 0, std::to_string(crowded) + " missed");
  }

  // ---- aging -----------------------------------------------------------------------------------------
  {
    // A stream that fills the processor (C = T = 1, priority 0) and one job of priority 5 at 0.
    std::vector<Job> jobs;
    for (int k = 0; k < 100; k++) jobs.push_back({0, static_cast<std::size_t>(k), static_cast<double>(k), 1.0, k + 1.0, 0, 0, no_server});
    jobs.push_back({1, 0, 0.0, 1.0, 1000.0, 5, 0, no_server});
    SimulationOptions o;
    o.policy = Policy::fixed_priority;
    const auto starved = simulate(jobs, 100, o);
    s.require("fixed priority without aging: the low job never runs beside a full stream", starved.outcomes.back().start < 0);
    o.aging_step = 1;
    const auto aged = simulate(jobs, 100, o);
    s.check("with aging (step 1): it runs when its level 5 has risen to the stream's 0, at t = 5", std::fabs(aged.outcomes.back().start - 5), 1e-9);
  }

  // ---- non-preemptive and global EDF --------------------------------------------------------------
  {
    Draw d(17);
    auto set = random_set(d, 5, 0.5);
    const double h = hyperperiod(set);
    SimulationOptions o;
    o.preemptive = false;
    const auto sim = simulate(periodic_jobs(set, h), 4 * h, o);
    std::vector<int> slices(sim.outcomes.size(), 0);
    for (const auto& sl : sim.slices) slices[sl.job]++;
    s.require("non-preemptive: every job that ran ran in one slice", std::all_of(slices.begin(), slices.end(), [](int k) { return k <= 1; }));
  }
  {
    Draw d(19);
    std::vector<Job> jobs;
    for (std::size_t k = 0; k < 60; k++) {
      const double r = static_cast<double>(d.below(40)), c = 1 + static_cast<double>(d.below(5));
      jobs.push_back({k, 0, r, c, r + c + static_cast<double>(d.below(10)), 0, 0, no_server});
    }
    SimulationOptions o;
    o.processors = 2;
    const auto sim = simulate(jobs, 1000, o);
    std::vector<double> times;
    for (const auto& sl : sim.slices) { times.push_back(sl.start); times.push_back(sl.end); }
    for (const auto& j : jobs) times.push_back(j.release);
    std::sort(times.begin(), times.end());
    times.erase(std::unique(times.begin(), times.end()), times.end());
    std::size_t too_many = 0, idle_with_work = 0;
    for (std::size_t k = 0; k + 1 < times.size(); k++) {
      const double mid = 0.5 * (times[k] + times[k + 1]);
      std::size_t on = 0;
      for (const auto& sl : sim.slices) on += sl.start < mid && mid < sl.end;
      std::size_t waiting = 0;
      for (std::size_t j = 0; j < jobs.size(); j++) {
        const auto& out = sim.outcomes[j];
        bool running = false;
        for (const auto& sl : sim.slices) running = running || (sl.job == j && sl.start < mid && mid < sl.end);
        if (jobs[j].release < mid && (!out.finished() || out.finish > mid) && !running) waiting++;
      }
      if (on > 2) too_many++;
      if (waiting > 0 && on < 2) idle_with_work++;
    }
    std::size_t unfinished = 0;
    for (const auto& out : sim.outcomes) unfinished += !out.finished();
    s.require("global EDF on two processors: never more than two at once", too_many == 0);
    s.require("and never a processor idle while a job waits", idle_with_work == 0, std::to_string(idle_with_work) + " intervals");
    s.require("and every job finished", unfinished == 0);
  }

  // ---- elastic compression -------------------------------------------------------------------------
  {
    // U = 1/4 + 2/5 + 3/6 = 1.15, compressed to 0.9; the third is rigid.
    std::vector<Elastic> set{{1, 4, 20, 1}, {2, 5, 25, 3}, {3, 6, 6, 0}};
    const auto p = compress(set, 0.9);
    double u = 0;
    bool within = true;
    if (p) for (std::size_t i = 0; i < set.size(); i++) { u += set[i].cost / (*p)[i]; within = within && (*p)[i] >= set[i].period - 1e-12 && (*p)[i] <= std::max(set[i].longest, set[i].period) + 1e-9; }
    s.require("elastic compression: a solution", p.has_value());
    s.check("its utilisation is the desired 0.9", std::fabs(u - 0.9), 1e-12);
    s.require("within the bounds, the rigid task unchanged", within && p && (*p)[2] == 6);
    // Each elastic task gave up utilisation in proportion to its coefficient: 0.25 → 0.1875, 0.4 → 0.2125.
    s.check("in proportion to the coefficients (1 : 3 of the excess 0.25)", std::fabs(1 / (*p)[0] - 0.1875) + std::fabs(2 / (*p)[1] - 0.2125), 1e-12);
    // A first task that may stretch only to 4.5 is clamped there and the second takes the rest.
    std::vector<Elastic> clamp{{1, 4, 4.5, 1}, {2, 5, 25, 1}, {3, 6, 6, 0}};
    const auto q = compress(clamp, 0.9);
    s.require("a task clamped at its longest period, the others compressed further", q && std::fabs((*q)[0] - 4.5) < 1e-12 && std::fabs(1 / 4.5 + 2 / (*q)[1] + 0.5 - 0.9) < 1e-12);
    s.require("an unreachable utilisation is refused", !compress(set, 0.55).has_value());
  }

  // ---- dispatch on threads ----------------------------------------------------------------------------
  {
    ThreadPool one(1);
    std::vector<LiveJob> jobs;
    Draw d(23);
    for (int k = 0; k < 30; k++) jobs.push_back({"j" + std::to_string(k), 1.0 + static_cast<double>(d.below(1000)), 0, 0, no_server, nullptr});
    const auto run = dispatch(one, jobs, Policy::edf);
    bool ordered = run.records.size() == jobs.size();
    for (std::size_t k = 1; k < run.records.size(); k++) ordered = ordered && jobs[run.records[k - 1].job].deadline <= jobs[run.records[k].job].deadline;
    s.require("one thread runs the jobs in EDF's order", ordered);
  }
  {
    ThreadPool four(4);
    std::vector<std::atomic<int>> runs(64);
    std::vector<LiveJob> jobs;
    for (std::size_t k = 0; k < 64; k++)
      jobs.push_back({"j", static_cast<double>(k % 7), static_cast<int>(k % 3), 0, k % 4 == 0 ? 0 : no_server, [&runs, k](std::size_t) {
                        runs[k]++;
                        volatile double x = 0;
                        for (int i = 0; i < 20000; i++) x = x + std::sqrt(static_cast<double>(i));
                      }});
    const auto run = dispatch(four, jobs, Policy::edf);
    bool once = true;
    for (auto& r : runs) once = once && r.load() == 1;
    s.require("four threads run each job once", once && run.records.size() == 64);
    bool apart = true;
    for (const auto& a : run.records)
      for (const auto& b : run.records)
        if (&a != &b && jobs[a.job].group == 0 && jobs[b.job].group == 0 && a.start < b.end && b.start < a.end) apart = false;
    s.require("the jobs of one group never overlap", apart);
    std::vector<LiveJob> past{{"overdue", 0.0, 0, 0, no_server, [](std::size_t) { volatile int x = 0; for (int i = 0; i < 1000; i++) x = x + i; }}};
    s.require("a job due before it starts is counted late", dispatch(four, past).late == 1);
  }

  return s.done();
}
