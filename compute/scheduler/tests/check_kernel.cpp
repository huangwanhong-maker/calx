// Checks of the kernel probe (kernel.hpp): it leaves the caller's scheduling class as it was; what it
// reports agrees with what the kernel then does (raising a thread to SCHED_FIFO succeeds exactly when
// the probe found it permitted, and a refusal changes nothing); a thread raised can be lowered again;
// the report says what the dispatcher's timing rests on. The checks hold on any kernel, permitted or
// not: they test the agreement, not the privilege.

#include <calx/scheduler/kernel.hpp>
#include <calx/verification/check.hpp>

#include <string>
#include <thread>

using namespace calx::scheduler;

int main() {
  calx::verification::Suite s;
  const std::string before = current_class();
  const KernelReport k = probe_kernel();
  s.require("kernel: the probe leaves the caller's class as it was", current_class() == before, before + " → " + current_class());
  s.require("kernel: the report says what the timing rests on", !k.about.empty(), k.about);
#if defined(__linux__)
  s.require("kernel: on Linux, SCHED_FIFO's priorities are an interval within 1–99",
            k.fifo_min >= 1 && k.fifo_max <= 99 && k.fifo_min <= k.fifo_max, std::to_string(k.fifo_min) + "–" + std::to_string(k.fifo_max));
  bool raised = false, lowered = false;
  std::string during, after;
  std::thread t([&] {
    raised = raise_to_fifo(k.fifo_min);
    during = current_class();
    lowered = lower_to_time_sharing();
    after = current_class();
  });
  t.join();
  s.require("kernel: raising a thread to SCHED_FIFO succeeds exactly when the probe found it permitted", raised == k.fifo_permitted,
            std::string(raised ? "raised" : "refused") + "; " + k.about);
  s.require("kernel: a raised thread runs under SCHED_FIFO, a refused one is unchanged",
            raised ? during == "SCHED_FIFO " + std::to_string(k.fifo_min) : during == "SCHED_OTHER", during);
  s.require("kernel: and it returns to time sharing", lowered && after == "SCHED_OTHER", after);
  s.require("kernel: a refusal is explained", k.fifo_permitted && k.deadline_permitted ? k.refusal.empty() : !k.refusal.empty(), k.refusal);
#else
  s.require("kernel: elsewhere than Linux nothing is raised", !raise_to_fifo(1) && !k.fifo_permitted);
#endif
  return s.done();
}
