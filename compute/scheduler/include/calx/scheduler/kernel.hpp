// The kernel beneath the dispatcher: what real-time scheduling the operating system offers this
// process, and its use where it is permitted (decision 0010).
//
// realtime.hpp's dispatcher orders jobs in user space and runs them on a pool of threads; the kernel
// decides when those threads run. On a general-purpose kernel under time sharing (Linux's
// SCHED_OTHER, Windows's normal priority) a ready worker can wait behind any other process, so
// lateness is measured, not bounded. Linux offers two real-time classes (sched(7)): SCHED_FIFO,
// fixed priorities 1–99 above every time-shared thread, run until they block or yield; and
// SCHED_DEADLINE, EDF with constant bandwidth servers in the kernel itself (a runtime, a deadline and
// a period per thread, admitted by the kernel's own utilisation test). Using either needs the
// privilege CAP_SYS_NICE or, for SCHED_FIFO, a non-zero RLIMIT_RTPRIO; a fully preemptible kernel
// (PREEMPT_RT, in the mainline since Linux 6.12 as a build option) bounds the kernel's own latency.
//
// probe_kernel() says which of these the process has, by asking the kernel and by trying each class
// on a short-lived thread of its own and leaving it at once, so that the caller's threads are never
// changed. raise_to_fifo() and lower_to_time_sharing() change the calling thread; with a ThreadPool,
// pool.run() applies them to every worker. On other systems the probe reports time sharing only and
// the changes refuse.
//
// Caveat. A SCHED_FIFO thread that never blocks starves the time-shared threads of its core: raise
// only workers that wait when idle (ThreadPool's do, after a bounded spin), and prefer
// SCHED_DEADLINE's reservations where the kernel permits them.

#pragma once

#include <cerrno>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <string>
#include <thread>

#if defined(__linux__)
#include <pthread.h>
#include <sched.h>
#include <sys/resource.h>
#include <sys/syscall.h>
#include <sys/utsname.h>
#include <unistd.h>
#endif

namespace calx::scheduler {

struct KernelReport {
  std::string system;              // "Linux", "Windows", …
  std::string release, version;    // uname's
  bool linux_kernel = false;
  bool wsl = false;                // a Linux kernel under Windows's subsystem
  bool preempt_rt = false;         // fully preemptible: /sys/kernel/realtime is 1, or PREEMPT_RT in the version
  bool fifo_permitted = false;     // SCHED_FIFO could be entered and left
  bool deadline_permitted = false; // SCHED_DEADLINE could be entered and left
  int fifo_min = 0, fifo_max = 0;  // its priorities
  long rtprio_limit = 0;           // RLIMIT_RTPRIO's soft limit (−1: unlimited)
  std::string refusal;             // why the classes were refused, when they were
  std::string about;               // one sentence: what the dispatcher's timing rests on
};

namespace kernel_detail {
#if defined(__linux__)
// struct sched_attr of sched_setattr(2), which glibc declared only lately.
struct SchedAttr {
  std::uint32_t size;
  std::uint32_t sched_policy;
  std::uint64_t sched_flags;
  std::int32_t sched_nice;
  std::uint32_t sched_priority;
  std::uint64_t sched_runtime, sched_deadline, sched_period;  // ns
};
constexpr std::uint32_t policy_deadline = 6;  // SCHED_DEADLINE

inline int set_attr(const SchedAttr& a) {
#if defined(SYS_sched_setattr)
  return int(syscall(SYS_sched_setattr, 0, &a, 0));
#else
  (void)a;
  errno = ENOSYS;
  return -1;
#endif
}

// Runs f on a new thread and returns its result: the probe never changes the caller's thread.
template <class F>
auto on_own_thread(F f) {
  decltype(f()) r{};
  std::thread t([&] { r = f(); });
  t.join();
  return r;
}
#endif
}  // namespace kernel_detail

inline KernelReport probe_kernel() {
  KernelReport k;
#if defined(__linux__)
  namespace d = kernel_detail;
  k.linux_kernel = true;
  utsname u{};
  if (uname(&u) == 0) {
    k.system = u.sysname;
    k.release = u.release;
    k.version = u.version;
  }
  k.wsl = k.release.find("microsoft") != std::string::npos || k.release.find("WSL") != std::string::npos;
  {
    std::ifstream rt("/sys/kernel/realtime");
    int flag = 0;
    if (rt >> flag && flag == 1) k.preempt_rt = true;
  }
  if (k.version.find("PREEMPT_RT") != std::string::npos) k.preempt_rt = true;
  k.fifo_min = sched_get_priority_min(SCHED_FIFO);
  k.fifo_max = sched_get_priority_max(SCHED_FIFO);
  rlimit lim{};
  if (getrlimit(RLIMIT_RTPRIO, &lim) == 0) k.rtprio_limit = lim.rlim_cur == RLIM_INFINITY ? -1 : long(lim.rlim_cur);
  // SCHED_FIFO at its lowest priority, then back.
  const int fifo_error = d::on_own_thread([&] {
    sched_param p{};
    p.sched_priority = k.fifo_min;
    const int e = pthread_setschedparam(pthread_self(), SCHED_FIFO, &p);
    if (e == 0) {
      sched_param o{};
      pthread_setschedparam(pthread_self(), SCHED_OTHER, &o);
    }
    return e;
  });
  k.fifo_permitted = fifo_error == 0;
  // SCHED_DEADLINE with a small reservation (1 ms every 10 ms), then back.
  const int deadline_error = d::on_own_thread([] {
    d::SchedAttr a{};
    a.size = sizeof a;
    a.sched_policy = d::policy_deadline;
    a.sched_runtime = 1'000'000;
    a.sched_deadline = a.sched_period = 10'000'000;
    if (d::set_attr(a) == 0) {
      d::SchedAttr o{};
      o.size = sizeof o;
      o.sched_policy = SCHED_OTHER;
      d::set_attr(o);
      return 0;
    }
    return errno;
  });
  k.deadline_permitted = deadline_error == 0;
  const auto why = [&](int e) {
    if (e == EPERM) return std::string("not permitted: it needs CAP_SYS_NICE") + (k.rtprio_limit == 0 ? ", and RLIMIT_RTPRIO is 0" : "");
    return std::string(std::strerror(e));
  };
  if (!k.fifo_permitted) k.refusal = "SCHED_FIFO " + why(fifo_error);
  if (!k.deadline_permitted) k.refusal += std::string(k.refusal.empty() ? "" : "; ") + "SCHED_DEADLINE " + why(deadline_error);
  if (k.deadline_permitted && k.fifo_permitted)
    k.about = "the kernel offers SCHED_DEADLINE and SCHED_FIFO to this process";
  else if (k.fifo_permitted)
    k.about = "the kernel offers SCHED_FIFO to this process (" + k.refusal + ")";
  else
    k.about = "the kernel's real-time classes are refused (" + k.refusal + "): the dispatcher's threads are time-shared, so lateness is measured, not bounded";
  if (!k.preempt_rt && k.fifo_permitted) k.about += "; the kernel is not fully preemptible (no PREEMPT_RT)";
#else
#if defined(_WIN32)
  k.system = "Windows";
#else
  k.system = "unknown";
#endif
  k.refusal = "not offered on this system by calx";
  k.about = "no real-time class is used on this system: the dispatcher's threads are time-shared, so lateness is measured, not bounded";
#endif
  return k;
}

// The calling thread to SCHED_FIFO at `priority`; false, with nothing changed, if refused.
inline bool raise_to_fifo(int priority) {
#if defined(__linux__)
  sched_param p{};
  p.sched_priority = priority;
  return pthread_setschedparam(pthread_self(), SCHED_FIFO, &p) == 0;
#else
  (void)priority;
  return false;
#endif
}

// The calling thread back to time sharing.
inline bool lower_to_time_sharing() {
#if defined(__linux__)
  sched_param p{};
  return pthread_setschedparam(pthread_self(), SCHED_OTHER, &p) == 0;
#else
  return true;
#endif
}

// The calling thread's class: "SCHED_FIFO 10", "SCHED_DEADLINE", "SCHED_OTHER", or "time sharing".
inline std::string current_class() {
#if defined(__linux__)
  int policy = 0;
  sched_param p{};
  if (pthread_getschedparam(pthread_self(), &policy, &p) != 0) return "unknown";
  if (policy == SCHED_FIFO) return "SCHED_FIFO " + std::to_string(p.sched_priority);
  if (policy == SCHED_RR) return "SCHED_RR " + std::to_string(p.sched_priority);
  if (policy == int(kernel_detail::policy_deadline)) return "SCHED_DEADLINE";
  return "SCHED_OTHER";
#else
  return "time sharing";
#endif
}

}  // namespace calx::scheduler
