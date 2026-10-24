// The target description of the CPU backend: what the host processor can do, and what this
// build was allowed to use (calx architecture §4: a backend declares its capabilities).
//
// Two descriptions are kept apart, because they differ whenever a program built for one
// machine runs on another:
//
//   host()            the processor the program runs on, read once from CPUID (x86-64) or
//                     assumed from the architecture (AArch64 always has NEON); its threads;
//   LANE_BYTES, ISA   the widest vector the compiler was allowed to emit for this build,
//                     from the predefined macros (__AVX512BW__, __AVX2__, __SSE2__, ...).
//
// A kernel written on calx::cpu::Lanes uses the build's width. A program may compare the two
// to say that it was built for less than the host offers (or, run elsewhere, for more).
//
// CPUID leaves (Intel 64 and IA-32 Architectures Software Developer's Manual, vol. 2A, CPUID;
// AMD64 Architecture Programmer's Manual, vol. 3, appendix E):
//   leaf 1:   EDX bit 26 SSE2; ECX bit 19 SSE4.1, 20 SSE4.2, 23 POPCNT, 27 OSXSAVE, 28 AVX;
//   leaf 7:   EBX bit 5 AVX2, 8 BMI2, 16 AVX-512F, 30 AVX-512BW, 31 AVX-512VL;
//   XGETBV(0) bits 1-2 (SSE and AVX state) and 5-7 (AVX-512 state), saved by the OS.
// AVX and AVX-512 count only when the operating system saves their registers.

#pragma once

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <string>
#include <thread>

#if defined(__x86_64__) || defined(__i386__)
#include <cpuid.h>
#endif

namespace calx::cpu {

// The widest vector, in bytes, that this build may emit, and the instruction set it names.
#if defined(__AVX512BW__)
inline constexpr std::size_t LANE_BYTES = 64;
inline constexpr const char* ISA = "avx512bw";
#elif defined(__AVX2__)
inline constexpr std::size_t LANE_BYTES = 32;
inline constexpr const char* ISA = "avx2";
#elif defined(__SSE4_1__)
inline constexpr std::size_t LANE_BYTES = 16;
inline constexpr const char* ISA = "sse4.1";
#elif defined(__SSE2__)
inline constexpr std::size_t LANE_BYTES = 16;
inline constexpr const char* ISA = "sse2";
#elif defined(__ARM_NEON) || defined(__aarch64__)
inline constexpr std::size_t LANE_BYTES = 16;
inline constexpr const char* ISA = "neon";
#else
inline constexpr std::size_t LANE_BYTES = 16;  // GCC emulates 16-byte vectors where there are none
inline constexpr const char* ISA = "generic";
#endif

struct Target {
  std::string architecture;  // "x86-64", "aarch64", or "other"
  std::string vendor;        // "GenuineIntel", "AuthenticAMD", ... (x86 only)
  std::string brand;         // the processor's name, as CPUID gives it (x86 only)
  unsigned threads = 1;      // hardware threads, as the standard library reports them
  bool sse2 = false, sse41 = false, sse42 = false, popcnt = false, avx = false, avx2 = false, bmi2 = false;
  bool avx512f = false, avx512bw = false, avx512vl = false, neon = false;

  // The widest vector the host could run, in bytes.
  std::size_t widest() const { return avx512bw ? 64 : avx2 ? 32 : (sse2 || neon) ? 16 : 0; }
  // The host's best instruction set, named as ISA names the build's.
  std::string best() const {
    if (avx512bw) return "avx512bw";
    if (avx2) return "avx2";
    if (sse41) return "sse4.1";
    if (sse2) return "sse2";
    if (neon) return "neon";
    return "generic";
  }
};

namespace detail {

#if defined(__x86_64__) || defined(__i386__)
inline std::uint64_t xgetbv0() {
  std::uint32_t lo = 0, hi = 0;
  __asm__ volatile("xgetbv" : "=a"(lo), "=d"(hi) : "c"(0));
  return (static_cast<std::uint64_t>(hi) << 32) | lo;
}
#endif

inline Target detect() {
  Target t;
  t.threads = std::max(1u, std::thread::hardware_concurrency());
#if defined(__x86_64__) || defined(__i386__)
  t.architecture = "x86-64";
  unsigned a = 0, b = 0, c = 0, d = 0;
  if (__get_cpuid(0, &a, &b, &c, &d)) {
    char v[13] = {};
    std::memcpy(v, &b, 4);
    std::memcpy(v + 4, &d, 4);
    std::memcpy(v + 8, &c, 4);
    t.vendor = v;
    const unsigned top = a;
    if (__get_cpuid(1, &a, &b, &c, &d)) {
      t.sse2 = d & (1u << 26);
      t.sse41 = c & (1u << 19);
      t.sse42 = c & (1u << 20);
      t.popcnt = c & (1u << 23);
      const bool osxsave = c & (1u << 27);
      const std::uint64_t xcr0 = osxsave ? xgetbv0() : 0;
      const bool ymm = (xcr0 & 0x6) == 0x6, zmm = (xcr0 & 0xE6) == 0xE6;
      t.avx = (c & (1u << 28)) && ymm;
      if (top >= 7 && __get_cpuid_count(7, 0, &a, &b, &c, &d)) {
        t.avx2 = t.avx && (b & (1u << 5));
        t.bmi2 = b & (1u << 8);
        t.avx512f = zmm && (b & (1u << 16));
        t.avx512bw = t.avx512f && (b & (1u << 30));
        t.avx512vl = t.avx512f && (b & (1u << 31));
      }
    }
  }
  if (__get_cpuid(0x80000000, &a, &b, &c, &d) && a >= 0x80000004) {
    char s[49] = {};
    for (unsigned leaf = 0; leaf < 3; leaf++) {
      __get_cpuid(0x80000002 + leaf, &a, &b, &c, &d);
      std::memcpy(s + 16 * leaf, &a, 4);
      std::memcpy(s + 16 * leaf + 4, &b, 4);
      std::memcpy(s + 16 * leaf + 8, &c, 4);
      std::memcpy(s + 16 * leaf + 12, &d, 4);
    }
    std::string brand(s);
    const auto first = brand.find_first_not_of(' '), last = brand.find_last_not_of(' ');
    t.brand = first == std::string::npos ? "" : brand.substr(first, last - first + 1);
  }
#elif defined(__aarch64__)
  t.architecture = "aarch64";
  t.neon = true;
#else
  t.architecture = "other";
#endif
  return t;
}

}  // namespace detail

// The host, detected once. It describes fixed hardware, so a function-local constant serves
// every caller (it is initialised thread-safely, and never changes).
inline const Target& host() {
  static const Target t = detail::detect();
  return t;
}

// One line for a report: "Intel(R) Core(TM) i7-11370H, 8 threads, best avx512bw; built for avx2".
inline std::string describe(const Target& t = host()) {
  std::string s = t.brand.empty() ? t.architecture : t.brand;
  s += ", " + std::to_string(t.threads) + " threads, best " + t.best() + "; built for " + ISA;
  return s;
}

}  // namespace calx::cpu
