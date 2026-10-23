// Checks of numerical results against references, reported line by line.
//
// A check compares a value (usually an error) with a limit; an assertion states a
// condition; a skipped check says why. The report is the one every check of The Science
// Infrastructures prints, whatever its language, so that one runner adds them up:
//
//   ok    name: 1.030e-2 unit (limit 0.02 unit)
//   FAIL  name: detail
//   skip  what
//
//   41 of 41 checks passed.
//
// The helpers measure errors as the checks need them: relative error, the largest
// difference of two sequences, and the order of convergence from errors at two steps.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <span>
#include <string>
#include <string_view>

namespace calx::verification {

// A number as JavaScript's toExponential(3) writes it: 1.030e-2, 5.000e+0.
inline std::string exponential(double v) {
  if (std::isnan(v)) return "NaN";
  if (std::isinf(v)) return v > 0 ? "Infinity" : "-Infinity";
  char buf[48];
  std::snprintf(buf, sizeof buf, "%.3e", v);
  std::string s(buf);
  const auto e = s.find('e');
  if (e == std::string::npos) return s;
  std::string mant = s.substr(0, e), exp = s.substr(e + 1);
  const char sign = exp[0];
  std::size_t i = 1;
  while (i + 1 < exp.size() && exp[i] == '0') i++;
  return mant + "e" + sign + exp.substr(i);
}

// The shortest form of a limit, as JavaScript prints a number: 0.02, 1e-12.
inline std::string shortest(double v) {
  char buf[48];
  for (int p = 1; p <= 17; p++) {
    std::snprintf(buf, sizeof buf, "%.*g", p, v);
    if (std::strtod(buf, nullptr) == v) break;
  }
  return buf;
}

class Suite {
 public:
  // value: an error or a deviation; it passes when finite and |value| <= limit.
  bool check(std::string_view name, double value, double limit, std::string_view unit = "") {
    count_++;
    const bool ok = std::isfinite(value) && std::fabs(value) <= limit;
    if (!ok) failures_++;
    const std::string u = unit.empty() ? "" : " " + std::string(unit);
    std::cout << (ok ? "ok  " : "FAIL") << "  " << name << ": " << exponential(value) << u << " (limit " << shortest(limit) << u << ")\n";
    return ok;
  }
  bool require(std::string_view name, bool condition, std::string_view detail = "") {
    count_++;
    if (!condition) failures_++;
    std::cout << (condition ? "ok  " : "FAIL") << "  " << name << (detail.empty() ? "" : ": ") << detail << "\n";
    return condition;
  }
  void skip(std::string_view what) { std::cout << "skip  " << what << "\n"; }
  // Prints the count; the exit status of the program.
  int done() {
    std::cout << "\n" << (count_ - failures_) << " of " << count_ << " checks passed.\n";
    std::cout.flush();
    return failures_ ? 1 : 0;
  }
  int failures() const { return failures_; }

 private:
  int count_ = 0, failures_ = 0;
};

// |a - b| / |b|, or |a - b| when b is 0.
inline double relative(double a, double b) { return b == 0 ? std::fabs(a - b) : std::fabs(a - b) / std::fabs(b); }

inline double max_difference(std::span<const double> a, std::span<const double> b) {
  double m = 0;
  for (std::size_t i = 0; i < std::min(a.size(), b.size()); i++) m = std::max(m, std::fabs(a[i] - b[i]));
  return a.size() == b.size() ? m : INFINITY;
}

// The observed order of convergence from errors e1 at step h and e2 at step h / r.
inline double order(double e1, double e2, double r = 2) { return std::log(e1 / e2) / std::log(r); }

}  // namespace calx::verification
