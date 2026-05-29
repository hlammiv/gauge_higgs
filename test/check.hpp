#pragma once
// Minimal test harness: counters + macros. Each test .cpp includes this, runs
// checks, and returns g_fail (0 == all passed).
#include <cstdio>
#include <cmath>
#include <string>

static int g_pass = 0;
static int g_fail = 0;

#define CHECK(cond, msg)                                                    \
  do {                                                                      \
    if (cond) { ++g_pass; }                                                 \
    else { ++g_fail; std::printf("  FAIL: %s  (%s:%d)\n", msg, __FILE__, __LINE__); } \
  } while (0)

#define CHECK_CLOSE(a, b, tol, msg)                                         \
  do {                                                                      \
    double _d = std::fabs(double(a) - double(b));                           \
    if (_d <= (tol)) { ++g_pass; }                                          \
    else { ++g_fail; std::printf("  FAIL: %s  |%.3e - %.3e| = %.3e > %.1e  (%s:%d)\n", \
                                 msg, double(a), double(b), _d, double(tol), __FILE__, __LINE__); } \
  } while (0)

inline int report(const char* name) {
  std::printf("[%s] %d passed, %d failed\n", name, g_pass, g_fail);
  return g_fail == 0 ? 0 : 1;
}
