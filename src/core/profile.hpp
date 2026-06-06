#pragma once
// Lightweight, OPT-IN profiling counters/timers. Compiled to NOTHING unless
// -DGH_PROFILE is set, so normal builds and the full test suite are byte-for-byte
// unaffected. Used to break down the GeneralRep large-irrep HMC hot path
// (fund_alg / A-build / dmat_expi / matvec) and the three MD force sweeps.
//
// Usage:
//   GH_PROF_COUNT(fast_D);              // bump a call counter
//   { GH_PROF_SCOPE(dmat_expi); ... }   // accumulate wall time of a scope
//   gh::prof::report();                 // print the table (call at driver exit)
//
// Timers are thread-local accumulators summed at report time; counters are atomic.
#ifdef GH_PROFILE
#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <mutex>
#include <string>
#include <vector>

namespace gh {
namespace prof {

constexpr int kMaxSlots = 64;

// Per-thread accumulators: ZERO shared state in the hot path (no atomics, no CAS),
// so the timers themselves do NOT serialize a heavily-threaded sweep. Each thread
// owns a private array; report() registers it once (under a lock) and sums at exit.
struct ThreadAcc {
  std::uint64_t calls[kMaxSlots] = {};
  double        secs [kMaxSlots] = {};
};

inline std::mutex& reg_mutex() { static std::mutex m; return m; }
inline std::vector<ThreadAcc*>& all_accs() { static std::vector<ThreadAcc*> v; return v; }
inline ThreadAcc& tls_acc() {
  static thread_local ThreadAcc acc;
  static thread_local bool registered = false;
  if (!registered) {
    std::lock_guard<std::mutex> lk(reg_mutex());
    all_accs().push_back(&acc);
    registered = true;
  }
  return acc;
}

inline std::vector<std::string>& names() {
  static std::vector<std::string> n(kMaxSlots); return n; }

// Register a named slot once; return its stable index. Called only at static-init
// of each instrumented call site (amortized), never in the steady-state hot loop.
inline int slot_id(const char* name) {
  static std::atomic<int> next{0};
  std::lock_guard<std::mutex> lk(reg_mutex());
  auto& nm = names();
  for (int i = 0; i < next.load(); ++i)
    if (nm[i] == name) return i;
  const int id = next.fetch_add(1);
  nm[id] = name;
  return id;
}

struct ScopeTimer {
  int id;
  ThreadAcc* acc;
  std::chrono::steady_clock::time_point t0;
  explicit ScopeTimer(int i) : id(i), acc(&tls_acc()),
                               t0(std::chrono::steady_clock::now()) {}
  ~ScopeTimer() {
    const double dt = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t0).count();
    acc->calls[id] += 1;       // thread-private: no contention
    acc->secs[id]  += dt;
  }
};

inline void report() {
  std::lock_guard<std::mutex> lk(reg_mutex());
  std::uint64_t calls[kMaxSlots] = {};
  double        secs [kMaxSlots] = {};
  for (auto* a : all_accs())
    for (int i = 0; i < kMaxSlots; ++i) { calls[i] += a->calls[i]; secs[i] += a->secs[i]; }
  auto& nm = names();
  std::printf("\n==== GH_PROFILE breakdown (timers thread-private; total_s summed over threads) ====\n");
  std::printf("%-26s %14s %14s %14s\n", "slot", "calls", "total_s", "us/call");
  for (int i = 0; i < kMaxSlots; ++i) {
    if (nm[i].empty()) continue;
    const std::uint64_t c = calls[i];
    const double tot = secs[i];
    const double uspc = c ? (tot / c) * 1e6 : 0.0;
    std::printf("%-26s %14llu %14.4f %14.4f\n",
                nm[i].c_str(), (unsigned long long)c, tot, uspc);
  }
  std::printf("====================================================================================\n");
}

}  // namespace prof
}  // namespace gh

// id is cached per call-site in a static, so slot_id's lock is paid once.
#define GH_PROF_SCOPE(NAME)                                            \
  static int _gh_prof_id_##NAME = ::gh::prof::slot_id(#NAME);          \
  ::gh::prof::ScopeTimer _gh_prof_scope_##NAME(_gh_prof_id_##NAME)
#define GH_PROF_COUNT(NAME)                                            \
  do { static int _gh_prof_cid_##NAME = ::gh::prof::slot_id(#NAME);    \
       ::gh::prof::tls_acc().calls[_gh_prof_cid_##NAME] += 1; } while (0)
#define GH_PROF_REPORT() ::gh::prof::report()

#else  // !GH_PROFILE  -- everything compiles away

#define GH_PROF_SCOPE(NAME)   do {} while (0)
#define GH_PROF_COUNT(NAME)   do {} while (0)
#define GH_PROF_REPORT()      do {} while (0)

#endif
