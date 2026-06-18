// Gate for the value-superop hot-path optimization (docs/su2_su3_llr_speedup_plan.md item 2):
// CasimirChannels::value_cached == value (uncached) to projector precision, and the
// MultiInvariantPotential value-cache self-check is ACTIVE. Also reports the speedup.
#include "action/scalar_invariants.hpp"
#include "rep/rep_general.hpp"
#include <cstdio>
#include <random>
#include <chrono>
#include <vector>

using namespace gh;
static int pass = 0, fail = 0;
static void check(bool c, const char* m) { if (c) ++pass; else { ++fail; std::printf("  FAIL: %s\n", m); } }

template <int N>
static void run(const std::vector<int>& rows, const char* tag) {
  GeneralRep<N> rep(rows);
  CasimirChannels<N> ch(rep);
  const int nc = ch.n_channels();
  // BT-style couplings: f_c = c (1..nc) just to exercise all channels with mu2.
  std::vector<Real> f(nc); for (int c = 0; c < nc; ++c) f[c] = 0.3 + 0.2 * c;
  const Real mu2 = 0.113;
  auto Q = ch.build_value_superop(f);

  std::mt19937_64 rng(12345); std::normal_distribution<Real> gd(0.0, 1.0);
  Real maxrel = 0.0;
  for (int t = 0; t < 200; ++t) {
    DVec phi(ch.d); for (int a = 0; a < ch.d; ++a) phi(a) = Complex(gd(rng), gd(rng));
    const Real vu = ch.value(phi, f, mu2), vc = ch.value_cached(phi, Q, mu2);
    const Real den = std::max<Real>(1e-30, std::fabs(vu));
    maxrel = std::max(maxrel, std::fabs(vc - vu) / den);
  }
  std::printf("  %s (d=%d, %d ch): max rel err value_cached vs value = %.2e\n", tag, ch.d, nc, maxrel);
  check(maxrel < 1e-6, tag);

  // the potential's own self-check must have ENABLED the cache
  MultiInvariantPotential<N> pot(ch, f, mu2);
  check(pot.value_cache_ok, "MultiInvariantPotential value cache enabled");

  // timing: cached vs uncached
  std::vector<DVec> phis; for (int t = 0; t < 2000; ++t) { DVec p(ch.d); for (int a = 0; a < ch.d; ++a) p(a) = Complex(gd(rng), gd(rng)); phis.push_back(p); }
  volatile Real sink = 0;
  auto t0 = std::chrono::high_resolution_clock::now();
  for (auto& p : phis) sink += ch.value(p, f, mu2);
  auto t1 = std::chrono::high_resolution_clock::now();
  for (auto& p : phis) sink += ch.value_cached(p, Q, mu2);
  auto t2 = std::chrono::high_resolution_clock::now();
  const double tu = std::chrono::duration<double>(t1 - t0).count();
  const double tc = std::chrono::duration<double>(t2 - t1).count();
  std::printf("    value() %.3fus  value_cached() %.3fus  -> %.1fx\n",
              1e6 * tu / phis.size(), 1e6 * tc / phis.size(), tu / tc);
  (void)sink;
}

int main() {
  std::printf("-- value-superop correctness + speedup --\n");
  run<2>({6}, "SU(2) spin-3 (2T)");
  run<2>({8}, "SU(2) spin-4 (2O)");
  run<2>({12}, "SU(2) spin-6 (2I)");
  std::printf("[test_valuecache] %d passed, %d failed\n", pass, fail);
  return fail ? 1 : 0;
}
