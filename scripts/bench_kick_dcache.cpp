// D-cache / fast-path kick micro-benchmark (task: enable per-link D-cache for spin-3/4).
// Build against a scratch src copy that adds the GH_FORCE_FAST knob to enable_fast_if_valid:
//   g++ -std=c++20 -O3 -march=native -funroll-loops -fopenmp -I<src> -o bench_kick scripts/bench_kick_dcache.cpp
//   OMP_NUM_THREADS=1 ./bench_kick <young 6|8|12> <L> <nkick> [kappa]
//   GH_FORCE_FAST=1 forces the GeneralRep fast/cacheable path on (still gated by the 1e-7 self-check).
// MEASURED (single thread, this box): spin-3{6} cache ON = 0.75x (SLOWER), spin-4{8} = 2.1-2.3x,
//   spin-6{12} already on. The fast path's single expm build only beats 3 tensor applies for d>=9.
//
// Micro-benchmark: time the HMC kick() for a given rep, frozen, cache via GH_FORCE_FAST.
// Mirrors the production frozen GaugeHiggsHMC<4,2> setup. Reports per-kick wall time.
//   ./bench_kick <rep-young> <L> <nkick> [kappa]
// GH_FORCE_FAST=1 forces the GeneralRep fast/cacheable path on (prototype knob).
#include "hmc/gauge_higgs_hmc.hpp"
#include "action/scalar_invariants.hpp"
#include "rep/rep_general.hpp"
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <vector>

using namespace gh;
using Clock = std::chrono::steady_clock;

int main(int argc, char** argv) {
  if (argc < 4) { std::fprintf(stderr, "usage: %s <young e.g.6> <L> <nkick> [kappa]\n", argv[0]); return 1; }
  std::vector<int> rows = { std::atoi(argv[1]) };
  const int L = std::atoi(argv[2]);
  const int nkick = std::atoi(argv[3]);
  const Real kappa = (argc > 4) ? std::atof(argv[4]) : 1.0;

  GeneralRep<2> rep(rows);
  std::array<int,4> ext{L,L,L,L};
  GaugeHiggsHMC<4,2> hmc(ext, rep, 12345);
  hmc.beta = 1.5; hmc.kappa = kappa; hmc.tau = 1.0; hmc.nmd = 10;
  hmc.frozen_phi = true;

  CasimirChannels<2> ch(rep);
  std::vector<Real> f(ch.n_channels(), 0.1);  // generic potential, channels active
  MultiInvariantPotential<2> pot(ch, f, 0.113);
  hmc.potential = &pot;

  hmc.U.hot(hmc.rng, 0.8);
  hmc.phi.gaussian(hmc.rng, 12345, rep.real, 0.3);
  hmc.normalize_phi();
  hmc.refresh_momenta();

  std::printf("# rep={%d} d=%d  use_fast=%d  use_link_cache=%d  L=%d kappa=%.3f\n",
              rows[0], rep.d, (int)rep.use_fast, (int)hmc.use_link_cache, L, kappa);

  // warmup
  for (int i = 0; i < 3; ++i) hmc.kick(0.01);
  auto t0 = Clock::now();
  for (int i = 0; i < nkick; ++i) hmc.kick(0.01);
  auto t1 = Clock::now();
  double sec = std::chrono::duration<double>(t1 - t0).count();
  std::printf("kicks=%d  total=%.4f s  per_kick=%.4f ms\n", nkick, sec, 1e3 * sec / nkick);
  // also: full traj timing (kick cost dominates; 2*nmd+1 kicks for Omelyan)
  return 0;
}
