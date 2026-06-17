// MD-force profiling driver for GaugeHiggsHMC<4,2> (frozen SU(2)+Higgs).
// Dumps the GH_PROFILE per-scope breakdown (kick -> gauge/matter/scalar -> rotate/fast_D/...).
// Vary ONE knob at a time: rep dim (rows), L (volume), nmd (per-MD-step force cost), frozen.
//
// Build:
//   g++ -std=c++20 -O3 -march=native -funroll-loops -fopenmp -Isrc \
//       -DGH_PROFILE -DNDIM=4 -DNCOL=2 -o /tmp/gh_md_profile scripts/gh_md_profile.cpp
// Run (single-thread for clean per-component attribution):
//   OMP_NUM_THREADS=1 /tmp/gh_md_profile <rows> <L> <ntraj> <nmd> [frozen=1] [seed=1]
// rep rows (SU(2) spin-j = Young {2j}): "1"=fund d=2, "6"=spin-3 d=7 (BT 2T),
//   "8"=spin-4 d=9 (BO 2O), "12"=spin-6 d=13 (BI 2I).
#include "hmc/gauge_higgs_hmc.hpp"
#include "action/scalar_invariants.hpp"
#include "rep/rep_general.hpp"
#include "core/profile.hpp"
#include <cstdio>
#include <cstdlib>
#include <string>
#include <chrono>
#include <memory>
using namespace gh;

int main(int argc, char** argv) {
  if (argc < 5) { std::fprintf(stderr,"usage: %s <rows> <L> <ntraj> <nmd> [frozen=1] [seed=1]\n", argv[0]); return 1; }
  std::vector<int> rows{ std::atoi(argv[1]) };
  const int L = std::atoi(argv[2]);
  const int ntraj = std::atoi(argv[3]);
  const int nmd = std::atoi(argv[4]);
  const bool frozen = (argc>5) ? std::atoi(argv[5]) : 1;
  const std::uint64_t seed = (argc>6) ? (std::uint64_t)std::atoll(argv[6]) : 1ULL;

  auto rep = std::make_unique<GeneralRep<2>>(rows, GeneralRep<2>::RealType::Complex);
  CasimirChannels<2> ch(*rep);
  std::vector<Real> f(ch.n_channels(), 1.0);
  MultiInvariantPotential<2> pot(ch, f, 0.113);

  std::array<int,4> ext{}; for(int m=0;m<4;++m) ext[m]=L;
  GaugeHiggsHMC<4,2> hmc(ext, *rep, seed);
  hmc.beta=1.5; hmc.kappa=1.0; hmc.tau=1.0; hmc.nmd=nmd; hmc.potential=&pot;
  hmc.frozen_phi=frozen;
  hmc.U.hot(hmc.rng,0.8); hmc.phi.gaussian(hmc.rng,12345,rep->real,0.3);
  if (frozen) hmc.normalize_phi();

  std::fprintf(stderr,"# rep d=%d use_fast=%d link_cacheable=%d n_channels=%d L=%d nmd=%d frozen=%d\n",
    rep->d, (int)((GeneralRep<2>*)rep.get())->use_fast, (int)rep->link_cacheable(), ch.n_channels(), L, nmd, frozen);

  // a few warmup/therm then time the measured block
  for (int t=0;t<3;++t) hmc.trajectory();
  auto t0 = std::chrono::steady_clock::now();
  for (int t=0;t<ntraj;++t) hmc.trajectory();
  auto t1 = std::chrono::steady_clock::now();
  double sec = std::chrono::duration<double>(t1-t0).count();
  std::fprintf(stderr,"# TIMING: %d trajs in %.3f s -> %.4f s/traj  accept=%.3f\n",
    ntraj, sec, sec/ntraj, hmc.acceptance());
  GH_PROF_REPORT();
  return 0;
}
