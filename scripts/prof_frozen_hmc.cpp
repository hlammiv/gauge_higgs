// Profiling harness for the frozen-Higgs SU(2) HMC. Runs GaugeHiggsHMC<4,2> with the
// multi-invariant potential (BT/BO/BI campaign config) and dumps the GH_PROFILE breakdown,
// isolating per-MD-step force cost. Build WITH -DGH_PROFILE; scratch/diagnostic only.
//   g++ -std=c++17 -O3 -march=native -funroll-loops -fopenmp -Isrc -DNDIM=4 -DNCOL=2 \
//       -DGH_PROFILE -o build/prof_frozen_hmc scripts/prof_frozen_hmc.cpp
//   OMP_NUM_THREADS=1 ./build/prof_frozen_hmc <rep_rows=6> <L=6> <ntraj=40> <nmd=20> <frozen=1> <seed=1>
#include "hmc/gauge_higgs_hmc.hpp"
#include "action/scalar_invariants.hpp"
#include "rep/rep_general.hpp"
#include "core/profile.hpp"
#include <cstdio>
#include <cstdlib>
#include <string>
#include <sstream>
#include <vector>
#include <chrono>

using namespace gh;

int main(int argc, char** argv) {
  std::string repspec = argc > 1 ? argv[1] : "6";   // SU(2) spin-j -> rows=[2j]; spin-3 = "6"
  int L      = argc > 2 ? std::atoi(argv[2]) : 6;
  int ntraj  = argc > 3 ? std::atoi(argv[3]) : 40;
  int nmd    = argc > 4 ? std::atoi(argv[4]) : 20;
  bool froz  = argc > 5 ? std::atoi(argv[5]) != 0 : true;
  std::uint64_t seed = argc > 6 ? std::atoll(argv[6]) : 1;

  std::vector<int> rows; std::stringstream ss(repspec); std::string t;
  while (std::getline(ss, t, ',')) if (!t.empty()) rows.push_back(std::atoi(t.c_str()));
  GeneralRep<2> rep(rows);

  CasimirChannels<2> ch(rep);
  std::vector<Real> f(ch.n_channels(), 0.0);
  for (int i = 0; i < (int)f.size(); ++i) f[i] = 0.1 + 0.05 * i;   // representative coupling vector
  MultiInvariantPotential<2> pot(ch, f, 1.0);

  std::array<int,4> ext{}; for (int m=0;m<4;++m) ext[m]=L;
  GaugeHiggsHMC<4,2> hmc(ext, rep, seed);
  hmc.beta=1.5; hmc.kappa=1.0; hmc.tau=1.0; hmc.nmd=nmd; hmc.potential=&pot;
  hmc.frozen_phi = froz;
  hmc.U.hot(hmc.rng, 0.8); hmc.phi.gaussian(hmc.rng, 12345, rep.real, 0.3);
  if (froz) hmc.normalize_phi();

  auto t0 = std::chrono::steady_clock::now();
  for (int i=0;i<ntraj;++i) hmc.trajectory();
  auto t1 = std::chrono::steady_clock::now();
  double wall = std::chrono::duration<double>(t1-t0).count();

  std::printf("# rep=%s d=%d L=%d ntraj=%d nmd=%d frozen=%d use_fast=%d link_cache=%d\n",
              rep.name().c_str(), rep.d, L, ntraj, nmd, (int)froz,
              (int)rep.use_fast, (int)hmc.use_link_cache);
  std::printf("# WALL %.4f s total, %.4f s/traj, accept=%.3f\n",
              wall, wall/ntraj, hmc.acceptance());
  GH_PROF_REPORT();
  return 0;
}
