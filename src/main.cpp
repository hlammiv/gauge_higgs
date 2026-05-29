// Driver for the arbitrary-D SU(N) gauge (+ Higgs, forthcoming) HMC.
// Compile-time D = NDIM, N = NCOL (set via Makefile). Pure-gauge run for now.
//   ./build/gh_hmc [L=8] [beta=2.3] [ntherm=50] [nmeas=100] [nmd=20] [tau=1.0] [seed=1]
#include "hmc/hmc.hpp"
#include <cstdio>
#include <cstdlib>
#include <string>

using namespace gh;

int main(int argc, char** argv) {
  auto argf = [&](int i, double def) { return i < argc ? std::atof(argv[i]) : def; };
  auto argi = [&](int i, long def)   { return i < argc ? std::atol(argv[i]) : def; };

  const int    Lext   = static_cast<int>(argi(1, 8));
  const Real   beta   = argf(2, 2.3);
  const int    ntherm = static_cast<int>(argi(3, 50));
  const int    nmeas  = static_cast<int>(argi(4, 100));
  const int    nmd    = static_cast<int>(argi(5, 20));
  const Real   tau    = argf(6, 1.0);
  const std::uint64_t seed = static_cast<std::uint64_t>(argi(7, 1));

  std::array<int, kDim> L{};
  for (int mu = 0; mu < kDim; ++mu) L[mu] = Lext;

  GaugeHMC<kDim, kN> hmc(L, seed);
  hmc.beta = beta; hmc.tau = tau; hmc.nmd = nmd;

  std::printf("# Pure-gauge HMC: D=%d  SU(%d)  L=%d^%d  vol=%lld  beta=%.4f  tau=%.3f  nmd=%d\n",
              kDim, kN, Lext, kDim, (long long)hmc.lat.vol, beta, tau, nmd);
  hmc.U.hot(hmc.rng, 0.8);

  for (int t = 0; t < ntherm; ++t) hmc.trajectory();
  std::printf("# thermalized %d trajectories; acceptance %.3f\n", ntherm, hmc.acceptance());

  hmc.traj_count = 0; hmc.accept_count = 0;
  double sP = 0.0, sP2 = 0.0, sExp = 0.0;
  for (int t = 0; t < nmeas; ++t) {
    hmc.trajectory();
    const double P = avg_plaquette<kDim, kN>(hmc.U);
    sP += P; sP2 += P * P; sExp += std::exp(-hmc.last_dH);
  }
  const double mP = sP / nmeas;
  const double varP = sP2 / nmeas - mP * mP;
  std::printf("# measured %d trajectories\n", nmeas);
  std::printf("plaquette       = %.6f +/- %.6f\n", mP, std::sqrt(std::max(0.0, varP) / nmeas));
  std::printf("acceptance      = %.4f\n", hmc.acceptance());
  std::printf("<exp(-dH)>      = %.5f  (should be ~1)\n", sExp / nmeas);
  std::printf("max ||U^dU-I||  = %.3e\n", hmc.U.max_unitarity_violation());
  return 0;
}
