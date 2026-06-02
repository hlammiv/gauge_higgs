// Per-(q,L) nmd PRE-TUNE driver for the U(1) charge-q campaign.
//
// WHY: leaving the per-point auto-tuner on during a campaign is a trap -- at a
// too-coarse start nmd, tune_nmd ramps nmd up (x8/iter, up to nmd_max) and burns
// 30-traj calibration bursts at those inflated step counts, so a single point can
// balloon into hours of calibration that isn't production. With the multi-timescale
// integrator (n_scalar>1) the gauge-timescale nmd needed for a target acceptance is
// nearly FLAT across the whole (beta,kappa) plane, so we tune it ONCE per (L,q) at a
// representative (stiff) grid corner and reuse it across that grid with autotune OFF.
//
// This driver thermalizes at one (beta,kappa) point, runs tune_nmd, prints the tuned
// integer nmd to STDOUT (one line, machine-readable); diagnostics go to STDERR.
//
//   ./build/u1_pretune <L> <beta> <kappa> <lambda> <q>
//                      [ntherm nmd_start tau seed n_scalar acc_lo acc_hi]
#include "u1/u1.hpp"
#include "u1/autotune.hpp"
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <array>

using namespace gh;

int main(int argc, char** argv) {
  if (argc < 6) {
    std::fprintf(stderr,
      "usage: %s <L> <beta> <kappa> <lambda> <q> "
      "[ntherm nmd_start tau seed n_scalar acc_lo acc_hi]\n", argv[0]);
    return 1;
  }
  auto af = [&](int i, double d) { return i < argc ? std::atof(argv[i]) : d; };
  auto ai = [&](int i, long d)   { return i < argc ? std::atol(argv[i]) : d; };
  const int    L         = static_cast<int>(ai(1, 8));
  const Real   beta      = af(2, 1.0);
  const Real   kappa     = af(3, 0.3);
  const Real   lambda    = af(4, 0.5);
  const int    q         = static_cast<int>(ai(5, 2));
  const int    ntherm    = static_cast<int>(ai(6, 200));
  const int    nmd_start = static_cast<int>(ai(7, 6));
  const Real   tau       = af(8, 1.0);
  const std::uint64_t seed = static_cast<std::uint64_t>(ai(9, 12345));
  const int    n_scalar  = static_cast<int>(ai(10, 6));
  const double acc_lo    = af(11, 0.70);   // slightly conservative band so the reused
  const double acc_hi    = af(12, 0.92);   // nmd is safe across the whole grid

  std::array<int, kDim> ext{}; for (int mu = 0; mu < kDim; ++mu) ext[mu] = L;
  u1::U1HMC<kDim> hmc(ext, seed);
  hmc.beta = beta; hmc.kappa = kappa; hmc.lambda = lambda; hmc.q = q;
  hmc.tau = tau; hmc.nmd = nmd_start; hmc.n_scalar = n_scalar;
  hmc.hot(0.8); hmc.cold_phi(0.5);
  for (int t = 0; t < ntherm; ++t) hmc.trajectory();

  const u1::TuneResult tr = u1::tune_nmd<kDim>(hmc, 30, acc_lo, acc_hi);
  std::fprintf(stderr,
    "# pretune L=%d^%d q=%d beta=%.4f kappa=%.4f lambda=%.4g n_scalar=%d "
    "-> nmd=%d (acc=%.3f in_band=%d iters=%d)\n",
    L, kDim, q, beta, kappa, lambda, n_scalar,
    hmc.nmd, tr.acceptance, (int)tr.in_band, tr.iters);
  std::fprintf(stdout, "%d\n", hmc.nmd);   // machine-readable: the tuned nmd
  return 0;
}
