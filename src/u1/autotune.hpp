#pragma once
// Per-point HMC step-count (nmd) auto-tuner for U1HMC.
//
// The MD step size is dt = tau/nmd. In the stiff deep-Higgs region (large kappa,
// and large charge q) the scalar force is large, so a fixed nmd makes dH blow up
// and the Metropolis acceptance collapse to ~0 (seen in the Stage-0 pilot at
// kappa=0.6). tune_nmd runs short calibration bursts and adjusts hmc.nmd until the
// acceptance lands in [acc_lo, acc_hi], capped at [nmd_min, nmd_max].
//
// The caller must THERMALIZE hmc first. The calibration trajectories are extra
// thermalization (they advance the chain); the tuner ONLY sets hmc.nmd and does
// NOT bias production -- the caller runs production with the tuned nmd after
// resetting the accept counters. Deterministic given the hmc RNG state.
#include "u1/u1.hpp"
#include <cmath>
#include <cstdint>

namespace gh {
namespace u1 {

struct TuneResult { int nmd; double acceptance; int iters; bool in_band; };

template <int D>
inline TuneResult tune_nmd(U1HMC<D>& hmc, int n_cal = 30, double acc_lo = 0.65,
                           double acc_hi = 0.92, int nmd_min = 4, int nmd_max = 400,
                           int max_iter = 10) {
  const double dH_target = 0.1;            // <dH> giving acceptance ~0.75 (mid band)
  if (n_cal < 1) n_cal = 1;
  if (hmc.nmd < nmd_min) hmc.nmd = nmd_min;
  if (hmc.nmd > nmd_max) hmc.nmd = nmd_max;
  double a = 0.0;
  int rounds = 0;
  for (int it = 0; it < max_iter; ++it) {
    // One calibration burst at the current nmd; measure acceptance via counter
    // deltas (do NOT reset hmc counters -- keeps the RNG stream advancing).
    const std::uint64_t a0 = hmc.accept_count, t0 = hmc.traj_count;
    double sumdH = 0.0;
    for (int i = 0; i < n_cal; ++i) { hmc.trajectory(); sumdH += hmc.last_dH; }
    ++rounds;
    const std::uint64_t dt = hmc.traj_count - t0;
    a = dt ? double(hmc.accept_count - a0) / double(dt) : 0.0;
    const double meanDH = sumdH / n_cal;
    if (a >= acc_lo && a <= acc_hi) break;                       // in band -> done

    int nmd_new = hmc.nmd;
    if (a < acc_lo) {
      // Steps too coarse -> raise nmd. For a 2nd-order integrator <dH> ~ nmd^{-4},
      // so to cut <dH> by F we scale nmd by F^{1/4}; floor the growth (>=1.3x) to
      // guarantee progress and cap it (<=8x) so a pathological burst can't jump
      // straight to nmd_max.
      double ratio = (meanDH > dH_target) ? meanDH / dH_target : 1.0;
      double factor = std::pow(ratio, 0.25);
      if (!(factor >= 1.3)) factor = 1.3;   // also catches NaN
      if (factor > 8.0) factor = 8.0;
      nmd_new = static_cast<int>(std::ceil(hmc.nmd * factor));
    } else {
      // Too fine (acceptance above the band) -> lower nmd for efficiency.
      nmd_new = static_cast<int>(std::floor(hmc.nmd * 0.8));
    }
    if (nmd_new < nmd_min) nmd_new = nmd_min;
    if (nmd_new > nmd_max) nmd_new = nmd_max;
    if (nmd_new == hmc.nmd) break;          // clamped/stalled -> cannot improve
    hmc.nmd = nmd_new;
  }
  return TuneResult{ hmc.nmd, a, rounds, (a >= acc_lo && a <= acc_hi) };
}

}  // namespace u1
}  // namespace gh
