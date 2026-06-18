#pragma once
// =====================================================================================
// LLR (density-of-states) wrapper around the validated GaugeHiggsHMC engine -- HARMONIC-RESTRAINT
// 1D-A scheme. Locates the precise first-order gauge-freezing beta_f(kappa) for SU(N)+frozen-Higgs
// (e.g. SU(2)->2T) where plain HMC CANNOT tunnel the strong bistability. The U(1) LLR (src/u1/u1_llr.hpp)
// is angle/heatbath-specific and does NOT generalize to SU(N)+general-rep; this REUSES the full HMC
// dynamics and adds only the energy accessor, a Gaussian restraint, and the Robbins-Monro a(A) solve.
//
// WHY A RESTRAINT, NOT A HARD WINDOW: a global HMC trajectory refreshes ALL link momenta and moves the
// whole field, so A jumps by >> any usable hard-window half-width -> 0 in-window acceptance and RM can't
// correct a1 (chicken-and-egg). A SMOOTH quadratic restraint V_w = (lw/2)(A-A0)^2 instead makes the gauge
// force self-correcting: beta_eff(U) = a1 + lw*(A(U)-A0) pulls A back toward A0 every kick, so trajectories
// stay near A0 and are accepted by the ordinary smooth dH. (Standard "LLR with Gaussian constraint";
// equivalent to the hard window in the small-width limit, RM still gives a1 -> d ln rho/dA at <A>=A0.)
//
//   A(U) = gauge_action<D,N>(U,1) = sum_pl (1-(1/N)ReTrU_pl)   [extensive, conjugate to beta]
// MATTER rides under its PHYSICAL dynamics at hmc.kappa (never tilted/constrained; A is phi-independent).
// beta_f = Maxwell equal-area level of the a1(e) curve, e=A/n_plaq (scripts/su2_llr_betaf.py).
// MERGE-SAFE: additive; touches nothing in the HMC engine (only sets hmc.beta around its own kicks).
// =====================================================================================
#include "hmc/gauge_higgs_hmc.hpp"
#include "action/gauge_wilson.hpp"
#include "action/scalar_higgs.hpp"
#include "measure/autocorr.hpp"     // tau_int / kSokalC for tau-adaptive K
#include "core/config.hpp"
#include <vector>
#include <algorithm>
#include <cmath>

namespace gh {

template <int D, int N>
struct LLRWindow {
  GaugeHiggsHMC<D, N>& hmc;
  std::uint64_t rmctr = 0;           // accept-RNG counter
  long traj_tot = 0, traj_acc = 0;   // HMC acceptance diagnostics

  explicit LLRWindow(GaugeHiggsHMC<D, N>& h) : hmc(h) {}

  // Extensive bistable gauge energy (conjugate to beta).
  Real A() const { return gauge_action<D, N>(hmc.U, 1.0); }

  Real S_matter() const {
    return hmc.potential ? scalar_action<D, N>(hmc.phi, hmc.U, *hmc.rep, hmc.kappa, *hmc.potential)
                         : scalar_action<D, N>(hmc.phi, hmc.U, *hmc.rep, hmc.kappa, hmc.lambda);
  }
  // Full restrained Hamiltonian: kinetic + a1*A + (lw/2)(A-A0)^2 + S_matter(kappa).
  Real Hrestr(Real a1, Real A0, Real lw) const {
    const Real Av = A();
    return hmc.P.kinetic() + hmc.pi.kinetic()
         + a1 * Av + 0.5 * lw * (Av - A0) * (Av - A0) + S_matter();
  }

  // Omelyan-2MN MD with the restraint folded into the gauge force: before EACH kick, set the engine's
  // beta to beta_eff = a1 + lw*(A-A0) so its gauge force = grad[a1*A + (lw/2)(A-A0)^2]. Matter force
  // uses hmc.kappa (physical) inside the same kick. Reproduces the engine's integrator structure.
  void md_restrained(Real a1, Real A0, Real lw) {
    const Real lam = kOmelyanLambda;
    const Real eps = hmc.tau / hmc.nmd;
    auto set_beta = [&] { hmc.beta = a1 + lw * (A() - A0); };
    set_beta(); hmc.kick(lam * eps);
    for (int i = 0; i < hmc.nmd; ++i) {
      hmc.drift(0.5 * eps);
      set_beta(); hmc.kick((1.0 - 2.0 * lam) * eps);
      hmc.drift(0.5 * eps);
      set_beta(); hmc.kick((i != hmc.nmd - 1 ? 2.0 * lam : lam) * eps);
    }
  }

  // One restrained HMC trajectory. Accept by the smooth restrained dH (no hard window rejection).
  bool restrained_traj(Real a1, Real A0, Real lw) {
    hmc.refresh_momenta();                       // keys off hmc.traj_count -> MUST advance it per traj
    std::vector<Cmat<N>> U_save = hmc.U.u;
    std::vector<Complex> phi_save = hmc.phi.data;
    const Real Hi = Hrestr(a1, A0, lw);
    md_restrained(a1, A0, lw);
    if (hmc.reunit_each_traj) hmc.U.reunitarize_all();
    const Real Hf = Hrestr(a1, A0, lw);
    const Real dH = Hf - Hi;
    ++traj_tot;
    const double r = hmc.rng.uniform(Rng::key(0x77, rmctr++));
    const bool acc = (dH <= 0.0) || (r < std::exp(-dH));
    if (acc) ++traj_acc; else { hmc.U.u = U_save; hmc.phi.data = phi_save; }
    ++hmc.traj_count;                            // advance the momentum-refresh RNG key (fresh momenta next traj)
    return acc;
  }

  // Seed: run nseed restrained trajectories to pull A toward A0 from wherever it is.
  void seed(Real a1, Real A0, Real lw, Real kappa, int nseed) {
    hmc.kappa = kappa;
    for (int i = 0; i < nseed; ++i) restrained_traj(a1, A0, lw);
  }

  // Robbins-Monro solve for a1(A0): pin <A>=A0. Step 12/(W^2 (m+1)). tau-ADAPTIVE K (item 1):
  // the RM block-average window K is sized PER CELL from this cell's autocorrelation -- but tau MUST be
  // measured DURING the RM at the converging a1, NOT from the seed (the seed runs at the off-branch
  // bracket a1, so the unstable back-bending cells look white-noise there and get fatally under-sampled,
  // biasing beta_f -- the gate caught exactly this). So: run a short full-Kmax WARMUP collecting the A
  // series at the (now near-converged) a1, estimate tau, then K_use = clamp(ceil(c*tau), Kfloor, Kmax)
  // for the remaining iterations. White-noise cells drop to Kfloor; back-bending cells keep ~Kmax.
  // K is the tracking window of CONTINUOUS RM -> trades variance for cost, does NOT bias the a1 fixed pt.
  Real rm_solve(Real A0, Real hw, Real lw, Real a1_0, Real kappa, int Kmax, int NRM,
                int Kfloor, bool adaptive,
                double* acc_hmc = nullptr, double* resid = nullptr, double* meanA_out = nullptr,
                int* Kused_out = nullptr) {
    hmc.kappa = kappa;
    const Real W = 2.0 * hw;
    Real a1 = a1_0;
    const int tail = std::max(1, NRM / 4);
    const int warm = adaptive ? std::max(2, NRM / 5) : NRM;        // full-Kmax warmup before adapting
    Real sumr = 0, sumA = 0; int nt = 0;
    const long t0 = traj_tot, a0 = traj_acc;
    int K_use = Kmax;
    std::vector<Real> Abuf;                                        // RM A-series (warmup) for tau
    for (int m = 0; m < NRM; ++m) {
      Real s = 0;
      for (int k = 0; k < K_use; ++k) { restrained_traj(a1, A0, lw); const Real Ak = A(); s += Ak;
                                        if (adaptive && m < warm) Abuf.push_back(Ak); }
      const Real meanA = s / K_use;
      a1 += (12.0 / (W * W * (m + 1))) * (meanA - A0);
      if (adaptive && m + 1 == warm) {                            // adapt K from the in-RM (on-branch) tau
        const Real tau = tau_int(Abuf);
        K_use = std::max(Kfloor, std::min(Kmax, static_cast<int>(std::ceil(kSokalC * tau))));
      }
      if (m + tail >= NRM) { sumr += (meanA - A0) * (meanA - A0); sumA += meanA; ++nt; }
    }
    if (acc_hmc) { const long dt = traj_tot - t0, da = traj_acc - a0; *acc_hmc = dt ? double(da) / double(dt) : 0.0; }
    if (resid)    *resid    = nt ? std::sqrt(sumr / nt) : 0.0;
    if (meanA_out) *meanA_out = nt ? sumA / nt : A();
    if (Kused_out) *Kused_out = K_use;
    return a1;
  }
};

}  // namespace gh
