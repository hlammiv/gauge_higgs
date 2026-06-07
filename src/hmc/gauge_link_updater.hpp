#pragma once
// =====================================================================================
// Tunneling-capable gauge-link updater for SU(N)+arbitrary-irrep-Higgs HMC.
//
// PROBLEM (see memory: discrete-betac-hmc-cannot-tunnel): at large kappa a center-blind
// Higgs (e.g. SU(2) spin-3, rep arg '6', d=7, frozen to 2T) drives a STRONG FIRST-ORDER
// gauge-freezing transition. The existing HMC is smooth, local, symplectic MD and CANNOT
// tunnel the extensive barrier: from a hot start the gauge stays trapped disordered while
// the equilibrium is the ordered branch. We need a Markov move that CAN cross the barrier.
//
// THE KEY SUBTLETY: classic SU(2) heatbath/overrelaxation need the action LINEAR in the
// fundamental link, S ~ Re Tr(U Staple). The Wilson term IS linear. But the matter hopping
// S_H = -2 kappa Re[phi^dag D^(R)(U) phi'] is in the spin-R rep: D^(R)(U) is a degree-(2j)
// POLYNOMIAL in U, so S_H is NONLINEAR in U and does not fit the heatbath conditional.
// At large kappa S_H pins each link near the discrete group; a naive gauge heatbath (which
// is essentially flat in the matter direction) proposes a link that the matter term then
// VETOES with probability -> 1. So a pure gauge heatbath canNOT tunnel where it matters.
//
// THE SYNTHESIS (this file): three INDIVIDUALLY-EXACT moves, composed. Each leaves the FULL
// joint Boltzmann weight exp(-S_g - S_H) invariant (detailed balance). Toggle with env;
// default behavior unchanged.
//
//   (A) CENTER FLIP  U -> z*U,  z a nontrivial center element (N=2: z=-I).
//       For a CENTER-BLIND rep [N-ality k(R) with k(R) mod N == 0] the matter term is
//       INVARIANT: D^(R)(zU) = z^{k(R)} D^(R)(U) = D^(R)(U), so dS_H == 0 EXACTLY. The move
//       is then a PURE-GAUGE Metropolis on the Wilson dS_g that the matter term can NEVER
//       veto. THIS is the natural barrier-crosser: it flips the link sign without paying any
//       matter energy. (Self-checked at runtime: aborts if the rep is not center-blind.)
//
//   (B) GAUGE HEATBATH (Kennedy-Pendleton) + MATTER METROPOLIS VETO  [N=2 only].
//       Draw U' from the EXACT Wilson conditional ~ exp[(beta/N) Re Tr(U V)], accept on the
//       matter change only (the gauge MH ratio cancels because the proposal is the exact
//       conditional). Excellent gauge decorrelation/mixing; correct at all kappa, effective
//       at small-to-moderate kappa. For N != 2 this move is skipped (no SU(N) heatbath here).
//
//   (C) FULL-ACTION MULTI-HIT METROPOLIS  U' = exp(i width g) U, g~Gaussian su(N).
//       Symmetric proposal -> plain Metropolis on the EXACT local dS_g + dS_H. No linearity
//       assumed: the nonlinear matter term is evaluated directly via rep->rotate(). The
//       robust exact fallback that tunnels at ANY kappa (slowly at very large kappa).
//
// All moves are single-link, local: only the Wilson staple of (s,mu) and the SINGLE forward
// hopping bond phi_s^dag D(U_mu(s)) phi_{s+mu} depend on U_mu(s). phi is held FIXED during a
// gauge sweep (a valid gauge-only conditional sub-kernel). Composes with HMC as an
// independent stationary kernel run BETWEEN trajectories -- no MD/reversibility interaction.
// reunitarize_all() still applies; the per-kick Dcache is rebuilt each kick so is unaffected.
// =====================================================================================
#include "hmc/gauge_higgs_hmc.hpp"
#include "action/gauge_wilson.hpp"
#include "core/linalg.hpp"
#include <cmath>
#include <cstdio>
#include <vector>

namespace gh {

// ---- Local single-link action deltas (EXACT; verified to ~1e-13 vs the global action) ----

// Wilson change for U(s,mu) -> Up (staple Sig precomputed, link-independent).
template <int D, int N>
inline Real link_dSg(const GaugeHiggsHMC<D, N>& h, const Cmat<N>& Uo, const Cmat<N>& Up,
                     const Cmat<N>& Sig) {
  Cmat<N> dU = Up; dU -= Uo;
  return -(h.beta / N) * reTrProd(dU, Sig);
}

// Matter (spin-R hopping) change for U(s,mu) -> Up. Only the single forward bond touches it.
template <int D, int N>
inline Real link_dSh(const GaugeHiggsHMC<D, N>& h, std::int64_t s, int mu,
                     const Cmat<N>& Uo, const Cmat<N>& Up) {
  const std::int64_t y = h.lat.neighbor_fwd(s, mu);
  const DVec px = h.phi.get(s), py = h.phi.get(y);
  return -2.0 * h.kappa * (dot(px, h.rep->rotate(Up, py)).real()
                         - dot(px, h.rep->rotate(Uo, py)).real());
}

// ---- (A) CENTER-FLIP SWEEP -----------------------------------------------------------------
// z = exp(2 pi i m / N) * I for some m in 1..N-1. For N=2 the only nontrivial center is -I.
// PRECONDITION: the rep must be center-blind, N-ality k(R) mod N == 0, so D^(R)(zU)=D^(R)(U)
// and dS_H == 0 exactly. We assert that and accept on the Wilson dS_g ALONE.
//   q(U->zU) = q(zU->U) since z^{-1}=z^{N-1} is also a center element (for N=2 z=z^{-1}=-I),
//   so the proposal is symmetric and Metropolis on dS_g is exact.
template <int D, int N>
inline double center_flip_sweep(GaugeHiggsHMC<D, N>& h, std::uint64_t sweepid) {
  // center-blind check (do it once -> cheap); if violated this move is INVALID, bail loudly.
  if (h.rep->nality % N != 0) {
    std::fprintf(stderr,
        "[gauge_link_updater] FATAL: center_flip on a NON-center-blind rep (N-ality %d mod %d != 0);"
        " dS_H is not 0 -> move would be wrong. Disable GH_GUPD center flips for this rep.\n",
        h.rep->nality, N);
    std::abort();
  }
  std::uint64_t acc = 0, tot = 0;
  for (std::int64_t s = 0; s < h.lat.vol; ++s)
    for (int mu = 0; mu < D; ++mu) {
      // pick a nontrivial center element z = w^m I, w = exp(2 pi i / N), m in 1..N-1.
      const int m = (N == 2) ? 1 : (1 + static_cast<int>(
          h.rng.uniform(Rng::key(0xCE71, sweepid, s, mu)) * (N - 1)));
      const Real ang = 2.0 * kPi * m / N;
      const Complex z(std::cos(ang), std::sin(ang));
      const Cmat<N> Uo = h.U(s, mu);
      Cmat<N> Up = Uo; Up *= z;                     // z * U (z is scalar*I)
      const Cmat<N> Sig = staple<D, N>(h.U, s, mu);
      const Real dSg = link_dSg<D, N>(h, Uo, Up, Sig);   // dS_H == 0 (center-blind)
      const double r = h.rng.uniform(Rng::key(0xCE72, sweepid, s, mu));
      if (dSg <= 0.0 || r < std::exp(-dSg)) { h.U(s, mu) = Up; ++acc; }
      ++tot;
    }
  return tot ? double(acc) / double(tot) : 0.0;
}

// ---- (B) SU(2) HEATBATH (Kennedy-Pendleton) + MATTER METROPOLIS VETO -----------------------
// Draw a0 = cos(theta) with measure sqrt(1-a0^2) exp(alpha a0) on [-1,1] (Kennedy-Pendleton
// a0-rejection), 3-vector uniform on S^2 with magnitude sqrt(1-a0^2). U' = X W^dag with
// W = V/sqrt(det V) in SU(2), alpha = 2 (beta/N) sqrt(det V). EXACT Wilson conditional, so the
// gauge MH ratio is 1 and cancels; accept on the matter change exp(-dS_H) only. N==2 only.
inline Cmat<2> su2_from_quat(Real a0, Real a1, Real a2, Real a3) {
  Cmat<2> U;
  U(0, 0) = Complex(a0, a3);  U(0, 1) = Complex(a2, a1);
  U(1, 0) = Complex(-a2, a1); U(1, 1) = Complex(a0, -a3);
  return U;
}
// Sample X in SU(2) ~ exp[alpha a0] measure (alpha = effective coupling). Returns X.
inline Cmat<2> su2_heatbath_draw(Real alpha, const Rng& rng, std::uint64_t key) {
  Real a0 = 0.0;
  std::uint64_t ctr = 0;
  if (alpha < 1e-12) {                 // alpha -> 0: uniform on S^3
    a0 = 2.0 * rng.uniform(Rng::key(key, ctr++)) - 1.0;
  } else {
    for (int it = 0; it < 1000; ++it) {
      Real r1 = rng.uniform(Rng::key(key, ctr++)); if (r1 < 1e-300) r1 = 1e-300;
      Real r2 = rng.uniform(Rng::key(key, ctr++)); if (r2 < 1e-300) r2 = 1e-300;
      Real r3 = rng.uniform(Rng::key(key, ctr++)); if (r3 < 1e-300) r3 = 1e-300;
      const Real cc = std::cos(2.0 * kPi * r2);
      const Real lam2 = -(1.0 / (2.0 * alpha)) * (std::log(r1) + cc * cc * std::log(r3));
      const Real rr = rng.uniform(Rng::key(key, ctr++));
      if (rr * rr <= 1.0 - lam2) { a0 = 1.0 - 2.0 * lam2; break; }
      if (it == 999) a0 = 1.0 - 2.0 * lam2;
    }
  }
  a0 = std::max(-1.0, std::min(1.0, a0));
  const Real rmag = std::sqrt(std::max(0.0, 1.0 - a0 * a0));
  const Real u  = 2.0 * rng.uniform(Rng::key(key, ctr++)) - 1.0;
  const Real ph = 2.0 * kPi * rng.uniform(Rng::key(key, ctr++));
  const Real st = std::sqrt(std::max(0.0, 1.0 - u * u));
  return su2_from_quat(a0, rmag * st * std::cos(ph), rmag * st * std::sin(ph), rmag * u);
}

// One heatbath+veto sweep. N==2 specialization checked at compile time via the caller.
template <int D>
inline double gauge_heatbath_veto_sweep(GaugeHiggsHMC<D, 2>& h, std::uint64_t sweepid) {
  std::uint64_t acc = 0, tot = 0;
  const Real coef = h.beta / 2.0;     // (beta/N), N=2
  for (std::int64_t s = 0; s < h.lat.vol; ++s)
    for (int mu = 0; mu < D; ++mu) {     // D directions (NOT N colors)
      const Cmat<2> V = staple<D, 2>(h.U, s, mu);
      const Real k = std::sqrt(std::max(std::abs(det(V)), 1e-300));
      Cmat<2> W = V; W *= Complex(1.0 / k, 0.0);          // W = V / |..| in SU(2)
      const Real alpha = 2.0 * coef * k;
      const Cmat<2> X  = su2_heatbath_draw(alpha, h.rng, Rng::key(0xAB31, sweepid, s, mu));
      Cmat<2> Up = X * W.dagger();                        // exact Wilson-conditional draw
      reunitarize<2>(Up);   // staple W=V/k carries ~1e-15 roundoff; reproject so the link stays
                            // exactly in SU(2) (otherwise the error compounds over many sweeps).
      const Cmat<2> Uo = h.U(s, mu);
      const Real dSh = link_dSh<D, 2>(h, s, mu, Uo, Up);  // matter veto ONLY (gauge ratio==1)
      const double r = h.rng.uniform(Rng::key(0xAB32, sweepid, s, mu));
      if (dSh <= 0.0 || r < std::exp(-dSh)) { h.U(s, mu) = Up; ++acc; }
      ++tot;
    }
  return tot ? double(acc) / double(tot) : 0.0;
}
// Generic-N fallback: no SU(N) heatbath implemented here (would need Cabibbo-Marinari).
template <int D, int N>
inline double gauge_heatbath_veto_sweep(GaugeHiggsHMC<D, N>& h, std::uint64_t sweepid) {
  (void)h; (void)sweepid; return -1.0;   // signal "not available"
}

// ---- (C) FULL-ACTION MULTI-HIT METROPOLIS (exact fallback, any kappa) ----------------------
// U' = exp(i width g) U, g zero-mean Gaussian su(N) (symmetric -> exact). Accept on the EXACT
// local dS_g + dS_H. n_hit hits per link random-walk the link across the matter pin at large
// kappa. (This is the same algorithm already on GaugeHiggsHMC::gauge_metropolis_sweep; kept
// here too so the orchestrator owns all three moves with consistent key prefixes.)
template <int D, int N>
inline double gauge_metro_multihit_sweep(GaugeHiggsHMC<D, N>& h, Real width, int n_hit,
                                         std::uint64_t sweepid) {
  std::uint64_t acc = 0, tot = 0;
  for (std::int64_t s = 0; s < h.lat.vol; ++s)
    for (int mu = 0; mu < D; ++mu) {
      const Cmat<N> Sig = staple<D, N>(h.U, s, mu);   // link-independent: reuse across hits
      const std::int64_t y = h.lat.neighbor_fwd(s, mu);
      const DVec px = h.phi.get(s), py = h.phi.get(y);
      for (int hh = 0; hh < n_hit; ++hh) {
        AlgVec<N> g{};
        for (int a = 0; a < n_gen<N>(); ++a)
          g[a] = width * h.rng.gauss(Rng::key(0xF131, sweepid, s, mu,
                                              static_cast<std::uint64_t>(a) * 64 + hh));
        const Cmat<N> Uo = h.U(s, mu);
        const Cmat<N> Up = expi<N>(alg_to_mat<N>(g)) * Uo;
        Cmat<N> dU = Up; dU -= Uo;
        const Real dSg = -(h.beta / N) * reTrProd(dU, Sig);
        const Real dSh = -2.0 * h.kappa * (dot(px, h.rep->rotate(Up, py)).real()
                                         - dot(px, h.rep->rotate(Uo, py)).real());
        const Real dS = dSg + dSh;
        const double r = h.rng.uniform(Rng::key(0xF132, sweepid, s, mu, hh));
        if (dS <= 0.0 || r < std::exp(-dS)) { h.U(s, mu) = Up; ++acc; }
        ++tot;
      }
    }
  return tot ? double(acc) / double(tot) : 0.0;
}

// ---- ORCHESTRATOR --------------------------------------------------------------------------
// Configuration for one composite gauge-link sweep-set, run between HMC trajectories.
struct GaugeUpdaterCfg {
  int   n_center  = 0;     // center-flip sweeps   (move A) -- the barrier crosser
  int   n_hb      = 0;     // heatbath+veto sweeps (move B) -- mixing (N=2 only)
  int   n_metro   = 0;     // full-action metro    (move C) -- exact fallback
  Real  metro_width = 0.5; // tuned-small width for move C
  int   metro_hits  = 20;  // hits/link for move C
  bool  enabled = false;
  // Parse "ncenter,nhb,nmetro[,width,hits]" (env GH_GUPD). Empty/unset -> disabled.
  static GaugeUpdaterCfg from_env(const char* s) {
    GaugeUpdaterCfg c;
    if (!s || !*s) return c;
    std::vector<Real> v; const char* p = s;
    while (*p) { v.push_back(std::atof(p)); while (*p && *p != ',') ++p; if (*p == ',') ++p; }
    if (v.size() >= 1) c.n_center = static_cast<int>(v[0]);
    if (v.size() >= 2) c.n_hb     = static_cast<int>(v[1]);
    if (v.size() >= 3) c.n_metro  = static_cast<int>(v[2]);
    if (v.size() >= 4) c.metro_width = v[3];
    if (v.size() >= 5) c.metro_hits  = static_cast<int>(v[4]);
    c.enabled = (c.n_center > 0 || c.n_hb > 0 || c.n_metro > 0);
    return c;
  }
};

// Run one set of composite gauge-link sweeps. `traj` keys the RNG streams so they never
// collide with HMC keys and are reproducible/decomposition-independent. Returns nothing;
// pass `verbose` accs out via the optional pointers if wanted.
template <int D, int N>
inline void run_gauge_update(GaugeHiggsHMC<D, N>& h, const GaugeUpdaterCfg& cfg,
                             std::uint64_t traj,
                             double* acc_center = nullptr, double* acc_hb = nullptr,
                             double* acc_metro = nullptr) {
  double ac = 0, ah = 0, am = 0; int kc = 0, kh = 0, km = 0;
  // (A) center flips first: cheap, exact, the move that can tunnel without a matter veto.
  for (int k = 0; k < cfg.n_center; ++k) {
    ac += center_flip_sweep<D, N>(h, Rng::key(0x6A01, traj, k)); ++kc;
  }
  // (B) heatbath+veto for mixing (N=2 only; -1.0 means unavailable -> ignore). Call WITHOUT
  // explicit template args so overload partial-ordering picks the SU(2) specialization for
  // N==2 (explicit <D,N> would force the generic no-op fallback that returns -1.0).
  for (int k = 0; k < cfg.n_hb; ++k) {
    double a = gauge_heatbath_veto_sweep(h, Rng::key(0x6A02, traj, k));
    if (a >= 0.0) { ah += a; ++kh; }
  }
  // (C) full-action multi-hit Metropolis fallback (exact, any kappa).
  for (int k = 0; k < cfg.n_metro; ++k) {
    am += gauge_metro_multihit_sweep<D, N>(h, cfg.metro_width, cfg.metro_hits,
                                           Rng::key(0x6A03, traj, k)); ++km;
  }
  if (acc_center) *acc_center = kc ? ac / kc : 0.0;
  if (acc_hb)     *acc_hb     = kh ? ah / kh : 0.0;
  if (acc_metro)  *acc_metro  = km ? am / km : 0.0;
}

}  // namespace gh
