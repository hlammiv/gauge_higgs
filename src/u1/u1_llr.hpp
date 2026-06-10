#pragma once
// =====================================================================================
// 2D LLR (Linear Logarithmic Relaxation / density-of-states) engine for the compact
// U(1) + charge-q Higgs model. Computes the JOINT density of states rho(A, B) so that a
// single LLR run reweights to the WHOLE (beta,kappa) phase plane (confined / Coulomb /
// Higgs lines + the triple point). Plain HMC and replica-tempering cannot cross the strong
// first-order gauge bistability; LLR's Robbins-Monro a(E) solve resolves it (validated on
// the compact-U(1) bulk transition by Langfeld-Lucini).
//
// MERGE-SAFETY: this header is STANDALONE and ADDITIVE. It only *reuses* (never modifies)
//   gauge_action / scalar_action / scalar_force / scan_obs functions used by u1_scan and the
//   1746-test suite. Every helper below (link_dA, link_dB, site_dB, link_staple, the
//   scalar-only kick/drift/refresh, the constrained update, RM solve, ridge tiling) is NEW.
//
// CONVENTIONS (authoritative; matches src/u1/scan_obs.hpp + src/u1/u1.hpp):
//   S = beta*A - kappa*B + S_pot(phi),   exp(-S) = exp(-beta*A + kappa*B - S_pot).
//     A = plaq_energy_sum<D>(th,lat) = sum_plaq (1-cos theta_plaq)         [EXTENSIVE, conj. beta]
//     B = hop_energy_sum<D>(phi,th,lat,q) = sum_{x,mu} 2Re[conj(phi_x)e^{iq theta}phi_{x+mu}]
//         (scalar_action holds -kappa*B, so -B is conjugate to kappa).
//   WINDOWED PAIR (S1,S2) = (A, E2=-B)  [scan_obs (E1,E2)=(A,-B), (l1,l2)=(beta,kappa)].
//   Cell weight = exp(-a1*S1 - a2*S2) = exp(-a1*A + a2*B).
//   Higgs HMC half samples exp(+a2*B - S_pot)  =>  kappa_eff = +a2  (POSITIVE).
//   (One-line invariant: window E2=-B, kappa_eff=+a2, gauge weight exp(-a1*dA + a2*dB).)
// =====================================================================================
#include "u1/u1.hpp"          // U1HMC machinery, plaq_angle, scalar_action/force
#include "u1/scan_obs.hpp"    // plaq_energy_sum (=A), hop_energy_sum (=B)
#include "core/geometry.hpp"
#include "core/rng.hpp"
#include "core/config.hpp"
#include <vector>
#include <complex>
#include <cmath>
#include <cstdint>
#include <algorithm>
#include <array>

namespace gh {
namespace u1 {

// -------------------------------------------------------------------------------------
// ADDITIVE LOCAL HELPERS (NOT in u1.hpp -- merge safety). All are FD-checkable: the delta
// of plaq_energy_sum / hop_energy_sum (full recompute before/after the move) equals the
// helper return to ~1e-9. See U1LLR::fd_self_check.
// -------------------------------------------------------------------------------------

// (1) Change in A when link (x,mu) is shifted theta -> theta+dtheta.
// The link (x,mu) enters the 2(D-1) plaquettes touching it (same set as add_gauge_force):
//   - forward plaq  P_fwd = plaq_angle(x; mu,nu)            with COEFFICIENT +1 on theta,
//   - backward plaq P_bwd = plaq_angle(neighbor_bwd(x,nu); mu,nu) with COEFFICIENT -1 on theta
//     (P_bwd's third term is -th[fwd(bwd(x,nu),mu)... = -th[x,mu]).
// A = sum_plaq (1 - cos P), so
//   dA = sum_{nu!=mu}[ (cos P_fwd - cos(P_fwd + dtheta)) + (cos P_bwd - cos(P_bwd - dtheta)) ].
template <int D>
Real link_dA(const std::vector<Real>& th, const Lattice<D>& lat,
             std::int64_t x, int mu, Real dtheta) {
  Real dA = 0.0;
  for (int nu = 0; nu < D; ++nu) {
    if (nu == mu) continue;
    const Real Pf = plaq_angle<D>(th, lat, x, mu, nu);
    const std::int64_t xbnu = lat.neighbor_bwd(x, nu);
    const Real Pb = plaq_angle<D>(th, lat, xbnu, mu, nu);
    dA += (std::cos(Pf) - std::cos(Pf + dtheta)) + (std::cos(Pb) - std::cos(Pb - dtheta));
  }
  return dA;
}

// (2) Change in B from a gauge move theta_mu(x) -> theta+dtheta. Exactly ONE bond:
//   y = neighbor_fwd(x,mu);  z = conj(phi[x])*phi[y]   (complex, theta-independent)
//   dB = 2Re[z*e^{iq(theta+dtheta)}] - 2Re[z*e^{iq theta}],   theta = th[x*D+mu].
template <int D>
Real link_dB(const std::vector<Complex>& phi, const std::vector<Real>& th, const Lattice<D>& lat,
             std::int64_t x, int mu, int q, Real dtheta) {
  const std::int64_t y = lat.neighbor_fwd(x, mu);
  const Complex z = std::conj(phi[x]) * phi[y];
  const Real theta = th[x * D + mu];
  return 2.0 * (z * std::polar(1.0, q * (theta + dtheta))).real()
       - 2.0 * (z * std::polar(1.0, q * theta)).real();
}

// (3) Change in B from a scalar move phi[x]: phi_old -> phi_new (A invariant). Through the
// 2D bonds touching x: D forward (x is the "left" site) + D backward (x is the "right" site).
//   dB = sum_mu[ 2Re[(conj(phi_new)-conj(phi_old)) e^{iq t_mu(x)} phi[fwd(x,mu)]]
//              + 2Re[conj(phi[bwd(x,mu)]) e^{iq t_mu(bwd(x,mu))} (phi_new-phi_old)] ].
template <int D>
Real site_dB(const Complex& phi_new, const Complex& phi_old,
             const std::vector<Complex>& phi, const std::vector<Real>& th, const Lattice<D>& lat,
             std::int64_t x, int q) {
  const Complex dconj = std::conj(phi_new) - std::conj(phi_old);
  const Complex dphi = phi_new - phi_old;
  Real dB = 0.0;
  for (int mu = 0; mu < D; ++mu) {
    const std::int64_t yf = lat.neighbor_fwd(x, mu);
    const std::int64_t yb = lat.neighbor_bwd(x, mu);
    const Complex pf = std::polar(1.0, q * th[x * D + mu]);          // e^{iq t_mu(x)}
    const Complex pb = std::polar(1.0, q * th[yb * D + mu]);         // e^{iq t_mu(x-mu)}
    dB += 2.0 * (dconj * pf * phi[yf]).real();
    dB += 2.0 * (std::conj(phi[yb]) * pb * dphi).real();
  }
  return dB;
}

// (4) Complex staple resultant for gauge overrelaxation. The theta-dependent part of A from
// the 2(D-1) plaquettes touching link (x,mu) is
//     f(theta) = sum_{nu!=mu}[ cos(theta + Xf_nu) + cos(theta - Xb_nu) ]  (theta-independent
//                                                                          Xf,Xb)
//   Xf_nu = P_fwd - theta  (P_fwd = plaq_angle(x; mu,nu) = +theta + Xf),
//   Xb_nu = P_bwd + theta  (P_bwd = plaq_angle(bwd(x,nu); mu,nu) = -theta + Xb, and
//                           cos(P_bwd)=cos(-theta+Xb)=cos(theta-Xb)).
// So f(theta) = Re[ e^{i theta} R ] = |R| cos(theta + arg R),  R = sum (e^{i Xf} + e^{-i Xb}).
// f is invariant under the reflection theta -> -theta - 2 arg(R). We return R; the caller
// computes phi_R = -arg(R) and the A-invariant reflected angle theta_new = 2*phi_R - theta
// (FD-checked A-invariant in fd_self_check).
template <int D>
Complex link_staple(const std::vector<Real>& th, const Lattice<D>& lat,
                    std::int64_t x, int mu) {
  const Real theta = th[x * D + mu];
  Complex R(0.0, 0.0);
  for (int nu = 0; nu < D; ++nu) {
    if (nu == mu) continue;
    const Real Pf = plaq_angle<D>(th, lat, x, mu, nu);             // = theta + Xf
    const std::int64_t xbnu = lat.neighbor_bwd(x, nu);
    const Real Pb = plaq_angle<D>(th, lat, xbnu, mu, nu);          // = -theta + Xb
    R += std::polar(1.0, Pf - theta);                              // e^{i Xf}
    R += std::polar(1.0, -(Pb + theta));                           // e^{-i Xb}  (Xb = Pb + theta)
  }
  return R;
}

// =====================================================================================
// 2D LLR engine. Owns its own lattice + fields (gauge links move by Metropolis only -- NO
// link momenta; the scalar sub-sector has HMC momenta). Keeps running scalars A_run,B_run
// mutated on every accepted local move; the windowed axis-2 variable is E2 = -B.
// =====================================================================================
template <int D>
struct U1LLR {
  Lattice<D> lat;
  std::vector<Real>    th;                 // link angles [vol*D]  (NO link momenta)
  std::vector<Complex> phi, pi, Fphi;      // scalar field, HMC momentum, force [vol]
  Rng rng;

  // running EXTENSIVE scalars (no /vol)
  Real A_run = 0.0, B_run = 0.0;

  // physics parameters
  int  q = 1;
  Real lambda = 0.5;

  // gauge local-Metropolis + overrelaxation tuning
  Real delta0 = 0.6;     // Metropolis proposal half-width (auto-tuned to ~50% accept)
  int  n_hit_g = 10;     // hits per link per Metropolis sweep
  int  n_over = 2;       // overrelaxation sweeps after Metropolis

  // scalar HMC (stiff matter) tuning. HEAVY defaults: the condensate needs strong phi
  // sampling per RM iteration or <B> under-condenses (a2 biased) -- validated 2026-06-10
  // that tau_s=1, n_md_s=20, n_site=8 give a2->kappa and <B> within ~3% (weak 0.5/8/1 -> ~27% low).
  Real tau_s = 1.0;      // scalar trajectory length
  int  n_md_s = 20;      // scalar MD fine steps

  // interleaved single-site scalar Metropolis (mixing rescue)
  int  n_site = 8;       // sweeps
  Real site_w = 0.3;     // proposal half-width (auto-tuned)

  // drift trackers (incremental A_run/B_run vs full recompute)
  Real g_driftA = 0.0, g_driftB = 0.0;

  // counters (diagnostics)
  std::uint64_t traj_s = 0;
  long acc_gauge = 0, hit_gauge = 0;
  long acc_hmc = 0, try_hmc = 0;
  long guard_reject = 0;
  long acc_site = 0, hit_site = 0;
  unsigned sweep_counter = 0;

  U1LLR(const std::array<int, D>& ext, std::uint64_t seed)
      : lat(ext), th(static_cast<std::size_t>(lat.vol) * D, 0.0),
        phi(lat.vol, Complex(1.0, 0.0)), pi(lat.vol), Fphi(lat.vol), rng(seed) {}

  // ---- recompute / resync ----
  Real recompute_A() const { return plaq_energy_sum<D>(th, lat); }
  Real recompute_B() const { return hop_energy_sum<D>(phi, th, lat, q); }
  void resync_E() {
    const Real dA = std::fabs(A_run - recompute_A());
    const Real dB = std::fabs(B_run - recompute_B());
    g_driftA = std::max(g_driftA, dA);
    g_driftB = std::max(g_driftB, dB);
    A_run = recompute_A();
    B_run = recompute_B();
  }

  // ---- snapshot / restore (must save BOTH th AND phi: forgetting phi desyncs B) ----
  struct Snap { std::vector<Real> th; std::vector<Complex> phi; Real A, B; };
  Snap snapshot() const { return Snap{th, phi, A_run, B_run}; }
  void restore(const Snap& s) { th = s.th; phi = s.phi; A_run = s.A; B_run = s.B; }

  // ---- start configurations ----
  void cold() {  // ordered gauge (th=0 -> small A) + cold phi (matches U1HMC::cold_phi)
    std::fill(th.begin(), th.end(), 0.0);
    for (auto& z : phi) z = Complex(1.0, 0.0);
    A_run = recompute_A(); B_run = recompute_B();
  }
  void hot() {   // random th + random phi (matches U1HMC::hot scaling)
    for (std::int64_t s = 0; s < lat.vol; ++s) {
      for (int mu = 0; mu < D; ++mu)
        th[s * D + mu] = kPi * (2.0 * rng.uniform(Rng::key(0x40, s, mu)) - 1.0);
      phi[s] = Complex(0.3 * rng.gauss(Rng::key(0x41, s, 0)), 0.3 * rng.gauss(Rng::key(0x41, s, 1)));
    }
    A_run = recompute_A(); B_run = recompute_B();
  }

  // ---- scalar-only HMC helpers (links frozen; A must NOT move -> no matter force on links) ----
  void refresh_scalar_momenta() {
    for (std::int64_t s = 0; s < lat.vol; ++s)
      pi[s] = Complex(rng.gauss(Rng::key(0x71, traj_s, s, 0)),
                      rng.gauss(Rng::key(0x71, traj_s, s, 1)));
  }
  Real kinetic_scalar() const { Real e = 0; for (const auto& z : pi) e += std::norm(z); return 0.5 * e; }
  void kick_scalar(Real eps, Real keff) {
    scalar_force<D>(phi, th, lat, q, keff, lambda, Fphi);   // NO add_matter_force (links frozen)
    #pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < pi.size(); ++i) pi[i] += eps * Fphi[i];
  }
  void drift_scalar(Real eps) {
    #pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < phi.size(); ++i) phi[i] += eps * pi[i];
  }
  // Omelyan-2MN over the scalar sector only (port of inner_fast restricted to phi/pi).
  void scalar_omelyan(Real keff) {
    const Real lam = kOmelyanLambda;
    const Real d = tau_s / n_md_s;
    kick_scalar(lam * d, keff);
    for (int j = 0; j < n_md_s; ++j) {
      drift_scalar(0.5 * d);
      kick_scalar((1.0 - 2.0 * lam) * d, keff);
      drift_scalar(0.5 * d);
      kick_scalar((j != n_md_s - 1 ? 2.0 * lam : lam) * d, keff);
    }
  }

  // -----------------------------------------------------------------------------------
  // HYBRID constrained-2d update (ported from update_constrained2d). Box in (A, E2=-B).
  // SCHEDULE (driveIn=false): 1 gauge Metropolis sweep (n_hit_g/link) + n_over overrelax
  //   sweeps + 1 guarded Higgs HMC trajectory + n_site single-site scalar Metropolis sweeps.
  // driveIn=true (seeding): gauge Metropolis drive-in sweep + n_site single-site drive-in
  //   sweep ONLY (cheap, monotone; no HMC, no overrelaxation).
  // INHERENTLY SERIAL over sites (A_run,B_run are global running sums).
  // -----------------------------------------------------------------------------------
  void update_constrained2d(Real a1, Real a2, Real A0, Real E2_0, Real hw1, Real hw2,
                            bool driveIn, long& hit, long& acc) {
    const unsigned sweep = sweep_counter++;
    const Real Alo = A0 - hw1, Ahi = A0 + hw1;
    const Real E2lo = E2_0 - hw2, E2hi = E2_0 + hw2;

    // ---- GAUGE local Metropolis (parity order: even then odd) ----
    long lhit = 0, lacc = 0;
    for (int par = 0; par < 2; ++par)
      for (std::int64_t x = 0; x < lat.vol; ++x) {
        if (lat.parity[x] != par) continue;
        for (int mu = 0; mu < D; ++mu) {
          for (int h = 0; h < n_hit_g; ++h) {
            const Real dtheta = delta0 * (2.0 * rng.uniform(Rng::key(0x6A, x, mu, sweep, h)) - 1.0);
            const Real dA = link_dA<D>(th, lat, x, mu, dtheta);
            const Real dB = link_dB<D>(phi, th, lat, x, mu, q, dtheta);
            const Real An = A_run + dA, Bn = B_run + dB, E2n = -Bn;
            const bool inbox_now = (A_run >= Alo && A_run <= Ahi && -B_run >= E2lo && -B_run <= E2hi);
            const bool inbox_new = (An >= Alo && An <= Ahi && E2n >= E2lo && E2n <= E2hi);
            bool accept = false;
            if (driveIn && !inbox_now) {
              const Real dc = std::hypot((A_run - A0) / hw1, (-B_run - E2_0) / hw2);
              const Real dn = std::hypot((An - A0) / hw1, (E2n - E2_0) / hw2);
              accept = (dn < dc);                                  // elliptic normalized drive-in
            } else if (inbox_new) {
              const Real w = std::exp(-a1 * dA + a2 * dB);         // tilt: exp(-a1 dS1 - a2 dS2)
              accept = (w > rng.uniform(Rng::key(0x6B, x, mu, sweep, h)));
            }
            ++lhit;
            if (accept) { th[x * D + mu] += dtheta; A_run = An; B_run = Bn; ++lacc; }
          }
        }
      }
    hit += lhit; acc += lacc;
    hit_gauge += lhit; acc_gauge += lacc;

    if (!driveIn) {
      // ---- GAUGE overrelaxation (A invariant -> box test effectively on B) ----
      for (int ov = 0; ov < n_over; ++ov)
        for (int par = 0; par < 2; ++par)
          for (std::int64_t x = 0; x < lat.vol; ++x) {
            if (lat.parity[x] != par) continue;
            for (int mu = 0; mu < D; ++mu) {
              const Complex R = link_staple<D>(th, lat, x, mu);
              if (std::abs(R) < 1e-14) continue;
              const Real phiR = -std::arg(R);                       // A-invariant reflection axis
              const Real theta = th[x * D + mu];
              const Real dtheta = (2.0 * phiR - theta) - theta;     // theta_new - theta
              const Real dB = link_dB<D>(phi, th, lat, x, mu, q, dtheta);
              const Real Bn = B_run + dB, E2n = -Bn;
              if (E2n >= E2lo && E2n <= E2hi) {                     // A invariant; only B can leave
                th[x * D + mu] = 2.0 * phiR - theta;
                B_run = Bn;
              }
            }
          }

      // ---- guarded Higgs HMC trajectory (samples exp(+a2*B - S_pot), kappa_eff=+a2) ----
      const Real keff = a2;
      const Real Bsave = B_run;
      std::vector<Complex> phi0 = phi;
      refresh_scalar_momenta();
      const Real Hi = kinetic_scalar() + scalar_action<D>(phi, th, lat, q, keff, lambda);
      scalar_omelyan(keff);
      const Real Bf = recompute_B();                                // links fixed -> A unchanged
      const Real Hf = kinetic_scalar() + scalar_action<D>(phi, th, lat, q, keff, lambda);
      ++try_hmc;
      bool acc_traj = false;
      if (-Bf < E2lo || -Bf > E2hi) {
        ++guard_reject;                                             // B-window guard
      } else {
        const Real dH = Hf - Hi;
        acc_traj = (dH <= 0.0) || (rng.uniform(Rng::key(0x73, traj_s)) < std::exp(-dH));
      }
      if (acc_traj) { B_run = Bf; ++acc_hmc; } else { phi = phi0; B_run = Bsave; }
      ++traj_s;
    }

    // ---- interleaved single-site scalar Metropolis (always respects the box exactly) ----
    const Real keff = a2;
    for (int ns = 0; ns < n_site; ++ns) {
      const unsigned ssw = sweep_counter++;
      for (std::int64_t x = 0; x < lat.vol; ++x) {
        const Complex phi_new = phi[x] + Complex(site_w * (2.0 * rng.uniform(Rng::key(0x75, x, ssw, 0)) - 1.0),
                                                 site_w * (2.0 * rng.uniform(Rng::key(0x75, x, ssw, 1)) - 1.0));
        const Real dB = site_dB<D>(phi_new, phi[x], phi, th, lat, x, q);
        const Real n2new = std::norm(phi_new), n2old = std::norm(phi[x]);
        const Real dSpot = (n2new + lambda * (n2new - 1.0) * (n2new - 1.0))
                         - (n2old + lambda * (n2old - 1.0) * (n2old - 1.0));
        const Real Bn = B_run + dB, E2n = -Bn;
        ++hit_site;
        bool accept = false;
        if (driveIn) {
          // drive-in: pull B toward window via the same elliptic metric (A unaffected by phi)
          const bool inbox_now = (-B_run >= E2lo && -B_run <= E2hi);
          if (!inbox_now) {
            const Real dc = std::fabs((-B_run - E2_0) / hw2);
            const Real dn = std::fabs((E2n - E2_0) / hw2);
            accept = (dn < dc);
          }
        } else {
          if (E2n >= E2lo && E2n <= E2hi) {
            const Real w = std::exp(keff * dB - dSpot);             // exp(+a2 dB - dS_pot)
            accept = (w > rng.uniform(Rng::key(0x76, x, ssw)));
          }
        }
        if (accept) { phi[x] = phi_new; B_run = Bn; ++acc_site; }
      }
    }
  }

  // ---- delta0 / site_w auto-tune toward ~50% accept (call between sweeps) ----
  void autotune_gauge() {
    if (hit_gauge <= 0) return;
    const double r = double(acc_gauge) / double(hit_gauge);
    if (r > 0.6) delta0 *= 1.1; else if (r < 0.4) delta0 *= 0.9;
    acc_gauge = hit_gauge = 0;
  }
  void autotune_site() {
    if (hit_site <= 0) return;
    const double r = double(acc_site) / double(hit_site);
    if (r > 0.6) site_w *= 1.1; else if (r < 0.4) site_w *= 0.9;
    acc_site = hit_site = 0;
  }

  // -----------------------------------------------------------------------------------
  // Drive (A,E2) into the cell box (greedy). Returns reachability. CAP from caller.
  // -----------------------------------------------------------------------------------
  bool seed_cell(Real A0, Real E2_0, Real hw1, Real hw2, Real a1, Real a2, unsigned CAP) {
    auto inbox = [&] {
      return A_run >= A0 - hw1 && A_run <= A0 + hw1 && -B_run >= E2_0 - hw2 && -B_run <= E2_0 + hw2;
    };
    long h = 0, a = 0; unsigned drive = 0;
    while (!inbox() && drive < CAP) {
      update_constrained2d(a1, a2, A0, E2_0, hw1, hw2, /*driveIn=*/true, h, a);
      ++drive;
    }
    resync_E();
    return inbox();
  }

  // -----------------------------------------------------------------------------------
  // Coupled Robbins-Monro (ported VERBATIM from rm_solve2d). W=2*hw; a += 12/(W^2 (m+1))(<.>-target).
  // Returns the FINAL iterate. Diagnostics fill the optional RMDiag.
  // -----------------------------------------------------------------------------------
  struct RMDiag {
    Real res_A = 0, res_E2 = 0;       // trailing RMS residual of (<A>-A0),(<E2>-E2_0)
    Real slope_a1 = 0, slope_a2 = 0;  // trailing slope of a1,a2 vs m
    Real acc_gauge = 0, acc_hmc = 0;  // trailing accept rates
    int  nwiden = 0;                  // adaptive-hw2 widenings used
  };

  void rm_solve2d(Real A0, Real E2_0, Real hw1, Real hw2, Real a1_0, Real a2_0,
                  unsigned K, unsigned NRM, Real& a1_out, Real& a2_out, RMDiag* diag = nullptr) {
    const Real W1 = 2.0 * hw1, W2 = 2.0 * hw2;
    Real a1 = a1_0, a2 = a2_0;
    const unsigned tail = std::max(1u, NRM / 4);
    Real sumr1 = 0, sumr2 = 0; int ntail = 0;
    std::vector<Real> a1hist, a2hist;
    long g_hit = 0, g_acc = 0, h_try = 0, h_acc = 0;
    for (unsigned m = 0; m < NRM; ++m) {
      Real s1 = 0, s2 = 0; long hit = 0, acc = 0;
      const long t0 = try_hmc, ta0 = acc_hmc;
      const long gh = hit_gauge, ga = acc_gauge;
      for (unsigned k = 0; k < K; ++k) {
        update_constrained2d(a1, a2, A0, E2_0, hw1, hw2, /*driveIn=*/false, hit, acc);
        s1 += A_run; s2 += (-B_run);
      }
      const Real meanA = s1 / K, meanE2 = s2 / K;
      a1 += (12.0 / (W1 * W1 * (m + 1))) * (meanA - A0);
      a2 += (12.0 / (W2 * W2 * (m + 1))) * (meanE2 - E2_0);
      resync_E();                                                   // kill FP drift
      if (m + tail >= NRM) {
        sumr1 += (meanA - A0) * (meanA - A0);
        sumr2 += (meanE2 - E2_0) * (meanE2 - E2_0);
        ++ntail;
        a1hist.push_back(a1); a2hist.push_back(a2);
        g_hit += (hit_gauge - gh); g_acc += (acc_gauge - ga);
        h_try += (try_hmc - t0); h_acc += (acc_hmc - ta0);
      }
    }
    a1_out = a1; a2_out = a2;
    if (diag) {
      diag->res_A = ntail ? std::sqrt(sumr1 / ntail) : 0.0;
      diag->res_E2 = ntail ? std::sqrt(sumr2 / ntail) : 0.0;
      const int nh = (int)a1hist.size();
      if (nh >= 2) {
        diag->slope_a1 = (a1hist.back() - a1hist.front()) / (nh - 1);
        diag->slope_a2 = (a2hist.back() - a2hist.front()) / (nh - 1);
      }
      diag->acc_gauge = g_hit ? double(g_acc) / double(g_hit) : 0.0;
      diag->acc_hmc = h_try ? double(h_acc) / double(h_try) : 0.0;
    }
  }

  // ===================================================================================
  // 1D-A LLR: window ONLY the gauge action A; the matter field phi runs under its
  // PHYSICAL dynamics at fixed kappa (never constrained, never tilted). Avoids the
  // phi-conditional-at-fixed-B problem of the 2D scheme (where windowing B forced a
  // wrong conditional). At a cell, RM tunes a1 (= d ln rho_kappa/dA) to pin <A>=A0;
  // a1 -> beta at the equilibrium A for (beta,kappa). The cell also measures the
  // PHYSICAL <B> at that (A,kappa) as a reweightable observable.
  // -----------------------------------------------------------------------------------
  void update_constrained1d(Real a1, Real kappa, Real A0, Real hw1, bool driveIn, long& hit, long& acc) {
    const unsigned sweep = sweep_counter++;
    const Real Alo = A0 - hw1, Ahi = A0 + hw1;
    // ---- GAUGE local Metropolis: tilt a1 on A, PHYSICAL +kappa on the hop dB, window A only ----
    long lhit = 0, lacc = 0;
    for (int par = 0; par < 2; ++par)
      for (std::int64_t x = 0; x < lat.vol; ++x) {
        if (lat.parity[x] != par) continue;
        for (int mu = 0; mu < D; ++mu)
          for (int h = 0; h < n_hit_g; ++h) {
            const Real dtheta = delta0 * (2.0 * rng.uniform(Rng::key(0x6A, x, mu, sweep, h)) - 1.0);
            const Real dA = link_dA<D>(th, lat, x, mu, dtheta);
            const Real dB = link_dB<D>(phi, th, lat, x, mu, q, dtheta);
            const Real An = A_run + dA;
            bool accept = false;
            if (driveIn && !(A_run >= Alo && A_run <= Ahi)) {
              accept = std::fabs(An - A0) < std::fabs(A_run - A0);     // drive A into the slab
            } else if (An >= Alo && An <= Ahi) {
              const Real w = std::exp(-a1 * dA + kappa * dB);          // tilt a1 on A; physical kappa on B
              accept = (w > rng.uniform(Rng::key(0x6B, x, mu, sweep, h)));
            }
            ++lhit;
            if (accept) { th[x * D + mu] += dtheta; A_run = An; B_run += dB; ++lacc; }
          }
      }
    hit += lhit; acc += lacc;
    hit_gauge += lhit; acc_gauge += lacc;

    if (!driveIn) {
      // ---- GAUGE overrelaxation (A invariant; B free -> always accept) ----
      for (int ov = 0; ov < n_over; ++ov)
        for (int par = 0; par < 2; ++par)
          for (std::int64_t x = 0; x < lat.vol; ++x) {
            if (lat.parity[x] != par) continue;
            for (int mu = 0; mu < D; ++mu) {
              const Complex R = link_staple<D>(th, lat, x, mu);
              if (std::abs(R) < 1e-14) continue;
              const Real phiR = -std::arg(R);
              const Real theta = th[x * D + mu];
              const Real dtheta = (2.0 * phiR - theta) - theta;
              B_run += link_dB<D>(phi, th, lat, x, mu, q, dtheta);     // A invariant; B free
              th[x * D + mu] = 2.0 * phiR - theta;
            }
          }
      // ---- PHYSICAL scalar HMC at kappa (no B-window guard: matter is unconstrained) ----
      const Real keff = kappa;
      std::vector<Complex> phi0 = phi;
      refresh_scalar_momenta();
      const Real Hi = kinetic_scalar() + scalar_action<D>(phi, th, lat, q, keff, lambda);
      scalar_omelyan(keff);
      const Real Hf = kinetic_scalar() + scalar_action<D>(phi, th, lat, q, keff, lambda);
      ++try_hmc;
      const Real dH = Hf - Hi;
      const bool acc_traj = (dH <= 0.0) || (rng.uniform(Rng::key(0x73, traj_s)) < std::exp(-dH));
      if (acc_traj) { ++acc_hmc; B_run = recompute_B(); } else { phi = phi0; }   // A unchanged either way
      ++traj_s;
    }
    // ---- PHYSICAL single-site scalar Metropolis at kappa (no B-window) ----
    const Real keff = kappa;
    for (int ns = 0; ns < n_site; ++ns) {
      const unsigned ssw = sweep_counter++;
      for (std::int64_t x = 0; x < lat.vol; ++x) {
        const Complex phi_new = phi[x] + Complex(site_w * (2.0 * rng.uniform(Rng::key(0x75, x, ssw, 0)) - 1.0),
                                                 site_w * (2.0 * rng.uniform(Rng::key(0x75, x, ssw, 1)) - 1.0));
        const Real dB = site_dB<D>(phi_new, phi[x], phi, th, lat, x, q);
        const Real n2new = std::norm(phi_new), n2old = std::norm(phi[x]);
        const Real dSpot = (n2new + lambda * (n2new - 1.0) * (n2new - 1.0))
                         - (n2old + lambda * (n2old - 1.0) * (n2old - 1.0));
        ++hit_site;
        const Real w = std::exp(keff * dB - dSpot);                    // physical exp(+kappa dB - dS_pot)
        if (w > rng.uniform(Rng::key(0x76, x, ssw))) { phi[x] = phi_new; B_run += dB; ++acc_site; }
      }
    }
  }

  // Drive A into [A0-hw1, A0+hw1] (1D slab; matter rides physically).
  bool seed_slab(Real A0, Real hw1, Real a1, Real kappa, unsigned CAP) {
    long h = 0, a = 0; unsigned drive = 0;
    while (!(A_run >= A0 - hw1 && A_run <= A0 + hw1) && drive < CAP) {
      update_constrained1d(a1, kappa, A0, hw1, /*driveIn=*/true, h, a);
      ++drive;
    }
    resync_E();
    return A_run >= A0 - hw1 && A_run <= A0 + hw1;
  }

  // RM on a1 only (pin <A>=A0); accumulate the PHYSICAL <B> at this cell as an observable.
  void rm_solve1d(Real A0, Real hw1, Real a1_0, Real kappa, unsigned K, unsigned NRM,
                  Real& a1_out, Real& Bmean_out, RMDiag* diag = nullptr) {
    const Real W1 = 2.0 * hw1;
    Real a1 = a1_0;
    const unsigned tail = std::max(1u, NRM / 4);
    Real sumr1 = 0, sumB = 0; int ntail = 0;
    std::vector<Real> a1hist;
    long g_hit = 0, g_acc = 0, h_try = 0, h_acc = 0;
    for (unsigned m = 0; m < NRM; ++m) {
      Real s1 = 0, sB = 0; long hit = 0, acc = 0;
      const long t0 = try_hmc, ta0 = acc_hmc, gh = hit_gauge, ga = acc_gauge;
      for (unsigned k = 0; k < K; ++k) {
        update_constrained1d(a1, kappa, A0, hw1, /*driveIn=*/false, hit, acc);
        s1 += A_run; sB += B_run;
      }
      const Real meanA = s1 / K;
      a1 += (12.0 / (W1 * W1 * (m + 1))) * (meanA - A0);
      resync_E();
      if (m + tail >= NRM) {
        sumr1 += (meanA - A0) * (meanA - A0);
        sumB += sB / K; ++ntail;
        a1hist.push_back(a1);
        g_hit += (hit_gauge - gh); g_acc += (acc_gauge - ga);
        h_try += (try_hmc - t0); h_acc += (acc_hmc - ta0);
      }
    }
    a1_out = a1; Bmean_out = ntail ? sumB / ntail : B_run;
    if (diag) {
      diag->res_A = ntail ? std::sqrt(sumr1 / ntail) : 0.0;
      diag->res_E2 = 0; diag->slope_a2 = 0; diag->nwiden = 0;
      const int nh = (int)a1hist.size();
      if (nh >= 2) diag->slope_a1 = (a1hist.back() - a1hist.front()) / (nh - 1);
      diag->acc_gauge = g_hit ? double(g_acc) / double(g_hit) : 0.0;
      diag->acc_hmc = h_try ? double(h_acc) / double(h_try) : 0.0;
    }
  }

  // -----------------------------------------------------------------------------------
  // FD self-check of the four local helpers (CI gate; STEP 1 of the smoke test).
  // Returns the max |helper - full-recompute-delta| over the sampled moves.
  // -----------------------------------------------------------------------------------
  Real fd_self_check(int nlinks = 50, int nsites = 50) {
    Real maxerr = 0.0;
    // randomize a config
    for (std::int64_t s = 0; s < lat.vol; ++s) {
      for (int mu = 0; mu < D; ++mu)
        th[s * D + mu] = 1.3 * (2.0 * rng.uniform(Rng::key(0xC0, s, mu)) - 1.0);
      phi[s] = Complex(rng.gauss(Rng::key(0xC1, s, 0)), rng.gauss(Rng::key(0xC1, s, 1)));
    }
    // link_dA + link_dB
    for (int t = 0; t < nlinks; ++t) {
      const std::int64_t x = std::int64_t(lat.vol * rng.uniform(Rng::key(0xC2, t, 0))) % lat.vol;
      const int mu = int(D * rng.uniform(Rng::key(0xC2, t, 1))) % D;
      const Real dtheta = 0.7 * (2.0 * rng.uniform(Rng::key(0xC2, t, 2)) - 1.0);
      const Real A_before = recompute_A(), B_before = recompute_B();
      const Real dA_h = link_dA<D>(th, lat, x, mu, dtheta);
      const Real dB_h = link_dB<D>(phi, th, lat, x, mu, q, dtheta);
      th[x * D + mu] += dtheta;
      const Real dA_f = recompute_A() - A_before, dB_f = recompute_B() - B_before;
      th[x * D + mu] -= dtheta;
      maxerr = std::max(maxerr, std::fabs(dA_h - dA_f));
      maxerr = std::max(maxerr, std::fabs(dB_h - dB_f));
      // overrelaxation A-invariance: theta -> 2 phiR - theta leaves A unchanged
      const Complex Rst = link_staple<D>(th, lat, x, mu);
      if (std::abs(Rst) > 1e-12) {
        const Real phiR = -std::arg(Rst);
        const Real dtr = (2.0 * phiR - th[x * D + mu]) - th[x * D + mu];
        maxerr = std::max(maxerr, std::fabs(link_dA<D>(th, lat, x, mu, dtr)));
      }
    }
    // site_dB
    for (int t = 0; t < nsites; ++t) {
      const std::int64_t x = std::int64_t(lat.vol * rng.uniform(Rng::key(0xC3, t, 0))) % lat.vol;
      const Complex phi_new = phi[x] + Complex(0.5 * (2.0 * rng.uniform(Rng::key(0xC3, t, 1)) - 1.0),
                                               0.5 * (2.0 * rng.uniform(Rng::key(0xC3, t, 2)) - 1.0));
      const Real B_before = recompute_B();
      const Real dB_h = site_dB<D>(phi_new, phi[x], phi, th, lat, x, q);
      const Complex save = phi[x]; phi[x] = phi_new;
      const Real dB_f = recompute_B() - B_before;
      phi[x] = save;
      maxerr = std::max(maxerr, std::fabs(dB_h - dB_f));
    }
    return maxerr;
  }
};

}  // namespace u1
}  // namespace gh
