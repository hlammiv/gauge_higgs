#pragma once
// Compact U(1) + charge-q FROZEN-LENGTH (|phi|=1) Higgs: local Monte Carlo.
//
// phi_x = e^{i chi_x} (unit modulus), so the on-site potential |phi|^2+lambda(|phi|^2-1)^2
// is constant and DROPS OUT. The action is
//   S = beta * sum_plaq (1 - cos t_pl)
//       - 2 kappa * sum_{x,mu} cos(chi_{x+mu} - chi_x + q theta_mu(x)).
//
// WHY THIS SAMPLER: the dynamical-modulus HMC/LLR/multicanonical stack stalls on the deep-Higgs
// first-order matter B-jump because the RADIAL mode is an unbounded, slow, radially-amplified
// coordinate (memory frozen-radial-route-large-kappa). Freezing |phi|=1 deletes that coordinate:
// the matter PHASE chi_x has an EXACT von Mises conditional p(chi_x) ∝ exp(2 kappa |R_x| cos(chi_x-arg R_x)),
// sampled by heatbath (Best-Fisher) + rejection-free overrelaxation, exact at ANY kappa.
//
// GAUGE MOVES, all built on the theta-independent staple resultant G (local_gauge_cos = Re[e^{i theta} G]):
//   - local Metropolis (continuous proposal);
//   - overrelaxation: reflect about the plaquette staple theta -> -2 arg G - theta (preserves the
//     plaquette action exactly) with a Metropolis accept on the matter term only;
//   - Z_q-sector HEATBATH: resample the link's center coset theta -> theta + 2 pi m/q (m=0..q-1) from
//     its exact conditional exp(beta Re[e^{i(theta+2pi m/q)} G]) -- the matter term is invariant under the
//     2pi/q shift so dB=0 exactly (abelian GH_GUPD; the "correlated" form of the Z_q link-jump).
//
// MULTICANONICAL (optional, ptr-gated): when `mucab` is set, a bias g(B) is added to the ACCEPT step of
// EVERY B-changing kernel (matter heatbath/overrelax, gauge Metropolis/overrelax) via exp(+(g(B+dB)-g(B)));
// the dB=0 Z_q heatbath is unaffected. Heatbath-as-proposal makes the matter accept reduce to exp(Δg) alone
// (the local Boltzmann cancels). This is the Berg-Neuhaus muca-in-ACCEPT on a BOUNDED domain |B|<=2*vol*D --
// the stable replacement for the muca-in-FORCE path that froze (memory llr-density-of-states-approach).
// B is tracked incrementally in B_run (resynced once per sweep).
//
// Conventions match u1.hpp/scan_obs.hpp EXACTLY: A=sum_plaq(1-cos) [conj. to beta]; B=sum 2cos(..) [S has -kappa*B].
#include "u1/u1.hpp"            // plaq_angle, avg_plaquette, gauge_action, wilson_loop, polyakov, MucaB
#include "core/geometry.hpp"
#include "core/rng.hpp"
#include "core/config.hpp"
#include <vector>
#include <array>
#include <cmath>
#include <complex>
#include <algorithm>

namespace gh {
namespace u1 {

// Best-Fisher-Sendra von Mises sampler, centred at 0, concentration a >= 0, returns angle in
// [-pi,pi]. Counter-based uniforms keyed by (kbase, attempt, slot) -> decomposition-independent.
inline Real draw_vonmises(Real a, const Rng& rng, std::uint64_t kbase) {
  if (a < 1e-8) return (rng.uniform(Rng::key(kbase, 0)) - 0.5) * 2.0 * kPi;   // a->0: uniform
  const Real tau = 1.0 + std::sqrt(1.0 + 4.0 * a * a);
  const Real rho = (tau - std::sqrt(2.0 * tau)) / (2.0 * a);
  const Real r   = (1.0 + rho * rho) / (2.0 * rho);
  for (int att = 1; att < 256; ++att) {
    const Real u1 = rng.uniform(Rng::key(kbase, att, 1));
    const Real u2 = rng.uniform(Rng::key(kbase, att, 2));
    const Real u3 = rng.uniform(Rng::key(kbase, att, 3));
    const Real z  = std::cos(kPi * u1);
    const Real f  = (1.0 + r * z) / (r + z);
    const Real c  = a * (r - f);
    if (c * (2.0 - c) - u2 > 0.0 || std::log(c / u2) + 1.0 - c >= 0.0) {
      const Real ff = std::min(1.0, std::max(-1.0, f));
      const Real ang = std::acos(ff);
      return (u3 > 0.5) ? ang : -ang;
    }
  }
  return 0.0;   // fallback (essentially never reached)
}

template <int D>
struct U1Frozen {
  Lattice<D> lat;
  std::vector<Real> th;     // link angles theta_mu(x)  [vol*D]
  std::vector<Real> chi;    // matter phase chi_x (phi = e^{i chi})  [vol]
  Rng rng;
  Real beta = 1.0, kappa = 0.3;
  int  q = 2;
  Real gauge_step = 0.5;    // link Metropolis proposal half-width
  int  n_or  = 3;           // matter overrelaxation sweeps per heatbath sweep
  int  n_gor = 1;           // gauge overrelaxation sweeps per Metropolis sweep
  bool zq_hop = true;       // perform the Z_q-sector heatbath sweep (q>1)
  bool matter_metro = false; // muca builds/runs: use a BROAD uniform-proposal matter Metropolis (so the
                             // bias can drive B across a first-order jump); a heatbath proposal stays trapped.
  std::uint64_t sweep_count = 0;
  std::uint64_t g_prop = 0, g_acc = 0;       // gauge Metropolis acceptance
  std::uint64_t gor_prop = 0, gor_acc = 0;   // gauge overrelaxation acceptance (matter accept)
  std::uint64_t z_prop = 0, z_flip = 0;      // Z_q heatbath: total / sectors with m!=0 chosen
  std::uint64_t m_prop = 0, m_acc = 0;       // matter heatbath/overrelax muca acceptance

  // ---- multicanonical bias in B (optional; nullptr -> unbiased) ----
  MucaB* mucab = nullptr;
  bool   muca_build = false;
  Real   B_run = 0.0;       // running global B (resynced each sweep when mucab is set)

  U1Frozen(const std::array<int, D>& ext, std::uint64_t seed)
      : lat(ext), th(static_cast<std::size_t>(lat.vol) * D, 0.0), chi(lat.vol, 0.0), rng(seed) {}

  void hot(Real sigma = 1.0) {
    for (std::int64_t s = 0; s < lat.vol; ++s) {
      for (int mu = 0; mu < D; ++mu) th[s * D + mu] = sigma * rng.gauss(Rng::key(1, s, mu));
      chi[s] = (rng.uniform(Rng::key(2, s)) - 0.5) * 2.0 * kPi;
    }
  }
  void cold() { std::fill(th.begin(), th.end(), 0.0); std::fill(chi.begin(), chi.end(), 0.0); }

  // ---- observables (conventions identical to scan_obs.hpp) ----
  Real avg_plaq() const { return avg_plaquette<D>(th, lat); }
  Real A() const { return gauge_action<D>(th, lat, 1.0); }     // extensive sum(1-cos t_pl)
  Real B() const {                                             // extensive sum 2 cos(..)
    Real s = 0.0;
    #pragma omp parallel for schedule(static) reduction(+:s)
    for (std::int64_t x = 0; x < lat.vol; ++x)
      for (int mu = 0; mu < D; ++mu) {
        const std::int64_t y = lat.neighbor_fwd(x, mu);
        s += 2.0 * std::cos(chi[y] - chi[x] + q * th[x * D + mu]);
      }
    return s;
  }
  Real link_energy() const { return B() / (2.0 * static_cast<Real>(lat.vol) * D); }  // <cos> per link
  Real cur_B() const { return mucab ? B_run : B(); }

  // Site matter resultant: S(chi_x) = -2 kappa Re[e^{-i chi_x} R_x] = -2 kappa |R_x| cos(chi_x - arg R_x).
  Complex matter_resultant(std::int64_t x) const {
    Complex R(0, 0);
    for (int mu = 0; mu < D; ++mu) {
      const std::int64_t yf = lat.neighbor_fwd(x, mu), yb = lat.neighbor_bwd(x, mu);
      R += std::polar(1.0, chi[yf] + q * th[x * D + mu]);       // forward link (x,mu)
      R += std::polar(1.0, chi[yb] - q * th[yb * D + mu]);      // backward link (x-mu,mu)
    }
    return R;
  }
  // The 2D links touching site x contribute 2 Re[e^{-i chi_x} R] to B.
  static Real site_B_contrib(Real chi_x, Complex R) {
    return 2.0 * (R.real() * std::cos(chi_x) + R.imag() * std::sin(chi_x));
  }

  // sum cos(plaq) over the 2(D-1) plaquettes containing link (x,mu), from the current th.
  Real local_gauge_cos(std::int64_t x, int mu) const {
    Real s = 0.0;
    for (int nu = 0; nu < D; ++nu) {
      if (nu == mu) continue;
      s += std::cos(plaq_angle<D>(th, lat, x, mu, nu));
      s += std::cos(plaq_angle<D>(th, lat, lat.neighbor_bwd(x, nu), mu, nu));
    }
    return s;
  }
  // Theta-INDEPENDENT staple resultant G: local_gauge_cos = Re[e^{i theta} G] for ANY theta on this link.
  Complex gauge_staple(std::int64_t x, int mu) const {
    Complex G(0, 0);
    const Real thm = th[x * D + mu];
    for (int nu = 0; nu < D; ++nu) {
      if (nu == mu) continue;
      const Real pf = plaq_angle<D>(th, lat, x, mu, nu);                        // = thm + s_fwd
      const Real pb = plaq_angle<D>(th, lat, lat.neighbor_bwd(x, nu), mu, nu);  // = -thm + s_bwd
      G += std::polar(1.0, pf - thm);        // e^{i s_fwd}
      G += std::polar(1.0, -(pb + thm));     // e^{-i s_bwd}
    }
    return G;
  }
  static Real plaq_cos_at(Complex G, Real theta) { return (std::polar(1.0, theta) * G).real(); }

  // Apply a muca-biased Metropolis accept for a B-changing move with unbiased action change dS and
  // hopping change dB. Returns true (accepted, B_run advanced) / false (rejected).
  bool muca_accept(Real dS, Real dB, std::uint64_t k) {
    Real ax = -dS;
    if (mucab) ax += Real(mucab->gval(double(B_run + dB)) - mucab->gval(double(B_run)));
    if (ax >= 0.0 || rng.uniform(k) < std::exp(ax)) { if (mucab) B_run += dB; return true; }
    return false;
  }

  // ---- one full sweep ----
  void sweep() {
    const std::uint64_t t = sweep_count;
    if (mucab) B_run = B();   // resync the running B once per sweep

    // 1. gauge link Metropolis. L(theta) = -beta Re[e^{i theta} G] - 2 kappa cos(q theta + psi).
    for (std::int64_t x = 0; x < lat.vol; ++x)
      for (int mu = 0; mu < D; ++mu) {
        const std::size_t idx = x * D + mu;
        const Complex G = gauge_staple(x, mu);
        const Real psi = chi[lat.neighbor_fwd(x, mu)] - chi[x];
        const Real th_old = th[idx];
        const Real u = rng.uniform(Rng::key(0xA1, t, x, mu, 0));
        const Real th_new = th_old + gauge_step * (2.0 * u - 1.0);
        const Real m_old = std::cos(q * th_old + psi), m_new = std::cos(q * th_new + psi);
        const Real dS = -beta * (plaq_cos_at(G, th_new) - plaq_cos_at(G, th_old))
                        - 2.0 * kappa * (m_new - m_old);
        const Real dB = 2.0 * (m_new - m_old);             // this link's B-contribution change
        ++g_prop;
        th[idx] = th_new;
        if (muca_accept(dS, dB, Rng::key(0xA2, t, x, mu, 0))) ++g_acc; else th[idx] = th_old;
      }

    // 2. gauge overrelaxation: reflect theta -> -2 arg G - theta (plaquette preserved), accept on matter.
    for (int it = 0; it < n_gor; ++it)
      for (std::int64_t x = 0; x < lat.vol; ++x)
        for (int mu = 0; mu < D; ++mu) {
          const std::size_t idx = x * D + mu;
          const Complex G = gauge_staple(x, mu);
          const Real psi = chi[lat.neighbor_fwd(x, mu)] - chi[x];
          const Real th_old = th[idx];
          const Real th_new = -2.0 * std::arg(G) - th_old;
          const Real m_old = std::cos(q * th_old + psi), m_new = std::cos(q * th_new + psi);
          const Real dS = -2.0 * kappa * (m_new - m_old);  // plaquette part is preserved by construction
          const Real dB = 2.0 * (m_new - m_old);
          ++gor_prop;
          th[idx] = th_new;
          if (muca_accept(dS, dB, Rng::key(0xA3, t, it, x * D + mu))) ++gor_acc; else th[idx] = th_old;
        }

    // 3. Z_q-sector heatbath: resample the center coset theta -> theta + 2 pi m/q from exp(beta Re[e^{i.}G]).
    //    The charge-q hopping is invariant under the 2pi/q shift -> dB=0, so the muca bias is transparent.
    if (zq_hop && q > 1)
      for (std::int64_t x = 0; x < lat.vol; ++x)
        for (int mu = 0; mu < D; ++mu) {
          const std::size_t idx = x * D + mu;
          const Complex G = gauge_staple(x, mu);
          const Real th0 = th[idx];
          Real lw[256]; Real wmax = -1e300;
          for (int m = 0; m < q; ++m) { lw[m] = beta * plaq_cos_at(G, th0 + 2.0 * kPi * m / q); wmax = std::max(wmax, lw[m]); }
          Real sum = 0.0; for (int m = 0; m < q; ++m) { lw[m] = std::exp(lw[m] - wmax); sum += lw[m]; }
          Real r = rng.uniform(Rng::key(0xB0, t, x, mu, 0)) * sum; int sel = 0;
          for (int m = 0; m < q; ++m) { r -= lw[m]; if (r <= 0.0) { sel = m; break; } }
          ++z_prop; if (sel != 0) { th[idx] = th0 + 2.0 * kPi * sel / q; ++z_flip; }
        }

    // 4. matter update.
    //    - muca with matter_metro: BROAD uniform-proposal Metropolis with the FULL action + bias, so the
    //      chain proposes both ordering AND disordering moves and the bias can drive B across a first-order
    //      jump (a heatbath proposal stays trapped: it only ever proposes the ordered direction -> can't lower B).
    //    - muca without matter_metro: heatbath proposal + accept exp(Δg) (exact: the local Boltzmann cancels).
    //    - unbiased: exact von Mises heatbath (optimal at fixed kappa).
    for (std::int64_t x = 0; x < lat.vol; ++x) {
      const Complex R = matter_resultant(x);
      if (mucab && matter_metro) {
        const Real chi_new = (rng.uniform(Rng::key(0xC0, t, x)) - 0.5) * 2.0 * kPi;
        const Real dB = site_B_contrib(chi_new, R) - site_B_contrib(chi[x], R);
        ++m_prop;                                                  // dS_local = -kappa*dB
        if (muca_accept(-kappa * dB, dB, Rng::key(0xC1, t, x))) { chi[x] = chi_new; ++m_acc; }
      } else {
        Real a = 2.0 * kappa * std::abs(R), alpha = std::arg(R);
        if (a < 0.0) { a = -a; alpha += kPi; }                     // kappa<0 safety
        const Real chi_new = alpha + draw_vonmises(a, rng, Rng::key(0xC0, t, x));
        if (mucab) {
          const Real dB = site_B_contrib(chi_new, R) - site_B_contrib(chi[x], R);
          ++m_prop;
          if (muca_accept(0.0, dB, Rng::key(0xC1, t, x))) { chi[x] = chi_new; ++m_acc; }
        } else {
          chi[x] = chi_new;
        }
      }
    }

    // 5. matter overrelaxation: chi -> 2 arg R - chi (action-preserving); with muca, accept exp(Δg).
    for (int it = 0; it < n_or; ++it)
      for (std::int64_t x = 0; x < lat.vol; ++x) {
        const Complex R = matter_resultant(x);
        const Real chi_new = 2.0 * std::arg(R) - chi[x];
        if (mucab) {
          const Real dB = site_B_contrib(chi_new, R) - site_B_contrib(chi[x], R);
          ++m_prop;
          if (muca_accept(0.0, dB, Rng::key(0xC2, t, it, x))) { chi[x] = chi_new; ++m_acc; }
        } else {
          chi[x] = chi_new;
        }
      }

    if (mucab && muca_build) mucab->wl_record(double(B_run));   // Wang-Landau record the current B
    ++sweep_count;
  }

  void thermalize(int nsweep, bool tune = true) {
    for (int i = 0; i < nsweep; ++i) {
      sweep();
      if (tune && (i % 16 == 15)) {                          // retune the Metropolis step toward ~0.5
        const double acc = g_prop ? double(g_acc) / double(g_prop) : 0.5;
        if (acc > 0.6) gauge_step *= 1.15; else if (acc < 0.4) gauge_step *= 0.85;
        gauge_step = std::min(Real(kPi), std::max(Real(0.05), gauge_step));
        g_prop = g_acc = 0;
      }
    }
  }
  double gauge_acc() const { return g_prop ? double(g_acc) / double(g_prop) : 0.0; }
  double gor_rate()  const { return gor_prop ? double(gor_acc) / double(gor_prop) : 0.0; }
  double zq_flip()   const { return z_prop ? double(z_flip) / double(z_prop) : 0.0; }
  double matter_acc() const { return m_prop ? double(m_acc) / double(m_prop) : 1.0; }
};

}  // namespace u1
}  // namespace gh
