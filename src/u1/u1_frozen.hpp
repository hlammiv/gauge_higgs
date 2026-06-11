#pragma once
// Compact U(1) + charge-q FROZEN-LENGTH (|phi|=1) Higgs: local Monte Carlo.
//
// phi_x = e^{i chi_x} (unit modulus), so the on-site potential |phi|^2+lambda(|phi|^2-1)^2
// is constant and DROPS OUT. The action is
//   S = beta * sum_plaq (1 - cos t_pl)
//       - 2 kappa * sum_{x,mu} cos(chi_{x+mu} - chi_x + q theta_mu(x)).
//
// WHY THIS SAMPLER: the dynamical-modulus HMC/LLR/multicanonical stack stalls on the
// deep-Higgs first-order matter B-jump because the RADIAL mode is an unbounded, slow,
// radially-amplified coordinate (memory frozen-radial-route-large-kappa). Freezing |phi|=1
// deletes that coordinate: the single matter PHASE chi_x has an EXACT von Mises conditional
//   p(chi_x) ∝ exp( 2 kappa |R_x| cos(chi_x - arg R_x) ),   R_x = sum over the 2D touching links,
// sampled rejection-bounded by heatbath (Best-Fisher) + rejection-free overrelaxation, exact
// at ANY kappa. Gauge links use local Metropolis + the EXACT Z_q link-jump
//   theta_mu(x) -> theta_mu(x) + 2 pi m / q   (m=1..q-1),
// which leaves the charge-q hopping cos(.. + q theta) invariant => dB = 0 exactly (the abelian
// twin of the SU(2) GH_GUPD center flip). This is the Bowler-1981 / Damgaard-Heller-1989 sampler
// class that mapped q=2,6 forty years ago; see docs/how_q2_q6_were_extracted.md.
//
// Conventions match u1.hpp / scan_obs.hpp EXACTLY: A = sum_plaq(1-cos t_pl) [extensive,
// conjugate to beta]; B = sum_{x,mu} 2 cos(chi_{x+mu}-chi_x+q theta) [extensive; S contains -kappa*B].
#include "u1/u1.hpp"            // plaq_angle, avg_plaquette, gauge_action, wilson_loop, polyakov
#include "core/geometry.hpp"
#include "core/rng.hpp"
#include "core/config.hpp"
#include <vector>
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
  int  n_or = 3;            // matter overrelaxation sweeps per heatbath sweep
  bool zq_hop = true;       // perform the exact Z_q link-jump sweep (q>1)
  std::uint64_t sweep_count = 0;
  std::uint64_t g_prop = 0, g_acc = 0;   // gauge Metropolis acceptance counters
  std::uint64_t z_prop = 0, z_acc = 0;   // Z_q link-jump acceptance counters

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

  // Site matter resultant: S(chi_x) = -2 kappa Re[e^{-i chi_x} R_x] = -2 kappa |R_x| cos(chi_x - arg R_x).
  // R_x sums the 2D links touching x (forward (x,mu) gives beta=chi_{x+mu}+q th; backward
  // (x-mu,mu) gives beta=chi_{x-mu}-q th_{x-mu}).
  Complex matter_resultant(std::int64_t x) const {
    Complex R(0, 0);
    for (int mu = 0; mu < D; ++mu) {
      const std::int64_t yf = lat.neighbor_fwd(x, mu), yb = lat.neighbor_bwd(x, mu);
      R += std::polar(1.0, chi[yf] + q * th[x * D + mu]);       // forward link (x,mu)
      R += std::polar(1.0, chi[yb] - q * th[yb * D + mu]);      // backward link (x-mu,mu)
    }
    return R;
  }

  // sum of cos(plaq) over the 2(D-1) plaquettes containing link (x,mu), from the current th.
  Real local_gauge_cos(std::int64_t x, int mu) const {
    Real s = 0.0;
    for (int nu = 0; nu < D; ++nu) {
      if (nu == mu) continue;
      s += std::cos(plaq_angle<D>(th, lat, x, mu, nu));
      s += std::cos(plaq_angle<D>(th, lat, lat.neighbor_bwd(x, nu), mu, nu));
    }
    return s;
  }

  // ---- one full sweep: gauge Metropolis -> Z_q hop -> matter heatbath -> matter overrelax ----
  void sweep() {
    const std::uint64_t t = sweep_count;

    // 1. gauge link Metropolis (serial local updates). Local action L(theta) =
    //    -beta * local_gauge_cos - 2 kappa cos(q theta + psi), psi = chi_{x+mu}-chi_x.
    for (std::int64_t x = 0; x < lat.vol; ++x)
      for (int mu = 0; mu < D; ++mu) {
        const std::size_t idx = x * D + mu;
        const Real psi = chi[lat.neighbor_fwd(x, mu)] - chi[x];
        const Real th_old = th[idx];
        const Real Lold = -beta * local_gauge_cos(x, mu) - 2.0 * kappa * std::cos(q * th_old + psi);
        const Real u = rng.uniform(Rng::key(0xA1, t, x, mu, 0));
        const Real th_new = th_old + gauge_step * (2.0 * u - 1.0);
        th[idx] = th_new;
        const Real Lnew = -beta * local_gauge_cos(x, mu) - 2.0 * kappa * std::cos(q * th_new + psi);
        const Real dS = Lnew - Lold;
        ++g_prop;
        const Real ua = rng.uniform(Rng::key(0xA2, t, x, mu, 0));
        if (dS <= 0.0 || ua < std::exp(-dS)) ++g_acc; else th[idx] = th_old;
      }

    // 2. exact Z_q link-jump theta -> theta + 2 pi m/q (q>1): matter invariant => dB=0, only plaq changes.
    if (zq_hop && q > 1)
      for (std::int64_t x = 0; x < lat.vol; ++x)
        for (int mu = 0; mu < D; ++mu) {
          const std::size_t idx = x * D + mu;
          int m = 1 + static_cast<int>(rng.uniform(Rng::key(0xB0, t, x, mu, 0)) * (q - 1));  // 1..q-1
          if (m > q - 1) m = q - 1;
          const Real th_old = th[idx];
          const Real c_old = local_gauge_cos(x, mu);
          th[idx] = th_old + 2.0 * kPi * m / q;
          const Real dS = -beta * (local_gauge_cos(x, mu) - c_old);   // matter term cancels exactly
          ++z_prop;
          const Real ua = rng.uniform(Rng::key(0xB1, t, x, mu, 0));
          if (dS <= 0.0 || ua < std::exp(-dS)) ++z_acc; else th[idx] = th_old;
        }

    // 3. matter heatbath: chi_x <- arg R_x + vonMises(0, 2 kappa |R_x|).
    for (std::int64_t x = 0; x < lat.vol; ++x) {
      const Complex R = matter_resultant(x);
      Real a = 2.0 * kappa * std::abs(R), alpha = std::arg(R);
      if (a < 0.0) { a = -a; alpha += kPi; }               // kappa<0 safety
      chi[x] = alpha + draw_vonmises(a, rng, Rng::key(0xC0, t, x));
    }

    // 4. matter overrelaxation: chi_x <- 2 arg R_x - chi_x (rejection-free, action-preserving).
    for (int it = 0; it < n_or; ++it)
      for (std::int64_t x = 0; x < lat.vol; ++x) {
        const Real alpha = std::arg(matter_resultant(x));
        chi[x] = 2.0 * alpha - chi[x];
      }

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
  double zq_acc()    const { return z_prop ? double(z_acc) / double(z_prop) : 0.0; }
};

}  // namespace u1
}  // namespace gh
