#pragma once
// Compact U(1) lattice gauge theory + charge-q complex Higgs (the abelian Fradkin-Shenker
// case: a condensing charge-q scalar breaks U(1) -> Z_q). Self-contained abelian analog of
// the SU(N) code, reusing the infrastructure (arbitrary-D Lattice, counter-based RNG, Stats).
//
// Links are phase ANGLES theta_mu(x) in R (u = e^{i theta}); gauge action beta sum(1-cos t_pl);
// charge-q scalar phi (1 complex component per site) with covariant hopping phi*_x u^q phi_{x+mu}.
// Gauge transform: theta_mu(x) -> theta_mu(x)+a(x)-a(x+mu), phi_x -> e^{i q a(x)} phi_x.
// All forces are finite-difference-validated in test_u1.cpp.
#include "core/geometry.hpp"
#include "core/rng.hpp"
#include "core/config.hpp"
#include <vector>
#include <cmath>
#include <complex>

namespace gh {
namespace u1 {

// ---- gauge sector (link angles theta[s*D+mu]) ----
template <int D>
Real plaq_angle(const std::vector<Real>& th, const Lattice<D>& lat, std::int64_t s, int mu, int nu) {
  const std::int64_t spmu = lat.neighbor_fwd(s, mu), spnu = lat.neighbor_fwd(s, nu);
  return th[s * D + mu] + th[spmu * D + nu] - th[spnu * D + mu] - th[s * D + nu];
}

template <int D>
Real avg_plaquette(const std::vector<Real>& th, const Lattice<D>& lat) {
  Real s = 0.0;
  #pragma omp parallel for schedule(static) reduction(+:s)
  for (std::int64_t x = 0; x < lat.vol; ++x)
    for (int mu = 0; mu < D; ++mu)
      for (int nu = mu + 1; nu < D; ++nu) s += std::cos(plaq_angle<D>(th, lat, x, mu, nu));
  return s / (static_cast<Real>(lat.n_plaq()));
}

template <int D>
Real gauge_action(const std::vector<Real>& th, const Lattice<D>& lat, Real beta) {
  Real s = 0.0;
  #pragma omp parallel for schedule(static) reduction(+:s)
  for (std::int64_t x = 0; x < lat.vol; ++x)
    for (int mu = 0; mu < D; ++mu)
      for (int nu = mu + 1; nu < D; ++nu) s += 1.0 - std::cos(plaq_angle<D>(th, lat, x, mu, nu));
  return beta * s;
}

// F_mu(x) = dS_g/dtheta_mu(x) = beta sum_{nu!=mu} [ sin(t_pl(x;mu,nu)) - sin(t_pl(x-nu;mu,nu)) ].
template <int D>
void add_gauge_force(const std::vector<Real>& th, const Lattice<D>& lat, Real beta, std::vector<Real>& F) {
  #pragma omp parallel for schedule(static)
  for (std::int64_t x = 0; x < lat.vol; ++x)
    for (int mu = 0; mu < D; ++mu) {
      Real f = 0.0;
      for (int nu = 0; nu < D; ++nu) {
        if (nu == mu) continue;
        const std::int64_t xmnu = lat.neighbor_bwd(x, nu);
        f += std::sin(plaq_angle<D>(th, lat, x, mu, nu)) - std::sin(plaq_angle<D>(th, lat, xmnu, mu, nu));
      }
      F[x * D + mu] += beta * f;
    }
}

// ---- charge-q scalar sector ----
// S_H = sum_x [ |phi|^2 + lambda(|phi|^2-1)^2 ] - kappa sum_{x,mu} 2 Re[ conj(phi_x) e^{i q theta} phi_{x+mu} ].
template <int D>
Real scalar_action(const std::vector<Complex>& phi, const std::vector<Real>& th, const Lattice<D>& lat,
                   int q, Real kappa, Real lambda) {
  Real S = 0.0;
  #pragma omp parallel for schedule(static) reduction(+:S)
  for (std::int64_t x = 0; x < lat.vol; ++x) {
    const Real n2 = std::norm(phi[x]);
    S += n2 + lambda * (n2 - 1.0) * (n2 - 1.0);
    for (int mu = 0; mu < D; ++mu) {
      const std::int64_t y = lat.neighbor_fwd(x, mu);
      const Complex ph = std::polar(1.0, q * th[x * D + mu]);  // u^q
      S -= 2.0 * kappa * (std::conj(phi[x]) * ph * phi[y]).real();
    }
  }
  return S;
}

// MD force on phi: F_phi = -2 dS_H/dphi^*
//   = -2 phi[1+2lambda(|phi|^2-1)] + 2 kappa sum_mu [ e^{i q t_mu(x)} phi_{x+mu} + e^{-i q t_mu(x-mu)} phi_{x-mu} ].
template <int D>
void scalar_force(const std::vector<Complex>& phi, const std::vector<Real>& th, const Lattice<D>& lat,
                  int q, Real kappa, Real lambda, std::vector<Complex>& Fphi) {
  #pragma omp parallel for schedule(static)
  for (std::int64_t x = 0; x < lat.vol; ++x) {
    const Real n2 = std::norm(phi[x]);
    Complex f = -2.0 * (1.0 + 2.0 * lambda * (n2 - 1.0)) * phi[x];
    for (int mu = 0; mu < D; ++mu) {
      const std::int64_t yf = lat.neighbor_fwd(x, mu), yb = lat.neighbor_bwd(x, mu);
      const Complex pf = std::polar(1.0, q * th[x * D + mu]);          // e^{i q t_mu(x)}
      const Complex pb = std::polar(1.0, -q * th[yb * D + mu]);        // e^{-i q t_mu(x-mu)}
      f += 2.0 * kappa * (pf * phi[yf] + pb * phi[yb]);
    }
    Fphi[x] = f;
  }
}

// Matter back-reaction on the links: dS_H/dtheta_mu(x) = 2 kappa q Im[ conj(phi_x) e^{i q theta} phi_{x+mu} ].
template <int D>
void add_matter_force(const std::vector<Complex>& phi, const std::vector<Real>& th, const Lattice<D>& lat,
                      int q, Real kappa, std::vector<Real>& F) {
  #pragma omp parallel for schedule(static)
  for (std::int64_t x = 0; x < lat.vol; ++x)
    for (int mu = 0; mu < D; ++mu) {
      const std::int64_t y = lat.neighbor_fwd(x, mu);
      const Complex ph = std::polar(1.0, q * th[x * D + mu]);
      F[x * D + mu] += 2.0 * kappa * q * (std::conj(phi[x]) * ph * phi[y]).imag();
    }
}

// ---- observables ----
template <int D>
Real higgs_length(const std::vector<Complex>& phi, const Lattice<D>& lat) {
  Real s = 0.0;
  #pragma omp parallel for schedule(static) reduction(+:s)
  for (std::int64_t x = 0; x < lat.vol; ++x) s += std::norm(phi[x]);
  return s / static_cast<Real>(lat.vol);
}

// Gauge-invariant charge-q hopping energy <Re conj(phi_x) e^{i q theta} phi_{x+mu}>.
template <int D>
Real link_energy(const std::vector<Complex>& phi, const std::vector<Real>& th, const Lattice<D>& lat, int q) {
  Real s = 0.0;
  #pragma omp parallel for schedule(static) reduction(+:s)
  for (std::int64_t x = 0; x < lat.vol; ++x)
    for (int mu = 0; mu < D; ++mu)
      s += (std::conj(phi[x]) * std::polar(1.0, q * th[x * D + mu]) * phi[lat.neighbor_fwd(x, mu)]).real();
  return s / (static_cast<Real>(lat.vol) * D);
}

// Charge-m rectangular Wilson loop <cos(m * theta_loop)>; a charge-m source is screened by
// the charge-q condensate iff q | m (the Z_q-neutral charges).
template <int D>
Real wilson_loop(const std::vector<Real>& th, const Lattice<D>& lat, int m, int Ra, int Rb) {
  Real sum = 0.0; std::int64_t count = 0;
  for (std::int64_t x = 0; x < lat.vol; ++x)
    for (int mu = 0; mu < D; ++mu)
      for (int nu = mu + 1; nu < D; ++nu) {
        Real loop = 0.0; std::int64_t s = x;
        for (int i = 0; i < Ra; ++i) { loop += th[s * D + mu]; s = lat.neighbor_fwd(s, mu); }
        for (int j = 0; j < Rb; ++j) { loop += th[s * D + nu]; s = lat.neighbor_fwd(s, nu); }
        for (int i = 0; i < Ra; ++i) { s = lat.neighbor_bwd(s, mu); loop -= th[s * D + mu]; }
        for (int j = 0; j < Rb; ++j) { s = lat.neighbor_bwd(s, nu); loop -= th[s * D + nu]; }
        sum += std::cos(m * loop); ++count;
      }
  return sum / static_cast<Real>(count);
}

// Charge-m Polyakov loop in direction tdir.
template <int D>
Real polyakov(const std::vector<Real>& th, const Lattice<D>& lat, int m, int tdir = D - 1) {
  const int Lt = lat.L[tdir];
  Real sum = 0.0; std::int64_t count = 0;
  std::array<int, D> xc{};
  for (std::int64_t s = 0; s < lat.vol; ++s) {
    lat.coords(s, xc); if (xc[tdir] != 0) continue;
    Real line = 0.0; std::int64_t cur = s;
    for (int t = 0; t < Lt; ++t) { line += th[cur * D + tdir]; cur = lat.neighbor_fwd(cur, tdir); }
    sum += std::cos(m * line); ++count;
  }
  return sum / static_cast<Real>(count);
}

enum class Integ { Leapfrog, Omelyan2MN };

// HMC for compact U(1) + charge-q Higgs.
template <int D>
struct U1HMC {
  Lattice<D> lat;
  std::vector<Real> th, p, Flink;       // link angles, momenta, force buffer  [vol*D]
  std::vector<Complex> phi, pi, Fphi;   // scalar, momentum, force buffer       [vol]
  Rng rng;
  Real beta = 1.0, kappa = 0.2, lambda = 0.5, tau = 1.0;
  int  q = 2, nmd = 20;
  int  n_scalar = 1;  // fast scalar sub-steps per coarse gauge step (1 == single timescale)
  Integ integ = Integ::Omelyan2MN;
  std::uint64_t traj_count = 0, accept_count = 0;
  Real last_dH = 0.0;

  U1HMC(const std::array<int, D>& ext, std::uint64_t seed)
      : lat(ext), th(static_cast<std::size_t>(lat.vol) * D, 0.0), p(th.size(), 0.0), Flink(th.size(), 0.0),
        phi(lat.vol, Complex(0, 0)), pi(lat.vol), Fphi(lat.vol), rng(seed) {}

  void hot(Real sigma = 1.0) {
    for (std::int64_t s = 0; s < lat.vol; ++s) {
      for (int mu = 0; mu < D; ++mu) th[s * D + mu] = sigma * rng.gauss(Rng::key(1, s, mu));
      phi[s] = Complex(rng.gauss(Rng::key(2, s, 0)), rng.gauss(Rng::key(2, s, 1))) * Complex(0.3, 0);
    }
  }
  void cold_phi(Real c = 1.0) { for (auto& z : phi) z = Complex(c, 0); }

  Real kinetic_link() const { Real e = 0; for (Real v : p) e += v * v; return 0.5 * e; }
  Real kinetic_scalar() const { Real e = 0; for (const auto& z : pi) e += std::norm(z); return 0.5 * e; }
  Real hamiltonian() {
    return kinetic_link() + kinetic_scalar()
         + gauge_action<D>(th, lat, beta) + scalar_action<D>(phi, th, lat, q, kappa, lambda);
  }

  // Combined single-timescale kick: gauge + matter -> p, scalar -> pi. Evaluates
  // the (slow) gauge force once, so it is tallied in kick_count_gauge for an
  // apples-to-apples gauge-force-evaluation count against the nested scheme.
  void kick(Real eps) {
    std::fill(Flink.begin(), Flink.end(), 0.0);
    add_gauge_force<D>(th, lat, beta, Flink);
    add_matter_force<D>(phi, th, lat, q, kappa, Flink);
    #pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < p.size(); ++i) p[i] -= eps * Flink[i];
    scalar_force<D>(phi, th, lat, q, kappa, lambda, Fphi);
    #pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < pi.size(); ++i) pi[i] += eps * Fphi[i];
    ++kick_count_gauge;
  }

  // ---- multi-timescale (Sexton-Weingarten) force split ----------------------
  // SLOW force = smooth gauge plaquette force only (on link momenta p). This is
  // the term whose evaluation we want to keep RARE on the coarse timescale.
  // FAST/STIFF force = the whole scalar sector: scalar_force on pi PLUS the matter
  // back-reaction add_matter_force on p. The high-kappa / high-q stiffness lives
  // entirely in this gauge-scalar coupling, so it is resolved on the fine grid.
  // kick_count_gauge counts SLOW (gauge) force evaluations for cost accounting.
  std::uint64_t kick_count_gauge = 0;

  // SLOW kick: gauge plaquette force only -> link momenta p.  p -= eps * F_gauge.
  void kick_gauge(Real eps) {
    std::fill(Flink.begin(), Flink.end(), 0.0);
    add_gauge_force<D>(th, lat, beta, Flink);
    #pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < p.size(); ++i) p[i] -= eps * Flink[i];
    ++kick_count_gauge;
  }

  // FAST kick: full scalar sector. matter back-reaction -> p, scalar EOM -> pi.
  void kick_fast(Real eps) {
    std::fill(Flink.begin(), Flink.end(), 0.0);
    add_matter_force<D>(phi, th, lat, q, kappa, Flink);
    #pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < p.size(); ++i) p[i] -= eps * Flink[i];
    scalar_force<D>(phi, th, lat, q, kappa, lambda, Fphi);
    #pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < pi.size(); ++i) pi[i] += eps * Fphi[i];
  }
  void drift(Real eps) {
    #pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < th.size(); ++i) th[i] += eps * p[i];
    #pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < phi.size(); ++i) phi[i] += eps * pi[i];
  }
  void md_evolve() {
    const Real eps = tau / nmd;
    if (integ == Integ::Leapfrog) {
      kick(0.5 * eps);
      for (int i = 0; i < nmd; ++i) { drift(eps); if (i != nmd - 1) kick(eps); }
      kick(0.5 * eps);
    } else {
      const Real lam = kOmelyanLambda;
      kick(lam * eps);
      for (int i = 0; i < nmd; ++i) { drift(0.5 * eps); kick((1 - 2 * lam) * eps); drift(0.5 * eps); kick((i != nmd - 1 ? 2 * lam : lam) * eps); }
    }
  }
  // ---- nested multi-timescale Sexton-Weingarten integrator -----------------
  // OUTER: Omelyan-2MN (or leapfrog) over the SLOW gauge force, nmd coarse steps.
  // INNER: each coarse "drift" interval H is replaced by a self-contained,
  // symmetric Omelyan-2MN of n_scalar fine sub-steps over the FAST scalar(+matter)
  // force together with the th/phi DRIFT. Because every inner integrator is itself
  // time-symmetric and the outer slow-kick schedule is symmetric, the whole map is
  // symmetric => reversible & symplectic for any n_scalar.
  //
  // The fast (matter+scalar) force depends on BOTH th and phi, and th is advanced
  // by the link-momentum drift; both fields therefore live in the inner drift, so
  // the gauge-scalar coupling is fully resolved on the fine grid while the slow
  // gauge plaquette force is touched only nmd times per trajectory.
  //
  // n_scalar == 1 special case: reproduces the existing single-timescale combined
  // Omelyan trajectory exactly (same operations, same coefficients).
  void md_evolve_mts() {
    if (n_scalar <= 1) { md_evolve(); return; }
    const Real eps = tau / nmd;                  // coarse (slow / gauge) step
    const Real lam = kOmelyanLambda;
    if (integ == Integ::Leapfrog) {
      kick_gauge(0.5 * eps);
      for (int i = 0; i < nmd; ++i) {
        inner_fast(eps);                         // fast Omelyan over the full coarse step
        if (i != nmd - 1) kick_gauge(eps);
      }
      kick_gauge(0.5 * eps);
    } else {
      kick_gauge(lam * eps);
      for (int i = 0; i < nmd; ++i) {
        inner_fast(0.5 * eps);
        kick_gauge((1.0 - 2.0 * lam) * eps);
        inner_fast(0.5 * eps);
        kick_gauge((i != nmd - 1 ? 2.0 * lam : lam) * eps);
      }
    }
  }

  // Symmetric Omelyan-2MN of n_scalar fine sub-steps over the FAST scalar(+matter)
  // force and the th/phi drift, integrating over a coarse sub-interval H.
  void inner_fast(Real H) {
    const Real lam = kOmelyanLambda;
    const int  m   = n_scalar;
    const Real d   = H / m;                      // fine step
    if (integ == Integ::Leapfrog) {
      kick_fast(0.5 * d);
      for (int j = 0; j < m; ++j) { drift(d); if (j != m - 1) kick_fast(d); }
      kick_fast(0.5 * d);
    } else {
      kick_fast(lam * d);
      for (int j = 0; j < m; ++j) {
        drift(0.5 * d);
        kick_fast((1.0 - 2.0 * lam) * d);
        drift(0.5 * d);
        kick_fast((j != m - 1 ? 2.0 * lam : lam) * d);
      }
    }
  }

  void refresh_momenta() {
    const std::uint64_t sg = Rng::key(0x51, traj_count), ss = Rng::key(0x52, traj_count);
    for (std::int64_t s = 0; s < lat.vol; ++s) {
      for (int mu = 0; mu < D; ++mu) p[s * D + mu] = rng.gauss(Rng::key(sg, s, mu));
      pi[s] = Complex(rng.gauss(Rng::key(ss, s, 0)), rng.gauss(Rng::key(ss, s, 1)));
    }
  }
  bool trajectory() {
    refresh_momenta();
    std::vector<Real> th0 = th; std::vector<Complex> phi0 = phi;
    const Real Hi = hamiltonian();
    if (n_scalar > 1) md_evolve_mts(); else md_evolve();
    const Real Hf = hamiltonian();
    last_dH = Hf - Hi; ++traj_count;
    const double r = rng.uniform(Rng::key(0x53, traj_count));
    bool acc = (last_dH <= 0.0) || (r < std::exp(-last_dH));
    if (acc) ++accept_count; else { th = th0; phi = phi0; }
    return acc;
  }
  double acceptance() const { return traj_count ? double(accept_count) / double(traj_count) : 0.0; }
};

}  // namespace u1
}  // namespace gh
