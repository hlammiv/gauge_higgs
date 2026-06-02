#pragma once
// Combined gauge + Higgs HMC for an arbitrary irrep R. Reuses the symplectic
// integrator structure of the pure-gauge engine; the kick now updates BOTH the
// link momentum (gauge force + matter back-reaction) and the scalar momentum.
// Nested multiple-timescale (gauge fine / Higgs coarse) is a planned extension;
// this is single-timescale leapfrog / Omelyan 2MN. See theory_notes §2, §4.
#include "action/gauge_wilson.hpp"
#include "action/scalar_higgs.hpp"
#include "hmc/hmc.hpp"  // for Integrator enum
#include <cmath>

namespace gh {

template <int D, int N>
struct GaugeHiggsHMC {
  Lattice<D>       lat;
  GaugeField<D, N> U;
  LinkMom<D, N>    P, Flink;      // link momentum + reusable link-force buffer
  CVecField<D>     phi, pi, Fphi; // scalar field, momentum, reusable scalar-force buffer
  const Representation<N>* rep = nullptr;
  Rng rng;
  // couplings
  Real beta = 2.3, kappa = 0.1, lambda = 0.5, tau = 1.0;
  // Optional pluggable on-site potential (e.g. the multi-invariant potential that locks
  // a discrete subgroup H). If null, the simple lambda(phi^dag phi-1)^2 quartic is used.
  const OnsitePotential<N>* potential = nullptr;
  int  nmd = 20;
  Integrator integ = Integrator::Omelyan2MN;
  bool reunit_each_traj = true;
  bool frozen_phi = false;   // |phi_x|=1 frozen-length (sphere) scalar HMC: geodesic drift +
                             // tangent-projected momentum, removing the radial |phi| mode so
                             // kappa probes the gauge-Higgs transition, not radial condensation.
  // stats
  std::uint64_t traj_count = 0, accept_count = 0;
  Real last_dH = 0.0;

  GaugeHiggsHMC(const std::array<int, D>& extents, const Representation<N>& R, std::uint64_t seed)
      : lat(extents), U(lat), P(lat), Flink(lat),
        phi(lat, R.d), pi(lat, R.d), Fphi(lat, R.d), rep(&R), rng(seed) {}

  // ---- frozen-length |phi_x|=1 helpers (per-site sphere; real rep -> S^{d-1}, complex -> S^{2d-1}) ----
  // Project pi onto the tangent space at phi (remove the radial part Re(phi^dag pi) phi), per site.
  void project_pi_tangent() {
    const int d = phi.d;
    #pragma omp parallel for schedule(static)
    for (std::int64_t s = 0; s < lat.vol; ++s) {
      const std::size_t o = static_cast<std::size_t>(s) * d;
      Real rp = 0.0;
      for (int k = 0; k < d; ++k) rp += (std::conj(phi.data[o + k]) * pi.data[o + k]).real();
      for (int k = 0; k < d; ++k) pi.data[o + k] -= Complex(rp, 0.0) * phi.data[o + k];
    }
  }
  // Normalize each site to |phi_x|=1 (call once at init for a frozen-length run).
  void normalize_phi() {
    const int d = phi.d;
    #pragma omp parallel for schedule(static)
    for (std::int64_t s = 0; s < lat.vol; ++s) {
      const std::size_t o = static_cast<std::size_t>(s) * d;
      Real n2 = 0.0; for (int k = 0; k < d; ++k) n2 += std::norm(phi.data[o + k]);
      if (n2 > 0.0) { const Real inv = 1.0 / std::sqrt(n2);
        for (int k = 0; k < d; ++k) phi.data[o + k] *= Complex(inv, 0.0); }
    }
  }
  // Exact great-circle (geodesic) drift of (phi,pi) on the unit sphere by time eps, per site:
  //   phi' = cos(r eps) phi + sin(r eps) pi/r ,  pi' = -r sin(r eps) phi + cos(r eps) pi ,  r=|pi|.
  // This is the exact constrained free flow: keeps |phi|=1 and pi tangent to machine precision,
  // and is time-reversible -> a valid symplectic scalar drift for the frozen-length HMC.
  void geodesic_drift(Real eps) {
    const int d = phi.d;
    #pragma omp parallel for schedule(static)
    for (std::int64_t s = 0; s < lat.vol; ++s) {
      const std::size_t o = static_cast<std::size_t>(s) * d;
      Real r2 = 0.0; for (int k = 0; k < d; ++k) r2 += std::norm(pi.data[o + k]);
      const Real r = std::sqrt(r2);
      if (r < 1e-300) continue;
      const Real c = std::cos(r * eps), sn = std::sin(r * eps);
      for (int k = 0; k < d; ++k) {
        const Complex f = phi.data[o + k], p = pi.data[o + k];
        phi.data[o + k] = Complex(c, 0.0) * f + Complex(sn / r, 0.0) * p;
        pi.data[o + k]  = Complex(-r * sn, 0.0) * f + Complex(c, 0.0) * p;
      }
    }
  }

  Real hamiltonian() const {
    return P.kinetic() + pi.kinetic()
         + gauge_action<D, N>(U, beta)
         + (potential ? scalar_action<D, N>(phi, U, *rep, kappa, *potential)
                      : scalar_action<D, N>(phi, U, *rep, kappa, lambda));
  }

  void kick(Real eps) {
    // link force: gauge staple + matter staple ; P -= eps F
    Flink.zero();
    add_gauge_force<D, N>(U, beta, Flink);
    add_matter_link_force<D, N>(phi, U, *rep, kappa, Flink);
    #pragma omp parallel for schedule(static)
    for (std::int64_t i = 0; i < static_cast<std::int64_t>(P.p.size()); ++i)
      for (int a = 0; a < n_gen<N>(); ++a) P.p[i][a] -= eps * Flink.p[i][a];
    // scalar force ; pi += eps F_phi
    if (potential) scalar_force<D, N>(phi, U, *rep, kappa, *potential, Fphi);
    else           scalar_force<D, N>(phi, U, *rep, kappa, lambda, Fphi);
    #pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < pi.data.size(); ++i) pi.data[i] += eps * Fphi.data[i];
    if (frozen_phi) project_pi_tangent();   // drop the radial force component (|phi|=1 constraint)
  }

  void drift(Real eps) {
    #pragma omp parallel for schedule(static)
    for (std::int64_t s = 0; s < lat.vol; ++s)
      for (int mu = 0; mu < D; ++mu) {
        Cmat<N> Pm = alg_to_mat<N>(P(s, mu));
        Pm *= Complex(eps, 0.0);
        U(s, mu) = expi<N>(Pm) * U(s, mu);
      }
    if (frozen_phi) { geodesic_drift(eps); return; }   // |phi|=1: great-circle drift on the sphere
    #pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < phi.data.size(); ++i) phi.data[i] += eps * pi.data[i];
  }

  void md_evolve() {
    const Real eps = tau / nmd;
    if (integ == Integrator::Leapfrog) {
      kick(0.5 * eps);
      for (int i = 0; i < nmd; ++i) { drift(eps); if (i != nmd - 1) kick(eps); }
      kick(0.5 * eps);
    } else {
      const Real lam = kOmelyanLambda;
      kick(lam * eps);
      for (int i = 0; i < nmd; ++i) {
        drift(0.5 * eps);
        kick((1.0 - 2.0 * lam) * eps);
        drift(0.5 * eps);
        kick((i != nmd - 1 ? 2.0 * lam : lam) * eps);
      }
    }
  }

  void refresh_momenta() {
    const std::uint64_t sg = Rng::key(0xC3, traj_count);
    const std::uint64_t ss = Rng::key(0xD4, traj_count);
    P.refresh(rng, sg);
    pi.gaussian(rng, ss, rep->real);
    if (frozen_phi) project_pi_tangent();   // Gaussian momentum restricted to the tangent space
  }

  bool trajectory() {
    refresh_momenta();
    std::vector<Cmat<N>> U_save  = U.u;
    std::vector<Complex> phi_save = phi.data;

    const Real H_i = hamiltonian();
    md_evolve();
    if (reunit_each_traj) U.reunitarize_all();
    const Real H_f = hamiltonian();

    last_dH = H_f - H_i;
    ++traj_count;
    const double r = rng.uniform(Rng::key(0xE5, traj_count));
    bool accept = (last_dH <= 0.0) || (r < std::exp(-last_dH));
    if (accept) ++accept_count;
    else { U.u = U_save; phi.data = phi_save; }
    return accept;
  }

  double acceptance() const { return traj_count ? double(accept_count) / double(traj_count) : 0.0; }
};

}  // namespace gh
