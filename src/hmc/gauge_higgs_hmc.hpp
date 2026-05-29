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
  // stats
  std::uint64_t traj_count = 0, accept_count = 0;
  Real last_dH = 0.0;

  GaugeHiggsHMC(const std::array<int, D>& extents, const Representation<N>& R, std::uint64_t seed)
      : lat(extents), U(lat), P(lat), Flink(lat),
        phi(lat, R.d), pi(lat, R.d), Fphi(lat, R.d), rep(&R), rng(seed) {}

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
  }

  void drift(Real eps) {
    #pragma omp parallel for schedule(static)
    for (std::int64_t s = 0; s < lat.vol; ++s)
      for (int mu = 0; mu < D; ++mu) {
        Cmat<N> Pm = alg_to_mat<N>(P(s, mu));
        Pm *= Complex(eps, 0.0);
        U(s, mu) = expi<N>(Pm) * U(s, mu);
      }
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
