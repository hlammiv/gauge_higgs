#pragma once
// HMC engine (pure gauge for now): symplectic integrators (leapfrog, Omelyan 2MN),
// momentum heatbath, Metropolis accept/reject, reversibility diagnostics.
// The gauge+Higgs extension reuses this integrator structure; see hmc/integrator notes.
#include "action/gauge_wilson.hpp"
#include <cmath>

namespace gh {

enum class Integrator { Leapfrog, Omelyan2MN };

template <int D, int N>
struct GaugeHMC {
  Lattice<D>      lat;
  GaugeField<D, N> U;
  LinkMom<D, N>    P;
  LinkMom<D, N>    F;   // reusable force buffer
  Rng              rng;
  // parameters
  Real    beta = 2.3;
  Real    tau  = 1.0;          // trajectory length
  int     nmd  = 20;           // MD steps per trajectory
  Integrator integ = Integrator::Omelyan2MN;
  bool    reunit_each_traj = true;
  // statistics
  std::uint64_t traj_count = 0, accept_count = 0;
  Real last_dH = 0.0;

  GaugeHMC(const std::array<int, D>& extents, std::uint64_t seed)
      : lat(extents), U(lat), P(lat), F(lat), rng(seed) {}

  Real hamiltonian() const { return P.kinetic() + gauge_action<D, N>(U, beta); }

  // P^a -= eps * F^a, with F freshly computed from the current U.
  void kick(Real eps) {
    F.zero();
    add_gauge_force<D, N>(U, beta, F);
    #pragma omp parallel for schedule(static)
    for (std::int64_t i = 0; i < static_cast<std::int64_t>(P.p.size()); ++i)
      for (int a = 0; a < n_gen<N>(); ++a) P.p[i][a] -= eps * F.p[i][a];
  }

  // U_mu(x) <- exp(i eps P_mu(x)) U_mu(x).
  void drift(Real eps) {
    #pragma omp parallel for schedule(static)
    for (std::int64_t s = 0; s < lat.vol; ++s)
      for (int mu = 0; mu < D; ++mu) {
        Cmat<N> Pm = alg_to_mat<N>(P(s, mu));
        Pm *= Complex(eps, 0.0);
        U(s, mu) = expi<N>(Pm) * U(s, mu);
      }
  }

  void md_evolve() {
    const Real eps = tau / nmd;
    if (integ == Integrator::Leapfrog) {
      kick(0.5 * eps);
      for (int i = 0; i < nmd; ++i) { drift(eps); if (i != nmd - 1) kick(eps); }
      kick(0.5 * eps);
    } else {  // Omelyan 2MN, merged kicks
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

  // One HMC trajectory. Returns true if accepted. `stream` seeds the momentum heatbath.
  bool trajectory() {
    const std::uint64_t stream = Rng::key(0xA1, traj_count);
    P.refresh(rng, stream);
    std::vector<Cmat<N>> U_save = U.u;  // for reject / reversibility

    const Real H_i = hamiltonian();
    md_evolve();
    if (reunit_each_traj) U.reunitarize_all();
    const Real H_f = hamiltonian();

    last_dH = H_f - H_i;
    ++traj_count;
    // Metropolis: accept with prob min(1, exp(-dH)).
    const double r = rng.uniform(Rng::key(0xB2, traj_count));
    bool accept = (last_dH <= 0.0) || (r < std::exp(-last_dH));
    if (accept) { ++accept_count; }
    else        { U.u = U_save; }
    return accept;
  }

  double acceptance() const {
    return traj_count ? double(accept_count) / double(traj_count) : 0.0;
  }
};

}  // namespace gh
