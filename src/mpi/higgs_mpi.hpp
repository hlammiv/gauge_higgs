#pragma once
// MPI domain-decomposed scalar (Higgs) sector + combined gauge+Higgs HMC, for an
// arbitrary irrep. Reuses DistLattice + halo_exchange_raw (the scalar field needs only
// face halos, a subset of the gauge halo) and the lattice-agnostic Representation
// interface. Validated like the gauge slice: decomposition-independence (np 1/2/4 give
// identical plaquette, L_link, L_phi via the global-index RNG).
#include "mpi/gauge_mpi.hpp"
#include "rep/representation.hpp"
#include "action/scalar_higgs.hpp"  // OnsitePotential, QuarticPotential

namespace gh {

// Per-site complex d-vector field over the distributed lattice (Higgs phi / momentum pi).
template <int D>
struct DistScalar {
  const DistLattice<D>* lat = nullptr;
  int d = 0;
  std::vector<Complex> data;  // [store*d]
  DistScalar(const DistLattice<D>& L, int dim)
      : lat(&L), d(dim), data(static_cast<std::size_t>(L.store) * dim, Complex(0, 0)) {}

  DVec get(long s) const { DVec v(d); const long o = s * d; for (int k = 0; k < d; ++k) v(k) = data[o + k]; return v; }
  void set(long s, const DVec& v) { const long o = s * d; for (int k = 0; k < d; ++k) data[o + k] = v(k); }

  void cold(Real c = 1.0) {
    std::fill(data.begin(), data.end(), Complex(0, 0));
    for (long i = 0; i < lat->vint; ++i) data[lat->interior[i] * d] = Complex(c, 0);
  }
  void gaussian(const Rng& rng, std::uint64_t stream, bool real_rep, Real sigma = 1.0) {
    for (long i = 0; i < lat->vint; ++i) {
      const long s = lat->interior[i], g = lat->gidx[i], o = s * d;
      for (int k = 0; k < d; ++k) {
        const Real re = sigma * rng.gauss(Rng::key(stream, g, 2 * k));
        const Real im = real_rep ? 0.0 : sigma * rng.gauss(Rng::key(stream, g, 2 * k + 1));
        data[o + k] = Complex(re, im);
      }
    }
  }
  Real kinetic() const {
    Real e = 0.0;
    for (long i = 0; i < lat->vint; ++i) { const long o = lat->interior[i] * d; for (int k = 0; k < d; ++k) e += std::norm(data[o + k]); }
    return 0.5 * mpi_allreduce_sum(e, lat->comm);
  }
};

template <int D>
void halo_exchange(DistScalar<D>& phi) { halo_exchange_raw<D>(*phi.lat, phi.data.data(), phi.d); }

// Global scalar action (halos must be current). On-site potential is pluggable.
template <int D, int N>
Real scalar_action_d(const DistScalar<D>& phi, const DistGauge<D, N>& U,
                     const Representation<N>& rep, Real kappa, const OnsitePotential<N>& pot) {
  const DistLattice<D>& L = *U.lat;
  Real local = 0.0;
  for (long i = 0; i < L.vint; ++i) {
    const long s = L.interior[i];
    const DVec px = phi.get(s);
    local += pot.value(px);
    for (int mu = 0; mu < D; ++mu) {
      const DVec Dpy = rep.rotate(U(s, mu), phi.get(s + L.stride[mu]));
      local -= 2.0 * kappa * dot(px, Dpy).real();
    }
  }
  return mpi_allreduce_sum(local, L.comm);
}
template <int D, int N>
Real scalar_action_d(const DistScalar<D>& phi, const DistGauge<D, N>& U,
                     const Representation<N>& rep, Real kappa, Real lambda) {
  QuarticPotential<N> q(lambda);
  return scalar_action_d<D, N>(phi, U, rep, kappa, q);
}

template <int D, int N>
void scalar_force_d(const DistScalar<D>& phi, const DistGauge<D, N>& U,
                    const Representation<N>& rep, Real kappa, const OnsitePotential<N>& pot,
                    DistScalar<D>& Fphi) {
  const DistLattice<D>& L = *U.lat;
  for (long i = 0; i < L.vint; ++i) {
    const long s = L.interior[i];
    const DVec px = phi.get(s);
    const DVec g = pot.dV_dphibar(px);
    DVec F(rep.d);
    for (int k = 0; k < rep.d; ++k) F(k) = -2.0 * g(k);
    for (int mu = 0; mu < D; ++mu) {
      const long yf = s + L.stride[mu], yb = s - L.stride[mu];
      const DVec Dpf = rep.rotate(U(s, mu), phi.get(yf));
      const DVec Dpb = rep.rotate_dag(U(yb, mu), phi.get(yb));
      for (int k = 0; k < rep.d; ++k) F(k) += 2.0 * kappa * (Dpf(k) + Dpb(k));
    }
    Fphi.set(s, F);
  }
}
template <int D, int N>
void scalar_force_d(const DistScalar<D>& phi, const DistGauge<D, N>& U,
                    const Representation<N>& rep, Real kappa, Real lambda, DistScalar<D>& Fphi) {
  QuarticPotential<N> q(lambda);
  scalar_force_d<D, N>(phi, U, rep, kappa, q, Fphi);
}

template <int D, int N>
void add_gauge_force_d(const DistGauge<D, N>& U, Real beta, std::vector<AlgVec<N>>& F) {
  const DistLattice<D>& L = *U.lat;
  const auto& gen = generators<N>();
  const Real c = beta / N;
  for (long i = 0; i < L.vint; ++i) {
    const long s = L.interior[i];
    for (int mu = 0; mu < D; ++mu) {
      const Cmat<N> Omega = U(s, mu) * staple_d<D, N>(U, s, mu);
      AlgVec<N>& Fa = F[s * D + mu];
      for (int a = 0; a < n_gen<N>(); ++a) Fa[a] += c * trProd(gen.T[a], Omega).imag();
    }
  }
}

template <int D, int N>
void add_matter_link_force_d(const DistScalar<D>& phi, const DistGauge<D, N>& U,
                             const Representation<N>& rep, Real kappa, std::vector<AlgVec<N>>& F) {
  const DistLattice<D>& L = *U.lat;
  for (long i = 0; i < L.vint; ++i) {
    const long s = L.interior[i];
    for (int mu = 0; mu < D; ++mu) {
      const AlgVec<N> g = rep.hop_link_g(U(s, mu), phi.get(s), phi.get(s + L.stride[mu]));
      AlgVec<N>& Fa = F[s * D + mu];
      for (int a = 0; a < n_gen<N>(); ++a) Fa[a] += -kappa * g[a];
    }
  }
}

template <int D, int N>
Real link_energy_d(const DistScalar<D>& phi, const DistGauge<D, N>& U, const Representation<N>& rep) {
  const DistLattice<D>& L = *U.lat;
  Real local = 0.0;
  for (long i = 0; i < L.vint; ++i) {
    const long s = L.interior[i];
    for (int mu = 0; mu < D; ++mu) local += dot(phi.get(s), rep.rotate(U(s, mu), phi.get(s + L.stride[mu]))).real();
  }
  return mpi_allreduce_sum(local, L.comm) / static_cast<Real>(L.n_link_global());
}

template <int D>
Real higgs_length_d(const DistScalar<D>& phi) {
  const DistLattice<D>& L = *phi.lat;
  Real local = 0.0;
  for (long i = 0; i < L.vint; ++i) local += phi.get(L.interior[i]).norm2();
  long gv = 1; for (int mu = 0; mu < D; ++mu) gv *= L.G[mu];
  return mpi_allreduce_sum(local, L.comm) / static_cast<Real>(gv);
}

// Combined gauge + Higgs HMC over the distributed lattice (arbitrary irrep).
template <int D, int N>
struct DistGaugeHiggsHMC {
  DistLattice<D> lat;
  DistGauge<D, N> U;
  std::vector<AlgVec<N>> P, Flink;
  DistScalar<D> phi, pi, Fphi;
  const Representation<N>* rep = nullptr;
  Rng rng;
  Real beta = 2.3, kappa = 0.1, lambda = 0.5, tau = 1.0;
  const OnsitePotential<N>* potential = nullptr;  // null => lambda quartic
  int  nmd = 20;
  MpiIntegrator integ = MpiIntegrator::Omelyan2MN;
  std::uint64_t traj_count = 0, accept_count = 0;
  Real last_dH = 0.0;

  DistGaugeHiggsHMC(const std::array<int, D>& G, const std::array<int, D>& pg,
                    const Representation<N>& R, std::uint64_t seed)
      : lat(G, pg), U(lat), P(static_cast<std::size_t>(lat.store) * D),
        Flink(static_cast<std::size_t>(lat.store) * D),
        phi(lat, R.d), pi(lat, R.d), Fphi(lat, R.d), rep(&R), rng(seed) {}

  Real link_kinetic() const {
    Real e = 0.0;
    for (long i = 0; i < lat.vint; ++i) { const long s = lat.interior[i]; for (int mu = 0; mu < D; ++mu) { const AlgVec<N>& v = P[s * D + mu]; for (int a = 0; a < n_gen<N>(); ++a) e += v[a] * v[a]; } }
    return 0.5 * mpi_allreduce_sum(e, lat.comm);
  }
  Real hamiltonian() {
    halo_exchange(U); halo_exchange(phi);
    return link_kinetic() + pi.kinetic() + gauge_action_d<D, N>(U, beta)
         + (potential ? scalar_action_d<D, N>(phi, U, *rep, kappa, *potential)
                      : scalar_action_d<D, N>(phi, U, *rep, kappa, lambda));
  }

  void kick(Real eps) {
    halo_exchange(U); halo_exchange(phi);
    std::fill(Flink.begin(), Flink.end(), AlgVec<N>{});
    add_gauge_force_d<D, N>(U, beta, Flink);
    add_matter_link_force_d<D, N>(phi, U, *rep, kappa, Flink);
    for (long i = 0; i < lat.vint; ++i) { const long s = lat.interior[i]; for (int mu = 0; mu < D; ++mu) { AlgVec<N>& p = P[s * D + mu]; const AlgVec<N>& f = Flink[s * D + mu]; for (int a = 0; a < n_gen<N>(); ++a) p[a] -= eps * f[a]; } }
    if (potential) scalar_force_d<D, N>(phi, U, *rep, kappa, *potential, Fphi);
    else           scalar_force_d<D, N>(phi, U, *rep, kappa, lambda, Fphi);
    for (long i = 0; i < lat.vint; ++i) { const long o = lat.interior[i] * rep->d; for (int k = 0; k < rep->d; ++k) pi.data[o + k] += eps * Fphi.data[o + k]; }
  }
  void drift(Real eps) {
    for (long i = 0; i < lat.vint; ++i) { const long s = lat.interior[i]; for (int mu = 0; mu < D; ++mu) { Cmat<N> Pm = alg_to_mat<N>(P[s * D + mu]); Pm *= Complex(eps, 0.0); U(s, mu) = expi<N>(Pm) * U(s, mu); } }
    for (long i = 0; i < lat.vint; ++i) { const long o = lat.interior[i] * rep->d; for (int k = 0; k < rep->d; ++k) phi.data[o + k] += eps * pi.data[o + k]; }
  }
  void md_evolve() {
    const Real eps = tau / nmd;
    if (integ == MpiIntegrator::Leapfrog) {
      kick(0.5 * eps);
      for (int i = 0; i < nmd; ++i) { drift(eps); if (i != nmd - 1) kick(eps); }
      kick(0.5 * eps);
    } else {
      const Real lam = kOmelyanLambda;
      kick(lam * eps);
      for (int i = 0; i < nmd; ++i) { drift(0.5 * eps); kick((1.0 - 2.0 * lam) * eps); drift(0.5 * eps); kick((i != nmd - 1 ? 2.0 * lam : lam) * eps); }
    }
  }
  bool trajectory() {
    const std::uint64_t sg = Rng::key(0xC3, traj_count), ss = Rng::key(0xD4, traj_count);
    for (long i = 0; i < lat.vint; ++i) { const long s = lat.interior[i], g = lat.gidx[i]; for (int mu = 0; mu < D; ++mu) { AlgVec<N>& v = P[s * D + mu]; for (int a = 0; a < n_gen<N>(); ++a) v[a] = rng.gauss(Rng::key(sg, g, mu, a)); } }
    pi.gaussian(rng, ss, rep->real);
    std::vector<Cmat<N>> Usave = U.u;
    std::vector<Complex> psave = phi.data;
    const Real Hi = hamiltonian();
    md_evolve();
    U.reunitarize_interior();
    const Real Hf = hamiltonian();
    last_dH = Hf - Hi; ++traj_count;
    const double r = rng.uniform(Rng::key(0xE5, traj_count));
    bool acc = (last_dH <= 0.0) || (r < std::exp(-last_dH));
    if (acc) ++accept_count; else { U.u = Usave; phi.data = psave; }
    return acc;
  }
  double acceptance() const { return traj_count ? double(accept_count) / double(traj_count) : 0.0; }
};

}  // namespace gh
