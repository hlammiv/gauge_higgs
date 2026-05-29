#pragma once
// MPI domain-decomposition layer (pure-gauge slice). Additive: reuses the math
// headers (linalg/algebra/sun/rng); the serial code is untouched. The lattice is
// split over a Cartesian process grid; each rank stores its local block padded by
// a width-1 halo (faces + edges + corners, filled by sequential directional
// exchange). The counter-based RNG is keyed on the GLOBAL site index, so a
// decomposed run reproduces the single-process result -> cheap correctness test.
//
// This proves the hard infrastructure (halo exchange, global reductions, global-RNG,
// decomposition-independent HMC). The scalar/arbitrary-irrep sector reuses the same
// DistLattice + halo_exchange machinery (the scalar field needs only face halos).
#include <mpi.h>
#include "group/sun.hpp"
#include "core/rng.hpp"
#include <array>
#include <vector>

namespace gh {

inline double mpi_allreduce_sum(double x, MPI_Comm comm) {
  double r = 0.0;
  MPI_Allreduce(&x, &r, 1, MPI_DOUBLE, MPI_SUM, comm);
  return r;
}

template <int D>
struct DistLattice {
  std::array<int, D>   G{};       // global extents
  std::array<int, D>   P{};       // process grid (ranks per dim)
  std::array<int, D>   pc{};      // this rank's process-grid coords
  std::array<int, D>   Lloc{};    // local interior extents
  std::array<int, D>   Lpad{};    // padded extents = Lloc + 2
  std::array<long, D>  stride{};  // padded strides (dim 0 fastest)
  std::array<long, D>  Gstride{}; // global serial strides (matches serial Lattice numbering)
  long store = 0, vint = 0;
  int  rank = 0, nranks = 1;
  MPI_Comm comm = MPI_COMM_NULL;
  std::array<int, D> nbr_plus{}, nbr_minus{};
  std::vector<long> interior;     // store indices of interior cells
  std::vector<long> gidx;         // global serial index of each interior cell (RNG key)

  DistLattice(const std::array<int, D>& global, const std::array<int, D>& procgrid)
      : G(global), P(procgrid) {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nranks);
    std::array<int, D> periods{}; periods.fill(1);
    MPI_Cart_create(MPI_COMM_WORLD, D, P.data(), periods.data(), 0, &comm);
    MPI_Cart_coords(comm, rank, D, pc.data());
    for (int mu = 0; mu < D; ++mu) { Lloc[mu] = G[mu] / P[mu]; Lpad[mu] = Lloc[mu] + 2; }
    store = 1; for (int mu = 0; mu < D; ++mu) { stride[mu] = store; store *= Lpad[mu]; }
    Gstride[0] = 1; for (int mu = 1; mu < D; ++mu) Gstride[mu] = Gstride[mu - 1] * G[mu - 1];
    vint = 1; for (int mu = 0; mu < D; ++mu) vint *= Lloc[mu];
    for (int mu = 0; mu < D; ++mu) { int s, d; MPI_Cart_shift(comm, mu, 1, &s, &d); nbr_minus[mu] = s; nbr_plus[mu] = d; }
    // enumerate interior cells (mixed radix over Lloc, dim 0 fastest)
    interior.reserve(vint); gidx.reserve(vint);
    std::array<int, D> ic{};
    for (long n = 0; n < vint; ++n) {
      long t = n;
      for (int mu = 0; mu < D; ++mu) { ic[mu] = static_cast<int>(t % Lloc[mu]); t /= Lloc[mu]; }
      long sidx = 0, g = 0;
      for (int mu = 0; mu < D; ++mu) {
        sidx += static_cast<long>(ic[mu] + 1) * stride[mu];
        g    += static_cast<long>(pc[mu] * Lloc[mu] + ic[mu]) * Gstride[mu];
      }
      interior.push_back(sidx); gidx.push_back(g);
    }
  }

  long n_plaq_global() const { long gv = 1; for (int mu = 0; mu < D; ++mu) gv *= G[mu]; return gv * (static_cast<long>(D) * (D - 1) / 2); }
  long n_link_global() const { long gv = 1; for (int mu = 0; mu < D; ++mu) gv *= G[mu]; return gv * D; }
};

template <int D, int N>
struct DistGauge {
  const DistLattice<D>* lat = nullptr;
  std::vector<Cmat<N>> u;  // [store*D + mu]
  explicit DistGauge(const DistLattice<D>& L) : lat(&L), u(static_cast<std::size_t>(L.store) * D) {}
  Cmat<N>&       operator()(long s, int mu)       { return u[s * D + mu]; }
  const Cmat<N>& operator()(long s, int mu) const { return u[s * D + mu]; }

  void hot(const Rng& rng, Real sigma = 0.8) {
    for (long i = 0; i < lat->vint; ++i) {
      const long s = lat->interior[i], g = lat->gidx[i];
      for (int mu = 0; mu < D; ++mu) {
        AlgVec<N> v{};
        for (int a = 0; a < n_gen<N>(); ++a) v[a] = sigma * rng.gauss(Rng::key(1, g, mu, a));
        (*this)(s, mu) = expi<N>(alg_to_mat<N>(v));
      }
    }
  }
  void reunitarize_interior() {
    for (long i = 0; i < lat->vint; ++i)
      for (int mu = 0; mu < D; ++mu) reunitarize<N>((*this)(lat->interior[i], mu));
  }
};

// Generic width-1 halo exchange (faces + edges + corners) via sequential directional
// swaps, for any field laid out as `ncomp` contiguous Complex per cell (cell s occupies
// data[s*ncomp .. s*ncomp+ncomp)). After all D directions, diagonals are filled. Used by
// both the gauge field (ncomp = D*N*N) and the scalar field (ncomp = d). The scalar only
// needs faces, but filling edges/corners too is harmless.
template <int D>
void halo_exchange_raw(const DistLattice<D>& L, Complex* data, long ncomp) {
  for (int mu = 0; mu < D; ++mu) {
    std::array<int, D> lo{}, hi{};
    long ncells = 1;
    for (int nu = 0; nu < D; ++nu) {
      if (nu == mu) { lo[nu] = hi[nu] = 0; continue; }
      if (nu < mu) { lo[nu] = 0; hi[nu] = L.Lpad[nu] - 1; }   // already-filled halos
      else         { lo[nu] = 1; hi[nu] = L.Lloc[nu]; }       // interior only
      ncells *= (hi[nu] - lo[nu] + 1);
    }
    std::vector<std::array<int, D>> tc; tc.reserve(ncells);
    std::array<int, D> c = lo;
    for (long k = 0; k < ncells; ++k) {
      tc.push_back(c);
      for (int nu = 0; nu < D; ++nu) { if (nu == mu) continue; if (++c[nu] <= hi[nu]) break; c[nu] = lo[nu]; }
    }
    std::vector<Complex> sbuf(ncells * ncomp), rbuf(ncells * ncomp);
    auto cell_index = [&](long k, int layer) {
      std::array<int, D> cc = tc[k]; cc[mu] = layer;
      long s = 0; for (int nu = 0; nu < D; ++nu) s += static_cast<long>(cc[nu]) * L.stride[nu];
      return s;
    };
    auto pack   = [&](int layer) { for (long k = 0; k < ncells; ++k) { long s = cell_index(k, layer); for (long e = 0; e < ncomp; ++e) sbuf[k * ncomp + e] = data[s * ncomp + e]; } };
    auto unpack = [&](int layer) { for (long k = 0; k < ncells; ++k) { long s = cell_index(k, layer); for (long e = 0; e < ncomp; ++e) data[s * ncomp + e] = rbuf[k * ncomp + e]; } };
    const int cnt = static_cast<int>(ncells * ncomp * 2);  // complex -> 2 doubles
    pack(L.Lloc[mu]);  // A: send last interior layer to +mu, recv lower halo from -mu
    MPI_Sendrecv(reinterpret_cast<double*>(sbuf.data()), cnt, MPI_DOUBLE, L.nbr_plus[mu], 100 + mu,
                 reinterpret_cast<double*>(rbuf.data()), cnt, MPI_DOUBLE, L.nbr_minus[mu], 100 + mu,
                 L.comm, MPI_STATUS_IGNORE);
    unpack(0);
    pack(1);           // B: send first interior layer to -mu, recv upper halo from +mu
    MPI_Sendrecv(reinterpret_cast<double*>(sbuf.data()), cnt, MPI_DOUBLE, L.nbr_minus[mu], 200 + mu,
                 reinterpret_cast<double*>(rbuf.data()), cnt, MPI_DOUBLE, L.nbr_plus[mu], 200 + mu,
                 L.comm, MPI_STATUS_IGNORE);
    unpack(L.Lloc[mu] + 1);
  }
}

// Gauge field halo: all D link matrices per cell are contiguous (vector<Cmat<N>>).
template <int D, int N>
void halo_exchange(DistGauge<D, N>& U) {
  halo_exchange_raw<D>(*U.lat, reinterpret_cast<Complex*>(U.u.data()), static_cast<long>(D) * N * N);
}

// Staple at interior link (s, mu) via padded-stride arithmetic (halos must be current).
template <int D, int N>
Cmat<N> staple_d(const DistGauge<D, N>& U, long s, int mu) {
  const DistLattice<D>& L = *U.lat;
  Cmat<N> sig{};
  const long s_pmu = s + L.stride[mu];
  for (int nu = 0; nu < D; ++nu) {
    if (nu == mu) continue;
    const long s_pnu = s + L.stride[nu], s_mnu = s - L.stride[nu], s_pmu_mnu = s_pmu - L.stride[nu];
    sig += U(s_pmu, nu) * U(s_pnu, mu).dagger() * U(s, nu).dagger();
    sig += U(s_pmu_mnu, nu).dagger() * U(s_mnu, mu).dagger() * U(s_mnu, nu);
  }
  return sig;
}

// Average plaquette (global, halos must be current).
template <int D, int N>
Real avg_plaquette_d(const DistGauge<D, N>& U) {
  const DistLattice<D>& L = *U.lat;
  Real local = 0.0;
  for (long i = 0; i < L.vint; ++i) {
    const long s = L.interior[i];
    for (int mu = 0; mu < D; ++mu) {
      const long s_pmu = s + L.stride[mu];
      for (int nu = mu + 1; nu < D; ++nu) {
        const long s_pnu = s + L.stride[nu];
        Cmat<N> Pl = U(s, mu) * U(s_pmu, nu) * U(s_pnu, mu).dagger() * U(s, nu).dagger();
        local += Pl.trace().real();
      }
    }
  }
  return mpi_allreduce_sum(local, L.comm) / (static_cast<Real>(N) * L.n_plaq_global());
}

// Global gauge action S_g = (beta/N) sum (N - Re Tr U_pl).
template <int D, int N>
Real gauge_action_d(const DistGauge<D, N>& U, Real beta) {
  const DistLattice<D>& L = *U.lat;
  Real localReTr = 0.0;
  for (long i = 0; i < L.vint; ++i) {
    const long s = L.interior[i];
    for (int mu = 0; mu < D; ++mu) {
      const long s_pmu = s + L.stride[mu];
      for (int nu = mu + 1; nu < D; ++nu) {
        const long s_pnu = s + L.stride[nu];
        Cmat<N> Pl = U(s, mu) * U(s_pmu, nu) * U(s_pnu, mu).dagger() * U(s, nu).dagger();
        localReTr += Pl.trace().real();
      }
    }
  }
  const Real globalReTr = mpi_allreduce_sum(localReTr, L.comm);
  return (beta / N) * (static_cast<Real>(N) * L.n_plaq_global() - globalReTr);
}

enum class MpiIntegrator { Leapfrog, Omelyan2MN };

// Domain-decomposed pure-gauge HMC. Momenta refreshed with the global-index RNG key
// (decomposition-independent); halos re-exchanged before every force/action eval;
// dH and the accept decision are global (Allreduce + same RNG key on all ranks).
template <int D, int N>
struct DistGaugeHMC {
  DistLattice<D> lat;
  DistGauge<D, N> U;
  std::vector<AlgVec<N>> P, F;  // [store*D]
  Rng rng;
  Real beta = 2.3, tau = 1.0;
  int  nmd = 20;
  MpiIntegrator integ = MpiIntegrator::Omelyan2MN;
  std::uint64_t traj_count = 0, accept_count = 0;
  Real last_dH = 0.0;

  DistGaugeHMC(const std::array<int, D>& G, const std::array<int, D>& pg, std::uint64_t seed)
      : lat(G, pg), U(lat), P(static_cast<std::size_t>(lat.store) * D),
        F(static_cast<std::size_t>(lat.store) * D), rng(seed) {}

  Real kinetic() const {
    Real e = 0.0;
    for (long i = 0; i < lat.vint; ++i) {
      const long s = lat.interior[i];
      for (int mu = 0; mu < D; ++mu) { const AlgVec<N>& v = P[s * D + mu]; for (int a = 0; a < n_gen<N>(); ++a) e += v[a] * v[a]; }
    }
    return 0.5 * mpi_allreduce_sum(e, lat.comm);
  }
  Real hamiltonian() { halo_exchange(U); return kinetic() + gauge_action_d<D, N>(U, beta); }

  void compute_force() {
    const auto& gen = generators<N>();
    const Real c = beta / N;
    for (long i = 0; i < lat.vint; ++i) {
      const long s = lat.interior[i];
      for (int mu = 0; mu < D; ++mu) {
        const Cmat<N> Omega = U(s, mu) * staple_d<D, N>(U, s, mu);
        AlgVec<N>& Fa = F[s * D + mu];
        for (int a = 0; a < n_gen<N>(); ++a) Fa[a] = c * trProd(gen.T[a], Omega).imag();
      }
    }
  }
  void kick(Real eps) {
    halo_exchange(U);
    compute_force();
    for (long i = 0; i < lat.vint; ++i) {
      const long s = lat.interior[i];
      for (int mu = 0; mu < D; ++mu) { AlgVec<N>& p = P[s * D + mu]; const AlgVec<N>& f = F[s * D + mu]; for (int a = 0; a < n_gen<N>(); ++a) p[a] -= eps * f[a]; }
    }
  }
  void drift(Real eps) {
    for (long i = 0; i < lat.vint; ++i) {
      const long s = lat.interior[i];
      for (int mu = 0; mu < D; ++mu) { Cmat<N> Pm = alg_to_mat<N>(P[s * D + mu]); Pm *= Complex(eps, 0.0); U(s, mu) = expi<N>(Pm) * U(s, mu); }
    }
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
    const std::uint64_t stream = Rng::key(0xA1, traj_count);
    for (long i = 0; i < lat.vint; ++i) {
      const long s = lat.interior[i], g = lat.gidx[i];
      for (int mu = 0; mu < D; ++mu) { AlgVec<N>& v = P[s * D + mu]; for (int a = 0; a < n_gen<N>(); ++a) v[a] = rng.gauss(Rng::key(stream, g, mu, a)); }
    }
    std::vector<Cmat<N>> save = U.u;
    const Real Hi = hamiltonian();
    md_evolve();
    U.reunitarize_interior();
    const Real Hf = hamiltonian();
    last_dH = Hf - Hi; ++traj_count;
    const double r = rng.uniform(Rng::key(0xB2, traj_count));  // same on all ranks
    bool acc = (last_dH <= 0.0) || (r < std::exp(-last_dH));
    if (acc) ++accept_count; else U.u = save;
    return acc;
  }
  double acceptance() const { return traj_count ? double(accept_count) / double(traj_count) : 0.0; }
};

}  // namespace gh
