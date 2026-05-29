#pragma once
// Wilson plaquette gauge action: staple, action value, average plaquette, and the
// HMC gauge force (component form, FD-validated). See theory_notes §1.1, §3.1.
#include "core/fields.hpp"

namespace gh {

// Sum of the 2(D-1) staples around link (s, mu):
//   Sigma = sum_{nu!=mu} [ U_nu(x+mu) U_mu(x+nu)^dag U_nu(x)^dag        (forward)
//                        + U_nu(x+mu-nu)^dag U_mu(x-nu)^dag U_nu(x-nu) ](backward)
template <int D, int N>
Cmat<N> staple(const GaugeField<D, N>& U, std::int64_t s, int mu) {
  const Lattice<D>& lat = *U.lat;
  Cmat<N> sig{};
  const std::int64_t s_pmu = lat.neighbor_fwd(s, mu);
  for (int nu = 0; nu < D; ++nu) {
    if (nu == mu) continue;
    const std::int64_t s_pnu  = lat.neighbor_fwd(s, nu);
    const std::int64_t s_mnu  = lat.neighbor_bwd(s, nu);
    const std::int64_t s_pmu_mnu = lat.neighbor_bwd(s_pmu, nu);
    // forward staple
    sig += U(s_pmu, nu) * U(s_pnu, mu).dagger() * U(s, nu).dagger();
    // backward staple
    sig += U(s_pmu_mnu, nu).dagger() * U(s_mnu, mu).dagger() * U(s_mnu, nu);
  }
  return sig;
}

// Total Wilson action S_g = (beta/N) sum_{x,mu<nu} Re Tr(1 - U_pl)
//                         = (beta/N) sum_pl (N - Re Tr U_pl).
template <int D, int N>
Real gauge_action(const GaugeField<D, N>& U, Real beta) {
  const Lattice<D>& lat = *U.lat;
  Real sumReTr = 0.0;
  #pragma omp parallel for schedule(static) reduction(+:sumReTr)
  for (std::int64_t s = 0; s < lat.vol; ++s)
    for (int mu = 0; mu < D; ++mu) {
      const std::int64_t s_pmu = lat.neighbor_fwd(s, mu);
      for (int nu = mu + 1; nu < D; ++nu) {
        const std::int64_t s_pnu = lat.neighbor_fwd(s, nu);
        // U_pl = U_mu(x) U_nu(x+mu) U_mu(x+nu)^dag U_nu(x)^dag
        Cmat<N> P = U(s, mu) * U(s_pmu, nu) * U(s_pnu, mu).dagger() * U(s, nu).dagger();
        sumReTr += P.trace().real();
      }
    }
  return (beta / N) * (static_cast<Real>(N) * lat.n_plaq() - sumReTr);
}

// Average plaquette P = <(1/N) Re Tr U_pl> in [.. , 1].
template <int D, int N>
Real avg_plaquette(const GaugeField<D, N>& U) {
  const Lattice<D>& lat = *U.lat;
  Real sumReTr = 0.0;
  #pragma omp parallel for schedule(static) reduction(+:sumReTr)
  for (std::int64_t s = 0; s < lat.vol; ++s)
    for (int mu = 0; mu < D; ++mu) {
      const std::int64_t s_pmu = lat.neighbor_fwd(s, mu);
      for (int nu = mu + 1; nu < D; ++nu) {
        const std::int64_t s_pnu = lat.neighbor_fwd(s, nu);
        Cmat<N> P = U(s, mu) * U(s_pmu, nu) * U(s_pnu, mu).dagger() * U(s, nu).dagger();
        sumReTr += P.trace().real();
      }
    }
  return sumReTr / (static_cast<Real>(N) * lat.n_plaq());
}

// Gauge force component F^a_mu(x) = (beta/N) Im Tr( T^a U_mu(x) Sigma_mu(x) ).
// EOM dP^a/dt = -F^a. Accumulates into `F` (added, not overwritten) so the same
// container can collect gauge + matter forces.
template <int D, int N>
void add_gauge_force(const GaugeField<D, N>& U, Real beta, LinkMom<D, N>& F) {
  const Lattice<D>& lat = *U.lat;
  const auto& gen = generators<N>();
  const Real c = beta / N;
  #pragma omp parallel for schedule(static)
  for (std::int64_t s = 0; s < lat.vol; ++s)
    for (int mu = 0; mu < D; ++mu) {
      const Cmat<N> Omega = U(s, mu) * staple<D, N>(U, s, mu);
      AlgVec<N>& Fa = F(s, mu);
      for (int a = 0; a < n_gen<N>(); ++a)
        Fa[a] += c * trProd(gen.T[a], Omega).imag();
    }
}

}  // namespace gh
