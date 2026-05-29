#pragma once
// Extensive per-configuration observables for the (beta,kappa) scan / reweighting.
//
// REWEIGHTING CONVENTION (reweight.hpp depends on this EXACTLY):
//   weight = exp(-S),  S = beta*A - kappa*B + (on-site scalar terms), where
//     A = sum_plaq (1 - cos theta_plaq)                              [EXTENSIVE]
//         = gauge_action(theta,lat,1.0);  conjugate variable is beta.
//     B = sum_{x,mu} 2 Re[ conj(phi_x) e^{i q theta_mu(x)} phi_{x+mu} ]  [EXTENSIVE]
//         the hopping sum; scalar_action contains the term -kappa*B, so the
//         variable conjugate to kappa is -B.
//   For reweight.hpp written as S = sum_i lambda_i E_i :
//     (E_1, E_2)      = (A, -B)
//     (lambda_1, lambda_2) = (beta, kappa).
//   plaq_energy_sum gives A, hop_energy_sum gives B; both EXTENSIVE (no /vol).
#include "u1/u1.hpp"
#include "core/geometry.hpp"
#include "core/config.hpp"
#include <vector>
#include <cmath>
#include <complex>

namespace gh {
namespace u1 {

// A = sum over unordered plaquettes {mu<nu}, each counted once, of (1 - cos theta_plaq).
// EXTENSIVE. Identically equal to gauge_action<D>(th, lat, 1.0); this is the variable
// conjugate to beta in S = beta*A - kappa*B + on-site.
template <int D>
Real plaq_energy_sum(const std::vector<Real>& th, const Lattice<D>& lat) {
  Real s = 0.0;
  #pragma omp parallel for schedule(static) reduction(+:s)
  for (std::int64_t x = 0; x < lat.vol; ++x)
    for (int mu = 0; mu < D; ++mu)
      for (int nu = mu + 1; nu < D; ++nu)
        s += 1.0 - std::cos(plaq_angle<D>(th, lat, x, mu, nu));
  return s;
}

// B = sum_{x,mu} 2 Re[ conj(phi_x) e^{i q theta_mu(x)} phi_{x+mu} ].
// EXTENSIVE hopping sum. scalar_action contains -kappa*B, so -B is conjugate to kappa.
// Related to the AVERAGED link_energy by  B = 2 * (vol*D) * link_energy<D>(...).
template <int D>
Real hop_energy_sum(const std::vector<Complex>& phi, const std::vector<Real>& th,
                    const Lattice<D>& lat, int q) {
  Real s = 0.0;
  #pragma omp parallel for schedule(static) reduction(+:s)
  for (std::int64_t x = 0; x < lat.vol; ++x)
    for (int mu = 0; mu < D; ++mu) {
      const std::int64_t y = lat.neighbor_fwd(x, mu);
      const Complex ph = std::polar(1.0, q * th[x * D + mu]);  // e^{i q theta_mu(x)}
      s += 2.0 * (std::conj(phi[x]) * ph * phi[y]).real();
    }
  return s;
}

}  // namespace u1
}  // namespace gh
