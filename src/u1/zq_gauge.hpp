#pragma once
// Pure Z_q gauge theory in D dimensions: links theta = 2 pi n/q (n in 0..q-1), Wilson action
// S = beta * sum_plaq (1 - cos theta_plaq). Updated by an EXACT Z_q heatbath per link (sample the sector
// from its q-way Boltzmann conditional). No matter, no continuous part.
//
// WHY: the kappa->infinity limit of the frozen U(1)+charge-q Higgs IS pure Z_q gauge theory -- the matter
// locks the link onto the residual Z_q (e^{i q theta}=const), so the within-sector continuous fluctuation is
// suppressed ~1/kappa and the link becomes a pure Z_q element at the same beta. So
//   <plaq>, sigma_1, ... of U(1)+q at (beta, kappa->inf)  ->  pure Z_q at beta.
// Demonstrating that convergence as kappa grows is the DIGITIZATION PROOF (the discrete gauge theory IS the
// large-kappa limit). The observables are functions of theta only -> reuse gauge_obs.hpp / monopole.hpp.
#include "u1/u1.hpp"          // plaq_angle, avg_plaquette, gauge_action
#include "core/geometry.hpp"
#include "core/rng.hpp"
#include "core/config.hpp"
#include <vector>
#include <array>
#include <complex>
#include <cmath>
#include <algorithm>

namespace gh {
namespace u1 {

template <int D>
struct ZqGauge {
  Lattice<D> lat;
  std::vector<Real> th;      // link angles = 2 pi n/q
  Rng rng;
  Real beta = 1.0;
  int  q = 2;
  std::uint64_t sweep_count = 0;

  ZqGauge(const std::array<int, D>& ext, std::uint64_t seed)
      : lat(ext), th(static_cast<std::size_t>(lat.vol) * D, 0.0), rng(seed) {}

  void hot() {
    for (std::int64_t s = 0; s < lat.vol; ++s)
      for (int mu = 0; mu < D; ++mu) {
        int m = static_cast<int>(rng.uniform(Rng::key(1, s, mu)) * q); if (m >= q) m = q - 1;
        th[s * D + mu] = 2.0 * kPi * m / q;
      }
  }
  void cold() { std::fill(th.begin(), th.end(), 0.0); }
  Real avg_plaq() const { return avg_plaquette<D>(th, lat); }
  Real A() const { return gauge_action<D>(th, lat, 1.0); }

  // theta-independent staple resultant: local_gauge_cos = Re[e^{i theta} G].
  Complex staple(std::int64_t x, int mu) const {
    Complex G(0, 0); const Real thm = th[x * D + mu];
    for (int nu = 0; nu < D; ++nu) {
      if (nu == mu) continue;
      const Real pf = plaq_angle<D>(th, lat, x, mu, nu);
      const Real pb = plaq_angle<D>(th, lat, lat.neighbor_bwd(x, nu), mu, nu);
      G += std::polar(1.0, pf - thm);
      G += std::polar(1.0, -(pb + thm));
    }
    return G;
  }

  // Exact Z_q heatbath: resample each link's sector m in 0..q-1 from P(m) ∝ exp(beta Re[e^{i 2pi m/q} G]).
  void sweep() {
    const std::uint64_t t = sweep_count;
    for (std::int64_t x = 0; x < lat.vol; ++x)
      for (int mu = 0; mu < D; ++mu) {
        const std::size_t idx = x * D + mu;
        const Complex G = staple(x, mu);
        Real lw[256]; Real wmax = -1e300;
        for (int m = 0; m < q; ++m) { lw[m] = beta * (std::polar(1.0, 2.0 * kPi * m / q) * G).real(); wmax = std::max(wmax, lw[m]); }
        Real sum = 0.0; for (int m = 0; m < q; ++m) { lw[m] = std::exp(lw[m] - wmax); sum += lw[m]; }
        Real r = rng.uniform(Rng::key(0xD0, t, x, mu, 0)) * sum; int sel = 0;
        for (int m = 0; m < q; ++m) { r -= lw[m]; if (r <= 0.0) { sel = m; break; } }
        th[idx] = 2.0 * kPi * sel / q;
      }
    ++sweep_count;
  }
  void thermalize(int n) { for (int i = 0; i < n; ++i) sweep(); }
};

}  // namespace u1
}  // namespace gh
