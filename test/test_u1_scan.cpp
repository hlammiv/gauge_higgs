// Tests for the (beta,kappa) scan extensive observables (scan_obs.hpp).
//
// The CONVENTION test is the critical deliverable: it guarantees the A,B written by the
// driver are EXACTLY the action's (beta,kappa)-conjugate energies so downstream
// reweighting (reweight.hpp) is exact:
//   weight=exp(-S), S=beta*A-kappa*B+on-site; A=sum_plaq(1-cos), B=hopping sum.
// No physics scan here -- only machine-precision algebraic identities + gauge invariance.
#include "check.hpp"
#include "u1/u1.hpp"
#include "u1/scan_obs.hpp"
#include <random>
#include <vector>

using namespace gh;
using namespace gh::u1;

template <int D>
static std::array<int, D> cube(int n) { std::array<int, D> L{}; for (int mu = 0; mu < D; ++mu) L[mu] = n; return L; }

// Independent re-implementation of B = sum_{x,mu} 2 Re[ conj(phi_x) e^{iq theta} phi_{x+mu} ].
template <int D>
static Real hop_sum_ref(const std::vector<Complex>& phi, const std::vector<Real>& th,
                        const Lattice<D>& lat, int q) {
  Real s = 0.0;
  for (std::int64_t x = 0; x < lat.vol; ++x)
    for (int mu = 0; mu < D; ++mu) {
      const std::int64_t y = lat.neighbor_fwd(x, mu);
      const Complex u = std::exp(Complex(0.0, static_cast<Real>(q) * th[x * D + mu]));
      s += 2.0 * (std::conj(phi[x]) * u * phi[y]).real();
    }
  return s;
}

int main() {
  std::mt19937_64 rng(20260529ULL); std::normal_distribution<Real> g(0, 1);
  Lattice<4> lat(cube<4>(4));
  std::vector<Real> th(lat.vol * 4); for (auto& t : th) t = 0.6 * g(rng);
  std::vector<Complex> phi(lat.vol); for (auto& z : phi) z = Complex(0.7 * g(rng), 0.7 * g(rng));

  // ---- 1. CONVENTION (machine precision) ----
  std::printf("-- convention: A,B are the exact (beta,kappa)-conjugates --\n");
  // A == gauge_action(.,.,1.0).
  CHECK_CLOSE(plaq_energy_sum<4>(th, lat), gauge_action<4>(th, lat, 1.0), 1e-12,
              "A = plaq_energy_sum == gauge_action(th,lat,1.0)");
  // A is exactly the beta-conjugate: gauge_action(.,.,beta) == beta*A for several beta.
  for (Real beta : {0.3, 1.0, 1.7, 2.5, 4.0}) {
    char m[96]; std::snprintf(m, sizeof m, "gauge_action(beta=%.2f) == beta*A", beta);
    CHECK_CLOSE(gauge_action<4>(th, lat, beta), beta * plaq_energy_sum<4>(th, lat), 1e-12, m);
  }
  // B == independent direct re-implementation, for q = 1,2,3.
  for (int q = 1; q <= 3; ++q) {
    char m[64]; std::snprintf(m, sizeof m, "B = hop_energy_sum == ref (q=%d)", q);
    CHECK_CLOSE(hop_energy_sum<4>(phi, th, lat, q), hop_sum_ref<4>(phi, th, lat, q), 1e-12, m);
  }

  // ---- 2. EXTENSIVITY ties to the averaged forms ----
  std::printf("-- extensivity ties --\n");
  const std::int64_t n_plaq  = lat.vol * 4 * (4 - 1) / 2;   // = lat.n_plaq()
  const std::int64_t n_links = lat.vol * 4;
  CHECK_CLOSE(static_cast<Real>(n_plaq), static_cast<Real>(lat.n_plaq()), 1e-12, "N_plaq = vol*D(D-1)/2");
  // A / N_plaq == 1 - <plaq>.
  CHECK_CLOSE(plaq_energy_sum<4>(th, lat) / static_cast<Real>(n_plaq),
              1.0 - avg_plaquette<4>(th, lat), 1e-12, "A/N_plaq == 1 - avg_plaquette");
  // B == 2 * N_links * link_energy  (B has the 2Re; link_energy is Re averaged over vol*D).
  for (int q = 1; q <= 3; ++q) {
    char m[80]; std::snprintf(m, sizeof m, "B == 2*N_links*link_energy (q=%d)", q);
    CHECK_CLOSE(hop_energy_sum<4>(phi, th, lat, q),
                2.0 * static_cast<Real>(n_links) * link_energy<4>(phi, th, lat, q), 1e-12, m);
  }

  // ---- 3. GAUGE INVARIANCE ----
  std::printf("-- gauge invariance --\n");
  for (int q = 1; q <= 3; ++q) {
    std::vector<Real> a(lat.vol); for (auto& v : a) v = g(rng);
    std::vector<Real> th2 = th;
    for (std::int64_t x = 0; x < lat.vol; ++x)
      for (int mu = 0; mu < 4; ++mu) th2[x * 4 + mu] += a[x] - a[lat.neighbor_fwd(x, mu)];
    std::vector<Complex> phi2(lat.vol);
    for (std::int64_t x = 0; x < lat.vol; ++x) phi2[x] = std::polar(1.0, static_cast<Real>(q) * a[x]) * phi[x];
    char m[80];
    std::snprintf(m, sizeof m, "gauge-inv A (q=%d)", q);
    CHECK_CLOSE(plaq_energy_sum<4>(th2, lat), plaq_energy_sum<4>(th, lat), 1e-10, m);
    std::snprintf(m, sizeof m, "gauge-inv B (q=%d)", q);
    CHECK_CLOSE(hop_energy_sum<4>(phi2, th2, lat, q), hop_energy_sum<4>(phi, th, lat, q), 1e-10, m);
  }

  return report("test_u1_scan");
}
