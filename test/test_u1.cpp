// Tests for the compact U(1) + charge-q Higgs module: cold limit, finite-difference
// checks of the gauge/scalar/matter forces, integrator reversibility, <exp(-dH)>~1,
// and charge-q gauge invariance of the observables.
#include "check.hpp"
#include "u1/u1.hpp"
#include <random>

using namespace gh;
using namespace gh::u1;

template <int D>
static std::array<int, D> cube(int n) { std::array<int, D> L{}; for (int mu = 0; mu < D; ++mu) L[mu] = n; return L; }

template <int D>
static void test_cold() {
  Lattice<D> lat(cube<D>(4));
  std::vector<Real> th(lat.vol * D, 0.0);
  Real pl = avg_plaquette<D>(th, lat), ac = gauge_action<D>(th, lat, 1.7);
  CHECK_CLOSE(pl, 1.0, 1e-13, "U(1) cold plaquette = 1");
  CHECK_CLOSE(ac, 0.0, 1e-12, "U(1) cold gauge action = 0");
}

// FD of the TOTAL link force (gauge + matter) vs d(S_g+S_H)/dtheta, and the gauge force alone.
template <int D>
static void test_link_force_fd(int q, std::uint64_t seed) {
  Lattice<D> lat(cube<D>(4));
  std::mt19937_64 rng(seed); std::normal_distribution<Real> g(0, 1);
  std::vector<Real> th(lat.vol * D);
  std::vector<Complex> phi(lat.vol);
  for (auto& t : th) t = 0.6 * g(rng);
  for (auto& z : phi) z = Complex(0.7 * g(rng), 0.7 * g(rng));
  const Real beta = 1.8, kappa = 0.35, lambda = 0.4;

  std::vector<Real> Fg(lat.vol * D, 0.0); add_gauge_force<D>(th, lat, beta, Fg);
  std::vector<Real> Ft(lat.vol * D, 0.0); add_gauge_force<D>(th, lat, beta, Ft); add_matter_force<D>(phi, th, lat, q, kappa, Ft);
  const Real eps = 1e-6; Real wg = 0, wt = 0;
  std::int64_t sites[] = {0, 5, lat.vol - 1};
  for (std::int64_t s : sites)
    for (int mu = 0; mu < D; ++mu) {
      const std::size_t idx = s * D + mu; const Real o = th[idx];
      th[idx] = o + eps; Real Sgp = gauge_action<D>(th, lat, beta); Real Stp = Sgp + scalar_action<D>(phi, th, lat, q, kappa, lambda);
      th[idx] = o - eps; Real Sgm = gauge_action<D>(th, lat, beta); Real Stm = Sgm + scalar_action<D>(phi, th, lat, q, kappa, lambda);
      th[idx] = o;
      wg = std::max(wg, std::fabs((Sgp - Sgm) / (2 * eps) - Fg[idx]));
      wt = std::max(wt, std::fabs((Stp - Stm) / (2 * eps) - Ft[idx]));
    }
  char m[64];
  std::snprintf(m, sizeof m, "U(1) gauge force FD (q=%d)", q);  CHECK_CLOSE(wg, 0.0, 1e-6, m);
  std::snprintf(m, sizeof m, "U(1) gauge+matter link force FD (q=%d)", q); CHECK_CLOSE(wt, 0.0, 1e-6, m);
}

// FD of the scalar force: -dS_H/dRe(phi)=Re(F), -dS_H/dIm(phi)=Im(F).
template <int D>
static void test_scalar_force_fd(int q, std::uint64_t seed) {
  Lattice<D> lat(cube<D>(4));
  std::mt19937_64 rng(seed); std::normal_distribution<Real> g(0, 1);
  std::vector<Real> th(lat.vol * D); for (auto& t : th) t = 0.6 * g(rng);
  std::vector<Complex> phi(lat.vol); for (auto& z : phi) z = Complex(0.7 * g(rng), 0.7 * g(rng));
  const Real kappa = 0.35, lambda = 0.4;
  std::vector<Complex> F(lat.vol); scalar_force<D>(phi, th, lat, q, kappa, lambda, F);
  const Real eps = 1e-6; Real w = 0;
  std::int64_t sites[] = {0, 9, lat.vol - 1};
  for (std::int64_t s : sites) {
    const Complex o = phi[s];
    phi[s] = o + Complex(eps, 0); Real Sp = scalar_action<D>(phi, th, lat, q, kappa, lambda);
    phi[s] = o - Complex(eps, 0); Real Sm = scalar_action<D>(phi, th, lat, q, kappa, lambda);
    phi[s] = o + Complex(0, eps); Real Tp = scalar_action<D>(phi, th, lat, q, kappa, lambda);
    phi[s] = o - Complex(0, eps); Real Tm = scalar_action<D>(phi, th, lat, q, kappa, lambda);
    phi[s] = o;
    w = std::max(w, std::fabs(-(Sp - Sm) / (2 * eps) - F[s].real()));
    w = std::max(w, std::fabs(-(Tp - Tm) / (2 * eps) - F[s].imag()));
  }
  char m[64]; std::snprintf(m, sizeof m, "U(1) scalar force FD (q=%d)", q); CHECK_CLOSE(w, 0.0, 1e-6, m);
}

template <int D>
static void test_gauge_invariance(int q, std::uint64_t seed) {
  Lattice<D> lat(cube<D>(4));
  std::mt19937_64 rng(seed); std::normal_distribution<Real> g(0, 1);
  std::vector<Real> th(lat.vol * D); for (auto& t : th) t = 0.6 * g(rng);
  std::vector<Complex> phi(lat.vol); for (auto& z : phi) z = Complex(0.7 * g(rng), 0.7 * g(rng));
  const Real p0 = avg_plaquette<D>(th, lat), le0 = link_energy<D>(phi, th, lat, q);
  const Real w0 = wilson_loop<D>(th, lat, q, 2, 2), pk0 = polyakov<D>(th, lat, q);
  // random gauge transform: theta_mu(x) += a(x)-a(x+mu); phi_x *= e^{i q a(x)}
  std::vector<Real> a(lat.vol); for (auto& x : a) x = g(rng);
  std::vector<Real> th2 = th;
  for (std::int64_t x = 0; x < lat.vol; ++x)
    for (int mu = 0; mu < D; ++mu) th2[x * D + mu] += a[x] - a[lat.neighbor_fwd(x, mu)];
  std::vector<Complex> phi2(lat.vol);
  for (std::int64_t x = 0; x < lat.vol; ++x) phi2[x] = std::polar(1.0, q * a[x]) * phi[x];
  char m[64];
  std::snprintf(m, sizeof m, "U(1) gauge-inv plaquette (q=%d)", q); CHECK_CLOSE(avg_plaquette<D>(th2, lat), p0, 1e-10, m);
  std::snprintf(m, sizeof m, "U(1) gauge-inv link energy (q=%d)", q); CHECK_CLOSE(link_energy<D>(phi2, th2, lat, q), le0, 1e-10, m);
  Real w2 = wilson_loop<D>(th2, lat, q, 2, 2), pk2 = polyakov<D>(th2, lat, q);
  std::snprintf(m, sizeof m, "U(1) gauge-inv Wilson loop (q=%d)", q); CHECK_CLOSE(w2, w0, 1e-10, m);
  std::snprintf(m, sizeof m, "U(1) gauge-inv Polyakov (q=%d)", q);     CHECK_CLOSE(pk2, pk0, 1e-10, m);
}

template <int D>
static void test_reversibility(int q) {
  U1HMC<D> hmc(cube<D>(4), 99);
  hmc.beta = 1.5; hmc.kappa = 0.25; hmc.lambda = 0.5; hmc.q = q; hmc.tau = 1.0; hmc.nmd = 12;
  hmc.hot(0.5); hmc.refresh_momenta();
  auto th0 = hmc.th; auto phi0 = hmc.phi;
  hmc.md_evolve();
  for (auto& v : hmc.p) v = -v; for (auto& z : hmc.pi) z = -z;
  hmc.md_evolve();
  Real w = 0;
  for (std::size_t i = 0; i < th0.size(); ++i) w = std::max(w, std::fabs(hmc.th[i] - th0[i]));
  for (std::size_t i = 0; i < phi0.size(); ++i) w = std::max(w, std::abs(hmc.phi[i] - phi0[i]));
  char m[48]; std::snprintf(m, sizeof m, "U(1) HMC reversibility (q=%d)", q); CHECK_CLOSE(w, 0.0, 1e-9, m);
}

template <int D>
static void test_expdH(int q, std::uint64_t seed) {
  U1HMC<D> hmc(cube<D>(4), seed);
  hmc.beta = 1.5; hmc.kappa = 0.2; hmc.lambda = 0.5; hmc.q = q; hmc.tau = 1.0; hmc.nmd = 20;
  hmc.hot(0.3); hmc.cold_phi(0.5);
  const int ntraj = 300; double s = 0;
  for (int t = 0; t < ntraj; ++t) { hmc.trajectory(); s += std::exp(-hmc.last_dH); }
  const double mean = s / ntraj;
  char m[96]; std::snprintf(m, sizeof m, "U(1) <exp(-dH)>~1 (q=%d) got %.4f acc %.2f", q, mean, hmc.acceptance());
  CHECK(std::fabs(mean - 1.0) < 0.12, m);
  CHECK(hmc.acceptance() > 0.5, "U(1) acceptance > 0.5");
}

int main() {
  std::printf("-- cold limit --\n"); test_cold<3>(); test_cold<4>();
  std::printf("-- force finite differences --\n");
  test_link_force_fd<3>(2, 11); test_link_force_fd<3>(3, 12); test_link_force_fd<4>(2, 13);
  test_scalar_force_fd<3>(2, 21); test_scalar_force_fd<3>(3, 22);
  std::printf("-- gauge invariance --\n");
  test_gauge_invariance<3>(2, 31); test_gauge_invariance<3>(3, 32);
  std::printf("-- reversibility --\n"); test_reversibility<3>(2); test_reversibility<3>(3);
  std::printf("-- <exp(-dH)> / acceptance --\n"); test_expdH<3>(2, 41); test_expdH<3>(3, 42);
  return report("test_u1");
}
