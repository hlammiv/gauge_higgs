// Tests for the DeGrand-Toussaint monopole density (compact U(1)), D=3 and D=4.
//
// We validate WITHOUT any HMC: cold configs (density 0), smooth gauge transforms of cold
// (no 2*pi wraps => 0), full gauge invariance of the density on a random config, the
// topological/current conservation laws on the periodic torus, and an explicit single
// Dirac monopole in D=3.  Detection is exact-integer arithmetic so tolerances are tight.
//
// NOTE: DeGrand-Toussaint counts are noisy below beta ~ 1; the density is a diagnostic of
// the (de)confining/Coulomb regimes, not a strict scaling order parameter.
#include "check.hpp"
#include "u1/monopole.hpp"
#include "u1/u1.hpp"
#include <random>
#include <vector>

using namespace gh;
using namespace gh::u1;

template <int D>
static std::array<int, D> cube(int n) { std::array<int, D> L{}; for (int mu = 0; mu < D; ++mu) L[mu] = n; return L; }

// ---- 1. Cold theta = 0 -> all n = 0 -> density exactly 0 ----
template <int D>
static void test_cold() {
  Lattice<D> lat(cube<D>(4));
  std::vector<Real> th(static_cast<std::size_t>(lat.vol) * D, 0.0);
  Real rho = monopole_density<D>(th, lat);
  char m[64]; std::snprintf(m, sizeof m, "monopole cold density==0 (D=%d)", D);
  CHECK_CLOSE(rho, 0.0, 0.0, m);
}

// ---- 2. Smooth gauge transform of cold: theta_mu(x) = a(x) - a(x+mu), |a| < 0.5 ----
// No plaquette can wrap 2*pi (each plaq angle is a sum of 4 such links, |.| < 2 < 2*pi
// before reduction, and is in fact identically 0 up to round-off), so density == 0.
template <int D>
static void test_smooth_gauge_of_cold(std::uint64_t seed) {
  Lattice<D> lat(cube<D>(4));
  std::mt19937_64 rng(seed); std::uniform_real_distribution<Real> u(-0.5, 0.5);
  std::vector<Real> a(lat.vol); for (auto& x : a) x = u(rng);
  std::vector<Real> th(static_cast<std::size_t>(lat.vol) * D, 0.0);
  for (std::int64_t x = 0; x < lat.vol; ++x)
    for (int mu = 0; mu < D; ++mu) th[x * D + mu] = a[x] - a[lat.neighbor_fwd(x, mu)];
  Real rho = monopole_density<D>(th, lat);
  char m[80]; std::snprintf(m, sizeof m, "monopole smooth-gauge-of-cold density==0 (D=%d)", D);
  CHECK_CLOSE(rho, 0.0, 0.0, m);
}

// ---- 3. Gauge invariance of the density on a random (wrapping) config ----
template <int D>
static void test_gauge_invariance(std::uint64_t seed) {
  Lattice<D> lat(cube<D>(4));
  std::mt19937_64 rng(seed); std::normal_distribution<Real> g(0, 1);
  // Larger angles so some plaquettes wrap and genuine monopoles exist.
  std::vector<Real> th(static_cast<std::size_t>(lat.vol) * D); for (auto& t : th) t = 2.5 * g(rng);
  const Real rho0 = monopole_density<D>(th, lat);
  // random gauge transform: theta_mu(x) += a(x) - a(x+mu)
  std::vector<Real> a(lat.vol); for (auto& x : a) x = 3.0 * g(rng);
  std::vector<Real> th2 = th;
  for (std::int64_t x = 0; x < lat.vol; ++x)
    for (int mu = 0; mu < D; ++mu) th2[x * D + mu] += a[x] - a[lat.neighbor_fwd(x, mu)];
  const Real rho1 = monopole_density<D>(th2, lat);
  char m[80];
  std::snprintf(m, sizeof m, "monopole density nonzero on random cfg (D=%d) rho=%.4f", D, rho0);
  CHECK(rho0 > 0.0, m);
  std::snprintf(m, sizeof m, "monopole density gauge invariant (D=%d)", D);
  CHECK_CLOSE(rho1, rho0, 1e-10, m);
}

// ---- 4a. D=3 topological conservation: sum over all cubes of signed charge == 0 ----
static void test_conservation_D3(std::uint64_t seed) {
  Lattice<3> lat(cube<3>(4));
  std::mt19937_64 rng(seed); std::normal_distribution<Real> g(0, 1);
  std::vector<Real> th(static_cast<std::size_t>(lat.vol) * 3); for (auto& t : th) t = 2.5 * g(rng);
  std::int64_t signed_sum = 0, abs_sum = 0;
  for (std::int64_t s = 0; s < lat.vol; ++s) {
    int mc = monopole_charge_D3<3>(th, lat, s);
    signed_sum += mc; abs_sum += std::abs(mc);
  }
  CHECK(abs_sum > 0, "monopole D=3 some charge present (sanity)");
  CHECK(signed_sum == 0, "monopole D=3 total signed charge == 0 (periodic torus)");
}

// ---- 4b. D=4 current conservation: discrete divergence == 0 at every dual site ----
static void test_conservation_D4(std::uint64_t seed) {
  Lattice<4> lat(cube<4>(4));
  std::mt19937_64 rng(seed); std::normal_distribution<Real> g(0, 1);
  std::vector<Real> th(static_cast<std::size_t>(lat.vol) * 4); for (auto& t : th) t = 2.5 * g(rng);
  // Precompute k_rho(s) for all sites/directions.
  std::vector<int> k(static_cast<std::size_t>(lat.vol) * 4);
  std::int64_t abs_sum = 0;
  for (std::int64_t s = 0; s < lat.vol; ++s)
    for (int rho = 0; rho < 4; ++rho) {
      int kr = monopole_current_D4<4>(th, lat, s, rho);
      k[s * 4 + rho] = kr; abs_sum += std::abs(kr);
    }
  // divergence sum_rho [ k_rho(s) - k_rho(s - rho^) ] must vanish at every site.
  std::int64_t max_abs_div = 0;
  for (std::int64_t s = 0; s < lat.vol; ++s) {
    int div = 0;
    for (int rho = 0; rho < 4; ++rho) {
      const std::int64_t sm = lat.neighbor_bwd(s, rho);
      div += k[s * 4 + rho] - k[sm * 4 + rho];
    }
    max_abs_div = std::max<std::int64_t>(max_abs_div, std::abs(div));
  }
  CHECK(abs_sum > 0, "monopole D=4 some current present (sanity)");
  char m[80]; std::snprintf(m, sizeof m, "monopole D=4 current conserved (max|div|=%lld)", (long long)max_abs_div);
  CHECK(max_abs_div == 0, m);
}

// ---- 5. Explicit single Dirac monopole in D=3 (Wu-Yang / 't Hooft construction) ----
// On a periodic torus a single isolated magnetic charge is impossible (Gauss: total charge
// must be 0).  The smallest exact, gauge-invariant object is therefore a monopole-antimonopole
// pair.  We place a unit + monopole at dual site r_p (the centre of one cube) and a unit -
// monopole at the adjacent dual site r_m (the centre of the neighbouring cube in the z=dir-2
// direction), and build the lattice link angles from the textbook continuum monopole vector
// potential in the Dirac-string-along-+z gauge:
//
//     A(r) = (1/2) * (-y, x, 0) * (-1 - z/|r|) / (x^2 + y^2)     [string along +z]
//
// with theta_{x,mu} = A_mu evaluated at the midpoint of the link.  The 1/2 normalizes the
// total magnetic flux out of a small surface around the charge to exactly 2*pi (DeGrand-
// Toussaint's convention), so each cube picks up a unit DT charge.  Placing the +/- pair on
// adjacent dual sites makes the two Dirac strings overlap above r_m and cancel, leaving a
// clean, fully-internal dipole: exactly one +1 cube and one -1 cube, zero everywhere else.
// The DT charge is what we then measure back -> exactly one +1 and one -1 cube, total 0,
// density = 2/V.  (This is the rigorous minimal monopole check possible on the torus.)
static void monopole_A(double x, double y, double z, double sign, double a[3]) {
  const double r = std::sqrt(x * x + y * y + z * z);
  const double rho2 = x * x + y * y;
  if (rho2 < 1e-12) { a[0] = a[1] = a[2] = 0.0; return; }
  const double pref = 0.5 * sign * (-1.0 - z / r) / rho2;  // 1/2 -> total flux 2*pi
  a[0] = -y * pref; a[1] = x * pref; a[2] = 0.0;
}

static void test_single_monopole_D3() {
  Lattice<3> lat(cube<3>(6));
  std::vector<Real> th(static_cast<std::size_t>(lat.vol) * 3, 0.0);
  // + monopole at centre of the cube based at (2,2,2); - monopole one step up in dir 2.
  const double rp[3] = {2.5, 2.5, 2.5};
  const double rm[3] = {2.5, 2.5, 3.5};
  std::array<int, 3> x{};
  for (std::int64_t s = 0; s < lat.vol; ++s) {
    lat.coords(s, x);
    for (int mu = 0; mu < 3; ++mu) {
      double mid[3] = {(double)x[0], (double)x[1], (double)x[2]}; mid[mu] += 0.5;
      double ap[3], am[3];
      monopole_A(mid[0] - rp[0], mid[1] - rp[1], mid[2] - rp[2], +1.0, ap);
      monopole_A(mid[0] - rm[0], mid[1] - rm[1], mid[2] - rm[2], -1.0, am);
      th[s * 3 + mu] = ap[mu] + am[mu];
    }
  }
  // Tally per-cube charges.
  int n_plus = 0, n_minus = 0;
  std::int64_t signed_sum = 0, abs_sum = 0;
  for (std::int64_t s = 0; s < lat.vol; ++s) {
    int mc = monopole_charge_D3<3>(th, lat, s);
    signed_sum += mc; abs_sum += std::abs(mc);
    if (mc > 0) n_plus += mc;
    if (mc < 0) n_minus += -mc;
  }
  CHECK(signed_sum == 0, "monopole D=3 single-pair total charge == 0");
  CHECK(abs_sum == 2, "monopole D=3 single-pair |charge| sum == 2 (one +/- pair)");
  CHECK(n_plus == 1 && n_minus == 1, "monopole D=3 exactly one +1 and one -1 cube, rest 0");
  // density = (number of |charge| units) / V = 2 / V.
  Real rho = monopole_density<3>(th, lat);
  char m[80]; std::snprintf(m, sizeof m, "monopole D=3 single-pair density == 2/V (V=%lld)", (long long)lat.vol);
  CHECK_CLOSE(rho, 2.0 / static_cast<Real>(lat.vol), 1e-13, m);
}

int main() {
  std::printf("-- cold limit --\n");            test_cold<3>(); test_cold<4>();
  std::printf("-- smooth gauge of cold --\n");  test_smooth_gauge_of_cold<3>(101); test_smooth_gauge_of_cold<4>(102);
  std::printf("-- gauge invariance --\n");       test_gauge_invariance<3>(201); test_gauge_invariance<4>(202);
  std::printf("-- conservation laws --\n");      test_conservation_D3(301); test_conservation_D4(302);
  std::printf("-- single Dirac monopole (D=3) --\n"); test_single_monopole_D3();
  return report("test_monopole");
}
