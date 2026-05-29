// Stage-0/Stage-2 tests for the pure-gauge sector: force finite-difference,
// plaquette sanity, integrator reversibility, energy-conservation scaling, and
// the <exp(-dH)> identity / acceptance for SU(2) and SU(3).
#include "check.hpp"
#include "hmc/hmc.hpp"
#include "measure/observables.hpp"
#include <numeric>

using namespace gh;

// Finite-difference check of the gauge force against dS/d(omega^a) for one link.
template <int D, int N>
static void test_force_fd(std::uint64_t seed) {
  std::array<int, D> L{};
  for (int mu = 0; mu < D; ++mu) L[mu] = 4;
  Lattice<D> lat(L);
  GaugeField<D, N> U(lat);
  Rng rng(seed);
  U.hot(rng, 0.6);

  LinkMom<D, N> F(lat);
  F.zero();
  const Real beta = 2.0;
  add_gauge_force<D, N>(U, beta, F);

  const auto& gen = generators<N>();
  const Real eps = 1e-5;
  Real worst = 0.0;
  // sample a handful of (site, mu, a) -- indices must be in-bounds for all D
  std::int64_t sites[] = {0, lat.vol / 4, lat.vol / 2, lat.vol - 1};
  for (std::int64_t s : sites)
    for (int mu = 0; mu < D; ++mu)
      for (int a = 0; a < n_gen<N>(); ++a) {
        Cmat<N> Ta = gen.T[a];
        Cmat<N> orig = U(s, mu);
        // +eps
        U(s, mu) = expi<N>(Ta * Complex(eps, 0.0)) * orig;
        Real Sp = gauge_action<D, N>(U, beta);
        // -eps
        U(s, mu) = expi<N>(Ta * Complex(-eps, 0.0)) * orig;
        Real Sm = gauge_action<D, N>(U, beta);
        U(s, mu) = orig;
        Real fd = (Sp - Sm) / (2.0 * eps);
        worst = std::max(worst, std::fabs(fd - F(s, mu)[a]));
      }
  char msg[64];
  std::snprintf(msg, sizeof msg, "gauge force FD match SU(%d) D=%d", N, D);
  CHECK_CLOSE(worst, 0.0, 1e-6, msg);
}

template <int D, int N>
static void test_cold_plaquette() {
  std::array<int, D> L{};
  for (int mu = 0; mu < D; ++mu) L[mu] = 4;
  Lattice<D> lat(L);
  GaugeField<D, N> U(lat);
  U.cold();
  Real plaq = avg_plaquette<D, N>(U);
  Real act  = gauge_action<D, N>(U, 2.5);
  CHECK_CLOSE(plaq, 1.0, 1e-13, "cold plaquette = 1");
  CHECK_CLOSE(act, 0.0, 1e-9, "cold action = 0");
}

// Integrator reversibility: evolve, flip P, evolve again -> return to start.
template <int D, int N>
static void test_reversibility() {
  std::array<int, D> L{};
  for (int mu = 0; mu < D; ++mu) L[mu] = 4;
  GaugeHMC<D, N> hmc(L, 999);
  hmc.beta = 2.3; hmc.tau = 1.0; hmc.nmd = 12;
  hmc.U.hot(hmc.rng, 0.5);
  hmc.P.refresh(hmc.rng, 7);
  hmc.reunit_each_traj = false;

  std::vector<Cmat<N>> U0 = hmc.U.u;
  hmc.md_evolve();
  for (auto& v : hmc.P.p) for (int a = 0; a < n_gen<N>(); ++a) v[a] = -v[a];  // flip P
  hmc.md_evolve();

  Real worst = 0.0;
  for (std::size_t i = 0; i < U0.size(); ++i) worst = std::max(worst, (hmc.U.u[i] - U0[i]).fnorm());
  char msg[64];
  std::snprintf(msg, sizeof msg, "reversibility SU(%d) D=%d", N, D);
  CHECK_CLOSE(worst, 0.0, 1e-9, msg);
}

// Energy-conservation scaling: |dH| should fall ~4x when nmd doubles (2nd order).
template <int D, int N>
static void test_dH_scaling() {
  std::array<int, D> L{};
  for (int mu = 0; mu < D; ++mu) L[mu] = 4;
  auto run = [&](int nmd) {
    GaugeHMC<D, N> hmc(L, 31415);
    hmc.beta = 2.3; hmc.tau = 1.0; hmc.nmd = nmd; hmc.reunit_each_traj = false;
    hmc.U.hot(hmc.rng, 0.4);
    hmc.P.refresh(hmc.rng, 5);
    Real Hi = hmc.hamiltonian();
    hmc.md_evolve();
    Real Hf = hmc.hamiltonian();
    return std::fabs(Hf - Hi);
  };
  Real d1 = run(8), d2 = run(16);
  Real ratio = d1 / d2;
  char msg[80];
  std::snprintf(msg, sizeof msg, "dH 2nd-order scaling SU(%d) D=%d (ratio %.2f)", N, D, ratio);
  CHECK(ratio > 3.0 && ratio < 5.0, msg);
}

// <exp(-dH)> = 1 and reasonable acceptance over a short run.
template <int D, int N>
static void test_exp_dH(std::uint64_t seed) {
  std::array<int, D> L{};
  for (int mu = 0; mu < D; ++mu) L[mu] = 4;
  GaugeHMC<D, N> hmc(L, seed);
  hmc.beta = 2.3; hmc.tau = 1.0; hmc.nmd = 16;
  hmc.U.hot(hmc.rng, 0.3);
  const int ntraj = 400;
  double sum_exp = 0.0;
  for (int t = 0; t < ntraj; ++t) {
    hmc.trajectory();
    sum_exp += std::exp(-hmc.last_dH);
  }
  double mean = sum_exp / ntraj;
  char msg[96];
  std::snprintf(msg, sizeof msg, "<exp(-dH)>~1 SU(%d) D=%d (got %.4f, acc %.2f)",
                N, D, mean, hmc.acceptance());
  CHECK(std::fabs(mean - 1.0) < 0.12, msg);
  CHECK(hmc.acceptance() > 0.55, "acceptance > 0.55");
}

static void test_observables() {
  // cold lattice: Polyakov loop = 1
  std::array<int, 3> L{4, 4, 4};
  Lattice<3> lat(L);
  GaugeField<3, 2> U(lat); U.cold();
  Real p0 = polyakov_loop<3, 2>(U);
  Real p1 = polyakov_loop<3, 2>(U, 0);
  CHECK_CLOSE(p0, 1.0, 1e-12, "cold Polyakov (last dir) = 1");
  CHECK_CLOSE(p1, 1.0, 1e-12, "cold Polyakov (dir 0) = 1");
  // Stats on a known series {1,2,3,4,5}: mean 3, var 2, chi(V=10)=20.
  Stats s; for (int v = 1; v <= 5; ++v) s.add(v);
  CHECK_CLOSE(s.mean(), 3.0, 1e-12, "Stats mean");
  CHECK_CLOSE(s.var(), 2.0, 1e-12, "Stats var");
  CHECK_CLOSE(s.susceptibility(10.0), 20.0, 1e-12, "Stats susceptibility = V*var");
}

int main() {
  std::printf("-- observables --\n");
  test_observables();
  std::printf("-- gauge force finite difference --\n");
  test_force_fd<2, 2>(11); test_force_fd<3, 2>(12);
  test_force_fd<4, 2>(13); test_force_fd<3, 3>(14);
  std::printf("-- cold plaquette --\n");
  test_cold_plaquette<3, 2>(); test_cold_plaquette<4, 3>();
  std::printf("-- reversibility --\n");
  test_reversibility<3, 2>(); test_reversibility<3, 3>();
  std::printf("-- dH 2nd-order scaling --\n");
  test_dH_scaling<3, 2>();
  std::printf("-- <exp(-dH)> and acceptance --\n");
  test_exp_dH<3, 2>(2024); test_exp_dH<3, 3>(2025);
  return report("test_gauge");
}
