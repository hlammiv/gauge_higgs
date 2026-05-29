// Tests for arbitrary-irrep Wilson loops + rep Polyakov loops:
// cold = 1; 1x1 fundamental loop == average plaquette; fundamental-rep loop ==
// fundamental loop; and gauge invariance under a random local gauge transform.
#include "check.hpp"
#include "hmc/hmc.hpp"
#include "measure/observables.hpp"
#include "rep/rep_fundamental.hpp"
#include "rep/rep_adjoint.hpp"
#include "rep/rep_general.hpp"
#include <random>

using namespace gh;

template <int D>
static std::array<int, D> cube(int n) { std::array<int, D> L{}; for (int mu = 0; mu < D; ++mu) L[mu] = n; return L; }

template <int N>
static Cmat<N> rand_su(std::mt19937_64& rng, Real s = 0.8) {
  std::normal_distribution<Real> g(0.0, s);
  AlgVec<N> v{};
  for (int a = 0; a < n_gen<N>(); ++a) v[a] = g(rng);
  return expi<N>(alg_to_mat<N>(v));
}

// U_mu(x) -> Omega(x) U_mu(x) Omega(x+mu)^dag (each link depends only on its own old value).
template <int D, int N>
static void gauge_transform(GaugeField<D, N>& U, const std::vector<Cmat<N>>& Om) {
  const Lattice<D>& lat = *U.lat;
  for (std::int64_t s = 0; s < lat.vol; ++s)
    for (int mu = 0; mu < D; ++mu)
      U(s, mu) = Om[s] * U(s, mu) * Om[lat.neighbor_fwd(s, mu)].dagger();
}

template <int D, int N>
static void test_cold() {
  Lattice<D> lat(cube<D>(4));
  GaugeField<D, N> U(lat); U.cold();
  FundamentalRep<N> fund; AdjointRep<N> adj;
  char m[64];
  int sizes[][2] = {{1, 1}, {2, 1}, {2, 2}, {3, 2}};
  for (auto& s : sizes) {
    Real w = wilson_loop_fund<D, N>(U, s[0], s[1]);
    std::snprintf(m, sizeof m, "cold W_fund(%d,%d)=1 SU(%d)", s[0], s[1], N);
    CHECK_CLOSE(w, 1.0, 1e-12, m);
  }
  Real wf = wilson_loop_rep<D, N>(U, fund, 2, 2);
  Real wa = wilson_loop_rep<D, N>(U, adj, 2, 2);
  Real pf = polyakov_loop_rep<D, N>(U, fund);
  Real pa = polyakov_loop_rep<D, N>(U, adj);
  CHECK_CLOSE(wf, 1.0, 1e-12, "cold W_rep(fund)=1");
  CHECK_CLOSE(wa, 1.0, 1e-12, "cold W_rep(adj)=1");
  CHECK_CLOSE(pf, 1.0, 1e-12, "cold Polyakov_rep(fund)=1");
  CHECK_CLOSE(pa, 1.0, 1e-12, "cold Polyakov_rep(adj)=1");
}

template <int D, int N>
static void test_consistency(std::uint64_t seed) {
  Lattice<D> lat(cube<D>(4));
  GaugeField<D, N> U(lat); Rng rng(seed); U.hot(rng, 0.7);
  FundamentalRep<N> fund;
  Real w11 = wilson_loop_fund<D, N>(U, 1, 1);
  Real plaq = avg_plaquette<D, N>(U);
  CHECK_CLOSE(w11, plaq, 1e-12, "1x1 fund loop == avg plaquette");
  Real wr = wilson_loop_rep<D, N>(U, fund, 2, 1);
  Real wfd = wilson_loop_fund<D, N>(U, 2, 1);
  CHECK_CLOSE(wr, wfd, 1e-12, "fund-rep loop == fund loop");
}

template <int D, int N>
static void test_gauge_invariance(std::uint64_t seed) {
  Lattice<D> lat(cube<D>(4));
  GaugeField<D, N> U(lat); Rng rng(seed); U.hot(rng, 0.7);
  FundamentalRep<N> fund; AdjointRep<N> adj;
  Real w21 = wilson_loop_fund<D, N>(U, 2, 1);
  Real w22 = wilson_loop_fund<D, N>(U, 2, 2);
  Real wa  = wilson_loop_rep<D, N>(U, adj, 2, 2);
  Real pa  = polyakov_loop_rep<D, N>(U, adj);
  // random local gauge transform
  std::mt19937_64 mt(seed + 99);
  std::vector<Cmat<N>> Om(lat.vol);
  for (std::int64_t s = 0; s < lat.vol; ++s) Om[s] = rand_su<N>(mt, 1.0);
  gauge_transform<D, N>(U, Om);
  Real w21b = wilson_loop_fund<D, N>(U, 2, 1);
  Real w22b = wilson_loop_fund<D, N>(U, 2, 2);
  Real wab  = wilson_loop_rep<D, N>(U, adj, 2, 2);
  Real pab  = polyakov_loop_rep<D, N>(U, adj);
  char m[64];
  std::snprintf(m, sizeof m, "gauge-inv W_fund(2,1) SU(%d)", N);
  CHECK_CLOSE(w21b, w21, 1e-10, m);
  CHECK_CLOSE(w22b, w22, 1e-10, "gauge-inv W_fund(2,2)");
  CHECK_CLOSE(wab, wa, 1e-10, "gauge-inv W_rep(adj)");
  CHECK_CLOSE(pab, pa, 1e-10, "gauge-inv Polyakov_rep(adj)");
}

int main() {
  std::printf("-- cold Wilson/Polyakov loops --\n");
  test_cold<3, 2>(); test_cold<4, 2>(); test_cold<3, 3>();
  std::printf("-- consistency (1x1==plaq, fund-rep==fund) --\n");
  test_consistency<3, 2>(11); test_consistency<3, 3>(12); test_consistency<4, 2>(13);
  std::printf("-- gauge invariance --\n");
  test_gauge_invariance<3, 2>(21); test_gauge_invariance<3, 3>(22);
  return report("test_wilson");
}
