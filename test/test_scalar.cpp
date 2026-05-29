// Stage-0/Stage-3/Stage-5 tests for representations + scalar Higgs sector:
// rep invariants, adjoint homomorphism/orthogonality, scalar-force and
// matter-link-force finite differences, combined gauge+Higgs reversibility and
// <exp(-dH)> for fundamental (complex) and adjoint (real) SU(2)/SU(3).
#include "check.hpp"
#include "hmc/gauge_higgs_hmc.hpp"
#include "rep/rep_fundamental.hpp"
#include "rep/rep_adjoint.hpp"
#include <random>

using namespace gh;

template <int D>
static std::array<int, D> cube(int n) { std::array<int, D> L{}; for (int mu = 0; mu < D; ++mu) L[mu] = n; return L; }

template <int N>
static void test_invariants() {
  FundamentalRep<N> fund;
  AdjointRep<N>     adj;
  const Real Tf = 0.5, C2f = (static_cast<Real>(N) * N - 1) / (2.0 * N);
  const Real Ta = N,   C2a = N;
  char m[80];
  std::snprintf(m, sizeof m, "fundamental invariants SU(%d)", N);
  CHECK_CLOSE(fund.invariant_violation(Tf, C2f), 0.0, 1e-10, m);
  std::snprintf(m, sizeof m, "adjoint invariants SU(%d) (T=%g,C2=%g)", N, Ta, C2a);
  CHECK_CLOSE(adj.invariant_violation(Ta, C2a), 0.0, 1e-10, m);
}

// Adjoint rep matrix must be a real orthogonal homomorphism.
template <int N>
static void test_adjoint_rep_matrix(std::uint64_t seed) {
  AdjointRep<N> adj;
  std::mt19937_64 rng(seed);
  std::normal_distribution<Real> g(0.0, 0.7);
  auto rand_su = [&]() {
    AlgVec<N> v{};
    for (int a = 0; a < n_gen<N>(); ++a) v[a] = g(rng);
    return expi<N>(alg_to_mat<N>(v));
  };
  const int n = n_gen<N>();
  for (int t = 0; t < 10; ++t) {
    Cmat<N> U1 = rand_su(), U2 = rand_su();
    DMat M1 = adj.rep_matrix(U1);
    // real + orthogonal: M M^T = I
    DMat MMt = M1 * M1.dagger();  // dagger == transpose for real M
    Real ovio = 0.0, imvio = 0.0;
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j) {
        ovio = std::max(ovio, std::abs(MMt(i, j) - Complex(i == j ? 1.0 : 0.0, 0.0)));
        imvio = std::max(imvio, std::fabs(M1(i, j).imag()));
      }
    CHECK_CLOSE(ovio, 0.0, 1e-9, "adjoint link orthogonal");
    CHECK_CLOSE(imvio, 0.0, 1e-12, "adjoint link real");
    // homomorphism: D(U1 U2) = D(U1) D(U2)
    DMat lhs = adj.rep_matrix(U1 * U2);
    DMat rhs = M1 * adj.rep_matrix(U2);
    Real hvio = 0.0;
    for (std::size_t i = 0; i < lhs.a.size(); ++i) hvio = std::max(hvio, std::abs(lhs.a[i] - rhs.a[i]));
    CHECK_CLOSE(hvio, 0.0, 1e-9, "adjoint homomorphism");
  }
}

// Finite-difference check of the scalar force: -dS_H/dRe(phi) = Re(F), -dS_H/dIm = Im(F).
template <int D, int N>
static void test_scalar_force_fd(Representation<N>& rep, std::uint64_t seed, const char* tag) {
  Lattice<D> lat(cube<D>(4));
  GaugeField<D, N> U(lat); Rng rng(seed); U.hot(rng, 0.5);
  CVecField<D> phi(lat, rep.d); phi.gaussian(rng, 99, rep.real, 0.7);
  const Real kappa = 0.3, lambda = 0.4;
  CVecField<D> F(lat, rep.d); scalar_force<D, N>(phi, U, rep, kappa, lambda, F);

  const Real eps = 1e-5; Real worst = 0.0;
  std::int64_t sites[] = {0, 3, 9, lat.vol - 1};
  for (std::int64_t s : sites)
    for (int k = 0; k < rep.d; ++k) {
      const std::size_t idx = static_cast<std::size_t>(s) * rep.d + k;
      const Complex orig = phi.data[idx];
      // d/dRe
      phi.data[idx] = orig + Complex(eps, 0); Real Sp = scalar_action<D, N>(phi, U, rep, kappa, lambda);
      phi.data[idx] = orig - Complex(eps, 0); Real Sm = scalar_action<D, N>(phi, U, rep, kappa, lambda);
      phi.data[idx] = orig;
      worst = std::max(worst, std::fabs(-(Sp - Sm) / (2 * eps) - F.data[idx].real()));
      if (!rep.real) {
        phi.data[idx] = orig + Complex(0, eps); Real Tp = scalar_action<D, N>(phi, U, rep, kappa, lambda);
        phi.data[idx] = orig - Complex(0, eps); Real Tm = scalar_action<D, N>(phi, U, rep, kappa, lambda);
        phi.data[idx] = orig;
        worst = std::max(worst, std::fabs(-(Tp - Tm) / (2 * eps) - F.data[idx].imag()));
      }
    }
  char m[96]; std::snprintf(m, sizeof m, "scalar force FD (%s)", tag);
  CHECK_CLOSE(worst, 0.0, 1e-5, m);
}

// Finite-difference check of the TOTAL link force (gauge + matter) vs d(S_g+S_H)/d(omega^a).
template <int D, int N>
static void test_link_force_fd(Representation<N>& rep, std::uint64_t seed, const char* tag) {
  Lattice<D> lat(cube<D>(4));
  GaugeField<D, N> U(lat); Rng rng(seed); U.hot(rng, 0.5);
  CVecField<D> phi(lat, rep.d); phi.gaussian(rng, 77, rep.real, 0.7);
  const Real beta = 2.0, kappa = 0.35, lambda = 0.4;
  LinkMom<D, N> F(lat); F.zero();
  add_gauge_force<D, N>(U, beta, F);
  add_matter_link_force<D, N>(phi, U, rep, kappa, F);

  const auto& gen = generators<N>();
  const Real eps = 1e-5; Real worst = 0.0;
  auto totalS = [&]() { return gauge_action<D, N>(U, beta) + scalar_action<D, N>(phi, U, rep, kappa, lambda); };
  std::int64_t sites[] = {0, 5, lat.vol - 1};
  for (std::int64_t s : sites)
    for (int mu = 0; mu < D; ++mu)
      for (int a = 0; a < n_gen<N>(); ++a) {
        const Cmat<N> orig = U(s, mu);
        U(s, mu) = expi<N>(gen.T[a] * Complex(eps, 0)) * orig;  Real Sp = totalS();
        U(s, mu) = expi<N>(gen.T[a] * Complex(-eps, 0)) * orig; Real Sm = totalS();
        U(s, mu) = orig;
        worst = std::max(worst, std::fabs((Sp - Sm) / (2 * eps) - F(s, mu)[a]));
      }
  char m[96]; std::snprintf(m, sizeof m, "link force (gauge+matter) FD (%s)", tag);
  CHECK_CLOSE(worst, 0.0, 1e-5, m);
}

template <int D, int N>
static void test_combined_reversibility(Representation<N>& rep, std::uint64_t seed, const char* tag) {
  GaugeHiggsHMC<D, N> hmc(cube<D>(4), rep, seed);
  hmc.beta = 2.3; hmc.kappa = 0.3; hmc.lambda = 0.5; hmc.tau = 1.0; hmc.nmd = 12;
  hmc.reunit_each_traj = false;
  hmc.U.hot(hmc.rng, 0.5);
  hmc.phi.gaussian(hmc.rng, 1, rep.real, 0.6);
  hmc.refresh_momenta();
  std::vector<Cmat<N>> U0 = hmc.U.u;
  std::vector<Complex> phi0 = hmc.phi.data;
  hmc.md_evolve();
  for (auto& v : hmc.P.p) for (int a = 0; a < n_gen<N>(); ++a) v[a] = -v[a];
  for (auto& z : hmc.pi.data) z = -z;
  hmc.md_evolve();
  Real worst = 0.0;
  for (std::size_t i = 0; i < U0.size(); ++i) worst = std::max(worst, (hmc.U.u[i] - U0[i]).fnorm());
  for (std::size_t i = 0; i < phi0.size(); ++i) worst = std::max(worst, std::abs(hmc.phi.data[i] - phi0[i]));
  char m[96]; std::snprintf(m, sizeof m, "combined reversibility (%s)", tag);
  CHECK_CLOSE(worst, 0.0, 1e-9, m);
}

template <int D, int N>
static void test_combined_expdH(Representation<N>& rep, std::uint64_t seed, const char* tag) {
  GaugeHiggsHMC<D, N> hmc(cube<D>(4), rep, seed);
  hmc.beta = 2.3; hmc.kappa = 0.2; hmc.lambda = 0.5; hmc.tau = 1.0; hmc.nmd = 20;
  hmc.U.hot(hmc.rng, 0.3);
  hmc.phi.cold(0.5);
  const int ntraj = 400; double s = 0.0;
  for (int t = 0; t < ntraj; ++t) { hmc.trajectory(); s += std::exp(-hmc.last_dH); }
  const double mean = s / ntraj;
  char m[128];
  std::snprintf(m, sizeof m, "<exp(-dH)>~1 combined (%s) got %.4f acc %.2f", tag, mean, hmc.acceptance());
  CHECK(std::fabs(mean - 1.0) < 0.15, m);
  CHECK(hmc.acceptance() > 0.5, "combined acceptance > 0.5");
}

int main() {
  std::printf("-- representation invariants --\n");
  test_invariants<2>(); test_invariants<3>();
  std::printf("-- adjoint rep matrix (orthogonal homomorphism) --\n");
  test_adjoint_rep_matrix<2>(101); test_adjoint_rep_matrix<3>(102);

  FundamentalRep<2> f2; AdjointRep<2> a2; FundamentalRep<3> f3; AdjointRep<3> a3;
  std::printf("-- scalar force finite difference --\n");
  test_scalar_force_fd<3, 2>(f2, 201, "SU(2) fund");
  test_scalar_force_fd<3, 2>(a2, 202, "SU(2) adj");
  test_scalar_force_fd<3, 3>(f3, 203, "SU(3) fund");
  test_scalar_force_fd<3, 3>(a3, 204, "SU(3) adj");
  std::printf("-- link force (gauge+matter) finite difference --\n");
  test_link_force_fd<3, 2>(f2, 301, "SU(2) fund");
  test_link_force_fd<3, 2>(a2, 302, "SU(2) adj");
  test_link_force_fd<3, 3>(f3, 303, "SU(3) fund");
  std::printf("-- combined reversibility --\n");
  test_combined_reversibility<3, 2>(f2, 401, "SU(2) fund");
  test_combined_reversibility<3, 2>(a2, 402, "SU(2) adj");
  std::printf("-- combined <exp(-dH)> / acceptance --\n");
  test_combined_expdH<3, 2>(f2, 501, "SU(2) fund");
  test_combined_expdH<3, 2>(a2, 502, "SU(2) adj");
  return report("test_scalar");
}
