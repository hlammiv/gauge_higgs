// Validation of the general highest-weight (tensor + Young symmetrizer) irrep
// engine: dimensions, N-ality, the algebra commutator identity [T^a,T^b]=if^abc T^c
// (the definitive proof it is a representation), Dynkin/Casimir invariants,
// basis-independent CHARACTER cross-checks against the explicit fundamental and
// adjoint reps, rep_matrix homomorphism, force finite differences, and combined HMC.
#include "check.hpp"
#include "hmc/gauge_higgs_hmc.hpp"
#include "rep/rep_fundamental.hpp"
#include "rep/rep_adjoint.hpp"
#include "rep/rep_general.hpp"
#include <random>

using namespace gh;

template <int N>
static Cmat<N> rand_su(std::mt19937_64& rng, Real s = 0.8) {
  std::normal_distribution<Real> g(0.0, s);
  AlgVec<N> v{};
  for (int a = 0; a < n_gen<N>(); ++a) v[a] = g(rng);
  return expi<N>(alg_to_mat<N>(v));
}

template <int D>
static std::array<int, D> cube(int nn) { std::array<int, D> L{}; for (int mu = 0; mu < D; ++mu) L[mu] = nn; return L; }

// [T^a_R, T^b_R] = i f^{abc} T^c_R  -- proves T^a_R is a representation of su(N).
template <int N>
static void test_commutator(const Representation<N>& rep, const char* tag) {
  const auto& f = structure_constants<N>();
  const int ng = n_gen<N>(), d = rep.d;
  Real worst = 0.0;
  for (int a = 0; a < ng; ++a)
    for (int b = 0; b < ng; ++b) {
      DMat comm = rep.T[a] * rep.T[b] - rep.T[b] * rep.T[a];
      DMat rhs(d, d);
      for (int c = 0; c < ng; ++c) {
        const Real fabc = f[(static_cast<std::size_t>(a) * ng + b) * ng + c];
        if (fabc != 0.0) rhs = rhs + rep.T[c] * Complex(0.0, fabc);  // i f^{abc} T^c
      }
      worst = std::max(worst, fnorm(comm - rhs));
    }
  char m[96]; std::snprintf(m, sizeof m, "commutator [T,T]=ifT (%s)", tag);
  CHECK_CLOSE(worst, 0.0, 1e-9, m);
}

// Compute the Dynkin index T(R) and Casimir C2(R) from the generators and check the
// rep is internally consistent; also check them against expected values.
template <int N>
static void test_invariants(const Representation<N>& rep, int expect_d, Real expect_T, Real expect_C2,
                            const char* tag) {
  char m[96];
  std::snprintf(m, sizeof m, "dim (%s)", tag);
  CHECK(rep.d == expect_d, m);
  // T(R) from Tr(T^0 T^0); C2 from d*C2 = T*(N^2-1).
  const Real T = (rep.T[0] * rep.T[0]).trace().real();
  const Real C2 = T * (static_cast<Real>(N) * N - 1) / rep.d;
  std::snprintf(m, sizeof m, "Dynkin index (%s) got %.4f", tag, T);
  CHECK_CLOSE(T, expect_T, 1e-7, m);
  std::snprintf(m, sizeof m, "Casimir (%s) got %.4f", tag, C2);
  CHECK_CLOSE(C2, expect_C2, 1e-7, m);
  // internal consistency: Tr(T^aT^b)=T delta, sum T^aT^a = C2 I
  std::snprintf(m, sizeof m, "invariant consistency (%s)", tag);
  CHECK_CLOSE(rep.invariant_violation(T, C2), 0.0, 1e-8, m);
}

// Character chi(U) = Tr D_R(U) must match a reference (basis-independent).
template <int N, class FRef>
static void test_character(const Representation<N>& rep, FRef ref, std::uint64_t seed, const char* tag) {
  std::mt19937_64 rng(seed);
  Real worst = 0.0;
  for (int t = 0; t < 8; ++t) {
    Cmat<N> U = rand_su<N>(rng);
    Complex chi = rep.rep_matrix(U).trace();
    Complex chiref = ref(U);
    worst = std::max(worst, std::abs(chi - chiref));
  }
  char m[96]; std::snprintf(m, sizeof m, "character match (%s)", tag);
  CHECK_CLOSE(worst, 0.0, 1e-9, m);
}

template <int N>
static void test_homomorphism(const Representation<N>& rep, std::uint64_t seed, const char* tag) {
  std::mt19937_64 rng(seed);
  Real worst = 0.0;
  for (int t = 0; t < 6; ++t) {
    Cmat<N> U1 = rand_su<N>(rng), U2 = rand_su<N>(rng);
    DMat lhs = rep.rep_matrix(U1 * U2);
    DMat rhs = rep.rep_matrix(U1) * rep.rep_matrix(U2);
    worst = std::max(worst, fnorm(lhs - rhs));
  }
  char m[96]; std::snprintf(m, sizeof m, "homomorphism (%s)", tag);
  CHECK_CLOSE(worst, 0.0, 1e-9, m);
}

// Scalar + link force finite differences for a general rep (reusing the validated
// rep-agnostic machinery, but exercising GeneralRep's T_R and rep_matrix).
template <int Dd, int N>
static void test_force_fd(Representation<N>& rep, std::uint64_t seed, const char* tag) {
  Lattice<Dd> lat(cube<Dd>(4));
  GaugeField<Dd, N> U(lat); Rng rng(seed); U.hot(rng, 0.5);
  CVecField<Dd> phi(lat, rep.d); phi.gaussian(rng, 55, rep.real, 0.7);
  const Real beta = 2.0, kappa = 0.3, lambda = 0.4;

  // scalar force
  CVecField<Dd> Fp(lat, rep.d); scalar_force<Dd, N>(phi, U, rep, kappa, lambda, Fp);
  const Real eps = 1e-5; Real ws = 0.0;
  std::int64_t sites[] = {0, 7, lat.vol - 1};
  for (std::int64_t s : sites)
    for (int k = 0; k < rep.d; ++k) {
      const std::size_t idx = static_cast<std::size_t>(s) * rep.d + k;
      const Complex orig = phi.data[idx];
      phi.data[idx] = orig + Complex(eps, 0); Real Sp = scalar_action<Dd, N>(phi, U, rep, kappa, lambda);
      phi.data[idx] = orig - Complex(eps, 0); Real Sm = scalar_action<Dd, N>(phi, U, rep, kappa, lambda);
      phi.data[idx] = orig;
      ws = std::max(ws, std::fabs(-(Sp - Sm) / (2 * eps) - Fp.data[idx].real()));
      if (!rep.real) {
        phi.data[idx] = orig + Complex(0, eps); Real Tp = scalar_action<Dd, N>(phi, U, rep, kappa, lambda);
        phi.data[idx] = orig - Complex(0, eps); Real Tm = scalar_action<Dd, N>(phi, U, rep, kappa, lambda);
        phi.data[idx] = orig;
        ws = std::max(ws, std::fabs(-(Tp - Tm) / (2 * eps) - Fp.data[idx].imag()));
      }
    }
  char m[96]; std::snprintf(m, sizeof m, "scalar force FD (%s)", tag);
  CHECK_CLOSE(ws, 0.0, 1e-5, m);

  // total link force
  LinkMom<Dd, N> F(lat); F.zero();
  add_gauge_force<Dd, N>(U, beta, F);
  add_matter_link_force<Dd, N>(phi, U, rep, kappa, F);
  const auto& gen = generators<N>();
  auto totalS = [&]() { return gauge_action<Dd, N>(U, beta) + scalar_action<Dd, N>(phi, U, rep, kappa, lambda); };
  Real wl = 0.0;
  std::int64_t ls[] = {0, lat.vol - 1};
  for (std::int64_t s : ls)
    for (int mu = 0; mu < Dd; ++mu)
      for (int a = 0; a < n_gen<N>(); ++a) {
        const Cmat<N> o = U(s, mu);
        U(s, mu) = expi<N>(gen.T[a] * Complex(eps, 0)) * o;  Real Sp = totalS();
        U(s, mu) = expi<N>(gen.T[a] * Complex(-eps, 0)) * o; Real Sm = totalS();
        U(s, mu) = o;
        wl = std::max(wl, std::fabs((Sp - Sm) / (2 * eps) - F(s, mu)[a]));
      }
  std::snprintf(m, sizeof m, "link force FD (%s)", tag);
  CHECK_CLOSE(wl, 0.0, 1e-5, m);
}

// The per-vector fast path (rotate/rotate_dag/hop_link_g) must match the full
// rep_matrix path bit-for-bit (same math, just no materialized matrix).
template <int N>
static void test_fastpath(const GeneralRep<N>& rep, std::uint64_t seed, const char* tag) {
  std::mt19937_64 rng(seed);
  std::normal_distribution<Real> g(0.0, 1.0);
  Real wr = 0, wrd = 0, wg = 0;
  const auto& f = structure_constants<N>();  // (unused) keep includes warm
  (void)f;
  for (int t = 0; t < 6; ++t) {
    Cmat<N> U = rand_su<N>(rng);
    DMat M = rep.rep_matrix(U);
    DVec phi(rep.d), psi(rep.d);
    for (int k = 0; k < rep.d; ++k) {
      phi(k) = Complex(g(rng), rep.real ? 0.0 : g(rng));
      psi(k) = Complex(g(rng), rep.real ? 0.0 : g(rng));
    }
    DVec r1 = rep.rotate(U, phi),      r2 = M * phi;
    DVec d1 = rep.rotate_dag(U, phi),  d2 = M.dagger() * phi;
    for (int k = 0; k < rep.d; ++k) {
      wr  = std::max(wr,  std::abs(r1(k) - r2(k)));
      wrd = std::max(wrd, std::abs(d1(k) - d2(k)));
    }
    AlgVec<N> ga = rep.hop_link_g(U, psi, phi);
    DVec Mphi = M * phi;
    for (int a = 0; a < n_gen<N>(); ++a) {
      const Real ref = -2.0 * dot(psi, rep.T[a] * Mphi).imag();
      wg = std::max(wg, std::fabs(ga[a] - ref));
    }
  }
  char m[96];
  std::snprintf(m, sizeof m, "fast-path rotate vs matrix (%s)", tag);     CHECK_CLOSE(wr,  0.0, 1e-11, m);
  std::snprintf(m, sizeof m, "fast-path rotate_dag vs matrix (%s)", tag); CHECK_CLOSE(wrd, 0.0, 1e-11, m);
  std::snprintf(m, sizeof m, "fast-path hop_link_g vs matrix (%s)", tag); CHECK_CLOSE(wg,  0.0, 1e-11, m);
}

template <int Dd, int N>
static void test_hmc(Representation<N>& rep, std::uint64_t seed, const char* tag) {
  GaugeHiggsHMC<Dd, N> hmc(cube<Dd>(4), rep, seed);
  hmc.beta = 2.3; hmc.kappa = 0.15; hmc.lambda = 0.5; hmc.tau = 1.0; hmc.nmd = 20;
  hmc.U.hot(hmc.rng, 0.3);
  hmc.phi.cold(0.4);
  const int ntraj = 300; double s = 0.0;
  for (int t = 0; t < ntraj; ++t) { hmc.trajectory(); s += std::exp(-hmc.last_dH); }
  const double mean = s / ntraj;
  char m[128];
  std::snprintf(m, sizeof m, "<exp(-dH)>~1 general (%s) got %.4f acc %.2f", tag, mean, hmc.acceptance());
  CHECK(std::fabs(mean - 1.0) < 0.15, m);
}

int main() {
  // ---- SU(2): fundamental [1], spin-1 [2], spin-3/2 [3] ----
  GeneralRep<2> g2_f({1});
  GeneralRep<2> g2_adj({2}, GeneralRep<2>::RealType::Real);  // spin 1 = adjoint, real
  GeneralRep<2> g2_32({3});
  std::printf("-- SU(2) dims: [1]=%d [2]=%d [3]=%d --\n", g2_f.d, g2_adj.d, g2_32.d);
  test_invariants<2>(g2_f,  2, 0.5, 0.75, "SU2[1]");
  test_invariants<2>(g2_adj, 3, 2.0, 2.0,  "SU2[2]");
  test_invariants<2>(g2_32, 4, 5.0, 3.75, "SU2[3]");
  test_commutator<2>(g2_f, "SU2[1]"); test_commutator<2>(g2_adj, "SU2[2]"); test_commutator<2>(g2_32, "SU2[3]");
  test_homomorphism<2>(g2_adj, 11, "SU2[2]"); test_homomorphism<2>(g2_32, 12, "SU2[3]");

  // ---- SU(3): 3bar [1,1], sextet [2], octet [2,1] ----
  GeneralRep<3> g3_3bar({1, 1});
  GeneralRep<3> g3_6({2});
  GeneralRep<3> g3_8({2, 1}, GeneralRep<3>::RealType::Real);  // adjoint, real
  std::printf("-- SU(3) dims: [1,1]=%d [2]=%d [2,1]=%d --\n", g3_3bar.d, g3_6.d, g3_8.d);
  test_invariants<3>(g3_3bar, 3, 0.5,  4.0 / 3,  "SU3[1,1]");
  test_invariants<3>(g3_6,    6, 2.5,  10.0 / 3, "SU3[2]");
  test_invariants<3>(g3_8,    8, 3.0,  3.0,      "SU3[2,1]");
  test_commutator<3>(g3_3bar, "SU3[1,1]"); test_commutator<3>(g3_6, "SU3[2]"); test_commutator<3>(g3_8, "SU3[2,1]");

  // N-ality
  CHECK(g3_3bar.nality == 2, "SU3[1,1] N-ality 2");
  CHECK(g3_6.nality == 2, "SU3[2] N-ality 2");
  CHECK(g3_8.nality == 0, "SU3[2,1] N-ality 0 (adjoint)");

  // ---- character cross-checks (basis-independent) ----
  std::printf("-- character cross-checks --\n");
  FundamentalRep<2> f2; AdjointRep<2> a2; AdjointRep<3> a3;
  test_character<2>(g2_f,   [](const Cmat<2>& U) { return U.trace(); },               21, "SU2[1] vs fund");
  test_character<2>(g2_adj, [&](const Cmat<2>& U) { return a2.rep_matrix(U).trace(); }, 22, "SU2[2] vs adj");
  test_character<3>(g3_3bar,[](const Cmat<3>& U) { return std::conj(U.trace()); },     23, "SU3[1,1] vs conj-fund");
  test_character<3>(g3_8,   [&](const Cmat<3>& U) { return a3.rep_matrix(U).trace(); }, 24, "SU3[2,1] vs adj");

  // ---- forces + HMC with general reps ----
  std::printf("-- fast-path (per-vector apply) consistency --\n");
  test_fastpath<2>(g2_adj, 61, "SU2[2]"); test_fastpath<2>(g2_32, 62, "SU2[3]");
  test_fastpath<3>(g3_6, 63, "SU3[2]");   test_fastpath<3>(g3_8, 64, "SU3[2,1]");
  std::printf("-- force finite differences (general reps) --\n");
  test_force_fd<3, 2>(g2_adj, 31, "SU2[2] real");
  test_force_fd<3, 3>(g3_6,   32, "SU3[2] complex");
  std::printf("-- combined HMC <exp(-dH)> (general reps) --\n");
  test_hmc<3, 2>(g2_adj, 41, "SU2[2] real");
  test_hmc<3, 3>(g3_6,   42, "SU3[2] complex");

  return report("test_general");
}
