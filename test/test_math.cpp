// Stage-0 unit tests: linear algebra, su(N) generators, exp, projections.
#include "check.hpp"
#include "group/sun.hpp"
#include <random>

using namespace gh;

template <int N>
static Cmat<N> random_herm_traceless(std::mt19937_64& rng) {
  std::normal_distribution<Real> g(0.0, 1.0);
  AlgVec<N> v{};
  for (int a = 0; a < n_gen<N>(); ++a) v[a] = g(rng);
  return alg_to_mat<N>(v);
}

template <int N>
static void test_generators() {
  const auto& gen = generators<N>();
  // Tr(T^a T^b) = 1/2 delta^{ab}
  for (int a = 0; a < n_gen<N>(); ++a)
    for (int b = 0; b < n_gen<N>(); ++b) {
      Complex t = trProd(gen.T[a], gen.T[b]);
      CHECK_CLOSE(t.real(), (a == b ? 0.5 : 0.0), 1e-12, "Tr(T^a T^b)=1/2 delta");
      CHECK_CLOSE(t.imag(), 0.0, 1e-12, "Tr(T^a T^b) real");
    }
  // Hermitian and traceless.
  for (int a = 0; a < n_gen<N>(); ++a) {
    Cmat<N> d = gen.T[a].dagger();
    Cmat<N> diff = gen.T[a] - d;
    CHECK_CLOSE(diff.fnorm(), 0.0, 1e-12, "T^a Hermitian");
    CHECK_CLOSE(gen.T[a].trace().real(), 0.0, 1e-12, "T^a traceless (re)");
    CHECK_CLOSE(gen.T[a].trace().imag(), 0.0, 1e-12, "T^a traceless (im)");
  }
}

template <int N>
static void test_exp_and_proj(std::mt19937_64& rng) {
  for (int trial = 0; trial < 20; ++trial) {
    Cmat<N> H = random_herm_traceless<N>(rng);
    Cmat<N> U = expi<N>(H);
    CHECK_CLOSE(unitarity_violation<N>(U), 0.0, 1e-10, "expi unitary");
    CHECK_CLOSE(std::abs(det<N>(U)), 1.0, 1e-10, "expi |det|=1");
    Complex dU = det<N>(U);
    CHECK_CLOSE(dU.real(), 1.0, 1e-9, "expi det=1 (re)");
    CHECK_CLOSE(dU.imag(), 0.0, 1e-9, "expi det=1 (im)");
    // alg_to_mat / mat_to_alg round trip on the Hermitian traceless H.
    AlgVec<N> v = mat_to_alg<N>(H);
    Cmat<N> H2 = alg_to_mat<N>(v);
    CHECK_CLOSE((H - H2).fnorm(), 0.0, 1e-12, "alg<->mat round trip");
    // proj_TA: anti-Hermitian + traceless.
    Cmat<N> M = random_herm_traceless<N>(rng);  // any matrix
    Cmat<N> A = proj_TA<N>(M);
    Cmat<N> sum = A + A.dagger();
    CHECK_CLOSE(sum.fnorm(), 0.0, 1e-12, "proj_TA anti-Hermitian");
    CHECK_CLOSE(std::abs(A.trace()), 0.0, 1e-12, "proj_TA traceless");
  }
}

static void test_su2_closed_vs_general(std::mt19937_64& rng) {
  // Compare SU(2) closed form to a brute-force series for the same H.
  std::normal_distribution<Real> g(0.0, 1.0);
  for (int trial = 0; trial < 20; ++trial) {
    AlgVec<2> v{ g(rng), g(rng), g(rng) };
    Cmat<2> H = alg_to_mat<2>(v);
    Cmat<2> U = expi<2>(H);  // closed form (specialization)
    // Brute-force exp(iH) via long Taylor with heavy scaling.
    int s = 8;
    Cmat<2> A = H * Complex(0.0, std::pow(0.5, s));
    Cmat<2> term = Cmat<2>::identity(), res = Cmat<2>::identity();
    for (int k = 1; k <= 30; ++k) { term = term * A; term *= Complex(1.0 / k, 0); res += term; }
    for (int i = 0; i < s; ++i) res = res * res;
    CHECK_CLOSE((U - res).fnorm(), 0.0, 1e-10, "SU(2) closed == series");
  }
}

template <int N>
static void test_structure_constants() {
  const auto& f = structure_constants<N>();
  const int n = n_gen<N>();
  // antisymmetry in first two indices
  for (int a = 0; a < n; ++a)
    for (int b = 0; b < n; ++b)
      for (int c = 0; c < n; ++c) {
        Real fabc = f[(a * n + b) * n + c];
        Real fbac = f[(b * n + a) * n + c];
        CHECK_CLOSE(fabc + fbac, 0.0, 1e-10, "f antisym in a,b");
      }
}

int main() {
  std::mt19937_64 rng(12345);
  std::printf("-- generators --\n");
  test_generators<2>(); test_generators<3>(); test_generators<4>();
  std::printf("-- exp / projections --\n");
  test_exp_and_proj<2>(rng); test_exp_and_proj<3>(rng); test_exp_and_proj<4>(rng);
  std::printf("-- SU(2) closed form --\n");
  test_su2_closed_vs_general(rng);
  std::printf("-- structure constants --\n");
  // SU(2): f^{012} = +1 (epsilon) with ordering (sx/2, sy/2, sz/2).
  {
    const auto& f = structure_constants<2>();
    const int n = 3;
    CHECK_CLOSE(f[(0 * n + 1) * n + 2], 1.0, 1e-10, "SU(2) f^{012}=1");
  }
  test_structure_constants<2>(); test_structure_constants<3>();
  return report("test_math");
}
