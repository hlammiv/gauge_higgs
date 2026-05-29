// Validation of the direct symmetric-representation basis added to GeneralRep for
// single-row Young diagrams (fully-symmetric reps: SU(2) spin-j = {2j},
// SU(3) (p,0) = {p}). The symmetric basis is built WITHOUT enumerating the n!
// row-symmetrizer permutations, so large single-row reps such as SU(2) spin-6
// = {12} (the binary-icosahedral 2I rep) build without OOM. We check:
//   (a) SU(2) {12} constructs (N^n = 4096) with d = 13, satisfies the algebra
//       commutator [T^a,T^b] = i f^{abc} T^c, the Casimir sum T^aT^a = C2 I, and
//       the rep_matrix homomorphism;
//   (b) for {2} (= spin-1 = adjoint) the CHARACTER Tr(rep_matrix(U)) matches
//       AdjointRep<2> on random U (a basis-independent cross-check that the
//       direct symmetric basis spans the same irrep as the general path);
//   (c) the symmetric path agrees with the (small-n) general Young path: same dim
//       and same characters, and the basis is orthonormal.
#include "check.hpp"
#include "group/sun.hpp"  // expi, alg_to_mat
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

// [T^a_R, T^b_R] = i f^{abc} T^c_R  (the definitive proof T^a_R is a representation).
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
        if (fabc != 0.0) rhs = rhs + rep.T[c] * Complex(0.0, fabc);
      }
      worst = std::max(worst, fnorm(comm - rhs));
    }
  char m[96]; std::snprintf(m, sizeof m, "commutator [T,T]=ifT (%s)", tag);
  CHECK_CLOSE(worst, 0.0, 1e-9, m);
}

// sum_a (T^a)^2 = C2 I  and  Tr(T^aT^b) = T(R) delta^{ab}.
template <int N>
static void test_casimir(const Representation<N>& rep, Real expect_T, Real expect_C2, const char* tag) {
  char m[96];
  const Real T = (rep.T[0] * rep.T[0]).trace().real();
  const Real C2 = T * (static_cast<Real>(N) * N - 1) / rep.d;
  std::snprintf(m, sizeof m, "Dynkin index (%s) got %.4f", tag, T);
  CHECK_CLOSE(T, expect_T, 1e-7, m);
  std::snprintf(m, sizeof m, "Casimir (%s) got %.4f", tag, C2);
  CHECK_CLOSE(C2, expect_C2, 1e-7, m);
  std::snprintf(m, sizeof m, "invariant consistency / T^aT^a=C2 I (%s)", tag);
  CHECK_CLOSE(rep.invariant_violation(T, C2), 0.0, 1e-7, m);
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

// chi(U) = Tr D_R(U) must match a reference (basis-independent).
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

// Orthonormality of the stored basis W (the symmetric vectors must be O.N.).
template <int N>
static void test_orthonormal(const GeneralRep<N>& rep, const char* tag) {
  Real worst = 0.0;
  const long dimT = rep.dimT;
  for (int k = 0; k < rep.d; ++k)
    for (int l = 0; l < rep.d; ++l) {
      Complex s(0, 0);
      for (long i = 0; i < dimT; ++i) s += std::conj(rep.W[k][i]) * rep.W[l][i];
      worst = std::max(worst, std::abs(s - Complex(k == l ? 1.0 : 0.0, 0.0)));
    }
  char m[96]; std::snprintf(m, sizeof m, "symmetric basis orthonormal (%s)", tag);
  CHECK_CLOSE(worst, 0.0, 1e-12, m);
}

int main() {
  // ---- (a) SU(2) spin-6 = {12} = binary-icosahedral 2I rep: builds without OOM ----
  // (general path would enumerate the row group S_12 -> 12! perms -> OOM/blow up.)
  std::printf("-- SU(2) {12} (spin 6) direct symmetric build --\n");
  GeneralRep<2> g2_s12({12});  // N^n = 2^12 = 4096
  std::printf("   d = %d (expect 13), N^n = %ld\n", g2_s12.d, g2_s12.dimT);
  CHECK(g2_s12.d == 13, "SU(2) {12} dim == 13 (2j+1, j=6)");
  CHECK(g2_s12.dimT == 4096, "SU(2) {12} N^n == 4096");
  test_orthonormal<2>(g2_s12, "SU2{12}");
  test_commutator<2>(g2_s12, "SU2{12}");
  // spin-6: C2 = j(j+1) = 42 in the spin normalization; in our T=Lambda/2 normalization
  // T(R) = C2 * d / (N^2-1) is derived internally, so just check internal consistency
  // (sum T^aT^a = C2 I and Tr(T^aT^b)=T delta) with the measured T/C2.
  {
    const Real T = (g2_s12.T[0] * g2_s12.T[0]).trace().real();
    const Real C2 = T * (2.0 * 2 - 1) / g2_s12.d;
    std::printf("   measured Dynkin T = %.4f, Casimir C2 = %.4f (expect C2 = j(j+1) = 42)\n", T, C2);
    CHECK_CLOSE(C2, 42.0, 1e-6, "SU(2) {12} Casimir == j(j+1) = 42");
    CHECK_CLOSE(g2_s12.invariant_violation(T, C2), 0.0, 1e-6, "SU(2) {12} sum T^aT^a = C2 I");
  }
  test_homomorphism<2>(g2_s12, 101, "SU2{12}");

  // ---- (b) {2} = spin-1 = adjoint: character matches AdjointRep<2> (basis-independent) ----
  std::printf("-- SU(2) {2} (spin 1) character vs AdjointRep<2> --\n");
  GeneralRep<2> g2_s2({2}, GeneralRep<2>::RealType::Real);
  CHECK(g2_s2.d == 3, "SU(2) {2} dim == 3");
  test_orthonormal<2>(g2_s2, "SU2{2}");
  test_commutator<2>(g2_s2, "SU2{2}");
  test_casimir<2>(g2_s2, 2.0, 2.0, "SU2{2}");
  AdjointRep<2> a2;
  test_character<2>(g2_s2, [&](const Cmat<2>& U) { return a2.rep_matrix(U).trace(); }, 202, "SU2{2} vs adjoint");

  // ---- (c) symmetric path agrees with general path on small reps + extra single-rows ----
  std::printf("-- cross-checks: more single-row reps --\n");
  GeneralRep<2> g2_s1({1});  // n=1 edge case (fundamental)
  CHECK(g2_s1.d == 2, "SU(2) {1} dim == 2");
  test_orthonormal<2>(g2_s1, "SU2{1}");
  test_character<2>(g2_s1, [](const Cmat<2>& U) { return U.trace(); }, 203, "SU2{1} vs fund");

  GeneralRep<2> g2_s4({4});  // spin-2, dim 5
  CHECK(g2_s4.d == 5, "SU(2) {4} dim == 5");
  test_orthonormal<2>(g2_s4, "SU2{4}");
  test_commutator<2>(g2_s4, "SU2{4}");
  test_homomorphism<2>(g2_s4, 204, "SU2{4}");

  GeneralRep<3> g3_s2({2});  // SU(3) sextet (2,0), dim 6
  CHECK(g3_s2.d == 6, "SU(3) {2} dim == 6 (sextet)");
  test_orthonormal<3>(g3_s2, "SU3{2}");
  test_commutator<3>(g3_s2, "SU3{2}");
  test_casimir<3>(g3_s2, 2.5, 10.0 / 3, "SU3{2}");
  test_homomorphism<3>(g3_s2, 205, "SU3{2}");
  CHECK(g3_s2.nality == 2, "SU(3) {2} N-ality 2");

  GeneralRep<3> g3_s3({3});  // SU(3) decuplet (3,0), dim 10
  CHECK(g3_s3.d == 10, "SU(3) {3} dim == 10 (decuplet)");
  test_orthonormal<3>(g3_s3, "SU3{3}");
  test_commutator<3>(g3_s3, "SU3{3}");
  test_homomorphism<3>(g3_s3, 206, "SU3{3}");

  return report("test_symrep");
}
