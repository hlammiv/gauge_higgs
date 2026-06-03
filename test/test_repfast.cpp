// Validate the GeneralRep exp(i w.T_R) fast path against the exact tensor path, for
// SU(2) (spin-3/4/6 = {6}/{8}/{12}) and SU(3) (to prove generality beyond N=2).
// Checks: (1) fundamental-log round-trip exp(i w.T_F)=U; (2) fast_D == tensor matrix;
// (3) rotate / rotate_dag / hop_link_g (the HMC hot path) == tensor reference.
#include "rep/rep_general.hpp"
#include "group/sun.hpp"
#include <cstdio>
#include <random>
#include <vector>

using namespace gh;

static int failures = 0;
static void check(const char* what, double err, double tol) {
  const bool ok = (err < tol);
  std::printf("    %-34s err=%.3e  %s\n", what, err, ok ? "ok" : "FAIL");
  if (!ok) ++failures;
}

template <int N>
static Cmat<N> rand_link(std::mt19937_64& rng) {
  std::uniform_real_distribution<double> U(-1.6, 1.6);
  AlgVec<N> v{};
  for (int a = 0; a < n_gen<N>(); ++a) v[a] = U(rng);
  return expi<N>(alg_to_mat<N>(v));
}

template <int N>
static double mat_err(const DMat& A, const DMat& B) { return fnorm(A - B); }

template <int N>
static void test_rep(std::vector<int> rows, const char* tag, std::mt19937_64& rng, double tol) {
  GeneralRep<N> R(rows);
  std::printf("  %s  d=%d  dimT=%ld  use_fast=%s\n", tag, R.d, R.dimT, R.use_fast ? "YES" : "no");
  std::normal_distribution<double> G(0.0, 1.0);
  double e_log = 0, e_mat = 0, e_rot = 0, e_rotd = 0, e_hop = 0;
  for (int trial = 0; trial < 20; ++trial) {
    const Cmat<N> Uf = rand_link<N>(rng);
    // (1) fundamental-log round-trip
    const Cmat<N> Ur = expi<N>(alg_to_mat<N>(fund_alg<N>(Uf)));
    e_log = std::max(e_log, (Ur - Uf).fnorm());
    // (2) fast matrix vs tensor matrix (independent of use_fast)
    const DMat Df = R.fast_D(Uf), Dt = R.rep_matrix_tensor(Uf);
    e_mat = std::max(e_mat, mat_err<N>(Df, Dt));
    // (3) hot-path ops vs tensor reference
    DVec phi(R.d), phiy(R.d);
    for (int i = 0; i < R.d; ++i) { phi(i) = Complex(G(rng), G(rng)); phiy(i) = Complex(G(rng), G(rng)); }
    const DVec ref  = Dt * phi;
    const DVec refd = Dt.dagger() * phi;
    DVec rr = R.rotate(Uf, phi), rrd = R.rotate_dag(Uf, phi);
    double a = 0, b = 0;
    for (int i = 0; i < R.d; ++i) { a += std::norm(rr(i) - ref(i)); b += std::norm(rrd(i) - refd(i)); }
    e_rot = std::max(e_rot, std::sqrt(a)); e_rotd = std::max(e_rotd, std::sqrt(b));
    // hop_link_g reference: g[a] = -2 Im( phi^dag T^a (D phiy) )
    const AlgVec<N> g = R.hop_link_g(Uf, phi, phiy);
    const DVec Dy = Dt * phiy;
    double h = 0;
    for (int aa = 0; aa < n_gen<N>(); ++aa) {
      const double gref = -2.0 * dot(phi, R.T[aa] * Dy).imag();
      h += (g[aa] - gref) * (g[aa] - gref);
    }
    e_hop = std::max(e_hop, std::sqrt(h));
  }
  check("fund-log round-trip exp(iw.T)=U", e_log, tol);
  check("fast_D vs tensor matrix", e_mat, tol);
  check("rotate vs tensor", e_rot, tol);
  check("rotate_dag vs tensor", e_rotd, tol);
  check("hop_link_g vs tensor", e_hop, tol);
}

int main() {
  std::mt19937_64 rng(20260603);
  std::printf("== SU(2): spin-3 {6}, spin-4 {8}, spin-6 {12} ==\n");
  test_rep<2>({6},  "2T spin-3 {6} ", rng, 1e-9);
  test_rep<2>({8},  "2O spin-4 {8} ", rng, 1e-9);
  test_rep<2>({12}, "2I spin-6 {12}", rng, 1e-9);
  std::printf("== SU(3): generality check (adjoint {1,1}, sym {2}, large {5}) ==\n");
  test_rep<3>({1, 1}, "SU(3) adj {1,1}", rng, 1e-8);
  test_rep<3>({2},    "SU(3) sym {2}  ", rng, 1e-8);
  test_rep<3>({5},    "SU(3) {5}      ", rng, 1e-8);
  std::printf("\n%s  (%d failures)\n", failures ? "FAILED" : "ALL PASS", failures);
  return failures ? 1 : 0;
}
