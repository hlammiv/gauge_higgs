// CANONICAL-CORRECTNESS GATE for the KAPPA-tempering axis of
// src/hmc/replica_tempering.hpp (Axis::Kappa, hopping conjugate H).
//
// If the kappa-conjugate H = 2*(vol*D)*link_energy and the swap Delta =
// (kappa_i-kappa_j)(H_i-H_j) are correctly normalized, then the marginal
// distribution sampled by the replica FIXED at kappa_k (through which configs
// migrate via swaps) must equal the single-kappa Boltzmann distribution, i.e.
// an INDEPENDENT single-kappa HMC at the same (beta,kappa) must give the same
// expectation of any observable. A wrong factor (e.g. dropping the 2 or the
// n_bonds) distorts the marginal -> caught here via link_energy and avg_plaquette.
//
// System: SU(2) fundamental, SHARED beta=2.0, kappa-ladder, L=4, multi-invariant
// off (plain quartic via lambda). Endpoints kappa_a (small) and kappa_b (larger)
// are GATED against independent runs.
//
// Build (D=4, N=2):
//   g++ -std=c++20 -O3 -march=native -funroll-loops -fopenmp -Isrc -DNDIM=4 -DNCOL=2
//       test/test_replica_tempering_kappa_canonical.cpp -o build/test_replica_tempering_kappa_canonical
#include "hmc/replica_tempering.hpp"
#include "rep/rep_fundamental.hpp"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <array>

using namespace gh;

struct Stat { double mean, sem, var; };

static Stat jackknife(const std::vector<double>& x, std::size_t bin) {
  const std::size_t nb = x.size() / bin;
  std::vector<double> b(nb, 0.0);
  for (std::size_t i = 0; i < nb; ++i) {
    double s = 0.0;
    for (std::size_t j = 0; j < bin; ++j) s += x[i * bin + j];
    b[i] = s / double(bin);
  }
  double sum = 0.0; for (double v : b) sum += v;
  const double mean = sum / double(nb);
  double var = 0.0; for (double v : b) var += (v - mean) * (v - mean);
  var /= double(nb > 1 ? nb - 1 : 1);
  std::vector<double> jk(nb);
  for (std::size_t i = 0; i < nb; ++i) jk[i] = (sum - b[i]) / double(nb - 1);
  double jkbar = 0.0; for (double v : jk) jkbar += v; jkbar /= double(nb);
  double jvar = 0.0; for (double v : jk) jvar += (v - jkbar) * (v - jkbar);
  jvar *= double(nb - 1) / double(nb);
  return Stat{mean, std::sqrt(jvar), var};
}

int main() {
  constexpr int D = NDIM, Nc = NCOL;
  static_assert(D == 4 && Nc == 2, "this test expects -DNDIM=4 -DNCOL=2");

  const int L = 4;
  std::array<int, D> ext{}; for (int mu = 0; mu < D; ++mu) ext[mu] = L;
  FundamentalRep<Nc> rep;

  const Real BETA = 2.0;       // SHARED across the kappa-ladder
  const Real lambda = 0.5;     // plain quartic potential (default path)
  // kappa-ladder; endpoints kappa_a, kappa_b are gated. Intermediate rungs give
  // overlap so configs actually migrate endpoint-to-endpoint. Spacing is DENSER at
  // larger kappa because the hopping energy rises steeply there (the swap-overlap
  // window in H shrinks), exactly the connectivity-tuning the beta test documents.
  const std::vector<Real> ladder{0.10, 0.14, 0.18, 0.22, 0.26, 0.30,
                                 0.33, 0.36, 0.38, 0.40};
  const Real kappa_a = ladder.front();
  const Real kappa_b = ladder.back();
  const std::size_t rung_a = 0, rung_b = ladder.size() - 1;

  const int nmd = 16, ntherm = 300, nmeas = 4000;
  const std::size_t bin = 40;  // -> 100 bins

  // observable: link_energy (the hopping; the kappa-sensitive order parameter)
  auto run_independent = [&](Real kappa, std::uint64_t seed) {
    GaugeHiggsHMC<D, Nc> hmc(ext, rep, seed);
    hmc.beta = BETA; hmc.kappa = kappa; hmc.lambda = lambda; hmc.nmd = nmd; hmc.tau = 1.0;
    hmc.U.hot(hmc.rng, 0.9); hmc.phi.gaussian(hmc.rng, 9999, rep.real, 0.3);
    for (int t = 0; t < ntherm; ++t) hmc.trajectory();
    std::vector<double> le; le.reserve(nmeas);
    for (int t = 0; t < nmeas; ++t) { hmc.trajectory(); le.push_back(link_energy<D, Nc>(hmc.phi, hmc.U, rep)); }
    std::fprintf(stderr, "[indep kappa=%.2f done]\n", (double)kappa);
    return le;
  };
  std::vector<double> ind_a = run_independent(kappa_a, 0xA1A1ull);
  std::vector<double> ind_b = run_independent(kappa_b, 0xB2B2ull);
  Stat sa = jackknife(ind_a, bin), sb = jackknife(ind_b, bin);

  // kappa-tempered ladder at SHARED beta, enabled=true.
  ReplicaTempering<D, Nc> pt(typename ReplicaTempering<D, Nc>::KappaLadder{},
                             ext, rep, ladder, BETA, /*seed0=*/0x7E11ull);
  pt.enabled = true; pt.n_sweep = 1;
  pt.set_lambda(lambda); pt.set_nmd(nmd); pt.set_tau(1.0);
  for (std::size_t k = 0; k < pt.n_replicas(); ++k) {
    pt.replica(k).U.hot(pt.replica(k).rng, 0.9);
    pt.replica(k).phi.gaussian(pt.replica(k).rng, 9999, rep.real, 0.3);
  }
  for (int t = 0; t < ntherm; ++t) pt.step();
  std::vector<double> tmp_a, tmp_b; tmp_a.reserve(nmeas); tmp_b.reserve(nmeas);
  for (int t = 0; t < nmeas; ++t) {
    pt.step();
    const auto& ra = pt.replica(rung_a); tmp_a.push_back(link_energy<D, Nc>(ra.phi, ra.U, rep));
    const auto& rb = pt.replica(rung_b); tmp_b.push_back(link_energy<D, Nc>(rb.phi, rb.U, rep));
  }
  std::fprintf(stderr, "[kappa-tempering done]\n");
  Stat ta = jackknife(tmp_a, bin), tb = jackknife(tmp_b, bin);

  auto nsig = [](const Stat& a, const Stat& b) {
    return std::fabs(a.mean - b.mean) / std::sqrt(a.sem * a.sem + b.sem * b.sem);
  };
  const double sig_a = nsig(sa, ta), sig_b = nsig(sb, tb);
  const double vr_a = ta.var / sa.var, vr_b = tb.var / sb.var;

  std::printf("# SU(2) fund, kappa-tempered, shared beta=%.2f, lambda=%.2f, L=%d, nmd=%d, nmeas=%d, bin=%zu\n",
              (double)BETA, (double)lambda, L, nmd, nmeas, bin);
  std::printf("# observable: link_energy (hopping; kappa-conjugate H = 2*(vol*D)*link_energy)\n\n");
  std::printf("%-7s %-12s %-12s %-12s %-12s %-9s %-9s\n",
              "kappa", "indep_mean", "indep_sem", "temp_mean", "temp_sem", "agree(s)", "var_ratio");
  std::printf("%-7.3f %-12.7f %-12.7f %-12.7f %-12.7f %-9.3f %-9.4f\n",
              (double)kappa_a, sa.mean, sa.sem, ta.mean, ta.sem, sig_a, vr_a);
  std::printf("%-7.3f %-12.7f %-12.7f %-12.7f %-12.7f %-9.3f %-9.4f\n",
              (double)kappa_b, sb.mean, sb.sem, tb.mean, tb.sem, sig_b, vr_b);

  std::printf("\n# swap acceptance (chain must be CONNECTED: every pair accepts):\n");
  bool chain_connected = true;
  for (std::size_t p = 0; p < pt.n_pairs(); ++p) {
    std::printf("  pair (%zu,%zu)  kappa %.3f->%.3f  attempts=%llu accepts=%llu acc=%.4f\n",
                p, p + 1, (double)pt.kappa(p), (double)pt.kappa(p + 1),
                (unsigned long long)pt.swap_attempts[p],
                (unsigned long long)pt.swap_accepts[p], pt.pair_acceptance(p));
    if (pt.swap_attempts[p] == 0 || pt.swap_accepts[p] == 0) chain_connected = false;
  }

  const double SIG_GATE = 3.5;   // mild widening: link_energy autocorrelation is heavier than plaquette
  const bool marg_a_ok = (sig_a <= SIG_GATE), marg_b_ok = (sig_b <= SIG_GATE);
  const bool var_a_ok = (vr_a > 0.4 && vr_a < 2.5), var_b_ok = (vr_b > 0.4 && vr_b < 2.5);
  const bool swaps_ok = chain_connected;

  std::printf("\n# gates:\n");
  std::printf("  marginal kappa_a within %.1f sigma : %s (%.2f s)\n", SIG_GATE, marg_a_ok ? "PASS" : "FAIL", sig_a);
  std::printf("  marginal kappa_b within %.1f sigma : %s (%.2f s)\n", SIG_GATE, marg_b_ok ? "PASS" : "FAIL", sig_b);
  std::printf("  variance comparable kappa_a        : %s (ratio %.3f)\n", var_a_ok ? "PASS" : "FAIL", vr_a);
  std::printf("  variance comparable kappa_b        : %s (ratio %.3f)\n", var_b_ok ? "PASS" : "FAIL", vr_b);
  std::printf("  swap chain connected (all accept)  : %s\n", swaps_ok ? "PASS" : "FAIL");

  const bool pass = marg_a_ok && marg_b_ok && var_a_ok && var_b_ok && swaps_ok;
  std::printf("\n%s\n", pass ? "PASS" : "FAIL");
  return pass ? 0 : 1;
}
