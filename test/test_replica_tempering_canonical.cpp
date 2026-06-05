// DECISIVE CANONICAL-CORRECTNESS GATE for src/hmc/replica_tempering.hpp.
//
// Parallel tempering is canonically correct iff, for EACH rung beta_k, the
// marginal distribution sampled by the replica fixed at beta_k (through which
// configurations migrate via swaps) is exactly the single-beta Boltzmann
// distribution pi_k(U) ~ exp(-S_g(U;beta_k)).  Equivalently: any observable's
// expectation under the tempered beta_k marginal must equal its expectation
// under an INDEPENDENT single-beta HMC at the same beta_k.  If the swap
// acceptance normalization (the P() energy or the sign of Delta) is wrong, the
// detailed-balance condition for the joint product measure is violated and the
// per-rung marginals are DISTORTED -- this test catches exactly that.
//
// Test system (cheap, fast-mixing): pure-gauge SU(2), rep adj, kappa=0, L=4.
// Two betas: beta_a = 1.5, beta_b = 2.0.
//
//   (1) INDEPENDENT single-beta HMC at beta_a and at beta_b; histogram
//       avg_plaquette -> jackknife mean & standard error, and variance.
//   (2) 2-replica ReplicaTempering({1.5,2.0}), enabled=TRUE, SAME total work;
//       histogram avg_plaquette of the replica FIXED at beta_a (and at beta_b).
//   (3) GATE: tempered mean vs independent mean must agree within ~3 sigma
//       (combined standard error), variances comparable.
//   (4) TOGGLE GATE: enabled=FALSE must reproduce the independent runs exactly
//       (no swaps => replica k is literally an independent single-beta stream;
//       with matched seeds the time series are bit-identical).
//
// Build (D=4, N=2):
//   g++ -std=c++20 -O3 -march=native -funroll-loops -fopenmp -Isrc -DNDIM=4 -DNCOL=2
//       test/test_replica_tempering_canonical.cpp -o build/test_replica_tempering_canonical
#include "hmc/replica_tempering.hpp"
#include "rep/rep_adjoint.hpp"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <array>

using namespace gh;

// --- jackknife mean & standard error of a time series, with binning to tame
//     autocorrelation (block average over `bin` consecutive samples first). ---
struct Stat { double mean, sem, var; std::size_t nbin; };

static Stat jackknife(const std::vector<double>& x, std::size_t bin) {
  const std::size_t nb = x.size() / bin;             // number of bins
  std::vector<double> b(nb, 0.0);
  for (std::size_t i = 0; i < nb; ++i) {
    double s = 0.0;
    for (std::size_t j = 0; j < bin; ++j) s += x[i * bin + j];
    b[i] = s / double(bin);
  }
  double sum = 0.0;
  for (double v : b) sum += v;
  const double mean = sum / double(nb);
  // raw variance of the (binned) samples, for the "variances comparable" check
  double var = 0.0;
  for (double v : b) var += (v - mean) * (v - mean);
  var /= double(nb > 1 ? nb - 1 : 1);
  // jackknife: leave-one-bin-out resampling
  std::vector<double> jk(nb);
  for (std::size_t i = 0; i < nb; ++i) jk[i] = (sum - b[i]) / double(nb - 1);
  double jkbar = 0.0; for (double v : jk) jkbar += v; jkbar /= double(nb);
  double jvar = 0.0;
  for (double v : jk) jvar += (v - jkbar) * (v - jkbar);
  jvar *= double(nb - 1) / double(nb);              // jackknife variance of the mean
  return Stat{mean, std::sqrt(jvar), var, nb};
}

int main() {
  constexpr int D = NDIM;   // 4
  constexpr int Nc = NCOL;  // 2
  static_assert(D == 4 && Nc == 2, "this test expects -DNDIM=4 -DNCOL=2");

  const int L = 4;
  std::array<int, D> ext{}; for (int mu = 0; mu < D; ++mu) ext[mu] = L;
  AdjointRep<Nc> rep;

  const Real beta_a = 1.5, beta_b = 2.0;   // the two betas we GATE (ladder endpoints)
  // A 2-rung {1.5,2.0} ladder gives ~0 swap acceptance on 4^4 (the plaquette
  // distributions don't overlap), so the swap would never fire and the
  // marginal test would be vacuous. We insert intermediate rungs so adjacent
  // pairs overlap and configurations actually MIGRATE from beta_a up to beta_b
  // and back. The endpoints stay exactly at beta_a and beta_b, so reading the
  // marginal off rung 0 and the top rung is still a clean canonical-correctness
  // gate -- but now it genuinely depends on the swap being correct, because the
  // configs sampled at each endpoint arrived there through accepted swaps.
  const std::vector<Real> ladder{1.5, 1.6, 1.7, 1.8, 1.9, 2.0};
  const std::size_t rung_a = 0;                       // beta_a = 1.5
  const std::size_t rung_b = ladder.size() - 1;       // beta_b = 2.0

  // Sampling budget. Pure-gauge SU(2) on 4^4 mixes fast; we use a generous
  // sample count so the standard errors are small enough to be a real gate.
  const int nmd     = 10;
  const int ntherm  = 200;
  const int nmeas   = 3000;    // measurements (trajectories) per stream
  const std::size_t bin = 30;  // binning to absorb autocorrelation (-> 100 bins)

  // ---------------------------------------------------------------------------
  // (1) INDEPENDENT single-beta HMC at beta_a and beta_b.
  // ---------------------------------------------------------------------------
  auto run_independent = [&](Real beta, std::uint64_t seed) {
    GaugeHiggsHMC<D, Nc> hmc(ext, rep, seed);
    hmc.beta = beta; hmc.kappa = 0.0; hmc.nmd = nmd; hmc.tau = 1.0;
    hmc.U.hot(hmc.rng, 0.9);
    for (int t = 0; t < ntherm; ++t) hmc.trajectory();
    std::vector<double> series; series.reserve(nmeas);
    for (int t = 0; t < nmeas; ++t) { hmc.trajectory(); series.push_back(avg_plaquette<D, Nc>(hmc.U)); }
    std::fprintf(stderr, "[indep beta=%.2f done]\n", (double)beta);
    return series;
  };

  std::vector<double> ind_a = run_independent(beta_a, 0xA1A1ull);
  std::vector<double> ind_b = run_independent(beta_b, 0xB2B2ull);
  Stat sa = jackknife(ind_a, bin);
  Stat sb = jackknife(ind_b, bin);

  // ---------------------------------------------------------------------------
  // (2) 2-replica tempering, enabled=TRUE, same total work; record the FIXED-
  //     beta marginal (replica 0 -> beta_a, replica 1 -> beta_b).
  // ---------------------------------------------------------------------------
  ReplicaTempering<D, Nc> pt(ext, rep, ladder, /*seed0=*/0x7E11ull);
  pt.enabled = true; pt.n_sweep = 1;
  pt.set_kappa(0.0); pt.set_nmd(nmd); pt.set_tau(1.0);
  for (std::size_t k = 0; k < pt.n_replicas(); ++k) pt.replica(k).U.hot(pt.replica(k).rng, 0.9);
  for (int t = 0; t < ntherm; ++t) pt.step();
  std::vector<double> tmp_a; tmp_a.reserve(nmeas);
  std::vector<double> tmp_b; tmp_b.reserve(nmeas);
  for (int t = 0; t < nmeas; ++t) {
    pt.step();
    tmp_a.push_back(pt.avg_plaq(rung_a));   // marginal at fixed beta_a (rung 0)
    tmp_b.push_back(pt.avg_plaq(rung_b));   // marginal at fixed beta_b (top rung)
  }
  std::fprintf(stderr, "[tempering done]\n");
  Stat ta = jackknife(tmp_a, bin);
  Stat tb = jackknife(tmp_b, bin);

  // sigma of disagreement (combined standard error of the two means)
  auto nsig = [](const Stat& a, const Stat& b) {
    return std::fabs(a.mean - b.mean) / std::sqrt(a.sem * a.sem + b.sem * b.sem);
  };
  const double sig_a = nsig(sa, ta);
  const double sig_b = nsig(sb, tb);
  // variance ratio (binned-sample variance); "comparable" => within ~2x either way
  const double vr_a = ta.var / sa.var;
  const double vr_b = tb.var / sb.var;

  std::printf("# pure-gauge SU(2) adj, kappa=0, L=%d, nmd=%d, nmeas=%d, bin=%zu\n",
              L, nmd, nmeas, bin);
  std::printf("# observable: avg_plaquette\n\n");
  std::printf("%-7s %-12s %-12s %-12s %-12s %-9s %-9s\n",
              "beta", "indep_mean", "indep_sem", "temp_mean", "temp_sem", "agree(s)", "var_ratio");
  std::printf("%-7.3f %-12.7f %-12.7f %-12.7f %-12.7f %-9.3f %-9.4f\n",
              (double)beta_a, sa.mean, sa.sem, ta.mean, ta.sem, sig_a, vr_a);
  std::printf("%-7.3f %-12.7f %-12.7f %-12.7f %-12.7f %-9.3f %-9.4f\n",
              (double)beta_b, sb.mean, sb.sem, tb.mean, tb.sem, sig_b, vr_b);

  std::printf("\n# swap acceptance (chain must be CONNECTED: every pair accepts):\n");
  bool chain_connected = true;
  for (std::size_t p = 0; p < pt.n_pairs(); ++p) {
    std::printf("  pair (%zu,%zu)  attempts=%llu  accepts=%llu  acc=%.4f\n",
                p, p + 1, (unsigned long long)pt.swap_attempts[p],
                (unsigned long long)pt.swap_accepts[p], pt.pair_acceptance(p));
    if (pt.swap_attempts[p] == 0 || pt.swap_accepts[p] == 0) chain_connected = false;
  }

  // ---------------------------------------------------------------------------
  // (4) TOGGLE GATE: enabled=FALSE => replica k is an independent single-beta
  //     stream. Seed each replica EXACTLY like a mirrored independent run at the
  //     same beta and demand the plaquette time series MATCH TO ROUND-OFF. (We
  //     allow a tiny tolerance, not bit-equality: avg_plaquette uses an OpenMP
  //     reduction whose summation order is not reproducible across separate
  //     runs, so two algebraically identical configs can differ at ~1e-15. The
  //     RNG streams and the dynamics are identical -- only the measurement's
  //     float accumulation order differs -- so a few-ULP tolerance is the
  //     correct no-op criterion.)
  // ---------------------------------------------------------------------------
  const std::uint64_t seedR0   = 0xC0FFEEull;
  const std::uint64_t seedStr  = 1000003ull;
  auto run_independent_short = [&](Real beta, std::uint64_t seed, int n) {
    GaugeHiggsHMC<D, Nc> hmc(ext, rep, seed);
    hmc.beta = beta; hmc.kappa = 0.0; hmc.nmd = nmd; hmc.tau = 1.0;
    hmc.U.hot(hmc.rng, 0.9);
    std::vector<double> s; s.reserve(n);
    for (int t = 0; t < n; ++t) { hmc.trajectory(); s.push_back(avg_plaquette<D, Nc>(hmc.U)); }
    return s;
  };
  const int ncheck = 120;
  // Reference: an independent stream for EACH rung, seeded exactly as the
  // toggle-off ReplicaTempering seeds that rung (seedR0 + k*seedStr, beta=ladder[k]).
  const std::size_t M = ladder.size();
  std::vector<std::vector<double>> ref(M);
  for (std::size_t k = 0; k < M; ++k)
    ref[k] = run_independent_short(ladder[k], seedR0 + k * seedStr, ncheck);

  ReplicaTempering<D, Nc> pt_off(ext, rep, ladder, /*seed0=*/seedR0, /*seed_stride=*/seedStr);
  pt_off.enabled = false; pt_off.n_sweep = 1;
  pt_off.set_kappa(0.0); pt_off.set_nmd(nmd); pt_off.set_tau(1.0);
  for (std::size_t k = 0; k < pt_off.n_replicas(); ++k) pt_off.replica(k).U.hot(pt_off.replica(k).rng, 0.9);
  std::vector<std::vector<double>> off(M);
  for (int t = 0; t < ncheck; ++t) {
    pt_off.step();
    for (std::size_t k = 0; k < M; ++k) off[k].push_back(pt_off.avg_plaq(k));
  }
  const double TOGGLE_TOL = 1e-12;   // >> the ~1e-15 OpenMP reduction reorder noise
  double max_abs_dev = 0.0;
  for (std::size_t k = 0; k < M; ++k)
    for (int t = 0; t < ncheck; ++t)
      max_abs_dev = std::max(max_abs_dev, std::fabs(off[k][t] - ref[k][t]));
  const bool toggle_matches = (max_abs_dev < TOGGLE_TOL);
  std::size_t off_attempts = 0;
  for (std::size_t p = 0; p < pt_off.n_pairs(); ++p) off_attempts += pt_off.swap_attempts[p];

  std::printf("\n# toggle-off gate (enabled=false): every rung == its mirrored independent run\n");
  std::printf("  swap_attempts=%zu (expect 0)\n", off_attempts);
  std::printf("  max|dev| over all %zu rungs x %d traj = %.3e  (tol %.0e) -> %s\n",
              M, ncheck, max_abs_dev, TOGGLE_TOL, toggle_matches ? "MATCH" : "MISMATCH");

  // ---------------------------------------------------------------------------
  // VERDICT
  // ---------------------------------------------------------------------------
  const double SIG_GATE = 3.0;
  const bool marg_a_ok = (sig_a <= SIG_GATE);
  const bool marg_b_ok = (sig_b <= SIG_GATE);
  const bool var_a_ok  = (vr_a > 0.5 && vr_a < 2.0);
  const bool var_b_ok  = (vr_b > 0.5 && vr_b < 2.0);
  const bool swaps_ok  = chain_connected;   // every adjacent pair accepted -> configs migrate endpoint-to-endpoint
  const bool toggle_ok = (off_attempts == 0) && toggle_matches;

  std::printf("\n# gates:\n");
  std::printf("  marginal beta_a within %.0f sigma : %s (%.2f s)\n", SIG_GATE, marg_a_ok ? "PASS" : "FAIL", sig_a);
  std::printf("  marginal beta_b within %.0f sigma : %s (%.2f s)\n", SIG_GATE, marg_b_ok ? "PASS" : "FAIL", sig_b);
  std::printf("  variance comparable beta_a       : %s (ratio %.3f)\n", var_a_ok ? "PASS" : "FAIL", vr_a);
  std::printf("  variance comparable beta_b       : %s (ratio %.3f)\n", var_b_ok ? "PASS" : "FAIL", vr_b);
  std::printf("  swap chain connected (all accept): %s\n", swaps_ok ? "PASS" : "FAIL");
  std::printf("  toggle-off no-op (to round-off)  : %s\n", toggle_ok ? "PASS" : "FAIL");

  const bool pass = marg_a_ok && marg_b_ok && var_a_ok && var_b_ok && swaps_ok && toggle_ok;
  std::printf("\n%s\n", pass ? "PASS" : "FAIL");
  return pass ? 0 : 1;
}
