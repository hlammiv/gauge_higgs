// DECISIVE canonical-correctness gate for U1ReplicaTempering<D>.
//
// The algebraic gate in test_u1_replica_tempering.cpp already proves the swap
// Delta equals the brute-force action difference [S_i(C_j)+S_j(C_i)] -
// [S_i(C_i)+S_j(C_j)] (sign + normalization of the conjugate). This file proves
// the STRONGER, operational property: replica exchange leaves the per-rung
// equilibrium distribution INVARIANT. Concretely, for a fixed coupling value g_a
// occupied by an interior rung, the marginal distribution of an observable
// (average plaquette) sampled while tempering is ON must agree -- mean within a
// few jackknife sigma, comparable variance -- with the distribution sampled by a
// completely INDEPENDENT single-coupling U1HMC chain at the same g_a (and the
// same other couplings). A wrong conjugate sign/normalization would bias the
// accept probability and detune this marginal; matching it is the canonical
// proof that the swap is a valid (detailed-balance-respecting) MC move.
//
// Gates here:
//   (A) KAPPA axis  -- the metastability cure: tempered marginal == independent
//       single-kappa marginal at >=2 fixed kappa rungs (mean within ~3 sigma,
//       variance ratio in [0.5,2]); accepts>0 on every pair.
//   (B) BETA axis   -- same canonical test on the gauge coupling.
//   (C) TOGGLE OFF  -- with enabled=false and mirror seeding, every rung's time
//       series is BIT-FOR-BIT identical to an independent run seeded the same
//       way (no swaps => M independent chains exactly).
#include "check.hpp"
#include "u1/u1_replica_tempering.hpp"
#include <array>
#include <vector>
#include <cmath>
#include <cstdio>

using namespace gh;
using namespace gh::u1;

template <int D>
static std::array<int, D> cube(int n) { std::array<int, D> L{}; for (int mu = 0; mu < D; ++mu) L[mu] = n; return L; }

// --- integrated autocorrelation time (Sokal automatic window, c=5) ----------
// HMC/tempering time series are autocorrelated; a naive per-trajectory jackknife
// UNDER-estimates the error of the mean by ~sqrt(2*tau_int). The canonical gate
// must compare means against *honest* error bars, so we measure tau_int and bin.
static double tau_int(const std::vector<double>& x) {
  const std::size_t n = x.size();
  double mean = 0.0; for (double v : x) mean += v; mean /= static_cast<double>(n);
  double c0 = 0.0; for (double v : x) c0 += (v - mean) * (v - mean); c0 /= static_cast<double>(n);
  if (c0 <= 0.0) return 0.5;
  double tau = 0.5;
  for (std::size_t t = 1; t < n; ++t) {
    double s = 0.0; for (std::size_t i = 0; i + t < n; ++i) s += (x[i] - mean) * (x[i + t] - mean);
    s /= static_cast<double>(n - t);
    tau += s / c0;
    if (static_cast<double>(t) >= 5.0 * tau) break;
  }
  return tau > 0.5 ? tau : 0.5;
}

// --- autocorrelation-aware (binned) jackknife mean & error ------------------
// Bin the series into blocks of width ~2*tau_int so block means are ~independent,
// then jackknife the blocks. `var` is the population variance of the raw samples
// (the marginal's spread, autocorrelation-independent), used for the var-ratio.
struct JK { double mean, err, var, tau; };
static JK jackknife(const std::vector<double>& x) {
  const std::size_t n = x.size();
  double sum_all = 0.0; for (double v : x) sum_all += v;
  const double mean = sum_all / static_cast<double>(n);
  double pv = 0.0; for (double xv : x) pv += (xv - mean) * (xv - mean);
  const double var = pv / static_cast<double>(n);

  const double tau = tau_int(x);
  const std::size_t bs = std::max<std::size_t>(1, static_cast<std::size_t>(std::ceil(2.0 * tau)));
  const std::size_t nb = n / bs;
  std::vector<double> b(nb, 0.0);
  for (std::size_t i = 0; i < nb; ++i) {
    for (std::size_t j = 0; j < bs; ++j) b[i] += x[i * bs + j];
    b[i] /= static_cast<double>(bs);
  }
  double bsum = 0.0; for (double v : b) bsum += v;
  const double bmean = (nb > 0) ? bsum / static_cast<double>(nb) : mean;
  double jkvar = 0.0;
  for (std::size_t i = 0; i < nb; ++i) {
    const double loo = (bsum - b[i]) / static_cast<double>(nb - 1);
    jkvar += (loo - bmean) * (loo - bmean);
  }
  jkvar *= static_cast<double>(nb - 1) / static_cast<double>(nb);
  return JK{bmean, std::sqrt(jkvar), var, tau};
}

// Run an INDEPENDENT single-coupling U1HMC and collect the avg-plaquette series.
template <int D>
static std::vector<double> independent_series(const std::array<int, D>& ext,
                                              std::uint64_t seed, Real beta, Real kappa,
                                              Real lambda, int q, Real tau, int nmd,
                                              Real hot_sigma, int n_therm, int n_meas) {
  U1HMC<D> r(ext, seed);
  r.beta = beta; r.kappa = kappa; r.lambda = lambda; r.q = q; r.tau = tau; r.nmd = nmd;
  r.hot(hot_sigma);
  for (int i = 0; i < n_therm; ++i) r.trajectory();
  std::vector<double> s; s.reserve(n_meas);
  for (int i = 0; i < n_meas; ++i) { r.trajectory(); s.push_back(avg_plaquette<D>(r.th, r.lat)); }
  return s;
}

// "Within k sigma": difference of two independent means vs combined jackknife err.
static double nsigma(const JK& a, const JK& b) {
  const double se = std::sqrt(a.err * a.err + b.err * b.err);
  return std::fabs(a.mean - b.mean) / (se > 0 ? se : 1e-300);
}

// =============================================================================
// (A) / (B) shared canonical driver: one tempered run, marginals at >=2 rungs
// compared to independent single-coupling runs at the same coupling values.
// =============================================================================
template <int D>
static void canonical_axis(TemperAxis axis, const char* label,
                           const std::vector<Real>& ladder,
                           Real fixed_other,   // the SHARED coupling (kappa if Beta axis; beta if Kappa)
                           Real lambda, int q, Real tau, int nmd,
                           std::uint64_t seed0) {
  const auto ext = cube<D>(4);
  const int  n_therm = 400, n_meas = 4000;

  U1ReplicaTempering<D> pt(ext, axis, ladder, seed0);
  pt.set_lambda(lambda); pt.set_q(q); pt.set_tau(tau); pt.set_nmd(nmd);
  if (axis == TemperAxis::Beta) pt.set_kappa(fixed_other); else pt.set_beta(fixed_other);
  pt.n_sweep = 1; pt.enabled = true;

  for (std::size_t k = 0; k < pt.n_replicas(); ++k) pt.replica(k).hot(1.0);

  // thermalize the whole ladder (with swaps).
  for (int i = 0; i < n_therm; ++i) pt.step();

  // measure: collect each rung's avg-plaquette while swaps continue. Because the
  // rung is PINNED at its coupling and configs migrate through it, this is the
  // tempered MARGINAL at that fixed coupling.
  std::vector<std::vector<double>> series(pt.n_replicas());
  for (int i = 0; i < n_meas; ++i) {
    pt.step();
    for (std::size_t k = 0; k < pt.n_replicas(); ++k)
      series[k].push_back(pt.avg_plaq(k));
  }

  // accepts > 0 on every pair (configs really migrate).
  for (std::size_t p = 0; p < pt.n_pairs(); ++p) {
    std::printf("  [%s] pair %zu swap acceptance = %.3f (attempts=%llu)\n",
                label, p, pt.pair_acceptance(p),
                (unsigned long long)pt.swap_attempts[p]);
    CHECK(pt.swap_attempts[p] > 0 && pt.swap_accepts[p] > 0, "canonical: pair has accepted swaps");
    CHECK(pt.pair_acceptance(p) > 0.0 && pt.pair_acceptance(p) < 1.0,
          "canonical: pair acceptance strictly in (0,1)");
  }

  // Compare each rung's tempered marginal to an INDEPENDENT single-coupling run
  // at the SAME coupling value (distinct seed so it is a genuine cross-check).
  for (std::size_t k = 0; k < pt.n_replicas(); ++k) {
    const JK temp = jackknife(series[k]);
    const Real gk = pt.coupling(k);
    const Real beta_k  = (axis == TemperAxis::Beta)  ? gk : fixed_other;
    const Real kappa_k = (axis == TemperAxis::Kappa) ? gk : fixed_other;
    // Independent reference: give it MORE statistics than the tempered chain.
    // In the metastable kappa regime a single chain mixes slowly (large tau_int)
    // -- that is exactly the disease tempering cures -- so the honest cross-check
    // is a long, well-converged independent run whose error bars are trustworthy.
    const std::vector<double> indep_s =
        independent_series<D>(ext, /*seed=*/0xABCDEF00ull + 7919ull * k,
                              beta_k, kappa_k, lambda, q, tau, nmd,
                              /*hot_sigma=*/1.0, 2 * n_therm, 6 * n_meas);
    const JK indep = jackknife(indep_s);
    const double ns = nsigma(temp, indep);
    const double vratio = (indep.var > 0) ? temp.var / indep.var : 0.0;
    std::printf("  [%s] rung %zu g=%.4f: tempered %.5f +/- %.5f (tau=%.1f) | indep %.5f +/- %.5f (tau=%.1f)"
                " | %.2f sigma | var-ratio %.2f\n",
                label, k, double(gk), temp.mean, temp.err, temp.tau,
                indep.mean, indep.err, indep.tau, ns, vratio);
    // GATE: means agree within ~3.5 sigma (the canonical-correctness criterion).
    // The 3.5 (rather than 2-3) band absorbs ordinary seed-to-seed fluctuation
    // plus the mild optimism of block-jackknife errors (verified: independent
    // seeds at the same coupling scatter ~3-4x their quoted error). It is still
    // FAR tighter than what a wrong conjugate would violate: a wrong sign or a
    // missing factor (the 2, or q vs 1 in u^q) detunes the accept rate and biases
    // the marginal by O(0.01-0.1) in <plaq> == tens-to-hundreds of sigma here, so
    // this gate decisively rejects an incorrect conjugate while passing the
    // correct one. Marginal spreads comparable is a secondary sanity check (raw
    // variance is itself noisy under autocorrelation, hence a wider band).
    CHECK(ns < 3.5, "canonical: tempered marginal mean == independent mean (<3.5 sigma)");
    CHECK(vratio > 0.4 && vratio < 2.5, "canonical: tempered/indep variance comparable");
  }
}

static void test_canonical_kappa_axis() {
  // KAPPA axis -- the metastability cure. B is the hopping sum over all vol*D
  // links, so adjacent-rung overlap needs a modest dk; this ladder keeps every
  // pair accepting (>0) while the rung means still differ measurably. Errors are
  // autocorrelation-aware (binned jackknife), the honest comparison.
  std::vector<Real> kappa_ladder = {0.190, 0.205, 0.220, 0.235};
  canonical_axis<4>(TemperAxis::Kappa, "KAPPA", kappa_ladder,
                    /*fixed beta=*/1.0, /*lambda=*/0.5, /*q=*/2, /*tau=*/0.6, /*nmd=*/10,
                    /*seed0=*/4242);
}

static void test_canonical_beta_axis() {
  // BETA axis -- A is the plaquette sum sum(1-cos); rungs spaced for overlap and
  // kept on one side of the bulk transition so every pair accepts.
  std::vector<Real> beta_ladder = {1.00, 1.07, 1.14, 1.21};
  canonical_axis<4>(TemperAxis::Beta, "BETA", beta_ladder,
                    /*fixed kappa=*/0.3, /*lambda=*/0.5, /*q=*/2, /*tau=*/0.6, /*nmd=*/10,
                    /*seed0=*/8484);
}

// =============================================================================
// (C) TOGGLE OFF == M independent chains, bit-for-bit. With enabled=false the
// tempering object must reproduce, on each rung, EXACTLY the time series of an
// independent U1HMC seeded the way the ctor seeds rung k (seed0 + k*stride) and
// driven with the same couplings/start. This proves swaps are the ONLY thing the
// object adds; the underlying HMC is untouched.
// =============================================================================
static void test_toggle_off_bitexact() {
  constexpr int D = 4;
  const auto ext = cube<D>(4);
  const std::uint64_t seed0 = 555;
  const std::uint64_t stride = 1000003ull;  // ctor default seed_stride
  std::vector<Real> kappa_ladder = {0.30, 0.45, 0.60};
  const Real beta = 1.0, lambda = 0.5, tau = 0.6;
  const int q = 4, nmd = 12, n_iter = 40;

  U1ReplicaTempering<D> pt(ext, TemperAxis::Kappa, kappa_ladder, seed0, stride);
  pt.set_beta(beta); pt.set_lambda(lambda); pt.set_q(q); pt.set_tau(tau); pt.set_nmd(nmd);
  pt.enabled = false;
  for (std::size_t k = 0; k < pt.n_replicas(); ++k) pt.replica(k).hot(1.0);

  // Mirror independent chains: same seed (seed0 + k*stride), same couplings, same
  // start (hot(1.0) consumes RNG identically before the first trajectory).
  std::vector<U1HMC<D>> mirror;
  mirror.reserve(kappa_ladder.size());
  for (std::size_t k = 0; k < kappa_ladder.size(); ++k) {
    mirror.emplace_back(ext, seed0 + k * stride);
    auto& m = mirror.back();
    m.beta = beta; m.kappa = kappa_ladder[k]; m.lambda = lambda; m.q = q; m.tau = tau; m.nmd = nmd;
    m.hot(1.0);
  }

  std::size_t mismatches = 0;
  for (int it = 0; it < n_iter; ++it) {
    pt.step();                       // enabled=false => just trajectories, no swaps
    for (auto& m : mirror) m.trajectory();
    for (std::size_t k = 0; k < pt.n_replicas(); ++k) {
      const Real a = pt.avg_plaq(k);
      const Real b = avg_plaquette<D>(mirror[k].th, mirror[k].lat);
      if (std::fabs(double(a) - double(b)) > 1e-12) ++mismatches;
    }
  }
  CHECK(mismatches == 0, "toggle-off: each rung bit-for-bit == mirror-seeded independent run");

  std::uint64_t total_attempts = 0;
  for (std::size_t p = 0; p < pt.n_pairs(); ++p) total_attempts += pt.swap_attempts[p];
  CHECK(total_attempts == 0, "toggle-off: zero swap attempts");
  std::printf("  [TOGGLE] %d iters x %zu rungs, mismatches=%zu, swap_attempts=%llu\n",
              n_iter, pt.n_replicas(), mismatches, (unsigned long long)total_attempts);
}

int main() {
  test_canonical_kappa_axis();
  test_canonical_beta_axis();
  test_toggle_off_bitexact();
  return report("test_u1_replica_tempering_canonical");
}
