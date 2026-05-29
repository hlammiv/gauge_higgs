// Tests for the autocorrelation-time / error-inflation module (src/measure/autocorr.hpp).
// PURE analysis -- NO HMC, NO lattice. Every series is a deterministic synthetic
// process (seeded std::mt19937_64) whose autocorrelation structure has a KNOWN closed
// form, against which we pin the estimators.
//
// The workhorse is the AR(1) (first-order autoregressive) process
//     x_{t+1} = a x_t + sqrt(1 - a^2) eta_t,   eta ~ N(0,1) white,
// the stationary Gaussian process with EXACT normalized autocorrelation
//     rho(t) = a^|t|,
// and hence EXACT integrated autocorrelation time
//     tau_int = 1/2 + sum_{t>=1} a^t = 1/2 + a/(1-a) = (1+a) / (2(1-a)).
// (a=0   -> tau_int = 1/2,  white noise;
//  a=0.5 -> tau_int = 1.5;
//  a=0.9 -> tau_int = 9.5,  slowly mixing -> looser tolerance & longer series.)
//
// Because these are STATISTICAL recoveries (not algebraic identities) we use loose,
// honest tolerances scaled by the analytic relative statistical error of tau_int,
// sigma(tau_int)/tau_int ~ sqrt(2 (2 W + 1) / N) (Madras-Sokal), W the window.
#include "check.hpp"
#include "measure/autocorr.hpp"
#include <random>
#include <vector>
#include <cmath>
#include <cstdio>

using namespace gh;

// Analytic integrated autocorrelation time of the AR(1) process with coefficient a.
static Real ar1_tau(Real a) { return (1.0 + a) / (2.0 * (1.0 - a)); }

// Generate a stationary AR(1) series of length N. The initial value is drawn from the
// stationary distribution N(0,1) so there is no burn-in transient to discard.
static std::vector<Real> make_ar1(Real a, int N, std::uint64_t seed) {
  std::mt19937_64 rng(seed);
  std::normal_distribution<Real> g(0.0, 1.0);
  std::vector<Real> x(static_cast<std::size_t>(N));
  const Real s = std::sqrt(1.0 - a * a);
  Real prev = g(rng);  // stationary start
  x[0] = prev;
  for (int t = 1; t < N; ++t) {
    prev = a * prev + s * g(rng);
    x[static_cast<std::size_t>(t)] = prev;
  }
  return x;
}

// ---------------------------------------------------------------- 1. AR(1) recovery
// For a in {0.0, 0.5, 0.9} recover tau_int and rho(1) from a long series.
static void test_ar1() {
  struct Case { Real a; int N; std::uint64_t seed; Real tau_rtol; Real rho_atol; };
  const Case cases[] = {
    {0.0, 200000, 0xA1u, 0.06, 0.02},   // white: tau_int = 0.5
    {0.5, 200000, 0xB2u, 0.06, 0.02},   // tau_int = 1.5
    {0.9, 600000, 0xC3u, 0.18, 0.02},   // tau_int = 9.5 (slow -> longer series, looser)
  };
  for (const Case& c : cases) {
    const std::vector<Real> x = make_ar1(c.a, c.N, c.seed);

    // rho(1) ~ a.  (cap the lag: we only inspect rho(0..2), and the full O(N^2)
    // function on a 6e5-long series would be needlessly expensive.)
    const std::vector<Real> rho = autocorr_function(x, 8);
    char m[96];
    std::snprintf(m, sizeof m, "AR(1) a=%.1f: rho(0) == 1", c.a);
    CHECK_CLOSE(rho[0], 1.0, 1e-12, m);
    std::snprintf(m, sizeof m, "AR(1) a=%.1f: rho(1) ~ a", c.a);
    CHECK_CLOSE(rho[1], c.a, c.rho_atol, m);
    if (c.a > 0.0) {  // rho(2) ~ a^2 for the genuinely correlated cases
      std::snprintf(m, sizeof m, "AR(1) a=%.1f: rho(2) ~ a^2", c.a);
      CHECK_CLOSE(rho[2], c.a * c.a, c.rho_atol, m);
    }

    // tau_int recovery to relative tolerance.
    int W = 0;
    const Real tau = tau_int(x, &W);
    const Real tau_exact = ar1_tau(c.a);
    std::snprintf(m, sizeof m, "AR(1) a=%.1f: tau_int ~ (1+a)/(2(1-a))=%.3f (got %.3f, W=%d)",
                  c.a, tau_exact, tau, W);
    CHECK_CLOSE(tau, tau_exact, c.tau_rtol * tau_exact, m);

    // The chosen Sokal window must satisfy the self-consistency rule W >= c*tau_int,
    // and must be comfortably larger than tau_int (so we are not truncating signal).
    std::snprintf(m, sizeof m, "AR(1) a=%.1f: Sokal window W=%d >= c*tau", c.a, W);
    CHECK(static_cast<Real>(W) >= kSokalC * tau - 1.0, m);
  }
}

// ---------------------------------------------------------------- 2. white noise
// rho=0 (a=0) -> tau_int -> 1/2 and autocorr_error ~ naive iid error (ratio ~ 1, i.e.
// sqrt(2*tau_int) ~ 1). Also N_eff ~ N.
static void test_white_noise() {
  const int N = 200000;
  const std::vector<Real> x = make_ar1(0.0, N, 0xD4u);

  int W = 0;
  const Real tau = tau_int(x, &W);
  CHECK_CLOSE(tau, 0.5, 0.04, "white: tau_int -> 1/2");

  const Real ae = autocorr_error(x);
  const Real ne = naive_error(x);
  CHECK(ne > 0.0, "white: naive error positive");
  // sqrt(2*tau_int) ~ sqrt(2*0.5) = 1 -> autocorr_error ~ naive_error (within ~few %).
  const Real ratio = ae / ne;
  CHECK_CLOSE(ratio, 1.0, 0.05, "white: autocorr_error/naive_error ~ 1 (no inflation)");

  // N_eff ~ N for white noise (2*tau_int ~ 1).
  const Real neff = effective_sample_size(x);
  CHECK_CLOSE(neff / static_cast<Real>(N), 1.0, 0.06, "white: N_eff ~ N");
}

// ---------------------------------------------------------------- 3. error inflation
// For AR(1) the autocorrelation-aware error is the naive iid error inflated by
// sqrt(2 tau_int). Verify the ratio matches sqrt(2 tau_int_exact), and that the
// blocking error plateaus to the SAME value as the block size grows past ~2 tau_int.
static void test_error_inflation() {
  const Real a = 0.9;
  const int N = 600000;
  const std::vector<Real> x = make_ar1(a, N, 0xE5u);

  const Real tau_exact = ar1_tau(a);                 // 9.5
  const Real inflate_exact = std::sqrt(2.0 * tau_exact);

  const Real ae = autocorr_error(x);
  const Real ne = naive_error(x);
  const Real ratio = ae / ne;
  // The measured ratio is sqrt(2 * tau_int_measured); since tau_int is recovered to
  // ~15% for a=0.9, the ratio (its sqrt) is good to ~8%.
  char m[120];
  std::snprintf(m, sizeof m, "AR(1) a=%.1f: autocorr/naive ratio ~ sqrt(2 tau)=%.3f (got %.3f)",
                a, inflate_exact, ratio);
  CHECK_CLOSE(ratio, inflate_exact, 0.10 * inflate_exact, m);

  // Internal consistency: autocorr_error == sqrt(2*tau_int_measured)*naive_error EXACTLY
  // (same Var, same tau) -- algebraic identity, machine precision.
  int W = 0;
  const Real tau_meas = tau_int(x, &W);
  CHECK_CLOSE(ratio, std::sqrt(2.0 * tau_meas), 1e-10,
              "autocorr_error/naive_error == sqrt(2 tau_int) exactly");

  // Blocking plateau: small blocks underestimate (still correlated), large blocks
  // (bs >> 2 tau_int ~ 19) plateau to the autocorrelation-aware error. Track the
  // blocking error across a geometric ladder of block sizes and assert it rises
  // monotonically into a plateau that agrees with autocorr_error.
  const int bss[] = {1, 4, 16, 64, 256, 1024, 4096};
  Real prev = 0.0;
  Real plateau = 0.0;
  for (int bs : bss) {
    const Real be = blocking_error(x, bs);
    std::snprintf(m, sizeof m, "AR(1) a=%.1f: blocking error rises with block size (bs=%d)", a, bs);
    // Non-decreasing up to small statistical wiggle as decorrelation kicks in.
    CHECK(be >= prev - 0.05 * be, m);
    prev = be;
    if (bs >= 256) plateau = be;  // record a large-block value
  }
  // The plateau (large-block blocking error) agrees with the spectral autocorr_error.
  std::snprintf(m, sizeof m, "AR(1) a=%.1f: blocking plateau ~ autocorr_error (%.3e vs %.3e)",
                a, plateau, ae);
  CHECK_CLOSE(plateau, ae, 0.15 * ae, m);

  // Small-block (bs=1) blocking error is just the naive iid error (no decorrelation),
  // so it is substantially SMALLER than the plateau: the inflation is real.
  const Real be1 = blocking_error(x, 1);
  CHECK(be1 < 0.6 * plateau, "AR(1): bs=1 blocking error << plateau (inflation is real)");
  // bs=1 blocking error equals the naive iid error up to the (N-1 vs N) Bessel factor.
  CHECK_CLOSE(be1, ne, 0.02 * ne, "AR(1): bs=1 blocking error == naive iid error");
}

// ---------------------------------------------------------------- 4. exact small cases
// A handful of algebraic identities on tiny hand-checkable inputs (machine precision).
static void test_exact_small() {
  // Constant series: zero variance -> rho(0)=1, rho(t>0)=0, tau_int=1/2, errors 0.
  std::vector<Real> cst(1000, 3.7);
  const std::vector<Real> rc = autocorr_function(cst);
  CHECK_CLOSE(rc[0], 1.0, 1e-15, "const: rho(0)==1");
  CHECK_CLOSE(rc[5], 0.0, 1e-15, "const: rho(t>0)==0");
  CHECK_CLOSE(tau_int(cst), 0.5, 1e-15, "const: tau_int==1/2");
  CHECK_CLOSE(autocorr_error(cst), 0.0, 1e-15, "const: autocorr_error==0");

  // series_mean and the population variance C[0] match a brute-force hand sum.
  std::vector<Real> v{1.0, 2.0, 3.0, 4.0};  // mean 2.5, var = ( (1.5^2+0.5^2)*2 )/4 = 1.25
  CHECK_CLOSE(series_mean(v), 2.5, 1e-15, "series_mean small exact");
  const std::vector<Real> C = autocov_function(v);
  CHECK_CLOSE(C[0], 1.25, 1e-15, "autocov C(0)==population variance exact");
  // C(1) by hand: (1/4) * [ (1-2.5)(2-2.5)+(2-2.5)(3-2.5)+(3-2.5)(4-2.5) ]
  //            = (1/4) * [ 0.75 + (-0.25) + 0.75 ] = (1/4)*1.25 = 0.3125
  CHECK_CLOSE(C[1], 0.3125, 1e-15, "autocov C(1) hand value exact");
  CHECK_CLOSE(autocorr_function(v)[1], 0.3125 / 1.25, 1e-15, "rho(1)=C(1)/C(0) exact");

  // Empty / length-1 guards: no crash, sane defaults.
  std::vector<Real> empt;
  CHECK(autocorr_function(empt).empty(), "empty: rho empty");
  CHECK_CLOSE(tau_int(empt), 0.5, 1e-15, "empty: tau_int default 0.5");
  CHECK_CLOSE(autocorr_error(empt), 0.0, 1e-15, "empty: autocorr_error 0");
  std::vector<Real> one{42.0};
  CHECK_CLOSE(tau_int(one), 0.5, 1e-15, "N=1: tau_int default 0.5");
  CHECK_CLOSE(autocorr_error(one), 0.0, 1e-15, "N=1: autocorr_error 0");
}

int main() {
  std::printf("-- AR(1) tau_int & rho recovery --\n");  test_ar1();
  std::printf("-- white noise (no inflation) --\n");    test_white_noise();
  std::printf("-- error inflation & blocking plateau --\n"); test_error_inflation();
  std::printf("-- exact small / guards --\n");          test_exact_small();
  return report("test_autocorr");
}
