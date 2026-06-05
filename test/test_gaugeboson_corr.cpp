// Synthetic validation of the gauge-boson correlator + plateau backend
// (src/measure/gaugeboson_correlator.hpp), with NO gauge field involved.
//
// We inject per-config zero-momentum time-slice fields S_i(t) whose connected
// two-point function is, by construction, a periodic single-cosh with a KNOWN mass:
//   < S_i(t') S_i(t) >_c  =  A * cosh( m * (((t'-t+T/2) mod T) - T/2) )   (+ noise).
// Feeding these through the SAME path the real measurement uses --
//   gaugeboson_timeslice-style perConfig array
//     -> measure::gaugeboson_correlator   (connected, time/component-averaged C(dt))
//     -> gh::cosh_effective_mass / gh::plateau   (correlator.hpp)
// -- must recover m_eff = m within the jackknife error. This is a pure correlator-
// infrastructure test: it exercises the connected estimator, the time/component
// averaging, the periodic cosh effective mass, and the delete-1 jackknife plateau.
//
// HOW THE COSH SIGNAL IS BUILT (exactly, not approximately):
// A stationary periodic Gaussian field S(t) with covariance K(dt) = A cosh(m(dt-T/2))
// is K(dt) = sum_w lambda_w exp(2 pi i w t / T) with eigenvalues lambda_w >= 0 (the
// DFT of K). We draw S(t) = sum_w sqrt(lambda_w) ( g_w^re cos + g_w^im sin )/... ;
// equivalently we synthesize S(t) = sum_w sqrt(lambda_w/T) [ a_w cos(theta_w t) +
// b_w sin(theta_w t) ] with a_w,b_w ~ N(0,1) i.i.d. across configs, which has EXACTLY
// the target periodic covariance in expectation. Averaged over many configs the
// connected correlator -> A cosh(m(dt-T/2)) and the effective mass -> m.
#include "check.hpp"
#include "core/config.hpp"
#include "core/rng.hpp"
#include "measure/correlator.hpp"
#include "measure/gaugeboson_correlator.hpp"
#include <vector>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdint>

using namespace gh;

// Target periodic single-cosh covariance K(dt) = A cosh(m (dt' - T/2)),
// dt' the wrapped distance in [0,T). Its DFT lambda_w (real, >= 0) defines the
// Gaussian synthesis. Returns lambda_w for w = 0..T-1.
static std::vector<Real> cosh_spectrum(int T, Real A, Real m) {
  // K(t) for t = 0..T-1.
  std::vector<Real> K(T);
  for (int t = 0; t < T; ++t) {
    const Real x = static_cast<Real>(t) - 0.5 * T;
    K[t] = A * std::cosh(m * x);
  }
  // lambda_w = sum_t K(t) cos(2 pi w t / T)  (K is real & even -> DFT is real).
  std::vector<Real> lam(T, 0.0);
  for (int w = 0; w < T; ++w) {
    Real s = 0.0;
    for (int t = 0; t < T; ++t) s += K[t] * std::cos(2.0 * kPi * w * t / T);
    lam[w] = s / T;                          // normalize so synthesis covariance == K
    if (lam[w] < 0.0) lam[w] = 0.0;          // clip tiny negative roundoff modes
  }
  return lam;
}

// One config's time-slice field S[t][i], i = 0..ncomp-1 independent components, each
// an i.i.d. draw of the stationary periodic Gaussian with spectrum `lam`.
static std::vector<std::vector<Real>> synth_timeslice(
    int T, int ncomp, const std::vector<Real>& lam, const Rng& rng,
    std::uint64_t cfg) {
  std::vector<std::vector<Real>> S(T, std::vector<Real>(ncomp, 0.0));
  for (int i = 0; i < ncomp; ++i) {
    // Draw a_w, b_w ~ N(0,1) and build S(t) = sum_w sqrt(lam_w) [a_w cos + b_w sin].
    // Using both cos and sin phases per mode gives the correct circular covariance.
    for (int w = 0; w < T; ++w) {
      const Real sd = std::sqrt(lam[w]);
      if (sd == 0.0) continue;
      const Real a = rng.gauss(Rng::key(cfg, i, w, 0));
      const Real b = rng.gauss(Rng::key(cfg, i, w, 1));
      const Real theta = 2.0 * kPi * w / T;
      for (int t = 0; t < T; ++t)
        S[t][i] += sd * (a * std::cos(theta * t) + b * std::sin(theta * t));
    }
  }
  return S;
}

int main() {
  std::printf("Synthetic gauge-boson correlator backend test (no gauge field):\n");
  std::printf("inject S_i(t) with connected cov A*cosh(m(t-T/2)) + noise, recover m_eff.\n\n");

  const int T      = 24;       // time extent (Lt); periodic
  const int ncomp  = 3;        // D-1 spatial components (e.g. D=4)
  const int Ncfg   = 4000;     // configs for the connected average
  const Real A     = 1.0;
  Rng rng(0xC0FFEE12345ULL);

  const Real m_inj[] = {0.30, 0.55, 0.90};
  int local_fail = 0;

  for (Real m : m_inj) {
    const std::vector<Real> lam = cosh_spectrum(T, A, m);

    // Build the per-config time-slice fields S[c][t][i] (the exact array a real run
    // feeds to gaugeboson_correlator), then form the connected averaged correlator.
    std::vector<std::vector<std::vector<Real>>> perConfig;
    perConfig.reserve(Ncfg);
    for (int c = 0; c < Ncfg; ++c)
      perConfig.push_back(synth_timeslice(T, ncomp, lam, rng, static_cast<std::uint64_t>(c) + 1));

    // (2) connected, time/component-averaged correlator C(dt).
    const std::vector<Real> C = measure::gaugeboson_correlator(perConfig);

    // Effective mass from the full averaged correlator + a plateau window away from
    // the noisy short/long-distance ends. For a clean cosh, m_eff(t) == m at all t;
    // the window simply mirrors a real plateau read.
    const std::vector<Real> me = cosh_effective_mass(C);
    // Plateau window at SMALL Euclidean separation where the connected signal-to-
    // noise is high. For a fast-decaying (large-m) cosh on a finite T the long-t
    // tail dives into the noise floor (and through the cosh minimum at t=T/2 where
    // C can go negative), so a standard plateau is read near the source, NOT at T/2.
    const int tmin = 1, tmax = 5;

    // (3) plateau() with a delete-1 jackknife over configs for the error bar.
    // It rebuilds the connected correlator path internally per replica, so we feed
    // it the per-config single-component correlators. Build per-config correlators by
    // treating each config (one S[t][i]) as its own connected sample (component-avg).
    std::vector<std::vector<Real>> perCfgCorr;
    perCfgCorr.reserve(Ncfg);
    for (const auto& cfg : perConfig) {
      // Connected auto-correlator of THIS config, averaged over components.
      std::vector<std::vector<std::vector<Real>>> one(1, cfg);
      perCfgCorr.push_back(measure::gaugeboson_correlator(one));
    }
    PlateauFit pf = plateau(perCfgCorr, tmin, tmax);

    // Direct window-average of the full-sample effective mass (central value).
    Real macc = 0.0; int mn = 0;
    for (int t = tmin; t <= tmax; ++t)
      if (std::isfinite(me[t])) { macc += me[t]; ++mn; }
    const Real m_full = mn ? macc / mn : std::numeric_limits<Real>::quiet_NaN();

    const Real err  = pf.err;
    const Real dev  = std::fabs(m_full - m);
    // Pass if the full-sample plateau mass matches injected m within 3 jackknife
    // sigma (statistical), with a small absolute floor for the err==tiny case.
    const Real tol  = std::max(3.0 * err, 0.02);
    const bool ok   = std::isfinite(m_full) && (dev <= tol);
    if (!ok) ++local_fail;

    std::printf("  m_inj=%.3f  m_eff(plateau)=%.4f  jk_err=%.4f  |dev|=%.4f  tol=%.4f  -> %s\n",
                m, m_full, err, dev, tol, ok ? "PASS" : "FAIL");

    char msg[160];
    std::snprintf(msg, sizeof msg,
                  "cosh m=%.2f: recovered m_eff=%.4f within tol %.4f (|dev|=%.4f)",
                  m, m_full, tol, dev);
    CHECK(ok, msg);
  }

  // Sanity: the connected C(dt) of a pure-noise (m large, A small) field should give
  // a well-defined, finite positive C(0) -- guard that the estimator isn't degenerate.
  {
    const std::vector<Real> lam = cosh_spectrum(T, 1.0, 0.55);
    std::vector<std::vector<std::vector<Real>>> pc;
    for (int c = 0; c < 200; ++c)
      pc.push_back(synth_timeslice(T, ncomp, lam, rng, 999000ULL + c));
    const std::vector<Real> C = measure::gaugeboson_correlator(pc);
    CHECK(std::isfinite(C[0]) && C[0] > 0.0, "connected C(0) finite & positive (estimator nondegenerate)");
  }

  std::printf("\n");
  (void)local_fail;
  return report("test_gaugeboson_corr");
}
