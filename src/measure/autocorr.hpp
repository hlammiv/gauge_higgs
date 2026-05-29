#pragma once
// Integrated autocorrelation time & autocorrelation-aware error bars for MCMC
// time series. PURE analysis: this knows nothing about the lattice, the gauge
// group, or the spacetime dimension -- feed it a 1D time series (one scalar
// observable measured once per Monte-Carlo step / trajectory) and it returns the
// normalized autocorrelation function, the integrated autocorrelation time, and
// the error on the mean inflated to honestly account for the autocorrelation.
//
// This complements Stats (observables.hpp): Stats::binned_error gives a blocking
// error; here we add the spectral (autocorrelation-function) estimate of the same
// thing plus the underlying tau_int that quantifies how many MC steps separate
// effectively independent samples -- exactly what is needed to convert a desired
// number of *independent* measurements into a trajectory budget for the phase-
// diagram campaign.
//
// DEFINITIONS (fixed; the tests pin them on analytic AR(1) / white-noise inputs)
//   * Autocovariance at lag t (population / biased 1/N normalization):
//         C(t) = (1/N) sum_{s=0}^{N-1-t} (x_s - xbar)(x_{s+t} - xbar),   t >= 0.
//     The biased 1/N (not 1/(N-t)) normalization is the standard choice in MCMC
//     error analysis (Madras-Sokal, Wolff): it tames the large-lag noise where
//     few pairs contribute and keeps the windowed sum well behaved.
//   * Normalized autocorrelation function:  rho(t) = C(t) / C(0),  rho(0) = 1.
//     For the AR(1) process x_{t+1} = a x_t + sqrt(1-a^2) eta_t one has the exact
//     rho(t) = a^|t|.
//   * Integrated autocorrelation time:
//         tau_int = 1/2 + sum_{t=1}^{W} rho(t),
//     truncated at an automatically chosen window W (see below). With the 1/2 and
//     the geometric rho(t)=a^t the closed form is tau_int = (1+a) / (2(1-a)); in
//     particular white noise (a=0) gives tau_int = 1/2.
//
// AUTOMATIC WINDOW (Madras-Sokal / Sokal self-consistent windowing)
//   Summing rho(t) to large t accumulates noise (each rho(t) has O(1/N) variance)
//   faster than signal once rho has decayed. Sokal's rule: choose the smallest
//   window M such that
//         M >= c * tau_int(M),                       (*)
//   where tau_int(M) = 1/2 + sum_{t=1}^{M} rho(t) is the running estimate and c is
//   a constant chosen so that for a single exponential autocorrelation the
//   systematic truncation bias is negligible. For rho(t) = exp(-t/tau) the bias of
//   stopping at M = c*tau is ~ exp(-c*tau/tau) = exp(-c); the conventional choice
//   c = 5..8 makes this < 1% while keeping the statistical noise small. We use
//       kSokalC = 6.0
//   (exp(-6) ~ 0.25%), the standard middle-of-the-road value (Wolff uses ~6 for
//   slowly decaying, ~4 for fast). The window is also capped at N-1.
#include "core/config.hpp"
#include <vector>
#include <cmath>
#include <cstdint>

namespace gh {

// Sokal self-consistent windowing constant c (see header comment). exp(-c) bounds
// the relative truncation bias for a single-exponential autocorrelation.
constexpr Real kSokalC = 6.0;

// Mean of a series (0 for an empty series).
inline Real series_mean(const std::vector<Real>& x) {
  if (x.empty()) return 0.0;
  Real s = 0.0;
  for (Real v : x) s += v;
  return s / static_cast<Real>(x.size());
}

// Mean square <x^2> of a series, used as the natural scale against which a variance
// is judged "genuinely zero" (no fluctuations) vs floating-point roundoff. For a
// literal constant series the population variance is mathematically zero but the
// finite-precision mean leaves O(eps) residuals -> a tiny positive C(0); we treat
// C(0) <= kDegenRelTol * <x^2> (plus a tiny absolute floor) as degenerate.
inline Real series_meansq(const std::vector<Real>& x) {
  if (x.empty()) return 0.0;
  Real s = 0.0;
  for (Real v : x) s += v * v;
  return s / static_cast<Real>(x.size());
}

// Relative tolerance below which the variance counts as floating-point noise rather
// than a real fluctuation (a constant series in finite precision lands ~1e-28 here).
constexpr Real kDegenRelTol = 1e-24;

// True if the series has no genuine fluctuations (constant up to roundoff): the
// variance C(0) is negligible compared to the data scale <x^2>.
inline bool is_degenerate(Real c0, Real meansq) {
  const Real floor_abs = 1e-300;                    // guard the all-zero series
  return c0 <= kDegenRelTol * meansq + floor_abs;
}

// Unnormalized autocovariance C(t) for lags t = 0..maxlag (clamped to N-1), with
// the biased 1/N normalization (see header). C[0] is the population variance.
// Cost is O(N * maxlag): pass a maxlag of a few hundred (>> the autocorrelation
// time) instead of the full N to keep long-series analysis cheap; the large-lag
// tail is pure noise and is discarded by the Sokal window anyway. maxlag < 0 means
// "all lags" (the full, O(N^2), function -- fine for the short series in tests).
// Returns an empty vector for an empty input; for a constant series C is all zeros.
inline std::vector<Real> autocov_function(const std::vector<Real>& x, int maxlag = -1) {
  const std::int64_t N = static_cast<std::int64_t>(x.size());
  if (N <= 0) return {};
  std::int64_t L = (maxlag < 0) ? (N - 1) : maxlag;
  if (L > N - 1) L = N - 1;
  std::vector<Real> C(static_cast<std::size_t>(L + 1), 0.0);
  const Real m = series_mean(x);
  for (std::int64_t t = 0; t <= L; ++t) {
    Real s = 0.0;
    for (std::int64_t s_idx = 0; s_idx + t < N; ++s_idx)
      s += (x[s_idx] - m) * (x[s_idx + t] - m);
    C[static_cast<std::size_t>(t)] = s / static_cast<Real>(N);  // biased 1/N
  }
  return C;
}

// Normalized autocorrelation function rho(t) = C(t)/C(0), rho(0) = 1, for lags
// 0..maxlag (clamped to N-1; maxlag<0 means all lags). If C(0) == 0 (constant
// series, no fluctuations) rho is rho(0)=1, rho(t>0)=0 (no memory) by convention.
inline std::vector<Real> autocorr_function(const std::vector<Real>& x, int maxlag = -1) {
  std::vector<Real> C = autocov_function(x, maxlag);
  if (C.empty()) return C;
  const Real c0 = C[0];
  if (is_degenerate(c0, series_meansq(x))) {  // constant / degenerate series
    std::vector<Real> rho(C.size(), 0.0);
    rho[0] = 1.0;
    return rho;
  }
  for (Real& v : C) v /= c0;
  return C;
}

// Default lag cap used by the windowed estimators (tau_int / errors) on long series:
// big enough to contain the Sokal window for any realistically slow MCMC chain, yet
// O(N * kMaxLagCap) stays cheap. If the self-consistent window is not reached within
// this cap (extremely slow mixing) the estimators sum to the cap and the returned
// window equals the cap -- a documented, conservative fallback.
constexpr int kMaxLagCap = 2000;

// Integrated autocorrelation time with the automatically chosen Sokal window.
// Optionally returns the chosen window W via *out_window. Returns 0.5 (the white-
// noise / iid value) for series too short to estimate (N < 2) or with no variance.
inline Real tau_int(const std::vector<Real>& x, int* out_window = nullptr) {
  const int N = static_cast<int>(x.size());
  if (N < 2) { if (out_window) *out_window = 0; return 0.5; }
  const int lagcap = (N - 1 < kMaxLagCap) ? (N - 1) : kMaxLagCap;
  const std::vector<Real> rho = autocorr_function(x, lagcap);
  if (rho.size() < 2 || rho[0] <= 0.0) { if (out_window) *out_window = 0; return 0.5; }

  // Running tau_int(M) = 1/2 + sum_{t=1}^{M} rho(t). Walk M up and stop at the
  // smallest window M with M >= c * tau_int(M) (Sokal). Cap M at the lag cap.
  const int Mmax = static_cast<int>(rho.size()) - 1;
  Real tau = 0.5;
  int W = Mmax;  // fall back to the full range if (*) is never satisfied
  for (int M = 1; M <= Mmax; ++M) {
    tau += rho[static_cast<std::size_t>(M)];
    if (static_cast<Real>(M) >= kSokalC * tau) { W = M; break; }
  }
  // tau now holds tau_int summed up to whatever M we stopped at. If we broke out
  // at the window that is exactly tau_int(W); if we ran to Mmax it is tau_int(Mmax).
  if (out_window) *out_window = W;
  // tau_int must be at least 1/2 (a single uncorrelated sample contributes the 1/2).
  return tau < 0.5 ? 0.5 : tau;
}

// Number of effectively independent samples N_eff = N / (2 tau_int).
inline Real effective_sample_size(const std::vector<Real>& x) {
  const int N = static_cast<int>(x.size());
  if (N <= 0) return 0.0;
  const Real tau = tau_int(x);
  return static_cast<Real>(N) / (2.0 * tau);
}

// Naive iid standard error on the mean, sqrt(Var/N) (population variance, 1/N).
// This ignores autocorrelation and underestimates the true error by sqrt(2 tau_int).
inline Real naive_error(const std::vector<Real>& x) {
  const int N = static_cast<int>(x.size());
  if (N < 2) return 0.0;
  const std::vector<Real> C = autocov_function(x, 0);  // only C[0] needed
  const Real var = C[0];  // population variance
  return std::sqrt(var / static_cast<Real>(N));
}

// Autocorrelation-aware standard error on the mean:
//     err = sqrt( 2 tau_int * Var / N ) = sqrt(2 tau_int) * naive_error.
// For white noise tau_int -> 1/2 so this reduces to the naive iid error.
inline Real autocorr_error(const std::vector<Real>& x) {
  const int N = static_cast<int>(x.size());
  if (N < 2) return 0.0;
  const std::vector<Real> C = autocov_function(x, 0);  // only C[0] needed
  const Real var = C[0];
  if (is_degenerate(var, series_meansq(x))) return 0.0;  // no fluctuations
  const Real tau = tau_int(x);
  return std::sqrt(2.0 * tau * var / static_cast<Real>(N));
}

// Blocking (batch-means) error on the mean using a given block size bs. Divides the
// series into N/bs non-overlapping blocks, averages within each block, and returns
// the standard error of those block means: sqrt( Var_blocks / nblocks ), with the
// sample (1/(nblocks-1)) variance. As bs grows past ~2 tau_int the block means
// decorrelate and this error rises to a plateau equal to autocorr_error.
inline Real blocking_error(const std::vector<Real>& x, int bs) {
  const int N = static_cast<int>(x.size());
  if (bs < 1) bs = 1;
  const int nblocks = N / bs;
  if (nblocks < 2) return 0.0;
  std::vector<Real> bm(static_cast<std::size_t>(nblocks), 0.0);
  for (int b = 0; b < nblocks; ++b) {
    Real s = 0.0;
    for (int i = b * bs; i < (b + 1) * bs; ++i) s += x[static_cast<std::size_t>(i)];
    bm[static_cast<std::size_t>(b)] = s / static_cast<Real>(bs);
  }
  Real mb = 0.0;
  for (Real v : bm) mb += v;
  mb /= static_cast<Real>(nblocks);
  Real v = 0.0;
  for (Real b : bm) v += (b - mb) * (b - mb);
  v /= static_cast<Real>(nblocks - 1);          // sample variance of block means
  return std::sqrt(v / static_cast<Real>(nblocks));  // standard error of the mean
}

}  // namespace gh
