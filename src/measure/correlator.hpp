#pragma once
// Generic time-correlator analysis: effective masses + jackknife plateau fit.
// Lattice-agnostic -- operates purely on real correlator arrays C(t) that are
// PERIODIC in (Euclidean) time with period T = C.size(). Reused for the photon
// mass (see u1/photon_mass.hpp) and later for general spectroscopy channels.
//
//   * cosh_effective_mass : m_eff(t) = arccosh( (C[t-1]+C[t+1]) / (2 C[t]) ),
//     the right estimator for a periodic single-cosh signal C(t)~cosh(m(t-T/2)).
//   * log_effective_mass  : m(t) = log( C[t] / C[t+1] ), the open/large-T estimator.
//   * plateau             : delete-1 jackknife mean & error of the cosh effective
//     mass averaged over a [tmin,tmax] window, taken across configurations.
#include "core/config.hpp"
#include <vector>
#include <cmath>
#include <limits>

namespace gh {

// Periodic cosh effective mass. For C(t) = A cosh(m (t - T/2)) the ratio
// (C[t-1]+C[t+1])/(2 C[t]) = cosh(m) exactly, so arccosh recovers m at every t.
// The arccosh argument must be >= 1; for noisy/zero data it can dip below 1 (or
// C[t] can vanish). We CLAMP arg to be >= 1 before arccosh when arg is in (0,1),
// and emit NaN when the ratio is non-finite or <= 0 (e.g. C[t]==0, sign flips).
// Callers / plateau() skip NaN entries. Output length == C.size(); entry t uses
// the periodic neighbors (t-1+T)%T and (t+1)%T so every t is defined.
inline std::vector<Real> cosh_effective_mass(const std::vector<Real>& C) {
  const int T = static_cast<int>(C.size());
  std::vector<Real> m(T, std::numeric_limits<Real>::quiet_NaN());
  if (T < 3) return m;
  for (int t = 0; t < T; ++t) {
    const Real cm = C[(t - 1 + T) % T];
    const Real cp = C[(t + 1) % T];
    const Real c0 = C[t];
    if (c0 == 0.0) { m[t] = std::numeric_limits<Real>::quiet_NaN(); continue; }
    Real arg = (cm + cp) / (2.0 * c0);
    if (!std::isfinite(arg) || arg <= 0.0) { m[t] = std::numeric_limits<Real>::quiet_NaN(); continue; }
    if (arg < 1.0) arg = 1.0;            // clamp tiny numerical dips below 1 -> m=0
    m[t] = std::acosh(arg);
  }
  return m;
}

// Open-boundary / large-T log effective mass m(t) = log( C[t] / C[t+1] ).
// NaN where the ratio is non-positive or non-finite. Length == C.size().
inline std::vector<Real> log_effective_mass(const std::vector<Real>& C) {
  const int T = static_cast<int>(C.size());
  std::vector<Real> m(T, std::numeric_limits<Real>::quiet_NaN());
  for (int t = 0; t < T; ++t) {
    const Real c0 = C[t];
    const Real c1 = C[(t + 1) % T];
    if (c1 == 0.0) { m[t] = std::numeric_limits<Real>::quiet_NaN(); continue; }
    const Real r = c0 / c1;
    m[t] = (std::isfinite(r) && r > 0.0) ? std::log(r) : std::numeric_limits<Real>::quiet_NaN();
  }
  return m;
}

struct PlateauFit {
  Real mass = 0.0;   // jackknife mean of the windowed cosh effective mass
  Real err  = 0.0;   // delete-1 jackknife standard error of that mean
};

// Plateau average of the cosh effective mass over t in [tmin,tmax] (inclusive),
// with a delete-1 jackknife mean & error across the K configurations.
//   samples[k] is one config's correlator C_k(t) (all the same length T).
// For each jackknife replica we (a) AVERAGE the K configs' correlators with one
// config deleted, (b) form the cosh effective mass of that averaged correlator,
// (c) average it over the window (skipping NaN). The jackknife mean is the mean
// of the K replica estimates; the error is the usual delete-1 jackknife sd
//   err = sqrt( (K-1)/K * sum_k (theta_k - theta_bar)^2 ).
// Building the effective mass from the AVERAGED correlator (not averaging
// per-config effective masses) is the standard, bias-correct procedure and makes
// K identical samples give err == 0 exactly. With K==1 there is no resampling so
// err == 0 and mass is just the single-config windowed effective mass.
inline PlateauFit plateau(const std::vector<std::vector<Real>>& samples, int tmin, int tmax) {
  PlateauFit out;
  const int K = static_cast<int>(samples.size());
  if (K == 0) return out;
  const int T = static_cast<int>(samples[0].size());
  if (T == 0 || tmin > tmax || tmin < 0 || tmax >= T) return out;

  // Window-average the cosh effective mass of a single correlator (skip NaN).
  auto windowed = [&](const std::vector<Real>& C) -> Real {
    const std::vector<Real> me = cosh_effective_mass(C);
    Real s = 0.0; int n = 0;
    for (int t = tmin; t <= tmax; ++t)
      if (std::isfinite(me[t])) { s += me[t]; ++n; }
    return n > 0 ? s / n : std::numeric_limits<Real>::quiet_NaN();
  };

  if (K == 1) { out.mass = windowed(samples[0]); out.err = 0.0; return out; }

  // Total sum of correlators (for fast delete-1 averages).
  std::vector<Real> tot(T, 0.0);
  for (const auto& Ck : samples)
    for (int t = 0; t < T; ++t) tot[t] += Ck[t];

  std::vector<Real> theta(K);
  for (int k = 0; k < K; ++k) {
    std::vector<Real> Cjk(T);
    for (int t = 0; t < T; ++t) Cjk[t] = (tot[t] - samples[k][t]) / (K - 1);  // delete config k
    theta[k] = windowed(Cjk);
  }
  Real tbar = 0.0; for (Real v : theta) tbar += v; tbar /= K;
  Real sw = 0.0; for (Real v : theta) sw += (v - tbar) * (v - tbar);
  out.mass = tbar;
  out.err  = std::sqrt((static_cast<Real>(K - 1) / K) * sw);
  return out;
}

}  // namespace gh
