#pragma once
// Jackknife + noise-guard layer for Creutz ratios / static potential.
//
// This header is ADDITIVE and is included ONLY by src/gh_string.cpp. It does NOT
// touch creutz.hpp: it reuses creutz_ratio / static_potential unchanged and adds
// (a) a blocked-sample container, (b) correlated jackknife error propagation, and
// (c) a noise guard that excludes loops that are non-positive or buried in their
// own statistical error.
//
// WHY JACKKNIFE (and not error-in-quadrature). A Creutz ratio
//   chi(R,T) = -log[ W(R,T) W(R-1,T-1) / ( W(R-1,T) W(R,T-1) ) ]
// is built from FOUR loop averages measured on the SAME configurations, so their
// fluctuations are strongly correlated. Adding the four individual loop errors in
// quadrature ignores that correlation and badly mis-estimates the error. Instead we
// keep, per measurement block b, the FULL loop grid W_b[R][T]; the leave-one-out
// mean grid Wbar^(-b) (mean of every block except b) is fed through the *same*
// estimator, giving chi_b; the jackknife mean and error
//   chi_jack = (1/Nb) sum_b chi_b ,
//   err      = sqrt[ (Nb-1)/Nb * sum_b (chi_b - chi_jack)^2 ]
// propagate the correlated four-loop error exactly.
//
// NOISE GUARD. The largest loops on a small lattice (e.g. W[3][3] ~ 1e-3) are pure
// noise: in a short run they can go negative, NaN-ing the log and producing a bogus
// negative "string tension". A loop W[R][T] is flagged UNRELIABLE if
//   W <= 0   OR   W < n_sigma * err(W).
// Any chi / V whose constituent loops include an unreliable loop is reported as
// EXCLUDED (never averaged into the plateau).
#include "measure/creutz.hpp"
#include "core/config.hpp"
#include <vector>
#include <cmath>
#include <limits>
#include <string>
#include <algorithm>

namespace gh {

using WGrid = std::vector<std::vector<Real>>;   // W[R][T], index 0 unused

// ----------------------------------------------------------------------------
// Blocked Wilson-loop samples: one full W[R][T] grid per block, plus the running
// mean grid and a per-loop error grid (std error on the mean across blocks).
// ----------------------------------------------------------------------------
struct CreutzJackknife {
  int Rmax = 0;
  std::vector<WGrid> blocks;        // blocks[b][R][T], one per measurement block

  explicit CreutzJackknife(int Rmax_) : Rmax(Rmax_) {}

  // Append a finished block (a full grid). Caller builds the block by averaging
  // n_block consecutive measurement trajectories into a (Rmax+1)x(Rmax+1) grid.
  void add_block(const WGrid& g) { blocks.push_back(g); }

  int Nb() const { return static_cast<int>(blocks.size()); }

  // Full-sample mean grid.
  WGrid mean_grid() const {
    WGrid m(Rmax + 1, std::vector<Real>(Rmax + 1, 0.0));
    const int nb = Nb();
    if (nb == 0) return m;
    for (const auto& g : blocks)
      for (int R = 0; R <= Rmax; ++R)
        for (int T = 0; T <= Rmax; ++T) m[R][T] += g[R][T];
    for (int R = 0; R <= Rmax; ++R)
      for (int T = 0; T <= Rmax; ++T) m[R][T] /= nb;
    return m;
  }

  // Std error on the mean of W[R][T] across blocks (block-jackknife-equivalent for
  // a plain mean): sqrt( Var_block / Nb ).
  WGrid err_grid() const {
    WGrid e(Rmax + 1, std::vector<Real>(Rmax + 1, 0.0));
    const int nb = Nb();
    if (nb < 2) return e;
    const WGrid m = mean_grid();
    for (int R = 0; R <= Rmax; ++R)
      for (int T = 0; T <= Rmax; ++T) {
        Real s2 = 0.0;
        for (const auto& g : blocks) { const Real d = g[R][T] - m[R][T]; s2 += d * d; }
        e[R][T] = std::sqrt(s2 / (static_cast<Real>(nb) * (nb - 1)));
      }
    return e;
  }

  // Leave-one-out mean grid Wbar^(-b): mean over all blocks except block b.
  WGrid loo_grid(int bskip) const {
    WGrid m(Rmax + 1, std::vector<Real>(Rmax + 1, 0.0));
    const int nb = Nb();
    if (nb < 2) return m;
    for (int b = 0; b < nb; ++b) {
      if (b == bskip) continue;
      for (int R = 0; R <= Rmax; ++R)
        for (int T = 0; T <= Rmax; ++T) m[R][T] += blocks[b][R][T];
    }
    for (int R = 0; R <= Rmax; ++R)
      for (int T = 0; T <= Rmax; ++T) m[R][T] /= (nb - 1);
    return m;
  }
};

// ----------------------------------------------------------------------------
// Noise guard.
// ----------------------------------------------------------------------------
// A loop is reliable iff it is strictly positive AND stands at least n_sigma above
// its own statistical error.
inline bool loop_reliable(Real w, Real err, Real n_sigma) {
  if (!(w > 0.0)) return false;            // also catches NaN
  if (!(err >= 0.0)) return false;
  if (w < n_sigma * err) return false;     // buried in noise
  return true;
}

// Reason string for an unreliable Creutz ratio's first offending constituent loop;
// empty string means all four constituents are reliable.
inline std::string creutz_excl_reason(const WGrid& Wm, const WGrid& We,
                                      int R, int T, Real n_sigma) {
  struct LP { int r, t; };
  const LP corners[4] = {{R, T}, {R - 1, T - 1}, {R - 1, T}, {R, T - 1}};
  for (const LP& c : corners) {
    const Real w = creutz_detail::at(Wm, c.r, c.t);
    const Real e = creutz_detail::at(We, c.r, c.t);
    if (!loop_reliable(w, e, n_sigma)) {
      char buf[96];
      if (!(w > 0.0))
        std::snprintf(buf, sizeof buf, "W[%d][%d]=%.2e <= 0", c.r, c.t, w);
      else
        std::snprintf(buf, sizeof buf, "W[%d][%d]=%.2e < %.0f*err(%.2e)",
                      c.r, c.t, w, n_sigma, e);
      return std::string(buf);
    }
  }
  return std::string();
}

// Reliability for a static potential point V(R,T) = -log(W[R][T]/W[R][T+1]):
// needs both W[R][T] and W[R][T+1] reliable.
inline std::string potential_excl_reason(const WGrid& Wm, const WGrid& We,
                                         int R, int T, Real n_sigma) {
  struct LP { int r, t; };
  const LP pts[2] = {{R, T}, {R, T + 1}};
  for (const LP& c : pts) {
    const Real w = creutz_detail::at(Wm, c.r, c.t);
    const Real e = creutz_detail::at(We, c.r, c.t);
    if (!loop_reliable(w, e, n_sigma)) {
      char buf[96];
      if (!(w > 0.0))
        std::snprintf(buf, sizeof buf, "W[%d][%d]=%.2e <= 0", c.r, c.t, w);
      else
        std::snprintf(buf, sizeof buf, "W[%d][%d]=%.2e < %.0f*err(%.2e)",
                      c.r, c.t, w, n_sigma, e);
      return std::string(buf);
    }
  }
  return std::string();
}

// ----------------------------------------------------------------------------
// Jackknife of a scalar estimator built from a W grid.
// ----------------------------------------------------------------------------
struct JackResult {
  Real value = std::numeric_limits<Real>::quiet_NaN();  // full-sample estimate
  Real error = std::numeric_limits<Real>::quiet_NaN();  // jackknife error
  bool ok = false;
};

// Generic jackknife: 'est' maps a W grid -> a scalar. The full-sample value is
// est(full mean grid); the error is the leave-one-out spread. If any leave-one-out
// estimate is non-finite the result is flagged not-ok (NaN error).
template <class Est>
inline JackResult jackknife(const CreutzJackknife& jk, Est&& est) {
  JackResult r;
  const int nb = jk.Nb();
  if (nb < 2) return r;
  const WGrid full = jk.mean_grid();
  r.value = est(full);
  if (!std::isfinite(r.value)) return r;

  std::vector<Real> qb;
  qb.reserve(nb);
  Real mean = 0.0;
  for (int b = 0; b < nb; ++b) {
    const Real q = est(jk.loo_grid(b));
    if (!std::isfinite(q)) return r;   // an unstable estimator -> no error
    qb.push_back(q);
    mean += q;
  }
  mean /= nb;
  Real s2 = 0.0;
  for (Real q : qb) { const Real d = q - mean; s2 += d * d; }
  r.error = std::sqrt((static_cast<Real>(nb - 1) / nb) * s2);
  r.ok = true;
  return r;
}

// On-diagonal Creutz ratio chi(R,R) with correlated jackknife error.
inline JackResult chi_diag_jack(const CreutzJackknife& jk, int R) {
  return jackknife(jk, [R](const WGrid& W) { return creutz_ratio(W, R, R); });
}

// Static potential V(R,T) with correlated jackknife error.
inline JackResult potential_jack(const CreutzJackknife& jk, int R, int T) {
  return jackknife(jk, [R, T](const WGrid& W) { return static_potential(W, R, T); });
}

// ----------------------------------------------------------------------------
// Plateau string tension over a hand-picked set of reliable R (chi(R,R)).
// Computed as the jackknife of the *average* of those chi(R,R) so that the
// correlated error across the surviving R is propagated correctly.
// ----------------------------------------------------------------------------
inline JackResult plateau_sigma_jack(const CreutzJackknife& jk,
                                     const std::vector<int>& reliableR) {
  if (reliableR.empty()) return JackResult{};
  return jackknife(jk, [&reliableR](const WGrid& W) {
    Real sum = 0.0; int n = 0;
    for (int R : reliableR) {
      const Real chi = creutz_ratio(W, R, R);
      if (!std::isfinite(chi)) return std::numeric_limits<Real>::quiet_NaN();
      sum += chi; ++n;
    }
    return n == 0 ? std::numeric_limits<Real>::quiet_NaN() : sum / n;
  });
}

// ----------------------------------------------------------------------------
// COMMON temporal extent for the V(R) linear fit.
//
// V(R) is only a clean static potential when every point is read off at the SAME
// temporal extent T: a line through points measured at INCONSISTENT T (e.g. R=1 at
// T=2 vs R=3 at T=1) is not a potential at all. So before fitting we choose ONE
// common T -- the largest temporal extent at which the static_potential estimator is
// reliable for EVERY R in the fit range -- and fit all points at that single T.
//
// static_potential(W,R,T) = -log(W[R][T]/W[R][T+1]) consumes the loops W[R][T] AND
// W[R][T+1]; a point at temporal extent T is "reliable for R" iff both of those loops
// pass the noise guard (potential_excl_reason empty). We scan T from the largest
// representable (Rmax-1, since T+1<=Rmax) downward and return the first T for which at
// least min_pts of the candidate R survive, together with the surviving R list.
// Returns T = -1 (and an empty R list) if no common T yields >= min_pts points.
// ----------------------------------------------------------------------------
struct CommonTFit {
  int T = -1;                 // chosen common temporal extent (-1 = none found)
  std::vector<int> R;         // R values reliable at that common T (the fit points)
};

inline CommonTFit common_T_potential(const WGrid& Wm, const WGrid& We,
                                     const std::vector<int>& candidateR,
                                     Real n_sigma, int Rmax, int min_pts = 2) {
  CommonTFit best;
  for (int T = Rmax - 1; T >= 1; --T) {     // need T+1 <= Rmax
    std::vector<int> surviving;
    for (int R : candidateR)
      if (potential_excl_reason(Wm, We, R, T, n_sigma).empty()) surviving.push_back(R);
    if (static_cast<int>(surviving.size()) >= min_pts) {
      best.T = T;
      best.R = std::move(surviving);
      return best;   // largest such T wins (scan is top-down)
    }
  }
  return best;       // none found
}

// ----------------------------------------------------------------------------
// Static-potential string tension at a SINGLE common temporal extent T. NOTE the
// creutz.hpp sign convention: the estimator static_potential(W,R,T) =
// -log(W[R][T]/W[R][T+1]) returns the NEGATED physical potential,
// -V_phys(R) = -(sigma R + p) for the area-law model, so its R-slope is -sigma. We fit
// a line to those estimator values across the supplied R (all at the SAME T) and
// return MINUS the slope, i.e. the physical string tension sigma_V (positive for
// confinement), matching the sign of chi(R,R). Correlated jackknife error (the whole
// fit is redone on every leave-one-out grid). Needs >= 2 R at the common T.
// ----------------------------------------------------------------------------
inline JackResult sigmaV_fit_jack(const CreutzJackknife& jk,
                                  const std::vector<int>& fitR, int commonT) {
  if (fitR.size() < 2 || commonT < 1) return JackResult{};
  return jackknife(jk, [&fitR, commonT](const WGrid& W) {
    // Ordinary least squares slope of the estimator -V_phys(R) vs R, all at T=commonT.
    Real sx = 0, sy = 0, sxx = 0, sxy = 0; int n = 0;
    for (int R : fitR) {
      const Real V = static_potential(W, R, commonT);
      if (!std::isfinite(V)) return std::numeric_limits<Real>::quiet_NaN();
      const Real x = R;
      sx += x; sy += V; sxx += x * x; sxy += x * V; ++n;
    }
    const Real denom = n * sxx - sx * sx;
    if (!(std::abs(denom) > 0.0)) return std::numeric_limits<Real>::quiet_NaN();
    const Real slope = (n * sxy - sx * sy) / denom;   // slope of -V_phys vs R = -sigma
    return -slope;                                    // physical sigma_V (>0 confining)
  });
}

}  // namespace gh
