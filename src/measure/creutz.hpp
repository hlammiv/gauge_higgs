#pragma once
// Static potential & Creutz-ratio post-processor.
//
// Group/dimension-agnostic: this is PURE post-processing on a grid of Wilson-loop
// expectation values W(R,T). It does not know about the lattice, the gauge group,
// or the spacetime dimension -- feed it the measured loop averages and it returns
// the static quark-antiquark potential V(R) and the Creutz ratio chi(R,T), the
// classic area-law (confinement / string tension) vs perimeter-law discriminator.
//
// INPUT GRID INDEXING
//   W is std::vector<std::vector<Real>> with W[R][T] == the expectation value of the
//   R x T rectangular Wilson loop, where R is the spatial extent and T the temporal
//   extent. Both R,T >= 1 must be present. Index 0 is unused (a 1x1 loop is W[1][1]);
//   callers typically size the grid as (Rmax+1) x (Tmax+1) and leave row/column 0
//   empty so that W[R][T] reads naturally. Every accessor bounds-checks against the
//   actual W.size() / W[R].size() and returns NaN on an out-of-range request.
//
// SIGN / DEFINITION CONVENTIONS (these are fixed; tests pin them to 1e-12)
//   * static_potential(W,R,T) = -log( W[R][T] / W[R][T+1] ).
//       This is the standard one-temporal-step effective-mass estimator of V(R):
//       for W(R,T) ~ exp(-V(R) T) at large T the ratio -> exp(-V(R)) so the estimator
//       -> V(R). For a model W = exp(-sigma R T - p(R+T) - c) it evaluates exactly to
//       -(sigma R + p); hence the discrete R-slope V(R+1,T)-V(R,T) = -sigma is exactly
//       the (negated) string tension and is independent of T.
//   * creutz_ratio(W,R,T) = -log( W[R][T] W[R-1][T-1] / ( W[R-1][T] W[R][T-1] ) ),
//       requires R,T >= 2. For a pure area law the perimeter and constant pieces cancel
//       identically and chi -> sigma; for a pure perimeter law chi -> 0; for a Coulombic
//       (1/R) potential chi is positive and decreasing in R. This is the cleanest local
//       estimator of the asymptotic string tension.
#include "core/config.hpp"
#include <vector>
#include <cmath>
#include <limits>

namespace gh {

// ----------------------------------------------------------------------------- helpers
namespace creutz_detail {

constexpr Real kNaN = std::numeric_limits<Real>::quiet_NaN();

// Bounds-checked read of W[R][T]; returns NaN if (R,T) is absent from the grid.
inline Real at(const std::vector<std::vector<Real>>& W, int R, int T) {
  if (R < 0 || static_cast<std::size_t>(R) >= W.size()) return kNaN;
  if (T < 0 || static_cast<std::size_t>(T) >= W[R].size()) return kNaN;
  return W[R][T];
}

}  // namespace creutz_detail

// Static-potential estimator V(R) from a one-step temporal ratio:
//   V(R) = -log( W[R][T] / W[R][T+1] ).
// Needs W[R][T] and W[R][T+1] present and strictly positive. Returns NaN otherwise.
inline Real static_potential(const std::vector<std::vector<Real>>& W, int R, int T) {
  const Real w0 = creutz_detail::at(W, R, T);
  const Real w1 = creutz_detail::at(W, R, T + 1);
  if (!(w0 > 0.0) || !(w1 > 0.0)) return creutz_detail::kNaN;  // also rejects NaN
  return -std::log(w0 / w1);
}

// Creutz ratio chi(R,T) -- the area-vs-perimeter discriminator:
//   chi(R,T) = -log( W[R][T] * W[R-1][T-1] / ( W[R-1][T] * W[R][T-1] ) ).
// Requires R,T >= 2 (it samples the 2x2 block of corners (R,T),(R-1,T-1),(R-1,T),(R,T-1)).
// Any missing / non-positive entry yields NaN. For a pure area law chi == sigma exactly
// (perimeter + constant cancel); for a pure perimeter law chi == 0 exactly.
inline Real creutz_ratio(const std::vector<std::vector<Real>>& W, int R, int T) {
  if (R < 2 || T < 2) return creutz_detail::kNaN;  // documented lower bound
  const Real a = creutz_detail::at(W, R,     T);
  const Real b = creutz_detail::at(W, R - 1, T - 1);
  const Real c = creutz_detail::at(W, R - 1, T);
  const Real d = creutz_detail::at(W, R,     T - 1);
  if (!(a > 0.0) || !(b > 0.0) || !(c > 0.0) || !(d > 0.0)) return creutz_detail::kNaN;
  return -std::log((a * b) / (c * d));
}

// Plateau (large-R,T) string tension: average of the on-diagonal Creutz ratios
// chi(R,R) for R in [Rmin, Rmax] (the most symmetric, fastest-converging estimators).
// Skips any (R,R) that returns NaN (missing/non-positive). Returns NaN if none usable.
// Rmax<0 means "use the largest R the grid supports".
inline Real string_tension_plateau(const std::vector<std::vector<Real>>& W,
                                    int Rmin = 2, int Rmax = -1) {
  if (Rmin < 2) Rmin = 2;
  if (Rmax < 0) Rmax = static_cast<int>(W.size()) - 1;  // largest representable R
  Real sum = 0.0;
  int n = 0;
  for (int R = Rmin; R <= Rmax; ++R) {
    const Real chi = creutz_ratio(W, R, R);
    if (std::isnan(chi)) continue;
    sum += chi;
    ++n;
  }
  return n == 0 ? creutz_detail::kNaN : sum / n;
}

// Perimeter-coefficient estimate p of W ~ exp(-sigma R T - p (R+T) - c). With the
// string tension sigma supplied (e.g. from string_tension_plateau), the static
// potential satisfies V(R,T) = -(sigma R + p), so p = -V(R,T) - sigma R, exactly and
// T-independent for the pure model. We average over the available (R,T) with R,T+1
// present. Returns NaN if no usable point exists.
inline Real perimeter_coeff(const std::vector<std::vector<Real>>& W, Real sigma,
                            int Rmin = 1, int Rmax = -1) {
  if (Rmin < 1) Rmin = 1;
  if (Rmax < 0) Rmax = static_cast<int>(W.size()) - 1;
  Real sum = 0.0;
  int n = 0;
  for (int R = Rmin; R <= Rmax; ++R) {
    // Use the smallest available T for this R that yields a valid V(R,T).
    for (std::size_t T = 1; T + 1 < W[R].size(); ++T) {
      const Real V = static_potential(W, R, static_cast<int>(T));
      if (std::isnan(V)) continue;
      sum += (-V - sigma * R);  // p = -V(R,T) - sigma R
      ++n;
      break;  // one (R,T) per R is enough; the model is T-independent
    }
  }
  return n == 0 ? creutz_detail::kNaN : sum / n;
}

}  // namespace gh
