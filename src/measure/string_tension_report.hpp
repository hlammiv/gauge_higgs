#pragma once
// Shared string-tension REPORTING block, extracted verbatim from gh_string.cpp so
// that both the single-(beta,kappa) driver (src/gh_string.cpp) and the tempered
// beta-ladder driver (src/gh_string_pt.cpp) print the IDENTICAL sigma analysis for
// each ensemble. This header is PURE printing/post-processing: it consumes a
// finished CreutzJackknife (the blocked Wilson-loop grid) plus a handful of
// already-measured ensemble scalars and emits, in order:
//   - ## ensemble indicators        (avg_plaquette, L_phi, L_link, polyakov_fund, acceptance)
//   - ## sanity: W[1][1] vs avg_plaquette
//   - ## blocking
//   - ## fundamental Wilson loop W[R][T] grid (with noise '*' flags)
//   - ## Creutz ratios chi(R,R)
//   - ## headline tension proxy: chi(2,2)
//   - ## string tension (Creutz plateau)
//   - ## static potential V_phys(R) at a common T
//   - ## string tension (potential fit)
//   - ## string-tension estimator comparison (finite-T/L systematic)
//
// The block is byte-for-byte what gh_string.cpp used to inline (verified by diff on
// a fixed run); see report_string_tension() below.
//
// ENSEMBLE INDICATORS are passed as OPTIONAL pointers so a driver can supply only
// the scalars it actually measured. gh_string supplies all four (plaq, Lphi, Llink,
// poly); the tempered driver supplies only plaq + poly. A null pointer suppresses
// that one indicator line. To stay byte-identical with gh_string the four lines are
// printed in the same fixed order whenever their pointer is non-null.
#include "measure/creutz.hpp"
#include "measure/creutz_jack.hpp"
#include "measure/observables.hpp"   // Stats
#include "core/config.hpp"
#include <cstdio>
#include <cmath>
#include <string>
#include <vector>

namespace gh {

// All the per-ensemble scalars the report needs that are NOT derivable from `jk`.
// w11 is the per-CONFIG grand mean of the 1x1 loop (NOT jk's block mean), used by
// the sanity check exactly as gh_string did.
struct StringReportInputs {
  Real n_sigma   = 3.0;
  int  n_block   = 1;     // measurement trajectories per jackknife block (for the blocking line)
  Real w11       = 0.0;   // per-config grand mean of W[1][1] (sanity vs avg_plaquette)
  Real acceptance = 0.0;  // HMC trajectory acceptance to print
  const Stats* plaq  = nullptr;  // avg_plaquette
  const Stats* Lphi  = nullptr;  // higgs length
  const Stats* Llink = nullptr;  // link energy
  const Stats* poly  = nullptr;  // fundamental polyakov loop
};

// Reproduces gh_string.cpp's sigma-reporting block exactly (the part from the W[R][T]
// grid through the estimator-spread comparison), prefixed by the ensemble-indicator,
// sanity and blocking lines. `label` (if non-empty) is printed as a header so a
// multi-rung driver can announce which beta rung this block belongs to; gh_string
// passes an empty label so its output is unchanged.
inline void report_string_tension(const char* label,
                                   const CreutzJackknife& jk,
                                   const StringReportInputs& in) {
  const int  Rmax    = jk.Rmax;
  const int  Rmin    = 2;   // Creutz ratios / plateau require R,T >= 2
  const Real n_sigma = in.n_sigma;

  if (label && label[0]) std::printf("\n%s\n", label);

  // Mean Wilson-loop grid W[R][T] and its per-loop error grid (block-jackknife).
  const WGrid W  = jk.mean_grid();
  const WGrid We = jk.err_grid();

  // ---- phase / cross-check indicators ----
  std::printf("\n## ensemble indicators\n");
  if (in.plaq)
    std::printf("avg_plaquette = %.8f +/- %.8f\n", in.plaq->mean(), in.plaq->binned_error());
  if (in.Lphi)
    std::printf("L_phi         = %.6f +/- %.6f\n", in.Lphi->mean(), in.Lphi->binned_error());
  if (in.Llink)
    std::printf("L_link        = %.6f +/- %.6f\n", in.Llink->mean(), in.Llink->binned_error());
  if (in.poly)
    std::printf("polyakov_fund = %.6f +/- %.6f\n", in.poly->mean(), in.poly->binned_error());
  std::printf("acceptance    = %.4f\n", in.acceptance);

  // ---- SANITY: the 1x1 Wilson loop IS the plaquette ----
  // Compare the per-config grand mean (not the mean-of-block-means, which can differ
  // by O(1e-15) when a trailing partial block is kept) to avoid a spurious MISMATCH.
  const Real plaq_mean = in.plaq ? in.plaq->mean() : W[1][1];
  const Real w11   = in.w11;
  const Real resid = w11 - plaq_mean;
  std::printf("\n## sanity: W[1][1] vs avg_plaquette\n");
  std::printf("W[1][1]       = %.12f\n", w11);
  std::printf("avg_plaquette = %.12f\n", plaq_mean);
  std::printf("residual      = %.3e   (%s)\n",
              resid, (std::abs(resid) < 1e-10) ? "OK (<1e-10)" : "MISMATCH");

  std::printf("\n## blocking: %d blocks of %d traj each (last may be partial); noise threshold n_sigma=%.2f\n",
              jk.Nb(), in.n_block, n_sigma);

  // ---- Wilson-loop grid (with per-loop block-jackknife errors) ----
  // A '*' marks a loop that fails the noise guard (W<=0 or W < n_sigma*err): such a
  // loop must not enter any Creutz ratio / potential.
  std::printf("\n## fundamental Wilson loop W[R][T] = <(1/N) Re Tr U_loop>  (value +/- err; * = below noise)\n");
  std::printf("# %-4s", "R\\T");
  for (int T = 1; T <= Rmax; ++T) std::printf(" %22d", T);
  std::printf("\n");
  for (int R = 1; R <= Rmax; ++R) {
    std::printf("  %-4d", R);
    for (int T = 1; T <= Rmax; ++T) {
      const bool ok = loop_reliable(W[R][T], We[R][T], n_sigma);
      char cell[40];
      std::snprintf(cell, sizeof cell, "%.6f+/-%.6f%s", W[R][T], We[R][T], ok ? " " : "*");
      std::printf(" %22s", cell);
    }
    std::printf("\n");
  }

  // ---- Creutz ratios chi(R,R) with jackknife errors + noise exclusion ----
  // Collect the RELIABLE R for the plateau as we go.
  std::printf("\n## Creutz ratios chi(R,R)  [on-diagonal; -> sigma at large R; jackknife err]\n");
  std::printf("# %-4s %-30s %s\n", "R", "chi(R,R) +/- err", "status");
  std::vector<int> reliableR;
  for (int R = Rmin; R <= Rmax; ++R) {
    const std::string reason = creutz_excl_reason(W, We, R, R, n_sigma);
    if (!reason.empty()) {
      std::printf("  %-4d %-30s excl (%s)\n", R, "--", reason.c_str());
      continue;
    }
    const JackResult cj = chi_diag_jack(jk, R);
    if (!cj.ok) {
      std::printf("  %-4d %-30s excl (jackknife unstable)\n", R, "--");
      continue;
    }
    char val[40]; std::snprintf(val, sizeof val, "%.6f +/- %.6f", cj.value, cj.error);
    std::printf("  %-4d %-30s reliable\n", R, val);
    reliableR.push_back(R);
  }

  // ---- HEADLINE: chi(2,2) -- smallest, best-determined loop, upper bound on sigma ----
  std::printf("\n## headline tension proxy: chi(2,2)\n");
  {
    const std::string reason = creutz_excl_reason(W, We, 2, 2, n_sigma);
    const JackResult c22 = chi_diag_jack(jk, 2);
    if (reason.empty() && c22.ok)
      std::printf("chi(2,2) = %.6f +/- %.6f   (UPPER bound on asymptotic sigma_fund)\n",
                  c22.value, c22.error);
    else
      std::printf("chi(2,2) = excl (%s)   [headline unavailable]\n",
                  reason.empty() ? "jackknife unstable" : reason.c_str());
  }

  // ---- PLATEAU string tension: jackknife mean of the SURVIVING chi(R,R) ----
  // Capture the Creutz-plateau result so the final block can quote the explicit
  // spread against sigma_V as a finite-T/L systematic.
  std::printf("\n## string tension (Creutz plateau)\n");
  JackResult sigma_creutz;          // the chi-plateau estimator (ok==false if none)
  bool sigma_creutz_is_bound = false;   // true when only chi(2,2) survives (upper bound)
  if (reliableR.empty()) {
    std::printf("sigma_fund: NO reliable chi(R,R) (all loops below noise) -- run longer / smaller R\n");
  } else if (reliableR.size() == 1 && reliableR[0] == 2) {
    sigma_creutz = chi_diag_jack(jk, 2);
    sigma_creutz_is_bound = true;
    std::printf("sigma_fund: only chi(2,2) reliable = %.6f +/- %.6f (upper bound)\n",
                sigma_creutz.value, sigma_creutz.error);
  } else {
    sigma_creutz = plateau_sigma_jack(jk, reliableR);
    std::printf("sigma_fund = %.6f +/- %.6f   (jackknife plateau of chi(R,R) over reliable R = {",
                sigma_creutz.value, sigma_creutz.error);
    for (std::size_t i = 0; i < reliableR.size(); ++i)
      std::printf("%s%d", i ? "," : "", reliableR[i]);
    std::printf("})\n");
  }

  // ---- static potential V(R) at a SINGLE COMMON T, + linear-fit sigma_V ----
  // V(R) is only a clean potential when every point is read off at the SAME temporal
  // extent T. We therefore (1) collect the candidate R that have ANY reliable T, then
  // (2) pick ONE common T = the largest temporal extent reliable for EVERY R in the
  // fit range, and (3) fit V_phys(R)=V0+sigma_V*R over those R, all at that single T.
  // Physical V_phys(R) = +log(W[R][T+1]/W[R][T]) = -static_potential(W,R,T).
  std::vector<int> candR;
  for (int R = 1; R <= Rmax; ++R)
    for (int T = Rmax - 1; T >= 1; --T)
      if (potential_excl_reason(W, We, R, T, n_sigma).empty()) { candR.push_back(R); break; }

  const CommonTFit cft = common_T_potential(W, We, candR, n_sigma, Rmax, /*min_pts=*/2);

  std::printf("\n## static potential V_phys(R) = log(W[R][T]/W[R][T+1])  (COMMON T; jackknife err)\n");
  if (cft.T < 0) {
    std::printf("# no common reliable T (each R reliable only at a different T) -- V(R) is not a clean potential\n");
  } else {
    std::printf("# common T = %d  (largest temporal extent reliable for every fitted R)\n", cft.T);
    std::printf("# %-4s %-4s %-30s %s\n", "R", "T", "V_phys(R,T) +/- err", "status");
    for (int R : cft.R) {
      const JackResult vj = potential_jack(jk, R, cft.T);
      if (!vj.ok) { std::printf("  %-4d %-4d %-30s excl (jackknife unstable)\n", R, cft.T, "--"); continue; }
      // potential_jack returns the creutz.hpp estimator (= -V_phys); flip sign for display.
      char val[40]; std::snprintf(val, sizeof val, "%.6f +/- %.6f", -vj.value, vj.error);
      std::printf("  %-4d %-4d %-30s reliable\n", R, cft.T, val);
    }
  }

  // ---- string tension from the V(R) fit (common T; non-positive-slope guarded) ----
  std::printf("\n## string tension (potential fit)\n");
  JackResult sigma_V;   // ok==false unless we obtain a significantly-positive slope
  if (cft.T < 0 || cft.R.size() < 2) {
    std::printf("sigma_V: unavailable (no common reliable T with >=2 R points)\n");
  } else {
    const JackResult svj = sigmaV_fit_jack(jk, cft.R, cft.T);
    if (!svj.ok) {
      std::printf("sigma_V: linear fit unstable (jackknife)\n");
    } else if (!(svj.value > 0.0) || svj.value < svj.error) {
      // GUARD: a non-positive (or not-significantly-positive) slope is NOT a tension.
      // A negative number must never be printed as a physics string tension.
      std::printf("sigma_V: no area-law signal (slope = %.6f +/- %.6f, not significantly > 0)"
                  "  [common T=%d, R={", svj.value, svj.error, cft.T);
      for (std::size_t i = 0; i < cft.R.size(); ++i) std::printf("%s%d", i ? "," : "", cft.R[i]);
      std::printf("}]\n");
    } else {
      sigma_V = svj;   // significantly-positive slope: a genuine tension estimate
      std::printf("sigma_V = %.6f +/- %.6f   (slope of V(R)=V0+sigma_V*R at common T=%d over R = {",
                  svj.value, svj.error, cft.T);
      for (std::size_t i = 0; i < cft.R.size(); ++i) std::printf("%s%d", i ? "," : "", cft.R[i]);
      std::printf("})\n");
    }
  }

  // ---- HONEST two-estimator comparison: spread = finite-T/L SYSTEMATIC ----
  // The jackknife errors above are STATISTICAL ONLY. The chi-plateau estimator and the
  // sigma_V fit converge to the same asymptotic sigma but differ at finite T/L; that
  // difference is a real systematic, NOT covered by the statistical errors, so we never
  // claim the two "agree within errors" -- we print their spread as the systematic.
  std::printf("\n## string-tension estimator comparison (finite-T/L systematic)\n");
  if (sigma_creutz.ok)
    std::printf("sigma_Creutz (chi plateau) = %.6f +/- %.6f%s\n",
                sigma_creutz.value, sigma_creutz.error,
                sigma_creutz_is_bound ? "   [chi(2,2) upper bound only]" : "");
  else
    std::printf("sigma_Creutz (chi plateau) = unavailable\n");
  if (sigma_V.ok)
    std::printf("sigma_V      (V(R) fit)    = %.6f +/- %.6f\n", sigma_V.value, sigma_V.error);
  else
    std::printf("sigma_V      (V(R) fit)    = unavailable (no significantly-positive slope at a common T)\n");
  if (sigma_creutz.ok && sigma_V.ok) {
    const Real d = std::abs(sigma_creutz.value - sigma_V.value);
    std::printf("estimator spread |sigma_Creutz - sigma_V| = %.6f"
                "  (finite-T/L systematic at this L; NOT covered by the statistical jackknife errors)\n", d);
  } else {
    std::printf("estimator spread: n/a (one estimator unavailable)\n");
  }
}

}  // namespace gh
