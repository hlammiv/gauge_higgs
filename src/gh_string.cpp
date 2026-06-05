// SU(2) -> H string-tension measurement driver, in D=4.
//
// PHYSICS. The fundamental Wilson loop W(R,T) is the static-source confinement
// probe. We accumulate the RxT loop grid over an SU(2)+Higgs HMC ensemble and
// post-process it into:
//   - the static quark-antiquark potential V(R)  (one temporal-step effective mass),
//   - the on-diagonal Creutz ratios chi(R,R)     (area-vs-perimeter discriminator),
//   - the fundamental string tension sigma_fund   (plateau of chi(R,R)),
//   - the fundamental Polyakov loop               (center-symmetry / deconfinement),
// plus the average plaquette as a cross-check. This mirrors src/screening.cpp
// (the D=3 SU(2)->2T probe-rep template) but here in D=4 with the FUNDAMENTAL loop
// and a string tension, and reuses src/hmc_higgs_multi.cpp's exact (rep, potential,
// freezing, hot-start) setup. The lattice is CUBIC L^4: the string tension needs a
// real spatial extent, so no time/space asymmetry is imposed here.
//
//   ./build/gh_string <rep> <L> <beta> <kappa> <mu2> <couplings|auto> [ntherm nmeas nmd tau seed]
//     <rep>        = fund | adj | <Young rows e.g. 6 (=SU(2) spin-3)>[:real]
//     <couplings>  = comma list f0,f1,... (one per C2 channel) or "auto" (all f_c=1).
//   Run with just <rep> to print the channel C2 values and the required #couplings.
//
//   GH_FROZEN env (any value) -> frozen-length |phi_x|=1 scalar (hmc.frozen_phi).
//
// Build (auto-discovered, default NDIM=4 NCOL=2):  make build/gh_string
#include "hmc/gauge_higgs_hmc.hpp"
#include "action/scalar_invariants.hpp"
#include "measure/observables.hpp"
#include "measure/creutz.hpp"
#include "measure/creutz_jack.hpp"
#include "rep/rep_fundamental.hpp"
#include "rep/rep_adjoint.hpp"
#include "rep/rep_general.hpp"
#include <cstdio>
#include <cstdlib>
#include <string>
#include <sstream>
#include <vector>
#include <memory>
#include <algorithm>

using namespace gh;

static_assert(kDim == 2 || kDim == 3 || kDim == 4, "gh_string is a Wilson-loop driver");
static_assert(kN == 2, "gh_string.cpp is the SU(2)->H string-tension driver: build with NCOL=2");

static std::unique_ptr<Representation<kN>> make_rep(std::string spec) {
  if (spec == "fund") return std::make_unique<FundamentalRep<kN>>();
  if (spec == "adj")  return std::make_unique<AdjointRep<kN>>();
  bool real = false;
  const auto colon = spec.find(':');
  if (colon != std::string::npos) { real = (spec.substr(colon + 1) == "real"); spec = spec.substr(0, colon); }
  std::vector<int> rows; std::stringstream ss(spec); std::string t;
  while (std::getline(ss, t, ',')) if (!t.empty()) rows.push_back(std::atoi(t.c_str()));
  if (rows.empty()) throw std::runtime_error("bad rep spec '" + spec + "'");
  return std::make_unique<GeneralRep<kN>>(rows, real ? GeneralRep<kN>::RealType::Real : GeneralRep<kN>::RealType::Complex);
}

int main(int argc, char** argv) {
  if (argc < 2) {
    std::fprintf(stderr,
      "usage: %s <rep> <L> <beta> <kappa> <mu2> <couplings|auto> [ntherm=80 nmeas=200 nmd=24 tau=1 seed=1]\n"
      "  run with just <rep> to print the channel C2 values and the required #couplings.\n", argv[0]);
    return 1;
  }
  std::unique_ptr<Representation<kN>> rep;
  try { rep = make_rep(argv[1]); }
  catch (const std::exception& e) { std::fprintf(stderr, "rep error: %s\n", e.what()); return 1; }

  CasimirChannels<kN> ch(*rep);
  std::printf("# rep=%s d=%d  %d quartic channels (C2): ", rep->name().c_str(), rep->d, ch.n_channels());
  for (Real c : ch.lambda) std::printf("%.4g ", c);
  std::printf("\n");
  if (argc < 7) { std::printf("# provide %d couplings f_c (comma list) as arg 6 (or 'auto').\n", ch.n_channels()); return 0; }

  auto argf = [&](int i, double d) { return i < argc ? std::atof(argv[i]) : d; };
  auto argi = [&](int i, long d)   { return i < argc ? std::atol(argv[i]) : d; };
  const int  Lext   = static_cast<int>(argi(2, 8));
  const Real beta   = argf(3, 2.3);
  const Real kappa  = argf(4, 0.2);
  const Real mu2    = argf(5, 1.0);
  const int  ntherm = static_cast<int>(argi(7, 80));
  const int  nmeas  = static_cast<int>(argi(8, 200));
  const int  nmd    = static_cast<int>(argi(9, 24));
  const Real tau    = argf(10, 1.0);
  const std::uint64_t seed = static_cast<std::uint64_t>(argi(11, 1));
  // Noise-guard threshold: a loop W is unreliable if W<=0 or W < n_sigma*err(W).
  // Configurable via GH_NSIGMA (default 3).
  const Real n_sigma = (std::getenv("GH_NSIGMA") ? std::atof(std::getenv("GH_NSIGMA")) : 3.0);

  std::vector<Real> f;
  if (std::string(argv[6]) == "auto") f.assign(ch.n_channels(), 1.0);
  else { std::stringstream ss(argv[6]); std::string t; while (std::getline(ss, t, ',')) if (!t.empty()) f.push_back(std::atof(t.c_str())); }
  if (static_cast<int>(f.size()) != ch.n_channels()) {
    std::fprintf(stderr, "error: need %d couplings, got %zu\n", ch.n_channels(), f.size()); return 1;
  }
  MultiInvariantPotential<kN> pot(ch, f, mu2);

  // CUBIC L^4 lattice (string tension needs a spatial extent -> no asymmetry).
  std::array<int, kDim> L{}; for (int mu = 0; mu < kDim; ++mu) L[mu] = Lext;
  GaugeHiggsHMC<kDim, kN> hmc(L, *rep, seed);
  hmc.beta = beta; hmc.kappa = kappa; hmc.tau = tau; hmc.nmd = nmd; hmc.potential = &pot;
  const bool frozen = (std::getenv("GH_FROZEN") != nullptr);   // |phi_x|=1 frozen-length scalar
  hmc.frozen_phi = frozen;

  // Wilson-loop grid extent: R,T = 1..Rmax with Rmax = min(L/2, 4).
  const int Rmax = std::min(Lext / 2, 4);
  const int Rmin = 2;   // Creutz ratios / plateau require R,T >= 2

  std::printf("# D=%d SU(%d) L=%d^%d (cubic)  beta=%.3f kappa=%.3f mu2=%.3f nmd=%d  potential=multi-invariant%s\n",
              kDim, kN, Lext, kDim, beta, kappa, mu2, nmd, frozen ? "  [FROZEN |phi|=1]" : "");
  std::printf("# Wilson-loop grid R,T = 1..%d  (Rmax = min(L/2,4)); string tension plateau over R in [%d,%d]\n",
              Rmax, Rmin, Rmax);
  std::printf("# ntherm=%d nmeas=%d tau=%.2f seed=%llu\n",
              ntherm, nmeas, tau, static_cast<unsigned long long>(seed));

  if (Rmax < 1) { std::fprintf(stderr, "error: L=%d too small for any Wilson loop\n", Lext); return 1; }

  // Start config (mirrors hmc_higgs_multi): hot by default; GH_COLD -> cold.
  if (std::getenv("GH_COLD")) { hmc.U.cold(); hmc.phi.cold(1.0); }
  else { hmc.U.hot(hmc.rng, 0.8); hmc.phi.gaussian(hmc.rng, 12345, rep->real, 0.3); }
  if (frozen) hmc.normalize_phi();   // project onto |phi_x|=1 before thermalizing
  for (int t = 0; t < ntherm; ++t) hmc.trajectory();

  hmc.traj_count = 0; hmc.accept_count = 0;
  Stats plaq, Lphi, Llink, poly;
  // Wilson-loop accumulators wl[R][T], R,T = 0..Rmax (index 0 unused, mirrors W[R][T]).
  std::vector<std::vector<Stats>> wl(Rmax + 1, std::vector<Stats>(Rmax + 1));

  // Blocked jackknife: aim for ~>=20 blocks; each block averages n_block
  // consecutive measurement trajectories into one full W[R][T] grid sample.
  const int target_blocks = 20;
  const int n_block  = std::max(1, nmeas / target_blocks);
  CreutzJackknife jk(Rmax);
  std::vector<std::vector<Real>> bsum(Rmax + 1, std::vector<Real>(Rmax + 1, 0.0));
  int bcount = 0;
  auto flush_block = [&]() {
    if (bcount == 0) return;
    WGrid g(Rmax + 1, std::vector<Real>(Rmax + 1, 0.0));
    for (int R = 1; R <= Rmax; ++R)
      for (int T = 1; T <= Rmax; ++T) g[R][T] = bsum[R][T] / bcount;
    jk.add_block(g);
    for (auto& row : bsum) std::fill(row.begin(), row.end(), 0.0);
    bcount = 0;
  };

  for (int t = 0; t < nmeas; ++t) {
    hmc.trajectory();
    plaq.add(avg_plaquette<kDim, kN>(hmc.U));
    Lphi.add(higgs_length<kDim>(hmc.phi));
    Llink.add(link_energy<kDim, kN>(hmc.phi, hmc.U, *rep));
    poly.add(polyakov_loop<kDim, kN>(hmc.U));
    for (int R = 1; R <= Rmax; ++R)
      for (int T = 1; T <= Rmax; ++T)
        wl[R][T].add(wilson_loop_fund<kDim, kN>(hmc.U, R, T));
    // accumulate this trajectory into the current block
    for (int R = 1; R <= Rmax; ++R)
      for (int T = 1; T <= Rmax; ++T)
        bsum[R][T] += wl[R][T].x.back();
    if (++bcount == n_block) flush_block();
  }
  flush_block();   // trailing partial block (kept; it is still an unbiased grid sample)

  // Mean Wilson-loop grid W[R][T] and its per-loop error grid (block-jackknife).
  const WGrid W  = jk.mean_grid();
  const WGrid We = jk.err_grid();

  // ---- phase / cross-check indicators ----
  std::printf("\n## ensemble indicators\n");
  std::printf("avg_plaquette = %.8f +/- %.8f\n", plaq.mean(), plaq.binned_error());
  std::printf("L_phi         = %.6f +/- %.6f\n", Lphi.mean(), Lphi.binned_error());
  std::printf("L_link        = %.6f +/- %.6f\n", Llink.mean(), Llink.binned_error());
  std::printf("polyakov_fund = %.6f +/- %.6f\n", poly.mean(), poly.binned_error());
  std::printf("acceptance    = %.4f\n", hmc.acceptance());

  // ---- SANITY: the 1x1 Wilson loop IS the plaquette ----
  // Compare the per-config grand mean (not the mean-of-block-means, which can differ
  // by O(1e-15) when a trailing partial block is kept) to avoid a spurious MISMATCH.
  const Real w11   = wl[1][1].mean();
  const Real resid = w11 - plaq.mean();
  std::printf("\n## sanity: W[1][1] vs avg_plaquette\n");
  std::printf("W[1][1]       = %.12f\n", w11);
  std::printf("avg_plaquette = %.12f\n", plaq.mean());
  std::printf("residual      = %.3e   (%s)\n",
              resid, (std::abs(resid) < 1e-10) ? "OK (<1e-10)" : "MISMATCH");

  std::printf("\n## blocking: %d blocks of %d traj each (last may be partial); noise threshold n_sigma=%.2f\n",
              jk.Nb(), n_block, n_sigma);

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
  return 0;
}
