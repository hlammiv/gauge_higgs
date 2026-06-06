// SU(2) -> H string-tension TEMPERED scan driver, in D=4.
//
// PHYSICS / METHOD. This is the parallel-tempered sibling of src/gh_string.cpp.
// Instead of one (beta,kappa) ensemble we run a whole beta-ladder
//   beta_0 < beta_1 < ... < beta_{M-1}     (kappa, mu2, couplings shared)
// of replicas at once via ReplicaTempering<D,N> (src/hmc/replica_tempering.hpp).
// After every replica advances n_sweep trajectories, adjacent rungs attempt config
// swaps (even/odd passes), so a configuration can random-walk up and down the ladder
// and decorrelate far faster than M independent runs -- especially across a
// confinement/Coulomb (or Higgs) transition where a single beta tunnels slowly.
//
// For EACH beta rung we accumulate the fundamental Wilson-loop grid W[R][T]
// (R,T=1..Rmax, Rmax=min(L/2,4)) into a per-rung blocked CreutzJackknife and the
// avg-plaquette + fundamental Polyakov loop, then hand each rung to the SHARED
// report_string_tension() (src/measure/string_tension_report.hpp) so every rung
// prints the IDENTICAL sigma analysis gh_string prints for a single ensemble. We
// then print the swap-acceptance table and whether the ladder is CONNECTED (every
// adjacent pair accepts with rate > 0; a zero-rate pair is a broken ladder rung).
//
//   ./build/gh_string_pt <rep> <L> <beta_list> <kappa> <mu2> <couplings|auto>
//                        [ntherm nmeas n_sweep nmd tau seed]
//     <rep>        = fund | adj | <Young rows e.g. 6 (=SU(2) spin-3)>[:real]
//     <beta_list>  = comma list  1.0,1.4,1.7,2.0,2.5   OR  min:max:step  1.0:2.5:0.3
//                    (SORTED ascending; must be STRICTLY increasing -- PT requires it)
//     <couplings>  = comma list f0,f1,... (one per C2 channel) or "auto" (all f_c=1).
//   Run with just <rep> to print the channel C2 values and the required #couplings.
//
//   GH_FROZEN env (any value) -> frozen-length |phi_x|=1 scalar on every replica.
//   GH_TEMPER env: tempering defaults ON; set 0/false/off/no to DISABLE swaps
//                  (-> M independent single-beta runs, the merge-safe no-op).
//   GH_NSIGMA env: loop noise-guard threshold (default 3).
//   GH_COLD env:   cold start (default hot).
//
// Build (auto-discovered, default NDIM=4 NCOL=2):  make build/gh_string_pt
#include "hmc/gauge_higgs_hmc.hpp"
#include "hmc/replica_tempering.hpp"
#include "action/scalar_invariants.hpp"
#include "measure/observables.hpp"
#include "measure/creutz.hpp"
#include "measure/creutz_jack.hpp"
#include "measure/string_tension_report.hpp"
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

static_assert(kDim == 2 || kDim == 3 || kDim == 4, "gh_string_pt is a Wilson-loop driver");
static_assert(kN == 2, "gh_string_pt.cpp is the SU(2)->H string-tension driver: build with NCOL=2");

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

// Parse a beta-ladder string: "1.0,1.4,1.7" (comma list) OR "min:max:step"
// (inclusive of max within fp tolerance). Result is SORTED ascending; duplicates are
// removed and the strict-increase requirement is enforced by ReplicaTempering's ctor.
static std::vector<Real> parse_beta_ladder(const std::string& s) {
  std::vector<Real> betas;
  if (s.find(':') != std::string::npos) {
    // min:max:step
    std::vector<Real> parts; std::stringstream ss(s); std::string t;
    while (std::getline(ss, t, ':')) if (!t.empty()) parts.push_back(std::atof(t.c_str()));
    if (parts.size() != 3) throw std::runtime_error("beta range needs min:max:step");
    const Real lo = parts[0], hi = parts[1], step = parts[2];
    if (!(step > 0.0)) throw std::runtime_error("beta range step must be > 0");
    for (Real b = lo; b <= hi + 1e-9 * std::max<Real>(1.0, std::abs(hi)); b += step)
      betas.push_back(b);
  } else {
    std::stringstream ss(s); std::string t;
    while (std::getline(ss, t, ',')) if (!t.empty()) betas.push_back(std::atof(t.c_str()));
  }
  std::sort(betas.begin(), betas.end());
  // Drop exact/near duplicates so the strict-increase ctor check does not trip on
  // a user-supplied repeated value; genuinely-too-close values still throw there.
  std::vector<Real> uniq;
  for (Real b : betas)
    if (uniq.empty() || b > uniq.back() + 1e-12) uniq.push_back(b);
  return uniq;
}

int main(int argc, char** argv) {
  if (argc < 2) {
    std::fprintf(stderr,
      "usage: %s <rep> <L> <beta_list> <kappa> <mu2> <couplings|auto> "
      "[ntherm=80 nmeas=200 n_sweep=1 nmd=24 tau=1 seed=1]\n"
      "  <beta_list> = comma list (1.0,1.4,1.7,2.0,2.5) OR min:max:step (1.0:2.5:0.3)\n"
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
  const int  Lext    = static_cast<int>(argi(2, 8));
  // arg 3 is the beta LIST (string); kappa/mu2 shift to args 4/5.
  std::vector<Real> betas;
  try { betas = parse_beta_ladder(argv[3]); }
  catch (const std::exception& e) { std::fprintf(stderr, "beta-list error: %s\n", e.what()); return 1; }
  if (betas.size() < 2) {
    std::fprintf(stderr, "error: tempering needs >=2 betas, got %zu (give a comma list or min:max:step)\n",
                 betas.size());
    return 1;
  }
  const Real kappa   = argf(4, 0.2);
  const Real mu2     = argf(5, 1.0);
  const int  ntherm  = static_cast<int>(argi(7, 80));
  const int  nmeas   = static_cast<int>(argi(8, 200));
  const int  n_sweep = static_cast<int>(argi(9, 1));
  const int  nmd     = static_cast<int>(argi(10, 24));
  const Real tau     = argf(11, 1.0);
  const std::uint64_t seed = static_cast<std::uint64_t>(argi(12, 1));
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

  // Build the tempered ladder. The ctor enforces strictly-increasing betas (we sorted
  // + de-duped above) and seeds replica k with seed + k*seed_stride.
  std::unique_ptr<ReplicaTempering<kDim, kN>> ptp;
  try { ptp = std::make_unique<ReplicaTempering<kDim, kN>>(L, *rep, betas, seed); }
  catch (const std::exception& e) { std::fprintf(stderr, "ladder error: %s\n", e.what()); return 1; }
  ReplicaTempering<kDim, kN>& pt = *ptp;

  // Shared couplings on every replica (kappa/potential/tau/nmd identical; beta differs).
  pt.set_kappa(kappa);
  pt.set_potential(&pot);
  pt.set_tau(tau);
  pt.set_nmd(nmd);
  pt.n_sweep = n_sweep;
  const bool frozen = (std::getenv("GH_FROZEN") != nullptr);
  if (frozen) pt.set_frozen_phi(true);

  // TOGGLE: tempering defaults ON; GH_TEMPER=0/off/false/no disables swaps. We do NOT
  // gate on the env unconditionally (env_enabled() is false when unset) -- the default
  // here is ON, and only an explicit falsy GH_TEMPER turns it off.
  pt.enabled = true;
  if (const char* e = std::getenv("GH_TEMPER")) {
    std::string v(e);
    std::transform(v.begin(), v.end(), v.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    if (v == "0" || v == "false" || v == "off" || v == "no") pt.enabled = false;
  }

  const int M = static_cast<int>(pt.n_replicas());
  const int Rmax = std::min(Lext / 2, 4);
  const int Rmin = 2;

  std::printf("# D=%d SU(%d) L=%d^%d (cubic)  kappa=%.3f mu2=%.3f nmd=%d  potential=multi-invariant%s\n",
              kDim, kN, Lext, kDim, kappa, mu2, nmd, frozen ? "  [FROZEN |phi|=1]" : "");
  std::printf("# TEMPERED beta-ladder (%d rungs): ", M);
  for (int k = 0; k < M; ++k) std::printf("%s%.4g", k ? "," : "", pt.beta(k));
  std::printf("\n");
  std::printf("# tempering swaps = %s   (GH_TEMPER toggle; OFF -> %d independent runs)\n",
              pt.enabled ? "ON" : "OFF", M);
  std::printf("# Wilson-loop grid R,T = 1..%d  (Rmax = min(L/2,4)); string tension plateau over R in [%d,%d]\n",
              Rmax, Rmin, Rmax);
  std::printf("# ntherm=%d nmeas=%d n_sweep=%d tau=%.2f seed=%llu\n",
              ntherm, nmeas, n_sweep, tau, static_cast<unsigned long long>(seed));

  if (Rmax < 1) { std::fprintf(stderr, "error: L=%d too small for any Wilson loop\n", Lext); return 1; }

  // -------- INITIALIZE each replica's fields (ctor does NOT start them) --------
  // Mirror gh_string's start: hot U + Gaussian phi (or cold via GH_COLD); normalize
  // phi onto |phi|=1 if frozen. Each replica uses its OWN rng so the streams differ.
  const bool cold = (std::getenv("GH_COLD") != nullptr);
  for (int k = 0; k < M; ++k) {
    GaugeHiggsHMC<kDim, kN>& r = pt.replica(k);
    if (cold) { r.U.cold(); r.phi.cold(1.0); }
    else      { r.U.hot(r.rng, 0.8); r.phi.gaussian(r.rng, 12345, rep->real, 0.3); }
    if (frozen) r.normalize_phi();   // project onto |phi_x|=1 before thermalizing
  }

  // -------- THERMALIZE: ntherm tempered steps (each = n_sweep traj/replica + swaps) --
  for (int t = 0; t < ntherm; ++t) pt.step();

  // Reset per-replica HMC accept counters so the reported acceptance covers only the
  // measurement phase (mirrors gh_string's traj_count/accept_count reset).
  for (int k = 0; k < M; ++k) { pt.replica(k).traj_count = 0; pt.replica(k).accept_count = 0; }

  // -------- MEASURE: nmeas tempered steps, accumulating per-rung observables --------
  // Per rung: avg-plaquette + fundamental Polyakov Stats, and a blocked CreutzJackknife
  // of the W[R][T] grid (~20 blocks like gh_string). We also keep wl[k][R][T] Stats so
  // the W[1][1] sanity uses the per-config grand mean exactly as gh_string does.
  const int target_blocks = 20;
  const int n_block = std::max(1, nmeas / target_blocks);

  std::vector<Stats> plaq(M), poly(M);
  std::vector<std::vector<std::vector<Stats>>> wl(
      M, std::vector<std::vector<Stats>>(Rmax + 1, std::vector<Stats>(Rmax + 1)));
  std::vector<CreutzJackknife> jks(M, CreutzJackknife(Rmax));
  // Per-rung running block buffer + block count.
  std::vector<std::vector<std::vector<Real>>> bsum(
      M, std::vector<std::vector<Real>>(Rmax + 1, std::vector<Real>(Rmax + 1, 0.0)));
  std::vector<int> bcount(M, 0);

  auto flush_block = [&](int k) {
    if (bcount[k] == 0) return;
    WGrid g(Rmax + 1, std::vector<Real>(Rmax + 1, 0.0));
    for (int R = 1; R <= Rmax; ++R)
      for (int T = 1; T <= Rmax; ++T) g[R][T] = bsum[k][R][T] / bcount[k];
    jks[k].add_block(g);
    for (auto& row : bsum[k]) std::fill(row.begin(), row.end(), 0.0);
    bcount[k] = 0;
  };

  for (int t = 0; t < nmeas; ++t) {
    pt.step();
    for (int k = 0; k < M; ++k) {
      const GaugeHiggsHMC<kDim, kN>& r = pt.replica(k);
      plaq[k].add(avg_plaquette<kDim, kN>(r.U));
      poly[k].add(polyakov_loop<kDim, kN>(r.U));
      for (int R = 1; R <= Rmax; ++R)
        for (int T = 1; T <= Rmax; ++T)
          wl[k][R][T].add(wilson_loop_fund<kDim, kN>(r.U, R, T));
      for (int R = 1; R <= Rmax; ++R)
        for (int T = 1; T <= Rmax; ++T)
          bsum[k][R][T] += wl[k][R][T].x.back();
      if (++bcount[k] == n_block) flush_block(k);
    }
  }
  for (int k = 0; k < M; ++k) flush_block(k);   // trailing partial block per rung

  // -------- PER-RUNG sigma report (identical analysis to gh_string) --------
  for (int k = 0; k < M; ++k) {
    char label[128];
    std::snprintf(label, sizeof label,
                  "########## rung %d/%d : beta = %.4g  (hmc_acc=%.4f) ##########",
                  k, M, pt.beta(k), pt.replica(k).acceptance());
    StringReportInputs in;
    in.n_sigma    = n_sigma;
    in.n_block    = n_block;
    in.w11        = wl[k][1][1].mean();
    in.acceptance = pt.replica(k).acceptance();
    in.plaq = &plaq[k];
    in.poly = &poly[k];
    report_string_tension(label, jks[k], in);
  }

  // -------- SWAP ACCEPTANCE table + ladder connectivity --------
  std::printf("\n########## tempering swap acceptance ##########\n");
  if (!pt.enabled) {
    std::printf("# swaps DISABLED (GH_TEMPER off): %d independent single-beta runs, no exchanges.\n", M);
  } else {
    std::printf("# adjacent-pair (k,k+1) config-swap acceptance over the measurement phase:\n");
    std::printf("# %-6s %-14s %-14s %s\n", "pair", "beta_k", "beta_{k+1}", "accept_rate");
    bool connected = true;
    for (std::size_t p = 0; p < pt.n_pairs(); ++p) {
      const double rate = pt.pair_acceptance(p);
      if (!(rate > 0.0)) connected = false;
      std::printf("  (%zu,%zu)  %-14.4g %-14.4g %.4f%s\n",
                  p, p + 1, pt.beta(p), pt.beta(p + 1), rate,
                  (rate > 0.0) ? "" : "   <-- BROKEN (no accepted swaps)");
    }
    std::printf("ladder connectivity: %s  (%s)\n",
                connected ? "CONNECTED" : "BROKEN",
                connected ? "every adjacent pair accepted >0 swaps"
                          : "at least one adjacent pair never swapped -- configs cannot traverse the ladder");
  }
  return 0;
}
