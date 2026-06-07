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
#include "hmc/gauge_link_updater.hpp"
#include "core/profile.hpp"
#include "action/scalar_invariants.hpp"
#include "measure/observables.hpp"
#include "measure/creutz.hpp"
#include "measure/creutz_jack.hpp"
#include "measure/string_tension_report.hpp"
#include "rep/rep_fundamental.hpp"
#include "rep/rep_adjoint.hpp"
#include "rep/rep_general.hpp"
#include <chrono>
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
  // Test hook: GH_NO_LINK_CACHE disables the per-link D^(R)(U) cache (per-call fast_D path).
  // Used only to verify the cache is pure memoization (cache-on vs cache-off bit-identical).
  if (std::getenv("GH_NO_LINK_CACHE")) hmc.use_link_cache = false;
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
  // Tunneling-capable gauge-link updater (center flip / heatbath / metro) run between
  // HMC trajectories. Default OFF (env GH_GUPD unset) -> behavior unchanged. Lets the
  // SU(2)->2T scan cross the first-order freezing barrier that local HMC cannot tunnel.
  //   GH_GUPD="ncenter,nhb,nmetro[,width,hits]"  (e.g. GH_GUPD=4,1,2,0.5,20)
  const GaugeUpdaterCfg gupd = GaugeUpdaterCfg::from_env(std::getenv("GH_GUPD"));
  for (int t = 0; t < ntherm; ++t) {
    hmc.trajectory();
    if (gupd.enabled) run_gauge_update<kDim, kN>(hmc, gupd, hmc.traj_count);
  }

#ifdef GH_PROFILE
  // Profiling-only isolation path: run PURE MD (nmeas trajectories' worth of
  // md_evolve, no hamiltonian/measurement/accept) so the GH_PROFILE breakdown
  // is uncontaminated by observable sweeps. Enabled with GH_PROF_MDONLY=1.
  if (std::getenv("GH_PROF_MDONLY")) {
    hmc.refresh_momenta();   // populate momenta once (geometry-only; not timed-critical)
    const auto _t0 = std::chrono::steady_clock::now();
    for (int t = 0; t < nmeas; ++t) hmc.md_evolve();
    const double _wall = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - _t0).count();
    const long _steps = static_cast<long>(nmeas) * hmc.nmd;
    std::printf("\n## GH_PROF_MDONLY: %d trajectories x nmd=%d = %ld MD steps\n",
                nmeas, hmc.nmd, _steps);
    std::printf("## total MD wall = %.4f s   per-MD-step = %.4f s\n",
                _wall, _steps ? _wall / _steps : 0.0);
    GH_PROF_REPORT();
    return 0;
  }
#endif

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
    if (gupd.enabled) run_gauge_update<kDim, kN>(hmc, gupd, hmc.traj_count);
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

  // ---- shared sigma-reporting block (also used by gh_string_pt) ----
  // The whole sigma analysis (ensemble indicators -> W[R][T] grid -> Creutz ratios ->
  // chi(2,2) headline -> plateau sigma -> common-T V(R) fit -> estimator spread) lives
  // in report_string_tension(); this driver supplies the ensemble scalars it measured.
  StringReportInputs in;
  in.n_sigma    = n_sigma;
  in.n_block    = n_block;
  in.w11        = wl[1][1].mean();   // per-config grand mean (sanity vs avg_plaquette)
  in.acceptance = hmc.acceptance();
  in.plaq = &plaq; in.Lphi = &Lphi; in.Llink = &Llink; in.poly = &poly;
  report_string_tension(/*label=*/"", jk, in);
  return 0;
}
