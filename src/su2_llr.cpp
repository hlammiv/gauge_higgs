// SU(N)+frozen-Higgs LLR (density-of-states) driver, 1D-A scheme -- locates the precise first-order
// gauge-freezing beta_f(kappa) that HMC cannot tunnel (wide hysteresis loop). REUSES GaugeHiggsHMC
// (src/hmc/llr_hmc.hpp) so all the validated dynamics carry over; this just tiles the bistable gauge
// energy A = sum_pl(1-(1/N)ReTrU_pl) into windows, RM-solves a1(A0)=d ln rho/dA per window, and prints
// the a1(e) curve (e=A/n_plaq=1-avg_plaq). Maxwell equal-area on a1(e) -> beta_f (scripts/su2_llr_betaf.py).
//
//   ./build/su2_llr <rep> <L> <kappa> <mu2> <couplings> [e_lo e_hi ncells hwfrac K NRM tau nmd ndrive seed]
//   FROZEN |phi|=1 via env GH_FROZEN=1 (recommended); GH_GUPD center-flip honored inside each trajectory
//   only if set (usually OFF for LLR -- the window already crosses the barrier).
#include "hmc/gauge_higgs_hmc.hpp"
#include "hmc/llr_hmc.hpp"
#include "action/scalar_invariants.hpp"
#include "measure/observables.hpp"
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
  if (argc < 6) {
    std::fprintf(stderr,
      "usage: %s <rep> <L> <kappa> <mu2> <couplings|auto> "
      "[e_lo=0.08 e_hi=0.46 ncells=20 hwfrac=0.65 K=40 NRM=40 tau=0.5 nmd=12 ndrive=60 seed=1]\n"
      "  e = A/n_plaq = 1 - avg_plaquette (action density); tile e in [e_lo,e_hi] across the bistability.\n",
      argv[0]);
    return 1;
  }
  std::unique_ptr<Representation<kN>> rep;
  try { rep = make_rep(argv[1]); }
  catch (const std::exception& e) { std::fprintf(stderr, "rep error: %s\n", e.what()); return 1; }

  CasimirChannels<kN> ch(*rep);
  auto argf = [&](int i, double d) { return i < argc ? std::atof(argv[i]) : d; };
  auto argi = [&](int i, long d)   { return i < argc ? std::atol(argv[i]) : d; };
  const int  Lext  = static_cast<int>(argi(2, 8));
  const Real kappa = argf(3, 1.0);
  const Real mu2   = argf(4, 0.113);
  // couplings (arg 6)
  std::vector<Real> f;
  if (std::string(argv[5]) == "auto") f.assign(ch.n_channels(), 1.0);
  else { std::stringstream ss(argv[5]); std::string t; while (std::getline(ss, t, ',')) if (!t.empty()) f.push_back(std::atof(t.c_str())); }
  if (static_cast<int>(f.size()) != ch.n_channels()) {
    std::fprintf(stderr, "error: need %d couplings, got %zu\n", ch.n_channels(), f.size()); return 1;
  }
  const Real e_lo   = argf(6, 0.08);
  const Real e_hi   = argf(7, 0.46);
  const int  ncells = static_cast<int>(argi(8, 20));
  const Real hwfrac = argf(9, 0.65);
  const int  K      = static_cast<int>(argi(10, 40));
  const int  NRM    = static_cast<int>(argi(11, 40));
  const Real tau    = argf(12, 0.5);
  const int  nmd    = static_cast<int>(argi(13, 20));     // cold start needs >=20 at kappa~6, more at 8
  const int  nseed  = static_cast<int>(argi(14, 40));
  const Real lwfrac = argf(15, 8.0);                       // restraint stiffness s: lw = s/hw^2
  const std::uint64_t seed = static_cast<std::uint64_t>(argi(16, 1));
  const int  Kfloor = static_cast<int>(argi(19, 8));       // tau-adaptive K floor (item 1)
  const bool adaptiveK = (std::getenv("GH_LLR_FIXEDK") == nullptr);  // GH_LLR_FIXEDK=1 -> fixed K (gate ref)

  MultiInvariantPotential<kN> pot(ch, f, mu2);
  std::array<int, kDim> L{}; for (int mu = 0; mu < kDim; ++mu) L[mu] = Lext;
  GaugeHiggsHMC<kDim, kN> hmc(L, *rep, seed);
  hmc.kappa = kappa; hmc.tau = tau; hmc.nmd = nmd; hmc.potential = &pot;
  const bool frozen = (std::getenv("GH_FROZEN") != nullptr);
  hmc.frozen_phi = frozen;

  const double nplaq = static_cast<double>(hmc.lat.n_plaq());
  std::printf("# SU(%d)->rep=%s d=%d LLR(1D-A)  D=%d L=%d^%d kappa=%.3f mu2=%.3f%s\n",
              kN, rep->name().c_str(), rep->d, kDim, Lext, kDim, kappa, mu2, frozen ? " [FROZEN |phi|=1]" : "");
  const Real de = (ncells > 1) ? (e_hi - e_lo) / (ncells - 1) : 0.0;
  const Real hw_e = hwfrac * (de > 0 ? de : (e_hi - e_lo));
  const Real hw = hw_e * nplaq;                  // cell half-width in extensive A
  const Real lw = lwfrac / (hw * hw);            // Gaussian restraint stiffness (lw>lnrho'' stabilizes loop)

  std::printf("# n_plaq=%.0f  e in [%.3f,%.3f] ncells=%d hwfrac=%.2f lwfrac=%.1f K=%d NRM=%d tau=%.2f nmd=%d nseed=%d seed=%llu\n",
              nplaq, e_lo, e_hi, ncells, hwfrac, lwfrac, K, NRM, tau, nmd, nseed, (unsigned long long)seed);
  std::printf("# columns: e0  A0  a1(beta_slope)  acc_hmc  resid_e  meanE  avg_plaq  L_link  side(0=cold/ordered,1=hot/disord)\n");
  std::printf("# TWO-SIDED annealing: cold pass ascends e from the ordered end, hot pass descends from the\n");
  std::printf("# disordered end; each cell seeded from its STABLE side so the restraint only interpolates\n");
  std::printf("# the unstable middle. Filter cells by small resid_e in the post-processor.\n");

  LLRWindow<kDim, kN> llr(hmc);
  const int nmid = ncells / 2;
  const int overlap = 2;                          // both passes penetrate the middle by this many cells

  struct Row { Real e0, A0, a1, accH, resid, meanA, avgpl, llink; int side; };
  std::vector<Row> rows;
  auto emit = [&](const Row& r) {                 // SAME divided format as the sorted footer
    std::printf("%.4f %.1f %.5f %.3f %.4g %.5f %.5f %.5f %d\n",
                r.e0, r.A0, r.a1, r.accH, r.resid / nplaq, r.meanA / nplaq, r.avgpl, r.llink, r.side);
  };
  auto run_cell = [&](int i, Real& a1g, int side) {
    const Real e0 = e_lo + i * de, A0 = e0 * nplaq;
    llr.seed(a1g, A0, lw, kappa, nseed);
    double accH = 0, resid = 0, meanA = 0; int Kused = 0;
    const Real a1 = llr.rm_solve(A0, hw, lw, a1g, kappa, K, NRM, Kfloor, adaptiveK,
                                 &accH, &resid, &meanA, &Kused);   // K=Kmax; tau-adaptive in-RM (item 1)
    a1g = a1;
    Row r{e0, A0, a1, (Real)accH, (Real)resid, (Real)meanA,
          avg_plaquette<kDim, kN>(hmc.U), link_energy<kDim, kN>(hmc.phi, hmc.U, *rep), side};
    rows.push_back(r);
    emit(r); std::fflush(stdout);                 // per-cell flush: crash/timeout recovery + live monitoring
  };

  // Per-pass starting a1 must keep beta_eff on the correct side of beta_f so the seed does NOT nucleate
  // the wrong phase: cold pass starts ABOVE beta_f (stays ordered), hot pass BELOW (stays disordered).
  // 3.0/1.8 bracket beta_f for the whole BT high-kappa range (~2.1-2.6); RM then relaxes a1 to a1(e).
  const Real a1_cold = argf(17, 3.0), a1_hot = argf(18, 1.8);
  // pass A: COLD (ordered) start, ascend e from e_lo to just past the middle
  hmc.U.cold(); hmc.phi.cold(1.0); if (frozen) hmc.normalize_phi();
  { Real a1g = a1_cold; for (int i = 0; i <= std::min(ncells - 1, nmid + overlap); ++i) run_cell(i, a1g, 0); }
  // pass B: HOT (disordered) start, descend e from e_hi to just past the middle
  hmc.U.hot(hmc.rng, 0.8); hmc.phi.gaussian(hmc.rng, 12345, rep->real, 0.3); if (frozen) hmc.normalize_phi();
  { Real a1g = a1_hot; for (int i = ncells - 1; i >= std::max(0, nmid - overlap); --i) run_cell(i, a1g, 1); }

  std::printf("# total_trajectories=%ld  (adaptiveK=%d)\n", llr.traj_tot, adaptiveK ? 1 : 0);
  // Redundant SORTED footer (human-readable; post-processor de-dups vs the streamed rows by round(e,5)).
  std::printf("# --- sorted footer (rows above are the per-cell streamed copy) ---\n");
  std::sort(rows.begin(), rows.end(), [](const Row& a, const Row& b) { return a.e0 < b.e0; });
  for (const Row& r : rows) emit(r);
  return 0;
}
