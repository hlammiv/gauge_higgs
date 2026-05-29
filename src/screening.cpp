// SU(2) -> 2T (binary tetrahedral) breaking on the lattice: the Wilson-loop
// SCREENING demonstration in probe irreps, in D=3.
//
// PHYSICS. A single SU(2) Higgs in the spin-3 (dim 7) irrep with the validated
// multi-invariant 2T-locking couplings condenses onto the 2T-singlet stratum,
// Higgsing SU(2) -> 2T. Once the residual gauge group is the discrete 2T, a
// static source in probe spin j' is SCREENED iff (spin-j')|_2T contains a
// 2T-singlet (the condensate can soak up the source's color). In D=3 the
// residual 2T gauge theory CONFINES its nontrivial irreps, so the screening
// selection rule is sharp: screened reps -> the rep-traced Wilson loop saturates
// to a roughly area-independent (perimeter-law) value; H-charged reps -> the loop
// keeps falling with area. (In 4D the discrete theory is deconfined, hence D=3.)
//
//   2T branching of the low SU(2) probe spins (draft Tab. su2subduction):
//     j'=1/2 {1} d2 : half-integer, FAITHFUL  -> NO 2T-singlet  (H-charged)
//     j'=1   {2} d3 : adjoint -> T (one 3-dim 2T irrep)  -> NO singlet (Goldstones)
//     j'=3/2 {3} d4 : half-integer            -> NO 2T-singlet  (H-charged)
//     j'=2   {4} d5 : -> E1 + E2 + T          -> NO 2T-singlet  (H-charged)
//     j'=3   {6} d7 : -> A0 + 2T              -> CONTAINS A0     (SCREENED)
//   => the screening table should single out j'=3 as the screened probe.
//
// This is a SMALL/CHEAP demonstrator: D=3, 6^3 or 8^3, ~tens of trajectories,
// a short kappa scan to locate the broken (Higgs) region, then the screening
// table at the deepest kappa. Build into BUILD=build_scr with NDIM=3 NCOL=2:
//   make BUILD=build_scr NDIM=3 NCOL=2 build_scr/screening
//
//   usage: ./build_scr/screening [L beta ntherm nmeas nmd tau seed]   (defaults below)
#include "hmc/gauge_higgs_hmc.hpp"
#include "action/scalar_invariants.hpp"
#include "action/gauge_wilson.hpp"
#include "action/scalar_higgs.hpp"
#include "measure/observables.hpp"
#include "rep/rep_general.hpp"
#include <cstdio>
#include <cstdlib>
#include <array>
#include <vector>
#include <memory>
#include <string>

using namespace gh;

static_assert(kDim == 3, "screening.cpp is the D=3 confining-test driver: build with NDIM=3");
static_assert(kN == 2, "screening.cpp is the SU(2)->2T driver: build with NCOL=2");

int main(int argc, char** argv) {
  auto argf = [&](int i, double d) { return i < argc ? std::atof(argv[i]) : d; };
  auto argi = [&](int i, long d)   { return i < argc ? std::atol(argv[i]) : d; };
  const int  Lext   = static_cast<int>(argi(1, 6));     // 6^3 default (cheap)
  const Real beta   = argf(2, 2.3);
  const int  ntherm = static_cast<int>(argi(3, 40));
  const int  nmeas  = static_cast<int>(argi(4, 50));
  const int  nmd    = static_cast<int>(argi(5, 20));
  const Real tau    = argf(6, 1.0);
  const std::uint64_t seed = static_cast<std::uint64_t>(argi(7, 1));

  // ---- Higgs irrep = SU(2) spin-3 = GeneralRep<2>({6}), d=7 ----
  GeneralRep<2> R({6});
  if (R.d != 7) { std::fprintf(stderr, "ERROR: spin-3 rep d=%d (expected 7)\n", R.d); return 1; }

  // ---- multi-invariant potential with the validated 2T-locking couplings ----
  CasimirChannels<2> ch(R);
  // Validated 2T couplings (docs/locking_couplings.md), in the channel-C2 order
  // 0,2,6,12,20,30,42:
  const std::vector<Real> f2T = {0.1287, 0.1548, 0.1835, 0.2399, 0.0056, 0.1745, 0.1130};
  const Real mu2 = 0.1134;
  std::printf("# Higgs rep = %s  d=%d   %d quartic channels (C2): ",
              R.name().c_str(), R.d, ch.n_channels());
  for (Real c : ch.lambda) std::printf("%.3g ", c);
  std::printf("\n");
  if (static_cast<int>(f2T.size()) != ch.n_channels()) {
    std::fprintf(stderr, "ERROR: %d channels but %zu 2T couplings\n",
                 ch.n_channels(), f2T.size());
    return 1;
  }
  MultiInvariantPotential<2> pot(ch, f2T, mu2);

  // ---- probe reps for the screening table (DO NOT use spin>=5 / rows>=10: OOM) ----
  struct Probe { const char* spin; std::unique_ptr<GeneralRep<2>> rep; bool has_singlet; };
  std::vector<Probe> probes;
  auto add_probe = [&](const char* spin, int row, bool has_singlet) {
    probes.push_back({spin, std::make_unique<GeneralRep<2>>(std::vector<int>{row}), has_singlet});
  };
  add_probe("1/2", 1, false);  // d2  faithful, H-charged
  add_probe("1",   2, false);  // d3  adjoint -> T, no singlet
  add_probe("3/2", 3, false);  // d4  H-charged
  add_probe("2",   4, false);  // d5  E1+E2+T, no singlet
  add_probe("3",   6, true);   // d7  A0+2T, CONTAINS 2T-singlet => screened

  const std::array<std::array<int,2>,3> sizes = {{ {1,1}, {2,2}, {3,3} }};

  std::printf("# D=%d SU(%d) L=%d^%d  beta=%.3f  mu2=%.4f  nmd=%d tau=%.2f seed=%llu\n",
              kDim, kN, Lext, kDim, beta, mu2, nmd, tau,
              static_cast<unsigned long long>(seed));
  std::printf("# 2T couplings f_c =");
  for (Real v : f2T) std::printf(" %.4f", v);
  std::printf("\n# ntherm=%d nmeas=%d per kappa\n", ntherm, nmeas);

  // ---- kappa scan: find the broken/Higgs region (L_link clearly nonzero) ----
  const std::vector<Real> kappas = {0.3, 0.5, 0.8};

  // Reusable lattice extents.
  std::array<int, kDim> L{}; for (int mu = 0; mu < kDim; ++mu) L[mu] = Lext;

  struct PhaseRow { Real kappa, plaq, Lphi, Llink; };
  std::vector<PhaseRow> phase_rows;
  Real best_kappa = kappas.back();

  // We run a phase scan first, then re-run the deepest kappa to take the screening
  // table (cheap; keeps the screening measurement isolated and reproducible).
  std::printf("\n## Phase scan (locating the broken/Higgs region)\n");
  std::printf("# %-7s %-10s %-10s %-10s %-8s\n", "kappa", "plaquette", "L_phi", "L_link", "accept");
  for (Real kappa : kappas) {
    GaugeHiggsHMC<kDim, kN> hmc(L, R, seed);
    hmc.beta = beta; hmc.kappa = kappa; hmc.tau = tau; hmc.nmd = nmd; hmc.potential = &pot;
    hmc.U.hot(hmc.rng, 0.8);
    hmc.phi.gaussian(hmc.rng, 12345, R.real, 0.3);
    for (int t = 0; t < ntherm; ++t) hmc.trajectory();
    hmc.traj_count = 0; hmc.accept_count = 0;
    Stats plaq, Lphi, Llink;
    for (int t = 0; t < nmeas; ++t) {
      hmc.trajectory();
      plaq.add(avg_plaquette<kDim, kN>(hmc.U));
      Lphi.add(higgs_length<kDim>(hmc.phi));
      Llink.add(link_energy<kDim, kN>(hmc.phi, hmc.U, R));
    }
    phase_rows.push_back({kappa, plaq.mean(), Lphi.mean(), Llink.mean()});
    std::printf("  %-7.3f %-10.5f %-10.5f %-10.5f %-8.4f\n",
                kappa, plaq.mean(), Lphi.mean(), Llink.mean(), hmc.acceptance());
  }
  // Deepest kappa = strongest Higgs region for the screening table.
  best_kappa = kappas.back();

  // ---- screening table at the deepest kappa ----
  std::printf("\n## Screening table at kappa=%.3f (deep Higgs region)\n", best_kappa);
  GaugeHiggsHMC<kDim, kN> hmc(L, R, seed + 1);
  hmc.beta = beta; hmc.kappa = best_kappa; hmc.tau = tau; hmc.nmd = nmd; hmc.potential = &pot;
  hmc.U.hot(hmc.rng, 0.8);
  hmc.phi.gaussian(hmc.rng, 23456, R.real, 0.3);
  for (int t = 0; t < ntherm; ++t) hmc.trajectory();

  hmc.traj_count = 0; hmc.accept_count = 0;
  Stats plaq, Lphi, Llink;
  // wl[probe][size] accumulators
  const int nP = static_cast<int>(probes.size());
  const int nS = static_cast<int>(sizes.size());
  std::vector<std::vector<Stats>> wl(nP, std::vector<Stats>(nS));
  for (int t = 0; t < nmeas; ++t) {
    hmc.trajectory();
    plaq.add(avg_plaquette<kDim, kN>(hmc.U));
    Lphi.add(higgs_length<kDim>(hmc.phi));
    Llink.add(link_energy<kDim, kN>(hmc.phi, hmc.U, R));
    for (int p = 0; p < nP; ++p)
      for (int s = 0; s < nS; ++s)
        wl[p][s].add(wilson_loop_rep<kDim, kN>(hmc.U, *probes[p].rep,
                                               sizes[s][0], sizes[s][1]));
  }

  std::printf("# phase indicators at kappa=%.3f:\n", best_kappa);
  std::printf("#   plaquette = %.6f +/- %.6f\n", plaq.mean(),  plaq.binned_error());
  std::printf("#   L_phi     = %.6f +/- %.6f\n", Lphi.mean(),  Lphi.binned_error());
  std::printf("#   L_link    = %.6f +/- %.6f   (clearly nonzero => broken/Higgs)\n",
              Llink.mean(), Llink.binned_error());
  std::printf("#   acceptance= %.4f\n", hmc.acceptance());

  std::printf("\n# Wilson-loop-rep < (1/d) Re Tr D^(R')(W) >  vs probe spin j' and loop size\n");
  std::printf("# %-6s %-4s %-4s %-19s %-19s %-19s\n",
              "spin", "d", "2T?", "W(1x1)", "W(2x2)", "W(3x3)");
  for (int p = 0; p < nP; ++p) {
    std::printf("  %-6s %-4d %-4s", probes[p].spin, probes[p].rep->d,
                probes[p].has_singlet ? "yes" : "no");
    for (int s = 0; s < nS; ++s)
      std::printf(" %9.6f+/-%-7.6f", wl[p][s].mean(), wl[p][s].binned_error());
    std::printf("\n");
  }
  std::printf("\n# screening signal: ratio W(3x3)/W(1x1) (closer to ~const/large => screened;\n");
  std::printf("#   strongly suppressed => H-charged / area-falloff)\n");
  std::printf("# %-6s %-4s %-12s\n", "spin", "2T?", "W(3x3)/W(1x1)");
  for (int p = 0; p < nP; ++p) {
    const Real w1 = wl[p][0].mean(), w3 = wl[p][nS-1].mean();
    std::printf("  %-6s %-4s %-12.5f\n", probes[p].spin,
                probes[p].has_singlet ? "yes" : "no",
                (w1 != 0.0) ? w3 / w1 : 0.0);
  }
  std::printf("\n# 2T?=yes (spin-3, contains the 2T-singlet) is the screened probe;\n");
  std::printf("# all 2T?=no probes are H-charged and should fall off more steeply.\n");
  return 0;
}
