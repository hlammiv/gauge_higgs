// Gauge + Higgs HMC with the GENERAL multi-invariant scalar potential
//   V = -mu2 (phi^dag phi) + sum_c f_c V_c[phi]
// (V_c = channel projection of phi phi^dag; tuning {f_c} aligns the VEV onto a chosen
// discrete subgroup H). Compile-time D=NDIM, N=NCOL.
//
//   ./build/hmc_higgs_multi <rep> <L> <beta> <kappa> <mu2> <couplings> [ntherm nmeas nmd tau seed]
//     <rep>        = fund | adj | <Young rows e.g. 6 (=SU(2) spin-3)>[:real]
//     <couplings>  = comma list f0,f1,... (one per channel, in the printed C2 order),
//                    or "auto" (all f_c=1, which reduces to (phi^dag phi)^2 -- a smoke test).
// Run with no couplings to first print the channel C2 values and the required count.
#include "hmc/gauge_higgs_hmc.hpp"
#include "action/scalar_invariants.hpp"
#include "measure/observables.hpp"
#include "measure/autocorr.hpp"
#include "rep/rep_fundamental.hpp"
#include "rep/rep_adjoint.hpp"
#include "rep/rep_general.hpp"
#include <cstdio>
#include <cstdlib>
#include <string>
#include <sstream>
#include <vector>
#include <fstream>
#include <memory>

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
  if (argc < 2) {
    std::fprintf(stderr,
      "usage: %s <rep> <L> <beta> <kappa> <mu2> <couplings|auto> [ntherm=80 nmeas=200 nmd=24 tau=1 seed=1 ts_outfile]\n"
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
  const std::string ts_out = (argc > 12) ? std::string(argv[12]) : std::string();  // per-traj time series

  std::vector<Real> f;
  if (std::string(argv[6]) == "auto") f.assign(ch.n_channels(), 1.0);
  else { std::stringstream ss(argv[6]); std::string t; while (std::getline(ss, t, ',')) if (!t.empty()) f.push_back(std::atof(t.c_str())); }
  if (static_cast<int>(f.size()) != ch.n_channels()) {
    std::fprintf(stderr, "error: need %d couplings, got %zu\n", ch.n_channels(), f.size()); return 1;
  }
  MultiInvariantPotential<kN> pot(ch, f, mu2);

  std::array<int, kDim> L{}; for (int mu = 0; mu < kDim; ++mu) L[mu] = Lext;
  GaugeHiggsHMC<kDim, kN> hmc(L, *rep, seed);
  hmc.beta = beta; hmc.kappa = kappa; hmc.tau = tau; hmc.nmd = nmd; hmc.potential = &pot;
  const bool frozen = (std::getenv("GH_FROZEN") != nullptr);   // |phi_x|=1 frozen-length scalar
  hmc.frozen_phi = frozen;
  std::printf("# D=%d SU(%d) L=%d^%d  beta=%.3f kappa=%.3f mu2=%.3f nmd=%d  potential=multi-invariant%s\n",
              kDim, kN, Lext, kDim, beta, kappa, mu2, nmd, frozen ? "  [FROZEN |phi|=1]" : "");

  // Start config: hot (disordered) by default; GH_COLD -> cold (ordered: identity links +
  // aligned phi). Hot-vs-cold hysteresis at fixed (beta,kappa) is the first-order test.
  if (std::getenv("GH_COLD")) { hmc.U.cold(); hmc.phi.cold(1.0); }
  else { hmc.U.hot(hmc.rng, 0.8); hmc.phi.gaussian(hmc.rng, 12345, rep->real, 0.3); }
  if (frozen) hmc.normalize_phi();   // project onto |phi_x|=1 before thermalizing
  for (int t = 0; t < ntherm; ++t) hmc.trajectory();

  hmc.traj_count = 0; hmc.accept_count = 0;
  Stats plaq, Lphi, Llink; double sExp = 0.0;
  const double V = static_cast<double>(hmc.lat.vol);
  std::vector<Real> ll_series, pl_series;            // per-trajectory series for tau_int
  ll_series.reserve(nmeas); pl_series.reserve(nmeas);
  for (int t = 0; t < nmeas; ++t) {
    hmc.trajectory();
    const Real pv = avg_plaquette<kDim, kN>(hmc.U);
    const Real lv = higgs_length<kDim>(hmc.phi);
    const Real kv = link_energy<kDim, kN>(hmc.phi, hmc.U, *rep);
    plaq.add(pv); Lphi.add(lv); Llink.add(kv);
    pl_series.push_back(pv); ll_series.push_back(kv);
    sExp += std::exp(-hmc.last_dH);
  }
  // Per-trajectory time series -> file (tau_int / blocking / thermalization checks).
  if (!ts_out.empty()) {
    std::ofstream f(ts_out);
    f << "# traj plaq L_link  (rep=" << rep->name() << " beta=" << beta
      << " kappa=" << kappa << " mu2=" << mu2 << " nmd=" << nmd << ")\n";
    for (int t = 0; t < nmeas; ++t) f << t << ' ' << pl_series[t] << ' ' << ll_series[t] << '\n';
  }
  int W_ll = 0; const Real tau_ll = tau_int(ll_series, &W_ll);
  const Real neff = effective_sample_size(ll_series);
  std::printf("plaquette   = %.6f +/- %.6f\n", plaq.mean(), plaq.binned_error());
  std::printf("L_phi       = %.6f +/- %.6f\n", Lphi.mean(), Lphi.binned_error());
  // L_link error is now the autocorrelation-aware (Madras-Sokal) estimate, not naive.
  std::printf("L_link      = %.6f +/- %.6f   chi_link=%.4f\n",
              Llink.mean(), autocorr_error(ll_series), Llink.susceptibility(V));
  std::printf("autocorr    = tau_int(L_link)=%.2f window=%d  N_eff=%.1f / %d meas\n",
              tau_ll, W_ll, neff, nmeas);
  std::printf("acceptance  = %.4f   <exp(-dH)>=%.5f\n", hmc.acceptance(), sExp / nmeas);
  return 0;
}
