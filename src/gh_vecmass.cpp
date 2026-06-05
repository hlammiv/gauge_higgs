// Gauge-boson (massive vector) mass measurement driver for a Higgsed SU(N) gauge
// theory. This is the runnable integration check that wires the new measurement
// module (src/measure/gaugeboson_op.hpp + gaugeboson_correlator.hpp) into an HMC
// stream and reports BOTH gauge-boson-mass estimators:
//
//   m_V(plateau) : connected, time/component-averaged zero-momentum timeslice
//                  correlator C(dt) of the electric operator O_{0i}, fed through
//                  cosh_effective_mass + plateau (delete-1 jackknife).
//   m_V(SF)      : static non-abelian magnetic structure factor + the
//                  R(p)^{-1} = a + b*phat2 pole fit (jackknife) -> m^2.
//
// MIRRORS src/hmc_higgs_multi.cpp (the SU(2)-Higgs multi-invariant HMC driver) but
// with ASYMMETRIC extents: L_t in the time direction (dir 0) and L_s in every
// spatial direction. The asymmetry (L_t > L_s) buys Euclidean-time separation for
// the timeslice correlator. Compile-time D=NDIM, N=NCOL (default 4, 2 -> SU(2)).
//
// SMOKE RUN only (a few hundred trajectories): this is an integration check that
// the module builds, runs, and yields FINITE m_V. For an SU(2) ADJOINT Higgs the
// adjoint condensate breaks SU(2)->U(1), so the surviving U(1) photon channel
// should come out LIGHT; do NOT over-interpret a short run. The deep-Higgs analytic
// ground truth (m_V^2 = (2/dim_G) C2 per SU(2) subgroup) lives in
// test/test_draft_gaugemass.cpp.
//
// usage: ./build/gh_vecmass <rep> <Ls> <Lt> <beta> <kappa> <mu2> <couplings|auto>
//                           [ntherm nmeas nmeasevery nmd tau seed]
//   <rep>       = fund | adj | <Young rows>[:real]   (default driver target: adj on SU(2))
//   <couplings> = comma list f_c (one per channel, printed-C2 order) or "auto"
//                 (all f_c=1 -> (phi^dag phi)^2, a smoke potential). Run with just
//                 <rep> to print the channel C2 values and the required #couplings.
#include "hmc/gauge_higgs_hmc.hpp"
#include "action/scalar_invariants.hpp"
#include "measure/observables.hpp"
#include "measure/correlator.hpp"
#include "measure/gaugeboson_op.hpp"
#include "measure/gaugeboson_correlator.hpp"
#include "rep/rep_fundamental.hpp"
#include "rep/rep_adjoint.hpp"
#include "rep/rep_general.hpp"
#include <cstdio>
#include <cstdlib>
#include <string>
#include <sstream>
#include <vector>
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

// Snapshot the current scalar field into the std::vector<DVec> the measurement
// module wants (one DVec per site, in lattice site order).
template <int D>
static std::vector<DVec> snapshot_phi(const CVecField<D>& phi) {
  std::vector<DVec> out;
  out.reserve(static_cast<std::size_t>(phi.lat->vol));
  for (std::int64_t s = 0; s < phi.lat->vol; ++s) out.push_back(phi.get(s));
  return out;
}

int main(int argc, char** argv) {
  std::setvbuf(stdout, nullptr, _IOLBF, 0);  // line-buffer stdout for live progress
  if (argc < 2) {
    std::fprintf(stderr,
      "usage: %s <rep> <Ls> <Lt> <beta> <kappa> <mu2> <couplings|auto>"
      " [ntherm=80 nmeas=200 nmeasevery=4 nmd=24 tau=1 seed=1]\n"
      "  run with just <rep> to print the channel C2 values and the required #couplings.\n",
      argv[0]);
    return 1;
  }
  std::unique_ptr<Representation<kN>> rep;
  try { rep = make_rep(argv[1]); }
  catch (const std::exception& e) { std::fprintf(stderr, "rep error: %s\n", e.what()); return 1; }

  CasimirChannels<kN> ch(*rep);
  std::printf("# rep=%s d=%d  %d quartic channels (C2): ", rep->name().c_str(), rep->d, ch.n_channels());
  for (Real c : ch.lambda) std::printf("%.4g ", c);
  std::printf("\n");
  if (argc < 8) { std::printf("# provide %d couplings f_c (comma list) as arg 7 (or 'auto').\n", ch.n_channels()); return 0; }

  auto argf = [&](int i, double d) { return i < argc ? std::atof(argv[i]) : d; };
  auto argi = [&](int i, long d)   { return i < argc ? std::atol(argv[i]) : d; };
  const int  Ls         = static_cast<int>(argi(2, 4));
  const int  Lt         = static_cast<int>(argi(3, 8));
  const Real beta       = argf(4, 2.3);
  const Real kappa      = argf(5, 1.5);   // deep-Higgs default
  const Real mu2        = argf(6, 1.0);
  const int  ntherm     = static_cast<int>(argi(8, 80));
  const int  nmeas      = static_cast<int>(argi(9, 200));
  const int  nmeasevery = static_cast<int>(argi(10, 4));   // trajectories between measurements
  const int  nmd        = static_cast<int>(argi(11, 24));
  const Real tau        = argf(12, 1.0);
  const std::uint64_t seed = static_cast<std::uint64_t>(argi(13, 1));

  std::vector<Real> f;
  if (std::string(argv[7]) == "auto") f.assign(ch.n_channels(), 1.0);
  else { std::stringstream ss(argv[7]); std::string t; while (std::getline(ss, t, ',')) if (!t.empty()) f.push_back(std::atof(t.c_str())); }
  if (static_cast<int>(f.size()) != ch.n_channels()) {
    std::fprintf(stderr, "error: need %d couplings, got %zu\n", ch.n_channels(), f.size()); return 1;
  }
  MultiInvariantPotential<kN> pot(ch, f, mu2);

  // ASYMMETRIC extents: dir 0 = time (Lt), dirs 1..D-1 = spatial (Ls).
  std::array<int, kDim> L{};
  L[measure::kGaugeBosonTimeDir] = Lt;
  for (int mu = 0; mu < kDim; ++mu) if (mu != measure::kGaugeBosonTimeDir) L[mu] = Ls;

  GaugeHiggsHMC<kDim, kN> hmc(L, *rep, seed);
  hmc.beta = beta; hmc.kappa = kappa; hmc.tau = tau; hmc.nmd = nmd; hmc.potential = &pot;

  std::printf("# D=%d SU(%d)  L_t=%d (time, dir 0) x L_s=%d^%d (spatial)  vol=%lld\n",
              kDim, kN, Lt, Ls, kDim - 1, static_cast<long long>(hmc.lat.vol));
  std::printf("# beta=%.3f kappa=%.3f mu2=%.3f nmd=%d tau=%.3g  ntherm=%d nmeas=%d every=%d traj  potential=multi-invariant\n",
              beta, kappa, mu2, nmd, tau, ntherm, nmeas, nmeasevery);

  // Start hot, thermalize.
  hmc.U.hot(hmc.rng, 0.8);
  hmc.phi.gaussian(hmc.rng, 12345, rep->real, 0.3);
  for (int t = 0; t < ntherm; ++t) hmc.trajectory();

  hmc.traj_count = 0; hmc.accept_count = 0;
  Stats plaq, Llink; double sExp = 0.0;
  const u1::PhotonMomenta<kDim> mom = u1::photon_momenta<kDim>(hmc.lat);

  // Per-config correlator samples (timeslice) + structure-factor samples (SF).
  std::vector<std::vector<Real>> corr_samples;   // [config][dt]
  std::vector<std::vector<Real>> sf_samples;      // [config][phat2 group]
  corr_samples.reserve(nmeas); sf_samples.reserve(nmeas);

  int nmeasured = 0;
  for (int t = 0; t < nmeas; ++t) {
    hmc.trajectory();
    sExp += std::exp(-hmc.last_dH);
    plaq.add(avg_plaquette<kDim, kN>(hmc.U));
    Llink.add(link_energy<kDim, kN>(hmc.phi, hmc.U, *rep));

    if ((t % nmeasevery) != 0) continue;   // measure every few trajectories
    const std::vector<DVec> phi_snap = snapshot_phi<kDim>(hmc.phi);

    // (a) zero-momentum timeslice field -> this config's connected correlator C(dt).
    auto S = measure::gaugeboson_timeslice_field<kDim, kN>(*rep, hmc.U, phi_snap, hmc.lat);
    corr_samples.push_back(measure::gaugeboson_correlator({S}));

    // (b) static non-abelian magnetic structure factor S(p) for the pole fit.
    sf_samples.push_back(measure::gaugeboson_structure_factor<kDim, kN>(*rep, hmc.U, phi_snap, hmc.lat, mom));
    ++nmeasured;
  }

  std::printf("\n# stream: acceptance=%.4f  <exp(-dH)>=%.5f  measured %d configs\n",
              hmc.acceptance(), sExp / nmeas, nmeasured);
  std::printf("plaquette = %.6f +/- %.6f\n", plaq.mean(), plaq.binned_error());
  std::printf("L_link    = %.6f +/- %.6f\n", Llink.mean(), Llink.binned_error());

  // ---- m_V(plateau): cosh effective mass of the connected timeslice correlator. ----
  // Plateau window: short Euclidean separation where the connected SNR is high (the
  // long-t tail dives through the cosh minimum at Lt/2 and is pure noise -- see the
  // header note in gaugeboson_correlator.hpp). Use [1, max(2, Lt/4)].
  Real mV_plat = 0.0, mV_plat_err = 0.0;
  if (!corr_samples.empty()) {
    const int Ltc  = static_cast<int>(corr_samples[0].size());
    const int tmin = 1;
    const int tmax = std::max(2, Ltc / 4);
    if (Ltc >= 3 && tmax > tmin && tmax < Ltc) {
      const PlateauFit pf = gh::plateau(corr_samples, tmin, tmax);
      mV_plat = pf.mass; mV_plat_err = pf.err;
    }
    // Mean correlator for visualization.
    std::vector<Real> Cbar(Ltc, 0.0);
    for (const auto& c : corr_samples) for (int dt = 0; dt < Ltc; ++dt) Cbar[dt] += c[dt];
    for (int dt = 0; dt < Ltc; ++dt) Cbar[dt] /= corr_samples.size();
    std::printf("\n# timeslice correlator Cbar(dt) (connected, time/component-averaged):\n");
    for (int dt = 0; dt < Ltc; ++dt) std::printf("#   C(%d) = % .6g\n", dt, Cbar[dt]);
    std::printf("==> m_V(plateau) [window t in [%d,%d]] = %.6g +/- %.4g\n",
                tmin, tmax, mV_plat, mV_plat_err);
  } else {
    std::printf("==> m_V(plateau): NO SAMPLES\n");
  }

  // ---- m_V(SF): R(p)^{-1} = a + b*phat2 pole fit of the structure factor. ----
  const int ng = mom.n_groups();
  const int nfit = (ng >= 3) ? 3 : ng;   // lowest 3 phat2 groups below the doubler edge
  const u1::PhotonMassFit fit = u1::photon_mass_fit<kDim>(sf_samples, mom, nfit);
  std::printf("\n# structure-factor R(phat2) table [R = <S(p)>/phat2 ~ Z/(phat2+m^2)]:\n");
  std::printf("#   %-12s %-14s %-14s %-12s %s\n", "phat2", "R", "R_err", "1/R", "in fit?");
  for (int g = 0; g < ng; ++g) {
    const Real invR = (g < static_cast<int>(fit.R.size()) && fit.R[g] != 0.0) ? 1.0 / fit.R[g] : 0.0;
    std::printf("#   %-12.6f %-14.6g %-14.3g %-12.6g %s\n",
                (g < static_cast<int>(fit.phat2.size()) ? fit.phat2[g] : 0.0),
                (g < static_cast<int>(fit.R.size()) ? fit.R[g] : 0.0),
                (g < static_cast<int>(fit.R_err.size()) ? fit.R_err[g] : 0.0),
                invR, (g < nfit) ? "yes" : "no");
  }
  std::printf("==> m_V(SF) [lowest %d groups] m2 = %.6g +/- %.4g   m_V = %.6g\n",
              nfit, fit.m2, fit.m2_err, fit.m_gamma);

  std::printf("\n# SUMMARY  m_V(plateau)=%.6g +/- %.4g   m_V(SF)=%.6g (m2=%.6g +/- %.4g)\n",
              mV_plat, mV_plat_err, fit.m_gamma, fit.m2, fit.m2_err);
  return 0;
}
