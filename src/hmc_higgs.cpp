// Gauge + Higgs HMC driver (arbitrary irrep). Compile-time D=NDIM, N=NCOL.
//   ./build/hmc_higgs <rep> <L> <beta> <kappa> <lambda> [ntherm] [nmeas] [nmd] [tau] [seed]
//     <rep> = fund | adj
// Prints plaquette, Higgs length <phi^dag phi>, gauge-invariant link energy
// L_link (the transition locator), acceptance, and <exp(-dH)>.
#include "hmc/gauge_higgs_hmc.hpp"
#include "measure/observables.hpp"
#include "rep/rep_fundamental.hpp"
#include "rep/rep_adjoint.hpp"
#include "rep/rep_general.hpp"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <string>
#include <vector>
#include <sstream>

using namespace gh;

// rep spec: "fund" | "adj" | "<r1,r2,...>" Young rows (general engine), optional
// ":real" suffix to declare a real-scalar rep, e.g. "2,1:real" (SU(3) adjoint).
static std::unique_ptr<Representation<kN>> make_rep(std::string spec) {
  if (spec == "fund") return std::make_unique<FundamentalRep<kN>>();
  if (spec == "adj")  return std::make_unique<AdjointRep<kN>>();
  bool real = false;
  const auto colon = spec.find(':');
  if (colon != std::string::npos) { real = (spec.substr(colon + 1) == "real"); spec = spec.substr(0, colon); }
  std::vector<int> rows; std::stringstream ss(spec); std::string tok;
  while (std::getline(ss, tok, ',')) if (!tok.empty()) rows.push_back(std::atoi(tok.c_str()));
  if (rows.empty()) throw std::runtime_error("bad rep spec '" + spec + "'");
  return std::make_unique<GeneralRep<kN>>(
      rows, real ? GeneralRep<kN>::RealType::Real : GeneralRep<kN>::RealType::Complex);
}

int main(int argc, char** argv) {
  if (argc < 6) {
    std::fprintf(stderr,
      "usage: %s <rep> <L> <beta> <kappa> <lambda> [ntherm=80] [nmeas=200] [nmd=24] [tau=1.0] [seed=1]\n"
      "  <rep> = fund | adj | <Young rows e.g. 2,1> [add ':real' for a real-scalar rep, e.g. 2,1:real]\n",
      argv[0]);
    return 1;
  }
  const char* repname = argv[1];
  auto argf = [&](int i, double def) { return i < argc ? std::atof(argv[i]) : def; };
  auto argi = [&](int i, long def)   { return i < argc ? std::atol(argv[i]) : def; };
  const int  Lext   = static_cast<int>(argi(2, 8));
  const Real beta   = argf(3, 2.3);
  const Real kappa  = argf(4, 0.2);
  const Real lambda = argf(5, 0.5);
  const int  ntherm = static_cast<int>(argi(6, 80));
  const int  nmeas  = static_cast<int>(argi(7, 200));
  const int  nmd    = static_cast<int>(argi(8, 24));
  const Real tau    = argf(9, 1.0);
  const std::uint64_t seed = static_cast<std::uint64_t>(argi(10, 1));

  std::unique_ptr<Representation<kN>> rep;
  try { rep = make_rep(repname); }
  catch (const std::exception& e) { std::fprintf(stderr, "rep error: %s\n", e.what()); return 1; }

  std::array<int, kDim> L{};
  for (int mu = 0; mu < kDim; ++mu) L[mu] = Lext;

  GaugeHiggsHMC<kDim, kN> hmc(L, *rep, seed);
  hmc.beta = beta; hmc.kappa = kappa; hmc.lambda = lambda; hmc.tau = tau; hmc.nmd = nmd;

  std::printf("# Gauge+Higgs HMC: D=%d SU(%d) rep=%s d=%d N-ality=%d L=%d^%d vol=%lld\n",
              kDim, kN, rep->name().c_str(), rep->d, rep->nality, Lext, kDim, (long long)hmc.lat.vol);
  std::printf("# beta=%.4f kappa=%.4f lambda=%.4f tau=%.3f nmd=%d  (kappa_c,free=1/(2D)=%.4f)\n",
              beta, kappa, lambda, tau, nmd, 1.0 / (2.0 * kDim));

  hmc.U.hot(hmc.rng, 0.8);
  hmc.phi.gaussian(hmc.rng, 12345, rep->real, 0.3);
  for (int t = 0; t < ntherm; ++t) hmc.trajectory();

  hmc.traj_count = 0; hmc.accept_count = 0;
  Stats plaq, Lphi, Llink, poly;
  double sExp = 0.0;
  const double V = static_cast<double>(hmc.lat.vol);
  for (int t = 0; t < nmeas; ++t) {
    hmc.trajectory();
    plaq.add(avg_plaquette<kDim, kN>(hmc.U));
    Lphi.add(higgs_length<kDim>(hmc.phi));
    Llink.add(link_energy<kDim, kN>(hmc.phi, hmc.U, *rep));
    poly.add(polyakov_loop<kDim, kN>(hmc.U));
    sExp += std::exp(-hmc.last_dH);
  }
  std::printf("plaquette   = %.6f +/- %.6f\n", plaq.mean(), plaq.binned_error());
  std::printf("L_phi       = %.6f +/- %.6f   (<phi^dag phi>)\n", Lphi.mean(), Lphi.binned_error());
  std::printf("L_link      = %.6f +/- %.6f   (gauge-inv hopping energy; transition locator)\n",
              Llink.mean(), Llink.binned_error());
  std::printf("Polyakov    = %.6f +/- %.6f\n", poly.mean(), poly.binned_error());
  std::printf("chi_link    = %.5f   (V*Var(L_link), susceptibility)\n", Llink.susceptibility(V));
  std::printf("Binder_link = %.5f   (U4 of L_link)\n", Llink.binder());
  std::printf("chi_phi     = %.5f   (V*Var(L_phi))\n", Lphi.susceptibility(V));
  std::printf("acceptance  = %.4f\n", hmc.acceptance());
  std::printf("<exp(-dH)>  = %.5f\n", sExp / nmeas);
  std::printf("max||UdU-I||= %.2e\n", hmc.U.max_unitarity_violation());
  return 0;
}
