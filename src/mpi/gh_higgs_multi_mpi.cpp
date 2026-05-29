// MPI domain-decomposed gauge + Higgs HMC with the GENERAL multi-invariant scalar potential.
//   V = -mu2 (phi^dag phi) + sum_c f_c V_c[phi]   (V_c = Casimir-channel projection of phi phi^dag)
// Models src/hmc_higgs_multi.cpp (serial) on top of the distributed DistGaugeHiggsHMC, whose
// hmc.potential (OnsitePotential) hook is wired through hamiltonian()/kick().
//
//   mpirun -np P ./build_mpi/gh_higgs_multi_mpi <rep> <L> <procgrid> <beta> <kappa> <mu2> <couplings|auto>
//                                               [ntherm=40] [nmeas=80] [nmd=20] [tau=1] [seed=1]
//     <rep>       = fund | adj | <Young rows e.g. 6 (=SU(2) spin-3)>[:real]
//     <procgrid>  = e.g. 2,1,1,1
//     <couplings> = comma list f0,f1,... (one per channel, printed C2 order), or "auto" (all f_c=1).
// Decomposition-independence (global-index RNG): -np 1/2/4 give identical plaquette, L_link, L_phi.
// Keep test runs small (few ranks, 1 thread).
#include "mpi/higgs_mpi.hpp"
#include "action/scalar_invariants.hpp"
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
  if (rows.empty()) throw std::runtime_error("bad rep spec");
  return std::make_unique<GeneralRep<kN>>(rows, real ? GeneralRep<kN>::RealType::Real : GeneralRep<kN>::RealType::Complex);
}

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  int rank = 0, nranks = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nranks);

  if (argc < 8) {
    if (rank == 0) std::fprintf(stderr,
      "usage: mpirun -np P %s <rep> <L> <procgrid> <beta> <kappa> <mu2> <couplings|auto> [ntherm=40] [nmeas=80] [nmd=20] [tau=1] [seed=1]\n"
      "  <couplings> = comma list f_c (one per channel, in the printed C2 order), or 'auto' (all f_c=1).\n",
      argv[0]);
    MPI_Finalize(); return 1;
  }
  const std::string repspec = argv[1];
  auto argf = [&](int i, double d) { return i < argc ? std::atof(argv[i]) : d; };
  auto argi = [&](int i, long d)   { return i < argc ? std::atol(argv[i]) : d; };
  const int  Lext   = static_cast<int>(argi(2, 8));
  const Real beta   = argf(4, 2.3);
  const Real kappa  = argf(5, 0.2);
  const Real mu2    = argf(6, 1.0);
  const int  ntherm = static_cast<int>(argi(8, 40));
  const int  nmeas  = static_cast<int>(argi(9, 80));
  const int  nmd    = static_cast<int>(argi(10, 20));
  const Real tau    = argf(11, 1.0);
  const std::uint64_t seed = static_cast<std::uint64_t>(argi(12, 1));

  std::array<int, kDim> G{}; for (int mu = 0; mu < kDim; ++mu) G[mu] = Lext;
  std::array<int, kDim> P{}; P.fill(1);
  { std::stringstream ss(argv[3]); std::string t; int mu = 0;
    while (std::getline(ss, t, ',') && mu < kDim) if (!t.empty()) P[mu++] = std::atoi(t.c_str()); }
  int prod = 1; for (int mu = 0; mu < kDim; ++mu) prod *= P[mu];
  bool ok = (prod == nranks);
  for (int mu = 0; mu < kDim; ++mu) ok = ok && (G[mu] % P[mu] == 0);
  if (!ok) { if (rank == 0) std::fprintf(stderr, "error: prod(procgrid)=%d must equal -np=%d and divide L=%d\n", prod, nranks, Lext); MPI_Finalize(); return 1; }

  std::unique_ptr<Representation<kN>> rep;
  try { rep = make_rep(repspec); }
  catch (const std::exception& e) { if (rank == 0) std::fprintf(stderr, "rep error: %s\n", e.what()); MPI_Finalize(); return 1; }

  CasimirChannels<kN> ch(*rep);
  if (rank == 0) {
    std::printf("# rep=%s d=%d  %d quartic channels (C2): ", rep->name().c_str(), rep->d, ch.n_channels());
    for (Real c : ch.lambda) std::printf("%.4g ", c);
    std::printf("\n");
  }

  // Couplings f_c: comma list (one per channel, printed C2 order) or "auto" (all 1).
  std::vector<Real> f;
  if (std::string(argv[7]) == "auto") f.assign(ch.n_channels(), 1.0);
  else { std::stringstream ss(argv[7]); std::string t; while (std::getline(ss, t, ',')) if (!t.empty()) f.push_back(std::atof(t.c_str())); }
  if (static_cast<int>(f.size()) != ch.n_channels()) {
    if (rank == 0) std::fprintf(stderr, "error: need %d couplings, got %zu\n", ch.n_channels(), f.size());
    MPI_Finalize(); return 1;
  }
  MultiInvariantPotential<kN> pot(ch, f, mu2);

  DistGaugeHiggsHMC<kDim, kN> hmc(G, P, *rep, seed);
  hmc.beta = beta; hmc.kappa = kappa; hmc.tau = tau; hmc.nmd = nmd; hmc.potential = &pot;
  hmc.U.hot(hmc.rng, 0.8);
  hmc.phi.gaussian(hmc.rng, 12345, rep->real, 0.3);

  halo_exchange(hmc.U); halo_exchange(hmc.phi);
  const Real pl0 = avg_plaquette_d<kDim, kN>(hmc.U);
  const Real ll0 = link_energy_d<kDim, kN>(hmc.phi, hmc.U, *rep);
  const Real lp0 = higgs_length_d<kDim>(hmc.phi);
  if (rank == 0) {
    std::printf("# MPI gauge+Higgs (multi-invariant): D=%d SU(%d) rep=%s d=%d L=%d^%d ranks=%d grid=", kDim, kN, rep->name().c_str(), rep->d, Lext, kDim, nranks);
    for (int mu = 0; mu < kDim; ++mu) std::printf("%d%s", P[mu], mu + 1 < kDim ? "x" : "");
    std::printf("  beta=%.3f kappa=%.3f mu2=%.3f nmd=%d\n", beta, kappa, mu2, nmd);
    std::printf("hot:  plaq=%.10f  L_link=%.10f  L_phi=%.10f  <- identical for any procgrid\n", pl0, ll0, lp0);
  }

  for (int t = 0; t < ntherm; ++t) hmc.trajectory();
  hmc.traj_count = 0; hmc.accept_count = 0;
  double sPl = 0, sLl = 0, sLp = 0, sExp = 0;
  for (int t = 0; t < nmeas; ++t) {
    hmc.trajectory();
    halo_exchange(hmc.U); halo_exchange(hmc.phi);
    sPl += avg_plaquette_d<kDim, kN>(hmc.U);
    sLl += link_energy_d<kDim, kN>(hmc.phi, hmc.U, *rep);
    sLp += higgs_length_d<kDim>(hmc.phi);
    sExp += std::exp(-hmc.last_dH);
  }
  if (rank == 0) {
    const double inv = 1.0 / nmeas;
    std::printf("meas: plaq=%.10f  L_link=%.10f  L_phi=%.10f  (%d+%d traj)\n", sPl * inv, sLl * inv, sLp * inv, ntherm, nmeas);
    std::printf("acceptance=%.4f  <exp(-dH)>=%.6f\n", hmc.acceptance(), sExp * inv);
  }
  MPI_Finalize();
  return 0;
}
