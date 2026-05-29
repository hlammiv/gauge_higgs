// MPI domain-decomposed pure-gauge HMC driver. Compile-time D=NDIM, N=NCOL.
//   mpirun -np <P> ./build/gh_hmc_mpi <L> <procgrid> <beta> [ntherm] [nmeas] [nmd] [tau] [seed]
//     <procgrid> = comma list, e.g. "2,1,1,1" (product must equal -np; each must divide L)
// Decomposition-independence (global-index RNG) means -np 1/2/4 give the SAME plaquette,
// which is how correctness is validated cheaply. Keep test runs small (few ranks, 1 thread).
#include "mpi/gauge_mpi.hpp"
#include <cstdio>
#include <cstdlib>
#include <string>
#include <sstream>
#include <vector>

using namespace gh;

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  int rank = 0, nranks = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nranks);

  if (argc < 4) {
    if (rank == 0) std::fprintf(stderr,
      "usage: mpirun -np P %s <L> <procgrid e.g. 2,1,1,1> <beta> [ntherm=20] [nmeas=50] [nmd=20] [tau=1.0] [seed=1]\n",
      argv[0]);
    MPI_Finalize(); return 1;
  }
  auto argf = [&](int i, double d) { return i < argc ? std::atof(argv[i]) : d; };
  auto argi = [&](int i, long d)   { return i < argc ? std::atol(argv[i]) : d; };
  const int  Lext   = static_cast<int>(argi(1, 8));
  const Real beta   = argf(3, 2.3);
  const int  ntherm = static_cast<int>(argi(4, 20));
  const int  nmeas  = static_cast<int>(argi(5, 50));
  const int  nmd    = static_cast<int>(argi(6, 20));
  const Real tau    = argf(7, 1.0);
  const std::uint64_t seed = static_cast<std::uint64_t>(argi(8, 1));

  std::array<int, kDim> G{}; for (int mu = 0; mu < kDim; ++mu) G[mu] = Lext;
  std::array<int, kDim> P{}; P.fill(1);
  { std::stringstream ss(argv[2]); std::string t; int mu = 0;
    while (std::getline(ss, t, ',') && mu < kDim) if (!t.empty()) P[mu++] = std::atoi(t.c_str()); }

  int prod = 1; for (int mu = 0; mu < kDim; ++mu) prod *= P[mu];
  bool ok = (prod == nranks);
  for (int mu = 0; mu < kDim; ++mu) ok = ok && (G[mu] % P[mu] == 0);
  if (!ok) {
    if (rank == 0) std::fprintf(stderr, "error: prod(procgrid)=%d must equal -np=%d and each P must divide L=%d\n",
                                prod, nranks, Lext);
    MPI_Finalize(); return 1;
  }

  DistGaugeHMC<kDim, kN> hmc(G, P, seed);
  hmc.beta = beta; hmc.tau = tau; hmc.nmd = nmd;
  hmc.U.hot(hmc.rng, 0.8);

  halo_exchange(hmc.U);
  const Real plaq_hot = avg_plaquette_d<kDim, kN>(hmc.U);
  if (rank == 0) {
    std::printf("# MPI pure-gauge HMC: D=%d SU(%d) L=%d^%d  ranks=%d  procgrid=", kDim, kN, Lext, kDim, nranks);
    for (int mu = 0; mu < kDim; ++mu) std::printf("%d%s", P[mu], mu + 1 < kDim ? "x" : "");
    std::printf("  beta=%.4f nmd=%d\n", beta, nmd);
    std::printf("plaquette(hot start)   = %.10f   <- must be identical for any procgrid (decomposition test)\n", plaq_hot);
  }

  for (int t = 0; t < ntherm; ++t) hmc.trajectory();
  hmc.traj_count = 0; hmc.accept_count = 0;
  double sExp = 0.0; Real sPl = 0.0;
  for (int t = 0; t < nmeas; ++t) {
    hmc.trajectory();
    halo_exchange(hmc.U);
    sPl  += avg_plaquette_d<kDim, kN>(hmc.U);
    sExp += std::exp(-hmc.last_dH);
  }
  if (rank == 0) {
    std::printf("plaquette(measured)    = %.10f   (after %d therm + %d meas trajectories)\n", sPl / nmeas, ntherm, nmeas);
    std::printf("acceptance             = %.4f\n", hmc.acceptance());
    std::printf("<exp(-dH)>             = %.6f\n", sExp / nmeas);
  }
  MPI_Finalize();
  return 0;
}
