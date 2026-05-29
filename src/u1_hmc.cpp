// Compact U(1) + charge-q Higgs HMC driver (abelian Fradkin-Shenker: charge-q Higgs
// condensation breaks U(1) -> Z_q). Compile-time D = NDIM.
//   ./build/u1_hmc <L> <beta> <kappa> <lambda> <q> [ntherm=80 nmeas=200 nmd=20 tau=1 seed=1]
// Prints plaquette, L_phi, charge-q link energy, charge-m Wilson loops / Polyakov, and the
// HMC diagnostics. A charge-m source is screened iff q | m (the residual Z_q is blind to it).
#include "u1/u1.hpp"
#include "measure/observables.hpp"   // Stats
#include <cstdio>
#include <cstdlib>

using namespace gh;

int main(int argc, char** argv) {
  if (argc < 6) {
    std::fprintf(stderr, "usage: %s <L> <beta> <kappa> <lambda> <q> [ntherm nmeas nmd tau seed]\n", argv[0]);
    return 1;
  }
  auto af = [&](int i, double d) { return i < argc ? std::atof(argv[i]) : d; };
  auto ai = [&](int i, long d)   { return i < argc ? std::atol(argv[i]) : d; };
  const int  L      = static_cast<int>(ai(1, 8));
  const Real beta   = af(2, 1.5);
  const Real kappa  = af(3, 0.2);
  const Real lambda = af(4, 0.5);
  const int  q      = static_cast<int>(ai(5, 2));
  const int  ntherm = static_cast<int>(ai(6, 80));
  const int  nmeas  = static_cast<int>(ai(7, 200));
  const int  nmd    = static_cast<int>(ai(8, 20));
  const Real tau    = af(9, 1.0);
  const std::uint64_t seed = static_cast<std::uint64_t>(ai(10, 1));

  std::array<int, kDim> ext{}; for (int mu = 0; mu < kDim; ++mu) ext[mu] = L;
  u1::U1HMC<kDim> hmc(ext, seed);
  hmc.beta = beta; hmc.kappa = kappa; hmc.lambda = lambda; hmc.q = q; hmc.tau = tau; hmc.nmd = nmd;
  std::printf("# U(1)+charge-%d Higgs HMC: D=%d L=%d^%d  beta=%.3f kappa=%.3f lambda=%.3f nmd=%d  (breaks U(1)->Z_%d)\n",
              q, kDim, L, kDim, beta, kappa, lambda, nmd, q);

  hmc.hot(0.8); hmc.cold_phi(0.5);
  for (int t = 0; t < ntherm; ++t) hmc.trajectory();
  hmc.traj_count = 0; hmc.accept_count = 0;

  Stats plaq, Lphi, Llink; std::vector<Stats> wl(q + 1);  // charge-1..q Wilson(2x2)
  double sExp = 0.0;
  for (int t = 0; t < nmeas; ++t) {
    hmc.trajectory();
    plaq.add(u1::avg_plaquette<kDim>(hmc.th, hmc.lat));
    Lphi.add(u1::higgs_length<kDim>(hmc.phi, hmc.lat));
    Llink.add(u1::link_energy<kDim>(hmc.phi, hmc.th, hmc.lat, q));
    for (int m = 1; m <= q; ++m) wl[m].add(u1::wilson_loop<kDim>(hmc.th, hmc.lat, m, 2, 2));
    sExp += std::exp(-hmc.last_dH);
  }
  std::printf("plaquette   = %.6f +/- %.6f\n", plaq.mean(), plaq.binned_error());
  std::printf("L_phi       = %.6f +/- %.6f\n", Lphi.mean(), Lphi.binned_error());
  std::printf("L_link(q=%d) = %.6f +/- %.6f\n", q, Llink.mean(), Llink.binned_error());
  for (int m = 1; m <= q; ++m)
    std::printf("W_%d(2x2)    = %.6f +/- %.6f   %s\n", m, wl[m].mean(), wl[m].binned_error(),
                (m % q == 0) ? "(Z_q-neutral: screened)" : "(Z_q-charged)");
  std::printf("acceptance  = %.4f   <exp(-dH)>=%.5f\n", hmc.acceptance(), sExp / nmeas);
  return 0;
}
