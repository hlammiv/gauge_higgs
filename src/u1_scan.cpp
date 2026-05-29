// Compact U(1) + charge-q Higgs (beta,kappa) GRID SCAN driver. Compile-time D = NDIM.
// For each grid node it runs the validated U1HMC, then writes the per-trajectory time
// series the reweighter (reweight.hpp) consumes, plus a per-node summary row.
//
//   ./build/u1_scan <L> <bmin> <bmax> <nb> <kmin> <kmax> <nk> <lambda> <q>
//                   [ntherm nmeas nmd tau base_seed measure_every outdir]
//
// REWEIGHTING CONVENTION (must match reweight.hpp EXACTLY):
//   weight = exp(-S),  S = beta*A - kappa*B + (on-site scalar terms), where
//     A = sum_plaq (1 - cos theta_plaq)                                  [EXTENSIVE]
//         = gauge_action(theta,lat,1.0);  conjugate to beta.
//     B = sum_{x,mu} 2 Re[ conj(phi_x) e^{i q theta_mu(x)} phi_{x+mu} ]  [EXTENSIVE]
//         scalar_action contains -kappa*B, so the variable conjugate to kappa is -B.
//   reweight.hpp form S = sum_i lambda_i E_i:
//     (E_1, E_2) = (A, -B),  (lambda_1, lambda_2) = (beta, kappa).
//   The time series stores A (plaq_energy_sum) and B (hop_energy_sum) per trajectory;
//   these are EXACTLY the (beta,kappa)-conjugate energies, so reweighting is exact.
#include "u1/u1.hpp"
#include "u1/scan_obs.hpp"
#include "u1/monopole.hpp"           // monopole_density<D> (DeGrand-Toussaint)
#include "u1/autotune.hpp"           // tune_nmd: per-point nmd auto-tune
#include "measure/observables.hpp"   // Stats
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>

using namespace gh;

int main(int argc, char** argv) {
  if (argc < 10) {
    std::fprintf(stderr,
      "usage: %s <L> <bmin> <bmax> <nb> <kmin> <kmax> <nk> <lambda> <q> "
      "[ntherm nmeas nmd tau base_seed measure_every outdir autotune]\n", argv[0]);
    return 1;
  }
  auto af = [&](int i, double d) { return i < argc ? std::atof(argv[i]) : d; };
  auto ai = [&](int i, long d)   { return i < argc ? std::atol(argv[i]) : d; };
  const int    L      = static_cast<int>(ai(1, 8));
  const Real   bmin   = af(2, 1.0);
  const Real   bmax   = af(3, 2.0);
  const int    nb     = static_cast<int>(ai(4, 5));
  const Real   kmin   = af(5, 0.0);
  const Real   kmax   = af(6, 0.4);
  const int    nk     = static_cast<int>(ai(7, 5));
  const Real   lambda = af(8, 0.5);
  const int    q      = static_cast<int>(ai(9, 2));
  const int    ntherm = static_cast<int>(ai(10, 80));
  const int    nmeas  = static_cast<int>(ai(11, 200));
  const int    nmd    = static_cast<int>(ai(12, 20));
  const Real   tau    = af(13, 1.0);
  const std::uint64_t base_seed = static_cast<std::uint64_t>(ai(14, 1));
  const int    measure_every = static_cast<int>(ai(15, 1));
  const std::string outdir = (16 < argc) ? std::string(argv[16]) : std::string("u1scan_out");
  // autotune (default ON): per-point, raise nmd from the given starting value until
  // acceptance is in band -- fixes the stiff deep-Higgs (large kappa/q) points.
  const int autotune = static_cast<int>(ai(17, 1));

  // Create the output directory if missing (idempotent; ignore "already exists").
  if (::mkdir(outdir.c_str(), 0755) != 0) { /* likely already exists -- proceed */ }

  std::array<int, kDim> ext{}; for (int mu = 0; mu < kDim; ++mu) ext[mu] = L;
  const Lattice<kDim> probe(ext);
  const std::int64_t n_links = probe.vol * kDim;
  const std::int64_t n_plaq  = probe.n_plaq();

  // summary.csv (one row per node); header documents columns + reweighting convention.
  const std::string sumpath = outdir + "/summary.csv";
  FILE* sf = std::fopen(sumpath.c_str(), "w");
  if (!sf) { std::fprintf(stderr, "ERROR: cannot open %s\n", sumpath.c_str()); return 1; }
  std::fprintf(sf,
    "# U(1)+charge-%d Higgs (beta,kappa) scan summary. D=%d L=%d^%d lambda=%.6g\n",
    q, kDim, L, kDim, lambda);
  std::fprintf(sf,
    "# Reweighting: weight=exp(-S), S=beta*A-kappa*B+on-site; A=sum_plaq(1-cos),"
    " B=sum_{x,mu}2Re[conj(phi)e^{iq theta}phi]; (E1,E2)=(A,-B),(l1,l2)=(beta,kappa).\n");
  // rho_M (DeGrand-Toussaint monopole density) is appended at the END of each row, so the
  // legacy column layout is untouched and old parsers keep working.
  std::fprintf(sf, "beta,kappa,L,q,plaq,plaq_err,Llink,Llink_err,phi2,phi2_err,acceptance,exp_mdH,rho_M,rho_M_err\n");

  std::fprintf(stderr, "# u1_scan: D=%d L=%d^%d q=%d lambda=%.4g grid=%dx%d (beta in [%.4g,%.4g], kappa in [%.4g,%.4g])\n",
               kDim, L, kDim, q, lambda, nb, nk, bmin, bmax, kmin, kmax);
  std::fprintf(stderr, "# ntherm=%d nmeas=%d nmd=%d(start) tau=%.4g base_seed=%llu measure_every=%d autotune=%d outdir=%s\n",
               ntherm, nmeas, nmd, tau, (unsigned long long)base_seed, measure_every, autotune, outdir.c_str());

  for (int ib = 0; ib < nb; ++ib) {
    const Real beta = (nb <= 1) ? bmin : bmin + (bmax - bmin) * ib / (nb - 1);
    for (int ik = 0; ik < nk; ++ik) {
      const Real kappa = (nk <= 1) ? kmin : kmin + (kmax - kmin) * ik / (nk - 1);
      const int node = ib * nk + ik;
      // Deterministic, node-independent seed: mix base_seed with the node index so the
      // scan is reproducible and each node's stream is decoupled from grid ordering.
      const std::uint64_t seed = Rng::key(base_seed, node);

      u1::U1HMC<kDim> hmc(ext, seed);
      hmc.beta = beta; hmc.kappa = kappa; hmc.lambda = lambda; hmc.q = q; hmc.tau = tau; hmc.nmd = nmd;
      hmc.hot(0.8); hmc.cold_phi(0.5);
      for (int t = 0; t < ntherm; ++t) hmc.trajectory();
      if (autotune) {
        const u1::TuneResult tr = u1::tune_nmd<kDim>(hmc);  // extra thermalization + sets hmc.nmd
        std::fprintf(stderr, "  node %d: autotuned nmd=%d (cal acc=%.3f in_band=%d)\n",
                     node, hmc.nmd, tr.acceptance, (int)tr.in_band);
      }
      hmc.traj_count = 0; hmc.accept_count = 0;

      // Per-node time series file.
      char fname[256];
      std::snprintf(fname, sizeof fname, "%s/ts_b%.6f_k%.6f.dat", outdir.c_str(), beta, kappa);
      FILE* tf = std::fopen(fname, "w");
      if (!tf) { std::fprintf(stderr, "ERROR: cannot open %s\n", fname); std::fclose(sf); return 1; }
      std::fprintf(tf,
        "# U(1)+charge-%d Higgs time series. D=%d L=%d^%d beta=%.6f kappa=%.6f lambda=%.6g q=%d\n",
        q, kDim, L, kDim, beta, kappa, lambda, q);
      std::fprintf(tf,
        "# Reweighting: weight=exp(-S), S=beta*A-kappa*B+on-site;"
        " A=plaq_energy_sum=sum_plaq(1-cos theta_pl) [conj to beta],"
        " B=hop_energy_sum=sum_{x,mu}2Re[conj(phi)e^{iq theta}phi] [-B conj to kappa].\n");
      std::fprintf(tf, "# columns: traj  A  B  avg_plaquette  higgs_length  link_energy  monopole_density\n");

      Stats plaq, Lphi, Llink, rhoM;
      double sExp = 0.0; int nrows = 0;
      for (int t = 0; t < nmeas; ++t) {
        for (int e = 0; e < measure_every; ++e) hmc.trajectory();
        const Real A   = u1::plaq_energy_sum<kDim>(hmc.th, hmc.lat);
        const Real B   = u1::hop_energy_sum<kDim>(hmc.phi, hmc.th, hmc.lat, q);
        const Real pl  = u1::avg_plaquette<kDim>(hmc.th, hmc.lat);
        const Real lp  = u1::higgs_length<kDim>(hmc.phi, hmc.lat);
        const Real le  = u1::link_energy<kDim>(hmc.phi, hmc.th, hmc.lat, q);
        const Real rho = u1::monopole_density<kDim>(hmc.th, hmc.lat);
        std::fprintf(tf, "%d %.15g %.15g %.15g %.15g %.15g %.15g\n", t, A, B, pl, lp, le, rho);
        plaq.add(pl); Lphi.add(lp); Llink.add(le); rhoM.add(rho);
        sExp += std::exp(-hmc.last_dH);
        ++nrows;
      }
      std::fclose(tf);

      const double expm = nrows ? sExp / nrows : 0.0;
      std::fprintf(sf, "%.6f,%.6f,%d,%d,%.8g,%.3g,%.8g,%.3g,%.8g,%.3g,%.4f,%.6g,%.8g,%.3g\n",
                   beta, kappa, L, q,
                   plaq.mean(), plaq.binned_error(),
                   Llink.mean(), Llink.binned_error(),
                   Lphi.mean(), Lphi.binned_error(),
                   hmc.acceptance(), expm,
                   rhoM.mean(), rhoM.binned_error());
      std::fflush(sf);

      std::fprintf(stderr,
        "  node %d/%d  beta=%.4f kappa=%.4f  acc=%.3f <exp(-dH)>=%.4f  plaq=%.4f Llink=%.4f Lphi=%.4f  -> %s\n",
        node + 1, nb * nk, beta, kappa, hmc.acceptance(), expm,
        plaq.mean(), Llink.mean(), Lphi.mean(), fname);
    }
  }
  std::fclose(sf);
  // Silence unused-variable warnings for the documentation-only extensive counts.
  (void)n_links; (void)n_plaq;
  std::fprintf(stderr, "# done: %d nodes -> %s\n", nb * nk, sumpath.c_str());
  return 0;
}
