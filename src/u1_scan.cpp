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
#include "u1/monopole.hpp"           // monopole_density<D> (DeGrand-Toussaint); also defines reduced_plaq_angle<D>
#include "u1/autotune.hpp"           // tune_nmd: per-point nmd auto-tune
// monopole.hpp already provides an identical reduced_plaq_angle<D>; tell photon_mass.hpp
// to reuse it instead of redefining it (avoids an in-TU template redefinition / ODR error).
#define GH_U1_HAVE_REDUCED_PLAQ_ANGLE
#include "u1/photon_mass.hpp"        // photon_timeslice_field<D>: transverse photon field per timeslice
#include "measure/observables.hpp"   // Stats
#include "measure/correlator.hpp"    // cosh_effective_mass, plateau (jackknife) -> m_gamma
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <string>
#include <vector>
#include <utility>
#include <algorithm>
#include <sys/stat.h>
#include <sys/types.h>

using namespace gh;

int main(int argc, char** argv) {
  if (argc < 10) {
    std::fprintf(stderr,
      "usage: %s <L> <bmin> <bmax> <nb> <kmin> <kmax> <nk> <lambda> <q> "
      "[ntherm nmeas nmd tau base_seed measure_every outdir autotune n_scalar]\n", argv[0]);
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
  // n_scalar (default 1 = single-timescale): scalar(+matter) sub-steps per gauge step
  // (multi-timescale Sexton-Weingarten). >1 resolves the stiff scalar sector cheaply
  // so the expensive gauge force stays on the coarse nmd timescale.
  const int n_scalar = static_cast<int>(ai(18, 1));

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
  // rho_M (DeGrand-Toussaint monopole density) and the photon observables are appended
  // at the END of each row (legacy column layout untouched). m_gamma is the Coulomb-
  // confinement LOCATOR: the contact-to-first-step log-ratio log(C(0)/|C(1)|) of the
  // transverse photon correlator -- SMALL when the photon is light (Coulomb), LARGE
  // when heavy (confined AND Higgs). It is robust on short Lt but is NOT a calibrated
  // mass (the heavy photon is sub-lattice on Lt=8; a precise m_gamma needs larger L_t).
  // m_gamma_cosh is the cosh-plateau effective mass, meaningful ONLY in the light/
  // Coulomb regime (it returns NaN/0 on the noisy heavy-phase tail -- do NOT use it as
  // the locator). The full correlator C(dt)+/-err is in cor_*.dat for re-extraction.
  std::fprintf(sf, "beta,kappa,L,q,plaq,plaq_err,Llink,Llink_err,phi2,phi2_err,acceptance,exp_mdH,rho_M,rho_M_err,m_gamma,m_gamma_err,m_gamma_cosh,m_gamma_cosh_err\n");

  std::fprintf(stderr, "# u1_scan: D=%d L=%d^%d q=%d lambda=%.4g grid=%dx%d (beta in [%.4g,%.4g], kappa in [%.4g,%.4g])\n",
               kDim, L, kDim, q, lambda, nb, nk, bmin, bmax, kmin, kmax);
  std::fprintf(stderr, "# ntherm=%d nmeas=%d nmd=%d(start) tau=%.4g base_seed=%llu measure_every=%d autotune=%d n_scalar=%d outdir=%s\n",
               ntherm, nmeas, nmd, tau, (unsigned long long)base_seed, measure_every, autotune, n_scalar, outdir.c_str());

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
      hmc.n_scalar = n_scalar;
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
      // Per-config transverse photon correlators C_k(dt). After the loop these feed
      // gh::plateau (delete-1 jackknife over configs) to extract the photon mass.
      std::vector<std::vector<Real>> photon_samples;
      photon_samples.reserve(nmeas);
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

        // Zero-spatial-momentum transverse photon field, recomputed for THIS config:
        // S[t][i] = sum_{x:x[0]=t} F_{0 i}(x), i=0..D-2 over the D-1 spatial dirs.
        // Build this config's own correlator C_k(dt) = (1/(Lt*ncomp)) *
        // sum_{t,i} S[(t+dt)%Lt][i]*S[t][i] for dt=0..Lt-1 (periodic in time).
        auto S = u1::photon_timeslice_field<kDim>(hmc.th, hmc.lat);
        const int Lt    = static_cast<int>(S.size());
        const int ncomp = Lt ? static_cast<int>(S[0].size()) : 0;  // = D-1
        const Real norm = static_cast<Real>(Lt) * static_cast<Real>(ncomp);
        std::vector<Real> Ck(Lt, 0.0);
        for (int dt = 0; dt < Lt; ++dt) {
          Real acc = 0.0;
          for (int tt = 0; tt < Lt; ++tt) {
            const int tp = (tt + dt) % Lt;
            for (int i = 0; i < ncomp; ++i) acc += S[tp][i] * S[tt][i];
          }
          Ck[dt] = (norm > 0.0) ? acc / norm : 0.0;
        }
        photon_samples.push_back(std::move(Ck));

        sExp += std::exp(-hmc.last_dH);
        ++nrows;
      }
      std::fclose(tf);

      // --- Transverse photon observables (Coulomb-confinement locator) -----------
      // PRIMARY m_gamma = contact-to-first-step log-ratio  log(C(0)/|C(1)|)  of the
      // mean correlator, with a delete-1 jackknife error over configs. SMALL when the
      // photon is light (Coulomb), LARGE when heavy (confined AND Higgs). Robust on
      // short Lt; NOT a calibrated mass. SECONDARY m_gamma_cosh = cosh-plateau eff.
      // mass over dt in [max(1,Lt/4), Lt/2] -- reliable ONLY in the light/Coulomb
      // regime (NaN/clamped on the noisy heavy-phase tail). The cosh estimator cannot
      // be the locator: it returns ~0 in the Higgs phase, which masquerades as massless.
      Real m_gamma = 0.0, m_gamma_err = 0.0;            // primary: first-step log-ratio
      Real m_gamma_cosh = 0.0, m_gamma_cosh_err = 0.0;  // secondary: cosh plateau
      const int Ltps = photon_samples.empty() ? 0 : static_cast<int>(photon_samples[0].size());
      const int Kc   = static_cast<int>(photon_samples.size());
      if (Ltps >= 2 && Kc >= 1) {
        // Mean correlator Cbar(dt) and its per-dt standard error (delete-1 jackknife
        // of a mean == SEM).
        std::vector<Real> tot(Ltps, 0.0), Cbar(Ltps, 0.0), Cerr(Ltps, 0.0);
        for (const auto& Ck : photon_samples)
          for (int dt = 0; dt < Ltps; ++dt) tot[dt] += Ck[dt];
        for (int dt = 0; dt < Ltps; ++dt) Cbar[dt] = tot[dt] / Kc;
        if (Kc >= 2) {
          for (const auto& Ck : photon_samples)
            for (int dt = 0; dt < Ltps; ++dt) { const Real d = Ck[dt] - Cbar[dt]; Cerr[dt] += d * d; }
          for (int dt = 0; dt < Ltps; ++dt)
            Cerr[dt] = std::sqrt(Cerr[dt] / (static_cast<Real>(Kc) * (Kc - 1)));
        }

        // First-step log-ratio discriminator + delete-1 jackknife. Floor |C(1)| so the
        // ratio stays finite (a vanishing tail -> large m_gamma, i.e. a heavy photon).
        const Real eps = 1e-300;
        auto first_step = [&](Real c0, Real c1) -> Real {
          const Real a = (c0 > 0.0) ? c0 : eps;
          const Real b = std::max(std::fabs(c1), eps);
          return std::log(a / b);
        };
        m_gamma = first_step(Cbar[0], Cbar[1]);
        if (Kc >= 2) {
          Real jbar = 0.0; std::vector<Real> jk(Kc);
          for (int k = 0; k < Kc; ++k) {
            const Real c0 = (tot[0] - photon_samples[k][0]) / (Kc - 1);
            const Real c1 = (tot[1] - photon_samples[k][1]) / (Kc - 1);
            jk[k] = first_step(c0, c1); jbar += jk[k];
          }
          jbar /= Kc;
          Real sw = 0.0; for (Real v : jk) sw += (v - jbar) * (v - jbar);
          m_gamma = jbar;
          m_gamma_err = std::sqrt((static_cast<Real>(Kc - 1) / Kc) * sw);
        }

        // Secondary cosh-plateau effective mass (light/Coulomb regime only).
        if (Ltps >= 3) {
          const int tmin = std::max(1, Ltps / 4);
          const int tmax = Ltps / 2;
          const gh::PlateauFit pf = gh::plateau(photon_samples, tmin, tmax);
          m_gamma_cosh = pf.mass; m_gamma_cosh_err = pf.err;
        }

        // Per-node correlator file: dt  Cbar(dt)  Cbar_err(dt). Storing the full
        // correlator + errors lets the analyzer re-extract ANY photon estimator with
        // errors WITHOUT re-running the (expensive) scan.
        char cname[256];
        std::snprintf(cname, sizeof cname, "%s/cor_b%.6f_k%.6f.dat", outdir.c_str(), beta, kappa);
        FILE* cf = std::fopen(cname, "w");
        if (cf) {
          std::fprintf(cf,
            "# U(1)+charge-%d transverse photon (zero-spatial-momentum F_{0i}) correlator."
            " D=%d L=%d^%d beta=%.6f kappa=%.6f lambda=%.6g q=%d\n",
            q, kDim, L, kDim, beta, kappa, lambda, q);
          std::fprintf(cf,
            "# Cbar(dt)=avg over %d configs of (1/(Lt*(D-1))) sum_{t,i} S_i(t+dt)S_i(t);"
            " Cbar_err=delete-1 jackknife (SEM). m_gamma=log(C0/|C1|)=%.8g+/-%.3g"
            " (Coulomb locator: small=light=Coulomb); m_gamma_cosh=%.8g+/-%.3g (light-only).\n",
            nrows, m_gamma, m_gamma_err, m_gamma_cosh, m_gamma_cosh_err);
          std::fprintf(cf, "# columns: dt  Cbar(dt)  Cbar_err(dt)\n");
          for (int dt = 0; dt < Ltps; ++dt)
            std::fprintf(cf, "%d %.15g %.15g\n", dt, Cbar[dt], Cerr[dt]);
          std::fclose(cf);
        }
      }

      const double expm = nrows ? sExp / nrows : 0.0;
      std::fprintf(sf, "%.6f,%.6f,%d,%d,%.8g,%.3g,%.8g,%.3g,%.8g,%.3g,%.4f,%.6g,%.8g,%.3g,%.8g,%.3g,%.8g,%.3g\n",
                   beta, kappa, L, q,
                   plaq.mean(), plaq.binned_error(),
                   Llink.mean(), Llink.binned_error(),
                   Lphi.mean(), Lphi.binned_error(),
                   hmc.acceptance(), expm,
                   rhoM.mean(), rhoM.binned_error(),
                   m_gamma, m_gamma_err, m_gamma_cosh, m_gamma_cosh_err);
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
