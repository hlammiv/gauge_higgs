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
#include "u1/u1_replica_tempering.hpp"  // U1ReplicaTempering<D>: kappa-axis parallel tempering (U1_TEMPER mode)
#include "u1/scan_obs.hpp"
#include "u1/monopole.hpp"           // monopole_density<D> (DeGrand-Toussaint); also defines reduced_plaq_angle<D>
#include "u1/autotune.hpp"           // tune_nmd: per-point nmd auto-tune
// monopole.hpp already provides an identical reduced_plaq_angle<D>; tell photon_mass.hpp
// to reuse it instead of redefining it (avoids an in-TU template redefinition / ODR error).
#define GH_U1_HAVE_REDUCED_PLAQ_ANGLE
#include "u1/photon_mass.hpp"        // photon_timeslice_field<D>: transverse photon field per timeslice
#include "measure/observables.hpp"   // Stats
#include "measure/correlator.hpp"    // cosh_effective_mass, plateau (jackknife) -> m_gamma
#include "measure/creutz.hpp"        // creutz_ratio (via creutz_jack.hpp)
#include "measure/creutz_jack.hpp"   // CreutzJackknife, chi_diag_jack, plateau_sigma_jack, creutz_excl_reason, WGrid
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <string>
#include <vector>
#include <utility>
#include <memory>
#include <algorithm>
#include <sys/stat.h>
#include <sys/types.h>

using namespace gh;

namespace {

// Full charge-1 AND charge-q Wilson-loop grids in ONE sweep over the lattice and all
// mu<nu planes (symmetrized, exactly as u1::wilson_loop). g1[R][T]=<cos(theta_loop)>,
// gq[R][T]=<cos(q*theta_loop)> for 1<=R,T<=Rmax (row/col 0 unused, mirrors W[R][T]).
// Ported verbatim from the validated POC src/u1_zqprobe.cpp.
template <int D>
void wilson_grids(const std::vector<Real>& th, const Lattice<D>& lat, int q, int Rmax,
                  WGrid& g1, WGrid& gq) {
  for (auto& r : g1) std::fill(r.begin(), r.end(), 0.0);
  for (auto& r : gq) std::fill(r.begin(), r.end(), 0.0);
  std::int64_t count = 0;
  #pragma omp parallel
  {
    WGrid l1(Rmax + 1, std::vector<Real>(Rmax + 1, 0.0));
    WGrid lq(Rmax + 1, std::vector<Real>(Rmax + 1, 0.0));
    std::int64_t lcount = 0;
    #pragma omp for schedule(static) nowait
    for (std::int64_t x = 0; x < lat.vol; ++x)
      for (int mu = 0; mu < D; ++mu)
        for (int nu = mu + 1; nu < D; ++nu) {
          for (int R = 1; R <= Rmax; ++R)
            for (int T = 1; T <= Rmax; ++T) {
              Real loop = 0.0; std::int64_t s = x;
              for (int i = 0; i < R; ++i) { loop += th[s * D + mu]; s = lat.neighbor_fwd(s, mu); }
              for (int j = 0; j < T; ++j) { loop += th[s * D + nu]; s = lat.neighbor_fwd(s, nu); }
              for (int i = 0; i < R; ++i) { s = lat.neighbor_bwd(s, mu); loop -= th[s * D + mu]; }
              for (int j = 0; j < T; ++j) { s = lat.neighbor_bwd(s, nu); loop -= th[s * D + nu]; }
              l1[R][T] += std::cos(loop);
              lq[R][T] += std::cos(q * loop);
              ++lcount;
            }
        }
    #pragma omp critical
    {
      for (int R = 1; R <= Rmax; ++R)
        for (int T = 1; T <= Rmax; ++T) { g1[R][T] += l1[R][T]; gq[R][T] += lq[R][T]; }
      count += lcount;
    }
  }
  // count tallies every (R,T,plane,site); per-(R,T) sample count = count / (Rmax*Rmax).
  const Real per = static_cast<Real>(count) / static_cast<Real>(Rmax * Rmax);
  if (per > 0.0)
    for (int R = 1; R <= Rmax; ++R)
      for (int T = 1; T <= Rmax; ++T) { g1[R][T] /= per; gq[R][T] /= per; }
}

// |P_n| = |(1/Vsp) sum_x exp(i n line_t(x))| for ONE config (the charge-n Polyakov
// modulus). <|P_1|> + Vspatial*Var(|P_1|) is the cheap finite-V-biased cross-check on
// the same charge-1 source. Ported verbatim from src/u1_zqprobe.cpp.
template <int D>
Real polyakov_abs(const std::vector<Real>& th, const Lattice<D>& lat, int n, int tdir = D - 1) {
  const int Lt = lat.L[tdir];
  Real re = 0.0, im = 0.0; std::int64_t count = 0;
  std::array<int, D> xc{};
  for (std::int64_t s = 0; s < lat.vol; ++s) {
    lat.coords(s, xc); if (xc[tdir] != 0) continue;
    Real line = 0.0; std::int64_t cur = s;
    for (int t = 0; t < Lt; ++t) { line += th[cur * D + tdir]; cur = lat.neighbor_fwd(cur, tdir); }
    re += std::cos(n * line); im += std::sin(n * line); ++count;
  }
  if (!count) return 0.0;
  re /= count; im /= count;
  return std::sqrt(re * re + im * im);
}

// ---------------------------------------------------------------------------
// Per-node measurement, refactored out of the main loop so it can be applied to
// EITHER an independent grid node (default path) OR each rung of a tempered
// kappa-ladder (U1_TEMPER). One MeasureCtx owns the per-rung accumulation state
// + the per-config ts file; the caller advances the chain (hmc.trajectory() for
// an independent node, pt.step() for a tempered ladder) and calls record() once
// per measurement, then finalize() to emit the cor file + summary row.
//
// EVERYTHING here is byte-for-byte the same computation as the original inline
// per-node block; only the plumbing (a struct instead of locals) changed.
struct MeasureCtx {
  // Identity / geometry (captured at construction, used in headers + summary row).
  Real beta, kappa, lambda;
  int  L, q, kDimv;
  int  Rmax, n_block;
  Real Vspatial;
  int  nmeas_hint;

  FILE* tf = nullptr;          // per-config time series (ts_*.dat)
  char  fname[256];            // ts path (echoed in the stderr progress line)

  Stats plaq, Lphi, Llink, rhoM, absP1;
  CreutzJackknife jk1, jkq;    // charge-1 PROBE + charge-q screened CONTROL
  WGrid b1, bq, w1, wq;        // current Creutz block + per-config scratch grids
  int   bcount = 0;
  std::vector<std::vector<Real>> photon_samples;  // per-config transverse photon C_k(dt)
  double sExp = 0.0; int nrows = 0;

  MeasureCtx(const std::array<int, kDim>& /*ext*/, const std::string& outdir,
             Real beta_, Real kappa_, Real lambda_, int L_, int q_,
             int Rmax_, int n_block_, Real Vspatial_, int nmeas_)
      : beta(beta_), kappa(kappa_), lambda(lambda_), L(L_), q(q_), kDimv(kDim),
        Rmax(Rmax_), n_block(n_block_), Vspatial(Vspatial_), nmeas_hint(nmeas_),
        jk1(Rmax_), jkq(Rmax_),
        b1(Rmax_ + 1, std::vector<Real>(Rmax_ + 1, 0.0)),
        bq(Rmax_ + 1, std::vector<Real>(Rmax_ + 1, 0.0)),
        w1(Rmax_ + 1, std::vector<Real>(Rmax_ + 1, 0.0)),
        wq(Rmax_ + 1, std::vector<Real>(Rmax_ + 1, 0.0)) {
    photon_samples.reserve(nmeas_);
    std::snprintf(fname, sizeof fname, "%s/ts_b%.6f_k%.6f.dat", outdir.c_str(), beta, kappa);
    tf = std::fopen(fname, "w");
    if (tf) {
      std::fprintf(tf,
        "# U(1)+charge-%d Higgs time series. D=%d L=%d^%d beta=%.6f kappa=%.6f lambda=%.6g q=%d\n",
        q, kDim, L, kDim, beta, kappa, lambda, q);
      std::fprintf(tf,
        "# Reweighting: weight=exp(-S), S=beta*A-kappa*B+on-site;"
        " A=plaq_energy_sum=sum_plaq(1-cos theta_pl) [conj to beta],"
        " B=hop_energy_sum=sum_{x,mu}2Re[conj(phi)e^{iq theta}phi] [-B conj to kappa].\n");
      std::fprintf(tf, "# columns: traj  A  B  avg_plaquette  higgs_length  link_energy  monopole_density\n");
    }
  }

  bool ok() const { return tf != nullptr; }

  void flush_probe() {
    if (!bcount) return;
    WGrid g1(Rmax + 1, std::vector<Real>(Rmax + 1, 0.0));
    WGrid gq(Rmax + 1, std::vector<Real>(Rmax + 1, 0.0));
    for (int R = 1; R <= Rmax; ++R)
      for (int T = 1; T <= Rmax; ++T) { g1[R][T] = b1[R][T] / bcount; gq[R][T] = bq[R][T] / bcount; }
    jk1.add_block(g1); jkq.add_block(gq);
    for (auto& r : b1) std::fill(r.begin(), r.end(), 0.0);
    for (auto& r : bq) std::fill(r.begin(), r.end(), 0.0);
    bcount = 0;
  }

  // Record ONE measurement from the chain's current config (after the caller has
  // advanced it). `t` is the time-series row index. Verbatim port of the inner
  // measurement body from the original per-node loop.
  void record(u1::U1HMC<kDim>& hmc, int t) {
    const Real A   = u1::plaq_energy_sum<kDim>(hmc.th, hmc.lat);
    const Real B   = u1::hop_energy_sum<kDim>(hmc.phi, hmc.th, hmc.lat, q);
    const Real pl  = u1::avg_plaquette<kDim>(hmc.th, hmc.lat);
    const Real lp  = u1::higgs_length<kDim>(hmc.phi, hmc.lat);
    const Real le  = u1::link_energy<kDim>(hmc.phi, hmc.th, hmc.lat, q);
    const Real rho = u1::monopole_density<kDim>(hmc.th, hmc.lat);
    std::fprintf(tf, "%d %.15g %.15g %.15g %.15g %.15g %.15g\n", t, A, B, pl, lp, le, rho);
    plaq.add(pl); Lphi.add(lp); Llink.add(le); rhoM.add(rho);

    absP1.add(polyakov_abs<kDim>(hmc.th, hmc.lat, 1));
    wilson_grids<kDim>(hmc.th, hmc.lat, q, Rmax, w1, wq);
    for (int R = 1; R <= Rmax; ++R)
      for (int T = 1; T <= Rmax; ++T) { b1[R][T] += w1[R][T]; bq[R][T] += wq[R][T]; }
    if (++bcount == n_block) flush_probe();

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

  // Close ts, extract photon observables (+cor file), Z_q probe summary, emit the
  // summary row to `sf`, and print the per-node stderr progress line. `acc_hmc`
  // is the measurement-phase HMC acceptance for this chain/rung. Verbatim port.
  void finalize(FILE* sf, const std::string& outdir, double acc_hmc,
                int node, int total) {
    if (tf) std::fclose(tf);
    flush_probe();  // commit any partial trailing block to the Creutz jackknife

    Real m_gamma = 0.0, m_gamma_err = 0.0;
    Real m_gamma_cosh = 0.0, m_gamma_cosh_err = 0.0;
    const int Ltps = photon_samples.empty() ? 0 : static_cast<int>(photon_samples[0].size());
    const int Kc   = static_cast<int>(photon_samples.size());
    if (Ltps >= 2 && Kc >= 1) {
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

      if (Ltps >= 3) {
        const int tmin = std::max(1, Ltps / 4);
        const int tmax = Ltps / 2;
        const gh::PlateauFit pf = gh::plateau(photon_samples, tmin, tmax);
        m_gamma_cosh = pf.mass; m_gamma_cosh_err = pf.err;
      }

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

    Real sigma1 = 0.0, sigma1_err = 0.0, sigmaq = 0.0, sigmaq_err = 0.0;
    const Real n_sigma = 2.0;
    const WGrid Wm1 = jk1.mean_grid(), We1 = jk1.err_grid();
    const WGrid Wmq = jkq.mean_grid(), Weq = jkq.err_grid();
    auto reliableR = [&](const CreutzJackknife& jk, const WGrid& Wm, const WGrid& We,
                         std::vector<int>& out) {
      for (int R = 2; R <= Rmax; ++R) {
        if (!creutz_excl_reason(Wm, We, R, R, n_sigma).empty()) continue;
        const JackResult c = chi_diag_jack(jk, R);
        if (!c.ok || !(c.value > 0.0)) continue;
        if (!(c.value >= n_sigma * c.error)) continue;
        out.push_back(R);
      }
    };
    std::vector<int> R1, Rq;
    reliableR(jk1, Wm1, We1, R1);
    reliableR(jkq, Wmq, Weq, Rq);
    const JackResult s1 = plateau_sigma_jack(jk1, R1);
    const JackResult sq = plateau_sigma_jack(jkq, Rq);
    if (s1.ok) { sigma1 = s1.value; sigma1_err = s1.error; }
    if (sq.ok) { sigmaq = sq.value; sigmaq_err = sq.error; }
    const Real absP1_mean = absP1.mean();
    const Real absP1_err  = absP1.binned_error();
    const Real chiP1      = absP1.susceptibility(Vspatial);

    const double expm = nrows ? sExp / nrows : 0.0;
    std::fprintf(sf, "%.6f,%.6f,%d,%d,%.8g,%.3g,%.8g,%.3g,%.8g,%.3g,%.4f,%.6g,%.8g,%.3g,%.8g,%.3g,%.8g,%.3g,%.8g,%.3g,%.8g,%.3g,%.8g,%.3g,%.8g\n",
                 beta, kappa, L, q,
                 plaq.mean(), plaq.binned_error(),
                 Llink.mean(), Llink.binned_error(),
                 Lphi.mean(), Lphi.binned_error(),
                 acc_hmc, expm,
                 rhoM.mean(), rhoM.binned_error(),
                 m_gamma, m_gamma_err, m_gamma_cosh, m_gamma_cosh_err,
                 sigma1, sigma1_err, sigmaq, sigmaq_err,
                 absP1_mean, absP1_err, chiP1);
    std::fflush(sf);

    std::fprintf(stderr,
      "  node %d/%d  beta=%.4f kappa=%.4f  acc=%.3f <exp(-dH)>=%.4f  plaq=%.4f Llink=%.4f Lphi=%.4f  -> %s\n",
      node + 1, total, beta, kappa, acc_hmc, expm,
      plaq.mean(), Llink.mean(), Lphi.mean(), fname);
  }
};

}  // namespace

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

  // Z_q-confinement probe geometry (ported from src/u1_zqprobe.cpp): charge-1/charge-q
  // Wilson grids out to Rmax = min(L/2,4); ~20 blocks for the correlated Creutz jackknife.
  int Rmax = std::min(L / 2, 4); if (Rmax < 2) Rmax = 2;
  const int target_blocks = 20;
  const int n_block = std::max(1, nmeas / target_blocks);
  // Spatial volume = transverse sites of the charge-1 Polyakov (one per x[tdir]=0 site).
  const Real Vspatial = static_cast<Real>(probe.vol) / static_cast<Real>(L);

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
  // APPENDED (Z_q-confinement probe; backward-compatible -- legacy columns above are
  // untouched): sigma1/sigma1_err = charge-1 string tension (chi(R,R) plateau over the
  // noise-reliable R) -- the PRIMARY order parameter that EXPOSES the deep-Higgs
  // Z_q-confined phase that rho_M is BLIND to (area law sigma1>0 when Z_q-confined for
  // q>=2). sigmaq/sigmaq_err = charge-q control (screened by the condensate -> ~0 deep
  // Higgs; q|q neutral). absP1/absP1_err = <|charge-1 Polyakov|> + its binned error;
  // chiP1 = Vspatial * Var(|P_1|) (cheap finite-V-biased cross-check on the charge-1
  // source). Diagonal chi(R,R) is orientation-symmetric so it is reliable even though
  // u1::wilson_loop symmetrizes over all mu<nu planes (incl spatial-spatial).
  std::fprintf(sf, "beta,kappa,L,q,plaq,plaq_err,Llink,Llink_err,phi2,phi2_err,acceptance,exp_mdH,rho_M,rho_M_err,m_gamma,m_gamma_err,m_gamma_cosh,m_gamma_cosh_err,sigma1,sigma1_err,sigmaq,sigmaq_err,absP1,absP1_err,chiP1\n");

  std::fprintf(stderr, "# u1_scan: D=%d L=%d^%d q=%d lambda=%.4g grid=%dx%d (beta in [%.4g,%.4g], kappa in [%.4g,%.4g])\n",
               kDim, L, kDim, q, lambda, nb, nk, bmin, bmax, kmin, kmax);
  std::fprintf(stderr, "# ntherm=%d nmeas=%d nmd=%d(start) tau=%.4g base_seed=%llu measure_every=%d autotune=%d n_scalar=%d outdir=%s\n",
               ntherm, nmeas, nmd, tau, (unsigned long long)base_seed, measure_every, autotune, n_scalar, outdir.c_str());

  // ---------------------------------------------------------------------------
  // TEMPERED-KAPPA MODE. U1_TEMPER set AND nk>=2: for EACH beta, run ONE
  // U1ReplicaTempering<kDim> kappa-ladder (the kappa grid as the ladder, shared
  // beta/lambda/q/tau/nmd/n_scalar) so a stuck deep-Higgs replica can tunnel by
  // swapping with neighbors. Every rung emits the SAME per-(beta,kappa) summary
  // rows + ts/cor files via MeasureCtx -- IDENTICAL observables, IDENTICAL CSV.
  // Default (unset) falls through to the byte-identical independent-grid path.
  const bool temper_mode = u1::U1ReplicaTempering<kDim>::env_enabled() && (nk >= 2);
  if (temper_mode) {
    // Kappa ladder = the kappa grid (strictly increasing, required by the module).
    std::vector<Real> kladder(nk);
    for (int ik = 0; ik < nk; ++ik) kladder[ik] = kmin + (kmax - kmin) * ik / (nk - 1);
    bool strictly_inc = true;
    for (int ik = 1; ik < nk; ++ik) if (!(kladder[ik] > kladder[ik - 1])) strictly_inc = false;
    if (!strictly_inc) {
      std::fprintf(stderr,
        "ERROR: U1_TEMPER kappa-ladder is not strictly increasing (need kmax>kmin, nk>=2)."
        " Got kmin=%.6g kmax=%.6g nk=%d.\n", kmin, kmax, nk);
      std::fclose(sf);
      return 1;
    }
    const bool cold = (std::getenv("U1_COLD") != nullptr);
    std::fprintf(stderr,
      "# U1_TEMPER: tempered-kappa mode. kappa-ladder (M=%d rungs) per beta:", nk);
    for (Real k : kladder) std::fprintf(stderr, " %.4g", k);
    std::fprintf(stderr, "  (init %s)\n", cold ? "COLD" : "HOT");

    for (int ib = 0; ib < nb; ++ib) {
      const Real beta = (nb <= 1) ? bmin : bmin + (bmax - bmin) * ib / (nb - 1);
      // One seed per beta-ladder; mix base_seed with ib so each ladder is decoupled.
      const std::uint64_t seed0 = Rng::key(base_seed, ib);
      u1::U1ReplicaTempering<kDim> pt(ext, u1::TemperAxis::Kappa, kladder, seed0);
      pt.set_beta(beta); pt.set_lambda(lambda); pt.set_q(q); pt.set_tau(tau);
      pt.set_nmd(nmd);   pt.set_n_scalar(n_scalar);
      pt.enabled = true; pt.n_sweep = measure_every;  // measure_every traj per rung per step
      const std::size_t M = pt.n_replicas();
      // Init each replica: HOT gauge (disordered) by default, COLD (theta=0) under
      // U1_COLD; scalar ordered in both -- mirrors the independent-grid start.
      for (std::size_t k = 0; k < M; ++k) {
        if (!cold) pt.replica(k).hot(0.8);
        pt.replica(k).cold_phi(0.5);
      }
      // Thermalize the whole ladder (swaps active during thermalization).
      for (int t = 0; t < ntherm; ++t) pt.step();
      // Reset HMC + swap counters so the reported acceptances are measurement-phase.
      for (std::size_t k = 0; k < M; ++k) { pt.replica(k).traj_count = 0; pt.replica(k).accept_count = 0; }
      for (std::size_t p = 0; p < pt.n_pairs(); ++p) { pt.swap_attempts[p] = 0; pt.swap_accepts[p] = 0; }

      // One MeasureCtx per rung; advance the ladder together (pt.step() runs all
      // replicas + swaps), recording each rung's current config every step.
      std::vector<std::unique_ptr<MeasureCtx>> ctx;
      ctx.reserve(M);
      for (std::size_t k = 0; k < M; ++k) {
        const Real kappa = pt.coupling(k);
        ctx.emplace_back(new MeasureCtx(ext, outdir, beta, kappa, lambda, L, q,
                                        Rmax, n_block, Vspatial, nmeas));
        if (!ctx.back()->ok()) {
          std::fprintf(stderr, "ERROR: cannot open ts file for beta=%.6f kappa=%.6f\n", beta, kappa);
          std::fclose(sf); return 1;
        }
      }
      for (int t = 0; t < nmeas; ++t) {
        pt.step();  // advances every replica measure_every trajectories + attempts swaps
        for (std::size_t k = 0; k < M; ++k) ctx[k]->record(pt.replica(k), t);
      }
      // Per-rung summary rows + cor files (node index = ladder-flattened grid node).
      for (std::size_t k = 0; k < M; ++k) {
        const int node = ib * nk + static_cast<int>(k);
        ctx[k]->finalize(sf, outdir, pt.hmc_acceptance(k), node, nb * nk);
      }

      // Swap diagnostics: per-pair acceptance + 'ladder connected?' line. WARN on
      // dead rungs (~0 acceptance => the ladder is too coarse at that pair, so the
      // tempering cannot transport configs across it).
      std::fprintf(stderr, "# beta=%.4f swap acceptances (pairs):", beta);
      bool connected = true;
      const double dead_thresh = 0.02;
      for (std::size_t p = 0; p < pt.n_pairs(); ++p) {
        const double a = pt.pair_acceptance(p);
        std::fprintf(stderr, " (%.3g,%.3g):%.3f", pt.coupling(p), pt.coupling(p + 1), a);
        if (a < dead_thresh) connected = false;
      }
      std::fprintf(stderr, "  -> ladder %s\n",
                   connected ? "CONNECTED" : "BROKEN (a pair has ~0 acceptance)");
      if (!connected) {
        for (std::size_t p = 0; p < pt.n_pairs(); ++p) {
          if (pt.pair_acceptance(p) < dead_thresh)
            std::fprintf(stderr,
              "  WARN: DEAD RUNG between kappa=%.4g and kappa=%.4g (swap acc=%.3f < %.2f):"
              " ladder too coarse here -- add intermediate kappa points.\n",
              pt.coupling(p), pt.coupling(p + 1), pt.pair_acceptance(p), dead_thresh);
        }
      }
    }
    std::fclose(sf);
    std::fprintf(stderr, "# done (tempered-kappa): %d nodes (%d beta x %d kappa-rungs) -> %s\n",
                 nb * nk, nb, nk, sumpath.c_str());
    return 0;
  }

  // ---------------------------------------------------------------------------
  // DEFAULT PATH (U1_TEMPER unset, or nk<2): independent single-chain grid nodes.
  // Byte-identical to the pre-tempering driver; measurement refactored into
  // MeasureCtx so the EXACT same code measures both grid nodes and ladder rungs.
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
      // Start: HOT (disordered gauge) by default; U1_COLD env -> COLD (theta=0, the ctor
      // default = ordered gauge) for hot/cold metastability bracketing. Scalar set ordered
      // in both. Mirrors gh_string's GH_COLD.
      if (!std::getenv("U1_COLD")) hmc.hot(0.8);
      hmc.cold_phi(0.5);
      for (int t = 0; t < ntherm; ++t) hmc.trajectory();
      if (autotune) {
        const u1::TuneResult tr = u1::tune_nmd<kDim>(hmc);  // extra thermalization + sets hmc.nmd
        std::fprintf(stderr, "  node %d: autotuned nmd=%d (cal acc=%.3f in_band=%d)\n",
                     node, hmc.nmd, tr.acceptance, (int)tr.in_band);
      }
      hmc.traj_count = 0; hmc.accept_count = 0;

      MeasureCtx ctx(ext, outdir, beta, kappa, lambda, L, q, Rmax, n_block, Vspatial, nmeas);
      if (!ctx.ok()) { std::fprintf(stderr, "ERROR: cannot open %s\n", ctx.fname); std::fclose(sf); return 1; }
      for (int t = 0; t < nmeas; ++t) {
        for (int e = 0; e < measure_every; ++e) hmc.trajectory();
        ctx.record(hmc, t);
      }
      ctx.finalize(sf, outdir, hmc.acceptance(), node, nb * nk);
    }
  }
  std::fclose(sf);
  // Silence unused-variable warnings for the documentation-only extensive counts.
  (void)n_links; (void)n_plaq;
  std::fprintf(stderr, "# done: %d nodes -> %s\n", nb * nk, sumpath.c_str());
  return 0;
}
