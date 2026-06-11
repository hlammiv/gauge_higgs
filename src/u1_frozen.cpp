// Compact U(1) + charge-q FROZEN-LENGTH (|phi|=1) Higgs local-MC driver. Compile-time D = NDIM.
//
// The frozen sampler (heatbath + overrelaxation on the matter phase, Metropolis + exact Z_q
// link-jump on the gauge links) deletes the dynamical-modulus deep-Higgs sampling wall and is the
// apples-to-apples match to Bowler-1981 / Damgaard-Heller-1989. Observables A,B use the SAME
// (beta,kappa)-conjugate convention as u1_scan.cpp / scan_obs.hpp (weight exp(-S), S=beta*A-kappa*B).
//
//   point: ./build/u1_frozen point <L> <beta> <kappa> <q> [nsweep ntherm seed n_or]
//          -> mean A, B, <plaq>, <cos>_link, gauge/Zq acceptance.
//
//   hyst:  ./build/u1_frozen hyst <L> <beta> <q> <kmin> <kmax> <nk> [nper ntherm0 seed n_or]
//          -> kappa-up then kappa-down sweep, carrying the config forward (NO re-thermalize):
//             prints <B>/link for both branches; a hysteresis gap = first-order line (the
//             Bowler/Damgaard-Heller thermal-cycling locator).
#include "u1/u1_frozen.hpp"
#include "u1/monopole.hpp"    // monopole_density<D> (rho_M); defines reduced_plaq_angle<D>
#include "u1/gauge_obs.hpp"   // wilson_grids, polyakov_abs (|P_n|), creutz_chi (sigma_1/sigma_q) -- DEFINITIVE discriminants
#define GH_U1_HAVE_REDUCED_PLAQ_ANGLE   // monopole.hpp already provided it -> photon_mass.hpp must not redefine
#include "u1/photon_structure.hpp"      // m_gamma via the static magnetic structure factor (needs L_s>=16)
#include "u1/zq_gauge.hpp"              // pure Z_q gauge theory = the kappa->inf MATCHING target (digitization proof)
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <string>

using namespace gh;
using namespace gh::u1;

namespace {
std::array<int, kDim> ext(int L) { std::array<int, kDim> e{}; for (int m = 0; m < kDim; ++m) e[m] = L; return e; }

int mode_point(int argc, char** argv) {
  auto af = [&](int i, double d){ return i < argc ? std::atof(argv[i]) : d; };
  auto ai = [&](int i, long d){ return i < argc ? std::atol(argv[i]) : d; };
  const int L = (int)ai(2, 4);
  const double beta = af(3, 1.0), kappa = af(4, 0.3);
  const int q = (int)ai(5, 2);
  const long nsweep = ai(6, 4000), ntherm = ai(7, 1000);
  const std::uint64_t seed = (std::uint64_t)ai(8, 7);
  const int n_or = (int)ai(9, 3);

  U1Frozen<kDim> s(ext(L), seed);
  s.beta = beta; s.kappa = kappa; s.q = q; s.n_or = n_or; s.hot(0.8);
  std::printf("# u1_frozen point: D=%d L=%d beta=%g kappa=%g q=%d  nsweep=%ld ntherm=%ld n_or=%d\n",
              kDim, L, beta, kappa, q, nsweep, ntherm, n_or);
  s.thermalize((int)ntherm);
  double mA = 0, mB = 0, mP = 0, mL = 0, mA2 = 0, mB2 = 0; long n = 0;
  // DEFINITIVE gauge-sector discriminants (docs/gauge_phase_observables_DEFINITIVE.md), measured on a stride.
  double mRho = 0, mP1 = 0, mPq = 0; long ng = 0;
  WGrid w1(3, std::vector<double>(3, 0.0)), wq(3, std::vector<double>(3, 0.0));
  WGrid acc1(3, std::vector<double>(3, 0.0)), accq(3, std::vector<double>(3, 0.0));
  for (long i = 0; i < nsweep; ++i) {
    s.sweep();
    const double A = s.A(), B = s.B();
    mA += A; mB += B; mA2 += A * A; mB2 += B * B; mP += s.avg_plaq(); mL += s.link_energy(); ++n;
    if (i % 4 == 0) {                                  // gauge discriminants on a stride (decorrelate + cheaper)
      mRho += monopole_density<kDim>(s.th, s.lat);
      mP1  += polyakov_abs<kDim>(s.th, s.lat, 1);
      mPq  += polyakov_abs<kDim>(s.th, s.lat, q);
      wilson_grids<kDim>(s.th, s.lat, q, 2, w1, wq);
      for (int R = 1; R <= 2; ++R) for (int T = 1; T <= 2; ++T) { acc1[R][T] += w1[R][T]; accq[R][T] += wq[R][T]; }
      ++ng;
    }
  }
  for (int R = 1; R <= 2; ++R) for (int T = 1; T <= 2; ++T) { acc1[R][T] /= ng; accq[R][T] /= ng; }
  const double sig1 = creutz_chi(acc1, 2), sigq = creutz_chi(accq, 2);   // sigma_1 (charge-1) + sigma_q control
  const double V = double(s.lat.vol);
  // susceptibilities = Var/vol (intensive); chi_link (conjugate to kappa) PEAKS at the Higgs transition,
  // chi_plaq (conjugate to beta) peaks at the gauge transition. The order parameter, NOT the hopping energy.
  const double chiB = (mB2 / n - (mB / n) * (mB / n)) / V;
  const double chiA = (mA2 / n - (mA / n) * (mA / n)) / V;
  std::printf("# <A>=%.4f  <B>=%.4f  <plaq>=%.5f  <cos>_link=%.5f\n", mA / n, mB / n, mP / n, mL / n);
  std::printf("# chi_plaq=%.5f  chi_link=%.5f   (Var/vol; peak => transition)\n", chiA, chiB);
  // IDENTIFY set: sigma_1>0 = charge-1 confined (Confined OR deep-Higgs-Z_q); sigma_q=control; |P1| = deconf order param.
  std::printf("# rho_M=%.5f  |P1|=%.5f  |Pq|=%.5f  sigma1=%.5f  sigmaq=%.5f   (sigma=chi(2,2); area>0=confine)\n",
              mRho / ng, mP1 / ng, mPq / ng, sig1, sigq);
  std::printf("# gauge_acc=%.2f gor_acc=%.2f zq_flip=%.2f matter_acc=%.2f step=%.3f\n",
              s.gauge_acc(), s.gor_rate(), s.zq_flip(), s.matter_acc(), s.gauge_step);
  return 0;
}

// Build the multicanonical weight g(B) by Wang-Landau so the frozen matter TUNNELS the deep-Higgs
// first-order B-jump (which plain heatbath cannot cross). Bounded domain B in [Bmin,Bmax]; bias enters
// the ACCEPT step of every B-changing kernel. Writes g(B) to <outbase>.g.
int mode_muca(int argc, char** argv) {
  auto af = [&](int i, double d){ return i < argc ? std::atof(argv[i]) : d; };
  auto ai = [&](int i, long d){ return i < argc ? std::atol(argv[i]) : d; };
  const int L = (int)ai(2, 4);
  const double beta = af(3, 1.0), kappa = af(4, 0.3);
  const int q = (int)ai(5, 2);
  const double Bmin = af(6, 0.0), Bmax = af(7, 2048.0); const int nbin = (int)ai(8, 64);
  const long maxsweep = ai(9, 400000), ntherm = ai(10, 2000);
  const std::uint64_t seed = (std::uint64_t)ai(11, 7);
  const int n_or = (int)ai(12, 3);
  const char* outbase = argc > 13 ? argv[13] : "u1frozen_muca";

  U1Frozen<kDim> s(ext(L), seed);
  s.beta = beta; s.kappa = kappa; s.q = q; s.n_or = n_or; s.hot(0.8);
  s.thermalize((int)ntherm);                       // unbiased thermalize (lands in ONE phase)

  // Convention: biased weight exp(-S + g(B)); accept uses +Δg (muca_accept); reweight uses exp(-g).
  // Wang-Landau converges g -> -ln rho_canonical by DECREMENTING g on each visit (this SUPPRESSES the
  // already-occupied phase; an INCREMENT would be positive feedback and pin the system in one phase).
  MucaB W(Bmin, Bmax, nbin);
  s.mucab = &W; s.muca_build = false; s.matter_metro = true;   // record WL here (correct sign); broad proposals
  const double lo = Bmin + 0.05 * (Bmax - Bmin), hi = Bmax - 0.05 * (Bmax - Bmin);
  const double Bmid = 0.5 * (Bmin + Bmax);
  std::printf("# u1_frozen muca: D=%d L=%d beta=%g kappa=%g q=%d  B=[%g,%g] nbin=%d\n",
              kDim, L, beta, kappa, q, Bmin, Bmax, nbin);
  long sweep = 0; int stage = 0; double bvlo = 1e18, bvhi = -1e18; bool tunneled = false;
  while (sweep < maxsweep && W.f > 0.01) {
    for (int k = 0; k < 2000 && sweep < maxsweep; ++k, ++sweep) {
      s.sweep();
      const double B = s.cur_B();
      const int bi = W.bin(B); W.g[bi] -= W.f; W.H[bi] += 1;     // WL: suppress the visited B-bin
      bvlo = std::min(bvlo, B); bvhi = std::max(bvhi, B);
      if (B < Bmid - 0.1 * (Bmax - Bmin)) tunneled |= (bvhi > Bmid + 0.1 * (Bmax - Bmin));
    }
    const double flat = W.flatness(lo, hi);
    std::printf("# sweep %ld stage %d f=%.4f flat=%.2f  B[%.0f,%.0f] tunneled=%d\n",
                sweep, stage, W.f, flat, bvlo, bvhi, (int)tunneled);
    if (flat > 0.8) { W.halve_f(); ++stage; }
  }
  char gpath[512]; std::snprintf(gpath, sizeof gpath, "%s.g", outbase); W.save(gpath);
  std::printf("# DONE: f=%.4f B[%.0f,%.0f] (mid=%.0f) TUNNELED=%s -> wrote %s\n",
              W.f, bvlo, bvhi, Bmid, tunneled ? "YES" : "NO", gpath);
  return 0;
}

// Build the multicanonical weight g(A) in the GAUGE ACTION A by Wang-Landau, so the gauge sector TUNNELS the
// large-beta Z_q (de)confinement ordering barrier (a single-link Z_q flip cannot cross it). Same machinery as
// mode_muca but on A; the matter stays unbiased (von Mises heatbath). Writes g(A) to <outbase>.gA.
int mode_mucaA(int argc, char** argv) {
  auto af = [&](int i, double d){ return i < argc ? std::atof(argv[i]) : d; };
  auto ai = [&](int i, long d){ return i < argc ? std::atol(argv[i]) : d; };
  const int L = (int)ai(2, 4);
  const double beta = af(3, 1.0), kappa = af(4, 0.3);
  const int q = (int)ai(5, 2);
  const double Amin = af(6, 0.0), Amax = af(7, 3072.0); const int nbin = (int)ai(8, 64);
  const long maxsweep = ai(9, 400000), ntherm = ai(10, 2000);
  const std::uint64_t seed = (std::uint64_t)ai(11, 7);
  const int n_or = (int)ai(12, 3);
  const char* outbase = argc > 13 ? argv[13] : "u1frozen_mucaA";

  U1Frozen<kDim> s(ext(L), seed);
  s.beta = beta; s.kappa = kappa; s.q = q; s.n_or = n_or; s.hot(0.8);
  s.thermalize((int)ntherm);                       // unbiased thermalize (fixes the gauge step; lands in ONE phase)

  // Convention: biased weight exp(-S + g(A)); accept uses +Δg; reweight uses exp(-g). WL DECREMENTS g on visit.
  MucaB WA(Amin, Amax, nbin);
  s.mucaA = &WA;
  const double lo = Amin + 0.05 * (Amax - Amin), hi = Amax - 0.05 * (Amax - Amin);
  const double Amid = 0.5 * (Amin + Amax);
  std::printf("# u1_frozen mucaA: D=%d L=%d beta=%g kappa=%g q=%d  A=[%g,%g] nbin=%d  (tunnels the large-beta Z_q barrier)\n",
              kDim, L, beta, kappa, q, Amin, Amax, nbin);
  long sweep = 0; int stage = 0; double avlo = 1e18, avhi = -1e18; bool tunneled = false;
  while (sweep < maxsweep && WA.f > 0.01) {
    for (int k = 0; k < 2000 && sweep < maxsweep; ++k, ++sweep) {
      s.sweep();
      const double A = s.cur_A();
      const int bi = WA.bin(A); WA.g[bi] -= WA.f; WA.H[bi] += 1;
      avlo = std::min(avlo, A); avhi = std::max(avhi, A);
      if (A < Amid - 0.1 * (Amax - Amin)) tunneled |= (avhi > Amid + 0.1 * (Amax - Amin));
    }
    const double flat = WA.flatness(lo, hi);
    std::printf("# sweep %ld stage %d f=%.4f flat=%.2f  A[%.0f,%.0f] tunneled=%d\n",
                sweep, stage, WA.f, flat, avlo, avhi, (int)tunneled);
    if (flat > 0.8) { WA.halve_f(); ++stage; }
  }
  char gpath[512]; std::snprintf(gpath, sizeof gpath, "%s.gA", outbase); WA.save(gpath);
  std::printf("# DONE: f=%.4f A[%.0f,%.0f] (mid=%.0f) TUNNELED=%s -> wrote %s\n",
              WA.f, avlo, avhi, Amid, tunneled ? "YES" : "NO", gpath);
  return 0;
}

// Measure <B>/link over nper sweeps WITHOUT re-thermalizing (carry the config forward).
double branch_point(U1Frozen<kDim>& s, double kappa, long nper) {
  s.kappa = kappa;
  double mL = 0; long n = 0;
  for (long i = 0; i < nper; ++i) { s.sweep(); mL += s.link_energy(); ++n; }
  return mL / n;
}

int mode_hyst(int argc, char** argv) {
  auto af = [&](int i, double d){ return i < argc ? std::atof(argv[i]) : d; };
  auto ai = [&](int i, long d){ return i < argc ? std::atol(argv[i]) : d; };
  const int L = (int)ai(2, 6);
  const double beta = af(3, 1.0);
  const int q = (int)ai(4, 2);
  const double kmin = af(5, 0.1), kmax = af(6, 0.6);
  const int nk = (int)ai(7, 11);
  const long nper = ai(8, 800), ntherm0 = ai(9, 1500);
  const std::uint64_t seed = (std::uint64_t)ai(10, 7);
  const int n_or = (int)ai(11, 3);

  U1Frozen<kDim> s(ext(L), seed);
  s.beta = beta; s.q = q; s.n_or = n_or;
  std::printf("# u1_frozen hyst: D=%d L=%d beta=%g q=%d  kappa[%g,%g]x%d  nper=%ld ntherm0=%ld n_or=%d\n",
              kDim, L, beta, q, kmin, kmax, nk, nper, ntherm0, n_or);
  std::printf("# col: kappa   <cos>_up   <cos>_down   (hysteresis gap => first-order)\n");

  const double dk = nk > 1 ? (kmax - kmin) / (nk - 1) : 0.0;
  std::vector<double> up(nk), dn(nk);
  // cold start at kmin (deep-confined branch coming UP)
  s.cold(); s.kappa = kmin; s.thermalize((int)ntherm0);
  for (int i = 0; i < nk; ++i) up[i] = branch_point(s, kmin + i * dk, nper);
  // continue DOWN from the high-kappa (Higgs) end, carrying the config
  for (int i = nk - 1; i >= 0; --i) dn[i] = branch_point(s, kmin + i * dk, nper);
  for (int i = 0; i < nk; ++i)
    std::printf("%.4f  %.5f  %.5f  %+.5f\n", kmin + i * dk, up[i], dn[i], dn[i] - up[i]);
  return 0;
}
// Photon mass m_gamma from the static magnetic structure factor on an ANISOTROPIC Ls^3 x Lt lattice
// (geometric anisotropy: large spatial L_s, modest L_t, same beta). m_gamma=0 ONLY in the Coulomb phase, so
// this is the observable that resolves the q>=5 intermediate-Coulomb WEDGE the gauge string tension cannot see.
int mode_pm(int argc, char** argv) {
  auto af = [&](int i, double d){ return i < argc ? std::atof(argv[i]) : d; };
  auto ai = [&](int i, long d){ return i < argc ? std::atol(argv[i]) : d; };
  const int Ls = (int)ai(2, 16), Lt = (int)ai(3, 8);
  const double beta = af(4, 1.0), kappa = af(5, 0.0);
  const int q = (int)ai(6, 2);
  const long nsweep = ai(7, 1500), ntherm = ai(8, 800);
  const int meas_every = (int)ai(9, 5);
  const std::uint64_t seed = (std::uint64_t)ai(10, 7);
  const int n_or = (int)ai(11, 3);

  std::array<int, kDim> e{}; e[0] = Lt; for (int m = 1; m < kDim; ++m) e[m] = Ls;   // dir 0=time(Lt), 1..=space(Ls)
  U1Frozen<kDim> s(e, seed);
  s.beta = beta; s.kappa = kappa; s.q = q; s.n_or = n_or; s.hot(0.8);
  s.thermalize((int)ntherm);
  const u1::PhotonMomenta<kDim> mom = u1::photon_momenta<kDim>(s.lat);
  const int ng = mom.n_groups();
  std::vector<std::vector<Real>> perConfig;
  for (long i = 0; i < nsweep; ++i) {
    s.sweep();
    if (i % meas_every == 0) perConfig.push_back(u1::photon_structure_factor<kDim>(s.th, s.lat, mom));
  }
  const int nfit = ng >= 4 ? 4 : ng;
  const u1::PhotonMassFit fit = u1::photon_mass_fit<kDim>(perConfig, mom, nfit);
  std::printf("# u1_frozen pm: Ls=%d Lt=%d beta=%g kappa=%g q=%d  nmeas=%zu ngroups=%d nfit=%d\n",
              Ls, Lt, beta, kappa, q, perConfig.size(), ng, nfit);
  std::printf("# m_gamma=%.5f  m2=%.6f +- %.6f  phat2_min=%.4f   (m_gamma~0 => COULOMB)\n",
              fit.m_gamma, fit.m2, fit.m2_err, ng ? mom.phat2[0] : 0.0);
  for (int g = 0; g < ng; ++g)
    std::printf("# phat2=%.4f  R=%.5f +- %.5f\n", mom.phat2[g],
                g < (int)fit.R.size() ? fit.R[g] : 0.0, g < (int)fit.R_err.size() ? fit.R_err[g] : 0.0);
  return 0;
}

// Pure Z_q gauge theory point run -- the kappa->inf matching target. Emits the SAME gauge discriminants as
// mode_point so the deep-Higgs U(1)+charge-q runs can be compared directly (digitization: U(1)+q(beta,kappa->inf)
// -> pure Z_q(beta)).
int mode_zq(int argc, char** argv) {
  auto af = [&](int i, double d){ return i < argc ? std::atof(argv[i]) : d; };
  auto ai = [&](int i, long d){ return i < argc ? std::atol(argv[i]) : d; };
  const int L = (int)ai(2, 8);
  const double beta = af(3, 1.0);
  const int q = (int)ai(4, 2);
  const long nsweep = ai(5, 4000), ntherm = ai(6, 1500);
  const std::uint64_t seed = (std::uint64_t)ai(7, 7);

  ZqGauge<kDim> s(ext(L), seed); s.beta = beta; s.q = q; s.hot();
  s.thermalize((int)ntherm);
  double mP = 0, mA = 0, mRho = 0, mP1 = 0, mPq = 0; long n = 0, ng = 0;
  WGrid w1(3, std::vector<double>(3, 0.0)), wq(3, std::vector<double>(3, 0.0));
  WGrid acc1(3, std::vector<double>(3, 0.0)), accq(3, std::vector<double>(3, 0.0));
  for (long i = 0; i < nsweep; ++i) {
    s.sweep();
    mP += s.avg_plaq(); mA += s.A(); ++n;
    if (i % 4 == 0) {
      mRho += monopole_density<kDim>(s.th, s.lat);
      mP1  += polyakov_abs<kDim>(s.th, s.lat, 1);
      mPq  += polyakov_abs<kDim>(s.th, s.lat, q);
      wilson_grids<kDim>(s.th, s.lat, q, 2, w1, wq);
      for (int R = 1; R <= 2; ++R) for (int T = 1; T <= 2; ++T) { acc1[R][T] += w1[R][T]; accq[R][T] += wq[R][T]; }
      ++ng;
    }
  }
  for (int R = 1; R <= 2; ++R) for (int T = 1; T <= 2; ++T) { acc1[R][T] /= ng; accq[R][T] /= ng; }
  const double sig1 = creutz_chi(acc1, 2), sigq = creutz_chi(accq, 2);
  std::printf("# u1_frozen zq (pure Z_q gauge): D=%d L=%d beta=%g q=%d  nsweep=%ld\n", kDim, L, beta, q, nsweep);
  std::printf("# <A>=%.4f  <plaq>=%.5f\n", mA / n, mP / n);
  std::printf("# rho_M=%.5f  |P1|=%.5f  |Pq|=%.5f  sigma1=%.5f  sigmaq=%.5f   (sigma=chi(2,2); area>0=confine)\n",
              mRho / ng, mP1 / ng, mPq / ng, sig1, sigq);
  return 0;
}

}  // namespace

int main(int argc, char** argv) {
  std::setvbuf(stdout, nullptr, _IOLBF, 0);
  if (argc < 2) {
    std::fprintf(stderr,
      "usage:\n"
      "  %s point <L> <beta> <kappa> <q> [nsweep ntherm seed n_or]\n"
      "  %s hyst  <L> <beta> <q> <kmin> <kmax> <nk> [nper ntherm0 seed n_or]\n"
      "  %s muca  <L> <beta> <kappa> <q> <Bmin> <Bmax> <nbin> [maxsweep ntherm seed n_or outbase]\n"
      "  %s mucaA <L> <beta> <kappa> <q> <Amin> <Amax> <nbin> [maxsweep ntherm seed n_or outbase]\n"
      "  %s pm    <Ls> <Lt> <beta> <kappa> <q> [nsweep ntherm meas_every seed n_or]   (m_gamma; Ls>=16)\n"
      "  %s zq    <L> <beta> <q> [nsweep ntherm seed]   (pure Z_q gauge; kappa->inf matching target)\n",
      argv[0], argv[0], argv[0], argv[0], argv[0], argv[0]);
    return 1;
  }
  if (!std::strcmp(argv[1], "point")) return mode_point(argc, argv);
  if (!std::strcmp(argv[1], "hyst"))  return mode_hyst(argc, argv);
  if (!std::strcmp(argv[1], "muca"))  return mode_muca(argc, argv);
  if (!std::strcmp(argv[1], "mucaA")) return mode_mucaA(argc, argv);
  if (!std::strcmp(argv[1], "pm"))    return mode_pm(argc, argv);
  if (!std::strcmp(argv[1], "zq"))    return mode_zq(argc, argv);
  std::fprintf(stderr, "unknown mode '%s' (use point|hyst|muca|mucaA|pm|zq)\n", argv[1]);
  return 1;
}
