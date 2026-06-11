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
  double mA = 0, mB = 0, mP = 0, mL = 0; long n = 0;
  for (long i = 0; i < nsweep; ++i) {
    s.sweep();
    mA += s.A(); mB += s.B(); mP += s.avg_plaq(); mL += s.link_energy(); ++n;
  }
  std::printf("# <A>=%.4f  <B>=%.4f  <plaq>=%.5f  <cos>_link=%.5f\n", mA / n, mB / n, mP / n, mL / n);
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
}  // namespace

int main(int argc, char** argv) {
  std::setvbuf(stdout, nullptr, _IOLBF, 0);
  if (argc < 2) {
    std::fprintf(stderr,
      "usage:\n"
      "  %s point <L> <beta> <kappa> <q> [nsweep ntherm seed n_or]\n"
      "  %s hyst  <L> <beta> <q> <kmin> <kmax> <nk> [nper ntherm0 seed n_or]\n"
      "  %s muca  <L> <beta> <kappa> <q> <Bmin> <Bmax> <nbin> [maxsweep ntherm seed n_or outbase]\n",
      argv[0], argv[0], argv[0]);
    return 1;
  }
  if (!std::strcmp(argv[1], "point")) return mode_point(argc, argv);
  if (!std::strcmp(argv[1], "hyst"))  return mode_hyst(argc, argv);
  if (!std::strcmp(argv[1], "muca"))  return mode_muca(argc, argv);
  std::fprintf(stderr, "unknown mode '%s' (use point|hyst|muca)\n", argv[1]);
  return 1;
}
