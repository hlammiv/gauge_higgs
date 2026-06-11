// =====================================================================================
// Multicanonical-in-matter driver for U(1)+charge-q Higgs. Builds the bias g(B) by Wang-Landau
// so the matter HMC TUNNELS the first-order symmetric<->Higgs B-jump (which plain HMC and the
// 2D-LLR matter-windowing cannot cross), then reweights to the canonical B density of states.
//
// modes:
//   build  L q lambda beta kappa Bmin Bmax nbin [maxsweep nmd tau seed]  -> WL g(B), report tunneling
//   (writes g(B) to <out>.g and the production B-histogram to <out>.hist for the Maxwell analysis)
// =====================================================================================
#include "u1/u1.hpp"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>

using namespace gh;

namespace {
std::array<int, kDim> ext(int L) { std::array<int, kDim> e{}; for (int m = 0; m < kDim; ++m) e[m] = L; return e; }
}

int main(int argc, char** argv) {
  std::setvbuf(stdout, nullptr, _IOLBF, 0);   // line-buffer so progress survives a kill/redirect
  if (argc < 9) {
    std::fprintf(stderr, "usage: %s L q lambda beta kappa Bmin Bmax nbin [maxsweep nmd tau seed outbase]\n", argv[0]);
    return 1;
  }
  auto af = [&](int i, double d){ return i < argc ? std::atof(argv[i]) : d; };
  auto ai = [&](int i, long d){ return i < argc ? std::atol(argv[i]) : d; };
  const int L = (int)ai(1, 4), Q = (int)ai(2, 2);
  const double lam = af(3, 0.5), beta = af(4, 1.0), kappa = af(5, 0.3);
  const double Bmin = af(6, 0), Bmax = af(7, 3000); const int nbin = (int)ai(8, 60);
  const long maxsweep = ai(9, 400000); const int nmd = (int)ai(10, 20); const double tau = af(11, 1.0);
  const std::uint64_t seed = (std::uint64_t)ai(12, 7);
  const char* outbase = argc > 13 ? argv[13] : "mucab";

  u1::U1HMC<kDim> hmc(ext(L), seed);
  hmc.beta = beta; hmc.kappa = kappa; hmc.lambda = lam; hmc.q = Q; hmc.tau = tau; hmc.nmd = nmd;
  hmc.hot(0.8); hmc.cold_phi(0.4);
  for (int t = 0; t < 300; ++t) hmc.trajectory();          // unbiased thermalize (lands in ONE phase)

  u1::MucaB W(Bmin, Bmax, nbin);
  hmc.mucab = &W; hmc.muca_build = true;
  const double band_lo = Bmin + 0.05 * (Bmax - Bmin), band_hi = Bmax - 0.05 * (Bmax - Bmin);
  const double Bmid = 0.5 * (Bmin + Bmax);

  std::printf("# muca-in-B build: L=%d q=%d beta=%g kappa=%g B=[%g,%g] nbin=%d\n", L, Q, beta, kappa, Bmin, Bmax, nbin);
  long sweep = 0; int stage = 0; double Bvis_lo = 1e18, Bvis_hi = -1e18; bool tunneled = false;
  while (sweep < maxsweep && W.f > 0.01) {
    for (int k = 0; k < 2000 && sweep < maxsweep; ++k, ++sweep) {
      hmc.trajectory();
      const double B = hmc.matter_B();
      Bvis_lo = std::min(Bvis_lo, B); Bvis_hi = std::max(Bvis_hi, B);
      if (B < Bmid - 0.1 * (Bmax - Bmin)) tunneled |= (Bvis_hi > Bmid + 0.1 * (Bmax - Bmin));
    }
    const double flat = W.flatness(band_lo, band_hi);
    std::printf("# sweep %ld stage %d f=%.4f flatness=%.2f  B visited [%.0f,%.0f] tunneled=%d\n",
                sweep, stage, W.f, flat, Bvis_lo, Bvis_hi, (int)tunneled);
    std::fflush(stdout);
    if (flat > 0.8) { W.halve_f(); ++stage; }
  }
  char gpath[512]; std::snprintf(gpath, sizeof gpath, "%s.g", outbase); W.save(gpath);
  std::printf("# DONE build: f=%.4f, B visited [%.0f,%.0f] (mid=%.0f). TUNNELED BOTH PHASES: %s\n",
              W.f, Bvis_lo, Bvis_hi, Bmid, tunneled ? "YES" : "NO");

  // ---- production: fixed g, reweight -> canonical B histogram (weights exp(-g(B))) ----
  hmc.muca_build = false; W.reset_hist();
  std::vector<double> Hw(nbin, 0.0);             // sum of exp(-g) per bin (reweighted density)
  const long nprod = 40000;
  for (long t = 0; t < nprod; ++t) {
    hmc.trajectory();
    const double B = hmc.matter_B();
    Hw[W.bin(B)] += std::exp(-W.gval(B));
  }
  char hp[512]; std::snprintf(hp, sizeof hp, "%s.hist", outbase);
  FILE* fp = std::fopen(hp, "w");
  std::fprintf(fp, "# B  reweighted_density  g(B)   (beta=%g kappa=%g)\n", beta, kappa);
  for (int i = 0; i < nbin; ++i) std::fprintf(fp, "%g %g %g\n", Bmin + (i + 0.5) * W.dB, Hw[i], W.g[i]);
  std::fclose(fp);
  std::printf("# wrote %s (g) and %s (reweighted canonical B density)\n", gpath, hp);
  return 0;
}
