// =====================================================================================
// Multicanonical bias in the MATTER hopping B = sum_{x,mu} 2 Re[conj(phi) e^{iq theta} phi'].
//
// WHY: the 2D-LLR windows B and HARD-CONSTRAINS the matter to a B-box; at high B the matter
// cannot be held there (guarded HMC acceptance collapses, res_E2 ~ window), so the deep-Higgs
// first-order region is unreachable (rigorously confirmed: 43/50 high-B cells hit nwiden=3 and
// still res_E2 ~= window). Instead, add a multicanonical weight g(B) to the action that the HMC
// FEELS, so the matter can TUNNEL the symmetric<->Higgs B-jump:
//
//   biased weight  exp(-S + g(B)),  S = beta*A - kappa*B + S_pot.
//   matter force   d/dphi[kappa*B - S_pot + g(B)] = (kappa + g'(B)) dB/dphi - dS_pot/dphi
//                  => use kappa_eff = kappa + g'(B) in scalar_force AND add_matter_force.
//   accept step    biased dH = dH_unbiased - (g(B_f) - g(B_i)).
//   reweight       <O>_canonical = <O exp(-g(B))>_bias / <exp(-g(B))>_bias.
//
// g(B) is built by WANG-LANDAU (g[bin] += f each visit; when the B-histogram is flat to a
// tolerance, f -> f/2; stop at small f) so the biased B-histogram is flat across the jump.
// Piecewise-linear g -> g'(B) is the bin slope (continuous enough for HMC).
// =====================================================================================
#pragma once
#include <vector>
#include <cmath>
#include <cstdio>
#include <algorithm>

namespace gh {
namespace u1 {

struct MucaB {
  double Bmin = 0, Bmax = 1, dB = 1;
  int nbin = 1;
  std::vector<double> g;       // multicanonical weight g(B) at bin centres
  std::vector<long>   H;       // visit histogram (for WL flatness)
  double f = 1.0;              // current WL increment (added to g per visit)

  MucaB() = default;
  MucaB(double lo, double hi, int n) { init(lo, hi, n); }
  void init(double lo, double hi, int n) {
    Bmin = lo; Bmax = hi; nbin = std::max(2, n); dB = (hi - lo) / nbin;
    g.assign(nbin, 0.0); H.assign(nbin, 0); f = 1.0;
  }
  int bin(double B) const {
    int i = (int)std::floor((B - Bmin) / dB);
    return std::min(nbin - 1, std::max(0, i));
  }
  // g(B): piecewise-linear interpolation between bin centres.
  double gval(double B) const {
    const double x = (B - Bmin) / dB - 0.5;     // bin-centre coordinate
    int i = (int)std::floor(x);
    if (i < 0) return g.front() + (g.size() > 1 ? (g[1] - g[0]) : 0.0) * (x - 0);
    if (i >= nbin - 1) return g.back() + (nbin > 1 ? (g[nbin - 1] - g[nbin - 2]) : 0.0) * (x - (nbin - 1));
    const double t = x - i;
    return g[i] * (1 - t) + g[i + 1] * t;
  }
  // g'(B): slope of the piecewise-linear g (per unit B).
  double gprime(double B) const {
    const double x = (B - Bmin) / dB - 0.5;
    int i = (int)std::floor(x);
    i = std::min(nbin - 2, std::max(0, i));
    return (g[i + 1] - g[i]) / dB;
  }
  // Wang-Landau update on the accepted config's B: bump g and histogram.
  void wl_record(double B) { int i = bin(B); g[i] += f; H[i] += 1; }
  // flatness of the histogram over the [lo,hi] B-band; returns min/mean ratio.
  double flatness(double lo, double hi) const {
    long mn = -1; double sum = 0; int cnt = 0;
    for (int i = 0; i < nbin; ++i) {
      const double Bc = Bmin + (i + 0.5) * dB;
      if (Bc < lo || Bc > hi) continue;
      if (mn < 0 || H[i] < mn) mn = H[i];
      sum += H[i]; ++cnt;
    }
    if (!cnt || sum <= 0) return 0.0;
    return (double)mn / (sum / cnt);
  }
  void reset_hist() { std::fill(H.begin(), H.end(), 0); }
  void halve_f() { f *= 0.5; reset_hist(); }

  bool save(const char* path) const {
    FILE* fp = std::fopen(path, "w"); if (!fp) return false;
    std::fprintf(fp, "# MucaB Bmin %g Bmax %g nbin %d f %g\n", Bmin, Bmax, nbin, f);
    for (int i = 0; i < nbin; ++i) std::fprintf(fp, "%g %g %ld\n", Bmin + (i + 0.5) * dB, g[i], H[i]);
    std::fclose(fp); return true;
  }
  bool load(const char* path) {
    FILE* fp = std::fopen(path, "r"); if (!fp) return false;
    char line[512]; double lo, hi, ff; int n;
    if (!std::fgets(line, sizeof line, fp) ||
        std::sscanf(line, "# MucaB Bmin %lf Bmax %lf nbin %d f %lf", &lo, &hi, &n, &ff) != 4) { std::fclose(fp); return false; }
    init(lo, hi, n); f = ff;
    for (int i = 0; i < nbin; ++i) {
      if (!std::fgets(line, sizeof line, fp)) break;
      double bc, gi; long hi2; if (std::sscanf(line, "%lf %lf %ld", &bc, &gi, &hi2) >= 2) { g[i] = gi; H[i] = (i < (int)H.size()) ? hi2 : 0; }
    }
    std::fclose(fp); return true;
  }
};

}  // namespace u1
}  // namespace gh
