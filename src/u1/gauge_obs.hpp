#pragma once
// Gauge-sector PHASE DISCRIMINANTS for U(1)+charge-q (functions of the link field theta only), so the
// frozen driver can emit the DEFINITIVE 3-phase fingerprint (docs/gauge_phase_observables_DEFINITIVE.md):
//   sigma_1 (charge-1 string tension; area>0 in Confined AND deep-Higgs-Z_q, the line {m_gamma,rho_M} miss),
//   sigma_q (screened charge-q CONTROL -> 0 while sigma_1>0  =>  residue is Z_q not trivial),
//   |P_1| (charge-1 Polyakov modulus; ~0 in Confined+deep-Z_q, >0 in Coulomb+deconfined-Higgs), |P_q| control.
// (rho_M from monopole.hpp; m_gamma from photon_structure.hpp needs L_s>=16, run separately.)
// wilson_grids / polyakov_abs ported verbatim from src/u1_scan.cpp (validated on the HMC path).
#include "core/geometry.hpp"
#include "core/config.hpp"
#include <vector>
#include <array>
#include <cmath>

namespace gh {
namespace u1 {

using WGrid = std::vector<std::vector<Real>>;

// charge-1 and charge-q rectangular Wilson loops g1[R][T]=<cos(loop)>, gq[R][T]=<cos(q*loop)>,
// 1<=R,T<=Rmax, symmetrized over all mu<nu planes and all sites (rows/col 0 unused).
template <int D>
void wilson_grids(const std::vector<Real>& th, const Lattice<D>& lat, int q, int Rmax, WGrid& g1, WGrid& gq) {
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
  const Real per = static_cast<Real>(count) / static_cast<Real>(Rmax * Rmax);
  if (per > 0.0)
    for (int R = 1; R <= Rmax; ++R)
      for (int T = 1; T <= Rmax; ++T) { g1[R][T] /= per; gq[R][T] /= per; }
}

// |P_n| = |(1/Vsp) sum_x exp(i n line_t(x))| for ONE config (charge-n Polyakov modulus). NOT polyakov(...,m)
// which returns <cos(line)> and averages to ~0 -- this magnitude is the deconfinement order parameter.
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

// Creutz ratio chi(R,R) = -ln[ W(R,R) W(R-1,R-1) / (W(R,R-1) W(R-1,R)) ]  -> string tension at scale R.
// Area law (confinement) => chi>0 and ~R-independent; perimeter (deconfined) => chi->0.
inline Real creutz_chi(const WGrid& W, int R) {
  if (R < 2 || (int)W.size() <= R) return std::nan("");
  const Real num = W[R][R] * W[R - 1][R - 1], den = W[R][R - 1] * W[R - 1][R];
  if (num <= 0.0 || den <= 0.0) return std::nan("");
  return -std::log(num / den);
}

}  // namespace u1
}  // namespace gh
