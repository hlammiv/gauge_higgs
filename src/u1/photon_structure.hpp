#pragma once
// Compact-U(1) photon mass from the STATIC (zero-frequency) magnetic structure
// factor at LOW NON-ZERO spatial momentum, fitting the transverse photon pole.
//
// WHY NOT photon_mass.hpp.  The zero-spatial-momentum, position-space transverse
// field-strength correlator C(dt) (u1/photon_mass.hpp) has three fatal problems:
//   (1) a flux-conservation sum rule sum_dt C(dt)=0 forces a negative tail, so
//       there is no clean exponential pole;
//   (2) the temporal plaquette F_{0i} straddles two time slices, smearing t;
//   (3) the p_t=pi temporal doubler contaminates the zero-frequency projection.
// The cure is to stay at p_0=0 by summing over ALL time slices (a static
// structure factor) and to probe LOW NON-ZERO spatial momentum, where the
// transverse photon propagator D_T(p) = Z/(phat2 + m^2) is directly visible.
//
// FIELD STRENGTH.  F_{ij}(x) is the reduced (principal-value) plaquette angle in
// a purely SPATIAL plane (i,j != time dir 0); it is gauge invariant (the
// plaquette angle is gauge invariant, so is its principal value). We reuse
// gh::u1::reduce_angle and the plaquette-angle accessor plaq_angle<D>.
//
// STATIC MAGNETIC STRUCTURE FACTOR.  For a spatial mode vector p (with p_0 = 0),
//   Ftilde_{ij}(p) = sum_{all sites x} exp(-i (p_1 x_1 + ... + p_{D-1} x_{D-1})) F_{ij}(x).
// Summing over x INCLUDES the time coordinate x_0, which is exactly what enforces
// p_0 = 0 (the static / zero-frequency projection). Per config and per momentum:
//   S(p) = sum_{spatial i<j} |Ftilde_{ij}(p)|^2.
// A direct O(V * n_momenta) DFT is used (no FFT needed at these sizes).
//
// LATTICE MOMENTUM.  For an integer spatial mode n = (n_1,...,n_{D-1}),
//   p_k = 2 pi n_k / L_k ,  phat_k = 2 sin(p_k / 2) ,  phat2 = sum_k phat_k^2.
// We enumerate the small phat2 values strictly below the doubler: each n_k in
// {0,1,2}, excluding n=0. The lowest distinct phat2 (for L large, equal extents)
// are (1,0,0) < (1,1,0) < (1,1,1) < (2,0,0) < ... We GROUP momenta with equal
// phat2 (cubic-symmetry copies) and AVERAGE S(p) over the group for statistics.
//
// POLE MODEL.  The transverse photon propagator D_T(p) = Z/(phat2 + m^2). The
// magnetic structure factor of the field strength measures the CURL of the gauge
// field, so <S(p)> is proportional to phat2 * D_T(p) = phat2 * Z/(phat2 + m^2)
// (one power of momentum from each of the two field-strength factors transverse
// to p, with the longitudinal part dropping out). Define
//   R(p) = <S(p)> / phat2  ~  Z/(phat2 + m^2),
// so the INVERSE is linear in phat2:
//   R(p)^{-1} = (phat2 + m^2)/Z = a + b*phat2 ,  a = m^2/Z ,  b = 1/Z .
// A linear least-squares fit of R^{-1} vs phat2 over the lowest few distinct
// phat2 (below the doubler) gives m^2 = a/b = intercept/slope and
// m_gamma = sqrt(max(0, m^2)). The overall V/normalization and the proportionality
// prefactor enter ONLY through Z (slope and intercept scale together), so m^2 is
// independent of them -- only the ratio a/b matters.
//
// ERRORS.  Per-config S(p) over the momentum set is stored so the fit can be
// redone under a delete-1 jackknife over configurations, giving m^2 +/- err.
#include "core/config.hpp"
#include "core/geometry.hpp"
#include "u1/u1.hpp"
// reduce_angle and reduced_plaq_angle<D> live in u1/photon_mass.hpp. A driver may
// already have pulled in monopole.hpp's identical reduced_plaq_angle<D> and set
// GH_U1_HAVE_REDUCED_PLAQ_ANGLE; including photon_mass.hpp here respects that guard.
#include "u1/photon_mass.hpp"
#include <vector>
#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <algorithm>
#include <limits>
#include <utility>

namespace gh {
namespace u1 {

// ---- momentum set ---------------------------------------------------------
// One group of cubic-symmetry-equivalent spatial momenta sharing a phat2 value.
// `modes[g]` is the list of integer spatial mode vectors (length D-1) in group g.
template <int D>
struct PhotonMomenta {
  std::vector<Real>                          phat2;  // distinct phat2, ascending
  std::vector<std::vector<std::array<int, D - 1>>> modes;  // modes[g] -> the n's
  int n_groups() const { return static_cast<int>(phat2.size()); }
  int n_momenta() const {
    int n = 0; for (const auto& g : modes) n += static_cast<int>(g.size()); return n;
  }
};

// Build the low-momentum set: spatial mode components n_k in {0,1,2}, excluding
// n = 0, keeping phat2 strictly at or below the (2,0,0) doubler-edge value so the
// fit window sits well below the lattice doubler. Momenta are grouped by (binned)
// equal phat2 and the groups are returned in ascending phat2 order.
template <int D>
PhotonMomenta<D> photon_momenta(const Lattice<D>& lat) {
  static_assert(D >= 2, "need at least one spatial direction");
  constexpr int Ds = D - 1;                          // spatial directions 1..D-1
  // phat for a single component n_k (direction d = k+1, since dir 0 is time).
  auto phat = [&](int k, int n) -> Real {
    const Real p = 2.0 * kPi * static_cast<Real>(n) / static_cast<Real>(lat.L[k + 1]);
    return 2.0 * std::sin(0.5 * p);
  };
  // Doubler-edge cap: the (2,0,0)-type value, i.e. phat^2 of a single n_k=2 mode.
  // Using direction 1's extent as the reference (extents are equal in validation).
  const Real cap = phat(0, 2) * phat(0, 2) + 1e-9;

  // Enumerate all n in {0,1,2}^Ds \ {0}, compute phat2, keep phat2 <= cap.
  struct Item { Real phat2; std::array<int, Ds> n; };
  std::vector<Item> items;
  std::array<int, Ds> n{};
  const long total = [&]{ long t = 1; for (int k = 0; k < Ds; ++k) t *= 3; return t; }();
  for (long idx = 0; idx < total; ++idx) {
    long r = idx;
    bool nonzero = false;
    for (int k = 0; k < Ds; ++k) { n[k] = static_cast<int>(r % 3); r /= 3; if (n[k]) nonzero = true; }
    if (!nonzero) continue;
    Real p2 = 0.0;
    for (int k = 0; k < Ds; ++k) { const Real ph = phat(k, n[k]); p2 += ph * ph; }
    if (p2 <= cap) items.push_back({p2, n});
  }
  std::sort(items.begin(), items.end(),
            [](const Item& a, const Item& b) { return a.phat2 < b.phat2; });

  PhotonMomenta<D> out;
  const Real tol = 1e-9;
  for (const auto& it : items) {
    if (out.phat2.empty() || std::fabs(it.phat2 - out.phat2.back()) > tol) {
      out.phat2.push_back(it.phat2);
      out.modes.emplace_back();
    }
    std::array<int, D - 1> nv{};
    for (int k = 0; k < Ds; ++k) nv[k] = it.n[k];
    out.modes.back().push_back(nv);
  }
  return out;
}

// ---- per-config structure factor ------------------------------------------
// Compute S(p), averaged over the cubic-symmetry copies in each phat2 group, for
// one configuration's link angles `th`. Output length == mom.n_groups(); entry g
// is < |Ftilde_{ij}(p)|^2 >_{spatial i<j, momenta in group g}.
//
// Ftilde_{ij}(p) = sum_x exp(-i sum_k p_k x_{k+1}) F_{ij}(x), summed over ALL sites
// x (so x_0 is summed -> p_0 = 0). p_k = 2 pi n_k / L_{k+1}. Direct DFT.
template <int D>
std::vector<Real> photon_structure_factor(const std::vector<Real>& th,
                                          const Lattice<D>& lat,
                                          const PhotonMomenta<D>& mom) {
  constexpr int Ds = D - 1;
  const int ng = mom.n_groups();
  std::vector<Real> S(ng, 0.0);
  if (ng == 0) return S;

  // Spatial planes (i,j), 1 <= i < j <= D-1.
  std::vector<std::pair<int, int>> planes;
  for (int i = 1; i < D; ++i)
    for (int j = i + 1; j < D; ++j) planes.emplace_back(i, j);
  const int npl = static_cast<int>(planes.size());

  // Flatten the momentum list so the (parallelizable) site loop accumulates each
  // (plane, momentum) Fourier sum once. flatM[m] = mode vector; flatG[m] = group.
  std::vector<std::array<int, Ds>> flatM;
  std::vector<int> flatG;
  for (int g = 0; g < ng; ++g)
    for (const auto& nv : mom.modes[g]) {
      std::array<int, Ds> a{};
      for (int k = 0; k < Ds; ++k) a[k] = nv[k];
      flatM.push_back(a); flatG.push_back(g);
    }
  const int nm = static_cast<int>(flatM.size());

  // Precompute per-direction angular frequency w_k = 2 pi / L_{k+1}; the phase of
  // site x for mode m is sum_k flatM[m][k] * w_k * x_{k+1}.
  std::array<Real, Ds> w{};
  for (int k = 0; k < Ds; ++k) w[k] = 2.0 * kPi / static_cast<Real>(lat.L[k + 1]);

  // Accumulators: re/im of Ftilde[plane][momentum]. Use a flat [npl*nm] buffer.
  const int nacc = npl * nm;
  std::vector<Real> re(nacc, 0.0), im(nacc, 0.0);

  #pragma omp parallel
  {
    std::vector<Real> lre(nacc, 0.0), lim(nacc, 0.0);
    std::vector<Real> Fpl(npl, 0.0);   // per-site reduced field strength, by plane
    std::array<int, D> xc{};
    #pragma omp for schedule(static) nowait
    for (std::int64_t s = 0; s < lat.vol; ++s) {
      lat.coords(s, xc);
      // Reduced spatial field strength for every plane at this site.
      // (Cache to avoid recomputing per momentum.)
      for (int pidx = 0; pidx < npl; ++pidx)
        Fpl[pidx] = reduced_plaq_angle<D>(th, lat, s, planes[pidx].first, planes[pidx].second);
      for (int m = 0; m < nm; ++m) {
        Real phase = 0.0;
        for (int k = 0; k < Ds; ++k) phase += static_cast<Real>(flatM[m][k]) * w[k] * xc[k + 1];
        const Real cph = std::cos(phase), sph = std::sin(phase);  // exp(-i phase): re=cos, im=-sin
        for (int pidx = 0; pidx < npl; ++pidx) {
          const Real f = Fpl[pidx];
          lre[pidx * nm + m] += f * cph;
          lim[pidx * nm + m] += -f * sph;
        }
      }
    }
    #pragma omp critical
    {
      for (int a = 0; a < nacc; ++a) { re[a] += lre[a]; im[a] += lim[a]; }
    }
  }

  // S(group) = average over (plane, momentum-in-group) of |Ftilde|^2.
  std::vector<Real> acc(ng, 0.0);
  std::vector<int>  cnt(ng, 0);
  for (int pidx = 0; pidx < npl; ++pidx)
    for (int m = 0; m < nm; ++m) {
      const int g = flatG[m];
      const int a = pidx * nm + m;
      acc[g] += re[a] * re[a] + im[a] * im[a];
      ++cnt[g];
    }
  for (int g = 0; g < ng; ++g) S[g] = cnt[g] ? acc[g] / cnt[g] : 0.0;
  return S;
}

// ---- pole fit + jackknife -------------------------------------------------
struct PhotonMassFit {
  Real m2 = 0.0;       // intercept/slope of the R^{-1} = a + b*phat2 fit
  Real m2_err = 0.0;   // delete-1 jackknife error over configs
  Real m_gamma = 0.0;  // sqrt(max(0, m2))
  // R(phat2) table from the full (all-config) mean, for visualization:
  std::vector<Real> phat2;   // distinct phat2 (== mom.phat2)
  std::vector<Real> R;       // R(p) = <S(p)>/phat2
  std::vector<Real> R_err;   // SEM of R(p) (delete-1 jackknife of a mean)
};

// Linear least-squares of y = a + b x over n points (returns m2 = a/b, or NaN if
// the slope is ~0). Equal weights (the cubic-copy averaging already balances stats).
inline bool linfit_m2(const std::vector<Real>& x, const std::vector<Real>& y,
                      Real& m2_out) {
  const int n = static_cast<int>(x.size());
  if (n < 2) return false;
  Real Sx = 0, Sy = 0, Sxx = 0, Sxy = 0;
  for (int i = 0; i < n; ++i) { Sx += x[i]; Sy += y[i]; Sxx += x[i] * x[i]; Sxy += x[i] * y[i]; }
  const Real den = n * Sxx - Sx * Sx;
  if (std::fabs(den) < 1e-300) return false;
  const Real b = (n * Sxy - Sx * Sy) / den;          // slope
  const Real a = (Sy - b * Sx) / n;                  // intercept
  if (std::fabs(b) < 1e-300) return false;
  m2_out = a / b;                                    // m^2 = (m^2/Z)/(1/Z)
  return true;
}

// Given the accumulated per-config structure factors perConfig[c][g] (g indexes
// the phat2 groups in `mom`), return the pole fit with a delete-1 jackknife error.
// `nfit` (default: use all groups) is how many lowest-phat2 groups enter the fit.
template <int D>
PhotonMassFit photon_mass_fit(const std::vector<std::vector<Real>>& perConfig,
                              const PhotonMomenta<D>& mom, int nfit = -1) {
  PhotonMassFit out;
  const int K  = static_cast<int>(perConfig.size());
  const int ng = mom.n_groups();
  if (K == 0 || ng == 0) return out;
  const int nf = (nfit <= 0 || nfit > ng) ? ng : nfit;

  out.phat2.assign(mom.phat2.begin(), mom.phat2.end());

  // Group totals (for fast delete-1 means).
  std::vector<Real> tot(ng, 0.0);
  for (const auto& s : perConfig)
    for (int g = 0; g < ng; ++g) tot[g] += s[g];

  // Full-sample mean S, R, and the fit.
  std::vector<Real> Sbar(ng), x(nf), y(nf);
  for (int g = 0; g < ng; ++g) Sbar[g] = tot[g] / K;
  out.R.assign(ng, 0.0);
  for (int g = 0; g < ng; ++g) out.R[g] = Sbar[g] / mom.phat2[g];
  for (int i = 0; i < nf; ++i) { x[i] = mom.phat2[i]; y[i] = 1.0 / out.R[i]; }
  Real m2_full = 0.0;
  const bool ok = linfit_m2(x, y, m2_full);
  out.m2 = ok ? m2_full : std::numeric_limits<Real>::quiet_NaN();

  // Delete-1 jackknife over configs for m2 and R errors.
  out.R_err.assign(ng, 0.0);
  if (K >= 2) {
    std::vector<Real> jk; jk.reserve(K);
    std::vector<std::vector<Real>> Rjk(ng, std::vector<Real>(K, 0.0));
    for (int k = 0; k < K; ++k) {
      std::vector<Real> Sj(ng), Rj(ng), xj(nf), yj(nf);
      for (int g = 0; g < ng; ++g) {
        Sj[g] = (tot[g] - perConfig[k][g]) / (K - 1);
        Rj[g] = Sj[g] / mom.phat2[g];
        Rjk[g][k] = Rj[g];
      }
      for (int i = 0; i < nf; ++i) { xj[i] = mom.phat2[i]; yj[i] = 1.0 / Rj[i]; }
      Real m2k = 0.0;
      if (linfit_m2(xj, yj, m2k)) jk.push_back(m2k);
    }
    if (jk.size() >= 2) {
      const int Kj = static_cast<int>(jk.size());
      Real jbar = 0.0; for (Real v : jk) jbar += v; jbar /= Kj;
      Real sw = 0.0; for (Real v : jk) sw += (v - jbar) * (v - jbar);
      // Report the jackknife-bias-aware central value as the full-sample estimate,
      // with the standard delete-1 jackknife standard error.
      out.m2_err = std::sqrt((static_cast<Real>(Kj - 1) / Kj) * sw);
    }
    // Per-group R error (delete-1 jackknife of a mean == SEM).
    for (int g = 0; g < ng; ++g) {
      Real rbar = 0.0; for (int k = 0; k < K; ++k) rbar += Rjk[g][k]; rbar /= K;
      Real sw = 0.0; for (int k = 0; k < K; ++k) sw += (Rjk[g][k] - rbar) * (Rjk[g][k] - rbar);
      out.R_err[g] = std::sqrt((static_cast<Real>(K - 1) / K) * sw);
    }
  }

  out.m_gamma = std::sqrt(std::max<Real>(0.0, out.m2));
  return out;
}

}  // namespace u1
}  // namespace gh
