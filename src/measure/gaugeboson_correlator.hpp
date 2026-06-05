#pragma once
// Zero-momentum time-slice correlator + static non-abelian structure factor for
// the gauge-boson (massive vector) of a Higgsed SU(N) gauge theory. This is the
// rep-general analogue of u1/photon_mass.hpp + u1/photon_structure.hpp: it mirrors
// those templates one-for-one, replacing the abelian reduced plaquette angle
//   F^{U(1)}_{mu nu}(x) = reduce_angle(theta_{mu nu}(x))   (a single Real)
// with the gauge-invariant, channel-correct scalar
//   O_{mu nu}(x) = sum_a n^a(x) F^a_{mu nu}(x) = dotAlg( n(x), F_{mu nu}(x) )
// from src/measure/gaugeboson_op.hpp, where n^a is the rep-R isospin density and
// F^a_{mu nu} the Hermitian su(N) plaquette field strength. Because O is already a
// gauge SINGLET Real per (site, plane), the entire abelian correlator backend --
// zero-momentum timeslice projection, connected time/component-averaged C(dt),
// static structure factor, and the R(p)^{-1} = a + b*phat2 pole fit -- carries over
// verbatim; only the per-(site,plane) field-strength value changes.
//
// TIME DIRECTION: mu_t = kGaugeBosonTimeDir = 0 (coordinate x[0]), as in the U(1)
// templates and gaugeboson_op.hpp.
//
// (1) gaugeboson_timeslice_field<D,N> : S[t][i] = sum_{x:x[0]=t} O_{0,i}(x),
//     i = 0..D-2 over the D-1 spatial directions. MIRRORS u1::photon_timeslice_field.
// (2) gaugeboson_correlator : the connected, time-translation- and component-
//     averaged C(dt) built from per-config timeslice fields. The U(1)
//     photon_correlator is rep-agnostic (it only sees a vector<vector<vector<Real>>>),
//     so it would be directly reusable; we clone it here so this header is
//     self-contained (no u1/ dependency) -- the math is identical. Feed C(dt) to
//     cosh_effective_mass + plateau from measure/correlator.hpp.
// (3) gaugeboson_structure_factor<D,N> : mirrors u1::photon_structure_factor but
//     replaces the scalar reduced-plaquette per-plane value with
//     dotAlg( n^a(x), F^a_{ij}(x) ) = O_{ij}(x) for the SPATIAL plane (i,j). The
//     momentum set (photon_momenta) and the pole fit (photon_mass_fit, the
//     R(p)^{-1} = a + b*phat2 fit + jackknife) are reused VERBATIM from
//     u1/photon_structure.hpp -- they are pure post-processing on per-config arrays.
#include "core/config.hpp"
#include "core/geometry.hpp"
#include "core/fields.hpp"
#include "core/linalg.hpp"
#include "group/algebra.hpp"
#include "rep/representation.hpp"
#include "measure/gaugeboson_op.hpp"
#include "measure/correlator.hpp"
// photon_momenta + photon_mass_fit (the lattice-momentum set and the R^{-1}=a+b*phat2
// pole fit with jackknife) are reused VERBATIM; they operate purely on per-config
// scalar arrays and need no abelian field. We pull them in and re-export below.
#include "u1/photon_structure.hpp"
#include <vector>
#include <array>
#include <cmath>
#include <cstdint>

namespace gh {
namespace measure {

// ---------------------------------------------------------------------------
// (1) Zero-spatial-momentum time-slice projection of the electric operator.
//     S[t][i] = sum_{x: x[0]=t} O_{0,i+offset}(x), i = 0..D-2 over spatial dirs.
//     MIRRORS u1::photon_timeslice_field (same loop structure, time dir 0); the
//     only change is the per-site value: O_{0i}(x) instead of reduced F_{0i}(x).
// ---------------------------------------------------------------------------
template <int D, int N>
std::vector<std::vector<Real>> gaugeboson_timeslice_field(
    const Representation<N>& rep, const GaugeField<D, N>& U,
    const std::vector<DVec>& phi, const Lattice<D>& lat) {
  const int mut = kGaugeBosonTimeDir;
  const int Lt  = lat.L[mut];
  std::vector<std::vector<Real>> S(Lt, std::vector<Real>(D - 1, 0.0));
  std::array<int, D> xc{};
  for (std::int64_t s = 0; s < lat.vol; ++s) {
    lat.coords(s, xc);
    const int t = xc[mut];
    const DVec& phi_s = phi[static_cast<std::size_t>(s)];
    int ci = 0;                                     // compact index over spatial dirs
    for (int i = 0; i < D; ++i) {
      if (i == mut) continue;
      S[t][ci] += op_electric<D, N>(rep, U, phi_s, s, i);  // O_{0 i}(x)
      ++ci;
    }
  }
  return S;
}

// ---------------------------------------------------------------------------
// (2) Connected, time-translation- and component-averaged gauge-boson correlator.
//     perConfig[c][t][i] = S_i(t) for config c (output of gaugeboson_timeslice_field).
//     C(dt) = (1/(Nc*Lt*(D-1))) sum_{c,t,i} S_i^{(c)}(t+dt) S_i^{(c)}(t)
//             - avg S * avg S    (connected; disconnected piece subtracted explicitly).
//     This is a verbatim clone of u1::photon_correlator (rep-agnostic: it only sees
//     a nested Real array). Feed C(dt) into cosh_effective_mass / plateau.
// ---------------------------------------------------------------------------
inline std::vector<Real> gaugeboson_correlator(
    const std::vector<std::vector<std::vector<Real>>>& perConfig) {
  const int Nc = static_cast<int>(perConfig.size());
  if (Nc == 0) return {};
  const int Lt = static_cast<int>(perConfig[0].size());
  if (Lt == 0) return {};
  const int ncomp = static_cast<int>(perConfig[0][0].size());

  // <S>: mean over all (config, time, component) samples.
  Real meanS = 0.0;
  std::int64_t nS = 0;
  for (const auto& cfg : perConfig)
    for (const auto& slice : cfg)
      for (Real v : slice) { meanS += v; ++nS; }
  meanS = nS ? meanS / nS : 0.0;

  std::vector<Real> C(Lt, 0.0);
  for (int dt = 0; dt < Lt; ++dt) {
    Real acc = 0.0;
    std::int64_t cnt = 0;
    for (const auto& cfg : perConfig)
      for (int t = 0; t < Lt; ++t) {
        const int tp = (t + dt) % Lt;               // periodic in time
        for (int i = 0; i < ncomp; ++i) {
          acc += cfg[tp][i] * cfg[t][i];
          ++cnt;
        }
      }
    C[dt] = (cnt ? acc / cnt : 0.0) - meanS * meanS;
  }
  return C;
}

// ---------------------------------------------------------------------------
// (3) Static (zero-frequency) non-abelian magnetic structure factor.
//     Mirrors u1::photon_structure_factor; per spatial plane (i,j) the abelian
//     reduced-plaquette value is replaced by the gauge-invariant operator
//       O_{ij}(x) = dotAlg( n^a(x), F^a_{ij}(x) ),
//     n^a the rep-R isospin density, F^a_{ij} the Hermitian su(N) field strength.
//     Ftilde_{ij}(p) = sum_x exp(-i sum_k p_k x_{k+1}) O_{ij}(x), summed over ALL
//     sites x (x_0 summed -> p_0 = 0). Output length == mom.n_groups(); entry g is
//     < |Ftilde_{ij}(p)|^2 >_{spatial planes i<j, momenta in group g}.
//
//     The momentum set type u1::PhotonMomenta<D> and builder u1::photon_momenta<D>
//     are reused (lattice momenta are theory-independent), as is the pole fit
//     u1::photon_mass_fit<D> (R(p)^{-1} = a + b*phat2). They are re-exported below.
// ---------------------------------------------------------------------------
template <int D, int N>
std::vector<Real> gaugeboson_structure_factor(const Representation<N>& rep,
                                              const GaugeField<D, N>& U,
                                              const std::vector<DVec>& phi,
                                              const Lattice<D>& lat,
                                              const u1::PhotonMomenta<D>& mom) {
  constexpr int Ds = D - 1;
  const int ng = mom.n_groups();
  std::vector<Real> S(ng, 0.0);
  if (ng == 0) return S;

  // Spatial planes (i,j), 1 <= i < j <= D-1 (excludes the time dir 0).
  std::vector<std::pair<int, int>> planes;
  for (int i = 1; i < D; ++i)
    for (int j = i + 1; j < D; ++j) planes.emplace_back(i, j);
  const int npl = static_cast<int>(planes.size());

  // Flatten the momentum list: flatM[m] = mode vector, flatG[m] = phat2 group.
  std::vector<std::array<int, Ds>> flatM;
  std::vector<int> flatG;
  for (int g = 0; g < ng; ++g)
    for (const auto& nv : mom.modes[g]) {
      std::array<int, Ds> a{};
      for (int k = 0; k < Ds; ++k) a[k] = nv[k];
      flatM.push_back(a); flatG.push_back(g);
    }
  const int nm = static_cast<int>(flatM.size());

  // Per-direction angular frequency w_k = 2 pi / L_{k+1} (dir 0 is time).
  std::array<Real, Ds> w{};
  for (int k = 0; k < Ds; ++k) w[k] = 2.0 * kPi / static_cast<Real>(lat.L[k + 1]);

  // Accumulators: re/im of Ftilde[plane][momentum] in a flat [npl*nm] buffer.
  const int nacc = npl * nm;
  std::vector<Real> re(nacc, 0.0), im(nacc, 0.0);

  #pragma omp parallel
  {
    std::vector<Real> lre(nacc, 0.0), lim(nacc, 0.0);
    std::vector<Real> Opl(npl, 0.0);   // per-site gauge-boson operator, by spatial plane
    std::array<int, D> xc{};
    #pragma omp for schedule(static) nowait
    for (std::int64_t s = 0; s < lat.vol; ++s) {
      lat.coords(s, xc);
      const DVec& phi_s = phi[static_cast<std::size_t>(s)];
      // Gauge-invariant operator O_{ij}(x) = dotAlg(n(x), F_{ij}(x)) for every
      // spatial plane at this site. (Cache once; reused across all momenta.)
      const AlgVec<N> n = isospin_density<N>(rep, phi_s);
      for (int pidx = 0; pidx < npl; ++pidx) {
        const AlgVec<N> F =
            fieldstrength_alg<D, N>(U, s, planes[pidx].first, planes[pidx].second);
        Opl[pidx] = dotAlg<N>(n, F);
      }
      for (int m = 0; m < nm; ++m) {
        Real phase = 0.0;
        for (int k = 0; k < Ds; ++k) phase += static_cast<Real>(flatM[m][k]) * w[k] * xc[k + 1];
        const Real cph = std::cos(phase), sph = std::sin(phase);  // exp(-i phase): re=cos, im=-sin
        for (int pidx = 0; pidx < npl; ++pidx) {
          const Real f = Opl[pidx];
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

// Re-export the reused, theory-independent momentum set + pole fit so callers can
// use them through the measure:: namespace without reaching into u1::.
template <int D>
using GaugeBosonMomenta = u1::PhotonMomenta<D>;
template <int D>
inline GaugeBosonMomenta<D> gaugeboson_momenta(const Lattice<D>& lat) {
  return u1::photon_momenta<D>(lat);          // identical lattice-momentum set
}
using GaugeBosonMassFit = u1::PhotonMassFit;
template <int D>
inline GaugeBosonMassFit gaugeboson_mass_fit(
    const std::vector<std::vector<Real>>& perConfig,
    const u1::PhotonMomenta<D>& mom, int nfit = -1) {
  return u1::photon_mass_fit<D>(perConfig, mom, nfit);  // R^{-1}=a+b*phat2 + jackknife
}

}  // namespace measure
}  // namespace gh
