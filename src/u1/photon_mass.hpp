#pragma once
// Compact-U(1) photon mass in D dimensions from the zero-spatial-momentum
// transverse field-strength correlator.
//
// TIME DIRECTION: mu_t = 0 (coordinate x[0]). A "timeslice" t collects all sites
// with x[0] == t; there are L[0] of them, periodic in t.
//
// PHYSICAL FIELD STRENGTH: the bare plaquette angle theta_{mu nu}(s) (see
// u1::plaq_angle) lives on the whole real line; the gauge-invariant, physical
// field strength is its REDUCED (principal-value) form in [-pi,pi),
//   F_{mu nu}(s) = theta_{mu nu}(s) - 2*pi*round(theta_{mu nu}(s)/(2*pi)),
// which subtracts the integer Dirac-string / monopole flux 2*pi*Z. F is gauge
// invariant: under theta_mu(x) -> theta_mu(x)+a(x)-a(x+mu) the plaquette angle is
// invariant, hence so is its principal value.
//
// ZERO-MOMENTUM PROJECTION: for each timeslice t and each spatial direction
// i != 0 we sum the reduced electric field F_{0 i} over the spatial sites of the
// slice, S_i(t) = sum_{x: x[0]=t} F_{0 i}(x). This projects onto zero spatial
// momentum, leaving a function of t only -- a transverse (i runs over the D-1
// spatial directions) interpolating operator for the photon.
//
// CORRELATOR: the connected, time-translation- and component-averaged
//   C(dt) = avg_i avg_t [ <S_i(t+dt) S_i(t)> - <S_i>^2 ],
// built from the per-config timeslice fields. In the Coulomb (massless) phase
// m_eff -> 0 (power-law decay); in the Higgs phase C(dt) shows a finite plateau
// m_eff ~ q*sqrt(kappa) (the photon eats the Goldstone and becomes massive).
#include "core/config.hpp"
#include "core/geometry.hpp"
#include "u1/u1.hpp"
#include <vector>
#include <cmath>
#include <cstdint>

namespace gh {
namespace u1 {

constexpr int kPhotonTimeDir = 0;  // mu_t -- the Euclidean time direction

// Principal value of an angle into [-pi, pi): subtract 2*pi*round(x/2pi).
// (std::round is half-away-from-zero, so the boundary +pi maps to -pi; this is a
// measure-zero case never hit by float HMC configs and physics-irrelevant.)
inline Real reduce_angle(Real x) {
  return x - 2.0 * kPi * std::round(x / (2.0 * kPi));
}

// Reduced (physical, principal-value) field strength F_{mu nu}(s) in [-pi,pi).
// NOTE: u1/monopole.hpp defines an identical reduced_plaq_angle<D> template. A driver
// that includes BOTH headers (e.g. u1_scan.cpp) would hit an in-TU template
// redefinition (ODR error). Such a driver may #define GH_U1_HAVE_REDUCED_PLAQ_ANGLE
// before including this header (after monopole.hpp) to suppress this definition and
// reuse monopole.hpp's identical one. Standalone users (e.g. test_photon_mass.cpp)
// leave the macro unset and get the definition here.
#ifndef GH_U1_HAVE_REDUCED_PLAQ_ANGLE
#define GH_U1_HAVE_REDUCED_PLAQ_ANGLE
template <int D>
inline Real reduced_plaq_angle(const std::vector<Real>& th, const Lattice<D>& lat,
                               std::int64_t s, int mu, int nu) {
  return reduce_angle(plaq_angle<D>(th, lat, s, mu, nu));
}
#endif

// Per-config zero-spatial-momentum projection of the reduced electric field.
// Returns timeslices[t][i] = S_i(t) = sum_{x: x[mu_t]=t} F_{mu_t, i+offset}(x),
// where i = 0..D-2 indexes the D-1 spatial directions (all dirs != mu_t). The
// outer length is L[mu_t]; the inner length is D-1.
template <int D>
std::vector<std::vector<Real>> photon_timeslice_field(const std::vector<Real>& th,
                                                      const Lattice<D>& lat) {
  const int mut = kPhotonTimeDir;
  const int Lt  = lat.L[mut];
  std::vector<std::vector<Real>> S(Lt, std::vector<Real>(D - 1, 0.0));
  std::array<int, D> xc{};
  for (std::int64_t s = 0; s < lat.vol; ++s) {
    lat.coords(s, xc);
    const int t = xc[mut];
    int ci = 0;                                     // compact index over spatial dirs
    for (int i = 0; i < D; ++i) {
      if (i == mut) continue;
      S[t][ci] += reduced_plaq_angle<D>(th, lat, s, mut, i);
      ++ci;
    }
  }
  return S;
}

// Connected, time-translation- and component-averaged photon correlator.
//   perConfig[c][t][i] = S_i(t) for configuration c (output of photon_timeslice_field).
// Returns C of length Lt with
//   C(dt) = (1/(Nc*Lt*(D-1))) sum_{c,t,i} S_i^{(c)}(t+dt) S_i^{(c)}(t)
//           - avg_{c,t,i} S_i^{(c)}(t) * avg_{c,t,i} S_i^{(c)}(t).
// The disconnected piece <S>^2 vanishes by parity for symmetric ensembles but is
// subtracted explicitly so the estimator is unbiased on finite, asymmetric data.
inline std::vector<Real> photon_correlator(
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

}  // namespace u1
}  // namespace gh
