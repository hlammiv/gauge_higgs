#pragma once
// Gauge-boson vector operator: a gauge-invariant, channel-correct (linear in A)
// interpolating field for the massive vector boson of a Higgsed gauge theory.
//
//   O_{mu nu}(x) = sum_a n^a(x) F^a_{mu nu}(x)
//
// with
//   n^a(x)       = Re[ phi(x)^dag T^a_R phi(x) ]              (rep-R isospin density)
//   F^a_{mu nu}  = mat_to_alg<N>( proj_TA<N>(U_{mu nu}(x)) * (-i) )  (Hermitian field strength)
//   U_{mu nu}(x) = U(s,mu) U(s+mu,nu) U(s+nu,mu)^dag U(s,nu)^dag     (plaquette)
//
// O contracts the rep-R isospin density (a "color" vector in su(N)*) against the
// plaquette field strength (a vector in su(N)) via the Killing form dotAlg, so it
// is a gauge SINGLET: under g(x) in SU(N) [phi -> D^(R)(g)phi, U_mu(x) -> g(x)U_mu(x)g(x+mu)^dag]
// both n^a and F^a rotate by the SAME adjoint matrix Ad(g) and the contraction is invariant.
// It is also LINEAR in the gauge potential A (one power of F), so it couples to the
// one-vector-boson state -- the channel the gauge-boson mass lives in.
//
// LOAD-BEARING DETAIL (do NOT regress): the bare mat_to_alg(proj_TA(U)) is
// IDENTICALLY ZERO. proj_TA returns an ANTI-Hermitian matrix, and mat_to_alg<N>(M) =
// 2 Re Tr(T^a M) projects onto the HERMITIAN part, whose real trace vanishes for an
// anti-Hermitian argument. Multiplying by Complex(0,-1) = -i turns the anti-Hermitian
// plaquette projection into the physical HERMITIAN field strength, which has nonzero
// algebra components. The * Complex(0,-1) factor below is MANDATORY.
//
// MASS: use the ELECTRIC components O_{0 i}, i = 1..D-1 (time direction mu_t = 0).
// These are the spin-1 / cubic-T1 vector that overlaps the massive photon.
//
// API: src/rep/representation.hpp (Representation<N>: d x d Hermitian T[a], rotate,
// dot on DVec), src/group/algebra.hpp (proj_TA, mat_to_alg, dotAlg, AlgVec, n_gen),
// src/core/fields.hpp (GaugeField<D,N>, lat.neighbor_fwd), src/core/linalg.hpp.
#include "rep/representation.hpp"
#include "core/fields.hpp"
#include "core/linalg.hpp"
#include "group/algebra.hpp"
#include <cstdint>

namespace gh {
namespace measure {

constexpr int kGaugeBosonTimeDir = 0;  // mu_t -- the Euclidean time direction

// Rep-R isospin density n^a(x): a su(N)*-valued "color direction" carried by the
// scalar at site x. It must (i) rotate by Ad(g) under a gauge transform (so the
// contraction dotAlg(n, F) is invariant) and (ii) be NONZERO on the actual HMC
// field.  The default is the moment-map bilinear
//     n^a(x) = Re[ phi(x)^dag T^a_R phi(x) ]
// which is the correct, nonzero isospin density for a COMPLEX rep (T^a_R Hermitian
// => phi^dag T^a phi real; we take Re explicitly).
//
// ADJOINT (real) SPECIAL CASE -- load-bearing, do NOT regress.  For the adjoint
// representation the scalar field components phi^a ARE the algebra coordinates
// (phi = phi^a T^a_fund, a real vector in su(N)).  The moment-map bilinear
// Re[phi^dag T^a phi] then VANISHES IDENTICALLY on the physical field: the adjoint
// generators are T^a_adj = -i f^{abc} (pure-imaginary, antisymmetric), so for a
// REAL adjoint field phi (the only kind the HMC produces -- rep.real => imag==0)
// phi^T (-i f^a) phi = 0 by antisymmetry, and Re/Im both give 0.  The op- and
// correlator-validation tests masked this by feeding a COMPLEX random vector to the
// real adjoint rep; on a genuine real adjoint condensate the operator is identically
// zero (verified).  The physically-correct 't Hooft adjoint-photon direction is the
// adjoint field itself:  n^a = phi^a.  It rotates by Ad(g) (phi -> Ad(g) phi under a
// gauge transform), so dotAlg(phi, F) is the gauge-invariant, linear-in-A vector
// (W/photon) operator.  We select it for the adjoint via (rep.real && N-ality 0 &&
// d == dim adjoint), which is exactly the adjoint rep.
template <int N>
AlgVec<N> isospin_density(const Representation<N>& rep, const DVec& phi) {
  AlgVec<N> n{};
  // Adjoint (real, N-ality 0, d = N^2-1): the field components are the algebra
  // direction; the bilinear moment map degenerates to zero, so use n^a = phi^a.
  if (rep.real && rep.nality == 0 && rep.d == n_gen<N>()) {
    for (int a = 0; a < n_gen<N>(); ++a) n[a] = phi(a).real();
    return n;
  }
  for (int a = 0; a < n_gen<N>(); ++a)
    n[a] = dot(phi, rep.T[a] * phi).real();  // x^dag y convention in dot(); Re of a real#
  return n;
}

// Plaquette U_{mu nu}(x) in the fundamental (SU(N)) at site s, plane (mu,nu):
//   U(s,mu) U(s+mu,nu) U(s+nu,mu)^dag U(s,nu)^dag.
template <int D, int N>
Cmat<N> plaquette(const GaugeField<D, N>& U, std::int64_t s, int mu, int nu) {
  const Lattice<D>& lat = *U.lat;
  const std::int64_t s_pmu = lat.neighbor_fwd(s, mu);
  const std::int64_t s_pnu = lat.neighbor_fwd(s, nu);
  return U(s, mu) * U(s_pmu, nu) * U(s_pnu, mu).dagger() * U(s, nu).dagger();
}

// Hermitian plaquette field strength as su(N) algebra components:
//   F^a_{mu nu}(x) = mat_to_alg<N>( proj_TA<N>(U_{mu nu}(x)) * Complex(0,-1) ).
// The Complex(0,-1) = -i factor is MANDATORY (see header note): it converts the
// anti-Hermitian TA projection into the physical Hermitian field strength.
template <int D, int N>
AlgVec<N> fieldstrength_alg(const GaugeField<D, N>& U, std::int64_t s, int mu, int nu) {
  const Cmat<N> Umunu = plaquette<D, N>(U, s, mu, nu);
  const Cmat<N> Fherm = proj_TA<N>(Umunu) * Complex(0.0, -1.0);  // -i [U_munu]_TA (Hermitian)
  return mat_to_alg<N>(Fherm);
}

// Electric component of the gauge-boson operator at site s, spatial direction i:
//   O_{0 i}(s) = dotAlg( n(s), F_{0 i}(s) ),  i = 1..D-1.
// `i` is the spatial direction index (1..D-1); the time direction is kGaugeBosonTimeDir = 0.
template <int D, int N>
Real op_electric(const Representation<N>& rep, const GaugeField<D, N>& U,
                 const DVec& phi_s, std::int64_t s, int i) {
  const AlgVec<N> n = isospin_density<N>(rep, phi_s);
  const AlgVec<N> F = fieldstrength_alg<D, N>(U, s, kGaugeBosonTimeDir, i);
  return dotAlg<N>(n, F);
}

// ---------------------------------------------------------------------------
// Secondary cross-check operator (FMS gauge-covariant hopping; optional, no GEVP):
//   O^FMS_mu(x) = sum_a n^a(x) K^a_mu(x),
//   K^a_mu(x)   = Im[ phi(x)^dag T^a_R ( D^(R)(U(s,mu)) phi(x+mu) - phi(x) ) ].
// Built from the SAME isospin density and a gauge-covariant hop; gauge invariant by
// the same Ad(g) cancellation. Used only as an independent overlap, never in the fit.
// ---------------------------------------------------------------------------
template <int D, int N>
AlgVec<N> hop_density(const Representation<N>& rep, const GaugeField<D, N>& U,
                      const DVec& phi_s, const DVec& phi_spmu, std::int64_t s, int mu) {
  const DVec hopped = rep.rotate(U(s, mu), phi_spmu);  // D^(R)(U(s,mu)) phi(x+mu)
  DVec diff(hopped.size());
  for (int k = 0; k < diff.size(); ++k) diff(k) = hopped(k) - phi_s(k);
  AlgVec<N> K{};
  for (int a = 0; a < n_gen<N>(); ++a)
    K[a] = dot(phi_s, rep.T[a] * diff).imag();
  return K;
}

template <int D, int N>
Real op_fms(const Representation<N>& rep, const GaugeField<D, N>& U,
            const DVec& phi_s, const DVec& phi_spmu, std::int64_t s, int mu) {
  const AlgVec<N> n = isospin_density<N>(rep, phi_s);
  const AlgVec<N> K = hop_density<D, N>(rep, U, phi_s, phi_spmu, s, mu);
  return dotAlg<N>(n, K);
}

}  // namespace measure
}  // namespace gh
