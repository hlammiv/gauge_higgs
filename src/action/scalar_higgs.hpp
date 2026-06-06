#pragma once
// Gauge-covariant scalar (Higgs) action in an arbitrary irrep R, plus the two MD
// forces it produces: the force on phi, and the matter back-reaction on the gauge
// links. All sign/factor choices are FD-validated in test_scalar.cpp.
// See theory_notes §1.2, §3.2, §3.3, conventions.md.
#include "core/scalar_field.hpp"
#include "core/fields.hpp"
#include "rep/representation.hpp"
#include "core/profile.hpp"  // GH_PROF_* (no-op unless -DGH_PROFILE)
#include <vector>

namespace gh {

// Pluggable on-site scalar potential V(phi). value() and dV/dphi^* (Wirtinger).
// QuarticPotential reproduces the original phi^dag phi + lambda(phi^dag phi - 1)^2.
// The general multi-invariant potential (MultiInvariantPotential, scalar_invariants.hpp)
// derives from this so the HMC can lock the VEV onto a chosen discrete subgroup H.
template <int N>
struct OnsitePotential {
  virtual ~OnsitePotential() = default;
  virtual Real value(const DVec& phi) const = 0;
  virtual DVec dV_dphibar(const DVec& phi) const = 0;  // dV/dphi^*
};

template <int N>
struct QuarticPotential : OnsitePotential<N> {
  Real lambda;
  explicit QuarticPotential(Real l) : lambda(l) {}
  Real value(const DVec& phi) const override {
    const Real n2 = phi.norm2();
    return n2 + lambda * (n2 - 1.0) * (n2 - 1.0);
  }
  DVec dV_dphibar(const DVec& phi) const override {
    const Real n2 = phi.norm2();
    const Real c = 1.0 + 2.0 * lambda * (n2 - 1.0);
    DVec g(phi.size());
    for (int k = 0; k < phi.size(); ++k) g(k) = Complex(c, 0) * phi(k);
    return g;
  }
};

// S_H = sum_x V(phi_x) - kappa sum_{x,mu} 2 Re[ phi_x^dag D^(R)(U_mu(x)) phi_{x+mu} ].
template <int D, int N>
Real scalar_action(const CVecField<D>& phi, const GaugeField<D, N>& U,
                   const Representation<N>& rep, Real kappa, const OnsitePotential<N>& pot) {
  const Lattice<D>& lat = *phi.lat;
  Real S = 0.0;
  #pragma omp parallel for schedule(static) reduction(+:S)
  for (std::int64_t s = 0; s < lat.vol; ++s) {
    const DVec px = phi.get(s);
    S += pot.value(px);
    for (int mu = 0; mu < D; ++mu) {
      const DVec Dpy = rep.rotate(U(s, mu), phi.get(lat.neighbor_fwd(s, mu)));
      S -= 2.0 * kappa * dot(px, Dpy).real();
    }
  }
  return S;
}
// Backward-compatible overload: lambda -> QuarticPotential.
template <int D, int N>
Real scalar_action(const CVecField<D>& phi, const GaugeField<D, N>& U,
                   const Representation<N>& rep, Real kappa, Real lambda) {
  QuarticPotential<N> q(lambda);
  return scalar_action<D, N>(phi, U, rep, kappa, q);
}

// MD force on phi: F_phi(x) = -2 dS_H/dphi_x^* (real-coordinate gradient as a complex
// vector; pi-kick is pi += dt F_phi).  On-site part = -2 dV/dphi^*; hopping part below.
// Dcache (optional): one D^(R)(U_link) per link [layout s*D+mu], built once per kick so the
// rep matrix is NOT recomputed (~3x) per link. Pure memoization -> bit-identical result.
template <int D, int N>
void scalar_force(const CVecField<D>& phi, const GaugeField<D, N>& U,
                  const Representation<N>& rep, Real kappa, const OnsitePotential<N>& pot,
                  CVecField<D>& Fphi, const std::vector<DMat>* Dcache = nullptr) {
  const Lattice<D>& lat = *phi.lat;
  #pragma omp parallel for schedule(static)
  for (std::int64_t s = 0; s < lat.vol; ++s) {
    const DVec px = phi.get(s);
    DVec g;
    { GH_PROF_SCOPE(scalar_force_onsite);
      g = pot.dV_dphibar(px); }
    DVec F(rep.d);
    for (int k = 0; k < rep.d; ++k) F(k) = -2.0 * g(k);  // on-site: -2 dV/dphi^*
    { GH_PROF_SCOPE(scalar_force_hopping);
    for (int mu = 0; mu < D; ++mu) {
      const std::int64_t yf = lat.neighbor_fwd(s, mu);
      const std::int64_t yb = lat.neighbor_bwd(s, mu);
      DVec Dpf, Dpb;
      if (Dcache) {
        Dpf = rep.rotate_cached    ((*Dcache)[s  * D + mu], phi.get(yf));   // D(U_mu(x)) phi_{x+mu}
        Dpb = rep.rotate_dag_cached((*Dcache)[yb * D + mu], phi.get(yb));   // D(U_mu(x-mu))^dag phi_{x-mu}
      } else {
        Dpf = rep.rotate    (U(s, mu),  phi.get(yf));
        Dpb = rep.rotate_dag(U(yb, mu), phi.get(yb));
      }
      for (int k = 0; k < rep.d; ++k) F(k) += 2.0 * kappa * (Dpf(k) + Dpb(k));
    } }
    Fphi.set(s, F);
  }
}
// Backward-compatible overload: lambda -> QuarticPotential.
template <int D, int N>
void scalar_force(const CVecField<D>& phi, const GaugeField<D, N>& U,
                  const Representation<N>& rep, Real kappa, Real lambda, CVecField<D>& Fphi,
                  const std::vector<DMat>* Dcache = nullptr) {
  QuarticPotential<N> q(lambda);
  scalar_force<D, N>(phi, U, rep, kappa, q, Fphi, Dcache);
}

// Matter back-reaction on the gauge links: dS_H/d(omega^a) for link (x,mu) = -kappa g^a,
// with g^a = rep.hop_link_g(...). Added into the link force F (alongside the gauge force);
// the link-momentum kick is  P -= dt * F.
template <int D, int N>
void add_matter_link_force(const CVecField<D>& phi, const GaugeField<D, N>& U,
                           const Representation<N>& rep, Real kappa, LinkMom<D, N>& F,
                           const std::vector<DMat>* Dcache = nullptr) {
  const Lattice<D>& lat = *phi.lat;
  #pragma omp parallel for schedule(static)
  for (std::int64_t s = 0; s < lat.vol; ++s)
    for (int mu = 0; mu < D; ++mu) {
      const std::int64_t y = lat.neighbor_fwd(s, mu);
      const AlgVec<N> g = Dcache
          ? rep.hop_link_g_cached((*Dcache)[s * D + mu], phi.get(s), phi.get(y))
          : rep.hop_link_g(U(s, mu), phi.get(s), phi.get(y));
      AlgVec<N>& Fa = F(s, mu);
      for (int a = 0; a < n_gen<N>(); ++a) Fa[a] += -kappa * g[a];
    }
}

// Gauge-invariant hopping energy (transition locator): L_link = <Re phi_x^dag D(U) phi_{x+mu}>
// averaged over all links.
template <int D, int N>
Real link_energy(const CVecField<D>& phi, const GaugeField<D, N>& U,
                 const Representation<N>& rep) {
  const Lattice<D>& lat = *phi.lat;
  Real s = 0.0;
  #pragma omp parallel for schedule(static) reduction(+:s)
  for (std::int64_t x = 0; x < lat.vol; ++x)
    for (int mu = 0; mu < D; ++mu) {
      const DVec Dpy = rep.rotate(U(x, mu), phi.get(lat.neighbor_fwd(x, mu)));
      s += dot(phi.get(x), Dpy).real();
    }
  return s / (static_cast<Real>(lat.vol) * D);
}

// Higgs length L_phi = <phi^dag phi>.
template <int D>
Real higgs_length(const CVecField<D>& phi) {
  const Lattice<D>& lat = *phi.lat;
  Real s = 0.0;
  #pragma omp parallel for schedule(static) reduction(+:s)
  for (std::int64_t x = 0; x < lat.vol; ++x) s += phi.get(x).norm2();
  return s / static_cast<Real>(lat.vol);
}

}  // namespace gh
