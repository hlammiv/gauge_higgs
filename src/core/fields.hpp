#pragma once
// Field containers over an arbitrary-D lattice. Links live on (site, mu); the
// link momentum is an su(N) algebra vector per (site, mu).
#include "core/geometry.hpp"
#include "core/rng.hpp"
#include "group/sun.hpp"
#include <vector>

namespace gh {

// Gauge field: U_mu(x) in SU(N), fundamental rep. Layout [site*D + mu].
template <int D, int N>
struct GaugeField {
  const Lattice<D>* lat = nullptr;
  std::vector<Cmat<N>> u;

  explicit GaugeField(const Lattice<D>& L) : lat(&L), u(static_cast<std::size_t>(L.vol) * D) {}

  Cmat<N>&       operator()(std::int64_t s, int mu)       { return u[s * D + mu]; }
  const Cmat<N>& operator()(std::int64_t s, int mu) const { return u[s * D + mu]; }

  void cold() { for (auto& m : u) m = Cmat<N>::identity(); }

  // Hot start: each link = exp(i sigma * random algebra element).
  void hot(const Rng& rng, Real sigma = 1.0) {
    #pragma omp parallel for schedule(static)
    for (std::int64_t s = 0; s < lat->vol; ++s)
      for (int mu = 0; mu < D; ++mu) {
        AlgVec<N> v{};
        for (int a = 0; a < n_gen<N>(); ++a)
          v[a] = sigma * rng.gauss(Rng::key(1, s, mu, a));
        (*this)(s, mu) = expi<N>(alg_to_mat<N>(v));
      }
  }

  void reunitarize_all() { for (auto& m : u) reunitarize<N>(m); }

  Real max_unitarity_violation() const {
    Real worst = 0.0;
    for (const auto& m : u) worst = std::max(worst, unitarity_violation<N>(m));
    return worst;
  }
};

// Link-momentum field: su(N) algebra vector per (site, mu).
template <int D, int N>
struct LinkMom {
  const Lattice<D>* lat = nullptr;
  std::vector<AlgVec<N>> p;

  explicit LinkMom(const Lattice<D>& L) : lat(&L), p(static_cast<std::size_t>(L.vol) * D) {}

  AlgVec<N>&       operator()(std::int64_t s, int mu)       { return p[s * D + mu]; }
  const AlgVec<N>& operator()(std::int64_t s, int mu) const { return p[s * D + mu]; }

  void zero() { for (auto& v : p) v.fill(0.0); }

  // Gaussian heatbath: P^a ~ N(0,1). `stream` distinguishes trajectories.
  void refresh(const Rng& rng, std::uint64_t stream) {
    #pragma omp parallel for schedule(static)
    for (std::int64_t s = 0; s < lat->vol; ++s)
      for (int mu = 0; mu < D; ++mu) {
        AlgVec<N>& v = (*this)(s, mu);
        for (int a = 0; a < n_gen<N>(); ++a) v[a] = rng.gauss(Rng::key(stream, s, mu, a));
      }
  }

  // Kinetic energy 1/2 sum (P^a)^2.
  Real kinetic() const {
    Real e = 0.0;
    for (const auto& v : p)
      for (int a = 0; a < n_gen<N>(); ++a) e += v[a] * v[a];
    return 0.5 * e;
  }
};

}  // namespace gh
