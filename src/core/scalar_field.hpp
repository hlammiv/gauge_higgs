#pragma once
// Per-site complex d-vector field (scalar Higgs phi, and its conjugate momentum pi).
// Stored flat [site*d + k]. For a REAL representation the field carries d real dof:
// the imaginary parts are held at zero (momentum heatbath draws only real Gaussians).
#include "core/geometry.hpp"
#include "core/rng.hpp"
#include "core/linalg.hpp"

namespace gh {

template <int D>
struct CVecField {
  const Lattice<D>* lat = nullptr;
  int d = 0;
  std::vector<Complex> data;  // [site*d + k]

  CVecField(const Lattice<D>& L, int dim)
      : lat(&L), d(dim), data(static_cast<std::size_t>(L.vol) * dim, Complex(0, 0)) {}

  DVec get(std::int64_t s) const {
    DVec v(d);
    const std::size_t off = static_cast<std::size_t>(s) * d;
    for (int k = 0; k < d; ++k) v(k) = data[off + k];
    return v;
  }
  void set(std::int64_t s, const DVec& v) {
    const std::size_t off = static_cast<std::size_t>(s) * d;
    for (int k = 0; k < d; ++k) data[off + k] = v(k);
  }

  void zero() { std::fill(data.begin(), data.end(), Complex(0, 0)); }

  // Uniform cold start: phi_k = c (real). With c=1, |phi|^2 = d.
  void cold(Real c = 1.0) {
    for (auto& z : data) z = Complex(0, 0);
    for (std::int64_t s = 0; s < lat->vol; ++s) data[static_cast<std::size_t>(s) * d] = Complex(c, 0);
  }

  // Hot/heatbath: Gaussian components. If `real_rep`, imaginary parts stay zero.
  void gaussian(const Rng& rng, std::uint64_t stream, bool real_rep, Real sigma = 1.0) {
    #pragma omp parallel for schedule(static)
    for (std::int64_t s = 0; s < lat->vol; ++s) {
      const std::size_t off = static_cast<std::size_t>(s) * d;
      for (int k = 0; k < d; ++k) {
        const Real re = sigma * rng.gauss(Rng::key(stream, s, 2 * k));
        const Real im = real_rep ? 0.0 : sigma * rng.gauss(Rng::key(stream, s, 2 * k + 1));
        data[off + k] = Complex(re, im);
      }
    }
  }

  // Kinetic energy 1/2 sum |pi_k|^2 (correct for both real and complex reps).
  Real kinetic() const {
    Real e = 0.0;
    for (const auto& z : data) e += std::norm(z);
    return 0.5 * e;
  }
};

}  // namespace gh
