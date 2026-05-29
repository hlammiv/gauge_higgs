#pragma once
// Adjoint representation: real, orthogonal link [U^adj]_{ab} = 2 Tr(T^a U T^b U^dag),
// no matrix exponential needed; generators (T^a_adj)_{bc} = -i f^{abc}. d = N^2-1,
// N-ality 0 (center-blind -> genuine Higgs/confinement transition can exist).
// See theory_notes §1.4, §6.2, memory: adjoint-link-fast-path.
#include "rep/representation.hpp"

namespace gh {

template <int N>
struct AdjointRep : Representation<N> {
  AdjointRep() {
    const int n = n_gen<N>();
    this->d = n; this->nality = 0; this->real = true;
    const auto& f = structure_constants<N>();
    this->T.resize(n);
    for (int a = 0; a < n; ++a) {
      DMat M(n, n);
      for (int b = 0; b < n; ++b)
        for (int c = 0; c < n; ++c)
          M(b, c) = Complex(0.0, -f[(static_cast<std::size_t>(a) * n + b) * n + c]);  // -i f^{abc}
      this->T[a] = M;
    }
  }
  std::string name() const override { return "adjoint"; }

  DMat rep_matrix(const Cmat<N>& U) const override {
    const int n = n_gen<N>();
    const auto& gen = generators<N>();
    const Cmat<N> Ud = U.dagger();
    DMat M(n, n);
    for (int a = 0; a < n; ++a) {
      const Cmat<N> TaU = gen.T[a] * U;  // T^a U
      for (int b = 0; b < n; ++b) {
        const Cmat<N> tmp = TaU * gen.T[b] * Ud;  // T^a U T^b U^dag
        M(a, b) = Complex(2.0 * tmp.trace().real(), 0.0);
      }
    }
    return M;
  }
};

}  // namespace gh
