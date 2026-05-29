#pragma once
// Fundamental representation: D^(R)(U) = U, T^a_R = T^a, d = N, N-ality 1.
#include "rep/representation.hpp"

namespace gh {

template <int N>
DMat cmat_to_dmat(const Cmat<N>& M) {
  DMat D(N, N);
  for (int i = 0; i < N; ++i)
    for (int j = 0; j < N; ++j) D(i, j) = M(i, j);
  return D;
}

template <int N>
struct FundamentalRep : Representation<N> {
  FundamentalRep() {
    this->d = N; this->nality = 1; this->real = false;
    const auto& gen = generators<N>();
    this->T.resize(n_gen<N>());
    for (int a = 0; a < n_gen<N>(); ++a) this->T[a] = cmat_to_dmat<N>(gen.T[a]);
  }
  std::string name() const override { return "fundamental"; }
  DMat rep_matrix(const Cmat<N>& U) const override { return cmat_to_dmat<N>(U); }
};

}  // namespace gh
