#pragma once
// su(N) Lie algebra: generalized Gell-Mann generators (Hermitian, Tr T^aT^b=1/2 d),
// algebra <-> matrix maps, TA projection, structure constants. See theory_notes §6.1.
#include "core/linalg.hpp"
#include <array>

namespace gh {

template <int N>
constexpr int n_gen() { return N * N - 1; }

// Algebra element as real components in the generalized Gell-Mann basis.
template <int N>
using AlgVec = std::array<Real, static_cast<std::size_t>(N * N - 1)>;

// The N^2-1 Hermitian generators T^a = Lambda^a / 2, ordered:
//   symmetric off-diagonal pairs, antisymmetric off-diagonal pairs, diagonal.
template <int N>
struct Generators {
  std::array<Cmat<N>, static_cast<std::size_t>(N * N - 1)> T{};
};

template <int N>
Generators<N> build_generators() {
  Generators<N> g{};
  int a = 0;
  // Symmetric: |j><k| + |k><j|, normalized so Tr(Lambda^2)=2 -> Lambda; T=Lambda/2.
  for (int j = 0; j < N; ++j)
    for (int k = j + 1; k < N; ++k) {
      Cmat<N>& T = g.T[a++];
      T(j, k) += Complex(0.5, 0.0);  // (1/2)(|j><k| + |k><j|)
      T(k, j) += Complex(0.5, 0.0);
    }
  // Antisymmetric: -i|j><k| + i|k><j|.
  for (int j = 0; j < N; ++j)
    for (int k = j + 1; k < N; ++k) {
      Cmat<N>& T = g.T[a++];
      T(j, k) += Complex(0.0, -0.5);
      T(k, j) += Complex(0.0, 0.5);
    }
  // Diagonal: l=1..N-1, Lambda^l = sqrt(2/(l(l+1))) (sum_{j<=l-1}|j><j| - l|l><l|) (0-indexed).
  for (int l = 1; l <= N - 1; ++l) {
    Cmat<N>& T = g.T[a++];
    const Real c = std::sqrt(2.0 / (static_cast<Real>(l) * (l + 1)));
    for (int j = 0; j < l; ++j) T(j, j) += Complex(0.5 * c, 0.0);
    T(l, l) += Complex(-0.5 * c * l, 0.0);
  }
  return g;
}

// Cached generator set (built once per N).
template <int N>
const Generators<N>& generators() {
  static const Generators<N> g = build_generators<N>();
  return g;
}

// Reconstruct the matrix M = sum_a v[a] T^a (Hermitian, traceless).
template <int N>
Cmat<N> alg_to_mat(const AlgVec<N>& v) {
  const auto& g = generators<N>();
  Cmat<N> M{};
  for (int a = 0; a < n_gen<N>(); ++a) {
    const Real va = v[a];
    if (va == 0.0) continue;
    const Cmat<N>& T = g.T[a];
    for (int i = 0; i < N; ++i)
      for (int j = 0; j < N; ++j) M(i, j) += va * T(i, j);
  }
  return M;
}

// Project a matrix onto algebra components: v[a] = 2 Re Tr(T^a M).
// For Hermitian traceless M = sum_b v[b] T^b this inverts alg_to_mat exactly.
template <int N>
AlgVec<N> mat_to_alg(const Cmat<N>& M) {
  const auto& g = generators<N>();
  AlgVec<N> v{};
  for (int a = 0; a < n_gen<N>(); ++a) v[a] = 2.0 * reTrProd(g.T[a], M);
  return v;
}

// Traceless anti-Hermitian projection: [M]_TA = 1/2(M - M^dag) - (1/2N)Tr(M - M^dag) I.
template <int N>
Cmat<N> proj_TA(const Cmat<N>& M) {
  Cmat<N> A = M.antiherm();           // 1/2 (M - M^dag)
  Complex tr = A.trace();             // = 1/2 Tr(M - M^dag)
  Complex shift = tr * Complex(1.0 / N, 0.0);
  for (int i = 0; i < N; ++i) A(i, i) -= shift;
  return A;
}

// AlgVec arithmetic helpers.
template <int N>
Real dotAlg(const AlgVec<N>& x, const AlgVec<N>& y) {
  Real s = 0.0;
  for (int a = 0; a < n_gen<N>(); ++a) s += x[a] * y[a];
  return s;
}

// Structure constants f^{abc} = -2 i Tr([T^a,T^b] T^c) (real). Cached.
template <int N>
const std::vector<Real>& structure_constants() {
  static const std::vector<Real> f = [] {
    const auto& g = generators<N>();
    const int n = n_gen<N>();
    std::vector<Real> f(static_cast<std::size_t>(n) * n * n, 0.0);
    for (int a = 0; a < n; ++a)
      for (int b = 0; b < n; ++b) {
        Cmat<N> comm = g.T[a] * g.T[b] - g.T[b] * g.T[a];  // [T^a,T^b]
        for (int c = 0; c < n; ++c) {
          Complex t = trProd(comm, g.T[c]);  // Tr([T^a,T^b] T^c) = i*y (purely imaginary)
          // f^{abc} = -2i Tr(..) = -2i(iy) = 2y = 2 Im(Tr(..)).
          f[(static_cast<std::size_t>(a) * n + b) * n + c] = 2.0 * t.imag();
        }
      }
    return f;
  }();
  return f;
}

}  // namespace gh
