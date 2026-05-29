#pragma once
// SU(N) group operations: matrix exponential exp(iH) (link drift), reunitarization,
// random group/algebra elements. See theory_notes §5.
#include "group/algebra.hpp"
#include <cmath>

namespace gh {

// exp(i H) for Hermitian H via scaling-and-squaring + Taylor. General N.
// (SU(2) closed form specialized below.) Result is unitary; if Tr H = 0 then det = 1.
template <int N>
Cmat<N> expi(const Cmat<N>& H) {
  const Real nrm = H.fnorm();
  int s = 0;
  Real scaled = nrm;
  while (scaled > 0.5) { scaled *= 0.5; ++s; }
  const Real inv2s = std::pow(0.5, s);
  // A = i H / 2^s
  Cmat<N> A = H * Complex(0.0, inv2s);
  Cmat<N> term = Cmat<N>::identity();
  Cmat<N> result = Cmat<N>::identity();
  for (int k = 1; k <= 18; ++k) {
    term = term * A;
    term *= Complex(1.0 / k, 0.0);
    result += term;
  }
  for (int i = 0; i < s; ++i) result = result * result;
  return result;
}

// SU(2) closed form: H = a·(sigma) traceless Hermitian, exp(iH) = cos r I + i sin(r)/r H,
// r = sqrt(Tr(H^2)/2). Exact, stays in SU(2).
template <>
inline Cmat<2> expi<2>(const Cmat<2>& H) {
  const Real r2 = 0.5 * trProd(H, H).real();   // Tr(H^2)/2 = a^2
  const Real r = std::sqrt(std::max(r2, 0.0));
  Real c, sinc;
  if (r < 1e-8) { c = 1.0 - 0.5 * r2; sinc = 1.0 - r2 / 6.0; }
  else          { c = std::cos(r);   sinc = std::sin(r) / r; }
  Cmat<2> U = Cmat<2>::identity();
  U *= Complex(c, 0.0);
  Cmat<2> iH = H * Complex(0.0, sinc);
  U += iH;
  return U;
}

// Reunitarize: modified Gram-Schmidt on rows, then rescale to det = 1.
template <int N>
void reunitarize(Cmat<N>& U) {
  for (int i = 0; i < N; ++i) {
    // Orthogonalize row i against previous rows.
    for (int k = 0; k < i; ++k) {
      Complex proj(0.0, 0.0);
      for (int j = 0; j < N; ++j) proj += std::conj(U(k, j)) * U(i, j);
      for (int j = 0; j < N; ++j) U(i, j) -= proj * U(k, j);
    }
    Real nrm = 0.0;
    for (int j = 0; j < N; ++j) nrm += std::norm(U(i, j));
    nrm = std::sqrt(nrm);
    const Complex inv(1.0 / nrm, 0.0);
    for (int j = 0; j < N; ++j) U(i, j) *= inv;
  }
  // Fix determinant phase to 1: divide row 0 by det.
  Complex d = det(U);
  const Real ad = std::abs(d);
  if (ad > 0.0) {
    const Complex phase = d / ad;          // unit-modulus determinant
    const Complex inv = std::conj(phase);  // 1/phase
    for (int j = 0; j < N; ++j) U(0, j) *= inv;
  }
}

// ||U^dag U - I|| Frobenius, unitarity diagnostic.
template <int N>
Real unitarity_violation(const Cmat<N>& U) {
  Cmat<N> d = U.dagger() * U;
  for (int i = 0; i < N; ++i) d(i, i) -= Complex(1.0, 0.0);
  return d.fnorm();
}

}  // namespace gh
