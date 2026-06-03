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

// --- Fundamental-link log: w with exp(i w.T_F) = U (inverse of expi(alg_to_mat(w))) ---
// Used by the GeneralRep fast path D^(R)(U) = exp(i w.T_R). Any valid branch of the log
// works: the rep is a homomorphism, so exp(i w.T_R) reproduces D^(R)(U) regardless of
// which w (mod 2pi shifts) maps to U in the fundamental.

// Small complex N x N inverse (Gauss-Jordan, partial pivot). General N.
template <int N>
Cmat<N> inverse(Cmat<N> A) {
  Cmat<N> Inv = Cmat<N>::identity();
  for (int col = 0; col < N; ++col) {
    int piv = col; Real best = std::abs(A(col, col));
    for (int r = col + 1; r < N; ++r) { Real v = std::abs(A(r, col)); if (v > best) { best = v; piv = r; } }
    if (piv != col)
      for (int j = 0; j < N; ++j) { std::swap(A(col, j), A(piv, j)); std::swap(Inv(col, j), Inv(piv, j)); }
    const Complex invd = Complex(1.0, 0.0) / A(col, col);
    for (int j = 0; j < N; ++j) { A(col, j) *= invd; Inv(col, j) *= invd; }
    for (int r = 0; r < N; ++r) if (r != col) {
      const Complex f = A(r, col);
      if (f == Complex(0, 0)) continue;
      for (int j = 0; j < N; ++j) { A(r, j) -= f * A(col, j); Inv(r, j) -= f * Inv(col, j); }
    }
  }
  return Inv;
}

// Principal matrix square root via Denman-Beavers (small N; quadratic convergence).
template <int N>
Cmat<N> sqrt_mat(const Cmat<N>& A) {
  Cmat<N> Y = A, Z = Cmat<N>::identity();
  for (int it = 0; it < 60; ++it) {
    Cmat<N> Yi = inverse<N>(Y), Zi = inverse<N>(Z);
    Cmat<N> Yn = (Y + Zi) * Complex(0.5, 0.0);
    Cmat<N> Zn = (Z + Yi) * Complex(0.5, 0.0);
    Y = Yn; Z = Zn;
    Cmat<N> r = Y * Y - A;
    if (r.fnorm() < 1e-14) break;
  }
  return Y;
}

// General N: X = -i log U via inverse-scaling-squaring + Mercator series; w = mat_to_alg(X).
template <int N>
AlgVec<N> fund_alg(const Cmat<N>& U) {
  const Cmat<N> I = Cmat<N>::identity();
  Cmat<N> V = U; int s = 0;
  while ((V - I).fnorm() > 0.3 && s < 50) { V = sqrt_mat<N>(V); ++s; }
  const Cmat<N> E = V - I;                       // log V = log(I+E) = sum_k (-1)^{k+1} E^k / k
  Cmat<N> Ek = E, logV = E;
  for (int k = 2; k <= 40; ++k) {
    Ek = Ek * E;
    logV += Ek * Complex(((k % 2) ? 1.0 : -1.0) / k, 0.0);
  }
  const Cmat<N> logU = logV * Complex(std::pow(2.0, s), 0.0);   // log U = 2^s log(U^{1/2^s})
  const Cmat<N> X = (logU * Complex(0.0, -1.0)).herm();         // -i log U (Hermitian); herm = denoise
  return mat_to_alg<N>(X);
}

// SU(2) closed form (exact, robust): U = cos r I + i (sin r / r) H, H = w.T traceless Hermitian.
template <>
inline AlgVec<2> fund_alg<2>(const Cmat<2>& U) {
  const Real c = 0.5 * U.trace().real();                       // = cos r
  const Real cc = std::max(-1.0, std::min(1.0, c));
  const Real r = std::acos(cc);
  const Real sinc = (r < 1e-8) ? 1.0 : std::sin(r) / r;        // sin r / r
  const Real invsinc = (std::fabs(sinc) > 1e-12) ? 1.0 / sinc : 0.0;
  // antiherm(U) = (U - U^dag)/2 = i sinc H  =>  H = -i antiherm(U) / sinc.
  const Cmat<2> H = U.antiherm() * Complex(0.0, -invsinc);
  return mat_to_alg<2>(H);
}

}  // namespace gh
