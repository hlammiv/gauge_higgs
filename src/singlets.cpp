// H-singlet multiplicity enumeration by CHARACTERS (no full rep build).
//
//   m_{sigma->tr} = (1/|H|) sum_{h in H} chi^sigma(h)
//
// where chi^sigma is the SU(N) group character of irrep sigma evaluated on the
// finite-subgroup element h (an SU(N) matrix from close_group). This counts how
// many copies of the trivial H-irrep sit inside sigma|_H -- i.e. the number of
// independent H-singlet directions ("J"/"rho'" channels of the draft). It is the
// orthogonality-relation projector onto the trivial rep, so it is an exact
// non-negative integer (up to FP noise) and needs ONLY the character (trace data),
// never the d x d irrep matrices. This stays cheap even for (4,4) (d=125).
//
//   SU(2): chi^j(h) = sin((2j+1) theta)/sin(theta), e^{+-i theta} the eigenvalues of h.
//   SU(3): chi^(p,q)(h) via the Weyl bialternant in the 3 eigenvalues of h.
//
// Convention (matches finite_subgroups.hpp / rep_general.hpp):
//   SU(2) spin j  = Young row {2j}            (so j=3 d7, j=4 d9, j=6 d13).
//   SU(3) (p,q)   = Young rows [p+q, q].
//
// Cross-checks vs draft-1.tex Tab. su2subduction / su3subduction (su3_table.tex):
// red entries there = irreps containing the H-singlet; we reproduce their
// multiplicities. See docs/draft_validation_singlets.md.
#include "rep/finite_subgroups.hpp"
#include <array>
#include <cmath>
#include <complex>
#include <cstdio>
#include <vector>

namespace gh {

// ---------------------------------------------------------------------------
// SU(2) character of spin j on a 2x2 SU(2) element h.
//   eigenvalues e^{+-i theta} with 2 cos theta = Re Tr h (Tr h is real for SU(2)).
//   chi^j = sin((2j+1) theta)/sin(theta); the theta->0 (and theta->pi) limits are
//   handled by the integer limit chi^j(+-1) = (+-1)^{2j} (2j+1).
// ---------------------------------------------------------------------------
inline Real su2_character(int twoj, const Cmat<2>& h) {
  const Real half_tr = 0.5 * h.trace().real();         // = cos theta
  const Real c = std::max(-1.0, std::min(1.0, half_tr));
  const Real theta = std::acos(c);
  const Real s = std::sin(theta);
  if (std::fabs(s) < 1e-12) {                          // theta ~ 0 or pi: integer limit
    const Real sign = (c > 0.0) ? 1.0 : ((twoj % 2 == 0) ? 1.0 : -1.0);
    return sign * (twoj + 1);
  }
  return std::sin((twoj + 1) * theta) / s;
}

// ---------------------------------------------------------------------------
// SU(3) character of (p,q) on a 3x3 SU(3) element h, via the Jacobi-Trudi
// identity (no eigenvalue extraction, no division-by-Vandermonde -> numerically
// stable even when h has (near-)degenerate eigenvalues, e.g. central elements):
//
//   chi_lambda = det( h_{lambda_i - i + j} )_{i,j=1..3},
//
// with partition lambda = (p+q, q, 0) and h_k the complete homogeneous symmetric
// polynomial of degree k in the eigenvalues (h_k = 0 for k<0, h_0 = 1). The h_k are
// generated from the power sums  p_k = Tr(h^k)  via Newton's identity
//   h_k = (1/k) sum_{i=1..k} p_i h_{k-i}.
// Only traces of matrix powers are used -- exact for unitary h, no cancellation.
// (This equals the Weyl bialternant chi but is its determinant-of-h_k form.)
// ---------------------------------------------------------------------------
inline Complex su3_character(int pp, int qq, const Cmat<3>& h) {
  const std::array<int, 3> lam{pp + qq, qq, 0};
  const int kmax = lam[0] + 2;                       // largest index lambda_1 - 1 + 3
  // power sums p_k = Tr(h^k), k = 1..kmax.
  std::vector<Complex> psum(kmax + 1, Complex(0, 0));
  Cmat<3> hp = Cmat<3>::identity();
  for (int k = 1; k <= kmax; ++k) { hp = hp * h; psum[k] = hp.trace(); }
  // complete homogeneous symmetric polynomials h_k via Newton's identity.
  std::vector<Complex> hk(kmax + 1, Complex(0, 0));
  hk[0] = Complex(1, 0);
  for (int k = 1; k <= kmax; ++k) {
    Complex s(0, 0);
    for (int i = 1; i <= k; ++i) s += psum[i] * hk[k - i];
    hk[k] = s / static_cast<Real>(k);
  }
  auto H = [&](int k) -> Complex {                   // h_k with the h_k=0 (k<0) convention
    if (k < 0) return Complex(0, 0);
    return hk[k];
  };
  // 3x3 Jacobi-Trudi determinant det( h_{lambda_i - i + j} ).
  auto M = [&](int i, int j) { return H(lam[i] - i + j); };  // i,j = 0..2 (i.e. row i+1, col j+1)
  return M(0,0)*(M(1,1)*M(2,2) - M(1,2)*M(2,1))
       - M(0,1)*(M(1,0)*M(2,2) - M(1,2)*M(2,0))
       + M(0,2)*(M(1,0)*M(2,1) - M(1,1)*M(2,0));
}

// SU(2): multiplicity of the trivial H-irrep inside spin j.
inline int su2_singlet_mult(int twoj, const std::vector<Cmat<2>>& H) {
  Complex acc(0, 0);
  for (const auto& h : H) acc += su2_character(twoj, h);
  return static_cast<int>(std::llround(acc.real() / H.size()));
}
// SU(3): multiplicity of the trivial H-irrep inside (p,q).
inline int su3_singlet_mult(int p, int q, const std::vector<Cmat<3>>& H) {
  Complex acc(0, 0);
  for (const auto& h : H) acc += su3_character(p, q, h);
  return static_cast<int>(std::llround(acc.real() / H.size()));
}
// raw (un-rounded) value, for FP-integrality reporting.
inline Real su2_singlet_raw(int twoj, const std::vector<Cmat<2>>& H) {
  Complex acc(0, 0); for (const auto& h : H) acc += su2_character(twoj, h);
  return acc.real() / H.size();
}
inline Real su3_singlet_raw(int p, int q, const std::vector<Cmat<3>>& H) {
  Complex acc(0, 0); for (const auto& h : H) acc += su3_character(p, q, h);
  return acc.real() / H.size();
}

}  // namespace gh

#ifndef SINGLETS_NO_MAIN
using namespace gh;

int main() {
  std::printf("=== H-singlet multiplicity m_{sigma->tr} = (1/|H|) sum_h chi^sigma(h) ===\n");
  std::printf("    (characters only; no full rep build)\n\n");

  // ---- SU(2) ----
  auto T = close_group<2>(gens_2T());
  auto O = close_group<2>(gens_2O());
  auto I = close_group<2>(gens_2I());
  std::printf("-- SU(2) -> 2T(|H|=%d) 2O(|H|=%d) 2I(|H|=%d) --\n",
              (int)T.size(), (int)O.size(), (int)I.size());
  std::printf("  j  dim |  2T  2O  2I\n");
  std::printf("  ---------+----------\n");
  for (int twoj = 0; twoj <= 20; twoj += 1) {  // j = twoj/2 ; integer + half-integer
    const int mT = su2_singlet_mult(twoj, T);
    const int mO = su2_singlet_mult(twoj, O);
    const int mI = su2_singlet_mult(twoj, I);
    if (twoj % 2 == 0)
      std::printf("  %-2d %3d | %3d %3d %3d%s\n", twoj / 2, twoj + 1, mT, mO, mI,
                  (mT || mO || mI) ? "   <- singlet" : "");
    else
      std::printf("  %d/2 %2d | %3d %3d %3d\n", twoj, twoj + 1, mT, mO, mI);
  }

  // ---- SU(3) ----
  auto S108  = close_group<3>(gens_Sigma108());
  auto S216  = close_group<3>(gens_Sigma216());
  auto S648  = close_group<3>(gens_Sigma648());
  auto S1080 = close_group<3>(gens_Sigma1080());
  std::printf("\n-- SU(3) -> Sigma(108)(|H|=%d) Sigma(216)(|H|=%d) "
              "Sigma(648)(|H|=%d) Sigma(1080)(|H|=%d) --\n",
              (int)S108.size(), (int)S216.size(), (int)S648.size(), (int)S1080.size());
  std::printf("  (p,q)  dim |  S108 S216 S648 S1080\n");
  std::printf("  -----------+----------------------\n");
  auto dim3 = [](int p, int q) { return (p + 1) * (q + 1) * (p + q + 2) / 2; };
  for (int s = 0; s <= 8; ++s)
    for (int p = 0; p <= s; ++p) {
      const int q = s - p;
      const int m108 = su3_singlet_mult(p, q, S108);
      const int m216 = su3_singlet_mult(p, q, S216);
      const int m648 = su3_singlet_mult(p, q, S648);
      const int m1080 = su3_singlet_mult(p, q, S1080);
      const bool any = m108 || m216 || m648 || m1080;
      std::printf("  (%d,%d) %4d | %4d %4d %4d %4d%s\n", p, q, dim3(p, q),
                  m108, m216, m648, m1080, any ? "   <- singlet" : "");
    }
  std::printf("\n(rows '<- singlet' = draft red entries; m=2 = two singlet directions)\n");
  return 0;
}
#endif
