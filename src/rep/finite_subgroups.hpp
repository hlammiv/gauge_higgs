#pragma once
// Fundamental-rep generators for the crystal-like finite subgroups used in the
// SU(N)->discrete-H breaking program, plus a brute-force group-closure routine.
// Generators verified (unitary, det=+1, closure orders) from Grimus-Ludl
// arXiv:1006.0098 (SU(3) Sigma series) and the standard binary-polyhedral / icosian
// construction (SU(2) 2T,2O,2I). See [[research-program]].
//   SU(2): 2T (24), 2O (48), 2I (120)        SU(3): Sigma(108),(216),(648),(1080)
#include "core/linalg.hpp"
#include <vector>

namespace gh {

// Brute-force closure: multiply by generators (BFS from identity) until saturation.
template <int N>
std::vector<Cmat<N>> close_group(const std::vector<Cmat<N>>& gens, int cap = 4000) {
  std::vector<Cmat<N>> G;
  auto known = [&](const Cmat<N>& M) {
    for (const auto& g : G) if ((g - M).fnorm() < 1e-7) return true;
    return false;
  };
  G.push_back(Cmat<N>::identity());
  std::vector<Cmat<N>> frontier = G;
  while (!frontier.empty()) {
    std::vector<Cmat<N>> next;
    for (const auto& g : frontier)
      for (const auto& s : gens) {
        Cmat<N> p = g * s;
        if (!known(p)) { G.push_back(p); next.push_back(p); if (static_cast<int>(G.size()) > cap) return G; }
      }
    frontier.swap(next);
  }
  return G;
}

namespace detail {
inline Cmat<2> mat2(Complex a, Complex b, Complex c, Complex d) { Cmat<2> M; M(0,0)=a; M(0,1)=b; M(1,0)=c; M(1,1)=d; return M; }
inline Cmat<3> mat3(std::initializer_list<Complex> e) { Cmat<3> M; int i=0; for (Complex z : e) M(i/3, i%3) = z, ++i; return M; }
const Complex I_(0,1);
const Complex w_(-0.5, 0.8660254037844387);          // e^{2pi i/3}
const Complex w2_(-0.5, -0.8660254037844387);         // e^{4pi i/3}
const double  r2_ = 0.7071067811865476;               // 1/sqrt2
const double  r3_ = 0.5773502691896258;               // 1/sqrt3
}  // namespace detail

// ---- SU(2) binary polyhedral ----
inline std::vector<Cmat<2>> gens_2T() {
  using namespace detail;
  Cmat<2> gi = mat2(I_, 0, 0, -I_);
  Cmat<2> gT = mat2(Complex(0.5,0.5), Complex(0.5,0.5), Complex(-0.5,0.5), Complex(0.5,-0.5));
  return {gi, gT};
}
inline std::vector<Cmat<2>> gens_2O() {
  using namespace detail;
  auto g = gens_2T();
  g.push_back(mat2(Complex(r2_,r2_), 0, 0, Complex(r2_,-r2_)));   // (1+i)/sqrt2
  return g;
}
inline std::vector<Cmat<2>> gens_2I() {
  using namespace detail;
  Cmat<2> gi = mat2(I_, 0, 0, -I_);
  // s = (phi + phi^{-1} i + j)/2 ; full precision so closure stays exactly in SU(2)
  const double ph2 = 0.8090169943749475;   // phi/2
  const double pi2 = 0.3090169943749475;   // phi^{-1}/2
  Cmat<2> s = mat2(Complex(ph2,pi2), 0.5, -0.5, Complex(ph2,-pi2));
  return {gi, s};
}
// Binary dihedral / quaternion group Q8 = preimage in SU(2) of the Klein four-group
// V4 = {1, R_x(pi), R_y(pi), R_z(pi)} subset SO(3). Generators i*sigma_1, i*sigma_2
// close to {+-I, +-i sigma_1, +-i sigma_2, +-i sigma_3} (order 8). Smallest crystal-like
// SU(2) subgroup -- the "real control" beneath 2T: a single REAL spin-2 (l=2) Higgs on
// the biaxial stratum (3 distinct quadrupole eigenvalues) breaks SU(2)->Q8; the uniaxial
// stratum instead leaves a continuous O(2). See [[discrete-observables-program]].
inline std::vector<Cmat<2>> gens_Q8() {
  using namespace detail;
  Cmat<2> a = mat2(0, I_, I_, 0);    // i*sigma_1 = [[0,i],[i,0]]
  Cmat<2> b = mat2(0, 1, -1, 0);     // i*sigma_2 = [[0,1],[-1,0]]
  return {a, b};
}

// ---- SU(3) Sigma series (Delta(27) blocks E,C + extra generators) ----
inline Cmat<3> su3_E() { using namespace detail; return mat3({0,1,0, 0,0,1, 1,0,0}); }
inline Cmat<3> su3_C() { using namespace detail; return mat3({1,0,0, 0,w_,0, 0,0,w2_}); }
inline Cmat<3> su3_V() {  // (1/(i sqrt3)) [[1,1,1],[1,w,w^2],[1,w^2,w]] ; prefactor -i/sqrt3 -> det=+1
  using namespace detail; Complex p(0, -r3_);
  return mat3({p*Complex(1,0), p*Complex(1,0), p*Complex(1,0),
               p*Complex(1,0), p*w_,          p*w2_,
               p*Complex(1,0), p*w2_,         p*w_});
}
inline Cmat<3> su3_X() {  // (1/(i sqrt3)) [[1,1,w^2],[1,w,w],[w,1,w]]
  using namespace detail; Complex p(0, -r3_);
  return mat3({p*Complex(1,0), p*Complex(1,0), p*w2_,
               p*Complex(1,0), p*w_,          p*w_,
               p*w_,           p*Complex(1,0), p*w_});
}
inline Cmat<3> su3_D() {  // diag(eps, eps, eps*w), eps=e^{4pi i/9}
  using namespace detail; Complex eps(0.17364817766693041, 0.984807753012208);
  return mat3({eps,0,0, 0,eps,0, 0,0,eps*w_});
}
inline Cmat<3> su3_A2() { using namespace detail; return mat3({1,0,0, 0,-1,0, 0,0,-1}); }
inline Cmat<3> su3_W() {  // (1/2)[[-1,mu2,mu1],[mu2,mu1,-1],[mu1,-1,mu2]], mu1=0.618034, mu2=-1.618034
  using namespace detail; double m1=0.618033988749895, m2=-1.618033988749895;
  return mat3({Complex(-0.5,0), Complex(0.5*m2,0), Complex(0.5*m1,0),
               Complex(0.5*m2,0), Complex(0.5*m1,0), Complex(-0.5,0),
               Complex(0.5*m1,0), Complex(-0.5,0), Complex(0.5*m2,0)});
}
inline Cmat<3> su3_F() { using namespace detail; return mat3({-1,0,0, 0,0,-w_, 0,-w2_,0}); }

inline std::vector<Cmat<3>> gens_Sigma108()  { return {su3_C(), su3_E(), su3_V()}; }
inline std::vector<Cmat<3>> gens_Sigma216()  { return {su3_C(), su3_E(), su3_V(), su3_X()}; }
inline std::vector<Cmat<3>> gens_Sigma648()  { return {su3_C(), su3_E(), su3_V(), su3_D()}; }
inline std::vector<Cmat<3>> gens_Sigma1080() { return {su3_A2(), su3_E(), su3_W(), su3_F()}; }

}  // namespace gh
