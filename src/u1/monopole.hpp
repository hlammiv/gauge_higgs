#pragma once
// DeGrand-Toussaint monopole density for compact U(1), templated for D=3 and D=4.
//
// Reference: T. DeGrand & D. Toussaint, "Topological excitations and Monte Carlo
// simulation of Abelian gauge theory," Phys. Rev. D 22, 2478 (1980).
//
// Idea.  Each plaquette angle theta_P = theta_{mu nu}(s) (from plaq_angle<D>) is split
// into a "physical" reduced angle theta_bar in [-pi,pi) and an integer number of 2*pi
// wraps -- the Dirac string n_{mu nu}(s):
//
//     theta_P = theta_bar + 2*pi * n_{mu nu}(s),   n_{mu nu}(s) = round(theta_P / 2*pi).
//
// (theta_bar is what is left after subtracting the wraps; theta_bar = theta_P - 2*pi*n.)
// The conserved magnetic charge is the lattice curl of n, i.e. the net Dirac string
// flux threading the surface of an elementary 3-cube.  Because n_{mu nu} is the round of
// a plaquette angle it is automatically antisymmetric, n_{nu mu}(s) = -n_{mu nu}(s)
// (plaq_angle<D>(s,nu,mu) = -plaq_angle<D>(s,mu,nu)), which is what makes the discrete
// curl below well defined and which guarantees current conservation (d^2 n = 0).
//
//   D=3:  one elementary 3-cube per site (dual site).  The monopole charge m(s) is the
//         oriented sum of n over the cube's 6 faces:
//             m(s) = sum_{(i,j,l) cyclic of (0,1,2)} [ n_{jl}(s+i^) - n_{jl}(s) ].
//         density = < sum_s |m(s)| > / V.
//
//   D=4:  there are D=4 distinct 3-cubes touching each site, each obtained by omitting
//         one direction rho.  The cube dual to rho carries the monopole CURRENT k_rho(s),
//         again the oriented 6-face sum of n over that 3-cube (the three directions
//         != rho).  density = < sum_{s,rho} |k_rho(s)| > / (4*V).
//
// Normalization.  In D=3 we divide the summed |m| by V (one cube per site).  In D=4 we
// divide the summed |k_rho| by 4*V = (number of directions) * (number of sites) -- i.e.
// the number of dual 3-cubes -- so the density is the mean |current| per dual cube and
// is directly comparable across D.
//
// Caveat.  DeGrand-Toussaint counts are NOISY (UV-dominated, lattice-artifact heavy)
// below beta ~ 1; treat the density as a qualitative diagnostic of the (de)confining /
// Coulomb regimes, NOT as a strict, scaling order parameter.
#include "u1/u1.hpp"        // plaq_angle<D> (read-only reuse)
#include "core/geometry.hpp"
#include "core/config.hpp"
#include <vector>
#include <cmath>
#include <array>
#include <cstdint>

namespace gh {
namespace u1 {

// Dirac string integer n_{mu nu}(s) = round(plaq_angle / 2*pi): the number of 2*pi wraps
// in the plaquette angle, equivalently (theta_P - theta_bar)/2*pi with theta_bar in
// [-pi,pi).  Antisymmetric in (mu,nu) because plaq_angle is.
template <int D>
int dirac_string(const std::vector<Real>& th, const Lattice<D>& lat,
                 std::int64_t s, int mu, int nu) {
  const Real tp = plaq_angle<D>(th, lat, s, mu, nu);
  return static_cast<int>(std::lround(tp / (2.0 * kPi)));
}

// Reduced ("physical") plaquette angle theta_bar in [-pi,pi): theta_P - 2*pi*n.
template <int D>
Real reduced_plaq_angle(const std::vector<Real>& th, const Lattice<D>& lat,
                        std::int64_t s, int mu, int nu) {
  const Real tp = plaq_angle<D>(th, lat, s, mu, nu);
  return tp - 2.0 * kPi * static_cast<Real>(std::lround(tp / (2.0 * kPi)));
}

// Oriented sum of the Dirac string n over the 6 faces of the elementary 3-cube based at
// site s and spanned by the three distinct directions (i,j,l).  This is the lattice curl
//   sum_cyclic [ n_{jl}(s + i^) - n_{jl}(s) ]
// taken over the cyclic permutations of (i,j,l).  For (i,j,l) a positively oriented cyclic
// triple this returns the magnetic charge / current piercing the cube.  It is an exact
// integer (sum of rounds).
template <int D>
int cube_charge(const std::vector<Real>& th, const Lattice<D>& lat,
                std::int64_t s, int i, int j, int l) {
  int m = 0;
  // three cyclic orientations (i;j,l), (j;l,i), (l;i,j)
  const int a[3] = {i, j, l};
  for (int c = 0; c < 3; ++c) {
    const int di = a[c], dj = a[(c + 1) % 3], dl = a[(c + 2) % 3];
    const std::int64_t sp = lat.neighbor_fwd(s, di);
    m += dirac_string<D>(th, lat, sp, dj, dl) - dirac_string<D>(th, lat, s, dj, dl);
  }
  return m;
}

// D=3: signed monopole charge in the single elementary 3-cube at site s (directions 0,1,2).
template <int D>
int monopole_charge_D3(const std::vector<Real>& th, const Lattice<D>& lat, std::int64_t s) {
  static_assert(D == 3, "monopole_charge_D3 is only defined for D=3");
  return cube_charge<D>(th, lat, s, 0, 1, 2);
}

// Levi-Civita symbol in 4D.
inline int eps4(int a, int b, int c, int d) {
  int p[4] = {a, b, c, d}, sgn = 1;
  for (int i = 0; i < 4; ++i)
    for (int j = i + 1; j < 4; ++j) {
      if (p[i] == p[j]) return 0;
      if (p[i] > p[j]) { int t = p[i]; p[i] = p[j]; p[j] = t; sgn = -sgn; }
    }
  return sgn;
}

// D=4: signed monopole current k_rho(s) through the 3-cube at site s dual to direction rho.
// It is the Hodge dual of the lattice curl of the Dirac string n (a 2-form), built with a
// BACKWARD lattice difference so that the current is conserved in the conventional backward
// form sum_rho [ k_rho(s) - k_rho(s - rho^) ] == 0:
//     k_rho(s) = 1/2 * sum_{sigma,mu,nu} eps_{rho sigma mu nu} [ n_{mu nu}(s) - n_{mu nu}(s-sigma^) ].
// The full antisymmetric eps sum (with the 1/2 absorbing the mu<->nu double count of the
// antisymmetric n) is what makes the dual current EXACTLY conserved (d^2 n = 0): the discrete
// divergence vanishes identically at every dual site.  The result is an integer (the (mu,nu)
// and (nu,mu) terms are equal, so the eps sum is even and the 1/2 is exact).  The choice of
// backward (vs forward) difference only relabels which dual site carries the current; it does
// not change |k_rho| (hence the density) and is fixed here to match the backward divergence.
template <int D>
int monopole_current_D4(const std::vector<Real>& th, const Lattice<D>& lat,
                        std::int64_t s, int rho) {
  static_assert(D == 4, "monopole_current_D4 is only defined for D=4");
  int acc = 0;
  for (int sigma = 0; sigma < 4; ++sigma) {
    if (sigma == rho) continue;
    const std::int64_t sm = lat.neighbor_bwd(s, sigma);
    for (int mu = 0; mu < 4; ++mu)
      for (int nu = 0; nu < 4; ++nu) {
        const int e = eps4(rho, sigma, mu, nu);
        if (e == 0) continue;
        acc += e * (dirac_string<D>(th, lat, s, mu, nu) - dirac_string<D>(th, lat, sm, mu, nu));
      }
  }
  return acc / 2;  // exact: (mu,nu) and (nu,mu) contribute equally
}

// Monopole density:
//   D=3: < sum_s |m(s)| > / V.
//   D=4: < sum_{s,rho} |k_rho(s)| > / (4*V).
template <int D>
Real monopole_density(const std::vector<Real>& th, const Lattice<D>& lat) {
  static_assert(D == 3 || D == 4, "monopole_density implemented for D=3 and D=4");
  std::int64_t total = 0;  // sum of |charge|/|current| (a non-negative integer)
  if constexpr (D == 3) {
    #pragma omp parallel for schedule(static) reduction(+:total)
    for (std::int64_t s = 0; s < lat.vol; ++s)
      total += std::abs(cube_charge<D>(th, lat, s, 0, 1, 2));
    return static_cast<Real>(total) / static_cast<Real>(lat.vol);
  } else {  // D == 4
    #pragma omp parallel for schedule(static) reduction(+:total)
    for (std::int64_t s = 0; s < lat.vol; ++s)
      for (int rho = 0; rho < 4; ++rho)
        total += std::abs(monopole_current_D4<D>(th, lat, s, rho));
    return static_cast<Real>(total) / (4.0 * static_cast<Real>(lat.vol));
  }
}

}  // namespace u1
}  // namespace gh
