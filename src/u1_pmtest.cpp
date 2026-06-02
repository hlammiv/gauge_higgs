// Compact-U(1) + charge-q Higgs PHOTON-MASS validation driver.
//
// Measures the gauge-invariant photon mass from the STATIC magnetic structure
// factor at low non-zero spatial momentum (see u1/photon_structure.hpp), fitting
// the transverse photon pole R(p)^{-1} = (phat2 + m^2)/Z = a + b*phat2 and
// reporting m^2 = a/b +/- (delete-1 jackknife error over configs).
//
//   ./build/u1_pmtest [Ls n_scalar nmd ntherm nmeas lambda tau base_seed which_point Lt]
//
// GEOMETRY.  Ls is the SPATIAL extent and Lt the temporal extent; the lattice is
// Ls^(D-1) x Lt (direction 0 = Euclidean time). Lt defaults to Ls (isotropic) when
// omitted or <=0. The photon resolution is set by Ls (p_hat2_min=(2 sin(pi/Ls))^2),
// NOT Lt, so a geometric anisotropy (large Ls, modest Lt) gives the same photon
// resolution as the isotropic Ls^D at lower cost -- this is the lever for getting
// m^2(Higgs)>0 resolved (L=8 cannot; needs Ls>=16). Watch finite-T: Lt too small
// shifts the phase boundaries (check Lt=8 vs 12).
//
// With no arguments it runs the three hard-coded validation points:
//   (A) deep Coulomb : beta=1.3 kappa=0.0 q=1  -> expect m2 ~ 0 (massless photon)
//   (B) confined     : beta=0.7 kappa=0.0 q=1  -> expect m2 > 0 (massive)
//   (C) Higgs        : beta=1.3 kappa=0.8 q=2  -> expect m2 > 0 (massive)
// PASS: m2(A) small / consistent with 0, m2(B)>0, m2(C)>0, with m2(A) << both.
//
// The HMC setup mirrors u1_scan.cpp: hot()/cold_phi() init, autotune OFF,
// trajectory() per measurement. Compile-time D = NDIM (validation uses NDIM=4).

// monopole.hpp is NOT included here, so there is no in-TU clash with
// reduced_plaq_angle<D>; photon_mass.hpp (pulled in by photon_structure.hpp)
// provides the single definition. (If a future driver also includes monopole.hpp,
// it would #define GH_U1_HAVE_REDUCED_PLAQ_ANGLE before photon_mass.hpp, exactly
// as u1_scan.cpp does.)
#include "u1/u1.hpp"
#include "u1/photon_structure.hpp"
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <array>
#include <string>

using namespace gh;

namespace {

struct Point { const char* name; Real beta; Real kappa; int q; const char* expect; };

void run_point(const Point& pt, int Ls, int Lt, int n_scalar, int nmd, int ntherm, int nmeas,
               Real lambda, Real tau, std::uint64_t seed) {
  // Direction 0 = Euclidean time (extent Lt); directions 1..D-1 = spatial (Ls).
  // The photon momentum set (photon_structure.hpp) keys off lat.L[k+1]=Ls and the
  // structure factor sums over x_0 to enforce p_0=0, so Lt enters only the
  // statistics / finite-T, not the momentum resolution -- hence Ls is the lever.
  std::array<int, kDim> ext{};
  ext[0] = Lt;
  for (int mu = 1; mu < kDim; ++mu) ext[mu] = Ls;

  u1::U1HMC<kDim> hmc(ext, seed);
  hmc.beta = pt.beta; hmc.kappa = pt.kappa; hmc.lambda = lambda; hmc.q = pt.q;
  hmc.tau = tau; hmc.nmd = nmd; hmc.n_scalar = n_scalar;
  hmc.hot(0.8); hmc.cold_phi(0.5);

  for (int t = 0; t < ntherm; ++t) hmc.trajectory();
  hmc.traj_count = 0; hmc.accept_count = 0;

  // Momentum set (computed once for this lattice).
  const u1::PhotonMomenta<kDim> mom = u1::photon_momenta<kDim>(hmc.lat);

  // Accumulate per-config structure factors S(p) (one row per measurement).
  std::vector<std::vector<Real>> perConfig;
  perConfig.reserve(nmeas);
  for (int t = 0; t < nmeas; ++t) {
    hmc.trajectory();
    perConfig.push_back(u1::photon_structure_factor<kDim>(hmc.th, hmc.lat, mom));
  }

  const int ng = mom.n_groups();
  // PRIMARY fit: the lowest THREE distinct phat2 = the (1,0,0),(1,1,0),(1,1,1)
  // modes, which are genuinely BELOW the doubler. The (2,0,0) group at phat2=2 is
  // the doubler EDGE: its R flattens (lattice curvature), so including it bends the
  // R^{-1} line and corrupts the intercept. We therefore exclude it from the
  // headline pole fit and report the all-group fit only as a diagnostic.
  const int nfit_primary = (ng >= 3) ? 3 : ng;
  const u1::PhotonMassFit fit = u1::photon_mass_fit<kDim>(perConfig, mom, nfit_primary);
  const u1::PhotonMassFit fit_all = u1::photon_mass_fit<kDim>(perConfig, mom, ng);

  std::printf("\n==== %-12s  beta=%.3f kappa=%.3f q=%d  (expect: %s) ====\n",
              pt.name, pt.beta, pt.kappa, pt.q, pt.expect);
  std::printf("  geom=%d^%d x %d (Ls=%d spatial, Lt=%d time)  ntherm=%d nmeas=%d nmd=%d n_scalar=%d tau=%.3g lambda=%.3g  acc=%.3f\n",
              Ls, kDim - 1, Lt, Ls, Lt, ntherm, nmeas, nmd, n_scalar, tau, lambda, hmc.acceptance());
  std::printf("  momentum groups: %d (n_momenta=%d; primary fit uses lowest %d below doubler)\n",
              ng, mom.n_momenta(), nfit_primary);
  std::printf("  R(phat2) table  [R = <S(p)>/phat2 ~ Z/(phat2+m^2)]:\n");
  std::printf("    %-12s %-16s %-16s %-12s %s\n", "phat2", "R", "R_err", "1/R", "(in primary fit?)");
  for (int g = 0; g < ng; ++g) {
    const Real invR = (fit.R[g] != 0.0) ? 1.0 / fit.R[g] : 0.0;
    std::printf("    %-12.6f %-16.6g %-16.3g %-12.6g %s\n",
                fit.phat2[g], fit.R[g], fit.R_err[g], invR,
                (g < nfit_primary) ? "yes" : "no (doubler edge)");
  }
  std::printf("  ==> [PRIMARY, lowest %d] m2 = %.6g +/- %.4g   m_gamma = %.6g\n",
              nfit_primary, fit.m2, fit.m2_err, fit.m_gamma);
  std::printf("  ==> [diagnostic, all %d] m2 = %.6g +/- %.4g   m_gamma = %.6g\n",
              ng, fit_all.m2, fit_all.m2_err, fit_all.m_gamma);
  std::fflush(stdout);  // land per-point results immediately (robust to early reap)
}

}  // namespace

int main(int argc, char** argv) {
  std::setvbuf(stdout, nullptr, _IOLBF, 0);  // line-buffer stdout for live progress
  auto ai = [&](int i, long d) { return i < argc ? std::atol(argv[i]) : d; };
  auto af = [&](int i, double d) { return i < argc ? std::atof(argv[i]) : d; };

  const int  Ls       = static_cast<int>(ai(1, 8));
  const int  n_scalar = static_cast<int>(ai(2, 6));
  const int  nmd      = static_cast<int>(ai(3, 8));
  const int  ntherm   = static_cast<int>(ai(4, 200));
  const int  nmeas    = static_cast<int>(ai(5, 400));
  const Real lambda   = af(6, 0.5);
  const Real tau      = af(7, 1.0);
  const std::uint64_t base_seed = static_cast<std::uint64_t>(ai(8, 12345));
  // which_point: -1 (default) runs all three; 0/1/2 runs only A/B/C (same seed as
  // in the full run, so a single-point run reproduces that point exactly).
  const int which_point = static_cast<int>(ai(9, -1));
  // Lt: temporal extent. Sentinel <=0 (or omitted) -> isotropic Lt=Ls. A smaller Lt
  // gives a geometric anisotropy Ls^(D-1) x Lt (cheap photon resolution at Ls>=16).
  int Lt = static_cast<int>(ai(10, 0));
  if (Lt <= 0) Lt = Ls;

  std::printf("# U(1) photon-mass validation (static magnetic structure factor pole fit)\n");
  std::printf("# D=%d  geom=%d^%d x %d (Ls=%d, Lt=%d)  n_scalar=%d  nmd=%d  ntherm=%d  nmeas=%d  lambda=%.3g  tau=%.3g  autotune=OFF\n",
              kDim, Ls, kDim - 1, Lt, Ls, Lt, n_scalar, nmd, ntherm, nmeas, lambda, tau);

  const Point points[] = {
    {"A:Coulomb", 1.3, 0.0, 1, "m2 ~ 0 (massless)"},
    {"B:confined", 0.7, 0.0, 1, "m2 > 0 (massive)"},
    {"C:Higgs",   1.3, 0.8, 2, "m2 > 0 (massive)"},
  };

  for (int i = 0; i < 3; ++i)
    if (which_point < 0 || which_point == i)
      run_point(points[i], Ls, Lt, n_scalar, nmd, ntherm, nmeas, lambda, tau,
                Rng::key(base_seed, i));

  return 0;
}
