// Bounded-from-below (BFB) checker for the SU(3) quartic potential -- asserts:
//   (a) all-positive f_c  -> BFB  (each V_c >= 0; trivially copositive)
//   (b) negative-dominant f_c -> NOT BFB (an unbounded direction is found; min < 0)
//   (c) the Sigma(108) locking couplings (docs/locking_couplings.md) -> BFB
// Plus a sanity check that the SU(3) (2,2) channel C2 list matches CasimirChannels.
// Pure group theory / sphere minimization on a 27-dim complex vector. See src/bfb_su3.cpp.
#include "check.hpp"
#define BFB_NO_MAIN                 // pull in check_bfb / V4 without the driver main
#include "../src/bfb_su3.cpp"
#include <cstdio>

using namespace gh;

int main() {
  GeneralRep<3> rep({4, 2});        // (2,2), d = 27, Sigma(108) Higgs irrep
  CasimirChannels<3> ch(rep);
  const int nc = ch.n_channels();

  // sanity: (2,2) gives the 9 channels 0,3,6,8,12,15,18,20,24 (docs/locking_couplings.md)
  CHECK(nc == 9, "(2,2) has 9 R x Rbar channels");
  const Real expect[] = {0, 3, 6, 8, 12, 15, 18, 20, 24};
  bool chan_ok = (nc == 9);
  for (int c = 0; c < nc && chan_ok; ++c)
    if (std::fabs(ch.lambda[c] - expect[c]) > 1e-3) chan_ok = false;
  CHECK(chan_ok, "(2,2) channel C2 list = 0,3,6,8,12,15,18,20,24");

  // (a) all-positive f -> BFB. With ALL f_c = 1, V_4 = sum_c ||P_c M||^2 = ||M||_F^2 =
  // (Phi^dag Phi)^2 = 1 identically on the unit sphere (sum_c P_c = identity), so the
  // minimum must equal 1 exactly -- a clean analytic anchor for the checker.
  {
    std::vector<Real> f(nc, 1.0);
    auto res = check_bfb(ch, f, /*nrestart=*/120, /*nstep=*/120);
    std::printf("  (a) all f=1:            vmin = %+.6e  bfb=%d\n", res.vmin, (int)res.bfb);
    CHECK(res.bfb, "(a) all-positive f_c is BFB");
    CHECK(res.vmin > -1e-6, "(a) all-positive f_c: V_4 min >= 0");
    CHECK_CLOSE(res.vmin, 1.0, 1e-6, "(a) all f=1: V_4 == ||M||_F^2 == 1 on unit sphere");
  }

  // (b) negative-dominant f -> NOT BFB: a strongly negative coupling on one channel
  // makes V_4 negative along that channel's eigen-direction; min < 0 (unbounded ray).
  {
    std::vector<Real> f(nc, 0.05);
    f[nc - 1] = -1.0;
    auto res = check_bfb(ch, f, /*nrestart=*/120, /*nstep=*/150);
    std::printf("  (b) f_last=-1:          vmin = %+.6e  bfb=%d\n", res.vmin, (int)res.bfb);
    CHECK(!res.bfb, "(b) negative-dominant f_c is NOT BFB");
    CHECK(res.vmin < -1e-3, "(b) negative-dominant f_c: unbounded direction (V_4 min < 0)");
  }

  // (c) Sigma(108) locking couplings -> BFB. These lock SU(3) -> Sigma(108) (a stable PSD
  // minimum); a fortiori the quartic part must be bounded below.
  {
    std::vector<Real> f = {0.2782, 0.0930, 0.1029, 0.1256, 0.0071, 0.0101, 0.0443, 0.2148, 0.1240};
    auto res = check_bfb(ch, f, /*nrestart=*/150, /*nstep=*/150);
    std::printf("  (c) Sigma(108) locking: vmin = %+.6e  bfb=%d\n", res.vmin, (int)res.bfb);
    CHECK(res.bfb, "(c) Sigma(108) locking couplings are BFB");
    CHECK(res.vmin > -1e-6, "(c) Sigma(108) locking: V_4 min >= 0");
  }

  // cross-check: the locking f also come out of find_stable_couplings -> BFB.
  {
    auto H = close_group<3>(gens_Sigma108());
    auto sv = singlet_vev<3>(rep, H);
    auto stab = find_stable_couplings<3>(ch, sv.vev, rep.real, 600, 12345);
    auto res = check_bfb(ch, stab.f, /*nrestart=*/150, /*nstep=*/150);
    std::printf("  (c') find_stable_couplings: stable=%d  vmin = %+.6e  bfb=%d\n",
                (int)stab.stable, res.vmin, (int)res.bfb);
    CHECK(stab.stable, "(c') find_stable_couplings found a stable Sigma(108) point");
    CHECK(res.bfb, "(c') find_stable_couplings f_c are BFB");
  }

  return report("test_bfb");
}
