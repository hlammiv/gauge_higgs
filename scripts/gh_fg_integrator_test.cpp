// Adversarial correctness gate for the ForceGradient2MNFG integrator branch.
// Checks: (1) combined reversibility of md_evolve under the FG branch (unfrozen + frozen),
//         (2) <exp(-dH)>~1 and acceptance under the FG branch.
#include "rep/rep_fundamental.hpp"
#include "rep/rep_general.hpp"
#include "action/scalar_invariants.hpp"
#include "hmc/gauge_higgs_hmc.hpp"
#include <cstdio>
#include <cmath>
using namespace gh;

template <int D>
static std::array<int, D> cube(int n) { std::array<int, D> L{}; for (int mu = 0; mu < D; ++mu) L[mu] = n; return L; }

template <int D, int N>
static int reversibility(Representation<N>& rep, std::uint64_t seed, bool frozen, const char* tag) {
  GaugeHiggsHMC<D, N> hmc(cube<D>(4), rep, seed);
  hmc.beta = 1.5; hmc.kappa = 1.0; hmc.lambda = 0.5; hmc.tau = 1.0; hmc.nmd = 8;
  hmc.integ = Integrator::ForceGradient2MNFG;
  hmc.reunit_each_traj = false;
  hmc.frozen_phi = frozen;
  hmc.U.hot(hmc.rng, 0.5);
  hmc.phi.gaussian(hmc.rng, 1, rep.real, 0.6);
  if (frozen) hmc.normalize_phi();
  hmc.refresh_momenta();
  std::vector<Cmat<N>> U0 = hmc.U.u;
  std::vector<Complex> phi0 = hmc.phi.data;
  hmc.md_evolve();
  for (auto& v : hmc.P.p) for (int a = 0; a < n_gen<N>(); ++a) v[a] = -v[a];
  for (auto& z : hmc.pi.data) z = -z;
  if (frozen) hmc.project_pi_tangent();
  hmc.md_evolve();
  Real worst = 0.0;
  for (std::size_t i = 0; i < U0.size(); ++i) worst = std::max(worst, (hmc.U.u[i] - U0[i]).fnorm());
  for (std::size_t i = 0; i < phi0.size(); ++i) worst = std::max(worst, std::abs(hmc.phi.data[i] - phi0[i]));
  bool ok = worst < 1e-9;
  std::printf("[%s] reversibility%s worst=%.2e %s\n", tag, frozen ? " (frozen)" : "", worst, ok ? "PASS" : "FAIL");
  return ok ? 0 : 1;
}

template <int D, int N>
static int expdh(Representation<N>& rep, std::uint64_t seed, bool frozen, int nmd, const char* tag) {
  GaugeHiggsHMC<D, N> hmc(cube<D>(4), rep, seed);
  hmc.beta = 1.5; hmc.kappa = 1.0; hmc.lambda = 0.5; hmc.tau = 1.0; hmc.nmd = nmd;
  hmc.integ = Integrator::ForceGradient2MNFG;
  hmc.frozen_phi = frozen;
  hmc.U.hot(hmc.rng, 0.3);
  hmc.phi.cold(0.5);
  if (frozen) hmc.normalize_phi();
  const int ntraj = 400; double s = 0.0;
  for (int t = 0; t < ntraj; ++t) { hmc.trajectory(); s += std::exp(-hmc.last_dH); }
  const double mean = s / ntraj;
  bool ok = std::fabs(mean - 1.0) < 0.15;
  std::printf("[%s] <exp(-dH)>%s=%.4f acc=%.2f nmd=%d %s\n", tag, frozen ? " (frozen)" : "", mean,
              hmc.acceptance(), nmd, ok ? "PASS" : "FAIL");
  return ok ? 0 : 1;
}

int main() {
  int fail = 0;
  FundamentalRep<2> f2;
  GeneralRep<2> g6({6});  // spin-3 (2T) = BT campaign rep
  std::printf("-- FG reversibility (must be machine-precision exact) --\n");
  fail += reversibility<3, 2>(f2, 401, false, "SU(2) fund");
  fail += reversibility<3, 2>(f2, 402, true,  "SU(2) fund");
  fail += reversibility<3, 2>(g6, 403, false, "spin-3 {6}");
  fail += reversibility<3, 2>(g6, 404, true,  "spin-3 {6}");
  std::printf("-- FG <exp(-dH)>~1 / acceptance (detailed balance) --\n");
  fail += expdh<3, 2>(f2, 501, false, 20, "SU(2) fund");
  fail += expdh<3, 2>(f2, 502, true,  20, "SU(2) fund");
  fail += expdh<3, 2>(g6, 503, true,  20, "spin-3 {6}");
  std::printf("== %s ==\n", fail ? "FAILURES PRESENT" : "ALL FG CORRECTNESS CHECKS PASS");
  return fail ? 1 : 0;
}
