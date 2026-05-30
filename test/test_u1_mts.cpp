// Tests for the nested multi-timescale (Sexton-Weingarten) MD integrator of the
// compact U(1) + charge-q Higgs HMC. The fast/stiff timescale carries the whole
// scalar sector (scalar_force on pi + matter back-reaction add_matter_force on p);
// the slow timescale carries only the smooth gauge plaquette force.
//   1. n_scalar=1 reproduces the single-timescale Omelyan trajectory exactly.
//   2. Reversibility of the nested integrator at a stiff point (q=2,3).
//   3. <exp(-dH)> ~ 1 and good acceptance at a stiff point.
//   4. Cost benefit: at equal gauge-force evaluations, multi-timescale acceptance
//      is clearly higher than single-timescale.
#include "check.hpp"
#include "u1/u1.hpp"

using namespace gh;
using namespace gh::u1;

template <int D>
static std::array<int, D> cube(int n) { std::array<int, D> L{}; for (int mu = 0; mu < D; ++mu) L[mu] = n; return L; }

// 1. n_scalar=1 must reproduce the single-timescale md_evolve to machine precision.
template <int D>
static void test_nscalar1_exact(int q) {
  // single-timescale reference
  U1HMC<D> a(cube<D>(4), 1234);
  a.beta = 1.0; a.kappa = 0.6; a.lambda = 0.5; a.q = q; a.tau = 1.0; a.nmd = 8; a.n_scalar = 1;
  a.hot(0.5); a.refresh_momenta();
  auto th0 = a.th; auto phi0 = a.phi; auto p0 = a.p; auto pi0 = a.pi;
  a.md_evolve();

  // multi-timescale path with n_scalar=1, identical start
  U1HMC<D> b(cube<D>(4), 1234);
  b.beta = 1.0; b.kappa = 0.6; b.lambda = 0.5; b.q = q; b.tau = 1.0; b.nmd = 8; b.n_scalar = 1;
  b.th = th0; b.phi = phi0; b.p = p0; b.pi = pi0;
  b.md_evolve_mts();

  Real w = 0;
  for (std::size_t i = 0; i < th0.size(); ++i) w = std::max(w, std::fabs(b.th[i] - a.th[i]));
  for (std::size_t i = 0; i < phi0.size(); ++i) w = std::max(w, std::abs(b.phi[i] - a.phi[i]));
  char m[80]; std::snprintf(m, sizeof m, "MTS n_scalar=1 == single-timescale (q=%d)", q);
  CHECK_CLOSE(w, 0.0, 1e-12, m);
}

// 2. Reversibility of the nested integrator at a stiff point.
template <int D>
static void test_reversibility_mts(int q) {
  U1HMC<D> hmc(cube<D>(4), 99);
  hmc.beta = 1.0; hmc.kappa = 0.6; hmc.lambda = 0.5; hmc.q = q; hmc.tau = 1.0; hmc.nmd = 8; hmc.n_scalar = 3;
  hmc.hot(0.5); hmc.refresh_momenta();
  auto th0 = hmc.th; auto phi0 = hmc.phi;
  hmc.md_evolve_mts();
  for (auto& v : hmc.p) v = -v;
  for (auto& z : hmc.pi) z = -z;
  hmc.md_evolve_mts();
  Real w = 0;
  for (std::size_t i = 0; i < th0.size(); ++i) w = std::max(w, std::fabs(hmc.th[i] - th0[i]));
  for (std::size_t i = 0; i < phi0.size(); ++i) w = std::max(w, std::abs(hmc.phi[i] - phi0[i]));
  char m[64]; std::snprintf(m, sizeof m, "MTS reversibility n_scalar=3 (q=%d)", q);
  CHECK_CLOSE(w, 0.0, 1e-9, m);
}

// 3. <exp(-dH)> ~ 1 and good acceptance at a stiff point with the nested integrator.
template <int D>
static void test_expdH_mts(int q, std::uint64_t seed) {
  U1HMC<D> hmc(cube<D>(4), seed);
  hmc.beta = 1.0; hmc.kappa = 0.6; hmc.lambda = 0.5; hmc.q = q; hmc.tau = 1.0; hmc.nmd = 8; hmc.n_scalar = 4;
  hmc.hot(0.3); hmc.cold_phi(0.5);
  const int ntraj = 300; double s = 0;
  for (int t = 0; t < ntraj; ++t) { hmc.trajectory(); s += std::exp(-hmc.last_dH); }
  const double mean = s / ntraj;
  char m[110];
  std::snprintf(m, sizeof m, "MTS <exp(-dH)>~1 (q=%d) got %.4f acc %.2f", q, mean, hmc.acceptance());
  CHECK(std::fabs(mean - 1.0) < 0.12, m);
  std::snprintf(m, sizeof m, "MTS acceptance > 0.6 (q=%d) acc %.2f", q, hmc.acceptance());
  CHECK(hmc.acceptance() > 0.6, m);
}

// 4. Cost benefit: SAME nmd (same number of slow gauge-force evals per trajectory),
//    multi-timescale acceptance must be clearly higher than single-timescale.
//    Demonstrated at a genuinely stiff point kappa=0.9 where the single-timescale
//    integrator collapses (the spec's kappa=0.6 is only mildly stiff on L=4: the
//    single-timescale path already sits near acc ~0.95 there, so the gap is small).
template <int D>
static void test_cost_benefit(int q, std::uint64_t seed) {
  const int ntraj = 300;
  const Real kappa = 0.9;
  // single-timescale
  U1HMC<D> s1(cube<D>(4), seed);
  s1.beta = 1.0; s1.kappa = kappa; s1.lambda = 0.5; s1.q = q; s1.tau = 1.0; s1.nmd = 8; s1.n_scalar = 1;
  s1.hot(0.3); s1.cold_phi(0.5);
  for (int t = 0; t < ntraj; ++t) s1.trajectory();
  // multi-timescale, SAME nmd (so the slow gauge force is evaluated the same number of times)
  U1HMC<D> s6(cube<D>(4), seed);
  s6.beta = 1.0; s6.kappa = kappa; s6.lambda = 0.5; s6.q = q; s6.tau = 1.0; s6.nmd = 8; s6.n_scalar = 6;
  s6.hot(0.3); s6.cold_phi(0.5);
  for (int t = 0; t < ntraj; ++t) s6.trajectory();

  const double g1 = double(s1.kick_count_gauge) / ntraj;   // gauge-force evals / traj
  const double g6 = double(s6.kick_count_gauge) / ntraj;
  std::printf("  [cost] q=%d kappa=%.1f nmd=8  single n_scalar=1: acc %.3f, gauge evals/traj %.1f\n",
              q, kappa, s1.acceptance(), g1);
  std::printf("  [cost] q=%d kappa=%.1f nmd=8  multi  n_scalar=6: acc %.3f, gauge evals/traj %.1f\n",
              q, kappa, s6.acceptance(), g6);
  char m[110];
  std::snprintf(m, sizeof m, "MTS equal gauge-evals/traj (q=%d) %.1f vs %.1f", q, g1, g6);
  CHECK_CLOSE(g6, g1, 1e-9, m);
  std::snprintf(m, sizeof m, "MTS acceptance higher than single (q=%d): %.3f > %.3f",
                q, s6.acceptance(), s1.acceptance());
  CHECK(s6.acceptance() > s1.acceptance() + 0.05, m);
}

int main() {
  std::printf("-- n_scalar=1 reproduces single timescale --\n");
  test_nscalar1_exact<3>(2); test_nscalar1_exact<3>(3); test_nscalar1_exact<4>(2);
  std::printf("-- nested integrator reversibility (stiff) --\n");
  test_reversibility_mts<3>(2); test_reversibility_mts<3>(3);
  std::printf("-- <exp(-dH)> / acceptance (stiff, nested) --\n");
  test_expdH_mts<3>(2, 41);
  std::printf("-- cost benefit at equal gauge-force evaluations --\n");
  test_cost_benefit<3>(2, 71);
  return report("test_u1_mts");
}
