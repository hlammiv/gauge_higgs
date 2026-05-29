// Mean-field strong-weak coupling transition: validate the saddle-point criticals
// (beta_c, alpha_B, v_B) for the SU(2) sector against the lattice-notes table
// (temporal gauge), and cross-check the w_BT(alpha) closed form against the direct
// group sum and the w_SU(2) Bessel form against the maximal-torus quadrature.
//
//   g++ -std=c++20 -O3 -march=native -fopenmp -Isrc -o build/test_meanfield test/test_meanfield.cpp
//   ./build/test_meanfield
#include "check.hpp"
#define MEANFIELD_NO_MAIN
#include "meanfield.cpp"   // header-only style: pull in Ensemble + saddle machinery
#include <cmath>

using namespace gh;

// closed forms
static Real w_su2_bessel(Real a) { return std::log(2.0 * std::cyl_bessel_i(1, a) / a); }
static Real w_bt_closed(Real a)  { return std::log(3.0 + 8.0 * std::cosh(a / 2.0) + std::cosh(a)) - std::log(12.0); }

int main() {
  // ---- ensembles ----
  Ensemble su2 = ensemble_su2();
  Ensemble bt  = ensemble_discrete<2>(close_group<2>(gens_2T()));
  Ensemble bo  = ensemble_discrete<2>(close_group<2>(gens_2O()));
  Ensemble bi  = ensemble_discrete<2>(close_group<2>(gens_2I()));

  // ---- group orders ----
  CHECK(bt.retr.size() == 24,  "|2T| = 24");
  CHECK(bo.retr.size() == 48,  "|2O| = 48");
  CHECK(bi.retr.size() == 120, "|2I| = 120");

  // ---- w_BT closed form vs group sum ----
  for (Real a : {0.5, 1.3, 2.0, 3.4883, 5.0, 7.73}) {
    char m[80]; std::snprintf(m, sizeof m, "w_BT closed == group-sum at a=%.3f", a);
    CHECK_CLOSE(bt.w(a), w_bt_closed(a), 1e-12, m);
  }
  // ---- w_SU(2) Bessel vs torus quadrature ----
  for (Real a : {0.5, 1.3, 2.4340, 4.0, 4.719}) {
    char m[80]; std::snprintf(m, sizeof m, "w_SU(2) Bessel == torus-quad at a=%.3f", a);
    CHECK_CLOSE(su2.w(a), w_su2_bessel(a), 1e-6, m);
  }
  // ---- w'(0) = 0 (symmetric), w(0) = 0 ----
  CHECK_CLOSE(su2.w(1e-6), 0.0, 1e-10, "w_SU(2)(0) = 0");
  // v_B = w'(alpha) -> 0 linearly as alpha->0 (symmetric ensemble: <ReTr>=0).
  CHECK_CLOSE(bt.wp(1e-8), 0.0, 1e-8, "w'_BT(0) = 0 (v_B->0 at alpha->0)");

  // ---- target criticals (temporal gauge), tol ~1e-3 ----
  struct T { const char* name; Ensemble* e; int D; Real amax; Real beta, alphaB, vB; };
  T tab[] = {
    {"SU(2) D3", &su2, 3, 25.0, 1.9569, 2.4340, 0.4982},
    {"SU(2) D4", &su2, 4, 25.0, 1.6817, 4.7190, 0.7043},
    {"BT    D3", &bt,  3, 40.0, 1.9394, 3.4883, 0.6387},
    {"BT    D4", &bt,  4, 40.0, 1.5374, 7.7309, 0.9260},
    {"BO    D3", &bo,  3, 40.0, 1.9566, 2.4528, 0.5010},
    {"BO    D4", &bo,  4, 40.0, 1.6738, 5.0523, 0.7303},
    {"BI    D3", &bi,  3, 40.0, 1.9569, 2.4339, 0.4982},
    {"BI    D4", &bi,  4, 40.0, 1.6817, 4.7201, 0.7044},
  };
  std::printf("  %-9s  %-21s  %-21s  %-21s\n", "case", "beta_c", "alpha_B", "v_B");
  for (auto& t : tab) {
    PlaqGeom g(t.D);
    Critical c = find_critical(*t.e, g, t.amax);
    char m[96];
    std::snprintf(m, sizeof m, "%s found", t.name); CHECK(c.found, m);
    std::printf("  %-9s  %8.5f (ref %7.4f)  %9.5f (ref %7.4f)  %8.5f (ref %7.4f)\n",
                t.name, c.beta_c, t.beta, c.alpha_B, t.alphaB, c.vB, t.vB);
    std::snprintf(m, sizeof m, "%s beta_c",  t.name); CHECK_CLOSE(c.beta_c,  t.beta,   1e-3, m);
    std::snprintf(m, sizeof m, "%s alpha_B", t.name); CHECK_CLOSE(c.alpha_B, t.alphaB, 1e-3, m);
    std::snprintf(m, sizeof m, "%s v_B",     t.name); CHECK_CLOSE(c.vB,      t.vB,     1e-3, m);
  }

  return report("test_meanfield");
}
