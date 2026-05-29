// Tests for the static-potential / Creutz-ratio post-processor (src/measure/creutz.hpp).
// PURE analytic post-processing -- no lattice, no RNG, no simulation. We feed synthetic
// W(R,T) grids with KNOWN closed forms and pin the exact cancellations to machine
// precision (area & perimeter laws), plus a qualitative monotonicity test for the
// Coulombic case.
#include "check.hpp"
#include "measure/creutz.hpp"
#include <vector>
#include <cmath>
#include <cstdio>
#include <limits>

using namespace gh;

static constexpr int kMax = 6;  // R,T range 1..6

// Build a (kMax+1)x(kMax+1) grid with W[R][T] = f(R,T); row/col 0 left as placeholder.
template <class F>
static std::vector<std::vector<Real>> make_grid(F f) {
  std::vector<std::vector<Real>> W(kMax + 1, std::vector<Real>(kMax + 1, 0.0));
  for (int R = 1; R <= kMax; ++R)
    for (int T = 1; T <= kMax; ++T) W[R][T] = f(R, T);
  return W;
}

// ---------------------------------------------------------------- 1. pure area law
// W = exp(-sigma R T - p(R+T) - c). Perimeter & constant cancel EXACTLY in chi, which
// must equal sigma to machine precision for ALL valid (R,T). The static-potential
// R-slope V(R+1,T)-V(R,T) = -sigma exactly and is independent of T.
static void test_area_law() {
  const Real sigma = 0.2, p = 0.1, c = 0.05;
  auto W = make_grid([&](int R, int T) {
    return std::exp(-sigma * R * T - p * (R + T) - c);
  });

  // chi(R,T) == sigma to 1e-12 for every R,T >= 2 (T+? not needed: chi uses only <=R,<=T).
  for (int R = 2; R <= kMax; ++R)
    for (int T = 2; T <= kMax; ++T) {
      char m[80]; std::snprintf(m, sizeof m, "area-law chi(%d,%d) == sigma", R, T);
      CHECK_CLOSE(creutz_ratio(W, R, T), sigma, 1e-12, m);
    }

  // static_potential closed form: V(R,T) = -(sigma R + p).
  for (int R = 1; R <= kMax; ++R)
    for (int T = 1; T < kMax; ++T) {  // needs W[R][T+1]
      char m[80]; std::snprintf(m, sizeof m, "area-law V(%d,%d) == -(sigma R + p)", R, T);
      CHECK_CLOSE(static_potential(W, R, T), -(sigma * R + p), 1e-12, m);
    }

  // R-slope of V recovers -sigma exactly and is T-independent.
  for (int R = 1; R < kMax; ++R)
    for (int T = 1; T < kMax; ++T) {
      const Real slope = static_potential(W, R + 1, T) - static_potential(W, R, T);
      char m[88]; std::snprintf(m, sizeof m, "area-law V-slope(%d,%d) == -sigma", R, T);
      CHECK_CLOSE(slope, -sigma, 1e-12, m);
    }

  // plateau string tension and recovered perimeter coefficient.
  CHECK_CLOSE(string_tension_plateau(W), sigma, 1e-12, "area-law plateau sigma");
  CHECK_CLOSE(perimeter_coeff(W, sigma), p, 1e-12, "area-law perimeter coeff p");
}

// ---------------------------------------------------------------- 2. pure perimeter law
// W = exp(-p(R+T) - c). No area term -> chi == 0 exactly; V(R) slope in R -> 0 exactly.
static void test_perimeter_law() {
  const Real p = 0.1, c = 0.05;
  auto W = make_grid([&](int R, int T) { return std::exp(-p * (R + T) - c); });

  for (int R = 2; R <= kMax; ++R)
    for (int T = 2; T <= kMax; ++T) {
      char m[80]; std::snprintf(m, sizeof m, "perimeter-law chi(%d,%d) == 0", R, T);
      CHECK_CLOSE(creutz_ratio(W, R, T), 0.0, 1e-12, m);
    }

  // V(R,T) = -p (constant in R), so the R-slope vanishes exactly.
  for (int R = 1; R < kMax; ++R)
    for (int T = 1; T < kMax; ++T) {
      const Real slope = static_potential(W, R + 1, T) - static_potential(W, R, T);
      char m[88]; std::snprintf(m, sizeof m, "perimeter-law V-slope(%d,%d) == 0", R, T);
      CHECK_CLOSE(slope, 0.0, 1e-12, m);
    }
  CHECK_CLOSE(string_tension_plateau(W), 0.0, 1e-12, "perimeter-law plateau sigma == 0");
}

// ---------------------------------------------------------------- 3. Coulombic law
// W = exp(-(a/R) T - p(R+T)). The perimeter piece cancels exactly in the 2x2 corner
// combination; the Coulomb piece survives and evaluates to the closed form
//   chi(R,T) = -a [ 1/(R-1) - 1/R ] = -a / (R (R-1)),  independent of T.
// With the FIXED sign convention chi = -log(...) the area law gives +sigma (confining,
// positive), so a Coulomb (screened, non-confining) potential necessarily gives the
// OPPOSITE sign: chi < 0, with a magnitude |chi| = a/(R(R-1)) that is strictly DECREASING
// in R -- the string-tension estimate falls off to zero, the hallmark of no confinement.
// We assert sign, the exact closed form (to 1e-12), and the strict monotone fall-off of
// the magnitude in R at fixed T.
static void test_coulomb_law() {
  const Real a = 0.5, p = 0.1;
  auto W = make_grid([&](int R, int T) {
    return std::exp(-(a / R) * T - p * (R + T));
  });

  for (int T = 2; T <= kMax; ++T) {
    Real prevMag = std::numeric_limits<Real>::infinity();
    for (int R = 2; R <= kMax; ++R) {
      const Real chi = creutz_ratio(W, R, T);
      char m[88];
      // Sign: non-confining Coulomb piece is negative under chi = -log(...).
      std::snprintf(m, sizeof m, "coulomb chi(%d,%d) < 0", R, T);
      CHECK(chi < 0.0, m);
      // Exact closed form -a/(R(R-1)), T-independent.
      std::snprintf(m, sizeof m, "coulomb chi(%d,%d) == -a/(R(R-1))", R, T);
      CHECK_CLOSE(chi, -a / (Real(R) * (R - 1)), 1e-12, m);
      // Magnitude (string-tension estimate) strictly decreasing in R.
      const Real mag = std::fabs(chi);
      std::snprintf(m, sizeof m, "coulomb |chi| decreasing in R at T=%d (R=%d)", T, R);
      CHECK(mag < prevMag, m);
      prevMag = mag;
    }
  }
}

// ---------------------------------------------------------------- 4. bounds / guards
// chi requires R,T>=2; accessors must return NaN (not crash) on out-of-range / R,T<2.
static void test_bounds() {
  auto W = make_grid([&](int R, int T) { return std::exp(-0.2 * R * T); });
  CHECK(std::isnan(creutz_ratio(W, 1, 3)), "chi(R<2) -> NaN");
  CHECK(std::isnan(creutz_ratio(W, 3, 1)), "chi(T<2) -> NaN");
  CHECK(std::isnan(creutz_ratio(W, 99, 3)), "chi(R out of range) -> NaN");
  CHECK(std::isnan(static_potential(W, kMax, kMax)), "V at top T (no T+1) -> NaN");
  CHECK(std::isnan(static_potential(W, 99, 1)), "V(R out of range) -> NaN");
  // A valid in-range query is finite.
  CHECK(std::isfinite(creutz_ratio(W, 3, 3)), "chi(3,3) finite");
  CHECK(std::isfinite(static_potential(W, 3, 3)), "V(3,3) finite");
}

int main() {
  std::printf("-- area law (exact) --\n");      test_area_law();
  std::printf("-- perimeter law (exact) --\n"); test_perimeter_law();
  std::printf("-- coulomb law (qualitative) --\n"); test_coulomb_law();
  std::printf("-- bounds / guards --\n");        test_bounds();
  return report("test_creutz");
}
