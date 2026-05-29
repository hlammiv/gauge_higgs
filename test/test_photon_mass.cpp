// Tests for the photon-mass / effective-mass / plateau module:
//   1. cosh & log effective-mass math exactness on analytic inputs
//   2. plateau() jackknife: identical samples -> err 0, scattered -> matches an
//      independent brute-force delete-1 jackknife reference
//   3. field-strength gauge invariance of the reduced plaquette angles, the
//      timeslice sums S_i(t), and the correlator C(dt)
//   4. cold sanity: theta=0 -> all field strengths 0 -> correlator all 0
// NO HMC; everything is a tiny constructed/cold/random config or synthetic data.
#include "check.hpp"
#include "measure/correlator.hpp"
#include "u1/photon_mass.hpp"
#include <random>
#include <vector>
#include <cmath>

using namespace gh;
using namespace gh::u1;

template <int D>
static std::array<int, D> cube(int n) { std::array<int, D> L{}; for (int mu = 0; mu < D; ++mu) L[mu] = n; return L; }

// ---- 1. effective-mass math exactness ----
static void test_effmass_math() {
  const int T = 16; const Real m = 0.7;
  std::vector<Real> C(T);
  for (int t = 0; t < T; ++t) C[t] = std::cosh(m * (t - T / 2.0));
  const std::vector<Real> me = cosh_effective_mass(C);
  // Recover m in the BULK. cosh(m(t-T/2)) is symmetric about T/2 but not strictly
  // periodic, so its argument sign flips across the periodic seam t=0 (neighbors
  // t=T-1 and t=1 sit on the same side, giving arg<1 there -- clamped to m_eff=0).
  // Every interior t=1..T-1 recovers m to machine precision.
  Real worst = 0.0;
  for (int t = 1; t < T; ++t) worst = std::max(worst, std::fabs(me[t] - m));
  CHECK_CLOSE(worst, 0.0, 1e-13, "cosh_eff recovers m=0.7 to machine precision (bulk)");
  CHECK_CLOSE(me[0], 0.0, 1e-12, "cosh_eff at seam t=0 clamps to 0");

  // Constant correlator -> arg=1 -> m_eff=0 to 1e-12.
  std::vector<Real> Cc(T, 3.3);
  const std::vector<Real> mec = cosh_effective_mass(Cc);
  Real worstc = 0.0;
  for (int t = 0; t < T; ++t) worstc = std::max(worstc, std::fabs(mec[t]));
  CHECK_CLOSE(worstc, 0.0, 1e-12, "cosh_eff(const)=0 to 1e-12");

  // log effective mass on a pure decaying exponential C(t)=exp(-m t) (open form):
  //   m(t)=log(C[t]/C[t+1])=m for t=0..T-2 (the wrap t=T-1 is not a clean single exp).
  std::vector<Real> Ce(T);
  for (int t = 0; t < T; ++t) Ce[t] = std::exp(-m * t);
  const std::vector<Real> mel = log_effective_mass(Ce);
  Real worstl = 0.0;
  for (int t = 0; t < T - 1; ++t) worstl = std::max(worstl, std::fabs(mel[t] - m));
  CHECK_CLOSE(worstl, 0.0, 1e-12, "log_eff recovers m on pure exp");
}

// Independent brute-force delete-1 jackknife reference for the windowed cosh
// effective mass (re-implements the estimator a different way as the hand check).
static PlateauFit ref_plateau(const std::vector<std::vector<Real>>& samples, int tmin, int tmax) {
  const int K = (int)samples.size();
  const int T = (int)samples[0].size();
  auto windowed = [&](const std::vector<Real>& C) {
    Real s = 0; int n = 0;
    for (int t = tmin; t <= tmax; ++t) {
      const Real arg = (C[(t - 1 + T) % T] + C[(t + 1) % T]) / (2.0 * C[t]);
      s += std::acosh(arg); ++n;
    }
    return s / n;
  };
  if (K == 1) return {windowed(samples[0]), 0.0};
  std::vector<Real> th(K);
  for (int k = 0; k < K; ++k) {
    std::vector<Real> Cjk(T, 0.0);
    for (int j = 0; j < K; ++j) if (j != k) for (int t = 0; t < T; ++t) Cjk[t] += samples[j][t];
    for (int t = 0; t < T; ++t) Cjk[t] /= (K - 1);
    th[k] = windowed(Cjk);
  }
  Real tbar = 0; for (Real v : th) tbar += v; tbar /= K;
  Real sw = 0; for (Real v : th) sw += (v - tbar) * (v - tbar);
  return {tbar, std::sqrt((Real)(K - 1) / K * sw)};
}

// ---- 2. plateau jackknife ----
static void test_plateau() {
  const int T = 16; const Real m = 0.55; const int tmin = 4, tmax = 12;

  // K identical known-mass cosh samples -> mass==m exactly, err==0.
  std::vector<Real> base(T);
  for (int t = 0; t < T; ++t) base[t] = std::cosh(m * (t - T / 2.0));
  std::vector<std::vector<Real>> ident(8, base);
  PlateauFit pid = plateau(ident, tmin, tmax);
  CHECK_CLOSE(pid.mass, m, 1e-10, "plateau identical cosh -> mass==m");
  CHECK_CLOSE(pid.err, 0.0, 1e-14, "plateau identical cosh -> err==0");

  // K=1 -> err 0 and mass == single-config windowed effective mass (== m here).
  PlateauFit p1 = plateau({base}, tmin, tmax);
  CHECK_CLOSE(p1.mass, m, 1e-10, "plateau K=1 -> single-config mass");
  CHECK_CLOSE(p1.err, 0.0, 1e-14, "plateau K=1 -> err==0");

  // Add small KNOWN scatter across samples: amplitude AND mass jitter so that the
  // delete-1 replicas genuinely differ; compare to the independent reference.
  std::mt19937_64 rng(12345); std::normal_distribution<Real> g(0, 1);
  const int K = 10;
  std::vector<std::vector<Real>> scat(K, std::vector<Real>(T));
  for (int k = 0; k < K; ++k) {
    const Real Ak = 1.0 + 0.05 * g(rng);    // amplitude scatter
    const Real mk = m + 0.03 * g(rng);      // mass scatter
    for (int t = 0; t < T; ++t) scat[k][t] = Ak * std::cosh(mk * (t - T / 2.0));
  }
  PlateauFit pj = plateau(scat, tmin, tmax);
  PlateauFit pr = ref_plateau(scat, tmin, tmax);
  CHECK_CLOSE(pj.mass, pr.mass, 1e-12, "plateau scatter mass == reference");
  CHECK_CLOSE(pj.err, pr.err, 1e-12, "plateau scatter err == hand jackknife");
  CHECK(pj.err > 1e-6, "plateau scatter err is genuinely nonzero");
}

// ---- 3. field-strength gauge invariance ----
template <int D>
static void test_gauge_invariance(std::uint64_t seed) {
  Lattice<D> lat(cube<D>(4));
  std::mt19937_64 rng(seed); std::normal_distribution<Real> g(0, 1);
  std::vector<Real> th(lat.vol * D); for (auto& t : th) t = 0.8 * g(rng);

  // Reduced plaquette angles before/after a random gauge transform.
  std::vector<Real> a(lat.vol); for (auto& x : a) x = g(rng);
  std::vector<Real> th2 = th;
  for (std::int64_t x = 0; x < lat.vol; ++x)
    for (int mu = 0; mu < D; ++mu) th2[x * D + mu] += a[x] - a[lat.neighbor_fwd(x, mu)];

  // (a) reduced field strengths invariant site-by-site, plane-by-plane.
  Real wF = 0.0;
  for (std::int64_t s = 0; s < lat.vol; ++s)
    for (int mu = 0; mu < D; ++mu)
      for (int nu = mu + 1; nu < D; ++nu)
        wF = std::max(wF, std::fabs(reduced_plaq_angle<D>(th, lat, s, mu, nu)
                                  - reduced_plaq_angle<D>(th2, lat, s, mu, nu)));
  CHECK_CLOSE(wF, 0.0, 1e-10, "reduced field strength gauge invariant");

  // (b) timeslice sums S_i(t) invariant.
  auto S1 = photon_timeslice_field<D>(th, lat);
  auto S2 = photon_timeslice_field<D>(th2, lat);
  Real wS = 0.0;
  for (std::size_t t = 0; t < S1.size(); ++t)
    for (std::size_t i = 0; i < S1[t].size(); ++i)
      wS = std::max(wS, std::fabs(S1[t][i] - S2[t][i]));
  CHECK_CLOSE(wS, 0.0, 1e-10, "timeslice S_i(t) gauge invariant");

  // (c) correlator C(dt) invariant.
  auto C1 = photon_correlator({S1});
  auto C2 = photon_correlator({S2});
  Real wC = 0.0;
  for (std::size_t dt = 0; dt < C1.size(); ++dt) wC = std::max(wC, std::fabs(C1[dt] - C2[dt]));
  CHECK_CLOSE(wC, 0.0, 1e-10, "correlator C(dt) gauge invariant");
}

// ---- 4. cold sanity: theta=0 -> F=0 -> S=0 -> C=0 ----
template <int D>
static void test_cold() {
  Lattice<D> lat(cube<D>(4));
  std::vector<Real> th(lat.vol * D, 0.0);
  // field strengths all zero
  Real wF = 0.0;
  for (std::int64_t s = 0; s < lat.vol; ++s)
    for (int mu = 0; mu < D; ++mu)
      for (int nu = mu + 1; nu < D; ++nu)
        wF = std::max(wF, std::fabs(reduced_plaq_angle<D>(th, lat, s, mu, nu)));
  CHECK_CLOSE(wF, 0.0, 1e-13, "cold reduced field strength = 0");

  auto S = photon_timeslice_field<D>(th, lat);
  Real wS = 0.0; for (auto& sl : S) for (Real v : sl) wS = std::max(wS, std::fabs(v));
  CHECK_CLOSE(wS, 0.0, 1e-13, "cold timeslice S_i(t) = 0");

  auto C = photon_correlator({S});
  Real wC = 0.0; for (Real v : C) wC = std::max(wC, std::fabs(v));
  CHECK_CLOSE(wC, 0.0, 1e-13, "cold correlator C(dt) = 0");
}

// Bonus: a synthetic monopole flux 2*pi in one plaquette is correctly subtracted
// by the principal value (physical F unchanged), proving F is the field strength
// mod the Dirac string -- the gauge-invariant content.
static void test_monopole_reduction() {
  // raw angle that is 0.3 + 2*pi should reduce to 0.3
  CHECK_CLOSE(reduce_angle(0.3 + 2.0 * kPi), 0.3, 1e-13, "reduce_angle strips +2pi");
  CHECK_CLOSE(reduce_angle(-0.4 - 4.0 * kPi), -0.4, 1e-13, "reduce_angle strips -4pi");
  CHECK(reduce_angle(kPi - 1e-9) > 0.0 && reduce_angle(kPi - 1e-9) < kPi, "reduce_angle in [-pi,pi)");
  CHECK_CLOSE(reduce_angle(kPi), -kPi, 1e-12, "reduce_angle(+pi) -> -pi (half-away boundary; measure-zero, physics-irrelevant)");
}

int main() {
  std::printf("-- effective-mass math --\n"); test_effmass_math();
  std::printf("-- plateau jackknife --\n");   test_plateau();
  std::printf("-- field-strength gauge invariance --\n");
  test_gauge_invariance<4>(101); test_gauge_invariance<3>(102);
  std::printf("-- cold sanity --\n"); test_cold<4>(); test_cold<3>();
  std::printf("-- monopole/principal-value --\n"); test_monopole_reduction();
  return report("test_photon_mass");
}
