// Tests for multi-coupling histogram reweighting (src/measure/reweight.hpp).
//
// PURE analytic post-processing -- no lattice, no HMC, no simulation. Every input
// is synthetic data drawn deterministically (seeded std::mt19937_64). We pin the
// estimator against KNOWN closed forms:
//
//   The Monte Carlo weight is exp(-S), S = sum_i lambda_i E_i, so reweighting a
//   reference ensemble at lambda0 to lambda tilts each sample by exp(-dlambda.E).
//   For Gaussian energies this is the exact Gaussian-tilting (Esscher) identity:
//     E ~ N(mu, Sigma)  =>  under exp(-dlambda.E) the tilted law is N(mu - Sigma.dlambda, Sigma).
//   Hence the reweighted means are exactly  <E>_tilted = mu - Sigma . dlambda,
//   and (1-D) <E^2>_tilted = (mu - dlambda sigma^2)^2 + sigma^2.
//
// Tests:
//   1. Single-histogram exactness, 1 coupling: E ~ N(mu,sigma^2), O=E and O=E^2.
//      Recover the tilted mean / second moment to ~ sigma/sqrt(N) MC tolerance.
//   2. 2-coupling reweighting: (E1,E2) jointly Gaussian; recover mu - Sigma.dlambda.
//   3. Multi-histogram WHAM: 3 reference ensembles at different lambda0 reweighted
//      to a common target agree with each other and with single-histogram within
//      errors; the f_k iteration converges (residual < tol, iters < cap).
//   4. Log-sum-exp stability: |dlambda*E| ~ 700+, finite & correct vs a tiny
//      exactly-computable two-sample reference; no inf/nan. Plus log_sum_exp/
//      weighted_mean_exp unit identities.
#include "check.hpp"
#include "measure/reweight.hpp"
#include <random>
#include <vector>
#include <array>
#include <cmath>
#include <cstdio>

using namespace gh;

// ---------------------------------------------------------------------------- 1
// Single-histogram exactness, ONE coupling. E ~ N(mu, sigma^2) drawn with a fixed
// seed. Reweight by dlambda. Gaussian tilting => <E>_tilt = mu - dlambda*sigma^2,
// <E^2>_tilt = (mu - dlambda*sigma^2)^2 + sigma^2. We assert recovery to a few MC
// sigmas. With N samples the statistical error on the tilted mean is ~ sigma/sqrt(N)
// (modestly inflated by the reweighting variance for the dlambda used here).
static void test_single_1coupling() {
  const Real mu = 1.3, sigma = 0.8;
  const int N = 4000000;
  std::mt19937_64 rng(0xC0FFEEu);
  std::normal_distribution<Real> gN(mu, sigma);

  std::vector<std::vector<Real>> E(1, std::vector<Real>(N));
  std::vector<Real> O(N), O2(N);
  for (int t = 0; t < N; ++t) {
    const Real e = gN(rng);
    E[0][t] = e;
    O[t]  = e;        // O = E
    O2[t] = e * e;    // O = E^2
  }

  const std::vector<Real> lambda0{0.0};  // reference: untilted
  // Several tilts including both signs; keep dlambda*sigma modest so the
  // reweighting variance stays controlled at this N.
  for (Real dl : {-0.6, -0.2, 0.0, 0.3, 0.7}) {
    const std::vector<Real> lambda{dl};
    const Real meanE  = single_histogram(E, O,  lambda0, lambda);
    const Real meanE2 = single_histogram(E, O2, lambda0, lambda);

    const Real exactE  = mu - dl * sigma * sigma;
    const Real exactE2 = exactE * exactE + sigma * sigma;

    // MC tolerance: base sigma/sqrt(N) inflated by the effective-sample-size loss
    // from reweighting, exp((dl*sigma)^2)/sqrt(N), times a safety factor. For the
    // untilted dl=0 case this is exact to ~ sigma/sqrt(N).
    const Real ess = std::exp((dl * sigma) * (dl * sigma));
    const Real tolE  = 12.0 * sigma * ess / std::sqrt((Real)N);
    const Real tolE2 = 12.0 * (sigma * sigma + 2.0 * std::fabs(exactE) * sigma) * ess
                         / std::sqrt((Real)N);

    char m[96];
    std::snprintf(m, sizeof m, "1-coupling <E>_tilt(dl=%.2f) == mu - dl*sig^2", dl);
    CHECK_CLOSE(meanE, exactE, tolE, m);
    std::snprintf(m, sizeof m, "1-coupling <E^2>_tilt(dl=%.2f) closed form", dl);
    CHECK_CLOSE(meanE2, exactE2, tolE2, m);
  }

  // dl=0 must reproduce the plain sample mean to MACHINE precision (no tilt at all).
  Real sm = 0.0; for (Real v : O) sm += v; sm /= N;
  const Real m0 = single_histogram(E, O, lambda0, lambda0);
  CHECK_CLOSE(m0, sm, 1e-10, "dl=0 reweight == plain sample mean (machine)");
}

// ---------------------------------------------------------------------------- 2
// Two-coupling reweighting. (E1,E2) jointly Gaussian with known mean mu and
// covariance Sigma (built from a Cholesky factor applied to iid normals). Tilting
// by dlambda=(dl1,dl2) gives the EXACT tilted means  <E>_tilt = mu - Sigma.dlambda.
static void test_two_coupling() {
  const Real mu1 = 0.4, mu2 = -0.7;
  // Covariance Sigma = L L^T with L lower-triangular:
  //   L = [[s1, 0],[rho*s2, sqrt(1-rho^2)*s2]]  => Sigma11=s1^2, Sigma22=s2^2,
  //   Sigma12 = rho*s1*s2.
  const Real s1 = 0.6, s2 = 0.9, rho = 0.5;
  const Real L11 = s1, L21 = rho * s2, L22 = std::sqrt(1.0 - rho * rho) * s2;
  const Real S11 = s1 * s1;
  const Real S12 = rho * s1 * s2;
  const Real S22 = s2 * s2;

  const int N = 6000000;
  std::mt19937_64 rng(0xABCDEF12u);
  std::normal_distribution<Real> g(0.0, 1.0);

  std::vector<std::vector<Real>> E(2, std::vector<Real>(N));
  for (int t = 0; t < N; ++t) {
    const Real z1 = g(rng), z2 = g(rng);
    E[0][t] = mu1 + L11 * z1;
    E[1][t] = mu2 + L21 * z1 + L22 * z2;
  }

  const std::vector<Real> lambda0{0.0, 0.0};
  const std::array<std::array<Real,2>,3> dls{{ {{0.3,-0.2}}, {{-0.4,0.25}}, {{0.15,0.35}} }};
  for (const auto& dl : dls) {
    const std::vector<Real> lambda{dl[0], dl[1]};
    const Real e1 = single_histogram(E, E[0], lambda0, lambda);
    const Real e2 = single_histogram(E, E[1], lambda0, lambda);

    // <E>_tilt = mu - Sigma . dlambda.
    const Real ex1 = mu1 - (S11 * dl[0] + S12 * dl[1]);
    const Real ex2 = mu2 - (S12 * dl[0] + S22 * dl[1]);

    // Effective-sample loss factor for a 2D tilt: exp( dl^T Sigma dl / 2 ).
    const Real quad = dl[0]*dl[0]*S11 + 2*dl[0]*dl[1]*S12 + dl[1]*dl[1]*S22;
    const Real ess = std::exp(0.5 * quad);
    const Real tol1 = 14.0 * s1 * ess / std::sqrt((Real)N);
    const Real tol2 = 14.0 * s2 * ess / std::sqrt((Real)N);

    char m[110];
    std::snprintf(m, sizeof m, "2-coupling <E1>_tilt(%.2f,%.2f) == mu1-(Sig.dl)1", dl[0], dl[1]);
    CHECK_CLOSE(e1, ex1, tol1, m);
    std::snprintf(m, sizeof m, "2-coupling <E2>_tilt(%.2f,%.2f) == mu2-(Sig.dl)2", dl[0], dl[1]);
    CHECK_CLOSE(e2, ex2, tol2, m);
  }
}

// ---------------------------------------------------------------------------- 3
// Multi-histogram (WHAM/MBAR) consistency. Three reference ensembles for a SINGLE
// coupling Gaussian model. Crucially the synthetic data must be self-consistent
// with the linear-action model: an ensemble sampled at coupling lambda0_k should
// have its energy distributed as N(mu - lambda0_k*sigma^2, sigma^2) -- exactly the
// Gaussian tilt of the lambda=0 law by lambda0_k. We build each ensemble that way,
// so that ALL three are reweightings of the same underlying exp(-(E-mu)^2/2sigma^2)
// reference; WHAM then has the correct generative model.
//
// We reweight to a common target lambda_T and assert:
//   (a) the WHAM combined estimate matches the exact tilted mean,
//   (b) each single-histogram estimate (from each ensemble alone) matches too,
//       hence agrees with WHAM, all within MC errors,
//   (c) the f_k iteration converged (residual < tol, iters < cap) and the recovered
//       free-energy DIFFERENCES match the analytic Gaussian free energy
//         f(lambda) = lambda*mu - 0.5*lambda^2*sigma^2  (up to a common constant),
//       since Z(lambda)/Z(0) = exp(0.5 lambda^2 sigma^2 - lambda mu) ... see below.
static void test_multi_histogram() {
  const Real mu = 0.5, sigma = 1.1;
  const Real s2 = sigma * sigma;
  const std::vector<Real> lam0s{-0.5, 0.0, 0.6};  // three reference couplings
  const int Nk = 800000;

  std::mt19937_64 rng(0x5EED1234u);

  std::vector<Ensemble> ens(lam0s.size());
  for (std::size_t k = 0; k < lam0s.size(); ++k) {
    const Real lk = lam0s[k];
    // Ensemble k is N(mu - lk*sigma^2, sigma^2): the exact tilt of the base law by lk.
    std::normal_distribution<Real> gk(mu - lk * s2, sigma);
    ens[k].lambda0 = {lk};
    ens[k].E.assign(1, std::vector<Real>(Nk));
    ens[k].O.assign(Nk, 0.0);
    for (int t = 0; t < Nk; ++t) {
      const Real e = gk(rng);
      ens[k].E[0][t] = e;
      ens[k].O[t]    = e;   // observable O = E
    }
  }

  const Real lamT = 0.25;                 // common target coupling
  const std::vector<Real> lambda{lamT};
  const Real exact = mu - lamT * s2;       // tilted mean at the target

  WhamSolution sol;
  const Real wham = multi_histogram(ens, lambda, 1e-12, 100000, &sol);

  // (c) convergence diagnostics.
  CHECK(sol.converged, "WHAM f_k iteration converged");
  CHECK(sol.residual < 1e-12, "WHAM residual < tol");
  CHECK(sol.iters > 0 && sol.iters < 100000, "WHAM iters within (0, cap)");
  CHECK(std::isfinite(wham), "WHAM target estimate finite");

  // Combined WHAM tolerance: pooled N, with reweighting inflation from the target.
  const Real Ntot = (Real)(lam0s.size() * (std::size_t)Nk);
  const Real tolW = 14.0 * sigma * std::exp((lamT * sigma) * (lamT * sigma))
                      / std::sqrt(Ntot);
  CHECK_CLOSE(wham, exact, tolW, "WHAM combined <E> == exact tilted mean");

  // (b) each single-histogram estimate agrees with the exact tilt and so with WHAM.
  for (std::size_t k = 0; k < lam0s.size(); ++k) {
    const Real sh = single_histogram(ens[k], lambda);
    const Real dl = lamT - lam0s[k];
    const Real tolk = 14.0 * sigma * std::exp((dl * sigma) * (dl * sigma))
                        / std::sqrt((Real)Nk);
    char m[110];
    std::snprintf(m, sizeof m, "single-hist ens[%zu] <E> == exact tilt", k);
    CHECK_CLOSE(sh, exact, tolk, m);
    // WHAM and the single-histogram estimate agree within the looser of the two errors.
    std::snprintf(m, sizeof m, "WHAM vs single-hist ens[%zu] agree", k);
    CHECK_CLOSE(wham, sh, tolk + tolW, m);
  }

  // (c continued) recovered free-energy differences vs analytic Gaussian.
  // With weight exp(-lambda E) and base law N(mu,sigma^2), Z(lambda)/Z(0) =
  //   E_0[ exp(-lambda E) ] = exp(-lambda mu + 0.5 lambda^2 sigma^2).
  // MBAR's f_k are defined so that exp(-f_k) ∝ Z(lambda0_k); i.e. f_k = -log Z(lambda0_k)
  // up to a common additive constant => f_k - f_0 = -[logZ(lam_k) - logZ(lam_0)]
  //   = -( (-lam_k mu + 0.5 lam_k^2 s2) - (-lam_0 mu + 0.5 lam_0^2 s2) ).
  // Our solver anchors f_0 := 0 (the FIRST ensemble, lam0s[0]).
  auto logZrel = [&](Real lam, Real lamRef) {
    const Real a = -lam    * mu + 0.5 * lam    * lam    * s2;
    const Real b = -lamRef * mu + 0.5 * lamRef * lamRef * s2;
    return a - b;  // logZ(lam) - logZ(lamRef)
  };
  CHECK_CLOSE(sol.f[0], 0.0, 1e-12, "WHAM f_0 anchored to 0");
  for (std::size_t k = 1; k < lam0s.size(); ++k) {
    const Real fk_exact = -logZrel(lam0s[k], lam0s[0]);   // f_k - f_0
    // Free-energy statistical error ~ a few / sqrt(Nk); allow a generous band.
    char m[96];
    std::snprintf(m, sizeof m, "WHAM f_%zu matches analytic Gaussian free energy", k);
    CHECK_CLOSE(sol.f[k], fk_exact, 0.02, m);
  }
}

// ---------------------------------------------------------------------------- 4
// Log-sum-exp numerical stability. Construct a case where dlambda*E reaches ~700+
// in magnitude (exp(700) ~ 1e304, near double overflow at ~1.8e308; exp(710)
// overflows). A naive estimator forming exp(+700) directly would overflow to inf;
// the log-sum-exp path must return the exact two-sample analytic answer.
//
// Reference: two samples E = {E_a, E_b}, observable O = {O_a, O_b}, tilt dl.
// Weights w_a = exp(-dl E_a), w_b = exp(-dl E_b). The ratio
//   <O> = (O_a w_a + O_b w_b)/(w_a + w_b)
// is invariant under multiplying both weights by a constant, so it equals
//   = (O_a + O_b r)/(1 + r),  r = exp(-dl (E_b - E_a)),
// which we compute from the SMALL difference dl*(E_b-E_a) exactly even though the
// individual exponents are huge. That closed form is the analytic target.
static void test_logsumexp_stability() {
  // Make individual exponents enormous: dl*E ~ 700.
  const Real dl = 7.0;
  const Real Ea = 100.0;   // -dl*Ea = -700
  const Real Eb = 99.0;    // -dl*Eb = -693  -> diff small & controlled
  const Real Oa = 2.0, Ob = 5.0;

  std::vector<std::vector<Real>> E(1, std::vector<Real>{Ea, Eb});
  std::vector<Real> O{Oa, Ob};
  const std::vector<Real> lambda0{0.0}, lambda{dl};

  const Real got = single_histogram(E, O, lambda0, lambda);

  // Exact closed form via the small difference (no large exponent ever formed).
  const Real r = std::exp(-dl * (Eb - Ea));   // exp(-7*(-1)) = exp(7)
  const Real exact = (Oa + Ob * r) / (1.0 + r);

  CHECK(std::isfinite(got), "stability: result finite (no inf/nan)");
  CHECK_CLOSE(got, exact, 1e-12, "stability: huge-exponent reweight == closed form");

  // Also exercise the opposite-sign extreme (large POSITIVE exponents would be the
  // overflow danger for a naive code): tilt the other way.
  const std::vector<Real> lambdaNeg{-dl};
  const Real got2 = single_histogram(E, O, lambda0, lambdaNeg);
  const Real r2 = std::exp(dl * (Eb - Ea));
  const Real exact2 = (Oa + Ob * r2) / (1.0 + r2);
  CHECK(std::isfinite(got2), "stability: opposite tilt finite");
  CHECK_CLOSE(got2, exact2, 1e-12, "stability: opposite-sign huge exponent == closed form");

  // Direct log_sum_exp identities.
  //  log_sum_exp({a}) == a; shift-invariance log_sum_exp({a+c,b+c}) = c + lse({a,b}).
  CHECK_CLOSE(log_sum_exp(std::vector<Real>{700.0}), 700.0, 1e-9,
              "log_sum_exp single element == that element (no overflow)");
  const Real base = log_sum_exp(std::vector<Real>{0.3, -1.1, 2.0});
  const Real shifted = log_sum_exp(std::vector<Real>{0.3 + 500.0, -1.1 + 500.0, 2.0 + 500.0});
  CHECK_CLOSE(shifted, base + 500.0, 1e-9, "log_sum_exp shift-invariance at +500");
  CHECK(std::isfinite(shifted), "log_sum_exp(+500 shift) finite");
  // Closed form for two equal exponents: log_sum_exp({x,x}) = x + log 2.
  CHECK_CLOSE(log_sum_exp(std::vector<Real>{350.0, 350.0}), 350.0 + std::log(2.0), 1e-9,
              "log_sum_exp({x,x}) == x + log2");
}

// ---------------------------------------------------------------------------- 5
// Jackknife helper sanity: mean matches single_histogram; identical-sample data
// gives zero error; the reported error is positive on scattered data and matches
// an independent delete-1 reference re-implemented here.
static void test_jackknife() {
  const Real mu = 0.2, sigma = 0.5;
  const int N = 2000;
  std::mt19937_64 rng(0x1357u);
  std::normal_distribution<Real> gN(mu, sigma);
  std::vector<std::vector<Real>> E(1, std::vector<Real>(N));
  std::vector<Real> O(N);
  for (int t = 0; t < N; ++t) { Real e = gN(rng); E[0][t] = e; O[t] = e; }

  const std::vector<Real> lambda0{0.0}, lambda{0.3};
  const ReweightEstimate est = reweighted_mean_jackknife(E, O, lambda0, lambda);
  const Real direct = single_histogram(E, O, lambda0, lambda);
  CHECK_CLOSE(est.mean, direct, 1e-12, "jackknife mean == single_histogram mean");
  CHECK(est.err > 0.0 && std::isfinite(est.err), "jackknife err positive & finite");

  // Independent brute-force delete-1 reference.
  std::vector<Real> a(N);
  for (int t = 0; t < N; ++t) a[t] = -(lambda[0] - lambda0[0]) * E[0][t];
  auto wmean = [&](int drop) {
    Real mx = -1e300; for (int t = 0; t < N; ++t) if (t != drop) mx = std::max(mx, a[t]);
    Real num = 0, den = 0;
    for (int t = 0; t < N; ++t) if (t != drop) { Real w = std::exp(a[t]-mx); num += O[t]*w; den += w; }
    return num/den;
  };
  std::vector<Real> th(N); for (int k = 0; k < N; ++k) th[k] = wmean(k);
  Real tbar = 0; for (Real v : th) tbar += v; tbar /= N;
  Real sw = 0; for (Real v : th) sw += (v-tbar)*(v-tbar);
  const Real refErr = std::sqrt((Real)(N-1)/N * sw);
  CHECK_CLOSE(est.err, refErr, 1e-10, "jackknife err == independent delete-1 reference");

  // n<2 -> zero error.
  std::vector<std::vector<Real>> E1(1, std::vector<Real>{0.7});
  std::vector<Real> O1{0.7};
  const ReweightEstimate e1 = reweighted_mean_jackknife(E1, O1, lambda0, lambda);
  CHECK_CLOSE(e1.err, 0.0, 1e-14, "jackknife n=1 -> err 0");
  CHECK_CLOSE(e1.mean, 0.7, 1e-14, "jackknife n=1 -> mean is the lone sample");
}

int main() {
  std::printf("-- single-histogram 1 coupling (Gaussian tilt) --\n"); test_single_1coupling();
  std::printf("-- single-histogram 2 couplings --\n");                 test_two_coupling();
  std::printf("-- multi-histogram WHAM/MBAR --\n");                     test_multi_histogram();
  std::printf("-- log-sum-exp stability --\n");                        test_logsumexp_stability();
  std::printf("-- jackknife helper --\n");                             test_jackknife();
  return report("test_reweight");
}
