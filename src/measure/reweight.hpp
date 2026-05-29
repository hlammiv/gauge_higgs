#pragma once
// Multi-coupling histogram reweighting (Ferrenberg-Swendsen / WHAM / MBAR).
//
// PURE post-processing of Monte Carlo time series -- no lattice, no group, no
// dimension. Given measurements made at one or several reference couplings, it
// estimates expectation values <O> at arbitrary TARGET couplings by reweighting,
// so a small set of simulated (beta, kappa, ...) points can trace a transition
// line / ridge and localize a triple point by interpolation.
//
// MODEL / CONVENTION (fixed; the tests pin the exact analytic identities)
//   The Monte Carlo weight of a configuration is exp(-S) with the action LINEAR
//   in the couplings,
//       S(config; lambda) = sum_i lambda_i E_i(config),
//   where lambda is the vector of couplings and E_i are the conjugate "energy"
//   observables. The CALLER fixes the signs so that this holds. For example for
//   the U(1) gauge+Higgs model used elsewhere in this code one may take
//       lambda = (beta, kappa, ...),
//       E_1 = sum_plaq (1 - cos theta_plaq),          (conjugate to beta)
//       E_2 = - sum_x 2 Re[ conj(phi) e^{i q theta} phi ],   (conjugate to kappa)
//   so that S = beta E_1 + kappa E_2 + ... . The module itself is action-agnostic:
//   it only ever sees numerical time series of the E_i and of the observable O.
//
//   Reweighting from a reference coupling lambda0 to a target lambda multiplies
//   each sample's relative weight by
//       exp( -(S(config;lambda) - S(config;lambda0)) )
//         = exp( -sum_i (lambda_i - lambda0_i) E_i(config) ).
//   Note the SIGN: increasing a coupling lambda_i suppresses configurations with
//   large conjugate energy E_i.
//
// NUMERICAL STABILITY
//   The reweighting exponents -sum_i (lambda_i - lambda0_i) E_i can be large and
//   positive, so exp() of them overflows. Every sum-of-exponentials below is done
//   with the LOG-SUM-EXP trick: shift by the maximum exponent, exponentiate the
//   (now <= 0) shifted exponents, sum, and add back the max in log space. We never
//   form exp() of a large positive number directly.
//
// CONTENTS
//   struct Ensemble                 -- one reference ensemble (lambda0, E matrix, O)
//   single_histogram(...)           -- one-ensemble Ferrenberg-Swendsen estimator
//   multi_histogram(...)            -- K-ensemble WHAM/MBAR estimator (self-consistent f_k)
//   reweighted_mean_jackknife(...)  -- single-histogram mean + delete-1 jackknife error
//   log_sum_exp(...)                -- the shared stable accumulator (exposed for tests)
//
// All arrays are plain std::vector; E is stored as E[i][t] (coupling-major,
// sample-minor) to match the "time series of the E_i" description.
#include "core/config.hpp"
#include <vector>
#include <cmath>
#include <cstddef>
#include <limits>
#include <algorithm>

namespace gh {

// ---------------------------------------------------------------------------
//  Numerically stable log-sum-exp:  log( sum_t exp(a_t) ).
//  Returns -inf for an empty input (the additive identity of a log-sum). The
//  shift by max(a) keeps every exponentiated term in (0,1], so no overflow.
// ---------------------------------------------------------------------------
inline Real log_sum_exp(const std::vector<Real>& a) {
  if (a.empty()) return -std::numeric_limits<Real>::infinity();
  Real m = -std::numeric_limits<Real>::infinity();
  for (Real v : a) m = std::max(m, v);
  if (!std::isfinite(m)) return m;  // all -inf (or a +inf) -> pass it through
  Real s = 0.0;
  for (Real v : a) s += std::exp(v - m);
  return m + std::log(s);
}

namespace reweight_detail {

constexpr Real kNaN = std::numeric_limits<Real>::quiet_NaN();

// Weighted-sum-exp ratio  ( sum_t o_t exp(a_t) ) / ( sum_t exp(a_t) )  done
// stably: shift by max(a) in BOTH sums so the common factor exp(max) cancels.
// o_t may be any sign. Returns NaN on empty input.
inline Real weighted_mean_exp(const std::vector<Real>& o, const std::vector<Real>& a) {
  const std::size_t n = a.size();
  if (n == 0) return kNaN;
  Real m = -std::numeric_limits<Real>::infinity();
  for (Real v : a) m = std::max(m, v);
  Real num = 0.0, den = 0.0;
  for (std::size_t t = 0; t < n; ++t) {
    const Real w = std::exp(a[t] - m);  // in (0,1]
    num += o[t] * w;
    den += w;
  }
  return num / den;
}

}  // namespace reweight_detail

// ---------------------------------------------------------------------------
//  One reference ensemble.
//    lambda0          : the couplings at which it was sampled (length n_couplings)
//    E[i][t]          : time series of conjugate energy E_i, i in [0,n_couplings),
//                       t in [0,n_samples). All rows MUST have the same length.
//    O[t]             : optional observable time series (same length as a row of E).
//                       If empty, single_histogram below can still reweight any of
//                       the E_i themselves (pass them in as O explicitly).
// ---------------------------------------------------------------------------
struct Ensemble {
  std::vector<Real> lambda0;
  std::vector<std::vector<Real>> E;
  std::vector<Real> O;

  std::size_t n_couplings() const { return lambda0.size(); }
  std::size_t n_samples() const { return E.empty() ? 0u : E[0].size(); }
};

// ---------------------------------------------------------------------------
//  SINGLE-HISTOGRAM (Ferrenberg-Swendsen) reweighting.
//
//  Estimate <O> at target coupling `lambda` from a single reference ensemble
//  sampled at `lambda0`:
//      a_t      = -sum_i (lambda_i - lambda0_i) E_i(t)      (log relative weight)
//      <O>_lam  = ( sum_t O_t e^{a_t} ) / ( sum_t e^{a_t} ).
//  Implemented with log-sum-exp (shift by max_t a_t).
//
//  `E`      : E[i][t], the reference energy time series.
//  `O`      : observable time series (same length as a row of E).
//  `lambda0`: reference couplings (length = E.size()).
//  `lambda` : target couplings (length = E.size()).
//  Returns NaN on shape mismatch or empty input.
// ---------------------------------------------------------------------------
inline Real single_histogram(const std::vector<std::vector<Real>>& E,
                             const std::vector<Real>& O,
                             const std::vector<Real>& lambda0,
                             const std::vector<Real>& lambda) {
  const std::size_t nc = E.size();
  if (nc == 0 || lambda0.size() != nc || lambda.size() != nc) return reweight_detail::kNaN;
  const std::size_t ns = E[0].size();
  if (ns == 0 || O.size() != ns) return reweight_detail::kNaN;
  for (std::size_t i = 0; i < nc; ++i)
    if (E[i].size() != ns) return reweight_detail::kNaN;

  // Per-sample log weight a_t = -sum_i dlambda_i E_i(t), dlambda = lambda - lambda0.
  std::vector<Real> a(ns, 0.0);
  for (std::size_t i = 0; i < nc; ++i) {
    const Real dl = lambda[i] - lambda0[i];
    if (dl == 0.0) continue;
    const std::vector<Real>& Ei = E[i];
    for (std::size_t t = 0; t < ns; ++t) a[t] -= dl * Ei[t];
  }
  return reweight_detail::weighted_mean_exp(O, a);
}

// Convenience overload taking an Ensemble (uses ensemble.O).
inline Real single_histogram(const Ensemble& ens, const std::vector<Real>& lambda) {
  return single_histogram(ens.E, ens.O, ens.lambda0, lambda);
}

// ---------------------------------------------------------------------------
//  SINGLE-HISTOGRAM mean + delete-1 jackknife error over samples.
//
//  Returns {mean, err}. The jackknife recomputes the reweighted mean on each
//  leave-one-out subset and forms the standard delete-1 jackknife error
//      err^2 = (n-1)/n * sum_k ( theta_k - theta_bar )^2 .
//  Each leave-one-out estimate is itself a stable weighted_mean_exp. With n<2
//  the error is 0. This is O(n^2) by design (transparent reference quality);
//  it is meant for modest n typical of reweighting studies.
// ---------------------------------------------------------------------------
struct ReweightEstimate {
  Real mean = reweight_detail::kNaN;
  Real err  = 0.0;
};

inline ReweightEstimate reweighted_mean_jackknife(const std::vector<std::vector<Real>>& E,
                                                  const std::vector<Real>& O,
                                                  const std::vector<Real>& lambda0,
                                                  const std::vector<Real>& lambda) {
  const std::size_t nc = E.size();
  if (nc == 0 || lambda0.size() != nc || lambda.size() != nc) return {};
  const std::size_t ns = E[0].size();
  if (ns == 0 || O.size() != ns) return {};
  for (std::size_t i = 0; i < nc; ++i)
    if (E[i].size() != ns) return {};

  // Full-sample log weights.
  std::vector<Real> a(ns, 0.0);
  for (std::size_t i = 0; i < nc; ++i) {
    const Real dl = lambda[i] - lambda0[i];
    if (dl == 0.0) continue;
    for (std::size_t t = 0; t < ns; ++t) a[t] -= dl * E[i][t];
  }

  ReweightEstimate out;
  out.mean = reweight_detail::weighted_mean_exp(O, a);
  if (ns < 2) { out.err = 0.0; return out; }

  // Delete-1 jackknife replicas.
  std::vector<Real> th(ns);
  std::vector<Real> aj(ns - 1), oj(ns - 1);
  for (std::size_t k = 0; k < ns; ++k) {
    std::size_t j = 0;
    for (std::size_t t = 0; t < ns; ++t) {
      if (t == k) continue;
      aj[j] = a[t];
      oj[j] = O[t];
      ++j;
    }
    th[k] = reweight_detail::weighted_mean_exp(oj, aj);
  }
  Real tbar = 0.0; for (Real v : th) tbar += v; tbar /= static_cast<Real>(ns);
  Real sw = 0.0;   for (Real v : th) sw += (v - tbar) * (v - tbar);
  out.err = std::sqrt(static_cast<Real>(ns - 1) / static_cast<Real>(ns) * sw);
  return out;
}

inline ReweightEstimate reweighted_mean_jackknife(const Ensemble& ens,
                                                  const std::vector<Real>& lambda) {
  return reweighted_mean_jackknife(ens.E, ens.O, ens.lambda0, lambda);
}

// ---------------------------------------------------------------------------
//  MULTI-HISTOGRAM (WHAM / MBAR).
//
//  Combine K reference ensembles sampled at couplings {lambda0_k}, with N_k
//  samples each, into a single optimal estimator. Index every sample by its
//  global position n (n = 1..M_total) regardless of which ensemble produced it;
//  let u_k(n) = sum_i lambda0_k[i] E_i(n) be the (negative log) reduced action
//  of sample n evaluated AT ensemble k's couplings, i.e. S = u_k for that k.
//
//  MBAR / WHAM self-consistency. The dimensionless free energies f_k (defined up
//  to an additive constant) satisfy the fixed-point equations
//      exp(-f_k) = sum_n  exp(-u_k(n)) / ( sum_l N_l exp(f_l - u_l(n)) ).      (1)
//  Taking logs and using log-sum-exp, with
//      g(n) = log_sum_exp_l [ log N_l + f_l - u_l(n) ]    (the per-sample
//             normalization / "mixture" denominator, common to all targets),
//  equation (1) becomes
//      f_k  <-  -log_sum_exp_n [ -u_k(n) - g(n) ].                            (1')
//  We iterate (1') to convergence, re-anchoring f_0 := 0 each sweep to fix the
//  gauge. CONVERGENCE CRITERION: stop when max_k |f_k^{new} - f_k^{old}| < tol
//  (after re-anchoring), or when `max_iter` sweeps are reached. The returned
//  WhamSolution reports the achieved `residual` and `iters` so callers can assert.
//
//  Reweighted expectation at an ARBITRARY target coupling lambda (reduced action
//  u_t(n) = sum_i lambda[i] E_i(n)) is the MBAR average
//      <O>_lambda = sum_n O(n) W(n) / sum_n W(n),
//      log W(n)   = -u_t(n) - g(n)                                            (2)
//  i.e. each sample is weighted by exp(-u_t(n)) divided by the same mixture
//  denominator exp(g(n)). Done with log-sum-exp via weighted_mean_exp on
//  a(n) = -u_t(n) - g(n).
//
//  Because g(n) is independent of the target, the f_k solve is done ONCE and any
//  number of targets are then cheap. multi_histogram() exposes both: solve_wham()
//  returns the solution (f_k, g(n), residual, iters); multi_histogram_eval()
//  reweights a target using a precomputed solution; multi_histogram() is the
//  one-shot convenience that solves then evaluates a single target.
// ---------------------------------------------------------------------------
struct WhamSolution {
  std::vector<Real> f;          // dimensionless free energies f_k, f_0 anchored to 0
  std::vector<Real> g;          // per-sample log mixture denominator g(n)
  std::vector<Real> O;          // flattened observable per global sample n
  std::vector<std::vector<Real>> E;  // flattened energies E[i][n] (coupling-major)
  Real residual = reweight_detail::kNaN;  // achieved max_k |df_k| at the last sweep
  int  iters    = 0;            // sweeps performed
  bool converged = false;       // residual < tol within max_iter

  std::size_t n_couplings() const { return E.size(); }
  std::size_t n_total() const { return O.size(); }
};

// Solve the MBAR/WHAM self-consistency for the free energies f_k.
//   `ensembles` : K reference ensembles. Every ensemble must share the same
//                 n_couplings; each must carry its O time series (used by eval).
//   `tol`       : convergence threshold on max_k |df_k|.
//   `max_iter`  : iteration cap.
inline WhamSolution solve_wham(const std::vector<Ensemble>& ensembles,
                               Real tol = 1e-10, int max_iter = 100000) {
  WhamSolution sol;
  const std::size_t K = ensembles.size();
  if (K == 0) return sol;
  const std::size_t nc = ensembles[0].n_couplings();
  if (nc == 0) return sol;

  // Flatten: global samples n = 0..M-1; remember each ensemble's N_k and the
  // per-ensemble reference couplings.
  std::vector<Real> Nk(K, 0.0);
  std::vector<std::vector<Real>> lam0(K);
  std::size_t M = 0;
  for (std::size_t k = 0; k < K; ++k) {
    if (ensembles[k].n_couplings() != nc) return sol;  // shape mismatch
    const std::size_t Nk_ = ensembles[k].n_samples();
    Nk[k] = static_cast<Real>(Nk_);
    lam0[k] = ensembles[k].lambda0;
    M += Nk_;
  }
  if (M == 0) return sol;

  sol.E.assign(nc, std::vector<Real>(M, 0.0));
  sol.O.assign(M, 0.0);
  {
    std::size_t off = 0;
    for (std::size_t k = 0; k < K; ++k) {
      const Ensemble& e = ensembles[k];
      const std::size_t Nk_ = e.n_samples();
      for (std::size_t i = 0; i < nc; ++i) {
        if (e.E[i].size() != Nk_) { sol = WhamSolution{}; return sol; }
        for (std::size_t t = 0; t < Nk_; ++t) sol.E[i][off + t] = e.E[i][t];
      }
      const bool haveO = e.O.size() == Nk_;
      for (std::size_t t = 0; t < Nk_; ++t) sol.O[off + t] = haveO ? e.O[t] : 0.0;
      off += Nk_;
    }
  }

  // Precompute u_k(n) = sum_i lambda0_k[i] E_i(n) for every ensemble k and sample n.
  // Stored as U[k][n].
  std::vector<std::vector<Real>> U(K, std::vector<Real>(M, 0.0));
  for (std::size_t k = 0; k < K; ++k)
    for (std::size_t i = 0; i < nc; ++i) {
      const Real l = lam0[k][i];
      if (l == 0.0) continue;
      for (std::size_t n = 0; n < M; ++n) U[k][n] += l * sol.E[i][n];
    }

  std::vector<Real> logNk(K);
  for (std::size_t k = 0; k < K; ++k)
    logNk[k] = (Nk[k] > 0.0) ? std::log(Nk[k]) : -std::numeric_limits<Real>::infinity();

  std::vector<Real> f(K, 0.0), fnew(K, 0.0);
  std::vector<Real> g(M, 0.0);
  std::vector<Real> tmpK(K), tmpM(M);

  Real residual = std::numeric_limits<Real>::infinity();
  int it = 0;
  for (; it < max_iter; ++it) {
    // g(n) = log_sum_exp_l [ logN_l + f_l - u_l(n) ].
    for (std::size_t n = 0; n < M; ++n) {
      for (std::size_t k = 0; k < K; ++k) tmpK[k] = logNk[k] + f[k] - U[k][n];
      g[n] = log_sum_exp(tmpK);
    }
    // f_k(new) = -log_sum_exp_n [ -u_k(n) - g(n) ].
    for (std::size_t k = 0; k < K; ++k) {
      for (std::size_t n = 0; n < M; ++n) tmpM[n] = -U[k][n] - g[n];
      fnew[k] = -log_sum_exp(tmpM);
    }
    // Re-anchor to fix the additive gauge: f_0 := 0.
    const Real shift = fnew[0];
    for (std::size_t k = 0; k < K; ++k) fnew[k] -= shift;
    // Residual after re-anchoring.
    residual = 0.0;
    for (std::size_t k = 0; k < K; ++k)
      residual = std::max(residual, std::fabs(fnew[k] - f[k]));
    f.swap(fnew);
    if (residual < tol) { ++it; break; }
  }

  // Final g(n) consistent with the converged f.
  for (std::size_t n = 0; n < M; ++n) {
    for (std::size_t k = 0; k < K; ++k) tmpK[k] = logNk[k] + f[k] - U[k][n];
    g[n] = log_sum_exp(tmpK);
  }

  sol.f = std::move(f);
  sol.g = std::move(g);
  sol.residual = residual;
  sol.iters = it;
  sol.converged = residual < tol;
  return sol;
}

// Reweight a target coupling using a precomputed WHAM solution.
//   <O>_lambda = ( sum_n O(n) exp(a_n) ) / ( sum_n exp(a_n) ),
//   a_n = -u_t(n) - g(n),  u_t(n) = sum_i lambda[i] E_i(n).
// `obs` overrides the flattened observable if non-empty (must be length n_total);
// otherwise the solution's stored O is used. Returns NaN on shape mismatch.
inline Real multi_histogram_eval(const WhamSolution& sol, const std::vector<Real>& lambda,
                                 const std::vector<Real>& obs = {}) {
  const std::size_t nc = sol.n_couplings();
  const std::size_t M = sol.n_total();
  if (nc == 0 || M == 0 || lambda.size() != nc) return reweight_detail::kNaN;
  const std::vector<Real>& O = obs.empty() ? sol.O : obs;
  if (O.size() != M) return reweight_detail::kNaN;

  std::vector<Real> a(M);
  for (std::size_t n = 0; n < M; ++n) a[n] = -sol.g[n];
  for (std::size_t i = 0; i < nc; ++i) {
    const Real l = lambda[i];
    if (l == 0.0) continue;
    for (std::size_t n = 0; n < M; ++n) a[n] -= l * sol.E[i][n];
  }
  return reweight_detail::weighted_mean_exp(O, a);
}

// One-shot WHAM/MBAR: solve f_k from the K ensembles, then reweight a single
// target coupling. `out_sol`, if non-null, receives the solution (f_k, residual,
// iters, converged) so the caller can assert convergence and reuse it.
inline Real multi_histogram(const std::vector<Ensemble>& ensembles,
                            const std::vector<Real>& lambda,
                            Real tol = 1e-10, int max_iter = 100000,
                            WhamSolution* out_sol = nullptr) {
  WhamSolution sol = solve_wham(ensembles, tol, max_iter);
  if (out_sol) *out_sol = sol;
  return multi_histogram_eval(sol, lambda);
}

}  // namespace gh
