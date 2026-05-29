// Bounded-from-below (BFB) conditions for the SU(3) quartic Higgs potential.
//
//   V_4(Phi) = sum_{rho'} f_{rho'} V_{rho'}(Phi),   V_{rho'}(Phi) = ||P_{rho'} M||^2 >= 0,
//
// with M = Phi Phi^dag (the d x d bilinear) and P_{rho'} the projector onto the
// rho'-th adjoint-Casimir channel of R (x) Rbar (our CasimirChannels<3>; the draft's
// "rho'" channels). Each V_{rho'} is a sum of squared magnitudes |Phi^dag Q_{rho'n} Phi|^2,
// hence individually >= 0; the channel value here is reTrDagProd(M, P_c M) = ||P_c M||^2.
//
// KEY POINT (this is the whole content of the SU(3) BFB question):
//   * For SU(2) integer-spin (real) reps the independent quartic invariants reduce to
//     QUADRATICS in the field bilinears, so BFB is trivial (a quadratic-form positivity).
//   * For SU(3) V_4 is genuinely QUARTIC in Phi: V_4 = sum_c f_c ||P_c M||^2 is a quadratic
//     form in the components of M, but M is constrained to the rank-1 cone M = Phi Phi^dag,
//     so BFB is a COPOSITIVITY condition on {f_c} over that cone, NOT plain positivity.
//   * all-f_c >= 0  =>  BFB  (each term >= 0): SUFFICIENT but NOT necessary. Physical minima
//     generally require some f_c < 0; how negative is bounded by the cross-channel structure
//     (sum_c P_c = 1, sum_c lambda_c P_c = Chat, etc.), which is exactly what the numeric
//     copositivity checker below probes.
//
// Because V_4 is homogeneous of degree 4 in Phi, V_4(t Phi) = t^4 V_4(Phi); so V_4 is
// bounded below (over all C^d) IFF min over the unit sphere |Phi| = 1 is >= 0. (A single
// unit vector with V_4 < 0 gives the unbounded ray t Phi, t -> infinity.) The checker
// minimizes V_4 on the sphere by many random restarts + projected-gradient descent and
// declares BFB iff that minimum is >= -eps.
//
// Build / run:
//   g++ -std=c++20 -O3 -march=native -fopenmp -Isrc -o build/bfb_su3 src/bfb_su3.cpp && ./build/bfb_su3
//
// See docs/draft_validation_bfb.md, memory: research-program / discrete-observables-program.
#include "action/scalar_invariants.hpp"
#include "action/vacuum_alignment.hpp"
#include "rep/rep_general.hpp"
#include "rep/finite_subgroups.hpp"
#include <random>
#include <cstdio>
#include <algorithm>

namespace gh {

struct BfbResult {
  Real vmin;        // minimum of V_4 over the unit sphere found across restarts
  bool bfb;         // vmin >= -eps
  DVec argmin;      // unit vector attaining vmin (an unbounded direction if !bfb)
};

// The combined weighted superoperator A_f = sum_c f_c P_c, applied to M, is a single
// polynomial in the adjoint-Casimir superoperator Chat:  A_f = sum_c f_c L_c(Chat),
// where L_c is the Lagrange interpolant with L_c(lambda_c')=delta_{cc'}. Hence
// A_f = poly(Chat) with poly(lambda_c) = f_c. We precompute the Krylov powers
// {M, Chat M, Chat^2 M, ...} ONCE (nc-1 Casimir applications total) and combine with the
// polynomial coefficients of poly -- instead of nc separate Lagrange chains (nc*(nc-1)
// Casimir applications). Big speedup for the sphere scan.
struct WeightedSuperOp {
  const CasimirChannels<3>* ch;
  std::vector<Real> coef;          // poly(x) = sum_k coef[k] x^k, poly(lambda_c) = f_c
  int npow;                        // number of Krylov powers = coef.size()

  WeightedSuperOp(const CasimirChannels<3>& c, const std::vector<Real>& f) : ch(&c) {
    const int nc = ch->n_channels();
    // Solve Vandermonde V a = f for poly coefficients (nc nodes lambda_c, degree nc-1).
    std::vector<Real> Vm(static_cast<std::size_t>(nc) * nc), a(f);
    for (int i = 0; i < nc; ++i) { Real p = 1.0; for (int j = 0; j < nc; ++j) { Vm[i * nc + j] = p; p *= ch->lambda[i]; } }
    // Gaussian elimination with partial pivoting on [Vm | a].
    for (int col = 0; col < nc; ++col) {
      int piv = col; Real best = std::fabs(Vm[col * nc + col]);
      for (int r = col + 1; r < nc; ++r) { Real v = std::fabs(Vm[r * nc + col]); if (v > best) { best = v; piv = r; } }
      if (piv != col) { for (int j = 0; j < nc; ++j) std::swap(Vm[piv * nc + j], Vm[col * nc + j]); std::swap(a[piv], a[col]); }
      Real inv = 1.0 / Vm[col * nc + col];
      for (int r = col + 1; r < nc; ++r) {
        Real fct = Vm[r * nc + col] * inv;
        for (int j = col; j < nc; ++j) Vm[r * nc + j] -= fct * Vm[col * nc + j];
        a[r] -= fct * a[col];
      }
    }
    for (int i = nc - 1; i >= 0; --i) { Real s = a[i]; for (int j = i + 1; j < nc; ++j) s -= Vm[i * nc + j] * a[j]; a[i] = s / Vm[i * nc + i]; }
    coef = a; npow = nc;
  }

  // A_f M = sum_k coef[k] Chat^k M.
  DMat apply(const DMat& M) const {
    DMat out(M.rows, M.cols);
    DMat pw = M;                                   // Chat^0 M
    for (int k = 0; k < npow; ++k) {
      if (coef[k] != 0.0) out = out + pw * Complex(coef[k], 0);
      if (k + 1 < npow) pw = ch->apply_casimir(pw);
    }
    return out;
  }
};

// Fused value + gradient at one point via the combined superoperator A_f = sum_c f_c P_c:
//   V_4 = <M, A_f M>            (quartic part, mu2 = 0; A_f is Hermitian, V_4 real)
//   g   = dV_4/dphibar = 2 (A_f M) phi
// One A_f-apply (npow-1 Casimir applications) gives both -- much cheaper than per-channel
// Lagrange chains. (mu2 = 0: BFB only cares about the quartic.)
inline Real V4_and_grad(const WeightedSuperOp& A, const DVec& phi, DVec* grad) {
  const int d = A.ch->d;
  DMat M = CasimirChannels<3>::outer(phi);
  DMat AM = A.apply(M);
  Real v = reTrDagProd(M, AM);
  if (grad) {
    DVec AMphi = AM * phi;
    *grad = DVec(d);
    for (int e = 0; e < d; ++e) (*grad)(e) = Complex(2.0, 0) * AMphi(e);
  }
  return v;
}

// V_4 over all channels (quartic part only; mu2 = 0). Homogeneous degree 4 in Phi.
inline Real V4(const WeightedSuperOp& A, const DVec& phi) {
  return V4_and_grad(A, phi, nullptr);
}

// Free complex gradient dV_4/dphibar at phi (mu2 = 0). Re<phi, g> = 2 V_4 (Euler, deg 4).
inline DVec gradV4(const WeightedSuperOp& A, const DVec& phi) {
  DVec g(A.ch->d);
  V4_and_grad(A, phi, &g);
  return g;
}

// One projected-gradient descent from a given start (unit vector). Returns the minimized
// V_4 and overwrites p with the minimizer. Riemannian gradient = free gradient minus its
// radial part; step along -rg and renormalize, with a backtracking line search. Stops on a
// flat gradient or vanishing relative improvement.
inline Real descend(const WeightedSuperOp& A, DVec& p, int nstep) {
  const int d = A.ch->d;
  auto normalize = [&](DVec& q) {
    Real nrm = std::sqrt(q.norm2());
    for (int k = 0; k < d; ++k) q(k) *= Complex(1.0 / nrm, 0);
  };
  Real v = V4(A, p);
  Real step = 0.25;
  int stall = 0;
  for (int it = 0; it < nstep; ++it) {
    DVec g = gradV4(A, p);                              // dV/dphibar
    Complex rad = dot(p, g);                            // radial (p^dag g), p unit
    DVec rg(d); Real gn2 = 0;
    for (int k = 0; k < d; ++k) { rg(k) = g(k) - rad * p(k); gn2 += std::norm(rg(k)); }
    // Converged: tangential gradient small relative to the local scale (V_4 ~ O(1) on
    // the sphere). Stops flat minima from burning the full nstep -- the dominant cost.
    if (gn2 < 1e-10 * (1.0 + v * v)) break;
    // Backtracking line search along -rg. Cap step growth modestly (a successful step
    // grows step only a little) so the next iteration rarely backtracks more than once.
    bool moved = false;
    for (int bt = 0; bt < 24; ++bt) {
      DVec q(d);
      for (int k = 0; k < d; ++k) q(k) = p(k) - Complex(step, 0) * rg(k);
      normalize(q);
      Real vq = V4(A, q);
      if (vq < v - 1e-15) {
        Real rel = (v - vq) / (1e-14 + std::fabs(v));
        p = q; v = vq; step *= 1.1; moved = true;
        if (rel < 1e-7) { if (++stall > 3) break; } else stall = 0;
        break;
      }
      step *= 0.4;
    }
    if (!moved) { step *= 0.4; if (++stall > 6 || step < 1e-13) break; }
  }
  return v;
}

// Minimize V_4 on the unit sphere |phi| = 1 (complex d-vector) by random-restart
// projected-gradient descent. V_4 is homogeneous degree 4, so bounded below over C^d IFF
// this sphere minimum is >= 0; BFB iff min >= -eps. The independent restarts are run in
// parallel (each its own RNG stream), reducing best vmin under a critical section.
inline BfbResult check_bfb(const CasimirChannels<3>& ch, const std::vector<Real>& f,
                           int nrestart = 200, int nstep = 120, std::uint64_t seed = 7,
                           Real eps = 1e-7) {
  const int d = ch.d;
  WeightedSuperOp A(ch, f);                             // A_f = sum_c f_c P_c, built once
  BfbResult best; best.vmin = 1e300; best.bfb = true; best.argmin = DVec(d);

  #pragma omp parallel
  {
    Real loc_v = 1e300; DVec loc_p(d);
    #pragma omp for schedule(dynamic) nowait
    for (int r = 0; r < nrestart; ++r) {
      std::mt19937_64 rng(seed + 0x9E3779B97F4A7C15ULL * static_cast<std::uint64_t>(r));
      std::normal_distribution<Real> gd(0.0, 1.0);
      DVec p(d);
      for (int k = 0; k < d; ++k) p(k) = Complex(gd(rng), gd(rng));
      Real nrm = std::sqrt(p.norm2());
      for (int k = 0; k < d; ++k) p(k) *= Complex(1.0 / nrm, 0);
      Real v = descend(A, p, nstep);
      if (v < loc_v) { loc_v = v; loc_p = p; }
    }
    #pragma omp critical
    if (loc_v < best.vmin) { best.vmin = loc_v; best.argmin = loc_p; }
  }
  best.bfb = (best.vmin >= -eps);
  return best;
}

}  // namespace gh

#ifndef BFB_NO_MAIN
using namespace gh;

static void print_channels(const CasimirChannels<3>& ch, const char* label) {
  std::printf("%s: %d channels  C2 = ", label, ch.n_channels());
  for (Real c : ch.lambda) std::printf("%.4g ", c);
  std::printf("\n");
}

int main() {
  std::printf("==== SU(3) bounded-from-below (BFB) analysis ====\n\n");
  std::printf("V_4(Phi) = sum_c f_c ||P_c M||^2,  M = Phi Phi^dag,  homogeneous deg-4.\n");
  std::printf("BFB  <=>  min_{|Phi|=1} V_4 >= 0  (copositivity of {f_c} on the rank-1 cone).\n");
  std::printf("all-f_c>=0 is SUFFICIENT but not necessary.\n\n");

  // (2,2) = Young rows {4,2}, d = 27 -- the Sigma(108) Higgs irrep.
  GeneralRep<3> rep({4, 2});
  CasimirChannels<3> ch(rep);
  const int nc = ch.n_channels();
  print_channels(ch, "(2,2)={4,2} d=27 [Sigma(108)]");
  std::printf("\n");

  // ---- 1. all-positive f -> BFB ----
  {
    std::vector<Real> f(nc, 1.0);
    auto res = check_bfb(ch, f);
    std::printf("[all f_c = 1]            vmin = %+.6e   BFB = %s\n",
                res.vmin, res.bfb ? "YES" : "NO");
  }

  // ---- 2. negative-dominant f -> NOT BFB ----
  {
    std::vector<Real> f(nc, 0.05);
    f[nc - 1] = -1.0;                       // strongly negative on the highest channel
    auto res = check_bfb(ch, f);
    std::printf("[f_last = -1, rest 0.05] vmin = %+.6e   BFB = %s (expect NO)\n",
                res.vmin, res.bfb ? "YES" : "NO");
  }

  // ---- 3. the Sigma(108) locking couplings (docs/locking_couplings.md) -> BFB ----
  {
    std::vector<Real> f = {0.2782, 0.0930, 0.1029, 0.1256, 0.0071, 0.0101, 0.0443, 0.2148, 0.1240};
    auto res = check_bfb(ch, f);
    std::printf("[Sigma(108) locking f]   vmin = %+.6e   BFB = %s (expect YES)\n",
                res.vmin, res.bfb ? "YES" : "NO");
  }

  // ---- 4. a slice of the BFB region: vary (f0, f_last) with the rest held at the
  //         locking values, to map where copositivity breaks. ----
  std::printf("\n-- BFB region slice: f = locking, but f[0]=x and f[last]=y --\n");
  std::vector<Real> base = {0.2782, 0.0930, 0.1029, 0.1256, 0.0071, 0.0101, 0.0443, 0.2148, 0.1240};
  std::printf("        y\\x");
  const Real xs[] = {-0.3, -0.1, 0.0, 0.1, 0.3, 0.6};
  const Real ys[] = {-0.6, -0.3, -0.1, 0.0, 0.1, 0.3};
  for (Real x : xs) std::printf("  %+5.2f", x);
  std::printf("\n");
  for (Real y : ys) {
    std::printf("  %+9.2f", y);
    for (Real x : xs) {
      std::vector<Real> f = base; f[0] = x; f[nc - 1] = y;
      auto res = check_bfb(ch, f, /*nrestart=*/24, /*nstep=*/100);
      std::printf("  %s   ", res.bfb ? " . " : " X ");   // '.' BFB, 'X' unbounded
    }
    std::printf("\n");
  }
  std::printf("(  .  = bounded below;   X  = unbounded direction found )\n");

  return 0;
}
#endif
