// Tree-level mean-field strong-weak coupling transition for SU(N) gauge theory
// frozen to a finite subgroup H (or kept continuous). Reproduces the lattice
// notes "Mean field method" + "tree level strong-weak coupling transition" and
// extends the criticals (beta_c, alpha_B, v_B) to SU(3) and the Sigma(108/216/
// 648/1080) crystal-like subgroups. Pure group theory + 1-D root finding; NO HMC.
//
//   w(alpha) = ln < exp( (alpha/N) Re Tr_fund U ) >            ("single-link" log)
//      Haar over continuous G, uniform over the |H| elements for discrete H.
//   Closed forms used as cross-checks:
//      w_SU(2)(a) = ln( 2 I_1(a)/a )
//      w_BT(a)    = ln( 3 + 8 cosh(a/2) + cosh(a) ) - ln 12
//
// Saddle point (TEMPORAL gauge, see notes):
//   N_P/N_l   = (D-2)/2                                        (plaquettes per link)
//   E_P(v)    = 1 - v^4 + (1-v^2) * 2/(D-2)                    (mean-field plaq energy)
//   v_B       = w'(alpha)                                      (mean-link saddle)
//   alpha     = -beta (N_P/N_l) E_P'(v_B)                      (stationarity)
//   s(v_B)    = w(alpha) - alpha w'(alpha)                     (Legendre transform of w)
//   F_P       = E_P(v_B) - (1/beta)(N_l/N_P) s(v_B)            (free energy per plaquette)
// beta_c is the first-order point where the nontrivial branch's F_P equals the
// trivial (alpha=0, v_B=0) value F_P^triv = E_P(0) = 1 + 2/(D-2).
//
// Build:
//   g++ -std=c++20 -O3 -march=native -fopenmp -Isrc -o build/meanfield src/meanfield.cpp
//   ./build/meanfield
#include "rep/finite_subgroups.hpp"
#include <cmath>
#include <cstdio>
#include <functional>
#include <string>
#include <vector>

namespace gh {

// ----------------------------------------------------------------------------
// An "ensemble" is a finite weighted set of (weight, ReTr_fund) nodes plus N.
// Both a discrete group (uniform weights) and a continuous-group quadrature on
// the maximal torus (Weyl-measure weights) reduce to this, so w/w'/w'' are all
// cheap weighted sums. The weights need not be normalized; w divides by sum w_k.
// ----------------------------------------------------------------------------
struct Ensemble {
  int N = 2;
  std::vector<Real> wt;     // weights
  std::vector<Real> retr;   // Re Tr_fund of each element / torus node
  Real wsum = 0.0;

  void finalize() { wsum = 0.0; for (Real w : wt) wsum += w; }

  // w(alpha) = ln( sum_k w_k exp((alpha/N) t_k) / sum_k w_k ).  Log-sum-exp safe.
  Real w(Real a) const {
    const Real c = a / N;
    Real mx = -1e300;
    for (Real t : retr) { Real e = c * t; if (e > mx) mx = e; }
    Real s = 0.0;
    for (std::size_t k = 0; k < retr.size(); ++k) s += wt[k] * std::exp(c * retr[k] - mx);
    return mx + std::log(s / wsum);
  }
  // w'(alpha) = (1/N) <t> under the tilted weights = the mean link v_B(alpha).
  Real wp(Real a) const {
    const Real c = a / N;
    Real mx = -1e300;
    for (Real t : retr) { Real e = c * t; if (e > mx) mx = e; }
    Real s = 0.0, st = 0.0;
    for (std::size_t k = 0; k < retr.size(); ++k) {
      Real p = wt[k] * std::exp(c * retr[k] - mx);
      s += p; st += p * retr[k];
    }
    return (st / s) / N;
  }
};

// Discrete subgroup H: uniform weights, ReTr of the fundamental N x N matrices.
template <int N>
Ensemble ensemble_discrete(const std::vector<Cmat<N>>& H) {
  Ensemble e; e.N = N;
  e.wt.reserve(H.size()); e.retr.reserve(H.size());
  for (const auto& U : H) { e.wt.push_back(1.0); e.retr.push_back(U.trace().real()); }
  e.finalize();
  return e;
}

// SU(2) continuous: torus theta in (0,pi), Haar weight (2/pi) sin^2 theta,
// ReTr = 2 cos theta. (Closed form ln(2 I_1(a)/a).)
inline Ensemble ensemble_su2(int M = 4000) {
  Ensemble e; e.N = 2;
  for (int i = 0; i < M; ++i) {
    Real th = (i + 0.5) / M * kPi;
    Real s = std::sin(th);
    e.wt.push_back(s * s);
    e.retr.push_back(2.0 * std::cos(th));
  }
  e.finalize();
  return e;
}

// SU(3) continuous: maximal torus diag(e^{i t1}, e^{i t2}, e^{i t3}), t3=-t1-t2,
// Weyl measure prod_{i<j} 4 sin^2((t_i - t_j)/2), ReTr = sum cos t_i.
inline Ensemble ensemble_su3(int M = 360) {
  Ensemble e; e.N = 3;
  const Real dt = 2.0 * kPi / M;
  for (int i = 0; i < M; ++i) {
    Real t1 = (i + 0.5) * dt;
    for (int j = 0; j < M; ++j) {
      Real t2 = (j + 0.5) * dt;
      Real t3 = -t1 - t2;
      Real s12 = std::sin((t1 - t2) / 2), s13 = std::sin((t1 - t3) / 2), s23 = std::sin((t2 - t3) / 2);
      Real meas = (4 * s12 * s12) * (4 * s13 * s13) * (4 * s23 * s23);
      e.wt.push_back(meas);
      e.retr.push_back(std::cos(t1) + std::cos(t2) + std::cos(t3));
    }
  }
  e.finalize();
  return e;
}

// ----------------------------------------------------------------------------
// Saddle-point machinery, temporal gauge.
// ----------------------------------------------------------------------------
struct PlaqGeom {
  int D;
  Real ratio;  // N_P / N_l = (D-2)/2
  explicit PlaqGeom(int d) : D(d), ratio((d - 2) / 2.0) {}
  Real EP(Real v)  const { return 1.0 - v * v * v * v + (1.0 - v * v) * 2.0 / (D - 2); }
  Real EPp(Real v) const { return -4.0 * v * v * v + (-2.0 * v) * 2.0 / (D - 2); }  // dE_P/dv
  Real FP_triv()   const { return EP(0.0); }  // s(0)=0, v_B(0)=0
};

struct Branch { Real alpha, vB, beta, s, FP; };

// One point of the nontrivial branch as a function of alpha.
inline Branch branch_at(const Ensemble& e, const PlaqGeom& g, Real alpha) {
  Branch b; b.alpha = alpha;
  b.vB    = e.wp(alpha);
  b.beta  = -alpha / (g.ratio * g.EPp(b.vB));
  b.s     = e.w(alpha) - alpha * b.vB;
  b.FP    = g.EP(b.vB) - (1.0 / b.beta) * (1.0 / g.ratio) * b.s;
  return b;
}

struct Critical { Real beta_c, alpha_B, vB; bool found; };

// First-order critical: scan alpha>0 for the point where F_P(branch) crosses
// F_P^triv, then Newton-refine the crossing in alpha.
inline Critical find_critical(const Ensemble& e, const PlaqGeom& g,
                              Real amax, Real da = 5e-4) {
  const Real Ftriv = g.FP_triv();
  Real prev_a = 0.0, prev_d = 0.0; bool have = false;
  Critical out{0, 0, 0, false};
  for (Real a = da; a < amax; a += da) {
    Real d = branch_at(e, g, a).FP - Ftriv;
    if (have && ((prev_d < 0 && d >= 0) || (prev_d > 0 && d <= 0))) {
      // bisection refine on [prev_a, a]
      Real lo = prev_a, hi = a, dlo = prev_d;
      for (int it = 0; it < 80; ++it) {
        Real mid = 0.5 * (lo + hi);
        Real dm = branch_at(e, g, mid).FP - Ftriv;
        if ((dlo < 0 && dm < 0) || (dlo > 0 && dm > 0)) { lo = mid; dlo = dm; } else hi = mid;
      }
      Real ac = 0.5 * (lo + hi);
      Branch b = branch_at(e, g, ac);
      out = {b.beta, b.alpha, b.vB, true};
      return out;
    }
    prev_a = a; prev_d = d; have = true;
  }
  return out;
}

}  // namespace gh

using namespace gh;

namespace {
struct Row { std::string name; const Ensemble* e; Real amax; };

void print_table(const std::vector<Row>& rows) {
  for (int D : {3, 4}) {
    PlaqGeom g(D);
    std::printf("\n  D = %d   (N_P/N_l = %.3f,  E_P^triv = %.4f)\n", D, g.ratio, g.FP_triv());
    std::printf("  %-14s  %10s  %10s  %10s   %12s\n",
                "group", "beta_c", "alpha_B", "v_B", "1-v_B");
    std::printf("  %s\n", std::string(64, '-').c_str());
    for (const auto& r : rows) {
      Critical c = find_critical(*r.e, g, r.amax);
      if (c.found)
        std::printf("  %-14s  %10.5f  %10.5f  %10.5f   %12.3e\n",
                    r.name.c_str(), c.beta_c, c.alpha_B, c.vB, 1.0 - c.vB);
      else
        std::printf("  %-14s  %10s\n", r.name.c_str(), "(no transition)");
    }
  }
}
}  // namespace

// The test (test_meanfield.cpp) #includes this file to reuse the Ensemble +
// saddle machinery, so guard the driver main() out of that translation unit.
#ifndef MEANFIELD_NO_MAIN
int main() {
  // --- SU(2) sector ensembles ---
  Ensemble su2 = ensemble_su2();
  Ensemble bt  = ensemble_discrete<2>(close_group<2>(gens_2T()));
  Ensemble bo  = ensemble_discrete<2>(close_group<2>(gens_2O()));
  Ensemble bi  = ensemble_discrete<2>(close_group<2>(gens_2I()));

  // --- SU(3) sector ensembles ---
  Ensemble su3   = ensemble_su3();
  Ensemble s108  = ensemble_discrete<3>(close_group<3>(gens_Sigma108()));
  Ensemble s216  = ensemble_discrete<3>(close_group<3>(gens_Sigma216()));
  Ensemble s648  = ensemble_discrete<3>(close_group<3>(gens_Sigma648()));
  Ensemble s1080 = ensemble_discrete<3>(close_group<3>(gens_Sigma1080()));

  std::printf("============================================================\n");
  std::printf(" Mean-field strong-weak criticals  (TEMPORAL gauge)\n");
  std::printf(" w(a)=ln<exp((a/N)ReTr U)>,  saddle E_P(v)=1-v^4+(1-v^2)2/(D-2)\n");
  std::printf("============================================================\n");
  std::printf(" |BT|=%zu |BO|=%zu |BI|=%zu | |S108|=%zu |S216|=%zu |S648|=%zu |S1080|=%zu\n",
              bt.retr.size(), bo.retr.size(), bi.retr.size(),
              s108.retr.size(), s216.retr.size(), s648.retr.size(), s1080.retr.size());

  std::printf("\n################  SU(2) SECTOR  ################\n");
  print_table({{"SU(2) cont", &su2, 25.0},
               {"2T (BT)",    &bt,  40.0},
               {"2O (BO)",    &bo,  40.0},
               {"2I (BI)",    &bi,  40.0}});

  std::printf("\n################  SU(3) SECTOR  ################\n");
  print_table({{"SU(3) cont", &su3,  40.0},
               {"Sigma(108)", &s108, 60.0},
               {"Sigma(216)", &s216, 60.0},
               {"Sigma(648)", &s648, 60.0},
               {"Sigma(1080)",&s1080,60.0}});

  // Frozenness signature: 1 - v_B(alpha) = 1 - w'(alpha) as alpha grows large.
  // Continuous G: power-law (SU(2): 1-v_B ~ 3/(2 alpha)). Discrete H: exponential,
  // governed by the gap between the largest ReTr (= N, identity) and the next
  // distinct value -- the link freezes onto the group element exponentially fast.
  std::printf("\n#### Discrete frozenness: 1 - v_B(alpha) = 1 - w'(alpha) ####\n");
  std::printf("   continuous -> power-law;  discrete -> exponential (link freezes).\n");
  std::printf("   %8s  %12s  %12s  %12s  %12s  %12s\n",
              "alpha", "SU(2)", "2T(BT)", "2I(BI)", "SU(3)", "S108");
  for (Real a : {4.0, 8.0, 12.0, 16.0, 20.0, 30.0}) {
    std::printf("   %8.1f  %12.4e  %12.4e  %12.4e  %12.4e  %12.4e\n",
                a, 1.0 - su2.wp(a), 1.0 - bt.wp(a), 1.0 - bi.wp(a),
                1.0 - su3.wp(a), 1.0 - s108.wp(a));
  }
  std::printf("   (SU(2) ~ 3/(2a): a=30 -> 4.96e-2 ~ 0.05;  BT ~ e^{-a/2}: a=20->30 drops ~150x)\n");
  return 0;
}
#endif  // MEANFIELD_NO_MAIN
