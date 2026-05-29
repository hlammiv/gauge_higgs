#pragma once
// Vacuum alignment: given a discrete subgroup H (as a list of SU(N) elements) and a
// Higgs irrep R, build the H-singlet VEV intrinsically (image of the group-averaging
// projector), then find quartic couplings f_c (for the multi-invariant potential) that
// make that VEV a STABLE minimum -- i.e. couplings that lock the gauge group to H.
// Validation: the Hessian must have exactly dim(SU(N))=N^2-1 zero modes (the broken-
// generator Goldstones, all eaten) and the remaining massive modes positive and grouped
// into H-irrep multiplets (draft Observation 2). See [[research-program]].
#include "action/scalar_invariants.hpp"
#include "rep/finite_subgroups.hpp"
#include <random>
#include <algorithm>

namespace gh {

// H-singlet VEV: a unit vector fixed by D^(R)(h) for all h in H. multiplicity = dim of
// the fixed subspace (= rank of the averaging projector P = (1/|H|) sum_h D^(R)(h)).
template <int N>
struct SingletVEV { DVec vev; int multiplicity; Real invariance; };

template <int N>
SingletVEV<N> singlet_vev(const Representation<N>& rep, const std::vector<Cmat<N>>& H) {
  const int d = rep.d;
  DMat P(d, d);
  for (const auto& h : H) P = P + rep.rep_matrix(h);
  P = P * Complex(1.0 / H.size(), 0);
  // rank of P (= multiplicity): orthonormalize images of basis vectors.
  std::vector<DVec> basis;
  for (int j = 0; j < d && (int)basis.size() < d; ++j) {
    DVec e(d); e(j) = Complex(1, 0);
    DVec v = P * e;
    for (const auto& b : basis) { Complex ov = dot(b, v); for (int k = 0; k < d; ++k) v(k) -= ov * b(k); }
    Real nrm = std::sqrt(v.norm2());
    if (nrm > 1e-8) { for (int k = 0; k < d; ++k) v(k) *= Complex(1.0 / nrm, 0); basis.push_back(v); }
  }
  SingletVEV<N> out;
  out.multiplicity = (int)basis.size();
  out.vev = basis.empty() ? DVec(d) : basis[0];
  // invariance: max_h ||D(h) vev - vev||
  Real worst = 0;
  for (const auto& h : H) {
    DVec r = rep.rotate(h, out.vev);
    Real e = 0; for (int k = 0; k < d; ++k) e += std::norm(r(k) - out.vev(k));
    worst = std::max(worst, std::sqrt(e));
  }
  out.invariance = worst;
  return out;
}

// Real gradient at phi. For a REAL rep the physical dof are the d real components
// (length d, dV/dRe phi_k = 2Re(dV/dphibar_k)); for a complex rep all 2d (Re,Im) dof.
template <int N>
std::vector<Real> real_grad(const CasimirChannels<N>& ch, const DVec& phi,
                            const std::vector<Real>& f, Real mu2, bool real_rep) {
  DVec g = ch.dV_dphibar(phi, f, mu2);
  if (real_rep) { std::vector<Real> r(ch.d); for (int k = 0; k < ch.d; ++k) r[k] = 2 * g(k).real(); return r; }
  std::vector<Real> r(2 * ch.d);
  for (int k = 0; k < ch.d; ++k) { r[2 * k] = 2 * g(k).real(); r[2 * k + 1] = 2 * g(k).imag(); }
  return r;
}

// Real Hessian (n x n, n = d for real rep else 2d) at phi via finite difference.
template <int N>
std::vector<Real> real_hessian(const CasimirChannels<N>& ch, const DVec& phi,
                               const std::vector<Real>& f, Real mu2, bool real_rep) {
  const int n = real_rep ? ch.d : 2 * ch.d;
  std::vector<Real> Hm(static_cast<std::size_t>(n) * n, 0.0);
  const Real eps = 1e-5;
  DVec p = phi;
  for (int j = 0; j < n; ++j) {
    const int k = real_rep ? j : j / 2; const bool im = (!real_rep) && (j % 2);
    const Complex o = p(k);
    p(k) = o + (im ? Complex(0, eps) : Complex(eps, 0)); auto gp = real_grad<N>(ch, p, f, mu2, real_rep);
    p(k) = o - (im ? Complex(0, eps) : Complex(eps, 0)); auto gm = real_grad<N>(ch, p, f, mu2, real_rep);
    p(k) = o;
    for (int i = 0; i < n; ++i) Hm[i * n + j] = (gp[i] - gm[i]) / (2 * eps);
  }
  return Hm;
}

// Radial-stationarity: set mu2 so the gradient along phi0 vanishes (phi0 unit).
template <int N>
Real radial_mu2(const CasimirChannels<N>& ch, const DVec& phi0, const std::vector<Real>& f) {
  // grad = -mu2 phi0 + sum f_c 2 P_c(M)phi0 ; component along phi0 = 0 => mu2 = Re<phi0, quartic-grad>
  std::vector<Real> fzero(f.size(), 0.0);
  DVec q = ch.dV_dphibar(phi0, f, 0.0);  // quartic part only (mu2=0)
  return dot(phi0, q).real() / phi0.norm2();
}

struct AlignResult { std::vector<Real> f; Real mu2; std::vector<Real> eig;
                     int nzero, nneg; Real min_massive; bool stable; };

// Search couplings f (simplex, f_c>=0) making phi0 a STABLE minimum, i.e. the Hessian is
// positive semidefinite (no negative modes): the zero modes are the Goldstones (broken
// generators, eaten; +1 global U(1) for complex reps) and the rest are positive masses.
// Hessian is LINEAR in f, so per-channel Hessians + radial mu2 are precomputed ONCE.
template <int N>
AlignResult find_stable_couplings(const CasimirChannels<N>& ch, const DVec& phi0,
                                  bool real_rep, int nsamples = 800, std::uint64_t seed = 1) {
  const int nc = ch.n_channels();
  const int n = real_rep ? ch.d : 2 * ch.d;
  std::vector<std::vector<Real>> Hc(nc);    // Hessian of V_c at phi0 (mu2=0)
  std::vector<Real> mu2c(nc);               // radial mu2 from channel c alone
  // Per-channel Hessians are independent -> parallelize (the heavy step, esp. d=64 Sigma(648)).
  #pragma omp parallel for schedule(dynamic)
  for (int c = 0; c < nc; ++c) {
    std::vector<Real> e(nc, 0.0); e[c] = 1.0;
    Hc[c] = real_hessian<N>(ch, phi0, e, 0.0, real_rep);
    mu2c[c] = radial_mu2<N>(ch, phi0, e);
  }
  std::mt19937_64 rng(seed);
  std::uniform_real_distribution<Real> U(0.0, 1.0);
  AlignResult best; best.stable = false; best.min_massive = -1e30;
  std::vector<Real> Hm(static_cast<std::size_t>(n) * n);
  const Real zt = 1e-5;  // zero-mode tolerance
  for (int s = 0; s < nsamples; ++s) {
    std::vector<Real> f(nc);
    Real sum = 0; for (auto& x : f) { x = U(rng); sum += x; } for (auto& x : f) x /= sum;
    Real mu2 = 0; for (int c = 0; c < nc; ++c) mu2 += f[c] * mu2c[c];
    std::fill(Hm.begin(), Hm.end(), 0.0);
    for (int c = 0; c < nc; ++c) { const Real fc = f[c]; for (std::size_t i = 0; i < Hm.size(); ++i) Hm[i] += fc * Hc[c][i]; }
    for (int i = 0; i < n; ++i) Hm[i * n + i] -= 2.0 * mu2;   // Hessian of -mu2|phi|^2
    std::vector<Real> ev = eig_sym(Hm, n);
    std::sort(ev.begin(), ev.end());
    int nneg = 0, nzero = 0; Real min_pos = 1e30;
    for (Real e : ev) { if (e < -zt) ++nneg; else if (e < zt) ++nzero; else min_pos = std::min(min_pos, e); }
    bool stable = (nneg == 0) && (min_pos < 1e29);
    // prefer stable; among those maximize the mass gap (min positive eigenvalue)
    bool better = stable && (!best.stable || min_pos > best.min_massive);
    if (better || (!best.stable && !stable && -nneg > -best.nneg)) {  // fallback: fewest negatives
      best.f = f; best.mu2 = mu2; best.eig = ev; best.nzero = nzero; best.nneg = nneg;
      best.min_massive = (min_pos < 1e29 ? min_pos : 0.0); best.stable = stable;
    }
  }
  return best;
}

}  // namespace gh
