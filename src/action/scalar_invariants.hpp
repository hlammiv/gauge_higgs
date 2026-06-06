#pragma once
// GENERAL multi-invariant quartic scalar potential, built from any Representation.
// The independent quartic G-invariants are the projections of the bilinear M=phi phi^dag
// onto the irreducible channels of R (x) Rbar (the d x d matrix space), which are the
// eigenspaces of the adjoint quadratic Casimir superoperator
//     Chat(M) = sum_a [T^a_R, [T^a_R, M]]   (eigenvalue C2(channel)).
// V[phi] = -mu2 (phi^dag phi) + sum_c f[c] V_c[phi],  V_c[phi] = Tr(M^dag P_c M),
// with P_c the Casimir eigenprojector (Lagrange polynomial in Chat). Tuning {f[c]}
// aligns the VEV onto a chosen H-stratum (vacuum alignment). Group-/rep-agnostic:
// the only input is the rep's generators T^a_R (which rep_general builds for any irrep).
// See research-program / discrete-observables-program; draft V_{rho'} construction.
#include "rep/representation.hpp"
#include "action/scalar_higgs.hpp"  // OnsitePotential base
#include <random>
#include <algorithm>

namespace gh {

template <int N>
struct CasimirChannels {
  const Representation<N>* rep = nullptr;
  int d = 0;
  std::vector<Real> lambda;   // distinct adjoint-Casimir eigenvalues = channel C2's

  explicit CasimirChannels(const Representation<N>& R) : rep(&R), d(R.d) {
    lambda = find_channels();
  }
  int n_channels() const { return static_cast<int>(lambda.size()); }

  // Chat(M) = sum_a [T^a,[T^a,M]].
  DMat apply_casimir(const DMat& M) const {
    DMat out(d, d);
    for (int a = 0; a < n_gen<N>(); ++a) {
      const DMat& T = rep->T[a];
      DMat c1 = T * M - M * T;          // [T^a, M]
      out = out + (T * c1 - c1 * T);    // [T^a, [T^a, M]]
    }
    return out;
  }

  // Project M onto channel c via the Lagrange polynomial in Chat.
  DMat apply_proj(int c, const DMat& M) const {
    DMat X = M;
    for (int cp = 0; cp < n_channels(); ++cp) {
      if (cp == c) continue;
      DMat CX = apply_casimir(X);
      X = (CX - X * Complex(lambda[cp], 0)) * Complex(1.0 / (lambda[c] - lambda[cp]), 0);
    }
    return X;
  }

  static DMat outer(const DVec& phi) {  // M_{ab} = phi_a conj(phi_b)
    const int d = phi.size();
    DMat M(d, d);
    for (int a = 0; a < d; ++a)
      for (int b = 0; b < d; ++b) M(a, b) = phi(a) * std::conj(phi(b));
    return M;
  }

  // V_c[phi] = Tr(M^dag P_c M) = <M, P_c M> = ||P_c M||^2 >= 0.
  Real channel_value(int c, const DVec& phi) const {
    DMat M = outer(phi);
    return reTrDagProd(M, apply_proj(c, M));
  }

  // Total V = -mu2 phi^dag phi + sum_c f[c] V_c. (Skip zero couplings -- big speedup when
  // few channels are active, e.g. the per-channel Hessian precompute uses f = e_c.)
  Real value(const DVec& phi, const std::vector<Real>& f, Real mu2) const {
    Real v = -mu2 * phi.norm2();
    DMat M = outer(phi);
    for (int c = 0; c < n_channels(); ++c) if (f[c] != 0.0) v += f[c] * reTrDagProd(M, apply_proj(c, M));
    return v;
  }

  // dV/dphibar (complex d-vector): -mu2 phi + sum_c f[c] * 2 P_c(M) phi.
  // (derivation: dV_c/dphibar = 2 P_c(M) phi, since P_c(M) is Hermitian.)
  DVec dV_dphibar(const DVec& phi, const std::vector<Real>& f, Real mu2) const {
    DMat M = outer(phi);
    DVec g(d);
    for (int e = 0; e < d; ++e) g(e) = Complex(-mu2, 0) * phi(e);
    for (int c = 0; c < n_channels(); ++c) {
      if (f[c] == 0.0) continue;           // skip inactive channels (big speedup)
      DMat PM = apply_proj(c, M);          // Hermitian
      DVec PMphi = PM * phi;
      for (int e = 0; e < d; ++e) g(e) += Complex(2.0 * f[c], 0) * PMphi(e);
    }
    return g;
  }

  // ---- phi-INDEPENDENT precompute (HMC hot path) ------------------------------
  // dV_dphibar's per-channel projector P_c is the SAME Lagrange polynomial in the
  // Casimir superoperator for every phi -- only M=phi phi^dag changes. The whole
  // sum_c f[c]*2*P_c is a fixed LINEAR superoperator A on the d x d matrix space.
  // We materialize A once as a (d^2) x (d^2) matrix Asuper by feeding the EXACT
  // apply_proj/apply_casimir code path the basis matrices E_{ij} -- so the cache
  // reproduces apply_proj to within roundoff (it REORDERS floating-point ops, so it
  // is NOT bit-identical; it matches the reference to the algorithm's full RELATIVE
  // precision -- the degree-(n_ch-1) Lagrange projector is itself roundoff-limited at
  // large d, e.g. ~2e-8 relative at d=13). dV then costs one d^2 x d^2 matvec + one
  // d x d matvec per site instead of O(n_ch^2 * n_gen * d^3).
  // Layout: vec index = i*d + j for matrix entry (i,j); Asuper is (d^2) x (d^2).
  std::vector<Complex> build_combined_superop(const std::vector<Real>& f) const {
    const int dd = d * d;
    std::vector<Complex> Asuper(static_cast<std::size_t>(dd) * dd, Complex(0, 0));
    for (int i = 0; i < d; ++i)
      for (int j = 0; j < d; ++j) {
        DMat E(d, d); E(i, j) = Complex(1, 0);     // basis matrix E_{ij}
        DMat AE(d, d);                              // AE = sum_c f[c]*2*P_c(E_{ij})
        for (int c = 0; c < n_channels(); ++c) {
          if (f[c] == 0.0) continue;
          AE = AE + apply_proj(c, E) * Complex(2.0 * f[c], 0);
        }
        const int col = i * d + j;
        for (int p = 0; p < d; ++p)
          for (int q = 0; q < d; ++q)
            Asuper[static_cast<std::size_t>(p * d + q) * dd + col] = AE(p, q);
      }
    return Asuper;
  }

  // Fast dV_dphibar using a prebuilt combined superoperator (build_combined_superop).
  // g = -mu2 phi + A(M) phi, with vec(A(M)) = Asuper * vec(M), M=phi phi^dag.
  DVec dV_dphibar_cached(const DVec& phi, const std::vector<Complex>& Asuper, Real mu2) const {
    const int dd = d * d;
    // vec(M)_{i*d+j} = phi_i conj(phi_j)
    std::vector<Complex> vM(static_cast<std::size_t>(dd));
    for (int i = 0; i < d; ++i) {
      const Complex pi = phi(i);
      for (int j = 0; j < d; ++j) vM[static_cast<std::size_t>(i) * d + j] = pi * std::conj(phi(j));
    }
    // AM = reshape(Asuper * vec(M))  (d x d)
    DMat AM(d, d);
    for (int r = 0; r < dd; ++r) {
      const Complex* row = Asuper.data() + static_cast<std::size_t>(r) * dd;
      Complex acc(0, 0);
      for (int k = 0; k < dd; ++k) acc += row[k] * vM[k];
      AM.a[r] = acc;
    }
    // g = -mu2 phi + AM * phi
    DVec g(d);
    for (int p = 0; p < d; ++p) {
      Complex acc = Complex(-mu2, 0) * phi(p);
      const Complex* arow = AM.a.data() + static_cast<std::size_t>(p) * d;
      for (int q = 0; q < d; ++q) acc += arow[q] * phi(q);
      g(p) = acc;
    }
    return g;
  }

 private:
  // Distinct eigenvalues of Chat via Lanczos in the d x d matrix space (Hermitian
  // operator, real inner product Re Tr(A^dag B)); generic Hermitian start sees every
  // channel, so the tridiagonal's eigenvalues are exactly the distinct C2's.
  std::vector<Real> find_channels() const {
    std::mt19937_64 rng(20260529ULL);
    std::normal_distribution<Real> gd(0.0, 1.0);
    DMat r(d, d);
    for (int a = 0; a < d; ++a)
      for (int b = 0; b < d; ++b) r(a, b) = Complex(gd(rng), gd(rng));
    DMat v0 = r + r.dagger();                       // Hermitian, generic
    v0 = v0 * Complex(1.0 / std::sqrt(reTrDagProd(v0, v0)), 0);

    std::vector<DMat> V{v0};
    std::vector<Real> alpha, beta;
    for (int k = 0; k < 60; ++k) {
      DMat w = apply_casimir(V[k]);
      Real a = reTrDagProd(V[k], w); alpha.push_back(a);
      w = w - V[k] * Complex(a, 0);
      if (k > 0) w = w - V[k - 1] * Complex(beta[k - 1], 0);
      // full reorthogonalization, twice (DGKS) -- a single pass loses orthogonality at
      // larger matrix dimension and Lanczos then emits ghost eigenvalues (spurious channels).
      for (int pass = 0; pass < 2; ++pass)
        for (const auto& u : V) { Real ov = reTrDagProd(u, w); w = w - u * Complex(ov, 0); }
      Real b = std::sqrt(std::max(0.0, reTrDagProd(w, w)));
      if (b < 1e-6) break;
      beta.push_back(b);
      V.push_back(w * Complex(1.0 / b, 0));
    }
    const int m = static_cast<int>(alpha.size());
    std::vector<Real> Tri(static_cast<std::size_t>(m) * m, 0.0);
    for (int i = 0; i < m; ++i) Tri[i * m + i] = alpha[i];
    for (int i = 0; i + 1 < m; ++i) { Tri[i * m + (i + 1)] = beta[i]; Tri[(i + 1) * m + i] = beta[i]; }
    std::vector<Real> ev = eig_sym(Tri, m);
    std::sort(ev.begin(), ev.end());
    std::vector<Real> distinct;
    for (Real e : ev) if (distinct.empty() || std::fabs(e - distinct.back()) > 1e-4) distinct.push_back(e);
    return distinct;
  }
};

// Pluggable multi-invariant potential for the HMC: V = -mu2 phi^dag phi + sum_c f[c] V_c.
// Tuning {f[c]} (per-channel couplings) aligns the VEV onto a chosen discrete subgroup H.
template <int N>
struct MultiInvariantPotential : OnsitePotential<N> {
  const CasimirChannels<N>* ch = nullptr;
  std::vector<Real> f;
  Real mu2 = 0.0;
  // Combined phi-INDEPENDENT superoperator A = sum_c f[c]*2*P_c, materialized once
  // (d^2 x d^2). dV_dphibar is the dominant MD cost; this hoists its entire
  // projector/Casimir rebuild out of the per-site inner loop. See build_combined_superop.
  // NOTE: built from f at construction -- f is fixed for the potential's lifetime
  // (all call sites pass f to the ctor and never mutate it afterward).
  std::vector<Complex> Asuper;
  MultiInvariantPotential(const CasimirChannels<N>& c, std::vector<Real> ff, Real m)
      : ch(&c), f(std::move(ff)), mu2(m), Asuper(c.build_combined_superop(f)) {}
  Real value(const DVec& phi) const override { return ch->value(phi, f, mu2); }
  DVec dV_dphibar(const DVec& phi) const override { return ch->dV_dphibar_cached(phi, Asuper, mu2); }
};

}  // namespace gh
