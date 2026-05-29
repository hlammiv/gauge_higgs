#pragma once
// Observables and time-series statistics for validation/physics.
//   - Polyakov loop (gauge sector)
//   - Stats: mean, variance, susceptibility chi = V*Var, Binder cumulant U4,
//     and a binned (block-jackknife) error that accounts for autocorrelation.
// Plaquette lives in gauge_wilson.hpp; L_phi / L_link in scalar_higgs.hpp.
#include "core/fields.hpp"
#include "rep/representation.hpp"
#include <vector>
#include <cmath>

namespace gh {

// Average Polyakov loop in direction tdir (default = last dim = "time"):
//   P = < (1/N) Re Tr [ prod_{t} U_tdir(x,t) ] >  averaged over the transverse sites.
// Screened by matter (not a strict order parameter for matter theories) but standard.
template <int D, int N>
Real polyakov_loop(const GaugeField<D, N>& U, int tdir = D - 1) {
  const Lattice<D>& lat = *U.lat;
  const int Lt = lat.L[tdir];
  Real sum = 0.0;
  std::int64_t count = 0;
  std::array<int, D> x{};
  for (std::int64_t s = 0; s < lat.vol; ++s) {
    lat.coords(s, x);
    if (x[tdir] != 0) continue;                 // one loop per transverse site
    Cmat<N> P = U(s, tdir);
    std::int64_t cur = lat.neighbor_fwd(s, tdir);
    for (int t = 1; t < Lt; ++t) { P = P * U(cur, tdir); cur = lat.neighbor_fwd(cur, tdir); }
    sum += P.trace().real();
    ++count;
  }
  return sum / (static_cast<Real>(count) * N);
}

// Path-ordered product of fundamental links around an Ra x Rb rectangle in the
// (mu,nu) plane based at site x (closes back to x). Ra=Rb=1 gives the plaquette.
template <int D, int N>
Cmat<N> wilson_rect(const GaugeField<D, N>& U, std::int64_t x, int mu, int nu, int Ra, int Rb) {
  const Lattice<D>& lat = *U.lat;
  Cmat<N> W = Cmat<N>::identity();
  std::int64_t s = x;
  for (int i = 0; i < Ra; ++i) { W = W * U(s, mu);          s = lat.neighbor_fwd(s, mu); }
  for (int j = 0; j < Rb; ++j) { W = W * U(s, nu);          s = lat.neighbor_fwd(s, nu); }
  for (int i = 0; i < Ra; ++i) { s = lat.neighbor_bwd(s, mu); W = W * U(s, mu).dagger(); }
  for (int j = 0; j < Rb; ++j) { s = lat.neighbor_bwd(s, nu); W = W * U(s, nu).dagger(); }
  return W;
}

// Average fundamental Wilson loop W(Ra,Rb) = <(1/N) Re Tr U_loop> over sites and mu<nu planes.
template <int D, int N>
Real wilson_loop_fund(const GaugeField<D, N>& U, int Ra, int Rb) {
  const Lattice<D>& lat = *U.lat;
  Real sum = 0.0; std::int64_t count = 0;
  for (std::int64_t x = 0; x < lat.vol; ++x)
    for (int mu = 0; mu < D; ++mu)
      for (int nu = mu + 1; nu < D; ++nu) {
        sum += wilson_rect<D, N>(U, x, mu, nu, Ra, Rb).trace().real();
        ++count;
      }
  return sum / (static_cast<Real>(count) * N);
}

// Average Wilson loop traced in an arbitrary probe irrep R':
//   W_R'(Ra,Rb) = < (1/d_R') Re Tr D^(R')(U_loop) > = < (1/d_R') Re chi_R'(U_loop) >.
// The Ra,Rb -> large area/perimeter behaviour over a set of probe reps R' fingerprints
// the residual symmetry: reps screened by the H-condensate -> perimeter law (-> const),
// unscreened -> area law. (Measurement only; rep_matrix is fine for occasional use.)
template <int D, int N>
Real wilson_loop_rep(const GaugeField<D, N>& U, const Representation<N>& rep, int Ra, int Rb) {
  const Lattice<D>& lat = *U.lat;
  Real sum = 0.0; std::int64_t count = 0;
  for (std::int64_t x = 0; x < lat.vol; ++x)
    for (int mu = 0; mu < D; ++mu)
      for (int nu = mu + 1; nu < D; ++nu) {
        const Cmat<N> loop = wilson_rect<D, N>(U, x, mu, nu, Ra, Rb);
        sum += rep.rep_matrix(loop).trace().real();
        ++count;
      }
  return sum / (static_cast<Real>(count) * rep.d);
}

// Polyakov loop traced in an arbitrary probe irrep R' (order parameter for the residual
// center symmetry): < (1/d_R') Re Tr D^(R')( prod_t U_tdir ) >.
template <int D, int N>
Real polyakov_loop_rep(const GaugeField<D, N>& U, const Representation<N>& rep, int tdir = D - 1) {
  const Lattice<D>& lat = *U.lat;
  const int Lt = lat.L[tdir];
  Real sum = 0.0; std::int64_t count = 0;
  std::array<int, D> xc{};
  for (std::int64_t s = 0; s < lat.vol; ++s) {
    lat.coords(s, xc);
    if (xc[tdir] != 0) continue;
    Cmat<N> P = U(s, tdir);
    std::int64_t cur = lat.neighbor_fwd(s, tdir);
    for (int t = 1; t < Lt; ++t) { P = P * U(cur, tdir); cur = lat.neighbor_fwd(cur, tdir); }
    sum += rep.rep_matrix(P).trace().real();
    ++count;
  }
  return sum / (static_cast<Real>(count) * rep.d);
}

// Time-series accumulator for an intensive observable O measured once per config.
struct Stats {
  std::vector<double> x;
  void add(double v) { x.push_back(v); }
  int  n() const { return static_cast<int>(x.size()); }

  double mean() const {
    double s = 0.0; for (double v : x) s += v; return x.empty() ? 0.0 : s / x.size();
  }
  double var() const {  // <O^2> - <O>^2  (population)
    double m = mean(), s2 = 0.0;
    for (double v : x) s2 += v * v;
    s2 = x.empty() ? 0.0 : s2 / x.size();
    return s2 - m * m;
  }
  // Susceptibility chi = V * Var(O).
  double susceptibility(double volume) const { return volume * var(); }
  // Binder cumulant U4 = 1 - <O^4> / (3 <O^2>^2), centered on the mean.
  double binder() const {
    double m = mean(), m2 = 0.0, m4 = 0.0;
    for (double v : x) { double d = v - m; m2 += d * d; m4 += d * d * d * d; }
    if (x.empty()) return 0.0;
    m2 /= x.size(); m4 /= x.size();
    return (m2 > 0.0) ? 1.0 - m4 / (3.0 * m2 * m2) : 0.0;
  }
  // Block-jackknife error on the mean using nbins blocks (autocorrelation-aware).
  double binned_error(int nbins = 16) const {
    const int N = n();
    if (N < 2 * nbins) nbins = std::max(2, N / 2);
    if (N < 2) return 0.0;
    const int bs = N / nbins;
    std::vector<double> bm;
    for (int b = 0; b < nbins; ++b) {
      double s = 0.0;
      for (int i = b * bs; i < (b + 1) * bs; ++i) s += x[i];
      bm.push_back(s / bs);
    }
    double mb = 0.0; for (double v : bm) mb += v; mb /= bm.size();
    double v = 0.0; for (double b : bm) v += (b - mb) * (b - mb);
    v /= (bm.size() - 1);
    return std::sqrt(v / bm.size());
  }
};

}  // namespace gh
