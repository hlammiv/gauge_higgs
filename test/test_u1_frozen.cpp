// Tests for the FROZEN-LENGTH (|phi|=1) compact U(1) + charge-q Higgs local sampler.
//   1. von Mises heatbath sampler reproduces <cos> = I1(a)/I0(a), <sin> = 0.
//   2. Z_q link-jump is EXACTLY dB=0 (the matter hopping is invariant).
//   3. matter resultant R_x predicts the global action change for a single-site chi move.
//   4. gauge local Metropolis dS matches the global S(theta) change for a single-link move.
//   5. gauge staple G: Re[e^{i theta} G] == local_gauge_cos, and the overrelax reflection preserves it.
//   6. kappa=0 pure-gauge <plaquette> agrees with the validated U1HMC (same compact-U(1) measure).
//   7. multicanonical g(B)=c*B in the accept reproduces an unbiased run at kappa+c (bias enters every kernel).
#include "check.hpp"
#include "u1/u1_frozen.hpp"
#include "u1/u1.hpp"
#include <cmath>

using namespace gh;
using namespace gh::u1;

template <int D>
static std::array<int, D> cube(int n) { std::array<int, D> L{}; for (int mu = 0; mu < D; ++mu) L[mu] = n; return L; }

template <int D>
static Real frozen_action(const U1Frozen<D>& s) { return s.beta * s.A() - s.kappa * s.B(); }

// ---- 1. von Mises sampler statistics ----
static void test_vonmises() {
  Rng rng(20260611);
  for (Real a : {0.5, 2.0, 5.0}) {
    const int N = 400000; double sc = 0, ss = 0;
    for (int i = 0; i < N; ++i) { const Real ang = draw_vonmises(a, rng, Rng::key(0x777, i)); sc += std::cos(ang); ss += std::sin(ang); }
    const double ref = std::cyl_bessel_i(1, a) / std::cyl_bessel_i(0, a);
    char m[96];
    std::snprintf(m, sizeof m, "vonMises <cos> a=%.1f got %.4f exp %.4f", a, sc / N, ref); CHECK_CLOSE(sc / N, ref, 0.01, m);
    std::snprintf(m, sizeof m, "vonMises <sin> a=%.1f got %.4f exp 0", a, ss / N);          CHECK_CLOSE(ss / N, 0.0, 0.01, m);
  }
}

// ---- 2. Z_q link-jump exactly preserves B ----
template <int D>
static void test_zq_hop_dB0(int q) {
  U1Frozen<D> s(cube<D>(4), 7); s.beta = 1.3; s.kappa = 0.6; s.q = q; s.hot(0.8);
  const Real B0 = s.B(); Real worst = 0.0;
  for (std::int64_t x = 0; x < s.lat.vol; ++x)
    for (int mu = 0; mu < D; ++mu)
      for (int m = 1; m < q; ++m) {
        const std::size_t idx = x * D + mu; const Real o = s.th[idx];
        s.th[idx] = o + 2.0 * kPi * m / q; worst = std::max(worst, std::fabs(s.B() - B0)); s.th[idx] = o;
      }
  char msg[64]; std::snprintf(msg, sizeof msg, "Z_q hop dB=0 (q=%d)", q); CHECK_CLOSE(worst, 0.0, 1e-11, msg);
}

// ---- 3. matter resultant predicts the single-site action change ----
template <int D>
static void test_matter_local(int q) {
  U1Frozen<D> s(cube<D>(4), 13); s.beta = 1.1; s.kappa = 0.45; s.q = q; s.hot(0.7);
  const Real delta = 0.123; Real worst = 0.0;
  for (std::int64_t x : {std::int64_t(0), std::int64_t(7), s.lat.vol - 1}) {
    const Complex R = s.matter_resultant(x);
    const Real a = 2.0 * s.kappa * std::abs(R), alpha = std::arg(R);
    const Real pred = -a * (std::cos(s.chi[x] + delta - alpha) - std::cos(s.chi[x] - alpha));
    const Real S0 = frozen_action<D>(s); const Real old = s.chi[x]; s.chi[x] = old + delta;
    const Real glob = frozen_action<D>(s) - S0; s.chi[x] = old;
    worst = std::max(worst, std::fabs(pred - glob));
  }
  char msg[72]; std::snprintf(msg, sizeof msg, "matter resultant vs global dS (q=%d)", q); CHECK_CLOSE(worst, 0.0, 1e-9, msg);
}

// ---- 4. gauge local Metropolis dS matches the global action change ----
template <int D>
static void test_gauge_local(int q) {
  U1Frozen<D> s(cube<D>(4), 17); s.beta = 1.6; s.kappa = 0.5; s.q = q; s.hot(0.7);
  const Real delta = 0.137; Real worst = 0.0;
  for (std::int64_t x : {std::int64_t(0), std::int64_t(9), s.lat.vol - 1})
    for (int mu = 0; mu < D; ++mu) {
      const std::size_t idx = x * D + mu;
      const Real psi = s.chi[s.lat.neighbor_fwd(x, mu)] - s.chi[x]; const Real o = s.th[idx];
      const Real Lold = -s.beta * s.local_gauge_cos(x, mu) - 2.0 * s.kappa * std::cos(q * o + psi);
      const Real S0 = frozen_action<D>(s); s.th[idx] = o + delta;
      const Real Lnew = -s.beta * s.local_gauge_cos(x, mu) - 2.0 * s.kappa * std::cos(q * (o + delta) + psi);
      const Real glob = frozen_action<D>(s) - S0; s.th[idx] = o;
      worst = std::max(worst, std::fabs((Lnew - Lold) - glob));
    }
  char msg[72]; std::snprintf(msg, sizeof msg, "gauge local dS vs global (q=%d)", q); CHECK_CLOSE(worst, 0.0, 1e-9, msg);
}

// ---- 5. gauge staple G: Re[e^{i theta} G] == local_gauge_cos; OR reflection preserves it ----
template <int D>
static void test_gauge_staple(int q) {
  U1Frozen<D> s(cube<D>(4), 23); s.beta = 1.4; s.kappa = 0.4; s.q = q; s.hot(0.9);
  Real wid = 0.0, wor = 0.0;
  for (std::int64_t x : {std::int64_t(0), std::int64_t(4), std::int64_t(11), s.lat.vol - 1})
    for (int mu = 0; mu < D; ++mu) {
      const std::size_t idx = x * D + mu; const Real th0 = s.th[idx];
      const Complex G = s.gauge_staple(x, mu);
      wid = std::max(wid, std::fabs(U1Frozen<D>::plaq_cos_at(G, th0) - s.local_gauge_cos(x, mu)));
      const Real thr = -2.0 * std::arg(G) - th0;   // overrelaxation reflection
      wor = std::max(wor, std::fabs(U1Frozen<D>::plaq_cos_at(G, thr) - U1Frozen<D>::plaq_cos_at(G, th0)));
    }
  char msg[80];
  std::snprintf(msg, sizeof msg, "staple Re[e^{it}G]=local_gauge_cos (q=%d)", q); CHECK_CLOSE(wid, 0.0, 1e-10, msg);
  std::snprintf(msg, sizeof msg, "gauge OR reflection preserves plaq (q=%d)", q); CHECK_CLOSE(wor, 0.0, 1e-10, msg);
}

// ---- 6. kappa=0: frozen pure-gauge <plaq> agrees with U1HMC ----
template <int D>
static void test_kappa0_vs_hmc() {
  const Real beta = 2.0;
  U1Frozen<D> s(cube<D>(4), 101); s.beta = beta; s.kappa = 0.0; s.q = 2; s.zq_hop = false; s.hot(0.8);
  s.thermalize(300);
  double pf = 0; int nf = 0; for (int i = 0; i < 600; ++i) { s.sweep(); pf += s.avg_plaq(); ++nf; } pf /= nf;
  U1HMC<D> h(cube<D>(4), 202); h.beta = beta; h.kappa = 0.0; h.lambda = 0.5; h.q = 2; h.tau = 1.0; h.nmd = 20;
  h.hot(0.8); h.cold_phi(0.5); for (int i = 0; i < 200; ++i) h.trajectory();
  double ph = 0; int nh = 0; for (int i = 0; i < 600; ++i) { h.trajectory(); ph += avg_plaquette<D>(h.th, h.lat); ++nh; } ph /= nh;
  char msg[96]; std::snprintf(msg, sizeof msg, "kappa=0 <plaq> frozen %.4f vs HMC %.4f", pf, ph); CHECK_CLOSE(pf, ph, 0.02, msg);
}

// ---- 7. muca g(B)=c*B reproduces an unbiased run at kappa+c (bias enters every B-changing kernel) ----
//   metro=false: heatbath-proposal matter accept; metro=true: broad-proposal matter Metropolis (the muca path).
template <int D>
static void test_muca_linear(int q, bool metro) {
  const Real beta = 1.0, kappa = 0.20, c = 0.10;
  U1Frozen<D> ref(cube<D>(4), 303); ref.beta = beta; ref.kappa = kappa + c; ref.q = q; ref.hot(0.8);
  ref.thermalize(400);
  double lr = 0; int nr = 0; for (int i = 0; i < 1500; ++i) { ref.sweep(); lr += ref.link_energy(); ++nr; } lr /= nr;

  U1Frozen<D> s(cube<D>(4), 404); s.beta = beta; s.kappa = kappa; s.q = q; s.matter_metro = metro; s.hot(0.8);
  MucaB W(-400.0, 2100.0, 250);
  for (int i = 0; i < W.nbin; ++i) W.g[i] = c * (W.Bmin + (i + 0.5) * W.dB);   // g(B) = c*B exactly
  s.mucab = &W; s.muca_build = false;
  s.thermalize(400);
  double ls = 0; int ns = 0; for (int i = 0; i < 1500; ++i) { s.sweep(); ls += s.link_energy(); ++ns; } ls /= ns;

  char msg[120];
  std::snprintf(msg, sizeof msg, "muca g=cB %s <cos> %.4f vs unbiased kappa+c %.4f (q=%d)",
                metro ? "metro" : "hb", ls, lr, q);
  CHECK_CLOSE(ls, lr, 0.02, msg);
}

int main() {
  std::printf("-- von Mises heatbath --\n"); test_vonmises();
  std::printf("-- Z_q link-jump dB=0 --\n"); test_zq_hop_dB0<4>(2); test_zq_hop_dB0<4>(4); test_zq_hop_dB0<3>(6);
  std::printf("-- matter resultant consistency --\n"); test_matter_local<4>(2); test_matter_local<3>(3);
  std::printf("-- gauge local dS consistency --\n"); test_gauge_local<4>(2); test_gauge_local<3>(3);
  std::printf("-- gauge staple + overrelaxation --\n"); test_gauge_staple<4>(2); test_gauge_staple<3>(3);
  std::printf("-- kappa=0 vs HMC --\n"); test_kappa0_vs_hmc<4>();
  std::printf("-- multicanonical-in-accept (linear g) --\n"); test_muca_linear<4>(2, false); test_muca_linear<3>(3, false); test_muca_linear<4>(2, true);
  return report("test_u1_frozen");
}
