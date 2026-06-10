// =====================================================================================
// 2D LLR driver for the compact U(1) + charge-q Higgs model. Computes the JOINT density of
// states rho(A, B) (windowed pair (S1,S2)=(A, E2=-B)) via the coupled Robbins-Monro solve
// in src/u1/u1_llr.hpp, so a single run reweights to the whole (beta,kappa) phase plane.
//
// Modes (mirror .llr_ref/llr2d.cpp):
//   (default / "full") band-around-ridge tiling, cold+hot meet-in-the-gap passes.
//   "presample"        plain U1HMC at several (beta,kappa) -> (A_mean,B_mean) for the ridge fit.
//   "seed"             sequential drive-in: dump per-cell configs + a manifest (level-2 //ism).
//   "solve"            RM-solve a [lo,hi) cell range from a manifest (parallel-safe workers).
//   "fdcheck"          CI smoke: helper FD self-check + one-cell RM convergence.
//
// Emits a "LLR2D ... E1=[Atop,Abot] step1=.. hw1=.. c0=.. c1=.. c2=.. step2=.. hw2=.. nperp=.."
// geometry header (byte-compatible token layout with the python grid_index parser) and
// "ANE2: A0 E2_0 a1 a2 sd1 sd2 nrep" lines, plus "RMDIAG:" sidecars.
//
// CONVENTION (matches scan_obs.hpp): window E2=-B, kappa_eff=+a2, gauge weight exp(-a1 dA + a2 dB).
// Reconstruction reweight is drop-in: lnrho - beta*A - kappa*E2 = lnrho - (beta*A - kappa*B).
// =====================================================================================
#include "u1/u1.hpp"
#include "u1/u1_llr.hpp"
#include "u1/scan_obs.hpp"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <string>
#include <vector>
#include <array>
#include <cmath>

using namespace gh;

namespace {

std::array<int, kDim> make_ext(int L) {
  std::array<int, kDim> e{}; for (int mu = 0; mu < kDim; ++mu) e[mu] = L; return e;
}

#define RIDGE(c0, c1, c2, A) ((c0) + (c1) * (A) + (c2) * (A) * (A))

// Apply common engine tuning from CLI to an U1LLR<kDim>.
struct LLRParams {
  int    L = 4, q = 1;
  Real   lambda = 0.5;
  Real   Atop = 0, Abot = 0, step1 = 1, hw1 = 1;
  Real   c0 = 0, c1 = 0, c2 = 0, step2 = 1, hw2 = 1;
  int    nperp = 3;
  Real   a0 = 1.0;
  std::uint64_t seed = 1;
  // RM budget: convergence (a1->beta, a2->kappa, res<<hw, acc_h~0.9) needs K~100 NRM~200 per cell
  // (validated single-cell 2026-06-10); the old 40/70 under-converges -> a2 stuck near init.
  unsigned K = 100, NRM = 200, R = 1;
  // engine knobs. GENTLE gauge mixing (n_hit_g=5, n_over=1) keeps the guarded scalar HMC healthy
  // (acc_h~0.9); the old 10/2 over-drives B -> scalar HMC guard-rejects -> a2 runs away.
  Real   delta0 = 0.6; int n_hit_g = 5, n_over = 1;
  // HEAVY matter sampling (validated): the phi condensate needs this or <B> under-condenses and
  // a2 is biased (~27% low at weak 0.5/8/1). tau_s=1, n_md_s=20, n_site=8 -> a2->kappa, <B> ~3%.
  Real   tau_s = 1.0;  int n_md_s = 20;
  int    n_site = 8;   Real site_w = 0.3;
  // Per-cell RM init ridges so a1~beta(A) and a2~kappa(A) start near the answer on BOTH axes
  // (the single shared a0 over-biases the kappa axis -> widening cascade). Set from the
  // presample fit via env U1_LLR_BRIDGE="d0 d1 d2" (beta vs A) and U1_LLR_KRIDGE="e0 e1 e2"
  // (kappa vs A); unset => fall back to a0 (old behavior).
  bool   bridge_set = false, kridge_set = false;
  Real   d0 = 0, d1 = 0, d2 = 0, e0 = 0, e1 = 0, e2 = 0;
};

// Per-cell RM seeds. a1_init ~ beta(A), a2_init ~ kappa(A); both default to a0 if no ridge given.
inline Real a1_init(const LLRParams& P, Real A) { return P.bridge_set ? (P.d0 + P.d1 * A + P.d2 * A * A) : P.a0; }
inline Real a2_init(const LLRParams& P, Real A) { return P.kridge_set ? (P.e0 + P.e1 * A + P.e2 * A * A) : P.a0; }

void parse_init_ridges(LLRParams& P) {
  if (const char* s = std::getenv("U1_LLR_BRIDGE")) {
    double a, b, c; if (std::sscanf(s, "%lf %lf %lf", &a, &b, &c) == 3) { P.d0 = a; P.d1 = b; P.d2 = c; P.bridge_set = true; }
  }
  if (const char* s = std::getenv("U1_LLR_KRIDGE")) {
    double a, b, c; if (std::sscanf(s, "%lf %lf %lf", &a, &b, &c) == 3) { P.e0 = a; P.e1 = b; P.e2 = c; P.kridge_set = true; }
  }
}

void apply_knobs(u1::U1LLR<kDim>& E, const LLRParams& P) {
  E.q = P.q; E.lambda = P.lambda;
  E.delta0 = P.delta0; E.n_hit_g = P.n_hit_g; E.n_over = P.n_over;
  E.tau_s = P.tau_s; E.n_md_s = P.n_md_s;
  E.n_site = P.n_site; E.site_w = P.site_w;
}

void emit_header(const LLRParams& P, int N1, int N2) {
  std::printf("LLR2D(grp,D,Nt,Nx): U1q%d %d %d %d | E1=[%g,%g] step1=%g hw1=%g | "
              "ridge c0=%g c1=%g c2=%g step2=%g hw2=%g nperp=%d | "
              "q=%d lambda=%g a0=%g seed=%llu K=%u NRM=%u R=%u | %dx%d cells\n",
              P.q, kDim, P.L, P.L, P.Atop, P.Abot, P.step1, P.hw1,
              P.c0, P.c1, P.c2, P.step2, P.hw2, P.nperp,
              P.q, P.lambda, P.a0, (unsigned long long)P.seed, P.K, P.NRM, P.R, N1, N2);
}

// Adaptive RM with hw2 widening on guarded-HMC accept collapse (contract: up to 3 widenings).
void rm_adaptive(u1::U1LLR<kDim>& E, Real A0, Real E2_0, Real hw1, Real hw2_in,
                 Real a1g, Real a2g, unsigned K, unsigned NRM,
                 Real& a1, Real& a2, Real& hw2_used, u1::U1LLR<kDim>::RMDiag& diag) {
  Real hw2 = hw2_in; Real tau0 = E.tau_s;
  int nwiden = 0;
  for (;;) {
    E.rm_solve2d(A0, E2_0, hw1, hw2, a1g, a2g, K, NRM, a1, a2, &diag);
    if (diag.acc_hmc >= 0.3 || nwiden >= 3) break;
    hw2 *= 1.5; E.tau_s *= 0.7; ++nwiden;
    E.seed_cell(A0, E2_0, hw1, hw2, a1, a2, 4000);
  }
  diag.nwiden = nwiden;
  hw2_used = hw2;
  E.tau_s = tau0;
}

// =====================================================================================
// Ridge tiling: cold (frozen->disordered) + hot (disordered->frozen) meet-in-the-gap.
// In (A,E2=-B): confined/disordered = LARGE A (=> i=0, Atop), Higgs/ordered = SMALL A
// (=> i=N1-1, Abot). Cold start = ordered (th=0 small A) walks frozen i=N1-1 -> i=0; hot
// start = random walks disordered i=0 -> frozen i=N1-1.
// =====================================================================================
struct Cell { Real A0, E2_0; std::vector<Real> a1s, a2s; std::vector<Real> hw2s; bool reached = false; };

int mode_full(const LLRParams& P) {
  const unsigned CAP = 4000;
  const int N1 = (int)((P.Atop - P.Abot) / P.step1 + 0.5) + 1;
  const int N2 = 2 * P.nperp + 1;
  emit_header(P, N1, N2);

  const std::int64_t vol = [&]{ std::int64_t v = 1; for (int mu = 0; mu < kDim; ++mu) v *= P.L; return v; }();
  std::printf("NPLAQ: %lld\nNLINKS: %lld\n",
              (long long)(vol * kDim * (kDim - 1) / 2), (long long)(vol * kDim));
  std::fflush(stdout);   // header out before the long seeding loop (else `> file` looks like a hang)

  std::vector<std::vector<Cell>> grid(N1, std::vector<Cell>(N2));
  for (int i = 0; i < N1; ++i)
    for (int j = 0; j < N2; ++j) {
      const Real A0 = P.Atop - i * P.step1;
      grid[i][j].A0 = A0;
      grid[i][j].E2_0 = RIDGE(P.c0, P.c1, P.c2, A0) + (j - P.nperp) * P.step2;
    }

  for (unsigned r = 0; r < P.R; ++r) {
    u1::U1LLR<kDim> E(make_ext(P.L), P.seed + r * 1000003u);
    apply_knobs(E, P);

    auto walk = [&](bool cold, int& reach) {
      if (cold) E.cold(); else E.hot();
      reach = 0;
      for (int s = 0; s < N1; ++s) {
        const int i = cold ? (N1 - 1 - s) : s;          // cold: frozen->disordered, hot: dis->frozen
        Cell& cc = grid[i][P.nperp];
        const Real g1 = a1_init(P, cc.A0), g2 = a2_init(P, cc.A0);   // per-cell beta(A), kappa(A) seeds
        if (!E.seed_cell(cc.A0, cc.E2_0, P.hw1, P.hw2, g1, g2, CAP)) continue;
        Real b1, b2, hw2u; u1::U1LLR<kDim>::RMDiag d;
        rm_adaptive(E, cc.A0, cc.E2_0, P.hw1, P.hw2, g1, g2, P.K, P.NRM, b1, b2, hw2u, d);
        cc.a1s.push_back(b1); cc.a2s.push_back(b2); cc.hw2s.push_back(hw2u); cc.reached = true; ++reach;
        std::printf("RMDIAG: %9.2f %9.2f %9.5f %9.5f %9.5f %9.5f %6.3f %6.3f %d\n",
                    cc.A0, cc.E2_0, d.res_A, d.res_E2, d.slope_a1, d.slope_a2, d.acc_gauge, d.acc_hmc, d.nwiden);
        std::fflush(stdout);
        auto ridgecfg = E.snapshot();
        for (int dir = -1; dir <= 1; dir += 2)
          for (int t = 1; t <= P.nperp; ++t) {
            const int j = P.nperp + dir * t;
            Cell& c = grid[i][j];
            E.restore(ridgecfg); E.resync_E();
            if (!E.seed_cell(c.A0, c.E2_0, P.hw1, P.hw2, b1, b2, CAP)) break;
            Real q1, q2, hw2v; u1::U1LLR<kDim>::RMDiag dj;
            rm_adaptive(E, c.A0, c.E2_0, P.hw1, P.hw2, b1, b2, P.K, P.NRM, q1, q2, hw2v, dj);
            c.a1s.push_back(q1); c.a2s.push_back(q2); c.hw2s.push_back(hw2v); c.reached = true;
            std::printf("RMDIAG: %9.2f %9.2f %9.5f %9.5f %9.5f %9.5f %6.3f %6.3f %d\n",
                        c.A0, c.E2_0, dj.res_A, dj.res_E2, dj.slope_a1, dj.slope_a2, dj.acc_gauge, dj.acc_hmc, dj.nwiden);
            std::fflush(stdout);
          }
        E.restore(ridgecfg); E.resync_E();
      }
    };

    int creach = 0, hreach = 0;
    walk(true, creach);
    walk(false, hreach);
    std::printf("repeat %u: cold ridge %d/%d, hot ridge %d/%d  (delta0=%.3f site_w=%.3f g_driftA=%.2e g_driftB=%.2e)\n",
                r, creach, N1, hreach, N1, E.delta0, E.site_w, E.g_driftA, E.g_driftB);
  }

  std::printf("\n# E1 E2 a1 a2 sd1 sd2 nrep\n");
  int nreached = 0;
  for (int i = 0; i < N1; ++i)
    for (int j = 0; j < N2; ++j) {
      Cell& c = grid[i][j];
      if (!c.reached) continue;
      ++nreached;
      const int n = (int)c.a1s.size();
      Real m1 = 0, m2 = 0; for (int k = 0; k < n; ++k) { m1 += c.a1s[k]; m2 += c.a2s[k]; } m1 /= n; m2 /= n;
      Real v1 = 0, v2 = 0; for (int k = 0; k < n; ++k) { v1 += (c.a1s[k]-m1)*(c.a1s[k]-m1); v2 += (c.a2s[k]-m2)*(c.a2s[k]-m2); }
      const Real sd1 = n > 1 ? std::sqrt(v1/(n-1)) : 0, sd2 = n > 1 ? std::sqrt(v2/(n-1)) : 0;
      std::printf("ANE2: %9.2f %9.2f %10.6f %10.6f %9.5f %9.5f %d\n", c.A0, c.E2_0, m1, m2, sd1, sd2, n);
    }
  std::printf("reached %d/%d cells\n", nreached, N1 * N2);
  return 0;
}

// =====================================================================================
// presample: plain U1HMC at npts (beta,kappa) along the expected ribbon -> (A_mean,B_mean).
// Emits raw points; the python ridge fit produces (c0,c1,c2,Atop,Abot).
// =====================================================================================
int mode_presample(int argc, char** argv) {
  // usage: presample L q lambda npts beta0 beta1 kappa0 kappa1 [ntherm nmeas nmd tau seed]
  if (argc < 9) {
    std::fprintf(stderr,
      "usage: %s presample L q lambda npts beta0 beta1 kappa0 kappa1 [ntherm nmeas nmd tau seed]\n",
      "u1_llr");
    return 1;
  }
  auto af = [&](int i, double d){ return i < argc ? std::atof(argv[i]) : d; };
  auto ai = [&](int i, long d){ return i < argc ? std::atol(argv[i]) : d; };
  const int L = (int)ai(1, 4), q = (int)ai(2, 1);
  const Real lambda = af(3, 0.5);
  const int npts = (int)ai(4, 8);
  const Real b0 = af(5, 0.8), b1 = af(6, 1.6), k0 = af(7, 0.5), k1 = af(8, 0.1);
  const int ntherm = (int)ai(9, 100), nmeas = (int)ai(10, 300), nmd = (int)ai(11, 20);
  const Real tau = af(12, 1.0);
  const std::uint64_t seed = (std::uint64_t)ai(13, 1);

  std::printf("# presample U(1)+q=%d L=%d^%d lambda=%g : (beta,kappa)->(A_mean,B_mean,E2=-B)\n",
              q, L, kDim, lambda);
  std::printf("# beta kappa A_mean B_mean E2_mean\n");
  for (int p = 0; p < npts; ++p) {
    const Real f = (npts > 1) ? Real(p) / (npts - 1) : 0.0;
    const Real beta = b0 + f * (b1 - b0);
    const Real kappa = k0 + f * (k1 - k0);
    u1::U1HMC<kDim> hmc(make_ext(L), seed + p);
    hmc.beta = beta; hmc.kappa = kappa; hmc.lambda = lambda; hmc.q = q; hmc.tau = tau; hmc.nmd = nmd;
    hmc.hot(0.8); hmc.cold_phi(0.5);
    for (int t = 0; t < ntherm; ++t) hmc.trajectory();
    Real sA = 0, sB = 0; int nm = 0;
    for (int t = 0; t < nmeas; ++t) {
      hmc.trajectory();
      sA += u1::plaq_energy_sum<kDim>(hmc.th, hmc.lat);
      sB += u1::hop_energy_sum<kDim>(hmc.phi, hmc.th, hmc.lat, q);
      ++nm;
    }
    const Real Am = sA / nm, Bm = sB / nm;
    std::printf("PRESAMPLE: %8.4f %8.4f %12.4f %12.4f %12.4f\n", beta, kappa, Am, Bm, -Bm);
  }
  return 0;
}

// =====================================================================================
// fdcheck: CI smoke (STEP 1 helper FD + STEP 2 one-cell RM convergence). D=2 L=4 q=1 in secs.
// =====================================================================================
int mode_fdcheck(int argc, char** argv) {
  auto ai = [&](int i, long d){ return i < argc ? std::atol(argv[i]) : d; };
  auto af = [&](int i, double d){ return i < argc ? std::atof(argv[i]) : d; };
  const int L = (int)ai(1, 4), q = (int)ai(2, 1);
  const Real lambda = af(3, 0.5);
  const std::uint64_t seed = (std::uint64_t)ai(4, 12345);

  u1::U1LLR<kDim> E(make_ext(L), seed);
  E.q = q; E.lambda = lambda;
  const Real fderr = E.fd_self_check(60, 60);
  std::printf("FDCHECK: D=%d L=%d q=%d  max|helper - recompute-delta| = %.3e  (pass < 1e-9: %s)\n",
              kDim, L, q, fderr, fderr < 1e-9 ? "YES" : "NO");

  // STEP 2: one-cell RM convergence at a target IN THE BULK of support (not on the A=0 / B=max
  // boundary, where <A>/<E2> cannot pin to the centre). Warm a config near (beta~1,kappa~0.3),
  // then centre the cell on a slightly-disordered (A,B) so the box has support on both sides.
  apply_knobs(E, LLRParams{});                            // default knobs
  E.q = q; E.lambda = lambda;
  E.hot();
  for (int w = 0; w < 30; ++w) { long h = 0, a = 0; E.update_constrained2d(1.0, 0.5, E.A_run, -E.B_run, 1e9, 1e9, false, h, a); }
  const Real A0 = E.A_run, E2_0 = -E.B_run;              // centre = the warmed (bulk) config
  const Real hw1 = std::max(2.0, 0.04 * std::fabs(A0) + 2.0);
  const Real hw2 = std::max(4.0, 0.04 * std::fabs(E2_0) + 4.0);
  if (!E.seed_cell(A0, E2_0, hw1, hw2, 1.0, 1.0, 4000)) {
    std::printf("FDCHECK: seed_cell FAILED (cell off-ribbon)\n");
    return (fderr < 1e-9) ? 0 : 1;
  }
  Real a1, a2; u1::U1LLR<kDim>::RMDiag d;
  E.rm_solve2d(A0, E2_0, hw1, hw2, 1.0, 1.0, 30, 50, a1, a2, &d);
  std::printf("FDCHECK: cell (A0=%.2f E2_0=%.2f hw1=%.2f hw2=%.2f) -> a1=%.4f a2=%.4f\n",
              A0, E2_0, hw1, hw2, a1, a2);
  std::printf("FDCHECK:  res_A=%.3e res_E2=%.3e slope_a1=%.3e slope_a2=%.3e acc_g=%.3f acc_hmc=%.3f\n",
              d.res_A, d.res_E2, d.slope_a1, d.slope_a2, d.acc_gauge, d.acc_hmc);
  std::printf("FDCHECK:  g_driftA=%.3e g_driftB=%.3e (vs hw1=%.2f hw2=%.2f)\n",
              E.g_driftA, E.g_driftB, hw1, hw2);
  const bool conv = (d.res_A < 0.25 * (2 * hw1)) && (d.res_E2 < 0.25 * (2 * hw2));
  std::printf("FDCHECK: RM converged in cell: %s\n", conv ? "YES" : "NO");
  return (fderr < 1e-9) ? 0 : 1;
}

// =====================================================================================
// seed / solve (level-2 parallelism). Configs serialized as raw th[] + phi[].
// =====================================================================================
bool save_cfg(const char* path, const u1::U1LLR<kDim>& E) {
  FILE* f = std::fopen(path, "wb");
  if (!f) return false;
  const std::size_t nl = E.th.size(), ns = E.phi.size();
  std::fwrite(E.th.data(), sizeof(Real), nl, f);
  std::fwrite(E.phi.data(), sizeof(Complex), ns, f);
  std::fclose(f);
  return true;
}
bool load_cfg(const char* path, u1::U1LLR<kDim>& E) {
  FILE* f = std::fopen(path, "rb");
  if (!f) return false;
  const std::size_t nl = E.th.size(), ns = E.phi.size();
  const std::size_t a = std::fread(E.th.data(), sizeof(Real), nl, f);
  const std::size_t b = std::fread(E.phi.data(), sizeof(Complex), ns, f);
  std::fclose(f);
  return a == nl && b == ns;
}

// Parse the geometry CLI shared by full/seed (returns false on usage error).
bool parse_geo(int argc, char** argv, int base, LLRParams& P) {
  // expects: L q lambda Atop Abot step1 hw1 c0 c1 c2 step2 hw2 nperp a0 seed [K NRM R]
  if (argc < base + 15) return false;
  auto af = [&](int i){ return std::atof(argv[i]); };
  auto ai = [&](int i){ return std::atol(argv[i]); };
  P.L = (int)ai(base + 0); P.q = (int)ai(base + 1); P.lambda = af(base + 2);
  P.Atop = af(base + 3); P.Abot = af(base + 4); P.step1 = af(base + 5); P.hw1 = af(base + 6);
  P.c0 = af(base + 7); P.c1 = af(base + 8); P.c2 = af(base + 9);
  P.step2 = af(base + 10); P.hw2 = af(base + 11); P.nperp = (int)ai(base + 12);
  P.a0 = af(base + 13); P.seed = (std::uint64_t)ai(base + 14);
  if (argc > base + 15) P.K = (unsigned)ai(base + 15);
  if (argc > base + 16) P.NRM = (unsigned)ai(base + 16);
  if (argc > base + 17) P.R = (unsigned)ai(base + 17);
  return true;
}

int mode_seed(int argc, char** argv) {
  // usage: seed L q lambda Atop Abot step1 hw1 c0 c1 c2 step2 hw2 nperp a0 seed cfgdir
  LLRParams P;
  if (!parse_geo(argc, argv, 1, P) || argc < 17) {
    std::fprintf(stderr, "usage: u1_llr seed L q lambda Atop Abot step1 hw1 c0 c1 c2 step2 hw2 nperp a0 seed cfgdir\n");
    return 1;
  }
  const char* cfgdir = argv[16];
  parse_init_ridges(P);   // per-cell beta(A)/kappa(A) drive-in seeds (env U1_LLR_BRIDGE/KRIDGE)
  const int N1 = (int)((P.Atop - P.Abot) / P.step1 + 0.5) + 1;
  const int N2 = 2 * P.nperp + 1;
  emit_header(P, N1, N2);

  std::vector<std::vector<Cell>> grid(N1, std::vector<Cell>(N2));
  std::vector<std::vector<bool>> saved(N1, std::vector<bool>(N2, false));
  for (int i = 0; i < N1; ++i)
    for (int j = 0; j < N2; ++j) {
      const Real A0 = P.Atop - i * P.step1;
      grid[i][j].A0 = A0; grid[i][j].E2_0 = RIDGE(P.c0, P.c1, P.c2, A0) + (j - P.nperp) * P.step2;
    }
  const unsigned CAP = 4000;
  u1::U1LLR<kDim> E(make_ext(P.L), P.seed);
  apply_knobs(E, P);
  int idx = 0;
  char path[1024];
  auto try_save = [&](int i, int j) {
    if (saved[i][j]) return;
    Cell& c = grid[i][j];
    if (!(E.A_run >= c.A0 - P.hw1 && E.A_run <= c.A0 + P.hw1 &&
          -E.B_run >= c.E2_0 - P.hw2 && -E.B_run <= c.E2_0 + P.hw2)) return;
    std::snprintf(path, sizeof(path), "%s/cell_%d_%d.cfg", cfgdir, i, j);
    if (!save_cfg(path, E)) return;
    saved[i][j] = true;
    std::printf("CELL: %d %d %d %.2f %.2f %.4f %.4f %s\n", idx++, i, j, c.A0, c.E2_0, P.hw1, P.hw2, path);
  };
  auto fan = [&](int i) {
    auto ridgecfg = E.snapshot();
    for (int dir = -1; dir <= 1; dir += 2)
      for (int t = 1; t <= P.nperp; ++t) {
        const int j = P.nperp + dir * t;
        Cell& c = grid[i][j];
        E.restore(ridgecfg); E.resync_E();
        if (!E.seed_cell(c.A0, c.E2_0, P.hw1, P.hw2, a1_init(P, c.A0), a2_init(P, c.A0), CAP)) break;
        try_save(i, j);
      }
    E.restore(ridgecfg); E.resync_E();
  };
  // cold pass (frozen -> disordered)
  E.cold();
  for (int s = 0; s < N1; ++s) {
    const int i = N1 - 1 - s;
    Cell& cc = grid[i][P.nperp];
    if (!E.seed_cell(cc.A0, cc.E2_0, P.hw1, P.hw2, a1_init(P, cc.A0), a2_init(P, cc.A0), CAP)) continue;
    try_save(i, P.nperp); fan(i);
  }
  // hot pass (disordered -> frozen)
  E.hot();
  for (int i = 0; i < N1; ++i) {
    Cell& cc = grid[i][P.nperp];
    if (!E.seed_cell(cc.A0, cc.E2_0, P.hw1, P.hw2, a1_init(P, cc.A0), a2_init(P, cc.A0), CAP)) continue;
    try_save(i, P.nperp); fan(i);
  }
  std::fprintf(stderr, "SEED: saved %d cells to %s\n", idx, cfgdir);
  return 0;
}

int mode_solve(int argc, char** argv) {
  // usage: solve L q lambda manifest lo hi a0 K NRM seed
  if (argc < 11) {
    std::fprintf(stderr, "usage: u1_llr solve L q lambda manifest lo hi a0 K NRM seed\n");
    return 1;
  }
  auto af = [&](int i){ return std::atof(argv[i]); };
  auto ai = [&](int i){ return std::atol(argv[i]); };
  LLRParams P;
  P.L = (int)ai(1); P.q = (int)ai(2); P.lambda = af(3);
  const char* manifest = argv[4];
  const int lo = (int)ai(5), hi = (int)ai(6);
  const Real a0 = af(7);
  P.a0 = a0;
  parse_init_ridges(P);   // per-cell beta(A)/kappa(A) RM seeds (env U1_LLR_BRIDGE/KRIDGE)
  const unsigned K = (unsigned)ai(8), NRM = (unsigned)ai(9);
  const std::uint64_t seed = (std::uint64_t)ai(10);

  FILE* mf = std::fopen(manifest, "r");
  if (!mf) { std::fprintf(stderr, "solve: cannot open manifest %s\n", manifest); return 1; }
  u1::U1LLR<kDim> E(make_ext(P.L), seed);
  E.q = P.q; E.lambda = P.lambda;
  char line[2048]; int done = 0;
  while (std::fgets(line, sizeof(line), mf)) {
    if (std::strncmp(line, "CELL:", 5) != 0) continue;
    int idx, i, j; double A0, E2_0, hw1, hw2; char path[1024];
    if (std::sscanf(line, "CELL: %d %d %d %lf %lf %lf %lf %1023s", &idx, &i, &j, &A0, &E2_0, &hw1, &hw2, path) != 8) continue;
    if (idx < lo || idx >= hi) continue;
    if (!load_cfg(path, E)) { std::fprintf(stderr, "solve: missing cfg %s\n", path); continue; }
    E.resync_E();
    Real b1, b2, hw2u; u1::U1LLR<kDim>::RMDiag d;
    rm_adaptive(E, A0, E2_0, hw1, hw2, a1_init(P, A0), a2_init(P, A0), K, NRM, b1, b2, hw2u, d);
    std::printf("ANE2: %9.2f %9.2f %10.6f %10.6f %9.5f %9.5f %d\n", A0, E2_0, b1, b2, 0.0, 0.0, 1);
    std::printf("RMDIAG: %9.2f %9.2f %9.5f %9.5f %9.5f %9.5f %6.3f %6.3f %d\n",
                A0, E2_0, d.res_A, d.res_E2, d.slope_a1, d.slope_a2, d.acc_gauge, d.acc_hmc, d.nwiden);
    std::fflush(stdout);
    ++done;
  }
  std::fclose(mf);
  std::fprintf(stderr, "SOLVE[%d,%d): %d cells\n", lo, hi, done);
  return 0;
}

// =====================================================================================
// slab: 1D-A LLR at FIXED kappa. Windows only the gauge action A; matter runs physically.
// Tiles A (cold+hot meet-in-the-gap), RM-solves a1(A), measures physical <B>(A).
// Emits "ANE1: A0 a1 Bmean sd1 sdB n" + "RMDIAG1: A0 res_A slope_a1 acc_g acc_h".
// =====================================================================================
int mode_slab(int argc, char** argv) {
  // usage: slab L q lambda kappa Atop Abot step1 hw1 a0 seed K NRM [delta0 n_hit_g n_over tau_s n_md_s n_site site_w]
  if (argc < 13) {
    std::fprintf(stderr, "usage: u1_llr slab L q lambda kappa Atop Abot step1 hw1 a0 seed K NRM "
                         "[delta0 n_hit_g n_over tau_s n_md_s n_site site_w]\n");
    return 1;
  }
  auto af = [&](int i, double d){ return i < argc ? std::atof(argv[i]) : d; };
  auto ai = [&](int i, long d){ return i < argc ? std::atol(argv[i]) : d; };
  LLRParams P;
  P.L = (int)ai(1, 4); P.q = (int)ai(2, 1); P.lambda = af(3, 0.5);
  const Real kappa = af(4, 0.1);
  P.Atop = af(5, 0); P.Abot = af(6, 0); P.step1 = af(7, 1); P.hw1 = af(8, 1);
  P.a0 = af(9, 1.0); P.seed = (std::uint64_t)ai(10, 1); P.K = (unsigned)ai(11, 100); P.NRM = (unsigned)ai(12, 200);
  P.delta0 = af(13, P.delta0); P.n_hit_g = (int)ai(14, P.n_hit_g); P.n_over = (int)ai(15, P.n_over);
  P.tau_s = af(16, P.tau_s); P.n_md_s = (int)ai(17, P.n_md_s); P.n_site = (int)ai(18, P.n_site); P.site_w = af(19, P.site_w);
  parse_init_ridges(P);                                   // a1~beta(A) seed via U1_LLR_BRIDGE

  const int N1 = (P.step1 > 0) ? (int)((P.Atop - P.Abot) / P.step1 + 0.5) + 1 : 1;
  std::printf("LLR1DA(grp,D,Nt,Nx): U1q%d %d %d %d | kappa=%g | A=[%g,%g] step1=%g hw1=%g | "
              "q=%d lambda=%g a0=%g seed=%llu K=%u NRM=%u | %d cells\n",
              P.q, kDim, P.L, P.L, kappa, P.Atop, P.Abot, P.step1, P.hw1,
              P.q, P.lambda, P.a0, (unsigned long long)P.seed, P.K, P.NRM, N1);
  const std::int64_t vol = [&]{ std::int64_t v = 1; for (int mu = 0; mu < kDim; ++mu) v *= P.L; return v; }();
  std::printf("NPLAQ: %lld\nNLINKS: %lld\n",
              (long long)(vol * kDim * (kDim - 1) / 2), (long long)(vol * kDim));
  std::fflush(stdout);

  std::vector<Real> A0s(N1);
  std::vector<std::vector<Real>> a1acc(N1), Bacc(N1);
  for (int i = 0; i < N1; ++i) A0s[i] = P.Atop - i * P.step1;

  const unsigned CAP = 4000;
  u1::U1LLR<kDim> E(make_ext(P.L), P.seed);
  apply_knobs(E, P);
  auto walk = [&](bool cold) {
    if (cold) E.cold(); else E.hot();
    for (int s = 0; s < N1; ++s) {
      const int i = cold ? (N1 - 1 - s) : s;
      const Real g1 = a1_init(P, A0s[i]);
      if (!E.seed_slab(A0s[i], P.hw1, g1, kappa, CAP)) continue;
      Real a1, Bm; u1::U1LLR<kDim>::RMDiag d;
      E.rm_solve1d(A0s[i], P.hw1, g1, kappa, P.K, P.NRM, a1, Bm, &d);
      a1acc[i].push_back(a1); Bacc[i].push_back(Bm);
      std::printf("RMDIAG1: %9.2f %9.5f %9.5f %6.3f %6.3f\n", A0s[i], d.res_A, d.slope_a1, d.acc_gauge, d.acc_hmc);
      std::fflush(stdout);
    }
  };
  walk(true); walk(false);

  std::printf("\n# A0 a1 Bmean sd1 sdB n   (kappa=%g)\n", kappa);
  int nreached = 0;
  for (int i = 0; i < N1; ++i) {
    const int n = (int)a1acc[i].size();
    if (!n) continue;
    ++nreached;
    Real m1 = 0, mB = 0; for (int k = 0; k < n; ++k) { m1 += a1acc[i][k]; mB += Bacc[i][k]; } m1 /= n; mB /= n;
    Real v1 = 0, vB = 0; for (int k = 0; k < n; ++k) { v1 += (a1acc[i][k]-m1)*(a1acc[i][k]-m1); vB += (Bacc[i][k]-mB)*(Bacc[i][k]-mB); }
    const Real sd1 = n > 1 ? std::sqrt(v1/(n-1)) : 0, sdB = n > 1 ? std::sqrt(vB/(n-1)) : 0;
    std::printf("ANE1: %9.2f %10.6f %12.4f %9.5f %9.4f %d\n", A0s[i], m1, mB, sd1, sdB, n);
  }
  std::printf("reached %d/%d cells\n", nreached, N1);
  return 0;
}

// =====================================================================================
// rect: RECTANGULAR (A,B) 2D tile -- for regions with a first-order B-jump (e.g. the
// q<=4 triple point) that the ridge-band cannot span. Rasters the grid boustrophedon
// with CARRY-FORWARD of both the config and (a1,a2), so each cell is seeded from its
// neighbour (short drive-in) -- this walks INTO the suppressed first-order gap gradually
// from the physical side. Emits ANE2/RMDIAG; reconstruct with u1_llr_rect_reconstruct.py.
// =====================================================================================
int mode_rect(int argc, char** argv) {
  // usage: rect L q lambda Atop Abot step1 hw1 Bmin Bmax step2 hw2 a0 seed K NRM [delta0 n_hit_g n_over tau_s n_md_s n_site site_w]
  if (argc < 16) {
    std::fprintf(stderr, "usage: u1_llr rect L q lambda Atop Abot step1 hw1 Bmin Bmax step2 hw2 a0 seed K NRM "
                         "[delta0 n_hit_g n_over tau_s n_md_s n_site site_w]\n");
    return 1;
  }
  auto af = [&](int i, double d){ return i < argc ? std::atof(argv[i]) : d; };
  auto ai = [&](int i, long d){ return i < argc ? std::atol(argv[i]) : d; };
  LLRParams P;
  P.L = (int)ai(1, 4); P.q = (int)ai(2, 1); P.lambda = af(3, 0.5);
  P.Atop = af(4, 0); P.Abot = af(5, 0); P.step1 = af(6, 1); P.hw1 = af(7, 1);
  const Real Bmin = af(8, 0), Bmax = af(9, 1); P.step2 = af(10, 1); P.hw2 = af(11, 1);
  P.a0 = af(12, 1.0); P.seed = (std::uint64_t)ai(13, 1); P.K = (unsigned)ai(14, 100); P.NRM = (unsigned)ai(15, 200);
  P.delta0 = af(16, P.delta0); P.n_hit_g = (int)ai(17, P.n_hit_g); P.n_over = (int)ai(18, P.n_over);
  P.tau_s = af(19, P.tau_s); P.n_md_s = (int)ai(20, P.n_md_s); P.n_site = (int)ai(21, P.n_site); P.site_w = af(22, P.site_w);
  parse_init_ridges(P);

  const int N1 = (int)((P.Atop - P.Abot) / P.step1 + 0.5) + 1;     // A columns (high A -> low A)
  const int N2 = (int)((Bmax - Bmin) / P.step2 + 0.5) + 1;          // B rows (low B -> high B)
  std::printf("LLR2DRECT(grp,D,Nt,Nx): U1q%d %d %d %d | A=[%g,%g] step1=%g hw1=%g | "
              "B=[%g,%g] step2=%g hw2=%g | q=%d lambda=%g K=%u NRM=%u | %dx%d cells\n",
              P.q, kDim, P.L, P.L, P.Atop, P.Abot, P.step1, P.hw1, Bmin, Bmax, P.step2, P.hw2,
              P.q, P.lambda, P.K, P.NRM, N1, N2);
  const std::int64_t vol = [&]{ std::int64_t v = 1; for (int mu = 0; mu < kDim; ++mu) v *= P.L; return v; }();
  std::printf("NPLAQ: %lld\nNLINKS: %lld\n",
              (long long)(vol * kDim * (kDim - 1) / 2), (long long)(vol * kDim));
  std::fflush(stdout);

  const unsigned CAP = 8000;
  std::vector<std::vector<int>> done(N1, std::vector<int>(N2, 0));   // reached flags (union of passes)
  // Two rasters meet in the first-order gap. cold() ~ (low A, HIGH B = Higgs corner); hot() ~
  // (high A, low B = confined corner). Each pass starts at its natural corner and carry-forwards
  // config + (a1,a2) cell-to-cell (boustrophedon), walking INTO the gap from its physical side.
  auto raster = [&](bool cold) -> int {
    u1::U1LLR<kDim> E(make_ext(P.L), P.seed + (cold ? 0u : 777u));
    apply_knobs(E, P);
    if (cold) E.cold(); else E.hot();
    Real a1 = P.a0, a2 = 0.0; bool have = false; int reach = 0;
    // cold: A low->high (ii: N1-1 -> 0), B high->low first column. hot: A high->low, B low->high.
    for (int ii = 0; ii < N1; ++ii) {
      const int i = cold ? (N1 - 1 - ii) : ii;
      const Real A0 = P.Atop - i * P.step1;
      for (int jj = 0; jj < N2; ++jj) {
        // boustrophedon in B; cold begins each block from the HIGH-B end, hot from LOW-B end
        int jB = (ii % 2 == 0) ? jj : (N2 - 1 - jj);                // local raster order
        const int j = cold ? (N2 - 1 - jB) : jB;                    // cold starts high B, hot low B
        if (done[i][j]) continue;                                   // other pass already solved it
        const Real B0 = Bmin + j * P.step2, E2_0 = -B0;
        const Real g1 = a1_init(P, A0);
        const Real s1 = have ? a1 : g1, s2 = have ? a2 : a2_init(P, A0);
        if (!E.seed_cell(A0, E2_0, P.hw1, P.hw2, s1, s2, CAP)) continue;  // short drive from carried cfg
        Real b1, b2; u1::U1LLR<kDim>::RMDiag d;
        E.rm_solve2d(A0, E2_0, P.hw1, P.hw2, s1, s2, P.K, P.NRM, b1, b2, &d);
        a1 = b1; a2 = b2; have = true;
        if (done[i][j]++ == 0) ++reach;                              // count first reach only
        std::printf("ANE2: %9.2f %9.2f %10.6f %10.6f %9.5f %9.5f 1\n", A0, E2_0, b1, b2, 0.0, 0.0);
        std::printf("RMDIAG: %9.2f %9.2f %9.5f %9.5f %9.5f %9.5f %6.3f %6.3f %d\n",
                    A0, E2_0, d.res_A, d.res_E2, d.slope_a1, d.slope_a2, d.acc_gauge, d.acc_hmc, d.nwiden);
        std::fflush(stdout);
      }
    }
    return reach;
  };
  const int rc = raster(true);
  const int rh = raster(false);
  int uni = 0; for (int i = 0; i < N1; ++i) for (int j = 0; j < N2; ++j) uni += (done[i][j] > 0);
  std::printf("reached %d/%d cells (cold pass %d, hot pass %d)\n", uni, N1 * N2, rc, rh);
  return 0;
}

// =====================================================================================
// rectseed: SEEDING half of a parallel rectangular tile. Same two-pass meet-in-the-gap
// raster as mode_rect, but instead of RM-solving each cell it DRIVES the config in and
// DUMPS it + a CELL manifest line. Reuse mode_solve (parallel over cells) for the RM, so
// the expensive K*NRM solve parallelizes -> K=100/NRM=200 becomes tractable.
// Emits the LLR2DRECT header + "CELL: idx i j A0 E2_0 hw1 hw2 path" lines; configs -> cfgdir.
// =====================================================================================
int mode_rect_seed(int argc, char** argv) {
  // usage: rectseed L q lambda Atop Abot step1 hw1 Bmin Bmax step2 hw2 a0 seed cfgdir [knobs...]
  if (argc < 15) {
    std::fprintf(stderr, "usage: u1_llr rectseed L q lambda Atop Abot step1 hw1 Bmin Bmax step2 hw2 a0 seed cfgdir "
                         "[delta0 n_hit_g n_over tau_s n_md_s n_site site_w]\n");
    return 1;
  }
  auto af = [&](int i, double d){ return i < argc ? std::atof(argv[i]) : d; };
  auto ai = [&](int i, long d){ return i < argc ? std::atol(argv[i]) : d; };
  LLRParams P;
  P.L = (int)ai(1, 4); P.q = (int)ai(2, 1); P.lambda = af(3, 0.5);
  P.Atop = af(4, 0); P.Abot = af(5, 0); P.step1 = af(6, 1); P.hw1 = af(7, 1);
  const Real Bmin = af(8, 0), Bmax = af(9, 1); P.step2 = af(10, 1); P.hw2 = af(11, 1);
  P.a0 = af(12, 1.0); P.seed = (std::uint64_t)ai(13, 1);
  const char* cfgdir = argv[14];
  P.delta0 = af(15, P.delta0); P.n_hit_g = (int)ai(16, P.n_hit_g); P.n_over = (int)ai(17, P.n_over);
  P.tau_s = af(18, P.tau_s); P.n_md_s = (int)ai(19, P.n_md_s); P.n_site = (int)ai(20, P.n_site); P.site_w = af(21, P.site_w);
  parse_init_ridges(P);

  const int N1 = (int)((P.Atop - P.Abot) / P.step1 + 0.5) + 1;
  const int N2 = (int)((Bmax - Bmin) / P.step2 + 0.5) + 1;
  std::printf("LLR2DRECT(grp,D,Nt,Nx): U1q%d %d %d %d | A=[%g,%g] step1=%g hw1=%g | "
              "B=[%g,%g] step2=%g hw2=%g | q=%d lambda=%g | %dx%d cells\n",
              P.q, kDim, P.L, P.L, P.Atop, P.Abot, P.step1, P.hw1, Bmin, Bmax, P.step2, P.hw2,
              P.q, P.lambda, N1, N2);
  const std::int64_t vol = [&]{ std::int64_t v = 1; for (int mu = 0; mu < kDim; ++mu) v *= P.L; return v; }();
  std::printf("NPLAQ: %lld\nNLINKS: %lld\n",
              (long long)(vol * kDim * (kDim - 1) / 2), (long long)(vol * kDim));
  std::fflush(stdout);

  const unsigned CAP = 8000;
  std::vector<std::vector<int>> done(N1, std::vector<int>(N2, 0));
  int idxc = 0;
  char path[1024];
  auto raster = [&](bool cold) {
    u1::U1LLR<kDim> E(make_ext(P.L), P.seed + (cold ? 0u : 777u));
    apply_knobs(E, P);
    if (cold) E.cold(); else E.hot();
    Real a1 = P.a0, a2 = 0.0; bool have = false;
    for (int ii = 0; ii < N1; ++ii) {
      const int i = cold ? (N1 - 1 - ii) : ii;
      const Real A0 = P.Atop - i * P.step1;
      for (int jj = 0; jj < N2; ++jj) {
        int jB = (ii % 2 == 0) ? jj : (N2 - 1 - jj);
        const int j = cold ? (N2 - 1 - jB) : jB;
        if (done[i][j]) continue;
        const Real B0 = Bmin + j * P.step2, E2_0 = -B0;
        const Real s1 = have ? a1 : a1_init(P, A0), s2 = have ? a2 : a2_init(P, A0);
        if (!E.seed_cell(A0, E2_0, P.hw1, P.hw2, s1, s2, CAP)) continue;
        a1 = s1; a2 = s2; have = true;                                // carry the drive config forward
        std::snprintf(path, sizeof(path), "%s/cell_%d_%d.cfg", cfgdir, i, j);
        if (!save_cfg(path, E)) continue;
        done[i][j] = 1;
        std::printf("CELL: %d %d %d %.2f %.2f %.4f %.4f %s\n", idxc++, i, j, A0, E2_0, P.hw1, P.hw2, path);
        std::fflush(stdout);
      }
    }
  };
  raster(true); raster(false);
  std::fprintf(stderr, "RECTSEED: %d cells dumped to %s\n", idxc, cfgdir);
  return 0;
}

void usage(const char* prog) {
  std::fprintf(stderr,
    "usage:\n"
    "  %s L q lambda Atop Abot step1 hw1 c0 c1 c2 step2 hw2 nperp a0 seed [K NRM R]\n"
    "        [delta0 n_hit_g n_over tau_s n_md_s n_site site_w]      (full ridge tiling)\n"
    "  %s presample L q lambda npts beta0 beta1 kappa0 kappa1 [ntherm nmeas nmd tau seed]\n"
    "  %s seed  L q lambda Atop Abot step1 hw1 c0 c1 c2 step2 hw2 nperp a0 seed cfgdir\n"
    "  %s solve L q lambda manifest lo hi a0 K NRM seed\n"
    "  %s fdcheck [L q lambda seed]                                  (CI smoke; build D=2 for speed)\n",
    prog, prog, prog, prog, prog);
}

}  // namespace

int main(int argc, char** argv) {
  if (argc > 1 && std::strcmp(argv[1], "presample") == 0) return mode_presample(argc - 1, argv + 1);
  if (argc > 1 && std::strcmp(argv[1], "seed") == 0)      return mode_seed(argc - 1, argv + 1);
  if (argc > 1 && std::strcmp(argv[1], "solve") == 0)     return mode_solve(argc - 1, argv + 1);
  if (argc > 1 && std::strcmp(argv[1], "fdcheck") == 0)   return mode_fdcheck(argc - 1, argv + 1);
  if (argc > 1 && std::strcmp(argv[1], "slab") == 0)      return mode_slab(argc - 1, argv + 1);
  if (argc > 1 && std::strcmp(argv[1], "rect") == 0)      return mode_rect(argc - 1, argv + 1);
  if (argc > 1 && std::strcmp(argv[1], "rectseed") == 0)  return mode_rect_seed(argc - 1, argv + 1);

  // default: full ridge tiling
  if (argc < 16) { usage(argv[0]); return 1; }
  LLRParams P;
  if (!parse_geo(argc, argv, 1, P)) { usage(argv[0]); return 1; }
  // optional engine knobs after [K NRM R] (positions 16,17,18 then 19..25)
  auto af = [&](int i, double d){ return i < argc ? std::atof(argv[i]) : d; };
  auto ai = [&](int i, long d){ return i < argc ? std::atol(argv[i]) : d; };
  P.delta0  = af(19, P.delta0);
  P.n_hit_g = (int)ai(20, P.n_hit_g);
  P.n_over  = (int)ai(21, P.n_over);
  P.tau_s   = af(22, P.tau_s);
  P.n_md_s  = (int)ai(23, P.n_md_s);
  P.n_site  = (int)ai(24, P.n_site);
  P.site_w  = af(25, P.site_w);
  parse_init_ridges(P);   // per-cell beta(A)/kappa(A) RM seeds (env U1_LLR_BRIDGE/KRIDGE)
  return mode_full(P);
}
