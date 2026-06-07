// MINIMAL proof-of-concept: charge-1 Wilson-loop string tension sigma_1 as the U(1)
// analog of SU(2) sigma_fund -- the order parameter that EXPOSES the deep-Higgs Z_q
// confined phase that the DeGrand-Toussaint monopole density rho_M is BLIND to.
//
// PHYSICS (the whole point).  A condensing charge-q scalar breaks U(1) -> Z_q. The
// charge-q condensate Higgs-SCREENS the U(1) monopoles, so rho_M collapses to ~0 by
// kappa~0.5 for ALL q -- it cannot see the residual Z_q gauge theory, which for q>=2
// STILL confines at small beta (via Z_q flux) all the way to kappa=inf. The charge-1
// Wilson loop carries Z_q charge 1 != 0 mod q, so it is UNSCREENED for q>=2 -> AREA law
// (sigma_1 > 0) in the Z_q-confined phase, perimeter law (sigma_1 = 0) when deconfined.
// For q=1 the charge-1 source IS the Higgs charge (1 == 0 mod 1, Z_1 trivial) -> always
// screened -> sigma_1 = 0: the clean NEGATIVE CONTROL.
//
// This is a standalone POC (NOT the full u1_scan integration). It reuses the existing
// charge-m Wilson loop machinery (u1.hpp) and the rep-agnostic Creutz/jackknife stack
// (measure/creutz*.hpp) verbatim, and reports the decisive table.
//
//   ./build/u1_zqprobe [L ntherm nmeas lambda base_seed Rmax n_scalar]
//
// Compile-time D = NDIM (use NDIM=4). NCOL is irrelevant for U(1).
#include "u1/u1.hpp"
#include "u1/monopole.hpp"            // monopole_density<D> (DeGrand-Toussaint rho_M)
#include "u1/autotune.hpp"            // tune_nmd
#include "measure/observables.hpp"    // Stats
#include "measure/creutz.hpp"         // creutz_ratio, string_tension_plateau
#include "measure/creutz_jack.hpp"    // CreutzJackknife, chi_diag_jack, plateau_sigma_jack
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <array>
#include <cmath>

using namespace gh;

namespace {

// Full charge-1 AND charge-q Wilson-loop grids in ONE sweep over the lattice and all
// mu<nu planes (symmetrized, exactly as u1::wilson_loop). g1[R][T]=<cos(theta_loop)>,
// gq[R][T]=<cos(q*theta_loop)> for 1<=R,T<=Rmax (row/col 0 unused, mirrors W[R][T]).
//
// EFFICIENCY: for each base site/plane, the oriented angle sum theta_loop(R,T) is built
// incrementally. We precompute the partial sums of the mu-bottom edge, mu-top edge, and
// the two nu side edges as functions of length, then loop_angle(R,T) = bottom(R) +
// rightside(R,T) - top(R,T) - leftside(T), reusing partials so the inner cost is O(Rmax^2)
// table lookups, not O(Rmax^2 * perimeter) re-walks. (Rmax<=4, so this stays trivial.)
template <int D>
void wilson_grids(const std::vector<Real>& th, const Lattice<D>& lat, int q, int Rmax,
                  WGrid& g1, WGrid& gq) {
  for (auto& r : g1) std::fill(r.begin(), r.end(), 0.0);
  for (auto& r : gq) std::fill(r.begin(), r.end(), 0.0);
  std::int64_t count = 0;
  #pragma omp parallel
  {
    WGrid l1(Rmax + 1, std::vector<Real>(Rmax + 1, 0.0));
    WGrid lq(Rmax + 1, std::vector<Real>(Rmax + 1, 0.0));
    std::int64_t lcount = 0;
    #pragma omp for schedule(static) nowait
    for (std::int64_t x = 0; x < lat.vol; ++x)
      for (int mu = 0; mu < D; ++mu)
        for (int nu = mu + 1; nu < D; ++nu) {
          // Direct incremental rectangle walk (clear, and Rmax tiny). For each R we
          // extend the bottom edge; for each T we extend the up/side; we recompute the
          // closing path lazily but cheaply. Simpler than caching all four edges and
          // identical numerically.
          for (int R = 1; R <= Rmax; ++R)
            for (int T = 1; T <= Rmax; ++T) {
              Real loop = 0.0; std::int64_t s = x;
              for (int i = 0; i < R; ++i) { loop += th[s * D + mu]; s = lat.neighbor_fwd(s, mu); }
              for (int j = 0; j < T; ++j) { loop += th[s * D + nu]; s = lat.neighbor_fwd(s, nu); }
              for (int i = 0; i < R; ++i) { s = lat.neighbor_bwd(s, mu); loop -= th[s * D + mu]; }
              for (int j = 0; j < T; ++j) { s = lat.neighbor_bwd(s, nu); loop -= th[s * D + nu]; }
              l1[R][T] += std::cos(loop);
              lq[R][T] += std::cos(q * loop);
              ++lcount;
            }
        }
    #pragma omp critical
    {
      for (int R = 1; R <= Rmax; ++R)
        for (int T = 1; T <= Rmax; ++T) { g1[R][T] += l1[R][T]; gq[R][T] += lq[R][T]; }
      count += lcount;
    }
  }
  // count tallies every (R,T,plane,site); per-(R,T) it is count/(Rmax*Rmax) but we tally
  // a uniform total, so normalize by the per-(R,T) sample count = count / (Rmax*Rmax).
  const Real per = static_cast<Real>(count) / static_cast<Real>(Rmax * Rmax);
  if (per > 0.0)
    for (int R = 1; R <= Rmax; ++R)
      for (int T = 1; T <= Rmax; ++T) { g1[R][T] /= per; gq[R][T] /= per; }
}

// <|P_n|> over configs and chi_P = Vspatial * Var(|P_n|): a cheap (finite-V-biased)
// cross-check on the same charge-n source. P_n = (1/Vsp) sum_x exp(i n line_t(x)).
template <int D>
Real polyakov_abs(const std::vector<Real>& th, const Lattice<D>& lat, int n, int tdir = D - 1) {
  const int Lt = lat.L[tdir];
  Real re = 0.0, im = 0.0; std::int64_t count = 0;
  std::array<int, D> xc{};
  for (std::int64_t s = 0; s < lat.vol; ++s) {
    lat.coords(s, xc); if (xc[tdir] != 0) continue;
    Real line = 0.0; std::int64_t cur = s;
    for (int t = 0; t < Lt; ++t) { line += th[cur * D + tdir]; cur = lat.neighbor_fwd(cur, tdir); }
    re += std::cos(n * line); im += std::sin(n * line); ++count;
  }
  if (!count) return 0.0;
  re /= count; im /= count;
  return std::sqrt(re * re + im * im);
}

struct Pt { int q; Real beta; Real kappa; };

}  // namespace

int main(int argc, char** argv) {
  auto af = [&](int i, double d) { return i < argc ? std::atof(argv[i]) : d; };
  auto ai = [&](int i, long d)   { return i < argc ? std::atol(argv[i]) : d; };
  const int    L       = static_cast<int>(ai(1, 8));
  const int    ntherm  = static_cast<int>(ai(2, 200));
  const int    nmeas   = static_cast<int>(ai(3, 800));
  const Real   lambda  = af(4, 0.5);
  const std::uint64_t base_seed = static_cast<std::uint64_t>(ai(5, 1));
  int          Rmax    = static_cast<int>(ai(6, std::min(L / 2, 4)));
  const int    n_scalar = static_cast<int>(ai(7, 6));
  if (Rmax < 2) Rmax = 2;

  std::array<int, kDim> ext{}; for (int mu = 0; mu < kDim; ++mu) ext[mu] = L;

  // Decisive grid: deep Higgs (kappa = 1.0, 1.5) scanned DOWN in beta, plus the kappa=0
  // U(1) sanity column. q=1 control, q=2 cleanest Z_2, q=8 receding/wedge.
  const std::vector<Real> kappas = {0.0, 1.0, 1.5};
  const std::vector<Real> betas  = {0.1, 0.2, 0.4, 0.6, 0.8};
  const std::vector<int>  qs     = {1, 2, 8};

  std::printf("# U(1)+charge-q Higgs: charge-1 Z_q-confinement probe (POC).\n");
  std::printf("# D=%d  L=%d^%d  lambda=%.3g  ntherm=%d  nmeas=%d  Rmax=%d  n_scalar=%d\n",
              kDim, L, kDim, lambda, ntherm, nmeas, Rmax, n_scalar);
  std::printf("# sigma_1 = charge-1 string tension (PRIMARY, area-vs-perimeter); chi11(2,2)\n");
  std::printf("#   = charge-1 Creutz ratio at (2,2); sigma_q = charge-q control (screened ->0\n");
  std::printf("#   deep Higgs); rho_M = DeGrand-Toussaint monopole density (BLIND to Z_q);\n");
  std::printf("#   |P1| = <|charge-1 Polyakov|> (cheap finite-V-biased cross-check).\n");
  std::printf("# n_block aims for ~20 blocks; jackknife errors are correlated.\n#\n");
  std::printf("%3s %6s %6s | %5s | %10s %9s | %10s | %10s %9s | %6s | %5s\n",
              "q", "beta", "kappa", "acc",
              "chi11(2,2)", "+/-", "sigma_1", "sigma_q", "+/-", "rho_M", "|P1|");
  std::printf("%s\n", "-----------------------------------------------------------------------------------------------------");

  const int target_blocks = 20;
  const int n_block = std::max(1, nmeas / target_blocks);

  for (int q : qs) {
    for (Real kappa : kappas) {
      for (Real beta : betas) {
        const std::uint64_t seed = Rng::key(base_seed, q, (std::int64_t)std::llround(beta * 1000 + kappa * 7));
        u1::U1HMC<kDim> hmc(ext, seed);
        hmc.beta = beta; hmc.kappa = kappa; hmc.lambda = lambda; hmc.q = q;
        hmc.tau = 1.0; hmc.nmd = 8; hmc.n_scalar = n_scalar;
        hmc.hot(0.8); hmc.cold_phi(0.5);
        for (int t = 0; t < ntherm; ++t) hmc.trajectory();
        const u1::TuneResult tr = u1::tune_nmd<kDim>(hmc);  // extra therm + sets nmd
        (void)tr;
        hmc.traj_count = 0; hmc.accept_count = 0;

        Stats rhoM, absP1;
        CreutzJackknife jk1(Rmax), jkq(Rmax);   // charge-1 and charge-q grids
        WGrid b1(Rmax + 1, std::vector<Real>(Rmax + 1, 0.0));
        WGrid bq(Rmax + 1, std::vector<Real>(Rmax + 1, 0.0));
        int bcount = 0;
        auto flush = [&]() {
          if (!bcount) return;
          WGrid g1(Rmax + 1, std::vector<Real>(Rmax + 1, 0.0));
          WGrid gq(Rmax + 1, std::vector<Real>(Rmax + 1, 0.0));
          for (int R = 1; R <= Rmax; ++R)
            for (int T = 1; T <= Rmax; ++T) {
              g1[R][T] = b1[R][T] / bcount; gq[R][T] = bq[R][T] / bcount;
            }
          jk1.add_block(g1); jkq.add_block(gq);
          for (auto& r : b1) std::fill(r.begin(), r.end(), 0.0);
          for (auto& r : bq) std::fill(r.begin(), r.end(), 0.0);
          bcount = 0;
        };

        WGrid w1(Rmax + 1, std::vector<Real>(Rmax + 1, 0.0));
        WGrid wq(Rmax + 1, std::vector<Real>(Rmax + 1, 0.0));
        for (int t = 0; t < nmeas; ++t) {
          hmc.trajectory();
          rhoM.add(u1::monopole_density<kDim>(hmc.th, hmc.lat));
          absP1.add(polyakov_abs<kDim>(hmc.th, hmc.lat, 1));
          wilson_grids<kDim>(hmc.th, hmc.lat, q, Rmax, w1, wq);
          for (int R = 1; R <= Rmax; ++R)
            for (int T = 1; T <= Rmax; ++T) { b1[R][T] += w1[R][T]; bq[R][T] += wq[R][T]; }
          if (++bcount == n_block) flush();
        }
        flush();

        // PRIMARY: chi(2,2) charge-1 Creutz ratio (orientation-symmetric, perimeter-clean).
        const JackResult c11 = chi_diag_jack(jk1, 2);
        // sigma_1 plateau over the NOISE-RELIABLE diagonal R only: a chi(R,R) is kept only
        // if its 2x2 constituent loops all stand >= n_sigma above their own jackknife error
        // (creutz_jack noise guard). On L=8, R=2 is solid; larger R is often pure noise at
        // small loops, so we let the guard decide rather than averaging garbage.
        const Real n_sigma = 2.0;
        const WGrid Wm1 = jk1.mean_grid(), We1 = jk1.err_grid();
        const WGrid Wmq = jkq.mean_grid(), Weq = jkq.err_grid();
        std::vector<int> R1, Rq;
        for (int R = 2; R <= Rmax; ++R) {
          if (creutz_excl_reason(Wm1, We1, R, R, n_sigma).empty()) R1.push_back(R);
          if (creutz_excl_reason(Wmq, Weq, R, R, n_sigma).empty()) Rq.push_back(R);
        }
        const JackResult s1 = plateau_sigma_jack(jk1, R1);
        const JackResult sq = plateau_sigma_jack(jkq, Rq);

        std::printf("%3d %6.2f %6.2f | %5.3f | %10.4f %9.4f | %10.4f | %10.4f %9.4f | %6.3f | %5.3f\n",
                    q, beta, kappa, hmc.acceptance(),
                    c11.value, c11.error,
                    s1.value,
                    sq.value, sq.error,
                    rhoM.mean(), absP1.mean());
        std::fflush(stdout);
      }
      std::printf("%s\n", "  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -");
    }
  }
  return 0;
}
