// Tests for the per-point nmd auto-tuner: it must RAISE nmd at a stiff (coarse)
// point until acceptance recovers, NOT blow up an already-fine point, respect the
// nmd_max cap, terminate, and be deterministic. Tiny L=4 HMC (no physics scan).
#include "check.hpp"
#include "u1/u1.hpp"
#include "u1/autotune.hpp"
#include <array>
#include <cstdio>

using namespace gh;
using namespace gh::u1;

template <int D>
static std::array<int, D> cube(int n) { std::array<int, D> L{}; for (int m = 0; m < D; ++m) L[m] = n; return L; }

template <int D>
static double measure_acc(U1HMC<D>& h, int n) {
  h.traj_count = 0; h.accept_count = 0;
  for (int i = 0; i < n; ++i) h.trajectory();
  return h.acceptance();
}

int main() {
  std::printf("-- stiff point: tuner raises nmd and recovers acceptance --\n");
  {
    // nmd=1, tau=2 => dt=2: catastrophically coarse, acceptance ~0 by construction.
    U1HMC<4> h(cube<4>(4), 12345);
    h.beta = 1.0; h.kappa = 0.6; h.q = 2; h.lambda = 0.5; h.tau = 2.0; h.nmd = 1;
    h.hot(0.5); h.cold_phi(0.5);
    for (int i = 0; i < 60; ++i) h.trajectory();
    const double acc_before = measure_acc(h, 40);
    const int nmd0 = h.nmd;
    TuneResult r = tune_nmd(h, 30, 0.65, 0.92, 4, 400, 10);
    const double acc_after = measure_acc(h, 80);
    char m[128];
    std::snprintf(m, sizeof m, "stiff: nmd %d->%d  acc %.3f->%.3f  (iters=%d band=%d)",
                  nmd0, h.nmd, acc_before, acc_after, r.iters, (int)r.in_band);
    std::printf("   %s\n", m);
    CHECK(acc_before < 0.65, "stiff: start really is low-acceptance");
    CHECK(h.nmd > nmd0, "stiff: tuner raised nmd");
    CHECK(acc_after >= 0.55 || h.nmd >= 400, "stiff: acceptance recovered (or capped)");
    CHECK(r.iters <= 10, "stiff: terminated within max_iter");
  }

  std::printf("-- easy point: tuner does not blow nmd up --\n");
  {
    U1HMC<4> h(cube<4>(4), 999);
    h.beta = 0.9; h.kappa = 0.1; h.q = 1; h.lambda = 0.5; h.tau = 1.0; h.nmd = 20;
    h.hot(0.5); h.cold_phi(0.5);
    for (int i = 0; i < 60; ++i) h.trajectory();
    TuneResult r = tune_nmd(h, 30, 0.65, 0.92, 4, 400, 10);
    const double acc_after = measure_acc(h, 80);
    char m[128];
    std::snprintf(m, sizeof m, "easy: nmd 20->%d  acc_after=%.3f  (iters=%d band=%d)", h.nmd, acc_after, r.iters, (int)r.in_band);
    std::printf("   %s\n", m);
    CHECK(h.nmd <= 40, "easy: nmd not blown up (<= 2x start)");
    CHECK(acc_after >= 0.50, "easy: acceptance stays sane after tune");
  }

  std::printf("-- cap respected + terminates --\n");
  {
    U1HMC<4> h(cube<4>(4), 7);
    h.beta = 1.0; h.kappa = 0.6; h.q = 2; h.lambda = 0.5; h.tau = 2.0; h.nmd = 1;
    h.hot(0.5); h.cold_phi(0.5);
    for (int i = 0; i < 40; ++i) h.trajectory();
    TuneResult r = tune_nmd(h, 25, 0.65, 0.92, 4, 12, 8);   // nmd_max = 12
    CHECK(h.nmd <= 12, "cap: nmd <= nmd_max=12");
    CHECK(r.iters <= 8, "cap: terminated within max_iter");
  }

  std::printf("-- determinism: same seed -> same tuned nmd --\n");
  {
    auto run = []() {
      U1HMC<4> h(cube<4>(4), 24680);
      h.beta = 1.0; h.kappa = 0.6; h.q = 2; h.lambda = 0.5; h.tau = 2.0; h.nmd = 1;
      h.hot(0.5); h.cold_phi(0.5);
      for (int i = 0; i < 60; ++i) h.trajectory();
      return tune_nmd(h, 30, 0.65, 0.92, 4, 400, 10).nmd;
    };
    const int a = run(), b = run();
    char m[64]; std::snprintf(m, sizeof m, "determinism: %d == %d", a, b);
    CHECK(a == b, m);
  }

  return report("test_autotune");
}
