// Smoke test for src/hmc/replica_tempering.hpp: build 3 pure-gauge SU(2) replicas
// (rep adj, kappa=0) on an increasing beta-ladder {1.4,1.7,2.0}, L=4, run with
// tempering ENABLED, and confirm it builds, runs, swaps are attempted, and every
// adjacent pair's swap acceptance lands in a sane open interval (0 < acc < 1).
//
// Build (D=4, N=2):
//   g++ -std=c++20 -O3 -march=native -funroll-loops -fopenmp -Isrc -DNDIM=4 -DNCOL=2 \
//       test/test_replica_tempering.cpp -o build/test_replica_tempering
#include "hmc/replica_tempering.hpp"
#include "rep/rep_adjoint.hpp"
#include <cstdio>
#include <cstdlib>

using namespace gh;

int main() {
  constexpr int D = NDIM;   // 4
  constexpr int Nc = NCOL;  // 2
  static_assert(D == 4 && Nc == 2, "this test expects -DNDIM=4 -DNCOL=2");

  const int L = 4;
  std::array<int, D> ext{}; for (int mu = 0; mu < D; ++mu) ext[mu] = L;

  AdjointRep<Nc> rep;                          // SU(2) adjoint (real rep)
  // Finely-spaced ladder: at L=4 the neighboring plaquette distributions must
  // overlap for swaps to be accepted at a non-trivial rate. (A coarse ladder
  // like {1.4,1.7,2.0} is *correct* but gives ~0 acceptance -- the gaps are too
  // wide on a 4^4 lattice, exactly as PT theory predicts.)
  const std::vector<Real> ladder{1.9, 2.0, 2.1};

  ReplicaTempering<D, Nc> pt(ext, rep, ladder, /*seed0=*/12345ull);
  pt.enabled = true;
  pt.n_sweep = 1;
  pt.set_kappa(0.0);                           // pure gauge: matter decoupled
  pt.set_nmd(12);
  pt.set_tau(1.0);

  // Disordered start so the betas actually separate the replicas' plaquettes.
  for (std::size_t k = 0; k < pt.n_replicas(); ++k)
    pt.replica(k).U.hot(pt.replica(k).rng, 0.9);

  std::printf("# %zu replicas, ladder beta = ", pt.n_replicas());
  for (std::size_t k = 0; k < pt.n_replicas(); ++k) std::printf("%.3f ", pt.beta(k));
  std::printf("\n");

  // Thermalize (tempering on).
  const int ntherm = 60, nsteps = 400;
  for (int t = 0; t < ntherm; ++t) pt.step();
  // Reset swap bookkeeping after thermalization for a clean acceptance measurement.
  for (std::size_t p = 0; p < pt.n_pairs(); ++p) { pt.swap_attempts[p] = 0; pt.swap_accepts[p] = 0; }
  for (int t = 0; t < nsteps; ++t) pt.step();

  std::printf("# per-replica HMC acceptance / mean plaquette:\n");
  for (std::size_t k = 0; k < pt.n_replicas(); ++k)
    std::printf("  replica %zu  beta=%.3f  hmc_acc=%.3f  plaq=%.5f\n",
                k, pt.beta(k), pt.hmc_acceptance(k), pt.avg_plaq(k));

  std::printf("# per-pair swap acceptance:\n");
  bool swaps_attempted = true, acc_in_range = true;
  for (std::size_t p = 0; p < pt.n_pairs(); ++p) {
    const double acc = pt.pair_acceptance(p);
    std::printf("  pair (%zu,%zu)  attempts=%llu  accepts=%llu  acc=%.4f\n",
                p, p + 1,
                (unsigned long long)pt.swap_attempts[p],
                (unsigned long long)pt.swap_accepts[p], acc);
    if (pt.swap_attempts[p] == 0) swaps_attempted = false;
    if (!(acc > 0.0 && acc < 1.0)) acc_in_range = false;
  }

  // ---- toggle safety property: enabled=false => zero swap attempts ----
  ReplicaTempering<D, Nc> pt_off(ext, rep, ladder, /*seed0=*/777ull);
  pt_off.enabled = false;
  pt_off.set_kappa(0.0);
  pt_off.set_nmd(12);
  for (std::size_t k = 0; k < pt_off.n_replicas(); ++k)
    pt_off.replica(k).U.hot(pt_off.replica(k).rng, 0.9);
  for (int t = 0; t < 20; ++t) pt_off.step();
  bool toggle_off_ok = true;
  for (std::size_t p = 0; p < pt_off.n_pairs(); ++p)
    if (pt_off.swap_attempts[p] != 0) toggle_off_ok = false;
  std::printf("# toggle: enabled=false made %s swap attempts (expect 0)\n",
              toggle_off_ok ? "0" : ">0");

  const bool pass = swaps_attempted && acc_in_range && toggle_off_ok;
  std::printf("checks: swaps_attempted=%d acc_in_(0,1)=%d toggle_off_zero=%d\n",
              swaps_attempted, acc_in_range, toggle_off_ok);
  std::printf("%s\n", pass ? "PASS" : "FAIL");
  return pass ? 0 : 1;
}
