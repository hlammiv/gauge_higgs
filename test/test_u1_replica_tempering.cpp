// Tests for U1ReplicaTempering<D>: the U(1)+charge-q analog of the SU(2)
// ReplicaTempering. Covers (1) the CANONICAL GATE -- the swap Delta computed
// from the derived conjugate equals the brute-force action difference
// [S_i(C_j)+S_j(C_i)] - [S_i(C_i)+S_j(C_j)] for BOTH axes, proving the
// conjugate normalization/sign; (2) a SMOKE test on a kappa-ladder (q=4, L=4,
// enabled) where swaps actually accept with per-pair acceptance in (0,1); and
// (3) the TOGGLE -- enabled=false => 0 swap attempts.
#include "check.hpp"
#include "u1/u1_replica_tempering.hpp"
#include <array>
#include <vector>

using namespace gh;
using namespace gh::u1;

template <int D>
static std::array<int, D> cube(int n) { std::array<int, D> L{}; for (int mu = 0; mu < D; ++mu) L[mu] = n; return L; }

// Total running action of a U1HMC replica (gauge + scalar), the quantity whose
// difference under a config swap defines the canonical Delta.
template <int D>
static Real total_action(const U1HMC<D>& r) {
  return gauge_action<D>(r.th, r.lat, r.beta)
       + scalar_action<D>(r.phi, r.th, r.lat, r.q, r.kappa, r.lambda);
}

// Brute-force swap Delta: actually exchange (th,phi), recompute the two actions,
// then swap back. Delta = S_after - S_before, with each replica keeping its OWN
// coupling. This is the ground-truth the conjugate formula must reproduce.
template <int D>
static Real bruteforce_delta(U1HMC<D>& a, U1HMC<D>& b) {
  const Real before = total_action<D>(a) + total_action<D>(b);
  std::swap(a.th, b.th); std::swap(a.phi, b.phi);
  const Real after = total_action<D>(a) + total_action<D>(b);
  std::swap(a.th, b.th); std::swap(a.phi, b.phi);  // restore
  return after - before;
}

// ---- (1) canonical gate: conjugate Delta == brute-force Delta, both axes -----
static void test_canonical_gate() {
  constexpr int D = 4;
  const auto ext = cube<D>(4);
  const int  q = 4;

  // ---- BETA axis ----
  {
    std::vector<Real> beta_ladder = {1.0, 1.5};
    U1ReplicaTempering<D> pt(ext, TemperAxis::Beta, beta_ladder, /*seed0=*/11);
    pt.set_kappa(0.8); pt.set_lambda(0.5); pt.set_q(q); pt.set_tau(0.5); pt.set_nmd(8);
    // Distinct, nontrivial configs in each replica.
    pt.replica(0).hot(0.7); pt.replica(1).hot(0.7);
    pt.replica(0).cold_phi(1.0);  // give them genuinely different fields

    const Real Ai = U1ReplicaTempering<D>::beta_conjugate(pt.replica(0));
    const Real Aj = U1ReplicaTempering<D>::beta_conjugate(pt.replica(1));
    const Real conj_delta = (pt.coupling(0) - pt.coupling(1)) * (Aj - Ai);
    const Real bf = bruteforce_delta<D>(pt.replica(0), pt.replica(1));
    CHECK_CLOSE(conj_delta, bf, 1e-9, "BETA-axis swap Delta == brute-force action diff");
  }

  // ---- KAPPA axis (the metastability cure) ----
  {
    std::vector<Real> kappa_ladder = {0.5, 1.2};
    U1ReplicaTempering<D> pt(ext, TemperAxis::Kappa, kappa_ladder, /*seed0=*/23);
    pt.set_beta(1.1); pt.set_lambda(0.5); pt.set_q(q); pt.set_tau(0.5); pt.set_nmd(8);
    pt.replica(0).hot(0.7); pt.replica(1).hot(0.9);

    const Real Bi = U1ReplicaTempering<D>::kappa_conjugate(pt.replica(0));
    const Real Bj = U1ReplicaTempering<D>::kappa_conjugate(pt.replica(1));
    const Real conj_delta = (pt.coupling(0) - pt.coupling(1)) * (Bi - Bj);
    const Real bf = bruteforce_delta<D>(pt.replica(0), pt.replica(1));
    CHECK_CLOSE(conj_delta, bf, 1e-9, "KAPPA-axis swap Delta == brute-force action diff");
  }
}

// ---- (2) smoke test: 3-replica kappa-ladder, q=4, L=4, swaps accept ----------
static void test_smoke_kappa_ladder() {
  constexpr int D = 4;
  const auto ext = cube<D>(4);

  // Kappa ladder chosen close enough that adjacent swaps actually accept. The
  // conjugate B is the hopping sum over ALL links (B ~ O(vol*D), here ~thousands
  // at q=4), so the swap Delta ~ (Delta kappa)*(B_i-B_j) is large unless the
  // rungs are FINELY spaced for the hopping-energy distributions to overlap; a
  // dk~0.004 spacing yields a healthy ~0.6 swap acceptance here.
  std::vector<Real> kappa_ladder = {0.400, 0.404, 0.408};
  U1ReplicaTempering<D> pt(ext, TemperAxis::Kappa, kappa_ladder, /*seed0=*/12345);
  pt.set_beta(1.0); pt.set_lambda(0.5); pt.set_q(4); pt.set_tau(0.6); pt.set_nmd(12);
  pt.n_sweep = 2;
  pt.enabled = true;

  CHECK(pt.n_replicas() == 3, "smoke: 3 replicas built");
  CHECK(pt.n_pairs() == 2, "smoke: 2 adjacent pairs");

  // Hot starts on all rungs; the ladder spacing (not the start) sets the swap
  // overlap, so per-pair acceptance is reproducible rather than basin-dependent.
  pt.replica(0).hot(1.0);
  pt.replica(1).hot(1.0);
  pt.replica(2).hot(1.0);

  for (int it = 0; it < 120; ++it) pt.step();

  bool ran = pt.replica(0).traj_count > 0;
  CHECK(ran, "smoke: replicas ran trajectories");

  // HMC acceptance should be healthy (this regime had 0.97-1.0).
  for (std::size_t k = 0; k < pt.n_replicas(); ++k) {
    const double a = pt.hmc_acceptance(k);
    CHECK(a > 0.3, "smoke: HMC acceptance healthy");
  }

  // Per-pair swap acceptance strictly in (0,1): swaps both happen and are
  // sometimes rejected, i.e. the Metropolis coin is doing real work.
  for (std::size_t p = 0; p < pt.n_pairs(); ++p) {
    const double a = pt.pair_acceptance(p);
    std::printf("  pair %zu swap acceptance = %.3f  (attempts=%llu)\n",
                p, a, (unsigned long long)pt.swap_attempts[p]);
    CHECK(pt.swap_attempts[p] > 0, "smoke: pair had swap attempts");
    CHECK(a > 0.0 && a < 1.0, "smoke: pair acceptance strictly in (0,1)");
  }

  // avg_plaq accessor returns a sane value (<cos> in (-1,1]).
  for (std::size_t k = 0; k < pt.n_replicas(); ++k) {
    const Real pl = pt.avg_plaq(k);
    CHECK(pl > -1.0001 && pl <= 1.0001, "smoke: avg_plaq in range");
  }
}

// ---- (3) toggle: enabled=false => zero swap attempts -------------------------
static void test_toggle_off() {
  constexpr int D = 4;
  const auto ext = cube<D>(4);
  std::vector<Real> kappa_ladder = {0.30, 0.45, 0.60};
  U1ReplicaTempering<D> pt(ext, TemperAxis::Kappa, kappa_ladder, /*seed0=*/777);
  pt.set_beta(1.0); pt.set_lambda(0.5); pt.set_q(4); pt.set_tau(0.6); pt.set_nmd(12);
  pt.enabled = false;

  pt.replica(0).hot(1.0); pt.replica(1).hot(1.0); pt.replica(2).hot(1.0);
  for (int it = 0; it < 20; ++it) pt.step();

  std::uint64_t total_attempts = 0;
  for (std::size_t p = 0; p < pt.n_pairs(); ++p) total_attempts += pt.swap_attempts[p];
  CHECK(total_attempts == 0, "toggle off: zero swap attempts");
  CHECK(pt.replica(0).traj_count > 0, "toggle off: replicas still ran HMC");
}

int main() {
  test_canonical_gate();
  test_smoke_kappa_ladder();
  test_toggle_off();
  return report("test_u1_replica_tempering");
}
