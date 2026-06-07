#pragma once
// Toggleable parallel tempering (replica exchange) for the compact U(1) +
// charge-q Higgs HMC -- the abelian analog of src/hmc/replica_tempering.hpp
// (the validated SU(2) ReplicaTempering<D,N>).
//
// MOTIVATION. At deep Higgs (kappa large, e.g. kappa=4, L=2) the U(1)+charge-q
// plaquette is STRONGLY METASTABLE for q>=2 (measured hot-vs-cold-start <plaq>
// gaps: q=4 -> 0.48, q=5 -> 0.60, q=8 -> 0.62) while HMC acceptance stays
// healthy (0.97-1.0): genuine metastability, the chain trapped in a
// start-dependent basin. q=1 is ergodic (gap 0). A single chain cannot sample
// the equilibrium distribution; replica exchange lets a stuck replica tunnel
// between basins by swapping with neighbors at adjacent couplings.
//
// LADDER. Replicas k=0..M-1 sit on a coupling LADDER along ONE axis (Beta or
// Kappa, selected at construction). The varying coupling increases strictly
// along the ladder; EVERY other coupling (kappa or beta, plus lambda, q, tau,
// nmd, n_scalar) is shared/fixed across all replicas. Each replica is its own
// U1HMC<D> with a distinct RNG seed.
//
// CYCLE (one call to step()):
//   1. every replica runs n_sweep independent HMC trajectories;
//   2. if enabled, attempt config swaps on adjacent pairs in two passes:
//        even pass: (0,1),(2,3),(4,5),...   then
//        odd  pass: (1,2),(3,4),(5,6),...
//      so on average every replica gets a swap proposal with each neighbor.
//
// =========================================================================
// SWAP CONJUGATES -- DERIVED FROM THE ACTUAL u1.hpp action functions.
// =========================================================================
// The total action evaluated by U1HMC::hamiltonian() is
//   S = gauge_action<D>(th,lat,beta) + scalar_action<D>(phi,th,lat,q,kappa,lambda)
// with, reading the code in src/u1/u1.hpp verbatim:
//
//   gauge_action  = beta * A,   A := sum_plaq (1 - cos theta_plaq)              (>=0)
//   scalar_action = onsite - kappa * B,
//                   B := sum_{x,mu} 2 Re[ conj(phi_x) e^{i q theta_mu(x)} phi_{x+mu} ]
//   onsite        = sum_x [ |phi|^2 + lambda(|phi|^2 - 1)^2 ]   (no beta/kappa dep.)
//
// So the conjugate energies (held BY THE CONFIG and exchanged with it) are:
//   BETA  axis:  conjugate = +dS/dbeta  = A   (the plaquette sum, sum(1-cos))
//   KAPPA axis:  conjugate = -dS/dkappa = B   (the gauge-invariant hopping sum)
// We compute A and B DIRECTLY from the action functions below (helpers
// beta_conjugate / kappa_conjugate) so the swap Delta carries the EXACT
// normalization of the running action -- no assumed prefactors.
//
// Metropolis weight for proposing to exchange the configs of adjacent replicas
// i (coupling g_i) and j (coupling g_j), accept w.p. min(1, exp(-Delta)):
//
//   Delta = [S_i(C_j) + S_j(C_i)] - [S_i(C_i) + S_j(C_j)].
//
//   BETA axis (only beta differs; onsite & -kappa*B identical => cancel):
//     S = beta*A + (const wrt swap)
//     Delta = (beta_i*A_j + beta_j*A_i) - (beta_i*A_i + beta_j*A_j)
//           = (beta_i - beta_j) * (A_j - A_i).
//
//   KAPPA axis (only kappa differs; beta*A & onsite identical => cancel; note
//   the SIGN: scalar_action contains MINUS kappa*B):
//     S = -kappa*B + (const wrt swap)
//     Delta = (-kappa_i*B_j - kappa_j*B_i) - (-kappa_i*B_i - kappa_j*B_j)
//           = (kappa_i - kappa_j) * (B_i - B_j).
//
// Why the other action pieces cancel: along the Beta axis only beta differs, so
// the scalar_action coefficient (kappa) and the onsite term are IDENTICAL in
// both replicas; the configs being exchanged carry their own scalar/onsite
// energy with them, so those contributions appear symmetrically in S_i(C_j) and
// S_j(C_i) and drop out of Delta. Symmetrically along the Kappa axis, beta*A and
// onsite are identical and cancel. Only the term whose coefficient DIFFERS
// between the two replicas survives in Delta. (Verified by the canonical-gate
// reduction above and by the smoke test.)
//
// Momenta p (link) and pi (scalar) are refreshed from a Gaussian heatbath at the
// start of every trajectory (U1HMC::refresh_momenta), so they need not be
// swapped -- the next trajectory redraws them regardless. On accept we exchange
// BOTH the link angles .th and the scalar field .phi.
//
// =========================================================================
// TOGGLE (safety property for merging).
// =========================================================================
// `enabled == false` => step() attempts NO swaps. The M replicas then evolve as
// M completely independent single-coupling HMC streams (each with its own seed),
// so their per-replica time series merge exactly as ordinary independent runs.
#include "u1/u1.hpp"
#include "core/geometry.hpp"
#include "core/rng.hpp"
#include "core/config.hpp"
#include <vector>
#include <array>
#include <cstdint>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <cctype>
#include <string>

namespace gh {
namespace u1 {

// Which coupling varies along the ladder. The OTHER coupling is shared/fixed.
enum class TemperAxis { Beta, Kappa };

template <int D>
struct U1ReplicaTempering {
  std::vector<U1HMC<D>> replicas;  // one per ladder rung
  std::vector<Real> ladder;        // varying coupling per rung (strictly increasing)
  TemperAxis axis = TemperAxis::Kappa;
  bool enabled = true;             // master swap toggle
  int  n_sweep = 1;                // HMC trajectories per replica per step()

  // Swap bookkeeping, indexed by adjacent pair p=(p,p+1), p=0..M-2.
  std::vector<std::uint64_t> swap_attempts;
  std::vector<std::uint64_t> swap_accepts;

  // Dedicated RNG stream for the swap Metropolis coin (independent of the HMC
  // accept/reject streams so tempering never correlates with trajectory MC).
  Rng           swap_rng;
  std::uint64_t swap_counter = 0;

  // Construct M replicas on the given coupling ladder along `axis`. All replicas
  // share lattice extents and (initially default) couplings; only the varying
  // coupling and the RNG seed differ. Use the shared setters below to fix the
  // OTHER couplings (and lambda/q/tau/nmd/n_scalar) before running. `seed0`
  // seeds replica 0; replica k is seeded seed0 + k*seed_stride (distinct
  // streams). The varying coupling is written into each replica here too.
  U1ReplicaTempering(const std::array<int, D>& extents,
                     TemperAxis axis_,
                     const std::vector<Real>& coupling_ladder,
                     std::uint64_t seed0,
                     std::uint64_t seed_stride = 1000003ull,
                     std::uint64_t swap_seed   = 0x5A7AB1ECull)
      : ladder(coupling_ladder), axis(axis_), swap_rng(swap_seed) {
    if (coupling_ladder.size() < 2)
      throw std::runtime_error("U1ReplicaTempering: need >=2 replicas");
    for (std::size_t k = 1; k < coupling_ladder.size(); ++k)
      if (!(coupling_ladder[k] > coupling_ladder[k - 1]))
        throw std::runtime_error("U1ReplicaTempering: coupling ladder must be strictly increasing");

    replicas.reserve(coupling_ladder.size());
    for (std::size_t k = 0; k < coupling_ladder.size(); ++k) {
      replicas.emplace_back(extents, seed0 + k * seed_stride);
      set_axis_coupling(replicas.back(), coupling_ladder[k]);
    }
    const std::size_t npair = coupling_ladder.size() - 1;
    swap_attempts.assign(npair, 0);
    swap_accepts.assign(npair, 0);
  }

  std::size_t n_replicas() const { return replicas.size(); }
  std::size_t n_pairs()    const { return replicas.empty() ? 0 : replicas.size() - 1; }

  // --- shared-coupling setters (apply identical couplings to every replica) ---
  // Setting the coupling that is the LADDER axis would clobber the ladder, so
  // those setters re-apply the per-rung ladder value afterwards to stay safe;
  // in normal use you only set the FIXED couplings here.
  void set_beta(Real beta) {
    for (auto& r : replicas) r.beta = beta;
    if (axis == TemperAxis::Beta) reapply_ladder();
  }
  void set_kappa(Real kappa) {
    for (auto& r : replicas) r.kappa = kappa;
    if (axis == TemperAxis::Kappa) reapply_ladder();
  }
  void set_lambda(Real lambda)   { for (auto& r : replicas) r.lambda = lambda; }
  void set_q(int q)              { for (auto& r : replicas) r.q = q; }
  void set_tau(Real tau)         { for (auto& r : replicas) r.tau = tau; }
  void set_nmd(int nmd)          { for (auto& r : replicas) r.nmd = nmd; }
  void set_n_scalar(int n)       { for (auto& r : replicas) r.n_scalar = n; }
  void set_integ(Integ integ)    { for (auto& r : replicas) r.integ = integ; }

  // ---- accessors ----------------------------------------------------------
  U1HMC<D>&       replica(std::size_t k)       { return replicas.at(k); }
  const U1HMC<D>& replica(std::size_t k) const { return replicas.at(k); }
  // Varying coupling (beta or kappa, per axis) of rung k.
  Real coupling(std::size_t k) const { return ladder.at(k); }

  // Average plaquette <cos theta_pl> in (.., 1] for replica k.
  Real avg_plaq(std::size_t k) const {
    const auto& r = replicas.at(k);
    return avg_plaquette<D>(r.th, r.lat);
  }

  // Empirical swap-acceptance rate for adjacent pair p=(p,p+1).
  double pair_acceptance(std::size_t p) const {
    return swap_attempts.at(p) ? double(swap_accepts.at(p)) / double(swap_attempts.at(p)) : 0.0;
  }
  // HMC trajectory acceptance for replica k.
  double hmc_acceptance(std::size_t k) const { return replicas.at(k).acceptance(); }

  // ---- the cycle ----------------------------------------------------------
  void step() {
    // 1. each replica advances n_sweep independent trajectories.
    for (auto& r : replicas)
      for (int s = 0; s < n_sweep; ++s) r.trajectory();

    // 2. swaps (only if enabled). Even pass then odd pass over adjacent pairs.
    if (!enabled) return;
    const std::size_t M = replicas.size();
    for (std::size_t start = 0; start < 2; ++start)            // 0 = even, 1 = odd
      for (std::size_t p = start; p + 1 < M; p += 2)
        attempt_swap(p);
  }

  // Propose exchanging the configs of adjacent replicas p and p+1.
  bool attempt_swap(std::size_t p) {
    U1HMC<D>& ri = replicas[p];
    U1HMC<D>& rj = replicas[p + 1];

    Real Delta;
    if (axis == TemperAxis::Beta) {
      // conjugate = A = sum_plaq(1-cos);  Delta = (beta_i-beta_j)(A_j-A_i).
      const Real Ai = beta_conjugate(ri);
      const Real Aj = beta_conjugate(rj);
      Delta = (ladder[p] - ladder[p + 1]) * (Aj - Ai);
    } else {
      // conjugate = B = hopping sum;  Delta = (kappa_i-kappa_j)(B_i-B_j).
      const Real Bi = kappa_conjugate(ri);
      const Real Bj = kappa_conjugate(rj);
      Delta = (ladder[p] - ladder[p + 1]) * (Bi - Bj);
    }

    ++swap_attempts[p];
    bool accept = (Delta <= 0.0);
    if (!accept) {
      const double u = swap_rng.uniform(Rng::key(0x57A9, ++swap_counter));
      accept = (u < std::exp(-Delta));
    } else {
      ++swap_counter;  // keep the coin stream marching even on auto-accept
    }
    if (accept) { exchange_configs(ri, rj); ++swap_accepts[p]; }
    return accept;
  }

  // ---- conjugate energies, computed DIRECTLY from the u1.hpp action --------
  // BETA conjugate  A = +dS/dbeta = gauge_action(beta=1) = sum_plaq(1-cos).
  static Real beta_conjugate(const U1HMC<D>& r) {
    return gauge_action<D>(r.th, r.lat, /*beta=*/1.0);
  }
  // KAPPA conjugate B = -dS/dkappa = sum_{x,mu} 2 Re[conj(phi) e^{i q th} phi'].
  // Read off scalar_action: with kappa=1 and lambda=0 the onsite reduces to
  // sum|phi|^2; subtracting that leaves exactly -(- B) = B. Equivalently we sum
  // the hopping directly to avoid any onsite contamination.
  static Real kappa_conjugate(const U1HMC<D>& r) {
    const auto& phi = r.phi; const auto& th = r.th; const auto& lat = r.lat;
    Real B = 0.0;
    #pragma omp parallel for schedule(static) reduction(+:B)
    for (std::int64_t x = 0; x < lat.vol; ++x)
      for (int mu = 0; mu < D; ++mu) {
        const std::int64_t y = lat.neighbor_fwd(x, mu);
        const Complex ph = std::polar(1.0, r.q * th[x * D + mu]);  // u^q
        B += 2.0 * (std::conj(phi[x]) * ph * phi[y]).real();
      }
    return B;
  }

private:
  void set_axis_coupling(U1HMC<D>& r, Real v) {
    if (axis == TemperAxis::Beta) r.beta = v; else r.kappa = v;
  }
  void reapply_ladder() {
    for (std::size_t k = 0; k < replicas.size(); ++k) set_axis_coupling(replicas[k], ladder[k]);
  }

  // Exchange the link angles th and scalar field phi between two replicas. The
  // lattice/D are identical, so the flat storage is layout-compatible and a
  // std::swap of the underlying vectors is the exchange. Momenta and force
  // buffers are NOT swapped: momenta are refreshed next trajectory, and the
  // force buffers are scratch recomputed from (th,phi) inside each trajectory.
  static void exchange_configs(U1HMC<D>& a, U1HMC<D>& b) {
    std::swap(a.th,  b.th);
    std::swap(a.phi, b.phi);
  }

public:
  // Driver-friendly env toggle. Returns false iff U1_TEMPER is unset or set to a
  // falsy value (0/false/off/no), true otherwise. A driver can do:
  //   pt.enabled = U1ReplicaTempering<D>::env_enabled();
  static bool env_enabled() {
    const char* e = std::getenv("U1_TEMPER");
    if (!e) return false;
    std::string v(e);
    std::transform(v.begin(), v.end(), v.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return !(v.empty() || v == "0" || v == "false" || v == "off" || v == "no");
  }
};

}  // namespace u1
}  // namespace gh
