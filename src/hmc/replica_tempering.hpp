#pragma once
// Toggleable parallel tempering in beta (the gauge coupling) for the SU(N)->H
// gauge-Higgs HMC. Replicas k=0..M-1 sit at a fixed, increasing beta-ladder
//   beta_0 < beta_1 < ... < beta_{M-1}
// and share EVERYTHING else (representation, kappa, potential, lattice, lambda,
// tau, nmd, integrator, frozen_phi). Each replica is its own GaugeHiggsHMC<D,N>
// with a distinct RNG seed.
//
// Cycle (one call to step()):
//   1. every replica runs n_sweep independent HMC trajectories;
//   2. if enabled, attempt config swaps on adjacent pairs in two passes:
//        even pass: (0,1),(2,3),(4,5),...   then
//        odd  pass: (1,2),(3,4),(5,6),...
//      so that on average every replica gets a swap proposal with each neighbor.
//
// -------------------------------------------------------------------------
// SWAP ACCEPTANCE -- normalization derived from the ACTUAL gauge action.
// -------------------------------------------------------------------------
// src/action/gauge_wilson.hpp defines
//   S_g(U;beta) = (beta/N) * ( N * n_plaq - sumReTr(U) ),
//   avg_plaquette(U) = sumReTr(U) / (N * n_plaq),
// where sumReTr = sum over oriented plaquettes (mu<nu) of Re Tr U_plaq, and
// n_plaq = vol * D(D-1)/2.
//
// We want the canonical PT form  S_g = -beta * P  with  dS_g/dbeta = -P, i.e.
//   dS_g/dbeta = (1/N)(N*n_plaq - sumReTr) = n_plaq - sumReTr/N  ==  -P
//   =>  P(U) = sumReTr/N - n_plaq = n_plaq * ( avg_plaquette(U) - 1 ).
// (Only the beta-dependent gauge term enters: kappa/lambda/mu2 and the scalar
//  hopping are IDENTICAL across replicas, so they cancel exactly in Delta. The
//  scalar gauge-covariant hopping does involve U, but it is multiplied by the
//  same kappa in both replicas and the *fields being exchanged carry their own
//  hopping energy with them*, so the matter action is invariant under the
//  config swap and contributes nothing to Delta. Only S_g, whose coefficient
//  beta differs between replicas, matters.)
//
// For adjacent replicas (i with beta_i, j with beta_j), proposing to exchange
// their configurations, the Metropolis weight is
//   Delta = [S_i(C_j) + S_j(C_i)] - [S_i(C_i) + S_j(C_j)]
//         = (-beta_i*P_j - beta_j*P_i) - (-beta_i*P_i - beta_j*P_j)
//         = (beta_i - beta_j) * (P_i - P_j),
// accept with probability min(1, exp(-Delta)). On accept we exchange the gauge
// field U and the scalar field phi between replicas i and j (each replica keeps
// its OWN fixed beta_k; it is the CONFIGS that migrate). The constant -n_plaq
// in P cancels in (P_i - P_j); we keep it in P() so the value is the true PT
// energy, but Delta is unaffected.
//
// Momenta P (link) and pi (scalar) are refreshed from a Gaussian heatbath at
// the start of every trajectory (GaugeHiggsHMC::refresh_momenta), so they need
// not be swapped -- the next trajectory will redraw them regardless.
//
// -------------------------------------------------------------------------
// TOGGLE (safety property for merging).
// -------------------------------------------------------------------------
// `enabled == false`  => step() attempts NO swaps. The M replicas then evolve as
// M completely independent single-beta HMC streams (each with its own seed), so
// their per-replica time series can be merged exactly as ordinary independent
// runs. A driver may also flip the toggle from the environment (GH_TEMPER):
//   ReplicaTempering<D,N>::env_enabled() returns true unless GH_TEMPER is unset,
//   "0", "false", "off", or "no" (case-insensitive).
#include "hmc/gauge_higgs_hmc.hpp"
#include "action/gauge_wilson.hpp"
#include "measure/observables.hpp"
#include "core/fields.hpp"
#include "core/scalar_field.hpp"
#include <vector>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <cctype>
#include <string>

namespace gh {

template <int D, int N>
struct ReplicaTempering {
  std::vector<GaugeHiggsHMC<D, N>> replicas;  // one per ladder rung
  std::vector<Real> betas;                    // fixed beta of each rung (increasing)
  bool enabled = true;                        // master swap toggle
  int  n_sweep = 1;                           // HMC trajectories per replica per step()

  // Swap bookkeeping, indexed by adjacent pair p = (p, p+1), p = 0..M-2.
  std::vector<std::uint64_t> swap_attempts;   // proposals for pair p
  std::vector<std::uint64_t> swap_accepts;    // accepted exchanges for pair p

  // A dedicated RNG stream for the swap Metropolis coin (independent of the HMC
  // accept/reject streams so tempering never correlates with trajectory MC).
  Rng     swap_rng;
  std::uint64_t swap_counter = 0;

  // Construct M replicas on the given beta-ladder. All replicas share the same
  // representation R, kappa, potential, lattice extents, lambda, tau, nmd,
  // integrator and frozen_phi; only beta and the RNG seed differ. `seed0` seeds
  // replica 0; replica k is seeded seed0 + k*seed_stride (distinct streams).
  ReplicaTempering(const std::array<int, D>& extents,
                   const Representation<N>& R,
                   const std::vector<Real>& beta_ladder,
                   std::uint64_t seed0,
                   std::uint64_t seed_stride = 1000003ull,
                   std::uint64_t swap_seed   = 0x5A7AB1ECull)  // distinct swap stream
      : betas(beta_ladder), swap_rng(swap_seed) {
    if (beta_ladder.size() < 2)
      throw std::runtime_error("ReplicaTempering: need >=2 replicas");
    for (std::size_t k = 1; k < beta_ladder.size(); ++k)
      if (!(beta_ladder[k] > beta_ladder[k - 1]))
        throw std::runtime_error("ReplicaTempering: beta-ladder must be strictly increasing");

    replicas.reserve(beta_ladder.size());
    for (std::size_t k = 0; k < beta_ladder.size(); ++k) {
      replicas.emplace_back(extents, R, seed0 + k * seed_stride);
      replicas.back().beta = beta_ladder[k];
    }
    const std::size_t npair = beta_ladder.size() - 1;
    swap_attempts.assign(npair, 0);
    swap_accepts.assign(npair, 0);
  }

  std::size_t n_replicas() const { return replicas.size(); }
  std::size_t n_pairs()    const { return replicas.empty() ? 0 : replicas.size() - 1; }

  // --- shared-coupling setters (apply identical couplings to every replica) ---
  void set_kappa(Real kappa)            { for (auto& r : replicas) r.kappa = kappa; }
  void set_lambda(Real lambda)          { for (auto& r : replicas) r.lambda = lambda; }
  void set_tau(Real tau)                { for (auto& r : replicas) r.tau = tau; }
  void set_nmd(int nmd)                 { for (auto& r : replicas) r.nmd = nmd; }
  void set_potential(const OnsitePotential<N>* pot) { for (auto& r : replicas) r.potential = pot; }
  void set_frozen_phi(bool f) {
    for (auto& r : replicas) { r.frozen_phi = f; if (f) r.normalize_phi(); }
  }
  void set_integrator(Integrator integ) { for (auto& r : replicas) r.integ = integ; }

  // ---- accessors ----------------------------------------------------------
  GaugeHiggsHMC<D, N>&       replica(std::size_t k)       { return replicas.at(k); }
  const GaugeHiggsHMC<D, N>& replica(std::size_t k) const { return replicas.at(k); }
  Real beta(std::size_t k) const { return betas.at(k); }

  // Average plaquette P_pl = <(1/N) Re Tr U_pl> in [.., 1] for replica k.
  Real avg_plaq(std::size_t k) const { return avg_plaquette<D, N>(replicas.at(k).U); }

  // PT energy P(U) = n_plaq*(avg_plaquette - 1) for replica k, with S_g = -beta*P.
  Real plaq_energy(std::size_t k) const {
    const auto& r = replicas.at(k);
    const Real nplaq = static_cast<Real>(r.lat.n_plaq());
    return nplaq * (avg_plaquette<D, N>(r.U) - 1.0);
  }

  // Empirical swap-acceptance rate for adjacent pair p=(p,p+1).
  double pair_acceptance(std::size_t p) const {
    return swap_attempts.at(p) ? double(swap_accepts.at(p)) / double(swap_attempts.at(p)) : 0.0;
  }
  // HMC trajectory acceptance for replica k (delegates to GaugeHiggsHMC).
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
    GaugeHiggsHMC<D, N>& ri = replicas[p];
    GaugeHiggsHMC<D, N>& rj = replicas[p + 1];

    // P(U) = n_plaq*(avg_plaquette-1); the constant cancels in (P_i - P_j) but
    // we form the true energies for transparency.
    const Real nplaq = static_cast<Real>(ri.lat.n_plaq());  // same lattice both replicas
    const Real Pi = nplaq * (avg_plaquette<D, N>(ri.U) - 1.0);
    const Real Pj = nplaq * (avg_plaquette<D, N>(rj.U) - 1.0);

    // Delta = (beta_i - beta_j)*(P_i - P_j); accept w.p. min(1, exp(-Delta)).
    const Real Delta = (betas[p] - betas[p + 1]) * (Pi - Pj);

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

private:
  // Exchange the gauge field U and scalar field phi between two replicas.
  // Lattice/rep/d are identical, so the flat storage is layout-compatible and a
  // std::swap of the underlying vectors is the exchange. Momenta and force
  // buffers are NOT swapped: momenta are refreshed next trajectory, and the
  // force buffers are scratch recomputed from (U,phi) inside each trajectory.
  static void exchange_configs(GaugeHiggsHMC<D, N>& a, GaugeHiggsHMC<D, N>& b) {
    std::swap(a.U.u,      b.U.u);       // GaugeField contents (std::vector<Cmat<N>>)
    std::swap(a.phi.data, b.phi.data);  // scalar field contents (std::vector<Complex>)
  }

public:
  // Driver-friendly env toggle. Returns false iff GH_TEMPER is unset or set to a
  // falsy value (0/false/off/no), true otherwise. A driver can do:
  //   pt.enabled = ReplicaTempering<D,N>::env_enabled();
  static bool env_enabled() {
    const char* e = std::getenv("GH_TEMPER");
    if (!e) return false;
    std::string v(e);
    std::transform(v.begin(), v.end(), v.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return !(v.empty() || v == "0" || v == "false" || v == "off" || v == "no");
  }
};

}  // namespace gh
