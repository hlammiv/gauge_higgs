// Gauge-invariance proof for the gauge-boson vector operator (src/measure/gaugeboson_op.hpp).
//
// We build a 4^4 lattice, a HOT SU(N) gauge field, and a random normalized rep-R
// scalar phi(x) per site. We compute all electric components O_{0 i}(x), i=1..D-1.
// Then we apply a RANDOM per-site SU(N) gauge transform g(x) = expi(alg_to_mat(rand)):
//     U_mu(x) -> g(x) U_mu(x) g(x+mu)^dag ,   phi(x) -> D^(R)(g(x)) phi(x) = rep.rotate(g, phi),
// recompute every O_{0 i}(x), and report the max |O_after - O_before| over all (x,i).
//
// The operator is a gauge SINGLET (n^a and F^a both rotate by Ad(g), the dotAlg
// contraction is invariant), so the difference must be ~roundoff (~4e-15). We PASS a
// rep iff err < 1e-9 (invariant) AND max|O| > 1e-8 (the operator is nontrivial -- a
// guard against the load-bearing "identically zero" bug if the -i factor were dropped).
//
// Reps: AdjointRep<2> (real, N-ality 0), GeneralRep<2>({4}) = SU(2) spin-2 (d=5),
// GeneralRep<3>({2}) = SU(3) symmetric (6, d=6).
#include "check.hpp"
#include "core/geometry.hpp"
#include "core/fields.hpp"
#include "core/rng.hpp"
#include "core/linalg.hpp"
#include "group/sun.hpp"
#include "group/algebra.hpp"
#include "rep/representation.hpp"
#include "rep/rep_adjoint.hpp"
#include "rep/rep_general.hpp"
#include "measure/gaugeboson_op.hpp"
#include <vector>
#include <cstdio>
#include <cmath>
#include <cstdint>

using namespace gh;

// A random unit-norm rep-R scalar field, one DVec(d) per site (deterministic in seed).
template <int N>
static std::vector<DVec> random_phi(const Representation<N>& rep, std::int64_t vol,
                                    const Rng& rng, std::uint64_t stream) {
  const int d = rep.d;
  std::vector<DVec> phi(static_cast<std::size_t>(vol), DVec(d));
  for (std::int64_t s = 0; s < vol; ++s) {
    DVec v(d);
    for (int k = 0; k < d; ++k) {
      const Real re = rng.gauss(Rng::key(stream, s, k, 0));
      const Real im = rng.gauss(Rng::key(stream, s, k, 1));
      v(k) = Complex(re, im);
    }
    const Real nrm = std::sqrt(v.norm2());
    const Real inv = nrm > 0 ? 1.0 / nrm : 0.0;
    for (int k = 0; k < d; ++k) v(k) = v(k) * Complex(inv, 0.0);
    phi[static_cast<std::size_t>(s)] = v;
  }
  return phi;
}

// A random SU(N) gauge transform g(x) per site, g(x) = expi(alg_to_mat(random AlgVec)).
template <int N>
static std::vector<Cmat<N>> random_gauge(std::int64_t vol, const Rng& rng,
                                         std::uint64_t stream) {
  std::vector<Cmat<N>> g(static_cast<std::size_t>(vol));
  for (std::int64_t s = 0; s < vol; ++s) {
    AlgVec<N> w{};
    for (int a = 0; a < n_gen<N>(); ++a) w[a] = rng.gauss(Rng::key(stream, s, a, 7));
    g[static_cast<std::size_t>(s)] = expi<N>(alg_to_mat<N>(w));
  }
  return g;
}

// Compute every electric component O_{0 i}(x), i=1..D-1, in a flat vector (x*(D-1)+(i-1)).
template <int D, int N>
static std::vector<Real> all_electric(const Representation<N>& rep,
                                      const GaugeField<D, N>& U,
                                      const std::vector<DVec>& phi) {
  const Lattice<D>& lat = *U.lat;
  std::vector<Real> O(static_cast<std::size_t>(lat.vol) * (D - 1), 0.0);
  for (std::int64_t s = 0; s < lat.vol; ++s)
    for (int i = 1; i < D; ++i)
      O[static_cast<std::size_t>(s) * (D - 1) + (i - 1)] =
          measure::op_electric<D, N>(rep, U, phi[static_cast<std::size_t>(s)], s, i);
  return O;
}

// Build the gauge-transformed gauge field and scalar field in place (new copies).
//   U'_mu(x) = g(x) U_mu(x) g(x+mu)^dag ,   phi'(x) = rep.rotate(g(x), phi(x)).
template <int D, int N>
static void apply_gauge(const Representation<N>& rep, const GaugeField<D, N>& U,
                        const std::vector<DVec>& phi, const std::vector<Cmat<N>>& g,
                        GaugeField<D, N>& Uout, std::vector<DVec>& phiout) {
  const Lattice<D>& lat = *U.lat;
  for (std::int64_t s = 0; s < lat.vol; ++s) {
    const Cmat<N>& gs = g[static_cast<std::size_t>(s)];
    for (int mu = 0; mu < D; ++mu) {
      const std::int64_t s_pmu = lat.neighbor_fwd(s, mu);
      Uout(s, mu) = gs * U(s, mu) * g[static_cast<std::size_t>(s_pmu)].dagger();
    }
    phiout[static_cast<std::size_t>(s)] = rep.rotate(gs, phi[static_cast<std::size_t>(s)]);
  }
}

template <int D, int N>
static void run_rep(const Representation<N>& rep, const char* tag,
                    const Lattice<D>& lat, std::uint64_t seed) {
  Rng rng(seed);

  // Hot gauge field + random unit scalar.
  GaugeField<D, N> U(lat);
  U.hot(rng, 1.0);
  std::vector<DVec> phi = random_phi<N>(rep, lat.vol, rng, /*stream=*/101);

  // Operator BEFORE the gauge transform.
  std::vector<Real> O_before = all_electric<D, N>(rep, U, phi);

  // Random per-site SU(N) gauge transform, applied to links and scalars.
  std::vector<Cmat<N>> g = random_gauge<N>(lat.vol, rng, /*stream=*/202);
  GaugeField<D, N> Ug(lat);
  std::vector<DVec> phig(static_cast<std::size_t>(lat.vol), DVec(rep.d));
  apply_gauge<D, N>(rep, U, phi, g, Ug, phig);

  // Operator AFTER the gauge transform.
  std::vector<Real> O_after = all_electric<D, N>(rep, Ug, phig);

  // Max invariance error and max |O| (operator nontriviality).
  Real err = 0.0, maxO = 0.0;
  for (std::size_t k = 0; k < O_before.size(); ++k) {
    err = std::max(err, std::fabs(O_after[k] - O_before[k]));
    maxO = std::max(maxO, std::fabs(O_before[k]));
  }

  const bool invariant = (err < 1e-9);
  const bool nontrivial = (maxO > 1e-8);
  std::printf("  %-22s d=%-2d  max|O|=%.3e  max|O_after-O_before|=%.3e  -> %s\n",
              tag, rep.d, maxO, err, (invariant && nontrivial) ? "PASS" : "FAIL");

  char m[160];
  std::snprintf(m, sizeof m, "%s: O_{0i} gauge invariant (err %.2e < 1e-9)", tag, err);
  CHECK(invariant, m);
  std::snprintf(m, sizeof m, "%s: O_{0i} nontrivial (max|O| %.2e > 1e-8, -i factor present)", tag, maxO);
  CHECK(nontrivial, m);
}

int main() {
  std::printf("Gauge-boson vector operator O_{0i} = sum_a n^a F^a_{0i}: gauge-invariance test\n");
  std::printf("(4^4 lattice, hot SU(N) links, random unit rep-R scalar, random per-site SU(N) gauge xfm)\n\n");

  const Lattice<4> lat({4, 4, 4, 4});

  // SU(2) adjoint (real, N-ality 0).
  { AdjointRep<2> rep;        run_rep<4, 2>(rep, "AdjointRep<2>", lat, 0xA11CE001ULL); }
  // SU(2) spin-2 = GeneralRep<2>({4}), d=5.
  { GeneralRep<2> rep({4});   run_rep<4, 2>(rep, "GeneralRep<2>({4}) spin2", lat, 0xB0B0B0B0ULL); }
  // SU(3) symmetric 6 = GeneralRep<3>({2}), d=6.
  { GeneralRep<3> rep({2});   run_rep<4, 3>(rep, "GeneralRep<3>({2})", lat, 0xC0FFEE99ULL); }

  std::printf("\n");
  return report("test_gaugeboson_op");
}
