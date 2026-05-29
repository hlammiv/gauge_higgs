// Tests for the general multi-invariant quartic potential (CasimirChannels):
// exact channel Casimir eigenvalues for SU(2) (= J(J+1)), projector completeness
// (sum_c P_c = I), and finite-difference check of dV/dphibar.
#include "check.hpp"
#include "action/scalar_invariants.hpp"
#include "hmc/gauge_higgs_hmc.hpp"
#include "rep/rep_general.hpp"
#include "rep/rep_adjoint.hpp"
#include <random>

using namespace gh;

template <int D>
static std::array<int, D> cube(int n) { std::array<int, D> L{}; for (int mu = 0; mu < D; ++mu) L[mu] = n; return L; }

// Validate the multi-invariant potential WIRED INTO the HMC: integrated scalar-force
// finite difference, combined-HMC reversibility, and <exp(-dH)>~1.
template <int Dd, int N>
static void test_hmc_multiinv(const Representation<N>& rep, const CasimirChannels<N>& ch,
                              std::uint64_t seed, const char* tag) {
  std::vector<Real> f(ch.n_channels());
  for (int c = 0; c < ch.n_channels(); ++c) f[c] = 0.5 + 0.2 * c;  // all >0 => bounded below
  MultiInvariantPotential<N> pot(ch, f, 1.0);
  const Real kappa = 0.3;

  // (1) integrated scalar-force FD: F = -2 dS_H/dphi^* with the multi-invariant potential.
  Lattice<Dd> lat(cube<Dd>(4));
  GaugeField<Dd, N> U(lat); Rng rng(seed); U.hot(rng, 0.5);
  CVecField<Dd> phi(lat, rep.d); phi.gaussian(rng, 55, rep.real, 0.6);
  CVecField<Dd> F(lat, rep.d); scalar_force<Dd, N>(phi, U, rep, kappa, pot, F);
  const Real eps = 1e-5; Real ws = 0.0;
  std::int64_t sites[] = {0, 7, lat.vol - 1};
  for (std::int64_t s : sites)
    for (int k = 0; k < rep.d; ++k) {
      const std::size_t idx = static_cast<std::size_t>(s) * rep.d + k;
      const Complex o = phi.data[idx];
      phi.data[idx] = o + Complex(eps, 0); Real Sp = scalar_action<Dd, N>(phi, U, rep, kappa, pot);
      phi.data[idx] = o - Complex(eps, 0); Real Sm = scalar_action<Dd, N>(phi, U, rep, kappa, pot);
      phi.data[idx] = o + Complex(0, eps); Real Tp = scalar_action<Dd, N>(phi, U, rep, kappa, pot);
      phi.data[idx] = o - Complex(0, eps); Real Tm = scalar_action<Dd, N>(phi, U, rep, kappa, pot);
      phi.data[idx] = o;
      ws = std::max(ws, std::fabs(-(Sp - Sm) / (2 * eps) - F.data[idx].real()));
      ws = std::max(ws, std::fabs(-(Tp - Tm) / (2 * eps) - F.data[idx].imag()));
    }
  char m[96]; std::snprintf(m, sizeof m, "multi-invariant scalar force FD (%s)", tag);
  CHECK_CLOSE(ws, 0.0, 1e-5, m);

  // (2) combined HMC reversibility with the potential plugged in.
  GaugeHiggsHMC<Dd, N> hmc(cube<Dd>(4), rep, seed + 1);
  hmc.beta = 2.3; hmc.kappa = kappa; hmc.tau = 1.0; hmc.nmd = 16; hmc.potential = &pot;
  hmc.reunit_each_traj = false;
  hmc.U.hot(hmc.rng, 0.4); hmc.phi.gaussian(hmc.rng, 1, rep.real, 0.4);
  hmc.refresh_momenta();
  std::vector<Cmat<N>> U0 = hmc.U.u; std::vector<Complex> phi0 = hmc.phi.data;
  hmc.md_evolve();
  for (auto& v : hmc.P.p) for (int a = 0; a < n_gen<N>(); ++a) v[a] = -v[a];
  for (auto& z : hmc.pi.data) z = -z;
  hmc.md_evolve();
  Real wr = 0.0;
  for (std::size_t i = 0; i < U0.size(); ++i) wr = std::max(wr, (hmc.U.u[i] - U0[i]).fnorm());
  for (std::size_t i = 0; i < phi0.size(); ++i) wr = std::max(wr, std::abs(hmc.phi.data[i] - phi0[i]));
  std::snprintf(m, sizeof m, "multi-invariant HMC reversibility (%s)", tag);
  CHECK_CLOSE(wr, 0.0, 1e-9, m);

  // (3) <exp(-dH)> ~ 1 over a short run.
  GaugeHiggsHMC<Dd, N> hmc2(cube<Dd>(4), rep, seed + 2);
  hmc2.beta = 2.3; hmc2.kappa = kappa; hmc2.tau = 1.0; hmc2.nmd = 24; hmc2.potential = &pot;
  hmc2.U.hot(hmc2.rng, 0.3); hmc2.phi.gaussian(hmc2.rng, 2, rep.real, 0.3);
  const int ntraj = 250; double s = 0.0;
  for (int t = 0; t < ntraj; ++t) { hmc2.trajectory(); s += std::exp(-hmc2.last_dH); }
  const double mean = s / ntraj;
  std::snprintf(m, sizeof m, "multi-invariant HMC <exp(-dH)>~1 (%s) got %.4f acc %.2f", tag, mean, hmc2.acceptance());
  CHECK(std::fabs(mean - 1.0) < 0.15, m);
}

// SU(2) spin-j (Young rows {2j}): channels of R x Rbar = spin-J, J=0..2j, C2(J)=J(J+1).
static void test_su2_channels() {
  { GeneralRep<2> r({6}); CasimirChannels<2> ch(r);   // {6} = spin-3 (BT rep), dim 7
    std::vector<Real> exp{0, 2, 6, 12, 20, 30, 42};
    char m[48]; std::snprintf(m, sizeof m, "SU(2) spin-3 (dim %d): %d channels", r.d, ch.n_channels());
    CHECK(r.d == 7 && ch.n_channels() == 7, m);
    Real worst = 0; for (int c = 0; c < ch.n_channels() && c < 7; ++c) worst = std::max(worst, std::fabs(ch.lambda[c] - exp[c]));
    CHECK_CLOSE(worst, 0.0, 1e-5, "SU(2) spin-3 channel C2 == J(J+1)");
  }
  { GeneralRep<2> r({8}); CasimirChannels<2> ch(r);   // {8} = spin-4 (BO rep), dim 9
    std::vector<Real> exp{0, 2, 6, 12, 20, 30, 42, 56, 72};
    CHECK(r.d == 9 && ch.n_channels() == 9, "SU(2) spin-4 (dim 9): 9 channels");
    Real worst = 0; for (int c = 0; c < ch.n_channels() && c < 9; ++c) worst = std::max(worst, std::fabs(ch.lambda[c] - exp[c]));
    CHECK_CLOSE(worst, 0.0, 1e-5, "SU(2) spin-4 channel C2 == J(J+1)");
  }
}

template <int N>
static void test_completeness(const Representation<N>& rep, const char* tag) {
  CasimirChannels<N> ch(rep);
  std::mt19937_64 rng(7);
  std::normal_distribution<Real> g(0.0, 1.0);
  DMat M(rep.d, rep.d);
  for (int a = 0; a < rep.d; ++a) for (int b = 0; b < rep.d; ++b) M(a, b) = Complex(g(rng), g(rng));
  DMat sum(rep.d, rep.d);
  for (int c = 0; c < ch.n_channels(); ++c) sum = sum + ch.apply_proj(c, M);
  char m[64]; std::snprintf(m, sizeof m, "completeness sum_c P_c = I (%s, %d chan)", tag, ch.n_channels());
  CHECK_CLOSE(fnorm(sum - M), 0.0, 1e-7, m);
}

// dV/dphibar vs finite difference: dV/dRe(phi_e)=2Re(g_e), dV/dIm(phi_e)=2Im(g_e).
template <int N>
static void test_gradient_fd(const Representation<N>& rep, std::uint64_t seed, const char* tag) {
  CasimirChannels<N> ch(rep);
  std::mt19937_64 rng(seed);
  std::normal_distribution<Real> gd(0.0, 0.6);
  DVec phi(rep.d);
  for (int e = 0; e < rep.d; ++e) phi(e) = Complex(gd(rng), gd(rng));  // generic complex
  std::vector<Real> f(ch.n_channels());
  for (auto& x : f) x = gd(rng);
  const Real mu2 = 0.7;
  DVec gvec = ch.dV_dphibar(phi, f, mu2);
  const Real eps = 1e-5; Real worst = 0;
  for (int e = 0; e < rep.d; ++e) {
    Complex o = phi(e);
    phi(e) = o + Complex(eps, 0); Real Vp = ch.value(phi, f, mu2);
    phi(e) = o - Complex(eps, 0); Real Vm = ch.value(phi, f, mu2);
    phi(e) = o + Complex(0, eps); Real Vpi = ch.value(phi, f, mu2);
    phi(e) = o - Complex(0, eps); Real Vmi = ch.value(phi, f, mu2);
    phi(e) = o;
    worst = std::max(worst, std::fabs((Vp - Vm) / (2 * eps) - 2.0 * gvec(e).real()));
    worst = std::max(worst, std::fabs((Vpi - Vmi) / (2 * eps) - 2.0 * gvec(e).imag()));
  }
  char m[64]; std::snprintf(m, sizeof m, "dV/dphibar FD (%s)", tag);
  CHECK_CLOSE(worst, 0.0, 1e-5, m);
}

int main() {
  std::printf("-- SU(2) channel Casimir eigenvalues --\n");
  test_su2_channels();
  std::printf("-- projector completeness --\n");
  { GeneralRep<2> r3({6}); test_completeness<2>(r3, "SU2 spin-3"); }
  { AdjointRep<3> a3;      test_completeness<3>(a3, "SU3 adj"); }
  { GeneralRep<3> r10({3, 0}); test_completeness<3>(r10, "SU3 10"); }
  std::printf("-- gradient finite difference --\n");
  { GeneralRep<2> r3({6}); test_gradient_fd<2>(r3, 101, "SU2 spin-3"); }
  { GeneralRep<3> r10({3, 0}); test_gradient_fd<3>(r10, 102, "SU3 10"); }
  std::printf("-- multi-invariant potential wired into HMC --\n");
  { GeneralRep<2> r({3}); CasimirChannels<2> ch(r); test_hmc_multiinv<3, 2>(r, ch, 301, "SU2 spin-3/2 d=4"); }
  return report("test_invariants");
}
