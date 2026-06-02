// Vacuum alignment: for each crystal-like subgroup H, build the H-singlet VEV, find
// stable multi-invariant couplings f_c that lock the gauge group to H, and report the
// spectrum. Validates: VEV is H-invariant; a stable minimum exists with exactly
// dim(SU(N))=N^2-1 Goldstone (zero) modes and positive massive modes (draft Obs 2).
// Run with arg "all" to also do the heavier SU(3) reps (Sigma 216/648/1080).
#include "check.hpp"
#include "action/vacuum_alignment.hpp"
#include "rep/rep_general.hpp"
#include <cstdio>

using namespace gh;

template <int N>
static void run(const Representation<N>& rep, const std::vector<Cmat<N>>& gens,
                const char* name, int nsamples) {
  auto H = close_group<N>(gens);
  auto sv = singlet_vev<N>(rep, H);
  CasimirChannels<N> ch(rep);
  auto res = find_stable_couplings<N>(ch, sv.vev, rep.real, nsamples, 12345);
  const int ngauge = N * N - 1;
  // zero modes = gauge Goldstones (broken generators, all eaten) + 1 global U(1) for a
  // complex rep + (multiplicity-1) flat MODULI when the H-singlet space is NOT unique.
  // Unique-singlet lockers (2T/2O/2I/Sigma, mult=1) add 0; Q8 in spin-2 has mult=2 -> one
  // biaxial-shape modulus -> a SOFT lock (the residual gauge group is still discrete Q8).
  const int expect_zero = (rep.real ? ngauge : ngauge + 1) + (sv.multiplicity - 1);

  char m[96];
  std::snprintf(m, sizeof m, "%s: VEV H-invariant", name);  CHECK(sv.invariance < 1e-8, m);
  std::snprintf(m, sizeof m, "%s: stable (PSD) locking couplings found", name); CHECK(res.stable, m);
  std::snprintf(m, sizeof m, "%s: zero modes = %d (expect %d Goldstones)", name, res.nzero, expect_zero);
  CHECK(res.nzero == expect_zero, m);

  std::printf("\n=== %s  |H|=%d  rep d=%d  singlet mult=%d  (inv %.1e) ===\n",
              name, (int)H.size(), rep.d, sv.multiplicity, sv.invariance);
  std::printf("  channels C2: "); for (Real c : ch.lambda) std::printf("%.4g ", c); std::printf("\n");
  std::printf("  f_c        : "); for (Real f : res.f) std::printf("%.4f ", f); std::printf("\n");
  std::printf("  mu^2 = %.4f   stable=%d   zero-modes=%d (expect %d)   negatives=%d   mass gap=%.4f\n",
              res.mu2, (int)res.stable, res.nzero, expect_zero, res.nneg, res.min_massive);
  // massive spectrum: positive eigenvalues grouped into degenerate multiplets
  std::vector<Real> pos; for (Real e : res.eig) if (e > 1e-5) pos.push_back(e);
  std::printf("  massive m^2 (multiplets): ");
  for (std::size_t i = 0; i < pos.size();) {
    std::size_t j = i; while (j < pos.size() && std::fabs(pos[j] - pos[i]) < 1e-3 * (1 + std::fabs(pos[i]))) ++j;
    std::printf("%.4f(x%zu) ", pos[i], j - i); i = j;
  }
  std::printf("\n");
}

int main(int argc, char** argv) {
  const bool all = (argc > 1 && std::string(argv[1]) == "all");
  // Treated as COMPLEX reps (the tensor basis is not the manifestly-real basis, so D is
  // complex unitary there): expected zero modes = (N^2-1) gauge Goldstones + 1 global U(1).
  std::printf("-- SU(2) binary polyhedral (2T spin-3, 2O spin-4) --\n");
  { GeneralRep<2> r({6});  run<2>(r, gens_2T(), "2T  (spin-3)", 800); }
  { GeneralRep<2> r({8});  run<2>(r, gens_2O(), "2O  (spin-4)", 800); }
  // Q8 control: spin-2 (l=2), COMPLEX like the other SU(2) entries (the GeneralRep tensor
  // basis is NOT manifestly real, so vacuum_alignment's real_rep=Re(phi) projection lands
  // on the wrong subspace and miscounts -- do NOT use :real here). The Q8-singlet space is
  // mult=2 (verified; character (1/8)[2(2j+1)+6(-1)^j]=2 at j=2): singlet_vev's orbit-rank
  // selector picks the BIAXIAL representative -> discrete Q8 gauge residual (all 3 W's
  // massive), but a biaxial-shape MODULUS survives -> expect_zero = 3 gauge + 1 U(1) + 1
  // modulus = 5. A SOFT lock (spin-2 is the natural-but-non-unique Q8 irrep).
  { GeneralRep<2> r({4}); run<2>(r, gens_Q8(), "Q8  (spin-2, soft/modulus)", 4000); }
  std::printf("\n-- SU(3) Sigma series --\n");
  { GeneralRep<3> r({4, 2}); run<3>(r, gens_Sigma108(), "Sigma(108) (2,2)", 600); }
  if (all) {  // heavier / larger coupling-search; 2I now buildable via the symmetric-rep basis
    { GeneralRep<2> r({12}); run<2>(r, gens_2I(), "2I  (spin-6)", 8000); }
    { GeneralRep<3> r({5, 1}); run<3>(r, gens_Sigma216(),  "Sigma(216) (4,1)", 600); }
    { GeneralRep<3> r({6});    run<3>(r, gens_Sigma1080(), "Sigma(1080) (6,0)", 600); }
    { GeneralRep<3> r({6, 3}); run<3>(r, gens_Sigma648(),  "Sigma(648) (3,3)", 400); }
  }
  return report("test_align");
}
