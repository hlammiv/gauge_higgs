// Gauge-boson mass matrix M_ab (draft: "Gauge bosons spectra" + "Enhanced symmetries").
// With phi_B the H-singlet VEV and T^(sigma) the irrep generators,
//     M_ab = Re( phi_B^dagger T^(sigma)_b T^(sigma)_a phi_B ).
// For each of the 7 crystal-like subgroups H (SU(2): 2T,2O,2I; SU(3): Sigma108/216/648/1080)
// we verify the draft's three claims:
//   (i)   M_ab is proportional to the identity (off-diagonal ~ 0, all diagonal entries equal):
//         all dim(G) gauge bosons are degenerate, M is rank-full -> the CONTINUOUS gauge
//         symmetry is COMPLETELY broken (no surviving continuous subgroup).
//   (ii)  Tr M = C2(sigma)  (quadratic Casimir of the Higgs irrep), since
//         Tr M = sum_a phi^dag (T^a)^2 phi = phi^dag (sum_a (T^a)^2) phi = C2(sigma)|phi|^2.
//   (iii) the common eigenvalue m_A^2 = (2/dim_G) kappa g^2 C2(sigma). In units kappa g^2 = 1
//         the physical gauge mass^2 is m_A^2 = 2 * (common diagonal of M) = (2/dim_G) C2.
//         For SU(3): m_A^2 = (1/12)(p^2+q^2+pq+3p+3q) with dim_G=8; for SU(2): m_A^2 prop to
//         C2(j)=j(j+1) with dim_G=3.
// Pure group theory / linear algebra: no HMC, no lattice. See [[research-program]],
// [[gauge-higgs-validation-benchmarks]].
#include "check.hpp"
#include "action/vacuum_alignment.hpp"     // singlet_vev
#include "rep/rep_general.hpp"             // GeneralRep
#include "rep/finite_subgroups.hpp"        // gens_2T/.., gens_Sigma..., close_group
#include <cstdio>
#include <vector>

using namespace gh;

// One subgroup's result, for the synthesis table.
struct Row {
  std::string name;
  int N, dimG, dimR, Hsize;
  Real C2_expected, trM, common_eig, mA2;     // mA2 = 2 * common_eig (kappa g^2 = 1)
  Real offdiag_max, diag_spread, invariance;
};

template <int N>
static Row run(const Representation<N>& rep, const std::vector<Cmat<N>>& gens,
               const char* name, Real C2_expected) {
  const int dimG = N * N - 1;
  const int d = rep.d;

  // H-singlet VEV phi_B (unit vector, image of the group-averaging projector).
  auto H = close_group<N>(gens);
  auto sv = singlet_vev<N>(rep, H);
  const DVec& phi = sv.vev;

  // Precompute T^a phi for all a (each is a d-vector).
  std::vector<DVec> Tphi(dimG);
  for (int a = 0; a < dimG; ++a) Tphi[a] = rep.T[a] * phi;

  // M_ab = Re( phi^dag T^(sigma)_b T^(sigma)_a phi ) = Re( (T_b phi)^dag (T_a phi) )
  // since T_b is Hermitian: phi^dag T_b (T_a phi) = (T_b phi)^dag (T_a phi).
  std::vector<Real> M(static_cast<std::size_t>(dimG) * dimG, 0.0);
  for (int a = 0; a < dimG; ++a)
    for (int b = 0; b < dimG; ++b)
      M[static_cast<std::size_t>(a) * dimG + b] = dot(Tphi[b], Tphi[a]).real();

  // (i) proportional to identity: max |off-diagonal| and spread of diagonal entries.
  Real offdiag_max = 0.0, diag_min = 1e300, diag_max = -1e300, trM = 0.0;
  for (int a = 0; a < dimG; ++a) {
    Real daa = M[static_cast<std::size_t>(a) * dimG + a];
    trM += daa;
    diag_min = std::min(diag_min, daa);
    diag_max = std::max(diag_max, daa);
    for (int b = 0; b < dimG; ++b)
      if (a != b) offdiag_max = std::max(offdiag_max, std::fabs(M[static_cast<std::size_t>(a) * dimG + b]));
  }
  Real common_eig = trM / dimG;            // if M = c I then c = Tr M / dimG
  Real diag_spread = diag_max - diag_min;

  // also confirm M is symmetric (real-symmetric => its eigenvalues are exactly the diagonal
  // if it is a multiple of I; we check via spread + off-diagonal, plus an eig_sym spread).
  std::vector<Real> Msym = M;              // eig_sym wants real-symmetric (M is symmetric)
  std::vector<Real> ev = eig_sym(Msym, dimG);
  Real ev_min = 1e300, ev_max = -1e300;
  for (Real e : ev) { ev_min = std::min(ev_min, e); ev_max = std::max(ev_max, e); }

  Row r;
  r.name = name; r.N = N; r.dimG = dimG; r.dimR = d; r.Hsize = (int)H.size();
  r.C2_expected = C2_expected; r.trM = trM; r.common_eig = common_eig;
  r.mA2 = 2.0 * common_eig;                // physical gauge mass^2 in kappa g^2 = 1 units
  r.offdiag_max = offdiag_max; r.diag_spread = diag_spread; r.invariance = sv.invariance;

  char m[160];
  // VEV really is an H-singlet.
  std::snprintf(m, sizeof m, "%s: VEV is H-invariant (singlet)", name);
  CHECK(sv.invariance < 1e-8, m);
  // (i-a) off-diagonal ~ 0.
  std::snprintf(m, sizeof m, "%s: M off-diagonal ~ 0 (%.2e)", name, offdiag_max);
  CHECK(offdiag_max < 1e-9, m);
  // (i-b) all diagonal entries equal.
  std::snprintf(m, sizeof m, "%s: M diagonal entries all equal (spread %.2e)", name, diag_spread);
  CHECK(diag_spread < 1e-9, m);
  // (i-c) eigenvalues all equal (rank-full, nonzero) => continuous symmetry completely broken.
  std::snprintf(m, sizeof m, "%s: M eigenvalues all equal -> M = c I rank-full (eig spread %.2e)",
                name, ev_max - ev_min);
  CHECK((ev_max - ev_min) < 1e-9, m);
  std::snprintf(m, sizeof m, "%s: M is positive (rank-full, no surviving continuous gauge symmetry)", name);
  CHECK(ev_min > 1e-9, m);
  // (ii) Tr M = C2(sigma).
  std::snprintf(m, sizeof m, "%s: Tr M = C2(sigma) = %.6f (got %.6f)", name, C2_expected, trM);
  CHECK_CLOSE(trM, C2_expected, 1e-7, m);
  // (iii) common eigenvalue = (1/dimG) C2  => m_A^2 = (2/dimG) C2.
  std::snprintf(m, sizeof m, "%s: common eig = C2/dimG = %.6f", name, C2_expected / dimG);
  CHECK_CLOSE(common_eig, C2_expected / dimG, 1e-7, m);
  std::snprintf(m, sizeof m, "%s: m_A^2 = (2/dimG) C2 = %.6f", name, 2.0 * C2_expected / dimG);
  CHECK_CLOSE(r.mA2, 2.0 * C2_expected / dimG, 1e-7, m);

  std::printf("=== %-18s N=%d |H|=%-4d dim_R=%-3d dim_G=%d  C2=%.4f ===\n",
              name, N, (int)H.size(), d, dimG, C2_expected);
  std::printf("    Tr M = %.6f (expect C2 %.6f) | off-diag %.2e | diag spread %.2e | eig spread %.2e\n",
              trM, C2_expected, offdiag_max, diag_spread, ev_max - ev_min);
  std::printf("    common eig (C2/dim_G) = %.6f | m_A^2 = 2*eig = %.6f | m_A = %.6f  (kappa g^2 = 1)\n\n",
              common_eig, r.mA2, std::sqrt(std::max(0.0, r.mA2)));
  return r;
}

int main() {
  std::printf("Gauge-boson mass matrix M_ab = Re(phi^dag T_b T_a phi), H-singlet VEV.\n");
  std::printf("(kappa g^2 = 1 units; m_A^2 = (2/dim_G) C2(sigma))\n\n");

  std::vector<Row> rows;

  // ----- SU(2) binary polyhedral. C2(j)=j(j+1). dim_G = 3. -----
  // 2T: spin-3 {6} d=7  C2 = 3*4 = 12
  { GeneralRep<2> r({6});  rows.push_back(run<2>(r, gens_2T(),  "BT  (2T, spin-3)", 3.0 * 4.0)); }
  // 2O: spin-4 {8} d=9  C2 = 4*5 = 20
  { GeneralRep<2> r({8});  rows.push_back(run<2>(r, gens_2O(),  "BO  (2O, spin-4)", 4.0 * 5.0)); }
  // 2I: spin-6 {12} d=13 C2 = 6*7 = 42
  { GeneralRep<2> r({12}); rows.push_back(run<2>(r, gens_2I(),  "BI  (2I, spin-6)", 6.0 * 7.0)); }

  // ----- SU(3) Sigma series. C2(p,q) = (1/3)(p^2+q^2+pq+3p+3q). dim_G = 8. -----
  auto c2_su3 = [](int p, int q) { return (p*p + q*q + p*q + 3*p + 3*q) / 3.0; };
  // Sigma108: (2,2) {4,2} d=27   bracket 24 -> C2 = 8
  { GeneralRep<3> r({4, 2}); rows.push_back(run<3>(r, gens_Sigma108(),  "Sigma108 (2,2)", c2_su3(2,2))); }
  // Sigma216: (4,1) {5,1} d=35   bracket 36 -> C2 = 12
  { GeneralRep<3> r({5, 1}); rows.push_back(run<3>(r, gens_Sigma216(),  "Sigma216 (4,1)", c2_su3(4,1))); }
  // Sigma648: (3,3) {6,3} d=64   bracket 45 -> C2 = 15
  { GeneralRep<3> r({6, 3}); rows.push_back(run<3>(r, gens_Sigma648(),  "Sigma648 (3,3)", c2_su3(3,3))); }
  // Sigma1080: (6,0) {6} d=28    bracket 54 -> C2 = 18
  { GeneralRep<3> r({6});    rows.push_back(run<3>(r, gens_Sigma1080(), "Sigma1080 (6,0)", c2_su3(6,0))); }

  // Cross-check the SU(3) draft brackets (p^2+q^2+pq+3p+3q) explicitly.
  struct PQ { const char* nm; int p, q, bracket; };
  PQ pqs[] = {{"(2,2)", 2, 2, 24}, {"(4,1)", 4, 1, 36}, {"(3,3)", 3, 3, 45}, {"(6,0)", 6, 0, 54}};
  for (auto& x : pqs) {
    int got = x.p*x.p + x.q*x.q + x.p*x.q + 3*x.p + 3*x.q;
    char m[96]; std::snprintf(m, sizeof m, "SU(3) bracket %s = %d", x.nm, x.bracket);
    CHECK(got == x.bracket, m);
  }

  // Synthesis table (for type-I/II comparison vs scalar m_H).
  std::printf("--- synthesis: gauge-boson mass m_A per subgroup (kappa g^2 = 1) ---\n");
  std::printf("%-18s %5s %5s %6s %10s %10s %10s\n", "subgroup", "dimG", "dimR", "C2", "m_A^2", "m_A", "|H|");
  for (const auto& r : rows)
    std::printf("%-18s %5d %5d %6.3f %10.5f %10.5f %6d\n",
                r.name.c_str(), r.dimG, r.dimR, r.C2_expected, r.mA2, std::sqrt(std::max(0.0, r.mA2)), r.Hsize);

  return report("test_draft_gaugemass");
}
