// Assertions for the character-based H-singlet multiplicity enumeration.
// Validates the draft's Tab. su2subduction / su3subduction (singlet-containing
// irreps, marked red there) via m_{sigma->tr} = (1/|H|) sum_h chi^sigma(h).
//
//   (i)   lowest unique singlets: 2T j=3, 2O j=4, 2I j=6 (each m=1);
//         Sigma(108) (2,2), Sigma(216) (4,1), Sigma(648) (3,3),
//         Sigma(1080) (6,0) -- each m=1.
//   (ii)  multiplicity-2 (two singlet directions): 2T j=6 (draft "2 A_0" -> m=2).
//   (iii) additional singlet-containing irreps from the notes:
//         Sigma(108)&Sigma(216): (7,1),(4,4),(1,7); Sigma(648): (4,4);
//         Sigma(1080): (4,4).
// Also: full cross-check of every m vs the draft tables (j<=10, p+q<=6) and the
// integrality of the raw projector value.
#define SINGLETS_NO_MAIN
#include "check.hpp"
#include "../src/singlets.cpp"
#include <cstdio>

using namespace gh;

int main() {
  auto T    = close_group<2>(gens_2T());
  auto O    = close_group<2>(gens_2O());
  auto I    = close_group<2>(gens_2I());
  auto S108 = close_group<3>(gens_Sigma108());
  auto S216 = close_group<3>(gens_Sigma216());
  auto S648 = close_group<3>(gens_Sigma648());
  auto S1080 = close_group<3>(gens_Sigma1080());

  // group orders (sanity: confirms close_group gives the right finite group).
  CHECK((int)T.size() == 24,    "|2T| = 24");
  CHECK((int)O.size() == 48,    "|2O| = 48");
  CHECK((int)I.size() == 120,   "|2I| = 120");
  CHECK((int)S108.size() == 108,  "|Sigma(108)| = 108");
  CHECK((int)S216.size() == 216,  "|Sigma(216)| = 216");
  CHECK((int)S648.size() == 648,  "|Sigma(648)| = 648");
  CHECK((int)S1080.size() == 1080, "|Sigma(1080)| = 1080");

  // m(trivial irrep) = 1 for every group (sin(theta)/sin(theta) limit, etc.).
  CHECK(su2_singlet_mult(0, T) == 1, "2T:  j=0 trivial m=1");
  CHECK(su2_singlet_mult(0, O) == 1, "2O:  j=0 trivial m=1");
  CHECK(su2_singlet_mult(0, I) == 1, "2I:  j=0 trivial m=1");
  CHECK(su3_singlet_mult(0, 0, S108)  == 1, "S108:  (0,0) trivial m=1");
  CHECK(su3_singlet_mult(0, 0, S216)  == 1, "S216:  (0,0) trivial m=1");
  CHECK(su3_singlet_mult(0, 0, S648)  == 1, "S648:  (0,0) trivial m=1");
  CHECK(su3_singlet_mult(0, 0, S1080) == 1, "S1080: (0,0) trivial m=1");

  // --- (i) lowest non-trivial unique singlets, each m=1 ---
  CHECK(su2_singlet_mult(6, T)  == 1, "(i) 2T: j=3 (d7)  m=1");   // twoj=6
  CHECK(su2_singlet_mult(8, O)  == 1, "(i) 2O: j=4 (d9)  m=1");   // twoj=8
  CHECK(su2_singlet_mult(12, I) == 1, "(i) 2I: j=6 (d13) m=1");   // twoj=12
  CHECK(su3_singlet_mult(2, 2, S108)  == 1, "(i) S108:  (2,2) m=1");
  CHECK(su3_singlet_mult(4, 1, S216)  == 1, "(i) S216:  (4,1) m=1");
  CHECK(su3_singlet_mult(3, 3, S648)  == 1, "(i) S648:  (3,3) m=1");
  CHECK(su3_singlet_mult(6, 0, S1080) == 1, "(i) S1080: (6,0) m=1");

  // these lowest irreps must be the FIRST non-trivial singlet for their group:
  // no smaller non-trivial j / (p,q) contains a singlet of that same group.
  CHECK(su2_singlet_mult(2, T) == 0 && su2_singlet_mult(4, T) == 0 &&
        su2_singlet_mult(5, T) == 0,                 "(i) 2T: no singlet for 0<j<3 (int j)");
  for (int tj = 1; tj <= 7; ++tj) if (tj != 8) CHECK(su2_singlet_mult(tj, O) == 0,
        "(i) 2O: no singlet for 0<2j<8");
  for (int tj = 1; tj <= 11; ++tj) CHECK(su2_singlet_mult(tj, I) == 0,
        "(i) 2I: no singlet for 0<2j<12");

  // --- (ii) multiplicity-2: 2T j=6 has TWO singlet directions ("2 A_0") ---
  CHECK(su2_singlet_mult(12, T) == 2, "(ii) 2T: j=6 (d13) m=2 (two A_0 directions)");

  // --- (iii) additional singlet-containing irreps flagged in the notes ---
  CHECK(su3_singlet_mult(7, 1, S108) >= 1, "(iii) S108: (7,1) contains singlet");
  CHECK(su3_singlet_mult(4, 4, S108) >= 1, "(iii) S108: (4,4) contains singlet");
  CHECK(su3_singlet_mult(1, 7, S108) >= 1, "(iii) S108: (1,7) contains singlet");
  CHECK(su3_singlet_mult(7, 1, S216) >= 1, "(iii) S216: (7,1) contains singlet");
  CHECK(su3_singlet_mult(4, 4, S216) >= 1, "(iii) S216: (4,4) contains singlet");
  CHECK(su3_singlet_mult(1, 7, S216) >= 1, "(iii) S216: (1,7) contains singlet");
  CHECK(su3_singlet_mult(4, 4, S648) >= 1, "(iii) S648: (4,4) contains singlet");
  CHECK(su3_singlet_mult(4, 4, S1080) >= 1, "(iii) S1080: (4,4) contains singlet");

  // --- full cross-check vs draft Tab. su2subduction (integer j; count of A_0/A_1) ---
  // 2T: A_0 multiplicities for j = 0..10.
  const int t2T[11]  = {1, 0, 0, 1, 1, 0, 2, 1, 1, 2, 2};
  // 2O: A_1 multiplicities for j = 0..10.
  const int t2O[11]  = {1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1};
  // 2I: A_0 multiplicities for j = 0..10.
  const int t2I[11]  = {1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1};
  for (int j = 0; j <= 10; ++j) {
    char m[64];
    std::snprintf(m, sizeof m, "su2subduction 2T j=%d : m=%d", j, t2T[j]);
    CHECK(su2_singlet_mult(2 * j, T) == t2T[j], m);
    std::snprintf(m, sizeof m, "su2subduction 2O j=%d : m=%d", j, t2O[j]);
    CHECK(su2_singlet_mult(2 * j, O) == t2O[j], m);
    std::snprintf(m, sizeof m, "su2subduction 2I j=%d : m=%d", j, t2I[j]);
    CHECK(su2_singlet_mult(2 * j, I) == t2I[j], m);
  }

  // --- full cross-check vs draft su3_table.tex (count of 1^(0); p+q<=6) ---
  // Indexed [p][q]; -1 = (p+q>6, not in the table). Multiplicity of trivial 1^(0).
  // From su3_table.tex: red 1^(0) entries.
  struct PQM { int p, q, m; };
  const std::vector<PQM> s108 = {
    {0,0,1},{2,2,1},{4,1,1},{1,4,1},{3,3,2},{6,0,2},{0,6,2}};
  const std::vector<PQM> s216 = {
    {0,0,1},{4,1,1},{1,4,1},{6,0,1},{3,3,1},{0,6,1}};
  const std::vector<PQM> s648 = {
    {0,0,1},{3,3,1}};
  const std::vector<PQM> s1080 = {
    {0,0,1},{6,0,1},{0,6,1}};
  auto verify = [&](const std::vector<PQM>& tab, const std::vector<Cmat<3>>& H,
                    const char* nm) {
    for (const auto& e : tab) {
      char m[80];
      std::snprintf(m, sizeof m, "su3subduction %s (%d,%d): m=%d", nm, e.p, e.q, e.m);
      CHECK(su3_singlet_mult(e.p, e.q, H) == e.m, m);
    }
  };
  verify(s108, S108, "S108");
  verify(s216, S216, "S216");
  verify(s648, S648, "S648");
  verify(s1080, S1080, "S1080");

  // every OTHER (p,q) with p+q<=6 NOT listed above must have m=0.
  auto in = [](const std::vector<PQM>& tab, int p, int q) {
    for (const auto& e : tab) if (e.p == p && e.q == q) return true; return false;
  };
  for (int sgrp = 0; sgrp < 4; ++sgrp) {
    const std::vector<PQM>& tab = (sgrp==0?s108:sgrp==1?s216:sgrp==2?s648:s1080);
    const std::vector<Cmat<3>>& H = (sgrp==0?S108:sgrp==1?S216:sgrp==2?S648:S1080);
    const char* nm = (sgrp==0?"S108":sgrp==1?"S216":sgrp==2?"S648":"S1080");
    for (int s = 0; s <= 6; ++s)
      for (int p = 0; p <= s; ++p) {
        const int q = s - p;
        if (in(tab, p, q)) continue;
        char m[80];
        std::snprintf(m, sizeof m, "su3subduction %s (%d,%d): NO singlet (m=0)", nm, p, q);
        CHECK(su3_singlet_mult(p, q, H) == 0, m);
      }
  }

  // --- integrality of the raw (un-rounded) projector value (FP sanity) ---
  for (int tj = 0; tj <= 20; ++tj) {
    Real r = su2_singlet_raw(tj, T);
    CHECK(std::fabs(r - std::llround(r)) < 1e-6, "2T raw m integral");
    r = su2_singlet_raw(tj, O);
    CHECK(std::fabs(r - std::llround(r)) < 1e-6, "2O raw m integral");
    r = su2_singlet_raw(tj, I);
    CHECK(std::fabs(r - std::llround(r)) < 1e-6, "2I raw m integral");
  }
  for (int s = 0; s <= 8; ++s)
    for (int p = 0; p <= s; ++p) {
      const int q = s - p;
      Real r = su3_singlet_raw(p, q, S108);
      CHECK(std::fabs(r - std::llround(r)) < 1e-5, "S108 raw m integral");
      r = su3_singlet_raw(p, q, S216);
      CHECK(std::fabs(r - std::llround(r)) < 1e-5, "S216 raw m integral");
      r = su3_singlet_raw(p, q, S648);
      CHECK(std::fabs(r - std::llround(r)) < 1e-5, "S648 raw m integral");
      r = su3_singlet_raw(p, q, S1080);
      CHECK(std::fabs(r - std::llround(r)) < 1e-5, "S1080 raw m integral");
    }

  return report("test_singlets");
}
