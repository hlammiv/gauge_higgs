// Closure validation for the crystal-like finite subgroups: each generator set must
// close to the correct order, with every element in SU(N) (unitary, det=+1).
#include "check.hpp"
#include "rep/finite_subgroups.hpp"
#include "group/sun.hpp"

using namespace gh;

template <int N>
static void check_group(const std::vector<Cmat<N>>& gens, int order, const char* name) {
  auto G = close_group<N>(gens);
  char m[64];
  std::snprintf(m, sizeof m, "%s closure order = %d (got %d)", name, order, (int)G.size());
  CHECK((int)G.size() == order, m);
  Real wu = 0, wd = 0;
  for (const auto& g : G) { wu = std::max(wu, unitarity_violation<N>(g)); wd = std::max(wd, std::abs(det<N>(g) - Complex(1,0))); }
  std::snprintf(m, sizeof m, "%s elements in SU(%d) (unit %.1e, det1 %.1e)", name, N, wu, wd);
  CHECK(wu < 1e-9 && wd < 1e-9, m);
}

int main() {
  std::printf("-- SU(2) binary polyhedral --\n");
  check_group<2>(gens_2T(), 24,  "2T");
  check_group<2>(gens_2O(), 48,  "2O");
  check_group<2>(gens_2I(), 120, "2I");
  std::printf("-- SU(3) Sigma series --\n");
  check_group<3>(gens_Sigma108(),  108,  "Sigma(108)");
  check_group<3>(gens_Sigma216(),  216,  "Sigma(216)");
  check_group<3>(gens_Sigma648(),  648,  "Sigma(648)");
  check_group<3>(gens_Sigma1080(), 1080, "Sigma(1080)");
  return report("test_groups");
}
