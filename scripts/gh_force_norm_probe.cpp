#include "hmc/gauge_higgs_hmc.hpp"
#include "action/scalar_invariants.hpp"
#include "rep/rep_general.hpp"
#include <cstdio>
#include <cmath>
using namespace gh;
static double norm_link(const LinkMom<4,2>& F){ double s=0; for(auto&v:F.p) for(int a=0;a<n_gen<2>();++a) s+=v[a]*v[a]; return std::sqrt(s);} 
static double norm_phi(const CVecField<4>& F){ double s=0; for(auto&z:F.data) s+=std::norm(z); return std::sqrt(s);} 
int main(int argc,char**argv){
  double kappa=atof(argv[1]); int L=atoi(argv[2]);
  std::vector<int> rows{6}; GeneralRep<2> rep(rows);
  CasimirChannels<2> ch(rep);
  std::vector<Real> f{0.1287,0.1548,0.1835,0.2399,0.0056,0.1745,0.1130};
  MultiInvariantPotential<2> pot(ch,f,0.113);
  std::array<int,4> ext{L,L,L,L};
  GaugeHiggsHMC<4,2> hmc(ext,rep,7);
  hmc.beta=1.5; hmc.kappa=kappa; hmc.tau=1.0; hmc.nmd=24; hmc.potential=&pot; hmc.frozen_phi=true;
  hmc.U.hot(hmc.rng,0.8); hmc.phi.gaussian(hmc.rng,12345,rep.real,0.3); hmc.normalize_phi();
  for(int t=0;t<40;++t) hmc.trajectory();
  // measure force norms at the current config
  LinkMom<4,2> Fg(hmc.lat), Fm(hmc.lat); CVecField<4> Fs(hmc.lat,rep.d);
  Fg.zero(); add_gauge_force<4,2>(hmc.U,hmc.beta,Fg);
  Fm.zero(); add_matter_link_force<4,2>(hmc.phi,hmc.U,rep,hmc.kappa,Fm,nullptr);
  scalar_force<4,2>(hmc.phi,hmc.U,rep,hmc.kappa,pot,Fs,nullptr);
  // frozen: project scalar force tangent for fair comparison? report raw + note
  printf("kappa=%.2f L=%d  |F_gauge|=%.4f  |F_matterlink|=%.4f  |F_scalar|=%.4f  ratio_scalar/gauge=%.2f  ratio_matterlink/gauge=%.2f\n",
         kappa,L,norm_link(Fg),norm_link(Fm),norm_phi(Fs),
         norm_phi(Fs)/norm_link(Fg), norm_link(Fm)/norm_link(Fg));
  return 0;
}
