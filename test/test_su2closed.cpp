// Adversarial exactness gate for the SU(2) closed-form D^(j)(U) drop-in:
//  (1) closed-form D == tensor ground-truth D bit-close, on many random links, spin 3/4/6
//  (2) D unitary, homomorphism D(U1 U2)=D(U1)D(U2)
//  (3) GeneralRep combined HMC reversibility on spin-3/4/6 (NOT in test_scalar.cpp)
//  (4) GeneralRep combined <exp(-dH)>~1 + acceptance on spin-3 (NOT in test_scalar.cpp)
//  (5) cached vs per-call rotate/rotate_dag/hop BIT-IDENTICAL (the noinline contract)
#include "check.hpp"
#include "hmc/gauge_higgs_hmc.hpp"
#include "rep/rep_general.hpp"
#include <random>
using namespace gh;

template <int D> static std::array<int,D> cube(int n){std::array<int,D> L{};for(int mu=0;mu<D;++mu)L[mu]=n;return L;}

static Cmat<2> rand_su2(std::mt19937_64& r){
  std::normal_distribution<Real> g(0.0,0.8); AlgVec<2> v{};
  for(int a=0;a<3;++a) v[a]=g(r); return expi<2>(alg_to_mat<2>(v));
}

static void test_matrix_equal(int twoj){
  GeneralRep<2> rep({twoj});
  char m[96]; std::snprintf(m,sizeof m,"spin-%g closed-form ENABLED (su2_closed)",twoj/2.0);
  CHECK(rep.su2_closed && rep.use_fast, m);
  std::mt19937_64 rng(12345+twoj);
  Real worstD=0, worstU=0, worstHom=0;
  for(int t=0;t<200;++t){
    Cmat<2> U1=rand_su2(rng), U2=rand_su2(rng);
    DMat Dc=rep.fast_D(U1), Dt=rep.rep_matrix_tensor(U1);
    worstD=std::max(worstD, fnorm(Dc-Dt));
    DMat UU=Dc*Dc.dagger();
    for(int i=0;i<rep.d;++i)for(int j=0;j<rep.d;++j)
      worstU=std::max(worstU,std::abs(UU(i,j)-Complex(i==j?1.0:0.0,0.0)));
    DMat lhs=rep.fast_D(U1*U2), rhs=Dc*rep.fast_D(U2);
    worstHom=std::max(worstHom, fnorm(lhs-rhs));
  }
  std::snprintf(m,sizeof m,"spin-%g closed-form vs tensor (200 links)",twoj/2.0);
  CHECK_CLOSE(worstD,0.0,1e-9,m);
  std::snprintf(m,sizeof m,"spin-%g closed-form unitary",twoj/2.0);
  CHECK_CLOSE(worstU,0.0,1e-9,m);
  std::snprintf(m,sizeof m,"spin-%g closed-form homomorphism",twoj/2.0);
  CHECK_CLOSE(worstHom,0.0,1e-9,m);
}

static void test_cached_bitidentical(int twoj){
  GeneralRep<2> rep({twoj}); std::mt19937_64 rng(777+twoj);
  std::normal_distribution<Real> g(0.0,0.7);
  Real worst=0;
  for(int t=0;t<50;++t){
    Cmat<2> U=rand_su2(rng);
    DVec px(rep.d), py(rep.d);
    for(int k=0;k<rep.d;++k){ px(k)=Complex(g(rng),g(rng)); py(k)=Complex(g(rng),g(rng)); }
    DMat Dc=rep.cache_rep_matrix(U);
    DVec a1=rep.rotate(U,py), a2=rep.rotate_cached(Dc,py);
    DVec b1=rep.rotate_dag(U,py), b2=rep.rotate_dag_cached(Dc,py);
    AlgVec<2> h1=rep.hop_link_g(U,px,py), h2=rep.hop_link_g_cached(Dc,px,py);
    for(int k=0;k<rep.d;++k){ worst=std::max(worst,std::abs(a1(k)-a2(k))); worst=std::max(worst,std::abs(b1(k)-b2(k))); }
    for(int a=0;a<3;++a) worst=std::max(worst,std::fabs(h1[a]-h2[a]));
  }
  char m[96]; std::snprintf(m,sizeof m,"spin-%g cached==per-call BIT-IDENTICAL",twoj/2.0);
  CHECK(worst==0.0, m);
}

template <int Dd> static void test_rev(int twoj,std::uint64_t seed){
  GeneralRep<2> rep({twoj});
  GaugeHiggsHMC<Dd,2> hmc(cube<Dd>(4),rep,seed);
  hmc.beta=2.3; hmc.kappa=0.3; hmc.lambda=0.5; hmc.tau=1.0; hmc.nmd=12; hmc.reunit_each_traj=false;
  hmc.U.hot(hmc.rng,0.5); hmc.phi.gaussian(hmc.rng,1,rep.real,0.6); hmc.refresh_momenta();
  std::vector<Cmat<2>> U0=hmc.U.u; std::vector<Complex> phi0=hmc.phi.data;
  hmc.md_evolve();
  for(auto& v:hmc.P.p) for(int a=0;a<3;++a) v[a]=-v[a];
  for(auto& z:hmc.pi.data) z=-z;
  hmc.md_evolve();
  Real worst=0;
  for(std::size_t i=0;i<U0.size();++i) worst=std::max(worst,(hmc.U.u[i]-U0[i]).fnorm());
  for(std::size_t i=0;i<phi0.size();++i) worst=std::max(worst,std::abs(hmc.phi.data[i]-phi0[i]));
  char m[96]; std::snprintf(m,sizeof m,"spin-%g combined reversibility",twoj/2.0);
  CHECK_CLOSE(worst,0.0,1e-9,m);
}

template <int Dd> static void test_expdH(int twoj,std::uint64_t seed){
  GeneralRep<2> rep({twoj});
  GaugeHiggsHMC<Dd,2> hmc(cube<Dd>(4),rep,seed);
  hmc.beta=2.3; hmc.kappa=0.2; hmc.lambda=0.5; hmc.tau=1.0; hmc.nmd=20;
  hmc.U.hot(hmc.rng,0.3); hmc.phi.cold(0.5);
  const int ntraj=300; double s=0;
  for(int t=0;t<ntraj;++t){ hmc.trajectory(); s+=std::exp(-hmc.last_dH); }
  double mean=s/ntraj;
  char m[128]; std::snprintf(m,sizeof m,"spin-%g <exp(-dH)>~1 got %.4f acc %.2f",twoj/2.0,mean,hmc.acceptance());
  CHECK(std::fabs(mean-1.0)<0.15,m);
}

int main(){
  std::printf("-- closed-form matrix == tensor --\n");
  test_matrix_equal(6); test_matrix_equal(8); test_matrix_equal(12);
  std::printf("-- cached == per-call bit-identical --\n");
  test_cached_bitidentical(6); test_cached_bitidentical(8); test_cached_bitidentical(12);
  std::printf("-- GeneralRep combined reversibility (spin 3/4/6) --\n");
  test_rev<3>(6,401); test_rev<3>(8,402); test_rev<3>(12,403);
  std::printf("-- GeneralRep <exp(-dH)> (spin 3) --\n");
  test_expdH<3>(6,501);
  return report("test_su2closed");
}
