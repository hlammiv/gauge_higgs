// Production-relevant microbenchmark: per-call cost of the THREE hot-path rep ops as
// the HMC actually calls them, comparing the CURRENT active path vs the proposed
// SU(2) closed-form D path (build D once + cached matvec/hop).
//
//  CURRENT (spin-3/4, use_fast=0): matrix-free apply_tensor in rotate/rotate_dag, and
//    apply_tensor inside hop_link_g. No d x d matrix is ever built.
//  CURRENT (spin-6, use_fast=1): fast_D = dmat_expi(...) built fresh inside rotate /
//    rotate_dag / hop_link_g (or once per link in the Dcache via cache_rep_matrix).
//  PROPOSED: build closed-form D^(j)(U) once per link, then rotate_cached / hop_cached.
//
// We measure a realistic per-link kick footprint: the scalar force does ONE rotate (fwd)
// + ONE rotate_dag (bwd), and the matter force does ONE hop_link_g, per link per kick.
// With the Dcache, D is built ONCE per link and reused by all three. So per link per kick:
//   current use_fast: 1 D build (cache) + 1 matvec + 1 dagger-matvec + 1 hop  (cache path)
//   current tensor  : 3 x apply_tensor (no cache; rotate + rotate_dag + hop each lift/apply/project)
//   proposed        : 1 closed-form D build + 1 matvec + 1 dagger-matvec + 1 hop
#include "rep/rep_general.hpp"
#include "group/sun.hpp"
#include "core/linalg.hpp"
#include <cstdio>
#include <chrono>
#include <vector>

using namespace gh;
using clk = std::chrono::steady_clock;

// closed-form S^n(U) builder (validated bit-close in bench_repD_su2.cpp), copied here.
static std::vector<Complex> binom_pow(Complex a, Complex c, int p) {
  std::vector<double> bc(p+1,1.0);
  for (int k=1;k<=p;++k) bc[k]=bc[k-1]*(p-k+1)/k;
  std::vector<Complex> apow(p+1), cpow(p+1), v(p+1);
  apow[0]=Complex(1,0); for(int k=1;k<=p;++k) apow[k]=apow[k-1]*a;
  cpow[0]=Complex(1,0); for(int k=1;k<=p;++k) cpow[k]=cpow[k-1]*c;
  for (int k=0;k<=p;++k) v[k]=bc[k]*apow[p-k]*cpow[k];
  return v;
}
static std::vector<int> ones_count_order(int n) {
  std::vector<int> t(n,0), order;
  while (true) {
    int ones=0; for(int x:t) ones+=x; order.push_back(ones);
    int j=n-1; while(j>=0 && t[j]==1) --j; if(j<0) break; ++t[j]; for(int k=j+1;k<n;++k) t[k]=t[j];
  }
  return order;
}
static DMat closed_form_D(const Cmat<2>& U, int n, int d, const std::vector<int>& qof) {
  Complex a=U(0,0), b=U(0,1), c=U(1,0), dd=U(1,1);
  auto factd=[&](int m){ double f=1; for(int i=2;i<=m;++i) f*=i; return f; };
  DMat D(d,d);
  for (int l=0;l<d;++l){
    int q=qof[l], p=n-q;
    std::vector<Complex> P0=binom_pow(a,c,p);
    std::vector<Complex> P1=binom_pow(b,dd,q);
    std::vector<Complex> prod(n+1, Complex(0,0));
    for (int i=0;i<=p;++i) for (int k=0;k<=q;++k) prod[i+k]+=P0[i]*P1[k];
    double Cl = factd(n)/(factd(q)*factd(p));
    for (int r=0;r<d;++r){
      double Cr = factd(n)/(factd(r)*factd(n-r));
      D(r,l) = prod[r]*std::sqrt(Cl/Cr);
    }
  }
  return D;
}

int main(int argc, char** argv){
  int twoj = argc>1 ? std::atoi(argv[1]) : 6;
  long iters = argc>2 ? std::atol(argv[2]) : 200000;
  std::vector<int> rows{twoj};
  GeneralRep<2> rep(rows);
  const int d = rep.d, n = twoj;
  std::printf("# SU(2) spin-%g d=%d n=%d use_fast=%d\n", twoj/2.0, d, n, (int)rep.use_fast);

  AlgVec<2> w{}; w[0]=0.7; w[1]=-0.4; w[2]=1.1;
  Cmat<2> U = expi<2>(alg_to_mat<2>(w));
  std::vector<int> qof = ones_count_order(n);

  // random-ish phi vectors
  DVec px(d), py(d);
  for (int k=0;k<d;++k){ px(k)=Complex(0.3*std::sin(1.0+k), -0.2*std::cos(0.5+k));
                          py(k)=Complex(-0.1*std::cos(2.0+k), 0.4*std::sin(0.7+k)); }

  // correctness: closed-form cached ops == tensor per-call ops
  DMat Dcf = closed_form_D(U, n, d, qof);
  DVec r_cf  = rep.rotate_cached(Dcf, py);
  DVec r_ten = rep.rotate(U, py);
  DVec rd_cf  = rep.rotate_dag_cached(Dcf, py);
  DVec rd_ten = rep.rotate_dag(U, py);
  AlgVec<2> h_cf  = rep.hop_link_g_cached(Dcf, px, py);
  AlgVec<2> h_ten = rep.hop_link_g(U, px, py);
  Real er=0,erd=0,eh=0;
  for(int k=0;k<d;++k){ er=std::max(er,std::abs(r_cf(k)-r_ten(k))); erd=std::max(erd,std::abs(rd_cf(k)-rd_ten(k))); }
  for(int a=0;a<3;++a) eh=std::max(eh,std::fabs(h_cf[a]-h_ten[a]));
  std::printf("# closed-form vs current rotate=%.2e rotate_dag=%.2e hop=%.2e\n", er, erd, eh);

  auto timeit=[&](const char* nm, auto fn){
    volatile double sink=0;
    auto t0=clk::now();
    for(long i=0;i<iters;++i) sink+=fn();
    auto t1=clk::now();
    double s=std::chrono::duration<double>(t1-t0).count();
    std::printf("%-34s %9.4f us/link-kick  (sink=%.3g)\n", nm, s/iters*1e6, (double)sink);
  };

  // CURRENT per-link-per-kick footprint (no Dcache): 1 rotate + 1 rotate_dag + 1 hop
  timeit("current  (rotate+rd+hop)", [&]{
    DVec a = rep.rotate(U, py);
    DVec bb = rep.rotate_dag(U, py);
    AlgVec<2> g = rep.hop_link_g(U, px, py);
    return a(0).real()+bb(0).real()+g[0];
  });
  // CURRENT with link cache (use_fast reps): 1 cache_rep_matrix + cached ops
  timeit("current+cache (Dexpi + 3 cached)", [&]{
    DMat D = rep.cache_rep_matrix(U);   // fast_D=dmat_expi for use_fast, else rep_matrix(U)
    DVec a = rep.rotate_cached(D, py);
    DVec bb = rep.rotate_dag_cached(D, py);
    AlgVec<2> g = rep.hop_link_g_cached(D, px, py);
    return a(0).real()+bb(0).real()+g[0];
  });
  // PROPOSED: 1 closed-form D build + 3 cached ops
  timeit("proposed (Dclosed + 3 cached)", [&]{
    DMat D = closed_form_D(U, n, d, qof);
    DVec a = rep.rotate_cached(D, py);
    DVec bb = rep.rotate_dag_cached(D, py);
    AlgVec<2> g = rep.hop_link_g_cached(D, px, py);
    return a(0).real()+bb(0).real()+g[0];
  });
  return 0;
}
