// Microbenchmark: cost of building D^(j)(U) for SU(2) spin-j by the three methods:
//   (A) tensor materialization  rep_matrix_tensor  (ground truth, O(d * dimT) = O(d*2^n))
//   (B) generic fast path       fast_D = dmat_expi(sum_a w_a T^a_R)  (O(d^3) expm)
//   (C) SU(2) closed-form Wigner-d in the irrep eigenbasis of T^3_R, then conjugated
//       back to the tensor basis.  (proposed cheap correct path)
// Also verifies (B)==(A) and (C)==(A) to machine precision (adversarial correctness).
#include "rep/rep_general.hpp"
#include "group/sun.hpp"
#include "core/linalg.hpp"
#include <cstdio>
#include <chrono>
#include <vector>
#include <complex>
#include <cmath>

using namespace gh;
using clk = std::chrono::steady_clock;

// Build the diagonalizing matrix V of T^3_R (the last su(2) generator a=2 in this basis is
// proportional to J_z up to convention). We just diagonalize T[2] (Hermitian) once.
// Then for SU(2), exp(i w.T_R) in that eigenbasis... but a generic U is not diagonal in T^3.
// Instead, the clean closed form: D^(j)(U) = sum over the SAME homomorphism but evaluated via
// the 2x2 U directly using the symmetric-tensor structure: D^(j)(U)_{kl} on the multiset basis
// IS a polynomial in the 2x2 entries. Because rep_general's single-row basis IS the symmetric
// (multiset) basis, we can build D directly as the induced action of U on degree-n symmetric
// tensors -- which is exactly rep_matrix_tensor but we can do it in O(d^2 * something) by
// recursion. Simplest correct O(d^2) builder: the symmetric power S^n(U) matrix element.
//
// S^n(U)_{(multiset beta),(multiset alpha)} = permanent-like sum. For SU(2) (N=2) each basis
// vector is a monomial x^p y^q (p+q=n). Under U: x->a x + c y, y-> b x + d y with
// U=[[a,b],[c,d]]. Then x^p y^q -> (a x+c y)^p (b x+d y)^q, and D_{(p',q'),(p,q)} is the
// coefficient of the normalized monomial. This is the standard SU(2) Wigner-D via generating
// function -- O(d) monomials, each expansion O(d^2) => O(d^3) but with TINY constant (just
// complex mul/add, no expm, no sqrt-matrix log). With normalization to match the orthonormal
// symmetric basis used in rep_general (norm = sqrt(C(n;mult)) = sqrt(p! q! ... )).

// poly: represent (a x + c y)^p as coefficients over x^{p-k} y^k, k=0..p (binomial).
static std::vector<Complex> binom_pow(Complex a, Complex c, int p) {
  // returns vec length p+1: coeff[k] = C(p,k) a^{p-k} c^{k}  (coefficient of x^{p-k} y^k)
  std::vector<Complex> v(p+1);
  // binomial coefficients
  std::vector<double> bc(p+1,1.0);
  for (int k=1;k<=p;++k) bc[k]=bc[k-1]*(p-k+1)/k;
  Complex ap = std::pow(a, p); // a^{p}, then divide by a per step (a may be 0) -> build iteratively
  // iterative powers to avoid pow on complex with zero base issues
  std::vector<Complex> apow(p+1), cpow(p+1);
  apow[0]=Complex(1,0); for(int k=1;k<=p;++k) apow[k]=apow[k-1]*a;
  cpow[0]=Complex(1,0); for(int k=1;k<=p;++k) cpow[k]=cpow[k-1]*c;
  for (int k=0;k<=p;++k) v[k]=bc[k]*apow[p-k]*cpow[k];
  return v;
}

// Build D^(j) in the orthonormal symmetric (multiset) basis used by GeneralRep single-row.
// Basis index l <-> (p=#of index-0, q=#of index-1) with p+q=n, ordered as rep_general builds:
// build_symmetric_basis enumerates non-decreasing tuples (0..0,..). Tuple with q ones, p zeros
// is the multiset; the l-th basis = that multiset. The ordering: t starts all-zero (q=0 -> p=n),
// increments -> visits multisets with increasing #ones? Let's just match by counting ones in the
// sorted tuple. We replicate the same odometer to get the l-> (#ones) mapping.
static std::vector<int> ones_count_order(int n) {
  std::vector<int> t(n,0), order;
  while (true) {
    int ones=0; for(int x:t) ones+=x; order.push_back(ones);
    int j=n-1; while(j>=0 && t[j]==1) --j; if(j<0) break; ++t[j]; for(int k=j+1;k<n;++k) t[k]=t[j];
  }
  return order; // order[l] = number of ones (=q) for basis vector l
}

int main(int argc, char** argv){
  int twoj = argc>1 ? std::atoi(argv[1]) : 12;  // rows=[2j]; spin-6 -> 12
  long iters = argc>2 ? std::atol(argv[2]) : 200000;
  std::vector<int> rows{twoj};
  GeneralRep<2> rep(rows);
  const int d = rep.d, n = twoj;
  std::printf("# SU(2) spin-%g  d=%d  n=%d dimT=%ld use_fast=%d\n",
              twoj/2.0, d, n, (long)1<<n, (int)rep.use_fast);

  // a generic SU(2) link
  AlgVec<2> w{}; w[0]=0.7; w[1]=-0.4; w[2]=1.1;
  Cmat<2> U = expi<2>(alg_to_mat<2>(w));
  // 2x2 entries U=[[a,b],[c,dd]]
  Complex a=U(0,0), b=U(0,1), c=U(1,0), dd=U(1,1);

  // (A) tensor ground truth
  DMat DA = rep.rep_matrix_tensor(U);
  // (B) generic fast path
  DMat DB = rep.fast_D(U);
  // (C) closed-form symmetric-power builder in the orthonormal multiset basis
  std::vector<int> qof = ones_count_order(n); // q for each column index
  auto factd=[&](int m){ double f=1; for(int i=2;i<=m;++i) f*=i; return f; };
  // norm of multiset basis vector with q ones, p zeros = sqrt( C(n,q) )  [since #distinct perms = n!/(p!q!), and build uses 1/sqrt(count) where count = #distinct perms]
  auto buildC=[&](Complex a,Complex b,Complex c,Complex dd){
    DMat D(d,d);
    // column l: monomial x^{p} y^{q} (p=n-q). image = (a x + c y)^p (b x + d y)^q  [x->col of U: x'=a x + c y? ]
    // Map: fundamental U acts on basis (e0,e1). symmetric power: e0^{p} e1^{q} -> (U e0)^p (U e1)^q
    // U e0 = a e0 + c e1 ; U e1 = b e0 + dd e1.
    for (int l=0;l<d;++l){
      int q=qof[l], p=n-q;
      std::vector<Complex> P0=binom_pow(a,c,p);   // (a e0 + c e1)^p -> coeff over e0^{p-i} e1^{i}
      std::vector<Complex> P1=binom_pow(b,dd,q);   // (b e0 + dd e1)^q
      // product polynomial: coeff of e0^{n-r} e1^{r}
      std::vector<Complex> prod(n+1, Complex(0,0));
      for (int i=0;i<=p;++i) for (int k=0;k<=q;++k) prod[i+k]+=P0[i]*P1[k];
      // prod[r] is coeff in the UNNORMALIZED monomial basis e0^{n-r} e1^{r}.
      // Convert to orthonormal multiset basis: unnormalized monomial e0^{n-r}e1^{r} =
      //   (1/ #perms_r?) ... The multiset basis vector W_r (orthonormal) = (1/sqrt(C(n,r))) * sum_perm.
      // The monomial-coefficient map to orthonormal basis row r' is prod[r'] * sqrt(C(n,r')) / sqrt(C(n,q))?
      // Use: orthonormal D_{r',l} = prod[r'] * sqrt( C(n,r') / C(n,q) )
      double Cl = factd(n)/(factd(q)*factd(p));
      for (int r=0;r<d;++r){
        double Cr = factd(n)/(factd(r)*factd(n-r));
        D(r,l) = prod[r]*std::sqrt(Cl/Cr);
      }
    }
    return D;
  };
  DMat DC = buildC(a,b,c,dd);

  std::printf("# ||DB-DA||=%.3e (generic fast vs tensor)\n", fnorm(DB-DA));
  std::printf("# ||DC-DA||=%.3e (SU2 closed-form vs tensor)\n", fnorm(DC-DA));

  // time each
  auto timeit=[&](const char* nm, auto fn){
    volatile double sink=0;
    auto t0=clk::now();
    for(long i=0;i<iters;++i){ DMat M=fn(); sink+=M.a[0].real(); }
    auto t1=clk::now();
    double s=std::chrono::duration<double>(t1-t0).count();
    std::printf("%-22s %10.4f s  %10.4f us/call  (sink=%.3g)\n", nm, s, s/iters*1e6, (double)sink);
  };
  timeit("tensor", [&]{ return rep.rep_matrix_tensor(U); });
  timeit("generic_fast(dmat_expi)", [&]{ return rep.fast_D(U); });
  timeit("SU2_closed_form", [&]{ return buildC(a,b,c,dd); });
  return 0;
}
