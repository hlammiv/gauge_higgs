#pragma once
// General arbitrary-irrep engine for SU(N), built by the tensor-power + Young-
// symmetrizer construction. An irrep is specified by a Young diagram (row lengths,
// <= N-1 rows) or equivalently Dynkin labels. The irrep space is the image of the
// Young symmetrizer c_lambda acting on V^{xn} (V = fundamental, n = #boxes); the
// rep generators are the fundamental generators acting as a derivation on the n
// tensor slots, restricted to that subspace. Because we tensor the actual
// generalized-Gell-Mann generators T^a_F, the resulting T^a_R live in the SAME
// basis as the fundamental, so the rep-agnostic HMC force machinery works unchanged.
//
//   - D^(R)(U) = (U^{xn}) restricted to the irrep subspace  (manifest homomorphism)
//   - T^a_R    = (sum_s 1..(T^a_F on slot s)..1) restricted     (-> [T^a_R,T^b_R]=if T^c_R)
//
// Correctness validated in test_general.cpp via the algebra commutator identity,
// dimension/Casimir, and character cross-checks against the explicit fundamental
// and adjoint reps. See theory_notes §6.2, memory: arbitrary-rep-generator-construction.
//
// Cost: O(N^n) intermediate (one-time build) and O(N^n) per rep_matrix call. Fine
// for small irreps; for large irreps switch rep_matrix to the cached exp(i w.T_R)
// path. N^n is capped (kMaxTensorDim) to stay safe.
#include "rep/rep_fundamental.hpp"  // generators<N>, cmat_to_dmat
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <string>

namespace gh {

template <int N>
struct GeneralRep : Representation<N> {
  enum class RealType { Complex, Real };
  static constexpr long kMaxTensorDim = 200000;

  std::vector<int> rows;                 // Young diagram row lengths
  int  n = 0;                            // number of boxes
  long dimT = 0;                         // N^n
  std::vector<long> powN;                // powN[k] = N^k
  std::vector<std::vector<Complex>> W;   // orthonormal basis of the irrep subspace in V^{xn}
  std::vector<Complex> Wflat;            // contiguous [l*dimT + i] copy of W for fast lift/project
  std::string repname;

  // Convert SU(N) Dynkin labels [a_1..a_{N-1}] to Young row lengths.
  static std::vector<int> rows_from_dynkin(const std::vector<int>& dynkin) {
    std::vector<int> r(dynkin.size(), 0);
    for (int i = 0; i < static_cast<int>(dynkin.size()); ++i)
      for (int k = i; k < static_cast<int>(dynkin.size()); ++k) r[i] += dynkin[k];
    while (!r.empty() && r.back() == 0) r.pop_back();
    return r;
  }

  explicit GeneralRep(std::vector<int> young_rows, RealType rt = RealType::Complex) {
    // ---- validate diagram ----
    while (!young_rows.empty() && young_rows.back() == 0) young_rows.pop_back();
    rows = young_rows;
    if (rows.empty()) throw std::runtime_error("GeneralRep: empty diagram (singlet not supported)");
    for (std::size_t i = 0; i + 1 < rows.size(); ++i)
      if (rows[i] < rows[i + 1]) throw std::runtime_error("GeneralRep: rows must be non-increasing");
    if (static_cast<int>(rows.size()) > N - 1)
      throw std::runtime_error("GeneralRep: > N-1 rows (strip full columns of N)");
    n = std::accumulate(rows.begin(), rows.end(), 0);

    powN.resize(n + 1); powN[0] = 1;
    for (int k = 1; k <= n; ++k) powN[k] = powN[k - 1] * N;
    dimT = powN[n];
    if (dimT > kMaxTensorDim)
      throw std::runtime_error("GeneralRep: N^n exceeds cap; irrep too large for tensor build");

    build_basis();
    build_generators();

    this->d = static_cast<int>(W.size());
    Wflat.resize(static_cast<std::size_t>(this->d) * dimT);
    for (int l = 0; l < this->d; ++l)
      std::copy(W[l].begin(), W[l].end(), Wflat.begin() + static_cast<std::size_t>(l) * dimT);
    this->nality = n % N;
    this->real = (rt == RealType::Real);
    repname = "general[";
    for (std::size_t i = 0; i < rows.size(); ++i) repname += (i ? "," : "") + std::to_string(rows[i]);
    repname += "]";
  }

  std::string name() const override { return repname; }

  // D^(R)(U) = <W_k| U^{xn} |W_l>. Materializes the full d x d matrix (used by tests/
  // characters). The hot-path operations below avoid building it.
  DMat rep_matrix(const Cmat<N>& U) const override {
    const int d = static_cast<int>(W.size());
    DMat D(d, d);
    for (int l = 0; l < d; ++l) {
      std::vector<Complex> v = W[l];
      for (int s = 0; s < n; ++s) v = apply_onebody(U, v, s);  // U on each slot -> U^{xn} W_l
      for (int k = 0; k < d; ++k) D(k, l) = inner(W[k], v);
    }
    return D;
  }

  // Fast path: apply D^(R)(U) (or its dagger) to ONE vector via lift -> U^{xn} -> project,
  // never forming the d x d matrix. Cost O(d*N^n + n*N^{n+1}) vs O(d^2*N^n) for the matrix.
  DVec rotate(const Cmat<N>& U, const DVec& phi) const override { return apply_tensor(U, phi, false); }
  DVec rotate_dag(const Cmat<N>& U, const DVec& phi) const override { return apply_tensor(U, phi, true); }
  AlgVec<N> hop_link_g(const Cmat<N>& U, const DVec& phi_x, const DVec& phi_y) const override {
    const DVec Dy = apply_tensor(U, phi_y, false);
    AlgVec<N> g{};
    for (int a = 0; a < n_gen<N>(); ++a) g[a] = -2.0 * dot(phi_x, this->T[a] * Dy).imag();
    return g;
  }

 private:
  // Apply the N x N matrix M to tensor slot s, src -> dst (dst overwritten). In-place
  // variant of apply_onebody (no allocation) for the hot path.
  void apply_onebody_into(const Cmat<N>& M, const Complex* src, Complex* dst, int s) const {
    std::fill(dst, dst + dimT, Complex(0, 0));
    const long blk = powN[n - 1 - s];
    for (long idx = 0; idx < dimT; ++idx) {
      const Complex c = src[idx];
      if (c == Complex(0, 0)) continue;
      const int digit = static_cast<int>((idx / blk) % N);
      const long base = idx - static_cast<long>(digit) * blk;
      for (int j = 0; j < N; ++j) dst[base + static_cast<long>(j) * blk] += M(j, digit) * c;
    }
  }

  DVec apply_tensor(const Cmat<N>& U, const DVec& phi, bool dag) const {
    const int dd = static_cast<int>(W.size());
    // thread-local ping-pong scratch (grown once per thread; no per-call allocation).
    static thread_local std::vector<Complex> sa, sb;
    if (static_cast<long>(sa.size()) < dimT) { sa.resize(dimT); sb.resize(dimT); }
    std::fill(sa.begin(), sa.begin() + dimT, Complex(0, 0));  // lift: sa = sum_l phi_l W_l
    for (int l = 0; l < dd; ++l) {
      const Complex c = phi(l);
      if (c == Complex(0, 0)) continue;
      const Complex* w = &Wflat[static_cast<std::size_t>(l) * dimT];
      for (long i = 0; i < dimT; ++i) sa[i] += c * w[i];
    }
    const Cmat<N> M = dag ? U.dagger() : U;                  // apply U^{xn} (or (U^dag)^{xn})
    Complex* cur = sa.data();
    Complex* oth = sb.data();
    for (int s = 0; s < n; ++s) { apply_onebody_into(M, cur, oth, s); std::swap(cur, oth); }
    DVec out(dd);                                            // project onto basis
    for (int k = 0; k < dd; ++k) {
      const Complex* w = &Wflat[static_cast<std::size_t>(k) * dimT];
      Complex acc(0, 0);
      for (long i = 0; i < dimT; ++i) acc += std::conj(w[i]) * cur[i];
      out(k) = acc;
    }
    return out;
  }

  // tuple <-> index (slot 0 most significant).
  long tuple_to_idx(const std::vector<int>& t) const {
    long idx = 0;
    for (int j = 0; j < n; ++j) idx += static_cast<long>(t[j]) * powN[n - 1 - j];
    return idx;
  }
  std::vector<int> idx_to_tuple(long idx) const {
    std::vector<int> t(n);
    for (int j = 0; j < n; ++j) { t[j] = static_cast<int>((idx / powN[n - 1 - j]) % N); }
    return t;
  }

  static Complex inner(const std::vector<Complex>& a, const std::vector<Complex>& b) {
    Complex s(0, 0);
    for (std::size_t i = 0; i < a.size(); ++i) s += std::conj(a[i]) * b[i];
    return s;
  }

  // Apply the N x N matrix M to tensor slot s of vector v (over V^{xn}).
  std::vector<Complex> apply_onebody(const Cmat<N>& M, const std::vector<Complex>& v, int s) const {
    std::vector<Complex> out(dimT, Complex(0, 0));
    const long blk = powN[n - 1 - s];
    for (long idx = 0; idx < dimT; ++idx) {
      const Complex c = v[idx];
      if (c == Complex(0, 0)) continue;
      const int digit = static_cast<int>((idx / blk) % N);
      const long base = idx - static_cast<long>(digit) * blk;
      for (int j = 0; j < N; ++j) out[base + static_cast<long>(j) * blk] += M(j, digit) * c;
    }
    return out;
  }

  // (sum_s T^a_F on slot s) v  -- the derivation action of generator a.
  std::vector<Complex> derivation(int a, const std::vector<Complex>& v) const {
    const Cmat<N>& Ta = generators<N>().T[a];
    std::vector<Complex> out(dimT, Complex(0, 0));
    for (int s = 0; s < n; ++s) {
      std::vector<Complex> tmp = apply_onebody(Ta, v, s);
      for (long i = 0; i < dimT; ++i) out[i] += tmp[i];
    }
    return out;
  }

  // sigma . t : (sigma t)[sigma[j]] = t[j].
  std::vector<int> apply_perm(const std::vector<int>& sigma, const std::vector<int>& t) const {
    std::vector<int> u(n);
    for (int j = 0; j < n; ++j) u[sigma[j]] = t[j];
    return u;
  }
  static int perm_sign(const std::vector<int>& p) {
    const int m = static_cast<int>(p.size());
    std::vector<char> seen(m, 0);
    int sign = 1;
    for (int i = 0; i < m; ++i) {
      if (seen[i]) continue;
      int len = 0, j = i;
      while (!seen[j]) { seen[j] = 1; j = p[j]; ++len; }
      if (len % 2 == 0) sign = -sign;  // even-length cycle is odd
    }
    return sign;
  }

  // Build a permutation group as the product of symmetric groups on disjoint blocks.
  // Returns full permutations (size n) and their signs.
  void build_group(const std::vector<std::vector<int>>& blocks,
                   std::vector<std::vector<int>>& perms, std::vector<int>& signs) const {
    std::vector<int> id(n); std::iota(id.begin(), id.end(), 0);
    perms = { id };
    for (const auto& blk : blocks) {
      std::vector<int> src = blk; std::sort(src.begin(), src.end());
      // all orderings (dst) of src
      std::vector<std::vector<int>> dsts;
      std::vector<int> cur = src;
      do { dsts.push_back(cur); } while (std::next_permutation(cur.begin(), cur.end()));
      std::vector<std::vector<int>> next;
      next.reserve(perms.size() * dsts.size());
      for (const auto& g : perms)
        for (const auto& dst : dsts) {
          std::vector<int> g2 = g;
          for (std::size_t k = 0; k < src.size(); ++k) g2[src[k]] = dst[k];
          next.push_back(std::move(g2));
        }
      perms.swap(next);
    }
    signs.resize(perms.size());
    for (std::size_t i = 0; i < perms.size(); ++i) signs[i] = perm_sign(perms[i]);
  }

  // Young symmetrizer c_lambda = b_lambda a_lambda applied to a standard tuple.
  std::vector<Complex> young_symmetrize(const std::vector<int>& t,
                                        const std::vector<std::vector<int>>& R,
                                        const std::vector<std::vector<int>>& Cperm,
                                        const std::vector<int>& Csign) const {
    // a_lambda: symmetrize over rows (coeff 1)
    std::vector<Complex> v1(dimT, Complex(0, 0));
    for (const auto& p : R) v1[tuple_to_idx(apply_perm(p, t))] += Complex(1, 0);
    // b_lambda: antisymmetrize over columns (coeff sgn)
    std::vector<Complex> res(dimT, Complex(0, 0));
    for (long idx = 0; idx < dimT; ++idx) {
      const Complex c = v1[idx];
      if (c == Complex(0, 0)) continue;
      const std::vector<int> ta = idx_to_tuple(idx);
      for (std::size_t q = 0; q < Cperm.size(); ++q)
        res[tuple_to_idx(apply_perm(Cperm[q], ta))] += Complex(Csign[q], 0) * c;
    }
    return res;
  }

  void build_basis() {
    // canonical tableau: positions 0..n-1 row by row; record row/column membership.
    std::vector<std::vector<int>> rowblocks(rows.size());
    int ncols = rows.empty() ? 0 : rows[0];
    std::vector<std::vector<int>> colblocks(ncols);
    int pos = 0;
    for (std::size_t i = 0; i < rows.size(); ++i)
      for (int j = 0; j < rows[i]; ++j) {
        rowblocks[i].push_back(pos);
        colblocks[j].push_back(pos);
        ++pos;
      }
    std::vector<std::vector<int>> R, Cperm; std::vector<int> Rsign, Csign;
    build_group(rowblocks, R, Rsign);
    build_group(colblocks, Cperm, Csign);

    // image of c_lambda over all standard tuples, orthonormalized.
    W.clear();
    for (long tc = 0; tc < dimT; ++tc) {
      std::vector<Complex> v = young_symmetrize(idx_to_tuple(tc), R, Cperm, Csign);
      for (const auto& w : W) {
        const Complex ov = inner(w, v);
        for (long i = 0; i < dimT; ++i) v[i] -= ov * w[i];
      }
      Real nrm = std::sqrt(std::max(0.0, inner(v, v).real()));
      if (nrm > 1e-9) {
        const Complex inv(1.0 / nrm, 0);
        for (auto& z : v) z *= inv;
        W.push_back(std::move(v));
      }
    }
  }

  void build_generators() {
    const int d = static_cast<int>(W.size());
    this->T.assign(n_gen<N>(), DMat(d, d));
    for (int a = 0; a < n_gen<N>(); ++a) {
      DMat Ta(d, d);
      for (int l = 0; l < d; ++l) {
        const std::vector<Complex> v = derivation(a, W[l]);
        for (int k = 0; k < d; ++k) Ta(k, l) = inner(W[k], v);
      }
      this->T[a] = Ta;
    }
  }
};

}  // namespace gh
