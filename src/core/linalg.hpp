#pragma once
// Small fixed-size complex linear algebra for the hot path (gauge links, su(N)
// generators), plus a lightweight dynamic complex matrix/vector for the
// arbitrary-irrep representation engine.
#include "core/config.hpp"
#include <array>
#include <vector>
#include <cmath>
#include <cassert>

namespace gh {

// ----------------------------------------------------------------------------
// Fixed-size N x N complex matrix, row-major (a[i*N+j]). Compile-time N.
// ----------------------------------------------------------------------------
template <int N>
struct Cmat {
  std::array<Complex, static_cast<std::size_t>(N) * N> a{};  // zero-initialized

  static constexpr int dim() { return N; }

  Complex&       operator()(int i, int j)       { return a[i * N + j]; }
  const Complex& operator()(int i, int j) const { return a[i * N + j]; }

  static Cmat zero() { return Cmat{}; }
  static Cmat identity() {
    Cmat m{};
    for (int i = 0; i < N; ++i) m(i, i) = Complex(1.0, 0.0);
    return m;
  }

  Cmat& operator+=(const Cmat& o) { for (std::size_t i = 0; i < a.size(); ++i) a[i] += o.a[i]; return *this; }
  Cmat& operator-=(const Cmat& o) { for (std::size_t i = 0; i < a.size(); ++i) a[i] -= o.a[i]; return *this; }
  Cmat& operator*=(Complex s)     { for (auto& z : a) z *= s; return *this; }

  Cmat dagger() const {
    Cmat r{};
    for (int i = 0; i < N; ++i)
      for (int j = 0; j < N; ++j) r(i, j) = std::conj((*this)(j, i));
    return r;
  }

  Complex trace() const {
    Complex t(0.0, 0.0);
    for (int i = 0; i < N; ++i) t += (*this)(i, i);
    return t;
  }

  // Hermitian / anti-Hermitian parts.
  Cmat herm()     const { Cmat d = dagger(); Cmat r = *this; r += d; r *= Complex(0.5, 0.0); return r; }
  Cmat antiherm() const { Cmat d = dagger(); Cmat r = *this; r -= d; r *= Complex(0.5, 0.0); return r; }

  // Frobenius norm.
  Real fnorm() const {
    Real s = 0.0;
    for (const auto& z : a) s += std::norm(z);
    return std::sqrt(s);
  }
};

template <int N>
Cmat<N> operator+(Cmat<N> A, const Cmat<N>& B) { A += B; return A; }
template <int N>
Cmat<N> operator-(Cmat<N> A, const Cmat<N>& B) { A -= B; return A; }
template <int N>
Cmat<N> operator*(const Cmat<N>& A, const Cmat<N>& B) {
  Cmat<N> C{};
  for (int i = 0; i < N; ++i)
    for (int k = 0; k < N; ++k) {
      const Complex aik = A(i, k);
      if (aik == Complex(0.0, 0.0)) continue;
      for (int j = 0; j < N; ++j) C(i, j) += aik * B(k, j);
    }
  return C;
}
template <int N>
Cmat<N> operator*(Cmat<N> A, Complex s) { A *= s; return A; }
template <int N>
Cmat<N> operator*(Complex s, Cmat<N> A) { A *= s; return A; }

// Real part of the trace of A*B without forming the product.
template <int N>
Real reTrProd(const Cmat<N>& A, const Cmat<N>& B) {
  Real s = 0.0;
  for (int i = 0; i < N; ++i)
    for (int k = 0; k < N; ++k) s += (A(i, k) * B(k, i)).real();
  return s;
}
// Full trace of A*B.
template <int N>
Complex trProd(const Cmat<N>& A, const Cmat<N>& B) {
  Complex s(0.0, 0.0);
  for (int i = 0; i < N; ++i)
    for (int k = 0; k < N; ++k) s += A(i, k) * B(k, i);
  return s;
}

// Determinant via Gaussian elimination with partial pivoting (small N).
template <int N>
Complex det(const Cmat<N>& M) {
  std::array<Complex, static_cast<std::size_t>(N) * N> m = M.a;
  auto at = [&](int i, int j) -> Complex& { return m[i * N + j]; };
  Complex d(1.0, 0.0);
  for (int col = 0; col < N; ++col) {
    int piv = col;
    Real best = std::abs(at(col, col));
    for (int r = col + 1; r < N; ++r) {
      Real v = std::abs(at(r, col));
      if (v > best) { best = v; piv = r; }
    }
    if (best == 0.0) return Complex(0.0, 0.0);
    if (piv != col) {
      for (int j = 0; j < N; ++j) std::swap(at(piv, j), at(col, j));
      d = -d;
    }
    d *= at(col, col);
    const Complex inv = Complex(1.0, 0.0) / at(col, col);
    for (int r = col + 1; r < N; ++r) {
      const Complex f = at(r, col) * inv;
      if (f == Complex(0.0, 0.0)) continue;
      for (int j = col; j < N; ++j) at(r, j) -= f * at(col, j);
    }
  }
  return d;
}

// ----------------------------------------------------------------------------
// Fixed-size complex N-vector (for fundamental-rep scalars / quick math).
// ----------------------------------------------------------------------------
template <int N>
struct Cvec {
  std::array<Complex, N> v{};
  Complex&       operator()(int i)       { return v[i]; }
  const Complex& operator()(int i) const { return v[i]; }
  static Cvec zero() { return Cvec{}; }
  Cvec& operator+=(const Cvec& o) { for (int i = 0; i < N; ++i) v[i] += o.v[i]; return *this; }
  Cvec& operator*=(Complex s)     { for (int i = 0; i < N; ++i) v[i] *= s; return *this; }
  Real norm2() const { Real s = 0; for (auto& z : v) s += std::norm(z); return s; }
};
template <int N>
Complex dot(const Cvec<N>& x, const Cvec<N>& y) {  // x† y
  Complex s(0.0, 0.0);
  for (int i = 0; i < N; ++i) s += std::conj(x.v[i]) * y.v[i];
  return s;
}
template <int N>
Cvec<N> operator*(const Cmat<N>& M, const Cvec<N>& x) {
  Cvec<N> y{};
  for (int i = 0; i < N; ++i)
    for (int j = 0; j < N; ++j) y.v[i] += M(i, j) * x.v[j];
  return y;
}

// ----------------------------------------------------------------------------
// Dynamic complex matrix / vector for the arbitrary-irrep representation engine
// (representation dimension d(R) is only known at runtime). Not hot-path; kept
// simple and correct. Row-major.
// ----------------------------------------------------------------------------
struct DMat {
  int rows = 0, cols = 0;
  std::vector<Complex> a;
  DMat() = default;
  DMat(int r, int c) : rows(r), cols(c), a(static_cast<std::size_t>(r) * c, Complex(0, 0)) {}
  Complex&       operator()(int i, int j)       { return a[static_cast<std::size_t>(i) * cols + j]; }
  const Complex& operator()(int i, int j) const { return a[static_cast<std::size_t>(i) * cols + j]; }
  static DMat identity(int n) { DMat m(n, n); for (int i = 0; i < n; ++i) m(i, i) = Complex(1, 0); return m; }
  DMat dagger() const {
    DMat r(cols, rows);
    for (int i = 0; i < rows; ++i)
      for (int j = 0; j < cols; ++j) r(j, i) = std::conj((*this)(i, j));
    return r;
  }
  Complex trace() const { Complex t(0, 0); for (int i = 0; i < std::min(rows, cols); ++i) t += (*this)(i, i); return t; }
};
inline DMat operator*(const DMat& A, const DMat& B) {
  assert(A.cols == B.rows);
  DMat C(A.rows, B.cols);
  for (int i = 0; i < A.rows; ++i)
    for (int k = 0; k < A.cols; ++k) {
      const Complex aik = A(i, k);
      if (aik == Complex(0, 0)) continue;
      for (int j = 0; j < B.cols; ++j) C(i, j) += aik * B(k, j);
    }
  return C;
}
inline DMat operator+(const DMat& A, const DMat& B) {
  DMat C(A.rows, A.cols);
  for (std::size_t i = 0; i < A.a.size(); ++i) C.a[i] = A.a[i] + B.a[i];
  return C;
}
inline DMat operator-(const DMat& A, const DMat& B) {
  DMat C(A.rows, A.cols);
  for (std::size_t i = 0; i < A.a.size(); ++i) C.a[i] = A.a[i] - B.a[i];
  return C;
}
inline DMat operator*(DMat A, Complex s) { for (auto& z : A.a) z *= s; return A; }
inline Real fnorm(const DMat& A) { Real s = 0; for (const auto& z : A.a) s += std::norm(z); return std::sqrt(s); }

// Re Tr(A^dag B) — real inner product on the space of matrices.
inline Real reTrDagProd(const DMat& A, const DMat& B) {
  Real s = 0.0;
  for (std::size_t i = 0; i < A.a.size(); ++i) s += (std::conj(A.a[i]) * B.a[i]).real();
  return s;
}

// Eigenvalues of a real-symmetric n x n matrix (row-major, modified in place) via
// cyclic Jacobi. Returns the n eigenvalues (unsorted). Small n only.
inline std::vector<Real> eig_sym(std::vector<Real>& A, int n) {
  std::vector<Real> off(0);
  auto at = [&](int i, int j) -> Real& { return A[i * n + j]; };
  for (int sweep = 0; sweep < 100; ++sweep) {
    Real offsum = 0.0;
    for (int p = 0; p < n; ++p) for (int q = p + 1; q < n; ++q) offsum += at(p, q) * at(p, q);
    if (offsum < 1e-28) break;
    for (int p = 0; p < n; ++p)
      for (int q = p + 1; q < n; ++q) {
        if (std::fabs(at(p, q)) < 1e-300) continue;
        const Real theta = (at(q, q) - at(p, p)) / (2.0 * at(p, q));
        const Real t = (theta >= 0 ? 1.0 : -1.0) / (std::fabs(theta) + std::sqrt(theta * theta + 1.0));
        const Real c = 1.0 / std::sqrt(t * t + 1.0), s = t * c;
        for (int k = 0; k < n; ++k) {
          const Real akp = at(k, p), akq = at(k, q);
          at(k, p) = c * akp - s * akq; at(k, q) = s * akp + c * akq;
        }
        for (int k = 0; k < n; ++k) {
          const Real apk = at(p, k), aqk = at(q, k);
          at(p, k) = c * apk - s * aqk; at(q, k) = s * apk + c * aqk;
        }
      }
  }
  std::vector<Real> ev(n);
  for (int i = 0; i < n; ++i) ev[i] = A[i * n + i];
  return ev;
}

struct DVec {
  std::vector<Complex> v;
  DVec() = default;
  explicit DVec(int n) : v(n, Complex(0, 0)) {}
  int        size() const { return static_cast<int>(v.size()); }
  Complex&       operator()(int i)       { return v[i]; }
  const Complex& operator()(int i) const { return v[i]; }
  Real norm2() const { Real s = 0; for (auto& z : v) s += std::norm(z); return s; }
};
inline Complex dot(const DVec& x, const DVec& y) {  // x† y
  Complex s(0, 0);
  for (int i = 0; i < x.size(); ++i) s += std::conj(x.v[i]) * y.v[i];
  return s;
}
inline DVec operator*(const DMat& M, const DVec& x) {
  DVec y(M.rows);
  for (int i = 0; i < M.rows; ++i)
    for (int j = 0; j < M.cols; ++j) y.v[i] += M(i, j) * x.v[j];
  return y;
}

}  // namespace gh
