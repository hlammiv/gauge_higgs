#pragma once
// Arbitrary-irrep representation interface (the HiRep principle): the HMC force
// machinery is rep-agnostic; a Representation only supplies the rep-R generators
// T^a_R (d x d) and the map D^(R)(U) from a fundamental link to its rep matrix.
// See theory_notes §6, memory: hirep-rep-agnostic-force-principle.
#include "group/algebra.hpp"
#include <string>
#include <vector>

namespace gh {

template <int N>
struct Representation {
  int  d = 0;        // representation dimension d(R)
  int  nality = 0;   // N-ality k(R): D^(R)(z 1) = z^k 1, z = exp(2pi i/N)
  bool real = false; // true => scalar carries d real dof (e.g. adjoint)
  std::vector<DMat> T;  // T^a_R, a = 0..N^2-2, each d x d (Hermitian)

  virtual ~Representation() = default;
  virtual std::string name() const = 0;
  // D^(R)(U): the fundamental link U promoted to rep R (d x d).
  virtual DMat rep_matrix(const Cmat<N>& U) const = 0;

  // Hot-path operations: apply D^(R)(U) (or its dagger) to a vector, and assemble
  // the matter-staple link force. The default materializes rep_matrix(U); reps with
  // a large/expensive rep matrix (GeneralRep) override these to apply directly to the
  // vector without ever building the d x d matrix.
  virtual DVec rotate(const Cmat<N>& U, const DVec& phi) const { return rep_matrix(U) * phi; }
  virtual DVec rotate_dag(const Cmat<N>& U, const DVec& phi) const {
    DMat M = rep_matrix(U);
    return M.dagger() * phi;
  }

  // For the matter-staple link force. Returns g^a = 2 Re[ phi_x^dag (i T^a_R) D^(R)(U) phi_y ].
  // Then dS_H/d(omega^a) = -kappa * g^a (see scalar_higgs.hpp).
  virtual AlgVec<N> hop_link_g(const Cmat<N>& U, const DVec& phi_x, const DVec& phi_y) const {
    const DMat DR = rep_matrix(U);
    const DVec Dy = DR * phi_y;                 // D^(R)(U) phi_y
    AlgVec<N> g{};
    for (int a = 0; a < n_gen<N>(); ++a) {
      const DVec t = T[a] * Dy;                 // T^a_R D^(R)(U) phi_y
      const Complex z = dot(phi_x, t);          // phi_x^dag T^a_R D^(R)(U) phi_y
      g[a] = -2.0 * z.imag();                   // 2 Re[ phi_x^dag (i T^a_R) ... ] = -2 Im(z)
    }
    return g;
  }

  // Runtime self-consistency: Tr(T^a_R T^b_R) = T(R) delta^{ab}, sum_a T^a_R T^a_R = C2(R) I.
  // Returns max deviation from those (given the expected Dynkin index and Casimir).
  Real invariant_violation(Real dynkin_T, Real casimir_C2) const {
    Real worst = 0.0;
    const int n = n_gen<N>();
    // Dynkin: Tr(T^a T^b)
    for (int a = 0; a < n; ++a)
      for (int b = 0; b < n; ++b) {
        DMat P = T[a] * T[b];
        Complex tr = P.trace();
        worst = std::max(worst, std::fabs(tr.real() - (a == b ? dynkin_T : 0.0)));
        worst = std::max(worst, std::fabs(tr.imag()));
      }
    // Casimir: sum_a (T^a)^2 = C2 I
    DMat C(d, d);
    for (int a = 0; a < n; ++a) C = C + (T[a] * T[a]);
    for (int i = 0; i < d; ++i)
      for (int j = 0; j < d; ++j)
        worst = std::max(worst, std::abs(C(i, j) - Complex(i == j ? casimir_C2 : 0.0, 0.0)));
    return worst;
  }
};

}  // namespace gh
