#pragma once
// Arbitrary-dimension lattice geometry: site<->coordinate maps, forward/backward
// neighbor tables, even-odd parity. Compile-time dimension D, runtime extents.
#include "core/config.hpp"
#include <array>
#include <vector>
#include <cstdint>

namespace gh {

template <int D>
struct Lattice {
  std::array<int, D>      L{};       // extents per direction
  std::int64_t            vol = 0;   // total number of sites
  std::array<std::int64_t, D> stride{};  // mixed-radix strides (dir 0 fastest)
  std::vector<std::int64_t> fwd;     // [site*D + mu] -> neighbor in +mu
  std::vector<std::int64_t> bwd;     // [site*D + mu] -> neighbor in -mu
  std::vector<int>          parity;  // [site] -> (sum coords) & 1

  explicit Lattice(const std::array<int, D>& extents) : L(extents) {
    vol = 1;
    for (int mu = 0; mu < D; ++mu) { stride[mu] = vol; vol *= L[mu]; }
    fwd.resize(static_cast<std::size_t>(vol) * D);
    bwd.resize(static_cast<std::size_t>(vol) * D);
    parity.resize(vol);
    std::array<int, D> x{};
    for (std::int64_t s = 0; s < vol; ++s) {
      coords(s, x);
      int par = 0;
      for (int mu = 0; mu < D; ++mu) par += x[mu];
      parity[s] = par & 1;
      for (int mu = 0; mu < D; ++mu) {
        std::array<int, D> y = x;
        y[mu] = (x[mu] + 1) % L[mu];      fwd[s * D + mu] = index(y);
        y[mu] = (x[mu] - 1 + L[mu]) % L[mu]; bwd[s * D + mu] = index(y);
      }
    }
  }

  void coords(std::int64_t s, std::array<int, D>& x) const {
    for (int mu = 0; mu < D; ++mu) { x[mu] = static_cast<int>(s % L[mu]); s /= L[mu]; }
  }
  std::int64_t index(const std::array<int, D>& x) const {
    std::int64_t s = 0;
    for (int mu = 0; mu < D; ++mu) s += static_cast<std::int64_t>(x[mu]) * stride[mu];
    return s;
  }
  std::int64_t neighbor_fwd(std::int64_t s, int mu) const { return fwd[s * D + mu]; }
  std::int64_t neighbor_bwd(std::int64_t s, int mu) const { return bwd[s * D + mu]; }

  // Number of plaquettes (oriented mu<nu) = vol * D(D-1)/2.
  std::int64_t n_plaq() const { return vol * (static_cast<std::int64_t>(D) * (D - 1) / 2); }
};

}  // namespace gh
