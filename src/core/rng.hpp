#pragma once
// Counter-based, stateless RNG (splitmix64 mixing). Draws are a pure function of
// a key tuple keyed on the GLOBAL site index, so results are identical regardless
// of MPI/OpenMP/SIMD decomposition (Random123/Philox philosophy). See architecture §3.
#include "core/config.hpp"
#include <cstdint>
#include <cmath>

namespace gh {

inline std::uint64_t splitmix64(std::uint64_t x) {
  x += 0x9E3779B97F4A7C15ULL;
  x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
  x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
  return x ^ (x >> 31);
}

struct Rng {
  std::uint64_t seed = 0x1234567890ABCDEFULL;

  explicit Rng(std::uint64_t s) : seed(s) {}

  // Combine several integer key fields into one 64-bit stream id.
  static std::uint64_t key(std::uint64_t a, std::uint64_t b = 0, std::uint64_t c = 0,
                           std::uint64_t d = 0, std::uint64_t e = 0) {
    std::uint64_t h = splitmix64(a);
    h = splitmix64(h ^ (b + 0x9E3779B97F4A7C15ULL));
    h = splitmix64(h ^ (c + 0x9E3779B97F4A7C15ULL));
    h = splitmix64(h ^ (d + 0x9E3779B97F4A7C15ULL));
    h = splitmix64(h ^ (e + 0x9E3779B97F4A7C15ULL));
    return h;
  }

  // Uniform double in (0,1). `k` is the full stream key.
  double uniform(std::uint64_t k) const {
    std::uint64_t h = splitmix64(seed ^ splitmix64(k));
    // 53-bit mantissa in (0,1): add 0.5 ULP to exclude 0.
    return ((h >> 11) + 0.5) * (1.0 / 9007199254740992.0);
  }

  // Standard normal via Box-Muller. Two uniforms from sub-keys of `k`.
  double gauss(std::uint64_t k) const {
    const double u1 = uniform(splitmix64(k * 2 + 1));
    const double u2 = uniform(splitmix64(k * 2 + 2));
    return std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * kPi * u2);
  }
};

}  // namespace gh
