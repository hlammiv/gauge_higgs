#pragma once
// Compile-time configuration. NDIM (spacetime dimension) and NCOL (number of
// colors N) are compile-time constants in the HiRep/HILA tradition: set them via
// -DNDIM=.. -DNCOL=.. at build time so the production binary is monomorphic and
// fast. The test suite instantiates the templates at several (D,N) directly.
#include <complex>

namespace gh {

using Real    = double;
using Complex = std::complex<Real>;

#ifndef NDIM
#define NDIM 4
#endif
#ifndef NCOL
#define NCOL 2
#endif

constexpr int  kDim = NDIM;   // default spacetime dimension for the driver
constexpr int  kN   = NCOL;   // default number of colors for the driver
constexpr Real kPi  = 3.14159265358979323846264338327950288;

// Omelyan 2MN minimum-norm integrator parameter (Omelyan-Mryglod-Folk).
constexpr Real kOmelyanLambda = 0.1931833275037836;

}  // namespace gh
