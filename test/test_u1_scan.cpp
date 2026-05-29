// Tests for the (beta,kappa) scan extensive observables (scan_obs.hpp).
//
// The CONVENTION test is the critical deliverable: it guarantees the A,B written by the
// driver are EXACTLY the action's (beta,kappa)-conjugate energies so downstream
// reweighting (reweight.hpp) is exact:
//   weight=exp(-S), S=beta*A-kappa*B+on-site; A=sum_plaq(1-cos), B=hopping sum.
// No physics scan here -- only machine-precision algebraic identities + gauge invariance.
#include "check.hpp"
#include "u1/u1.hpp"
#include "u1/scan_obs.hpp"
#include "u1/monopole.hpp"   // monopole_density<D> (the new ts column)
#include <random>
#include <vector>
#include <cstdio>
#include <cstring>
#include <string>
#include <sstream>
#include <cmath>

using namespace gh;
using namespace gh::u1;

template <int D>
static std::array<int, D> cube(int n) { std::array<int, D> L{}; for (int mu = 0; mu < D; ++mu) L[mu] = n; return L; }

// Independent re-implementation of B = sum_{x,mu} 2 Re[ conj(phi_x) e^{iq theta} phi_{x+mu} ].
template <int D>
static Real hop_sum_ref(const std::vector<Complex>& phi, const std::vector<Real>& th,
                        const Lattice<D>& lat, int q) {
  Real s = 0.0;
  for (std::int64_t x = 0; x < lat.vol; ++x)
    for (int mu = 0; mu < D; ++mu) {
      const std::int64_t y = lat.neighbor_fwd(x, mu);
      const Complex u = std::exp(Complex(0.0, static_cast<Real>(q) * th[x * D + mu]));
      s += 2.0 * (std::conj(phi[x]) * u * phi[y]).real();
    }
  return s;
}

int main() {
  std::mt19937_64 rng(20260529ULL); std::normal_distribution<Real> g(0, 1);
  Lattice<4> lat(cube<4>(4));
  std::vector<Real> th(lat.vol * 4); for (auto& t : th) t = 0.6 * g(rng);
  std::vector<Complex> phi(lat.vol); for (auto& z : phi) z = Complex(0.7 * g(rng), 0.7 * g(rng));

  // ---- 1. CONVENTION (machine precision) ----
  std::printf("-- convention: A,B are the exact (beta,kappa)-conjugates --\n");
  // A == gauge_action(.,.,1.0).
  CHECK_CLOSE(plaq_energy_sum<4>(th, lat), gauge_action<4>(th, lat, 1.0), 1e-12,
              "A = plaq_energy_sum == gauge_action(th,lat,1.0)");
  // A is exactly the beta-conjugate: gauge_action(.,.,beta) == beta*A for several beta.
  for (Real beta : {0.3, 1.0, 1.7, 2.5, 4.0}) {
    char m[96]; std::snprintf(m, sizeof m, "gauge_action(beta=%.2f) == beta*A", beta);
    CHECK_CLOSE(gauge_action<4>(th, lat, beta), beta * plaq_energy_sum<4>(th, lat), 1e-12, m);
  }
  // B == independent direct re-implementation, for q = 1,2,3.
  for (int q = 1; q <= 3; ++q) {
    char m[64]; std::snprintf(m, sizeof m, "B = hop_energy_sum == ref (q=%d)", q);
    CHECK_CLOSE(hop_energy_sum<4>(phi, th, lat, q), hop_sum_ref<4>(phi, th, lat, q), 1e-12, m);
  }

  // ---- 2. EXTENSIVITY ties to the averaged forms ----
  std::printf("-- extensivity ties --\n");
  const std::int64_t n_plaq  = lat.vol * 4 * (4 - 1) / 2;   // = lat.n_plaq()
  const std::int64_t n_links = lat.vol * 4;
  CHECK_CLOSE(static_cast<Real>(n_plaq), static_cast<Real>(lat.n_plaq()), 1e-12, "N_plaq = vol*D(D-1)/2");
  // A / N_plaq == 1 - <plaq>.
  CHECK_CLOSE(plaq_energy_sum<4>(th, lat) / static_cast<Real>(n_plaq),
              1.0 - avg_plaquette<4>(th, lat), 1e-12, "A/N_plaq == 1 - avg_plaquette");
  // B == 2 * N_links * link_energy  (B has the 2Re; link_energy is Re averaged over vol*D).
  for (int q = 1; q <= 3; ++q) {
    char m[80]; std::snprintf(m, sizeof m, "B == 2*N_links*link_energy (q=%d)", q);
    CHECK_CLOSE(hop_energy_sum<4>(phi, th, lat, q),
                2.0 * static_cast<Real>(n_links) * link_energy<4>(phi, th, lat, q), 1e-12, m);
  }

  // ---- 3. GAUGE INVARIANCE ----
  std::printf("-- gauge invariance --\n");
  for (int q = 1; q <= 3; ++q) {
    std::vector<Real> a(lat.vol); for (auto& v : a) v = g(rng);
    std::vector<Real> th2 = th;
    for (std::int64_t x = 0; x < lat.vol; ++x)
      for (int mu = 0; mu < 4; ++mu) th2[x * 4 + mu] += a[x] - a[lat.neighbor_fwd(x, mu)];
    std::vector<Complex> phi2(lat.vol);
    for (std::int64_t x = 0; x < lat.vol; ++x) phi2[x] = std::polar(1.0, static_cast<Real>(q) * a[x]) * phi[x];
    char m[80];
    std::snprintf(m, sizeof m, "gauge-inv A (q=%d)", q);
    CHECK_CLOSE(plaq_energy_sum<4>(th2, lat), plaq_energy_sum<4>(th, lat), 1e-10, m);
    std::snprintf(m, sizeof m, "gauge-inv B (q=%d)", q);
    CHECK_CLOSE(hop_energy_sum<4>(phi2, th2, lat, q), hop_energy_sum<4>(phi, th, lat, q), 1e-10, m);
  }

  // ---- 4. DRIVER SMOKE: ts carries a finite monopole_density column ----
  // Mirror u1_scan.cpp's per-trajectory measurement loop on a TINY 4^4 lattice (no physics
  // scan): run a few trajectories, write the ts file with the SAME header + row layout as
  // the driver, then re-parse the '# columns:' line (header-driven) and confirm a finite
  // monopole_density column is present in every data row.
  std::printf("-- driver smoke: ts monopole_density column --\n");
  {
    const char* tspath = "/tmp/test_u1_scan_mono_smoke.dat";
    const int q = 2, nmeas = 6;
    u1::U1HMC<4> hmc(cube<4>(4), 12345ULL);
    hmc.beta = 1.0; hmc.kappa = 0.2; hmc.lambda = 0.5; hmc.q = q; hmc.tau = 1.0; hmc.nmd = 8;
    hmc.hot(0.8); hmc.cold_phi(0.5);
    for (int t = 0; t < 5; ++t) hmc.trajectory();   // tiny thermalization

    FILE* tf = std::fopen(tspath, "w");
    CHECK(tf != nullptr, "smoke: opened ts file for writing");
    if (tf) {
      // EXACT header + column line the driver writes (must contain monopole_density).
      std::fprintf(tf, "# U(1)+charge-%d Higgs time series (test smoke).\n", q);
      std::fprintf(tf, "# columns: traj  A  B  avg_plaquette  higgs_length  link_energy  monopole_density\n");
      for (int t = 0; t < nmeas; ++t) {
        hmc.trajectory();
        const Real A   = u1::plaq_energy_sum<4>(hmc.th, hmc.lat);
        const Real B   = u1::hop_energy_sum<4>(hmc.phi, hmc.th, hmc.lat, q);
        const Real pl  = u1::avg_plaquette<4>(hmc.th, hmc.lat);
        const Real lp  = u1::higgs_length<4>(hmc.phi, hmc.lat);
        const Real le  = u1::link_energy<4>(hmc.phi, hmc.th, hmc.lat, q);
        const Real rho = u1::monopole_density<4>(hmc.th, hmc.lat);
        std::fprintf(tf, "%d %.15g %.15g %.15g %.15g %.15g %.15g\n", t, A, B, pl, lp, le, rho);
      }
      std::fclose(tf);
    }

    // Re-parse: find the '# columns:' line, map NAME -> index, locate monopole_density,
    // then verify every data row has that many columns and a finite value there.
    FILE* rf = std::fopen(tspath, "r");
    CHECK(rf != nullptr, "smoke: reopened ts file for reading");
    if (rf) {
      char line[1024];
      int mono_idx = -1, ncols = 0;
      const char* tag = "# columns:";
      while (std::fgets(line, sizeof line, rf)) {
        if (line[0] == '#') {
          if (std::strncmp(line, tag, std::strlen(tag)) == 0) {
            // Tokenize names after the tag; record index of monopole_density.
            std::istringstream is(std::string(line + std::strlen(tag)));
            std::string name; int idx = 0;
            while (is >> name) { if (name == "monopole_density") mono_idx = idx; ++idx; }
            ncols = idx;
          }
          continue;
        }
        break;  // first data row (already in 'line')
      }
      CHECK(mono_idx >= 0, "smoke: '# columns:' line contains monopole_density");
      CHECK(ncols == 7, "smoke: header lists 7 columns (5 legacy + traj + monopole_density)");

      // 'line' currently holds the first data row; validate it and the rest.
      int nrows = 0, ok_rows = 0;
      do {
        if (line[0] == '#' || line[0] == '\n' || line[0] == '\0') continue;
        std::istringstream ds(line);
        std::vector<double> vals; double v;
        while (ds >> v) vals.push_back(v);
        ++nrows;
        const bool good = (static_cast<int>(vals.size()) == ncols) &&
                          (mono_idx < static_cast<int>(vals.size())) &&
                          std::isfinite(vals[mono_idx]) && (vals[mono_idx] >= 0.0);
        if (good) ++ok_rows;
      } while (std::fgets(line, sizeof line, rf));
      std::fclose(rf);
      CHECK(nrows == nmeas, "smoke: read all data rows");
      CHECK(ok_rows == nrows, "smoke: every data row has a finite, non-negative monopole_density");
    }
  }

  return report("test_u1_scan");
}
