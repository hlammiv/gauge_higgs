#!/usr/bin/env bash
# Profile / benchmark GaugeHiggsHMC<4,2> (SU(2) frozen-Higgs, hmc_higgs_multi).
# Builds a GH_PROFILE driver that prints the per-component MD-force timer table
# (kick / add_gauge_force / add_matter_link_force / scalar_force / rotate / fast_D ...)
# and a plain driver for clean wall-clock. Use to attribute per-traj cost and test
# the two levers: (1) integrator nmd vs acceptance, (2) OpenMP thread scaling.
#
# Usage:
#   scripts/gh_hmc_profile.sh build            # build /tmp/ghprof/{driver_prof,driver_plain}
#   scripts/gh_hmc_profile.sh prof  <rep L beta kappa mu2 coup ntherm nmeas nmd>   # timer table
#   scripts/gh_hmc_profile.sh omp   <rep L beta kappa mu2 coup ntherm nmeas nmd>   # OMP=1,2,4,8 wall
#   scripts/gh_hmc_profile.sh nmd   <rep L beta kappa mu2 coup>                    # accept vs nmd
# rep: 2(fund) | 6(spin3=BT) | 8(spin4) | 12(spin6=2I). coup = comma f_c list or 'auto'.
set -u
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SCRATCH=/tmp/ghprof
build(){
  mkdir -p "$SCRATCH"
  cp "$REPO/src/hmc_higgs_multi.cpp" "$SCRATCH/driver_prof.cpp"
  python3 - "$SCRATCH/driver_prof.cpp" <<'PY'
import sys; p=sys.argv[1]; s=open(p).read()
s=s.replace('#include "rep/rep_general.hpp"','#include "rep/rep_general.hpp"\n#include "core/profile.hpp"',1)
s=s.replace('hmc.acceptance(), sExp / nmeas);\n  return 0;','hmc.acceptance(), sExp / nmeas);\n  GH_PROF_REPORT();\n  return 0;')
open(p,'w').write(s)
PY
  g++ -std=c++17 -O3 -march=native -funroll-loops -fopenmp -DNDIM=4 -DNCOL=2 -I"$REPO/src" -o "$SCRATCH/driver_plain" "$SCRATCH/driver_prof.cpp"
  g++ -std=c++17 -O3 -march=native -funroll-loops -fopenmp -DNDIM=4 -DNCOL=2 -DGH_PROFILE -I"$REPO/src" -o "$SCRATCH/driver_prof" "$SCRATCH/driver_prof.cpp"
  echo "built $SCRATCH/driver_prof $SCRATCH/driver_plain"
}
cmd="${1:-build}"; shift || true
case "$cmd" in
  build) build ;;
  prof)  GH_FROZEN=1 OMP_NUM_THREADS=1 "$SCRATCH/driver_prof" "$@" 1 2>&1 | grep -E "plaq|accept|====|slot|rotate|kick|gauge|matter|scalar|fast_D|build_Dcache|dmat" ;;
  omp)   for T in 1 2 4 8; do S=$(date +%s.%N); GH_FROZEN=1 OMP_NUM_THREADS=$T "$SCRATCH/driver_plain" "$@" 1 >/dev/null 2>&1; E=$(date +%s.%N); echo "OMP=$T wall=$(echo "$E-$S"|bc)s"; done ;;
  nmd)   for n in 4 8 16 32; do echo -n "nmd=$n "; GH_FROZEN=1 OMP_NUM_THREADS=4 "$SCRATCH/driver_plain" "$@" 8 40 $n 1 1 2>&1 | grep accept; done ;;
  *) echo "unknown: $cmd"; exit 1 ;;
esac
