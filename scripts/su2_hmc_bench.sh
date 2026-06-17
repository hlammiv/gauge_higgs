#!/usr/bin/env bash
# su2_hmc_bench.sh -- profile the SU(2) frozen-Higgs HMC (GaugeHiggsHMC<4,2>, driver
# hmc_higgs_multi) and attribute per-trajectory cost to MD-force components.
#
# Two builds, from an UNMODIFIED copy of src/hmc_higgs_multi.cpp:
#   * clean    : plain -O3, used for honest wall/user-CPU timing of cost(nmd), cost(L).
#   * profiled : -DGH_PROFILE, the src/core/profile.hpp scope timers fire and the driver
#                (with one added GH_PROF_REPORT() before its final return) prints a
#                per-component breakdown table. Ratios are robust to machine load.
#
# Cost model:  cost_per_traj(nmd) = a + b*nmd  with b = ONE Omelyan kick (= one MD force
# eval: build_Dcache + add_gauge_force + add_matter_link_force + scalar_force). The
# Omelyan integrator does (2*nmd+1) kicks and 2*nmd drifts per trajectory; a is the
# per-traj fixed cost (2 hamiltonian() action evals + momentum refresh + save/restore).
#
# RIGOR: measure USER CPU (%U), not elapsed -- this box is shared/loaded, so elapsed is
# contended. Always OMP_NUM_THREADS=1 (single-thread per-component cost) unless -t given.
# Keep jobs small (L<=8, few trajs); the production machines are busy.
#
# This script does NOT modify production src/. It builds into a scratch dir.
#
# Usage:
#   scripts/su2_hmc_bench.sh build           # build clean + profiled scratch drivers
#   scripts/su2_hmc_bench.sh nmd   [L] [rep]  # cost(nmd) at fixed L,rep  -> fit a+b*nmd
#   scripts/su2_hmc_bench.sh rep   [L] [nmd]  # cost vs rep dim d (fund,6,8,12)
#   scripts/su2_hmc_bench.sh vol   [rep][nmd] # cost vs L (4,6,8)
#   scripts/su2_hmc_bench.sh prof  [L][rep][nmd]   # one profiled component-breakdown run
#   scripts/su2_hmc_bench.sh frozen [L][rep][nmd]  # frozen vs unfrozen comparison
#   scripts/su2_hmc_bench.sh all              # the full campaign (slow)
#
# Env: SCRATCH (build dir, default /tmp/su2bench), NMEAS (trajs/timing run, default 20),
#      THREADS (OMP, default 1).
set -u
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SCRATCH="${SCRATCH:-/tmp/su2bench}"
NMEAS="${NMEAS:-20}"
THREADS="${THREADS:-1}"
# spin-3 (2T/BT) couplings from the campaign example; rep `6` needs 7 of them. Reps with
# a different channel count fall back to "auto" (all f_c=1) -- cost is what we measure,
# not physics, so the coupling values do not matter for timing.
COUP_SPIN3="0.1287,0.1548,0.1835,0.2399,0.0056,0.1745,0.1130"
FROZEN="${FROZEN:-1}"   # 1 -> GH_FROZEN set

build() {
  mkdir -p "$SCRATCH"
  cp "$REPO/src/hmc_higgs_multi.cpp" "$SCRATCH/drv.cpp"
  # add the profile report before main's final return (idempotent)
  if ! grep -q GH_PROF_REPORT "$SCRATCH/drv.cpp"; then
    perl -0pi -e 's/(\n\s*)return 0;\n\}\s*$/$1GH_PROF_REPORT();$1return 0;\n}\n/' "$SCRATCH/drv.cpp"
  fi
  local F="-std=c++20 -O3 -march=native -funroll-loops -fopenmp -I$REPO/src -DNDIM=4 -DNCOL=2"
  g++ $F            -o "$SCRATCH/hmc_clean" "$SCRATCH/drv.cpp" && echo "built $SCRATCH/hmc_clean"
  g++ $F -DGH_PROFILE -o "$SCRATCH/hmc_prof"  "$SCRATCH/drv.cpp" && echo "built $SCRATCH/hmc_prof"
}

# couplings string for a given rep: spin-3 uses the campaign values, else "auto".
coup_for() { [ "$1" = "6" ] && echo "$COUP_SPIN3" || echo "auto"; }
frozen_env() { [ "$FROZEN" = "1" ] && echo "GH_FROZEN=1" || echo ""; }

# user-CPU seconds for nmeas trajectories (ntherm=0): prints "user_s"
time_run() { # rep L beta kappa mu2 nmd
  local rep="$1" L="$2" beta="$3" kappa="$4" mu2="$5" nmd="$6"
  local coup; coup="$(coup_for "$rep")"
  /usr/bin/time -f "%U" env $(frozen_env) OMP_NUM_THREADS="$THREADS" \
      "$SCRATCH/hmc_clean" "$rep" "$L" "$beta" "$kappa" "$mu2" "$coup" 0 "$NMEAS" "$nmd" 1 1 \
      >/dev/null 2>"$SCRATCH/.t"; tail -1 "$SCRATCH/.t"
}

cmd_nmd() {
  local L="${1:-8}" rep="${2:-6}"
  echo "# cost(nmd): rep=$rep L=$L frozen=$FROZEN nmeas=$NMEAS  -> user_s/traj"
  echo "# nmd  user_s  s_per_traj"
  for nmd in 4 8 16 32; do
    local u; u="$(time_run "$rep" "$L" 1.5 1.0 0.113 "$nmd")"
    awk -v n="$nmd" -v u="$u" -v m="$NMEAS" 'BEGIN{printf "%d %.2f %.4f\n", n, u, u/m}'
  done
}

cmd_rep() {
  local L="${1:-6}" nmd="${2:-8}"
  echo "# cost vs rep: L=$L nmd=$nmd frozen=$FROZEN nmeas=$NMEAS  -> user_s/traj"
  echo "# rep  d  user_s  s_per_traj"
  for rep in fund 6 8 12; do
    local u; u="$(time_run "$rep" "$L" 1.5 1.0 0.113 "$nmd")"
    local d; d="$("$SCRATCH/hmc_clean" "$rep" 2>/dev/null | sed -n 's/.* d=\([0-9]*\).*/\1/p' | head -1)"
    awk -v r="$rep" -v d="$d" -v u="$u" -v m="$NMEAS" 'BEGIN{printf "%s %s %.2f %.4f\n", r, d, u, u/m}'
  done
}

cmd_vol() {
  local rep="${1:-6}" nmd="${2:-8}"
  echo "# cost vs L: rep=$rep nmd=$nmd frozen=$FROZEN nmeas=$NMEAS  -> user_s/traj"
  echo "# L  vol  user_s  s_per_traj  s_per_traj_per_site"
  for L in 4 6 8; do
    local u; u="$(time_run "$rep" "$L" 1.5 1.0 0.113 "$nmd")"
    awk -v L="$L" -v u="$u" -v m="$NMEAS" \
        'BEGIN{v=L*L*L*L; printf "%d %d %.2f %.4f %.3e\n", L, v, u, u/m, (u/m)/v}'
  done
}

cmd_prof() {
  local L="${1:-8}" rep="${2:-6}" nmd="${3:-10}" nm="${NMEAS_PROF:-12}"
  local coup; coup="$(coup_for "$rep")"
  echo "# profiled breakdown: rep=$rep L=$L nmd=$nmd frozen=$FROZEN nmeas=$nm"
  env $(frozen_env) OMP_NUM_THREADS="$THREADS" \
      "$SCRATCH/hmc_prof" "$rep" "$L" 1.5 1.0 0.113 "$coup" 0 "$nm" "$nmd" 1 1 2>&1 | \
      sed -n '/GH_PROFILE breakdown/,/====/p'
}

cmd_frozen() {
  local L="${1:-8}" rep="${2:-6}" nmd="${3:-8}"
  echo "# frozen vs unfrozen: rep=$rep L=$L nmd=$nmd nmeas=$NMEAS"
  for fz in 1 0; do
    FROZEN="$fz"; local u; u="$(time_run "$rep" "$L" 1.5 1.0 0.113 "$nmd")"
    echo "frozen=$fz  user_s=$u  s_per_traj=$(awk -v u="$u" -v m="$NMEAS" 'BEGIN{printf "%.4f",u/m}')"
  done
}

case "${1:-help}" in
  build)  build ;;
  nmd)    shift; cmd_nmd  "$@" ;;
  rep)    shift; cmd_rep  "$@" ;;
  vol)    shift; cmd_vol  "$@" ;;
  prof)   shift; cmd_prof "$@" ;;
  frozen) shift; cmd_frozen "$@" ;;
  all)    build; cmd_nmd 8 6; cmd_rep 6 8; cmd_vol 6 8; cmd_prof 8 6 10; cmd_frozen 8 6 8 ;;
  *) sed -n '2,40p' "$0" ;;
esac
