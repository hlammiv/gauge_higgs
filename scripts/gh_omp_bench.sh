#!/usr/bin/env bash
# OpenMP strong-scaling benchmark for the frozen SU(2) gauge-Higgs HMC (hmc_higgs_multi).
# Holds (rep, L, beta, kappa, mu2, couplings, nmd, ntherm, nmeas, seed) fixed and varies
# ONLY OMP_NUM_THREADS, reporting wall-clock s/traj. The per-traj cost is ~100% MD force
# evaluation (linear in nmd), all three forces + build_Dcache carry #pragma omp parallel for,
# so this measures the realizable intra-job (latency) speedup of running at OMP=k.
#
# Usage: scripts/gh_omp_bench.sh <rep> <L> <nmd> <ntherm> <nmeas> "<thread list>"
#   e.g. scripts/gh_omp_bench.sh 6  8 10 4 20 "1 2 4 8"     # BT spin-3 d=7
#        scripts/gh_omp_bench.sh 12 4 8  4 20 "1 2 4 8"     # BI 2I spin-6 d=13
set -euo pipefail
BIN=./build/hmc_higgs_multi
REP=${1:-6}; L=${2:-8}; NMD=${3:-10}; NTH=${4:-4}; NMEAS=${5:-20}; THREADS=${6:-"1 2 4 8"}
BETA=1.5; KAPPA=1.0; MU2=0.113
# Determine #quartic channels for this rep from the driver's header line (generic
# non-degenerate coupling list of that length keeps the potential non-trivial).
NCH=$($BIN "$REP" 2>&1 | grep -oP '\K[0-9]+(?= quartic channels)' | head -1)
COUP=$(python3 -c "import sys;n=int(sys.argv[1]);print(','.join(['%.4f'%(0.10+0.013*i) for i in range(n)]))" "$NCH")
echo "# rep=$REP L=$L nmd=$NMD ntherm=$NTH nmeas=$NMEAS  beta=$BETA kappa=$KAPPA mu2=$MU2  nchan=$NCH"
echo "# OMP   s/traj   speedup-vs-OMP1"
base=""
for T in $THREADS; do
  t0=$(date +%s.%N)
  OMP_NUM_THREADS=$T GH_FROZEN=1 $BIN "$REP" "$L" "$BETA" "$KAPPA" "$MU2" "$COUP" "$NTH" "$NMEAS" "$NMD" 1 1 >/dev/null 2>&1
  t1=$(date +%s.%N)
  ntot=$((NTH + NMEAS))
  spt=$(python3 -c "print((($t1)-($t0))/$ntot)")
  if [ -z "$base" ]; then base=$spt; fi
  spd=$(python3 -c "print('%.2f'%($base/$spt))")
  printf "  %-4s %8.4f   %sx\n" "$T" "$spt" "$spd"
done
